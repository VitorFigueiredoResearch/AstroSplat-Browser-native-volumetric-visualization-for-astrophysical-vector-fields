#!/usr/bin/env python3
"""
solar_to_splats.py
══════════════════════════════════════════════════════════════════════════════
SDO/HMI Active Region — Field Line Structure JSON Converter
─────────────────────────────────────────────────────────────────────────────
PIPELINE:
  Mode A · JSOC/DRMS   — sunpy + drms → JSOC HMI magnetogram → FITS → JSON
  Mode B · HTTP FITS   — direct urllib retrieval from a known public FITS URL
  Mode C · SYNTHETIC   — analytically generated bipolar active region

PHYSICS (Potential Field Extrapolation):
  The photospheric line-of-sight magnetic field B_z is analytically
  continued into the low corona assuming a current-free (potential) field.
  The governing equations are:
      B = -∇φ,   ∇²φ = 0  (Laplace equation)
  The Fourier-domain solution (Green's function approach, Alissandrakis 1981)
  yields at z = 0 (photospheric boundary):
      φ̂(kx, ky) = B̂_z(kx, ky) / k
      B̂_x(k)    = -i · (kx/k) · B̂_z(k)
      B̂_y(k)    = -i · (ky/k) · B̂_z(k)
  where k = √(kx² + ky²) is the horizontal wavenumber magnitude.
  Field components are recovered by inverse FFT (Nakagawa & Raadu 1972).

OUTPUT:
  solar_active_region.json — compatible with field_visualizer.html
  Each record: { x, y, z, Bx, By, Bz, intensity, gray }
  The first record carries metadata only (source = "meta").
══════════════════════════════════════════════════════════════════════════════
USAGE:
  python solar_to_splats.py                    # cascading auto mode
  python solar_to_splats.py --mode synthetic   # guaranteed offline result
  python solar_to_splats.py --mode real        # requires sunpy + JSOC email
  python solar_to_splats.py --mode http        # public FITS, no authentication
  python solar_to_splats.py --ar 13590         # NOAA active region number
  python solar_to_splats.py --date 2024-03-14  # observation date

DEPENDENCIES (Mode A):
  pip install sunpy drms astropy scipy --break-system-packages
  Email registration: http://jsoc.stanford.edu/ajax/register_email.html
"""

import argparse
import gc
import json
import math
import os
import struct
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

import numpy as np

# ════════════════════════════════════════════════════════════════════════════
#  CONFIGURATION
# ════════════════════════════════════════════════════════════════════════════
JSOC_EMAIL    = "your.email@example.com"   # JSOC registered email address
NOAA_AR       = 13590                       # Default NOAA Active Region number
DATE_STR      = "2024-01-15"               # Observation date [ISO 8601]
GRID_SIZE     = 60                          # Output grid resolution: 60×60 = 3,600 points
MAX_POINTS    = 5000                        # Absolute point count ceiling
Z_SCALE       = 0.012                       # Bz [G] → scene Z-coordinate scaling factor
OUTPUT_FILE   = "solar_active_region.json"

# Known public FITS endpoints (SDO/HMI science archive mirrors).
# Authentication is not required for these URLs.
PUBLIC_FITS_URLS = [
    "https://hesperia.gsfc.nasa.gov/sftp/pub/hmi/magnetogram/hmi_m_2024_01_15.fits",
    "https://jsoc.stanford.edu/data/hmi/fits/2014/01/01/HMI20140101_000000_45s.M.fits",
]

RESET = "\033[0m"; RED  = "\033[91m"; YEL  = "\033[93m"
GREEN = "\033[92m"; CYAN = "\033[96m"; BOLD = "\033[1m"
DIM   = "\033[2m"


# ════════════════════════════════════════════════════════════════════════════
#  MINIMAL FITS READER  (stdlib only — astropy not required for Modes B/C)
# ════════════════════════════════════════════════════════════════════════════
class MiniFITS:
    """
    Partial FITS reader supporting single-extension image HDUs.

    BITPIX values supported: 8, 16, 32, 64, -32, -64.
    BSCALE and BZERO scaling are applied. Data are returned as float32.
    This implementation is sufficient for HMI line-of-sight magnetogram files.

    References
    ----------
    Wells, D.C. et al. (1981). FITS: a flexible image transport system.
    A&A Supplement Series, 44, 363–370.
    """
    DTYPE_MAP = {8: ">u1", 16: ">i2", 32: ">i4", 64: ">i8",
                 -32: ">f4", -64: ">f8"}

    def __init__(self, path_or_bytes):
        if isinstance(path_or_bytes, (str, Path)):
            with open(path_or_bytes, "rb") as f:
                raw = f.read()
        else:
            raw = path_or_bytes
        self.header, self.data = self._parse(raw)

    def _parse(self, raw: bytes):
        BLOCK = 2880
        pos   = 0
        header = {}
        while True:
            block = raw[pos:pos+BLOCK]
            pos  += BLOCK
            for i in range(36):
                card = block[i*80:(i+1)*80].decode("ascii", errors="replace")
                key  = card[:8].strip()
                if key == "END":
                    break
                if "=" in card:
                    k, v = card.split("=", 1)
                    v = v.split("/")[0].strip().strip("'").strip()
                    header[k.strip()] = v
            if "END" in block.decode("ascii", errors="replace"):
                break

        bitpix = int(header.get("BITPIX", -32))
        naxis  = int(header.get("NAXIS",   2))
        axes   = [int(header.get(f"NAXIS{i+1}", 1)) for i in range(naxis)]
        bscale = float(header.get("BSCALE", 1.0))
        bzero  = float(header.get("BZERO",  0.0))

        dtype  = np.dtype(self.DTYPE_MAP.get(bitpix, ">f4"))
        n_elem = 1
        for a in axes: n_elem *= a
        data_raw = np.frombuffer(raw[pos:pos + n_elem * dtype.itemsize], dtype=dtype)
        data     = data_raw.astype(np.float32) * bscale + bzero

        if naxis >= 2:
            data = data.reshape(axes[1], axes[0])   # (NAXIS2, NAXIS1) → (rows, cols)

        return header, data


# ════════════════════════════════════════════════════════════════════════════
#  POTENTIAL FIELD EXTRAPOLATION  (FFT / Green's function method)
# ════════════════════════════════════════════════════════════════════════════
def potential_field_photosphere(bz_2d: np.ndarray, pixel_size: float = 1.0):
    """
    Derive horizontal field components at the photospheric boundary from
    the observed vertical field B_z using the FFT-based Green's function
    approach (Alissandrakis 1981; Nakagawa & Raadu 1972).

    The potential field satisfies B = -∇φ with ∇²φ = 0. The Fourier-domain
    solution at the lower boundary (z = 0) is:
        φ̂(kx, ky) = B̂_z(kx, ky) / k
        B̂_x       = -i · kx · φ̂
        B̂_y       = -i · ky · φ̂

    The DC wavenumber component (k = 0) is regularised to zero, consistent
    with the assumption of zero net magnetic flux over the domain.

    Parameters
    ----------
    bz_2d : np.ndarray, shape (ny, nx)
        Photospheric B_z field [G].
    pixel_size : float
        Physical pixel scale used to compute wavenumbers [arbitrary units].

    Returns
    -------
    bx, by : np.ndarray
        Horizontal field components at z = 0 [same units as bz_2d].
    """
    ny, nx = bz_2d.shape

    bz_hat = np.fft.fft2(bz_2d)

    # Horizontal wavenumber arrays [rad / pixel_size]
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=pixel_size)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=pixel_size)
    KX, KY = np.meshgrid(kx, ky)
    K      = np.sqrt(KX**2 + KY**2)
    K[0, 0] = 1e-12            # regularise DC component; net flux assumed zero

    # Horizontal components in Fourier space
    bx_hat = -1j * (KX / K) * bz_hat
    by_hat = -1j * (KY / K) * bz_hat

    bx = np.real(np.fft.ifft2(bx_hat))
    by = np.real(np.fft.ifft2(by_hat))

    return bx, by


def bz_to_z_position(bz_2d: np.ndarray, scale: float = Z_SCALE) -> np.ndarray:
    """
    Map photospheric B_z to a 3D height coordinate above the photosphere.

    Positive (outward) polarity is mapped to positive Z (upward);
    negative polarity to negative Z (subsurface, for visual representation).
    A tanh saturation function prevents extreme field strengths from
    producing numerically unreasonable scene coordinates.

    Parameters
    ----------
    bz_2d : np.ndarray
        Photospheric vertical field component [G].
    scale : float
        Linear scaling factor [scene units / G].

    Returns
    -------
    np.ndarray : Z-coordinate array [scene units].
    """
    bz_norm = bz_2d / (np.percentile(np.abs(bz_2d[np.isfinite(bz_2d)]), 99) + 1e-9)
    return np.tanh(bz_norm * 1.8) * (GRID_SIZE * scale * 4.0)


# ════════════════════════════════════════════════════════════════════════════
#  DATA PROCESSING UTILITIES
# ════════════════════════════════════════════════════════════════════════════
def downsample_region(data: np.ndarray, grid: int = GRID_SIZE) -> np.ndarray:
    """
    Identify the most magnetically active sub-region and bin-average it
    to the target grid resolution.

    The activity criterion is the total squared field strength (ΣB²) within
    a sliding window of size grid×grid. A coarse 16-step scan is applied
    in both dimensions to reduce computational cost. NaN and Inf pixels
    (common in FITS magnetograms at array boundaries) are replaced with zero.

    Parameters
    ----------
    data : np.ndarray
        Full-resolution input magnetogram [G].
    grid : int
        Target output resolution (grid × grid pixels).

    Returns
    -------
    np.ndarray, shape (grid, grid) : Downsampled sub-region [G].
    """
    ny, nx = data.shape
    data   = np.where(np.isfinite(data), data, 0.0)

    best_score       = -1
    best_iy, best_ix = 0, 0

    step_y = max(1, (ny - grid) // 16)
    step_x = max(1, (nx - grid) // 16)
    for iy in range(0, max(1, ny - grid), step_y):
        for ix in range(0, max(1, nx - grid), step_x):
            patch = data[iy:iy+grid, ix:ix+grid]
            score = float(np.sum(patch**2))
            if score > best_score:
                best_score       = score
                best_iy, best_ix = iy, ix

    cropped = data[best_iy:best_iy+grid, best_ix:best_ix+grid]

    if cropped.shape == (grid, grid):
        return cropped

    # Bin-average when the cropped region exceeds the target dimensions
    factor_y = cropped.shape[0] / grid
    factor_x = cropped.shape[1] / grid
    result   = np.zeros((grid, grid), dtype=np.float32)
    for iy in range(grid):
        for ix in range(grid):
            y0 = int(iy * factor_y);  y1 = int((iy+1) * factor_y) or y0+1
            x0 = int(ix * factor_x);  x1 = int((ix+1) * factor_x) or x0+1
            result[iy, ix] = cropped[y0:y1, x0:x1].mean()
    return result


def gray_from_intensity(intensity: np.ndarray, gamma: float = 0.55) -> np.ndarray:
    """
    Map absolute field intensity to an 8-bit grayscale luminance index.

    Intensity values are normalised to [0, 1] using the 99th percentile
    as the saturation threshold, then gamma-compressed with exponent γ = 0.55
    to improve perceptual uniformity across the dynamic range.

    Parameters
    ----------
    intensity : np.ndarray
        Absolute field magnitude |B| [G].
    gamma : float
        Compression exponent; values < 1 expand the mid-range.

    Returns
    -------
    np.ndarray (int) : Luminance indices in [0, 255].
    """
    b_max = np.percentile(intensity[np.isfinite(intensity)], 99) + 1e-9
    norm  = np.clip(intensity / b_max, 0, 1)
    return (np.power(norm, gamma) * 255).astype(int)


def build_json_records(bz_grid: np.ndarray,
                       bx_grid: np.ndarray,
                       by_grid: np.ndarray,
                       source_label: str = "unknown") -> list:
    """
    Construct field line records from 2D field component arrays.

    Each grid pixel is converted to one 3D record. Spatial coordinates
    are normalised to [-5, +5] scene units. Field vectors are normalised
    by their respective component maxima to produce dimensionless direction
    cosines. Points with total field magnitude < 0.5 G are excluded
    as below-threshold dead pixels.

    Parameters
    ----------
    bz_grid, bx_grid, by_grid : np.ndarray, shape (GRID_SIZE, GRID_SIZE)
        Vertical and horizontal field components [G].
    source_label : str
        Provenance identifier embedded in each record.

    Returns
    -------
    list[dict] : Field line records, each containing x, y, z, Bx, By, Bz,
                 intensity [G], gray [0–255], and source.
    """
    ny, nx  = bz_grid.shape
    z_grid  = bz_to_z_position(bz_grid)
    B_mag   = np.sqrt(bx_grid**2 + by_grid**2 + bz_grid**2)
    gray_g  = gray_from_intensity(np.abs(bz_grid))    # grayscale driven by |Bz|

    bxy_max = max(np.abs(bx_grid).max(), np.abs(by_grid).max(), 1.0)
    bz_max  = max(np.abs(bz_grid).max(), 1.0)

    records = []
    for iy in range(ny):
        for ix in range(nx):
            bz_val = float(bz_grid[iy, ix])
            bx_val = float(bx_grid[iy, ix])
            by_val = float(by_grid[iy, ix])
            z_val  = float(z_grid[iy, ix])
            mag    = float(B_mag[iy, ix])
            g      = int(gray_g[iy, ix])

            # Spatial coordinates normalised to scene domain [-5, +5]
            x_norm = (ix / (nx - 1)) * 2.0 - 1.0
            y_norm = (iy / (ny - 1)) * 2.0 - 1.0

            # Sub-threshold pixels are excluded from the output dataset
            if mag < 0.5:
                continue

            records.append({
                "x":         round(x_norm * 5.0,     4),
                "y":         round(y_norm * 5.0,     4),
                "z":         round(z_val,             4),
                "Bx":        round(bx_val / bxy_max,  5),
                "By":        round(by_val / bxy_max,  5),
                "Bz":        round(bz_val / bz_max,   5),
                "intensity": round(mag,                2),
                "gray":      g,
                "source":    source_label,
            })

    # Records exceeding the point ceiling are discarded in descending order
    # of field intensity, retaining the most active field concentrations.
    if len(records) > MAX_POINTS:
        records.sort(key=lambda r: -r["intensity"])
        records = records[:MAX_POINTS]

    return records


# ════════════════════════════════════════════════════════════════════════════
#  MODE A — JSOC/DRMS REAL DATA (requires sunpy + drms + JSOC email)
# ════════════════════════════════════════════════════════════════════════════
def mode_real(noaa_ar: int, date_str: str) -> list | None:
    """
    Retrieve an HMI SHARP (Space-weather HMI Active Region Patch) cutout
    from the JSOC data archive for the specified NOAA active region.

    The hmi.sharp_cea_720s series provides pre-extracted active region patches
    in cylindrical equal-area (CEA) projection, avoiding the need to
    download or process a full 4096×4096 disk magnetogram. The radial
    field component Br is used as the photospheric B_z boundary condition.

    Dependencies: pip install sunpy drms astropy scipy
    """
    print(f"\n{CYAN}[Mode A] Initiating JSOC query via drms...{RESET}")

    try:
        import drms                                          # type: ignore
        from astropy.io import fits as astrofits            # type: ignore
    except ImportError:
        print(f"  {YEL}Required packages not available (sunpy, drms, astropy).{RESET}")
        print(f"  {DIM}Install with: pip install sunpy drms astropy scipy{RESET}")
        return None

    if "your.email" in JSOC_EMAIL:
        print(f"  {YEL}JSOC_EMAIL has not been configured.{RESET}")
        print(f"  {DIM}Registration: http://jsoc.stanford.edu/ajax/register_email.html{RESET}")
        return None

    try:
        print(f"  Connecting to JSOC (registered address: {JSOC_EMAIL})...")
        client = drms.Client(email=JSOC_EMAIL, verbose=False)

        t_str = date_str.replace("-", ".") + "_00:00:00_TAI"
        ds    = f"hmi.sharp_cea_720s[][{t_str}]{{Bp,Bt,Br}}"

        print(f"  Querying series: {ds}")
        keys, segs = client.query(ds, key=["T_REC", "NOAA_AR", "HARPNUM"],
                                  seg=["Bp", "Br"])

        if keys.empty:
            print(f"  {YEL}No records returned for AR{noaa_ar} on {date_str}.{RESET}")
            return None

        mask = keys["NOAA_AR"].astype(str).str.contains(str(noaa_ar), na=False)
        if not mask.any():
            print(f"  {YEL}AR{noaa_ar} absent from results; first available record used.{RESET}")
            mask = pd.Series([True] + [False]*(len(keys)-1), index=keys.index)

        row     = keys[mask].iloc[0]
        seg_row = segs[mask].iloc[0]
        harpnum = row["HARPNUM"]
        print(f"  SHARP #{harpnum} identified for NOAA AR{noaa_ar}.")

        # Radial (Br) component is retrieved as the B_z boundary condition
        br_url = "http://jsoc.stanford.edu" + seg_row["Br"]
        print(f"  Retrieving Br segment: {br_url}")
        with urllib.request.urlopen(br_url, timeout=60) as resp:
            fits_bytes = resp.read()

        hdu   = astrofits.open(fits_bytes)[0]
        bz_2d = hdu.data.astype(np.float32)
        print(f"  {GREEN}Retrieval complete: {bz_2d.shape} pixels.{RESET}")

        bz_ds        = downsample_region(bz_2d, GRID_SIZE)
        bx_ds, by_ds = potential_field_photosphere(bz_ds)

        label = f"HMI_SHARP_AR{noaa_ar}_{date_str}"
        return build_json_records(bz_ds, bx_ds, by_ds, source_label=label)

    except Exception as e:
        print(f"  {RED}Mode A failed: {e}{RESET}")
        return None


# ════════════════════════════════════════════════════════════════════════════
#  MODE B — PUBLIC HTTP FITS RETRIEVAL (no authentication required)
# ════════════════════════════════════════════════════════════════════════════
def mode_http() -> list | None:
    """
    Retrieve a publicly accessible HMI FITS file from a known archive mirror.

    Failure of all configured URLs causes a silent transition to Mode C.
    The internal MiniFITS parser is used; astropy is not required.
    """
    print(f"\n{CYAN}[Mode B] Attempting retrieval from public FITS archive...{RESET}")

    for url in PUBLIC_FITS_URLS:
        try:
            print(f"  → {url[:70]}...")
            req = urllib.request.Request(
                url, headers={"User-Agent": "SolarFieldConverter/1.0"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                fits_bytes = resp.read()

            print(f"  {GREEN}Retrieved {len(fits_bytes)/1024:.0f} KB.{RESET}")
            fits_obj = MiniFITS(fits_bytes)
            bz_2d    = fits_obj.data

            if bz_2d is None or bz_2d.size == 0:
                raise ValueError("FITS data array is empty.")

            print(f"  Input dimensions: {bz_2d.shape}.  "
                  f"Bz range: [{bz_2d.min():.0f}, {bz_2d.max():.0f}] G.")
            bz_ds        = downsample_region(bz_2d, GRID_SIZE)
            bx_ds, by_ds = potential_field_photosphere(bz_ds)

            return build_json_records(bz_ds, bx_ds, by_ds, source_label="HMI_HTTP")

        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
            print(f"  {YEL}Retrieval failed ({type(e).__name__}). Proceeding to next URL.{RESET}")
        except Exception as e:
            print(f"  {YEL}Parse error: {e}{RESET}")

    print(f"  {RED}All configured HTTP sources are unreachable.{RESET}")
    return None


# ════════════════════════════════════════════════════════════════════════════
#  MODE C — SYNTHETIC ACTIVE REGION  (offline, deterministic)
# ════════════════════════════════════════════════════════════════════════════
def generate_synthetic_active_region(seed: int = 42) -> np.ndarray:
    """
    Generate a synthetic bipolar active region magnetogram.

    The model is structured around observational properties of NOAA active
    regions at solar maximum (van Driel-Gesztelyi & Green 2015):

    Components
    ----------
    Primary bipole:
        Two Gaussian flux concentrations representing the leading and following
        sunspot polarities. The leading polarity (positive in the northern
        hemisphere) is modelled as compact (σ ~ 0.6–0.85 scene units) with a
        peak of 1200–1800 G. The following polarity is broader (σ ~ 0.75–1.0)
        and weaker (900–1400 G), consistent with the Joy's Law asymmetry.
    Satellite flux:
        Three to six minor Gaussian concentrations distributed around the
        primary bipole, representing parasitic polarity regions and ephemeral
        flux emergence.
    Polarity Inversion Line (PIL) shear:
        A sigmoid-modulated contribution along the PIL introduces a non-potential
        appearance without requiring a full MHD simulation.
    Diffuse background:
        Broad low-amplitude halo (~40 G) representing the dispersed background
        network field.
    Instrument noise:
        Gaussian noise with σ = 15 G, consistent with the HMI magnetogram
        noise floor under standard observing conditions.
    Apodisation:
        A circular attenuation mask suppresses edge artefacts, consistent with
        the reduced signal reliability near magnetogram boundaries.

    Parameters
    ----------
    seed : int
        Random number generator seed for reproducibility.

    Returns
    -------
    np.ndarray, shape (GRID_SIZE, GRID_SIZE) : Synthetic B_z [G].
    """
    rng    = np.random.default_rng(seed)
    g      = GRID_SIZE
    x      = np.linspace(-4.0, 4.0, g)
    y      = np.linspace(-4.0, 4.0, g)
    XX, YY = np.meshgrid(x, y)

    bz = np.zeros((g, g), dtype=np.float32)

    # ── Primary bipole ────────────────────────────────────────────────────
    pos_cx,   pos_cy  = -1.35,  0.10
    pos_peak          =  rng.uniform(1200, 1800)    # [G]
    pos_sigma         =  rng.uniform(0.60, 0.85)

    neg_cx,   neg_cy  =  1.50, -0.05
    neg_peak          = -rng.uniform(900, 1400)     # [G]
    neg_sigma         =  rng.uniform(0.75, 1.00)    # broader, consistent with flux dispersal

    bz += pos_peak * np.exp(-((XX - pos_cx)**2 + (YY - pos_cy)**2) / (2*pos_sigma**2))
    bz += neg_peak * np.exp(-((XX - neg_cx)**2 + (YY - neg_cy)**2) / (2*neg_sigma**2))

    # ── Satellite flux concentrations ─────────────────────────────────────
    n_satellites = rng.integers(3, 7)
    for _ in range(n_satellites):
        sc_x    = rng.uniform(-3.0,  3.0)
        sc_y    = rng.uniform(-2.5,  2.5)
        sc_peak = rng.choice([-1, 1]) * rng.uniform(80, 350)    # [G]
        sc_sig  = rng.uniform(0.20, 0.50)
        bz     += sc_peak * np.exp(-((XX - sc_x)**2 + (YY - sc_y)**2) / (2*sc_sig**2))

    # ── PIL shear contribution ────────────────────────────────────────────
    pil_shear = rng.uniform(0.08, 0.25)
    pil_warp  = np.tanh((XX + pil_shear * YY) * 1.5) * rng.uniform(30, 80)
    bz       += pil_warp * np.exp(-(XX**2 + YY**2) / 8.0)

    # ── Diffuse background halo ───────────────────────────────────────────
    bz += rng.uniform(-40, 40) * np.exp(-(XX**2 + YY**2) / 12.0)

    # ── Instrument noise (~15 G RMS, HMI LOS specification) ──────────────
    bz += rng.normal(0, 15.0, (g, g)).astype(np.float32)

    # ── Boundary apodisation (circular mask) ─────────────────────────────
    r2   = (XX / 4.2)**2 + (YY / 4.2)**2
    apod = np.clip(1.0 - r2, 0, 1)**0.5
    bz  *= apod

    return bz


def mode_synthetic(ar_number: int = 0) -> list:
    """Execute the data transformation pipeline using a synthetic active region."""
    print(f"\n{CYAN}[Mode C] Generating synthetic active region dataset...{RESET}")

    seed  = (ar_number % 1000) if ar_number else 42
    bz_2d = generate_synthetic_active_region(seed=seed)

    b_abs = np.abs(bz_2d)
    print(f"  Grid resolution:  {GRID_SIZE}×{GRID_SIZE} = {GRID_SIZE**2} points")
    print(f"  Bz range:         [{bz_2d.min():.0f}, {bz_2d.max():.0f}] G")
    print(f"  |B| mean / p99:   {b_abs.mean():.1f} G / {np.percentile(b_abs,99):.0f} G")
    print(f"  Net flux:         {bz_2d.sum():.0f} G·px  (bipolar domain; expected ~0)")

    bx, by = potential_field_photosphere(bz_2d)
    print(f"  {GREEN}Potential field extrapolation complete "
          f"(FFT Green's function, Alissandrakis 1981).{RESET}")

    label = f"SYNTH_AR{'%05d' % ar_number if ar_number else 'GEN'}"
    return build_json_records(bz_2d, bx, by, source_label=label)


# ════════════════════════════════════════════════════════════════════════════
#  OUTPUT
# ════════════════════════════════════════════════════════════════════════════
def save_json(records: list, out_path: Path, meta: dict):
    """
    Serialise field line records to JSON.

    A single metadata record is prepended to the array. Metadata records
    carry source = "meta" and are filtered by the visualizer at load time.
    All numpy scalar types are cast to native Python types before serialisation.
    """
    output = [{
        "source":    "meta",
        "total_pts": len(records),
        "grid":      f"{GRID_SIZE}x{GRID_SIZE}",
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        **meta
    }] + records

    def _jsonify(obj):
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, (np.integer,)):  return int(obj)
        if isinstance(obj, np.ndarray):     return obj.tolist()
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable.")

    raw = json.dumps(output, separators=(",", ":"), default=_jsonify)
    out_path.write_text(raw, encoding="utf-8")

    kb = len(raw.encode()) / 1024
    print(f"\n{GREEN}{BOLD}Written: {out_path.name}  ({kb:.0f} KB,  {len(records)} records){RESET}")
    print(f"   Dataset is compatible with field_visualizer.html.\n")


# ════════════════════════════════════════════════════════════════════════════
#  VALIDATION REPORT
# ════════════════════════════════════════════════════════════════════════════
def print_field_stats(records: list):
    """Print a summary of dataset statistics for validation."""
    if not records:
        print(f"{RED}No records produced.{RESET}")
        return

    intensities = [r["intensity"] for r in records]
    grays       = [r["gray"]      for r in records]
    zvals       = [r["z"]         for r in records]
    sources     = set(r.get("source", "?") for r in records)

    print(f"\n{BOLD}{'─'*52}")
    print(f"  DATASET VALIDATION REPORT")
    print(f"{'─'*52}{RESET}")
    print(f"  Records    : {len(records)}")
    print(f"  Source     : {', '.join(sources)}")
    print(f"  |B| range  : [{min(intensities):.1f}, {max(intensities):.1f}] G  "
          f"(mean {sum(intensities)/len(intensities):.1f} G)")
    print(f"  Gray range : [{min(grays)}, {max(grays)}]  (expected: 0–255)")
    print(f"  Z range    : [{min(zvals):.2f}, {max(zvals):.2f}] scene units")
    print(f"{'─'*52}")


# ════════════════════════════════════════════════════════════════════════════
#  ENTRY POINT
# ════════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description="SDO/HMI Active Region — Field Line Structure JSON Converter",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("--mode",  choices=["auto", "real", "http", "synthetic"],
                        default="auto",
                        help="Data source mode (auto cascades: real → http → synthetic)")
    parser.add_argument("--ar",    type=int,  default=NOAA_AR,    help="NOAA active region number")
    parser.add_argument("--date",  type=str,  default=DATE_STR,   help="Observation date [ISO 8601]")
    parser.add_argument("--email", type=str,  default=JSOC_EMAIL, help="JSOC registered email")
    parser.add_argument("--out",   type=str,  default=OUTPUT_FILE, help="Output JSON path")
    args = parser.parse_args()

    globals()["JSOC_EMAIL"] = args.email

    print(f"\n{BOLD}{CYAN}{'═'*52}")
    print(f"  SDO/HMI ACTIVE REGION — FIELD LINE JSON CONVERTER")
    print(f"{'═'*52}{RESET}")
    print(f"  Target:   NOAA AR{args.ar}  |  Date: {args.date}")
    print(f"  Grid:     {GRID_SIZE}×{GRID_SIZE} = {GRID_SIZE**2} pts  (ceiling: {MAX_POINTS})")
    print(f"  Output:   {args.out}")

    records   = None
    mode_used = "none"

    if args.mode in ("auto", "real"):
        records = mode_real(args.ar, args.date)
        if records: mode_used = "real"

    if records is None and args.mode in ("auto", "http"):
        records = mode_http()
        if records: mode_used = "http"

    if records is None and args.mode in ("auto", "synthetic", "http"):
        records   = mode_synthetic(ar_number=args.ar)
        mode_used = "synthetic"

    if records is None:
        print(f"\n{RED}All data source modes failed. Exiting.{RESET}")
        sys.exit(1)

    print_field_stats(records)

    meta = {
        "mode":          mode_used,
        "noaa_ar":       args.ar,
        "date":          args.date,
        "extrapolation": "potential_field_FFT_Alissandrakis1981_Nakagawa-Raadu1972",
        "z_mapping":     "tanh(Bz/B_p99) * scale",
        "grayscale":     "gamma_0.55(|Bz|/B_p99)",
        "grid":          f"{GRID_SIZE}x{GRID_SIZE}",
        "max_points":    MAX_POINTS,
    }

    out_path = Path(args.out)
    save_json(records, out_path, meta)

    gc.collect()
    print(f"{DIM}Resource usage may be reviewed with: python memory_guard.py{RESET}\n")


if __name__ == "__main__":
    main()
