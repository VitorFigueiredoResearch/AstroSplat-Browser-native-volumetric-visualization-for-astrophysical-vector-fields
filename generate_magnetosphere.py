"""
generate_magnetosphere.py
═════════════════════════════════════════════════════════════════════════════
Saturn Polar Magnetic Flux Tube — Synthetic Field Line Dataset Generator
─────────────────────────────────────────────────────────────────────────────
Synthetic field line data are generated for a single magnetic flux tube
located at Saturn's geographic pole. The dataset is produced at reduced
fidelity, strictly capped at 5,000 sample points, to support operation
within a memory-constrained prototype visualizer environment.

Physical model:
  A centred magnetic dipole (axis aligned to Z) is evaluated along eight
  discrete field line filaments, each traced from the photospheric level
  upward through the corona. Filament trajectories are subject to a
  seeded stochastic displacement (Gaussian random walk, σ = NOISE_SCALE m)
  and a soft radial convergence term that causes field lines to tighten
  toward the flux tube axis at higher altitudes, consistent with
  conservation of magnetic flux in a divergence-free field.

Output: magnetosphere_data.json
  Array of records, each containing spatial coordinates (x, y, z),
  magnetic field vector components (Bx, By, Bz) in arbitrary units,
  field magnitude (intensity), and a grayscale luminance index (gray, 0–255).
"""

import json
import math
import random

# ── Point count constraint ───────────────────────────────────────────────────
MAX_POINTS = 5000

# ── Flux tube geometry parameters ────────────────────────────────────────────
TUBE_HEIGHT       = 10.0     # Total Z extent of the flux tube [scene units]
TUBE_RADIUS       = 0.6      # Maximum radial spread of the flux bundle [scene units]
N_STRANDS         = 8        # Number of discrete field line filaments
POINTS_PER_STRAND = 500      # Sample points per filament: 8 × 500 = 4,000 total
NOISE_SCALE       = 0.08     # Standard deviation of stochastic displacement per step

# ── Saturn dipole field parameters ───────────────────────────────────────────
# Saturn's axisymmetric dipole is tilted < 1° from the rotation axis
# (Russell & Dougherty 2010). The value below is exaggerated for visual
# discriminability in the prototype visualizer.
DIPOLE_TILT_DEG   = 3.0      # Dipole tilt angle relative to rotation axis [degrees]
FIELD_STRENGTH    = 20_000.0 # Dipole moment scaling constant [arbitrary units, nT-equivalent]


def dipole_field(x: float, y: float, z: float) -> tuple[float, float, float]:
    """
    Evaluate the magnetic field of a centred axial dipole at position (x, y, z).

    The magnetic dipole field is given by:
        B = (μ₀ m / 4π r⁵) · [3xz, 3yz, 2z² - x² - y²]
    where m is the dipole moment aligned to the Z axis.

    A regularisation term is applied to r⁵ to prevent singularity at r → 0.

    Parameters
    ----------
    x, y, z : float
        Cartesian coordinates [scene units].

    Returns
    -------
    (Bx, By, Bz) : tuple[float, float, float]
        Magnetic field vector components [arbitrary units].
    """
    r2 = x * x + y * y + z * z
    r  = math.sqrt(r2) if r2 > 1e-6 else 1e-3
    r5 = r ** 5 + 1e-6          # regularisation to prevent division by zero

    scale = FIELD_STRENGTH / r5
    bx = scale * 3.0 * x * z
    by = scale * 3.0 * y * z
    bz = scale * (2.0 * z * z - x * x - y * y)
    return bx, by, bz


def intensity_to_gray(intensity: float, i_min: float, i_max: float) -> int:
    """
    Map a scalar field magnitude to an 8-bit grayscale luminance index.

    Intensity values are linearly normalised to [0, 1] over the observed
    range [i_min, i_max], then gamma-compressed with exponent γ = 0.7
    to expand the perceptual mid-range before quantisation to [0, 255].

    Parameters
    ----------
    intensity : float
        Field magnitude value to be mapped.
    i_min, i_max : float
        Observed minimum and maximum of the intensity distribution.

    Returns
    -------
    int : Grayscale index in [0, 255].
    """
    span = i_max - i_min if i_max > i_min else 1.0
    normalised = (intensity - i_min) / span
    normalised = math.pow(max(0.0, min(1.0, normalised)), 0.7)
    return int(round(normalised * 255))


def generate_filament(filament_idx: int, n_filaments: int, n_points: int) -> list[dict]:
    """
    Trace a single field line filament through the flux tube volume.

    The filament is initialised at the base of the flux tube at an
    angular offset determined by its index within the bundle, then
    advanced upward in Z by a fixed step. At each step, position is
    updated by a soft radial restoring force toward the converging
    tube axis plus a stochastic displacement drawn from N(0, NOISE_SCALE²).
    The dipole field is evaluated analytically at each point.

    Parameters
    ----------
    filament_idx : int
        Zero-based index of this filament within the bundle.
    n_filaments : int
        Total number of filaments in the flux tube.
    n_points : int
        Number of sample points to generate along this filament.

    Returns
    -------
    list[dict] : Field line records with keys x, y, z, Bx, By, Bz, intensity.
    """
    rng = random.Random(filament_idx * 7919)    # deterministic seed per filament

    # Angular position within the flux bundle cross-section
    phi0 = (2.0 * math.pi * filament_idx) / n_filaments

    # Initial radial offset; distributes filaments across the tube cross-section
    r0 = TUBE_RADIUS * rng.uniform(0.15, 0.85)

    points = []
    # Filaments are initialised at the lower boundary of the flux tube
    cx, cy = r0 * math.cos(phi0), r0 * math.sin(phi0)
    cz = -TUBE_HEIGHT / 2.0

    step_z = TUBE_HEIGHT / n_points    # fixed vertical step per integration point

    # Dipole tilt introduces a small X-directed offset proportional to height
    tilt_rad = math.radians(DIPOLE_TILT_DEG)
    tilt_x   = math.sin(tilt_rad)

    for i in range(n_points):
        frac  = i / n_points           # normalised height parameter in [0, 1)

        # Radial convergence: filaments contract toward the axis at higher altitudes,
        # consistent with the divergence-free constraint ∇·B = 0
        r_now    = r0 * (1.0 - 0.4 * frac)
        target_x = r_now * math.cos(phi0) + tilt_x * cz * 0.1
        target_y = r_now * math.sin(phi0)

        # Euler integration step: restoring force + stochastic displacement
        cx += (target_x - cx) * 0.08 + rng.gauss(0, NOISE_SCALE)
        cy += (target_y - cy) * 0.08 + rng.gauss(0, NOISE_SCALE)
        cz += step_z

        # Dipole field evaluated at the updated position
        bx, by, bz = dipole_field(cx, cy, cz)
        mag = math.sqrt(bx * bx + by * by + bz * bz)

        points.append({
            "x":  round(cx, 5),
            "y":  round(cy, 5),
            "z":  round(cz, 5),
            "Bx": round(bx, 4),
            "By": round(by, 4),
            "Bz": round(bz, 4),
            "intensity": round(mag, 4),
        })

    return points


def main():
    print("Generating Saturn polar flux tube dataset...")
    all_points: list[dict] = []

    for s in range(N_STRANDS):
        filament = generate_filament(s, N_STRANDS, POINTS_PER_STRAND)
        all_points.extend(filament)
        print(f"  Filament {s+1}/{N_STRANDS}: {len(filament)} points  "
              f"(cumulative: {len(all_points)})")

    # ── Point count constraint (hardware safety) ─────────────────────────────
    if len(all_points) > MAX_POINTS:
        all_points = all_points[:MAX_POINTS]
        print(f"Note: dataset truncated to {MAX_POINTS} points per hardware constraint.")

    print(f"Dataset complete: {len(all_points)} sample points.")

    # ── Grayscale encoding ────────────────────────────────────────────────────
    # Intensity values are normalised over the full dataset range and mapped
    # to an 8-bit luminance index with gamma compression (γ = 0.7).
    intensities = [p["intensity"] for p in all_points]
    i_min, i_max = min(intensities), max(intensities)

    for p in all_points:
        p["gray"] = intensity_to_gray(p["intensity"], i_min, i_max)

    # ── JSON output ───────────────────────────────────────────────────────────
    out_path = "magnetosphere_data.json"
    with open(out_path, "w") as f:
        json.dump(all_points, f, separators=(",", ":"))   # compact encoding

    size_kb = len(json.dumps(all_points, separators=(",", ":"))) / 1024
    print(f"Written: {out_path}  ({size_kb:.1f} KB)")
    print("Dataset ready. Load this file in field_visualizer.html to render the flux tube.")


if __name__ == "__main__":
    main()
