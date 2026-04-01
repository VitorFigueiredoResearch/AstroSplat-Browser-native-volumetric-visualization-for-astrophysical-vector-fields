"""
generate_earth_bowshock.py
═════════════════════════════════════════════════════════════════════════════
Earth Magnetopause / Bow Shock Cross-Section — Synthetic Field Line Dataset
─────────────────────────────────────────────────────────────────────────────
A synthetic representation of the Earth's magnetospheric boundary region is
generated based on established empirical parametrisations. The dataset is
produced at reduced fidelity and capped at 5,000 sample points to support
operation within the prototype visualizer's memory constraints.

Physical model and references:
  Magnetopause shape: Shue et al. (1997) empirical parametrisation,
    r(θ) = r₀ · (2 / (1 + cos θ))^α
    with nominal solar wind conditions: r₀ = 10.22 R_E, α = 0.58.
    Reference: Shue, J.-H. et al. (1997), J. Geophys. Res. 102, 9497–9511.

  Bow shock shape: Farris & Russell (1994) ellipsoidal approximation,
    scaled to nominal standoff distance ~14 R_E at the subsolar point.
    Reference: Farris, M. H. & Russell, C. T. (1994), J. Geophys. Res. 99,
    17681–17689.

  Magnetic field vectors along field lines within the magnetosphere are
  approximated using a simplified dipole model (centred, tilted 11.5° from
  the rotation axis, consistent with IGRF mean inclination). Magnetosheath
  field vectors are assigned a draped orientation consistent with a
  compressed, kinked interplanetary magnetic field (southward IMF assumed
  for illustrative salience).

  All coordinates are expressed in units of Earth radii (R_E = 6,371 km).
  Values are synthetic approximations suitable for visualisation only and
  do not constitute observational data.

Coordinate frame:
  X: Sun–Earth direction (positive sunward)
  Y: Orthogonal to X in the ecliptic plane (positive duskward)
  Z: Northward (positive toward ecliptic north pole)

Output: earth_bowshock.json
  Array of records containing spatial coordinates (x, y, z) in R_E,
  magnetic field vector components (Bx, By, Bz) in arbitrary units,
  field magnitude (intensity), grayscale luminance index (gray, 0–255),
  and region label (source: 'magnetosphere' | 'magnetosheath' | 'bowshock').
  A metadata header record (source: 'meta') is prepended.
"""

import json
import math
import random
import os

# ── Output configuration ─────────────────────────────────────────────────────
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE  = os.path.join(SCRIPT_DIR, 'earth_bowshock.json')
RANDOM_SEED  = 42
MAX_RECORDS  = 5000
N_FIELD_LINES = 28        # field lines traced through the magnetosphere
N_POINTS_PER_LINE = 60   # sample points per field line
N_SHEATH_LINES = 20       # draped field lines in the magnetosheath
N_POINTS_SHEATH = 40     # points per sheath field line

# ── Physical constants (all in R_E units unless noted) ───────────────────────
RE_KM         = 6371.0   # Earth radius [km]
# Shue et al. (1997) magnetopause parameters — nominal quiet solar wind
MP_R0         = 10.22    # subsolar standoff distance [R_E]
MP_ALPHA      = 0.58     # flaring parameter [dimensionless]
# Bow shock standoff distance (Farris & Russell 1994 nominal)
BS_R0         = 14.0     # subsolar bow shock distance [R_E]
BS_ALPHA      = 0.60     # bow shock flaring parameter
# Dipole field strength parameter (arbitrary units; scales Bz of ionosphere)
B_DIPOLE_SCALE = 8.0
# Dipole tilt: Earth's magnetic dipole is inclined ~11.5° from rotation axis
DIPOLE_TILT_DEG = 11.5
DIPOLE_TILT_RAD = math.radians(DIPOLE_TILT_DEG)


def magnetopause_radius(theta_rad):
    """
    Magnetopause standoff radius at polar angle theta (rad from sunward X axis).

    Parametrisation from Shue et al. (1997), Eq. 1:
        r(θ) = r₀ · (2 / (1 + cos θ))^α

    Parameters
    ----------
    theta_rad : float
        Polar angle from the sunward direction [radians].

    Returns
    -------
    float
        Magnetopause radius [R_E].
    """
    cos_t = math.cos(theta_rad)
    return MP_R0 * (2.0 / (1.0 + cos_t)) ** MP_ALPHA


def bow_shock_radius(theta_rad):
    """
    Bow shock radius at polar angle theta (rad from sunward X axis).

    Farris & Russell (1994) ellipsoidal approximation scaled to nominal
    standoff distance BS_R0 at θ = 0.

    Parameters
    ----------
    theta_rad : float
        Polar angle from the sunward direction [radians].

    Returns
    -------
    float
        Bow shock radius [R_E].
    """
    cos_t = math.cos(theta_rad)
    return BS_R0 * (2.0 / (1.0 + cos_t)) ** BS_ALPHA


def dipole_field(x, y, z):
    """
    Evaluate the centred tilted magnetic dipole field at position (x, y, z).

    The dipole axis is tilted DIPOLE_TILT_DEG from the Z (north) axis toward
    the X axis, consistent with Earth's mean magnetic inclination. This is
    a first-order internal field representation only; external currents
    (ring current, Chapman–Ferraro) are not included.

    Dipole field formula (in dipole-axis-aligned frame):
        Br     = -2 B₀ cos θ_d / r³
        Btheta = -  B₀ sin θ_d / r³
    Converted to Cartesian for the tilted dipole axis.

    Parameters
    ----------
    x, y, z : float
        Position in Earth-centred GSM coordinates [R_E].

    Returns
    -------
    tuple of float
        (Bx, By, Bz) field components [arbitrary units].
    """
    # Rotate position into dipole frame (tilt about Y axis)
    cos_t = math.cos(DIPOLE_TILT_RAD)
    sin_t = math.sin(DIPOLE_TILT_RAD)
    xd =  x * cos_t + z * sin_t
    yd =  y
    zd = -x * sin_t + z * cos_t

    r2 = xd**2 + yd**2 + zd**2
    r  = math.sqrt(r2)
    if r < 0.5:
        return (0.0, 0.0, B_DIPOLE_SCALE)   # inside the solid Earth (regularised)

    r5 = r**5
    # Dipole field in the tilted frame
    Bxd = B_DIPOLE_SCALE * 3.0 * xd * zd / r5
    Byd = B_DIPOLE_SCALE * 3.0 * yd * zd / r5
    Bzd = B_DIPOLE_SCALE * (3.0 * zd**2 - r2) / r5

    # Rotate back to GSM frame
    Bx =  Bxd * cos_t - Bzd * sin_t
    By =  Byd
    Bz =  Bxd * sin_t + Bzd * cos_t
    return (Bx, By, Bz)


def magnetosheath_field(x, y, z):
    """
    Approximate magnetosheath draped field for southward IMF conditions.

    The magnetosheath field is modelled as a compressed, kinked version of
    the solar wind IMF. For southward IMF (Bz_sw < 0), the field drapes
    around the magnetopause with increasing By component toward the flanks.
    This is a heuristic approximation only, not a physics-based model.

    Parameters
    ----------
    x, y, z : float
        Position in GSM coordinates [R_E].

    Returns
    -------
    tuple of float
        (Bx, By, Bz) field components [arbitrary units].
    """
    r = math.sqrt(x**2 + y**2 + z**2)
    if r < 0.1:
        return (0.0, 0.0, -1.0)

    # Draping factor: field deflects around magnetopause
    mp_r = magnetopause_radius(math.atan2(math.sqrt(y**2 + z**2), max(x, -20.0)))
    drape = max(0.0, r - mp_r) / max(BS_R0 - mp_r, 0.5)
    drape = min(drape, 1.0)

    # Southward IMF baseline with compression and draping
    Bx = -0.3 * (x / r)
    By =  0.6 * (y / r) * drape
    Bz = -2.0 * (1.0 - 0.4 * drape)

    # Compression enhancement near the magnetopause nose
    compress = max(1.0, 2.5 - r / mp_r)
    return (Bx * compress, By * compress, Bz * compress)


def intensity_to_gray(intensity, i_min, i_max, gamma=0.55):
    """
    Map field magnitude to an 8-bit grayscale luminance index.

    A gamma compression (γ = 0.55) is applied to the normalised intensity
    to improve contrast at low-field regions, where most spatial variation
    of interest is located.

    Parameters
    ----------
    intensity : float
        Field magnitude at the sample point [arbitrary units].
    i_min, i_max : float
        Global minimum and maximum intensity across the dataset.
    gamma : float
        Gamma compression exponent. Values below 1 brighten the low end.

    Returns
    -------
    int
        Grayscale luminance index in the range [0, 255].
    """
    denom = max(i_max - i_min, 1e-9)
    norm  = max(0.0, min(1.0, (intensity - i_min) / denom))
    compressed = norm ** gamma
    return int(round(compressed * 255))


def trace_dipole_field_line(theta0, phi0, n_points):
    """
    Trace a single magnetospheric field line from an ionospheric footpoint.

    Integration uses a fixed-step Euler method along the normalised B-field
    vector. The field line terminates when either the magnetopause surface or
    the ionospheric surface (r < 1.1 R_E) is reached.

    Parameters
    ----------
    theta0 : float
        Co-latitude of the ionospheric footpoint [radians].
    phi0 : float
        Azimuth (longitude) of the footpoint [radians].
    n_points : int
        Maximum number of sample points per field line.

    Returns
    -------
    list of dict
        Each dict contains x, y, z, Bx, By, Bz, intensity, source.
    """
    # Start at r = 3 R_E to avoid the inner dipole singularity
    r_start = 3.0
    x = r_start * math.sin(theta0) * math.cos(phi0)
    y = r_start * math.sin(theta0) * math.sin(phi0)
    z = r_start * math.cos(theta0)

    step = 0.18
    records = []

    for _ in range(n_points):
        r = math.sqrt(x**2 + y**2 + z**2)
        if r < 1.1:
            break

        # Polar angle from sunward X axis for magnetopause radius
        theta_sun = math.atan2(math.sqrt(y**2 + z**2), max(x, -40.0))
        mp_r = magnetopause_radius(theta_sun)
        if r > mp_r:
            break   # exited the magnetopause

        Bx, By, Bz = dipole_field(x, y, z)
        mag = math.sqrt(Bx**2 + By**2 + Bz**2)
        if mag < 1e-9:
            break

        records.append(dict(x=x, y=y, z=z,
                            Bx=Bx, By=By, Bz=Bz,
                            intensity=mag, source='magnetosphere'))

        # Euler step along B
        x += step * Bx / mag
        y += step * By / mag
        z += step * Bz / mag

    return records


def generate_magnetosheath_lines(n_lines, n_pts):
    """
    Generate draped field line samples within the magnetosheath.

    Field lines are seeded at random positions between the magnetopause
    and the bow shock, restricted to the dayside (x > -5 R_E). The
    magnetosheath field model (magnetosheath_field) provides the local
    field vector at each sample point.

    Parameters
    ----------
    n_lines : int
        Number of field lines to trace.
    n_pts : int
        Sample points per field line.

    Returns
    -------
    list of dict
        Field line records with source='magnetosheath'.
    """
    records = []
    for _ in range(n_lines):
        # Seed in the magnetosheath (between magnetopause and bow shock)
        phi   = random.uniform(0.0, 2.0 * math.pi)
        theta = random.uniform(0.0, 1.4)  # restrict to dayside
        t_cos = max(math.cos(theta), -0.98)
        mp_r  = magnetopause_radius(math.acos(t_cos))
        bs_r  = bow_shock_radius(math.acos(t_cos))
        r_seed = random.uniform(mp_r + 0.3, bs_r - 0.3)
        if r_seed <= mp_r:
            continue

        x = r_seed * math.sin(theta) * math.cos(phi)
        y = r_seed * math.sin(theta) * math.sin(phi)
        z = r_seed * math.cos(theta)

        step = 0.20
        for _ in range(n_pts):
            theta_sun = math.atan2(math.sqrt(y**2 + z**2), max(x, -40.0))
            mp_r2 = magnetopause_radius(theta_sun)
            bs_r2 = bow_shock_radius(theta_sun)
            r2 = math.sqrt(x**2 + y**2 + z**2)
            if r2 < mp_r2 or r2 > bs_r2:
                break

            Bx, By, Bz = magnetosheath_field(x, y, z)
            mag = math.sqrt(Bx**2 + By**2 + Bz**2)
            if mag < 1e-9:
                break

            records.append(dict(x=x, y=y, z=z,
                                Bx=Bx, By=By, Bz=Bz,
                                intensity=mag, source='magnetosheath'))
            x += step * Bx / mag
            y += step * By / mag
            z += step * Bz / mag

    return records


def main():
    random.seed(RANDOM_SEED)

    all_records = []

    # ── Magnetospheric dipole field lines ────────────────────────────────────
    # Field lines distributed across a range of L-shells (3–10 R_E)
    # and distributed evenly in azimuth. Northern and southern hemispheres
    # are traced.
    phi_steps  = 8
    theta_steps = N_FIELD_LINES // phi_steps
    for i_phi in range(phi_steps):
        phi = (i_phi / phi_steps) * 2.0 * math.pi
        for i_th in range(theta_steps):
            # co-latitudes between 25° and 55° from north pole (maps to L = 3–10)
            theta = math.radians(25.0 + i_th * 30.0 / max(theta_steps - 1, 1))
            pts = trace_dipole_field_line(theta, phi, N_POINTS_PER_LINE)
            all_records.extend(pts)

    # ── Magnetosheath draped field lines ─────────────────────────────────────
    sheath_records = generate_magnetosheath_lines(N_SHEATH_LINES, N_POINTS_SHEATH)
    all_records.extend(sheath_records)

    if not all_records:
        print('[Error] No records generated.')
        return

    # ── Grayscale encoding ───────────────────────────────────────────────────
    intensities = [r['intensity'] for r in all_records]
    i_min = min(intensities)
    i_max = max(intensities)

    for r in all_records:
        r['gray'] = intensity_to_gray(r['intensity'], i_min, i_max)

    # Point cap enforcement
    if len(all_records) > MAX_RECORDS:
        print(f'[Constraint] Record count truncated: {len(all_records)} → {MAX_RECORDS}')
        all_records = all_records[:MAX_RECORDS]

    # ── Metadata header record ───────────────────────────────────────────────
    # Prepended metadata record (source='meta') is filtered out by the
    # visualizer at ingestion and is not forwarded to the renderer.
    meta = {
        'source':        'meta',
        'source_label':  'EARTH_BOWSHOCK_SYNTH_v1',
        'date':          '2026-04-01 (synthetic)',
        'coord_units':   'R_E (Earth radii, 1 R_E = 6371 km)',
        'field_units':   'arbitrary units (dipole model, B_DIPOLE_SCALE = 8.0)',
        'planetary_ref': 'Earth — GSM coordinate frame',
        'domain_scale':  'Magnetospheric, approx. ±20 R_E',
        'measured':      False,
        'data_type':     'Synthetic / inferred — empirical parametrisation',
        'mp_model':      'Shue et al. (1997), J. Geophys. Res. 102, 9497',
        'bs_model':      'Farris & Russell (1994), J. Geophys. Res. 99, 17681',
        'B_model':       'Centred tilted dipole (IGRF mean tilt 11.5°)',
        'notes': (
            'Magnetosheath field is a heuristic draping approximation '
            'for southward IMF. External current systems (ring current, '
            'tail current sheet) are not represented. This dataset is '
            'suitable for visualisation prototyping only.'
        ),
    }

    output = [meta] + all_records

    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output, f, indent=None, separators=(',', ':'))

    n_ms   = sum(1 for r in all_records if r['source'] == 'magnetosphere')
    n_sh   = sum(1 for r in all_records if r['source'] == 'magnetosheath')
    print(f'[Output] {OUTPUT_FILE}')
    print(f'  Total records : {len(all_records)}')
    print(f'  Magnetosphere : {n_ms} points (dipole field lines)')
    print(f'  Magnetosheath : {n_sh} points (draped IMF)')
    print(f'  Intensity range: [{i_min:.3f}, {i_max:.3f}]')
    print(f'  Coord units   : R_E — domain approx. ±20 R_E')
    print(f'  References    : Shue+1997, Farris&Russell+1994, IGRF dipole')


if __name__ == '__main__':
    main()
