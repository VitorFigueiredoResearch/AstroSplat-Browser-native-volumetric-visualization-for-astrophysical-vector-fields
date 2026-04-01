"""
Saturn Polar Magnetic Flux Tube — "Magnetic Hair" Data Generator
================================================================
Generates a synthetic dataset of a single magnetic flux tube at Saturn's pole.
Strictly capped at 5,000 points for hardware safety.

Output: magnetic_hair.json
"""

import json
import math
import random

# ── Safety cap ──────────────────────────────────────────────────────────────
MAX_POINTS = 5000

# ── Flux tube parameters ─────────────────────────────────────────────────────
TUBE_HEIGHT      = 10.0     # Z extent of the tube (aligned to Saturn's pole axis)
TUBE_RADIUS      = 0.6      # Radial spread of the flux tube
N_STRANDS        = 8        # Sub-filaments within the tube
POINTS_PER_STRAND = 500     # 8 × 500 = 4,000 points — safely under 5 k
NOISE_SCALE      = 0.08     # Slight random walk per step (keeps it looking organic)

# ── Saturn dipole field approximation ───────────────────────────────────────
# At the pole (θ→0), B is nearly axial. We tilt slightly and add radial twist.
DIPOLE_TILT_DEG  = 3.0      # Saturn's dipole is almost perfectly aligned (~< 1°, exaggerated for visual)
FIELD_STRENGTH   = 20_000.0 # Nanotesla-like scale (arbitrary units)


def dipole_field(x: float, y: float, z: float) -> tuple[float, float, float]:
    """
    Simplified analytic dipole centred at origin, axis along Z.
    B = (B0/r^5) * [3xz, 3yz, 2z^2 - x^2 - y^2]  (dipole formula)
    Clamped to avoid singularities at r→0.
    """
    r2 = x * x + y * y + z * z
    r  = math.sqrt(r2) if r2 > 1e-6 else 1e-3
    r5 = r ** 5 + 1e-6          # avoid ÷0

    scale = FIELD_STRENGTH / r5
    bx = scale * 3.0 * x * z
    by = scale * 3.0 * y * z
    bz = scale * (2.0 * z * z - x * x - y * y)
    return bx, by, bz


def intensity_to_gray(intensity: float, i_min: float, i_max: float) -> int:
    """Map intensity to [0, 255] grayscale (Step 2 requirement)."""
    span = i_max - i_min if i_max > i_min else 1.0
    normalized = (intensity - i_min) / span
    # Gamma-correct slightly so mid-range isn't too washed out
    normalized = math.pow(max(0.0, min(1.0, normalized)), 0.7)
    return int(round(normalized * 255))


def generate_strand(strand_idx: int, n_strands: int, n_points: int) -> list[dict]:
    """Generate one filament of the flux tube."""
    rng = random.Random(strand_idx * 7919)          # reproducible per strand

    # Angular offset for this strand
    phi0 = (2.0 * math.pi * strand_idx) / n_strands

    # Slight initial radial offset so filaments don't all start at centre
    r0 = TUBE_RADIUS * rng.uniform(0.15, 0.85)

    points = []
    # Start at the bottom of the tube (near Saturn's surface)
    cx, cy = r0 * math.cos(phi0), r0 * math.sin(phi0)
    cz = -TUBE_HEIGHT / 2.0

    step_z = TUBE_HEIGHT / n_points          # climb along Z each step

    # Tilt vector (slight dipole misalignment)
    tilt_rad = math.radians(DIPOLE_TILT_DEG)
    tilt_x   = math.sin(tilt_rad)

    for i in range(n_points):
        # Random walk in XY, drift toward axis at upper heights
        frac  = i / n_points                 # 0→1 from bottom to top
        # Filaments converge slightly toward the axis as they rise
        r_now = r0 * (1.0 - 0.4 * frac)
        target_x = r_now * math.cos(phi0) + tilt_x * cz * 0.1
        target_y = r_now * math.sin(phi0)

        # Soft pull toward target + small noise
        cx += (target_x - cx) * 0.08 + rng.gauss(0, NOISE_SCALE)
        cy += (target_y - cy) * 0.08 + rng.gauss(0, NOISE_SCALE)
        cz += step_z

        # Evaluate dipole field at this position
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
    print("🌀  Generating Saturn magnetic flux tube dataset...")
    all_points: list[dict] = []

    for s in range(N_STRANDS):
        strand = generate_strand(s, N_STRANDS, POINTS_PER_STRAND)
        all_points.extend(strand)
        print(f"  Strand {s+1}/{N_STRANDS} → {len(strand)} pts  (total so far: {len(all_points)})")

    # ── Hard safety cap ──────────────────────────────────────────────────────
    if len(all_points) > MAX_POINTS:
        all_points = all_points[:MAX_POINTS]
        print(f"⚠️  Capped at {MAX_POINTS} points for hardware safety.")

    print(f"✅  Total points: {len(all_points)}")

    # ── Normalise intensities and add grayscale (Step 2) ────────────────────
    intensities = [p["intensity"] for p in all_points]
    i_min, i_max = min(intensities), max(intensities)

    for p in all_points:
        p["gray"] = intensity_to_gray(p["intensity"], i_min, i_max)

    # ── Write JSON ───────────────────────────────────────────────────────────
    out_path = "magnetic_hair.json"
    with open(out_path, "w") as f:
        json.dump(all_points, f, separators=(",", ":"))   # compact, no whitespace

    size_kb = len(json.dumps(all_points, separators=(",", ":"))) / 1024
    print(f"💾  Written → {out_path}  ({size_kb:.1f} KB)")
    print("Done. Open the HTML viewer and load this file to see your magnetic hair! 🚀")


if __name__ == "__main__":
    main()
