# AstroSplat-Browser-native-volumetric-visualization-for-astrophysical-vector-fields — Prototype v2.0
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19377621.svg)](https://doi.org/10.5281/zenodo.19377621)

A browser-native, memory-constrained prototype for three-dimensional
visualisation of synthetic and observational planetary magnetic field
line datasets. All computation runs in-browser via Three.js (r128);
no server, no build step, no external dependencies beyond a CDN import.

---

## Files

| File | Purpose |
|---|---|
| `field_visualizer.html` | Main visualizer — open directly in a browser |
| `generate_magnetosphere.py` | Generates the Saturn polar flux tube dataset |
| `solar_to_splats.py` | Fetches or synthesises a solar active region magnetogram dataset |
| `generate_earth_bowshock.py` | Generates the Earth magnetopause / bow shock dataset |
| `memory_guard.py` | JSON compression and RAM monitoring utility |
| `magnetosphere_data.json` | Saturn polar flux tube — 4,000 points, 8 filaments |
| `magnetosphere_data_mini.json` | Compressed variant — 3-decimal float precision |
| `solar_active_region.json` | Synthetic solar active region — 3,600 points (60×60 grid) |
| `earth_bowshock.json` | Synthetic Earth magnetopause + bow shock — 1,456 points |

---

## Quick Start

1. Open `field_visualizer.html` in any modern browser (Chrome, Firefox, Edge).
2. Drag a dataset JSON file onto the drop zone, or click to browse.
3. Use the left-hand control panel to adjust rendering and overlays.

No installation required.

---

## Capabilities

### Rendering

Two render modes are available via the **Basic / Enhanced** toggle:

- **Basic** — `MeshPhongMaterial` with additive blending. Each field line
  point is represented as a low-polygon ellipsoid aligned to the local B-field
  direction. Colour is intensity-mapped grayscale with a blue-tinted floor.

- **Enhanced** — The same geometry, extended via `onBeforeCompile` shader
  injection to add: (a) view-dependent alpha, which reduces opacity when
  the ellipsoid is viewed end-on along its B-field axis, improving depth
  legibility; (b) a Fresnel-like rim glow term that reinforces the tubular
  geometry at grazing angles. These effects are perceptual aids only and do
  not represent physical emission or radiative transfer properties.

Rendering is rate-limited to 30 FPS and uses a hard ceiling of 5,000 points
to maintain stable performance on integrated GPUs.

### Reference Overlays

- **Reference grid** — a 20×20 scene-unit grid with labelled X/Y/Z axes
  and tick marks. Units correspond to the coordinate frame of the loaded
  dataset (see metadata). Toggle via the control panel.

- **Saturn body (schematic)** — a wireframe oblate spheroid representing
  Saturn's planetary body at correct oblateness (f = 0.09796, Archinal et al.
  2018), sized at a schematic 0.55 scene units equatorial radius. This is
  not to scale against field data unless the dataset uses correctly scaled
  physical coordinates. Includes a schematic ring system indicator.

### Spatial Domain Clipping

The **Domain Settings** panel provides a bounding box filter with four
presets:

| Preset | Domain [scene units] | Intended use |
|---|---|---|
| Saturn | ±15 | Saturn magnetospheric dataset |
| Solar AR | ±4 | Solar active region (60×60 grid) |
| Wide | ±30 | Large-scale or high-altitude datasets |
| None | (disabled) | No clipping — all points rendered |

Custom min/max values can be entered manually. Points outside the domain
are excluded from the renderer at each rebuild; the source data are not
modified. A warning is printed to the browser console if the domain
truncates a significant fraction of the loaded dataset.

### Multi-Dataset Merge

The **Dataset Merge** section controls how successive file loads interact:

- **Replace** (default) — discards all existing datasets and renders the
  newly loaded file alone.

- **Append** — retains all previously loaded datasets and renders them
  simultaneously. Each dataset is assigned a distinct source colour from a
  five-colour palette, blended with the point's intensity-based grayscale
  at 40% weight. Datasets retain their individual identities in the
  **Dataset Info** panel.

When Append mode is active with multiple datasets, a coordinate-proximity
deduplication pass is applied at render time. Points within `threshold`
scene units of an already-accepted point (first-encountered, by load order)
are discarded. The default threshold is 0.05 scene units. This prevents
GPU overdraw artifacts at dataset boundaries without modifying the source
records.

The deduplication step is a spatial filter only; it does not remove
physically distinct points that happen to be co-located in different
coordinate frames. Datasets recorded in different physical unit systems
must be normalised to a common frame before being merged for the
deduplication to have physical meaning.

### Metadata Display

When a dataset JSON includes a metadata header record
(`"source": "meta"`), its fields are displayed in the **Dataset Info**
panel. The following fields are recognised:

| Field | Description |
|---|---|
| `source_label` | Short identifier string for the dataset |
| `date` | Acquisition or generation date |
| `coord_units` | Physical units of the x, y, z coordinates |
| `field_units` | Physical units of Bx, By, Bz |
| `planetary_ref` | Planetary body and coordinate frame |
| `domain_scale` | Approximate spatial extent of the dataset |
| `measured` | Boolean — `true` if derived from real observations |
| `data_type` | Free-text description of the data provenance |

Fields not present in the JSON are displayed as "Unknown" or "Unspecified".
The metadata record is filtered from the render set at ingestion; it does
not contribute a visible point to the scene.

---

## Dataset JSON Schema

Each record in a compliant dataset JSON array must contain the following
fields. Additional fields are preserved but not used by the renderer.

```json
{
  "x":         0.0,
  "y":         0.0,
  "z":         0.0,
  "Bx":        0.0,
  "By":        0.0,
  "Bz":        0.0,
  "intensity": 0.0,
  "gray":      128
}
```

- `x`, `y`, `z`: position in the dataset's coordinate frame. Units are
  specified in the metadata record; the renderer treats them as scene units
  without conversion.
- `Bx`, `By`, `Bz`: magnetic field vector components in the same frame.
  The renderer normalises these internally to compute ellipsoid orientation.
- `intensity`: scalar field magnitude, used only for metadata statistics.
- `gray`: 8-bit grayscale luminance index [0, 255], used for per-point
  colouring. If absent, defaults to 128 (mid-grey).

A metadata header record with `"source": "meta"` may optionally be placed
as the first element of the array.

---

## Scale Handling

Coordinates are loaded and rendered without unit conversion. The
visualizer treats all values as dimensionless scene units. Physical
scale correspondence requires the user to confirm that the dataset's
coordinate frame matches the reference frame implied by any overlay:

- The Saturn reference body overlay is schematic only. It is sized at
  0.55 scene units equatorial radius regardless of the loaded dataset's
  physical scale. It conveys body shape and oblateness, not scale.

- The reference grid is labelled in scene units. When the Saturn polar
  flux tube dataset (`magnetosphere_data.json`) is loaded, one scene unit
  corresponds roughly to one Saturn radius (informal convention of the
  generator script); this is not enforced or verified by the visualizer.

- The Earth bow shock dataset (`earth_bowshock.json`) uses R_E (Earth radii)
  as the coordinate unit. The domain preset "Saturn" (±15) is appropriate
  for this dataset as well, given its ±16 R_E spatial extent.

---

## Included Datasets

### `magnetosphere_data.json` — Saturn Polar Flux Tube

- **Source**: fully synthetic (centred dipole Euler integration)
- **Generator**: `generate_magnetosphere.py`
- **Points**: 4,000 (8 filaments × 500 points each)
- **Physical model**: centred dipole field, Saturn axial alignment
- **Coordinate units**: informal Saturn radii (scene units)
- **Use**: primary prototype dataset; tests single-file load and
  Saturn-appropriate render settings

### `solar_active_region.json` — Synthetic Solar Active Region

- **Source**: synthetic (Gaussian bipole + potential field extrapolation,
  Nakagawa & Raadu 1972; Alissandrakis 1981)
- **Generator**: `solar_to_splats.py --mode synthetic`
- **Points**: 3,600 (60×60 photospheric grid)
- **Physical model**: potential field photosphere (FFT Green's function);
  Bz from magnetogram, Bx/By extrapolated
- **Coordinate units**: pixel-scale scene units
- **Use**: tests non-dipole field topology; Solar AR domain preset applies

### `earth_bowshock.json` — Earth Magnetopause / Bow Shock

- **Source**: synthetic (empirical parametrisation)
- **Generator**: `generate_earth_bowshock.py`
- **Points**: 1,456 (magnetospheric dipole lines + magnetosheath draped field)
- **Physical model**: Shue et al. (1997) magnetopause, Farris & Russell (1994)
  bow shock, centred tilted dipole (IGRF mean tilt 11.5°)
- **Coordinate units**: R_E (Earth radii)
- **Use**: second example dataset for append-mode testing; provides a
  structurally different field topology for multi-dataset merge

---

## Limitations

The following limitations apply to the current prototype and should be
considered when interpreting visualisation outputs:

1. **All datasets are synthetic.** No observational magnetometer or
   remote-sensing data are used. Field line geometry reflects the
   mathematical model, not any spacecraft measurement.

2. **Ellipsoidal splats are a proxy only.** Each field line sample point is
   represented as a low-polygon ellipsoid aligned to the local B-vector.
   This is a visualisation convention, not a physically meaningful geometric
   primitive for magnetic structures.

3. **Enhanced rendering is perceptual, not physical.** View-dependent alpha
   and the Fresnel rim glow are visual depth cues; they do not represent
   the emissive, absorptive, or scattering properties of the modelled plasma.

4. **No field-line connectivity is enforced.** Points are rendered
   independently. The apparent continuity of field line structures arises
   from the spatial density of sample points, not from explicit polyline
   geometry.

5. **Potential field extrapolation (solar dataset) is an approximation.**
   The Nakagawa–Raadu Green's function method recovers a current-free
   field above the photosphere. Real coronal fields include electric
   currents not captured by this model.

6. **Deduplication is spatial, not physical.** The coordinate-proximity
   deduplication in append mode operates on scene-unit coordinates. It has
   no physical basis when datasets use different coordinate systems or units.

7. **Memory ceiling is enforced but not validated against GPU limits.**
   The 5,000-point ceiling reduces overdraw load, but GPU memory limits
   vary by device. The in-browser heap monitor (`MEM` badge) tracks JS
   heap only; VRAM is not monitored.

8. **Saturn reference body scale is schematic.** The wireframe spheroid is
   sized for visual reference only. Its scene-unit size does not correspond
   to any specific physical radius in the loaded datasets.

---

## Future Development Recommendations

The following directions are suggested as next steps, ordered by scientific
legibility impact:

1. **Polyline field line rendering**: replace independent splats with
   connected `THREE.TubeGeometry` or `THREE.Line3` chains per traced field
   line. This would correctly convey field line connectivity and direction.

2. **Physical unit normalisation**: add a metadata-driven scale factor so
   datasets in different physical units (km, R_E, R_Saturn) can be
   co-displayed with correct relative scale.

3. **Measured data ingestion**: implement a loader for NASA CDAWeb or
   SPDF heliospheric ASCII files so that real spacecraft magnetometer
   records can be loaded alongside synthetic datasets.

4. **Inferred structure overlay**: the "Show inferred structures" toggle
   is currently a disabled placeholder. A future release should populate
   it with structures derived from field-line tracing through the loaded
   potential field, tagged as `data_type: 'inferred'` and rendered at
   reduced opacity to distinguish them from the source data.

5. **Coordinate frame selector**: add a UI dropdown to specify the
   physical coordinate frame of each loaded dataset (GSM, GSE, SIII, HGS),
   enabling the visualizer to apply rotation transforms before rendering.

## Citation

If you use AstroSplat in research, teaching, outreach, or derivative work,
please cite the repository release:

Figueiredo, V. M. F. (2026). *AstroSplat: Browser-native volumetric
visualization for astrophysical vector fields* (v0.1.0) [Computer software].
Zenodo. https://doi.org/10.5281/zenodo.19377621
