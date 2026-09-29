# Caliber Validation v2: mesh cross-sectional equivalent diameter

## Scope and motivation

This is geometric/synthetic validation, **not clinical validation**. V1 found grid-dependent EDT error, particularly for 1 mm cylinders, mixed bias signs, and inconsistent benefits from endpoint trimming. V2 compares a different geometric estimand on the same masks and recovered sample positions, then adds modest grid-phase and orientation sweeps. The new method is experimental; the production Case 001 estimator, frontend/API measurements, source NIfTI, and archived v1 files have not been changed.

## Reproduce

From the repository root:

```bash
source .venv/bin/activate
python3 scripts/run_caliber_validation_v2.py
python3 -m pytest -q
```

The unchanged v1 command still works. To verify it without overwriting the archived v1 results:

```bash
python3 scripts/run_caliber_validation.py --output validation/outputs/caliber_v2/v1_reproduction
```

That isolated reproduction produced `results.csv` and `profiles.csv` exactly equal to the archived v1 CSVs. Its dependency/source metadata can differ because the source tree now also contains experimental v2 modules. All original v1 files remain byte-identical.

V2 output defaults to `validation/outputs/caliber_v2/`; the runner rejects the v1 output directory, its descendants, and destinations outside the validation-output tree. Re-runs overwrite only the selected validation output directory's artifacts. There is no random seed or variable timestamp because the matrix is deterministic. `summary.json` records dependency versions, source hashes, parameters and counting conventions.

## Methods

### Unchanged EDT baseline

The existing Lee skeletonization, vascular graph, ordered main-path extraction, and `build_diameter_profile` are reused unchanged. Raw EDT diameter is twice distance to the nearest background voxel center, with physical voxel spacing. V2's core EDT metrics reproduce all 52 archived v1 runs. No Gaussian diameter smoothing is applied to either estimator in the head-to-head results.

### Experimental cross-sectional method

`neurovasc/geometry/cross_sectional_caliber.py` accepts a `Volume` and an ordered `MainCenterlinePath`. It builds the existing binary 0.5 marching-cubes surface using `mask_to_mesh`; that function applies the NIfTI affine to produce world coordinates in mm. It constructs the mesh once per phantom, not once per sample.

For every sample:

1. Estimate a normalized local tangent in physical coordinates.
2. Choose the Cartesian axis least aligned with the tangent; form orthonormal in-plane vectors by cross products. This avoids near-parallel basis construction.
3. Intersect the surface with the plane through the sample, normal to the tangent, using [PyVista's plane slicing](https://docs.pyvista.org/api/core/_autosummary/pyvista.dataobjectfilters.slice).
4. Merge intersection vertices with an **absolute 1e-6 mm tolerance**. PyVista documents absolute versus relative tolerances in [mesh cleaning](https://docs.pyvista.org/api/core/_autosummary/pyvista.polydatafilters.clean).
5. Reconstruct the contour from unique undirected line segments. Require every vertex to have degree two and exactly one closed loop with at least three vertices.
6. Project that loop into the plane basis and calculate polygon area by the shoelace formula, in mm². Require a finite positive area and that the sample origin lies inside the contour.
7. Compute `D_eq = 2 * sqrt(area_mm2 / pi)`.

This is not EDT under another name, nor is it an analytic-truth calculation. Area depends only on the voxel-derived mesh and recovered path. It measures the 0.5 isosurface of a binary mask, which still has voxelization and marching-cubes bias. It does not recover a subvoxel boundary that was never present in the segmentation.

### Tangent definition

Tangents use a local least-squares fit of physical xyz coordinates against cumulative physical path arc length. The explicit window is **6 mm total**, centered on each sample (±3 mm along the recovered path). Windows are clipped at endpoints, making those fits one-sided. If the window contains fewer than two samples, the two nearest arc-length samples are used. The fitted vector is normalized; zero-length, nonfinite, or directionless fits are invalid.

The full path is used for tangent fitting **before trimming**. V2 then asks the unchanged production profile function which sample coordinates survive each trim setting and selects those exact sections. Thus the two trim variants do not silently use different tangents for the same point. The existing Gaussian caliber smoothing is not used for tangent estimation.

The fixed 6 mm window reduces voxel-step direction noise but introduces a scale choice and can average over sharp bends. Fits use equal sample weights, not continuous arc-length quadrature. Tests check straight direction, reversal, endpoint behavior, analytic-arc angle error, normalization, and affine rotation equivariance.

### Failure handling

Invalid sections retain their sample slot and reason. The Python estimator returns NaN area/diameter and `valid=False`; CSV uses blank numeric values and explicit status. Failure reasons include invalid tangents/origins, empty intersections, open/branched contours, multiple contours, degenerate area, outside-contour origins, and mesh-construction failure. Multiple loops are conservatively rejected, never summed. No EDT fallback or interpolation hides missing cross-sectional estimates.

A shared centerline failure is different: neither method has any sample positions to evaluate. Both run rows are retained with `status=upstream_failed`, zero attempted samples, a null valid fraction and a diagnostic reason. The matched pair is unscorable. These failures are counted separately from section invalidity.

## Experiment matrix and units of counting

All phantoms are 30 mm long. V1 analytic voxel-center inclusion, padding, coordinate-based ground truth, flat caps and binary mask construction remain unchanged in the core cohort.

| Cohort | Added geometries | Definition |
|---|---:|---|
| Core v1 | 26 | Straight 1–5 mm; oblique/curved 2,3,4 mm; tapers 4→2 and 3→1.5 mm; both spacings |
| Grid phase | 20 | Straight 1–5 mm, both spacings, additional (+0.25,+0.25,+0.25) and (+0.5,+0.5,+0.5) voxel grid offsets |
| Orientation | 12 | Oblique 2,3,4 mm, both spacings, additional normalized directions (1,2,1) and (2,1,3) |

The centered phase and normalized (1,1,1) direction are already present in the core cohort. The grid moves relative to the fixed analytic vessel; the affine contains that shift, and truth remains evaluated at the physical sample coordinates.

Spacings are `(0.5,0.5,0.5)` and `(0.46875,0.46875,0.8)` mm. Curved phantoms retain the v1 24 mm bend radius and x-z arc. Taper truth is evaluated from physical position along the known axis, not normalized recovered path distance.

There are **58 geometries × 2 trim settings = 116 paired experiments**, or **232 estimator run rows**. Of these, **92 pairs are scorable** and **24 pairs fail upstream** (12 geometries, each with two trim settings). Paired trims are not independent replicates. No partial-volume experiment is included in v2.

## Metrics and fair comparison

`results.csv` retains all requested raw diameter summaries, MAE/RMSE/max absolute error, bias, absolute mean bias, absolute/signed relative percentage errors, valid/invalid counts, valid fraction, and taper correlation when defined. Phase, orientation, spacing, cohort and trim identify each run.

Per-estimator metrics use that estimator's valid samples. `matched_comparison.csv` instead calculates MAE/RMSE/bias differences on the **intersection of valid sample positions**, with common coverage explicitly recorded. Negative `cross_sectional_mae_minus_edt_mae` favors cross-sectional. A MAE difference of at most **0.01 mm** is a tie. No winner is assigned to an unscorable pair. In this run every attempted section was valid, so own-support and common-support metrics use identical samples.

Aggregate MAE is the equal-run mean of run MAEs. Aggregate RMSE is `sqrt(mean(run MSE))`; it is not the mean of run RMSEs. Both trim settings receive equal weight. This weighting is applied only to scorable runs, with failed counts always shown. Mean differences across geometries are descriptive; deterministic phantoms do not supply clinical confidence intervals.

## Main results

| Evaluation set | EDT MAE (mm) | Cross-sectional MAE (mm) | EDT RMSE (mm) | Cross-sectional RMSE (mm) |
|---|---:|---:|---:|---:|
| Core v1, all 52 pairs | 0.201283 | 0.099387 | 0.252097 | 0.121852 |
| Expanded matrix, 92 scorable pairs | 0.229955 | 0.080268 | 0.303408 | 0.102973 |

Expanded-matrix signed bias is **-0.039117 mm EDT** and **-0.021461 mm cross-sectional**. Core EDT retains its v1 positive bias **+0.093299 mm**; core cross-sectional bias is **-0.010751 mm**. Changing the experiment mix changes the overall sign: neither method has a universal bias correction.

Cross-sectional wins **78**, EDT wins **8**, and **6** comparisons tie. **24** remain unscorable. The mean matched MAE delta is **-0.149688 mm**; aggregate RMSE delta is **-0.200434 mm**. Oblique geometry has the largest mean improvement, curved geometry the smallest (both average deltas favor cross-sectional). These labels describe mean differences, not guaranteed performance of every vessel.

### 1. Does it reduce the 1 mm vessel error?

For centered, trimmed straight cylinders:

| Spacing | EDT mean (mm) / error (%) | Cross-sectional mean (mm) / error (%) |
|---|---:|---:|
| Isotropic | 1.414214 / 41.4214% | 1.196827 / 19.6827% |
| Case-like anisotropic | 1.325825 / 32.5825% | 1.122025 / 12.2025% |

It reduces the large centered-grid error but does not eliminate it. The largest cross-sectional run percentage error remains **19.6827%** on the centered 1 mm isotropic cylinder (both trims have essentially identical error). At quarter-voxel isotropic phase, EDT happens to measure **1.000000 mm** while cross-sectional gives **0.892062 mm** (10.7938% error). At half-voxel phase the skeleton is empty. These contrasts show why a single aligned phantom cannot establish small-vessel accuracy.

### 2. Does it reduce orientation dependence?

Mean MAE across 2,3,4 mm oblique diameters at trim=2 mm:

| Spacing | Direction | EDT MAE (mm) | Cross-sectional MAE (mm) |
|---|---|---:|---:|
| Isotropic | (1,1,1) | 0.225266 | 0.118159 |
| Isotropic | (1,2,1) | 0.346479 | 0.072480 |
| Isotropic | (2,1,3) | 0.277732 | 0.052395 |
| Anisotropic | (1,1,1) | 0.264187 | 0.063418 |
| Anisotropic | (1,2,1) | 0.350519 | 0.059072 |
| Anisotropic | (2,1,3) | 0.229746 | 0.043460 |

Across these three direction-level means, the MAE range shrinks from **0.1212 to 0.0658 mm** isotropically and **0.1208 to 0.0200 mm** anisotropically. This is reduced variation for the tested directions, not rotational invariance. Phase and direction are not fully crossed in v2.

### 3. Does it perform better on anisotropic spacing?

Expanded scorable matrix: isotropic MAE **0.214197 → 0.081989 mm**, anisotropic **0.247146 → 0.078389 mm** (EDT → cross-sectional).

For the complete core cohort, avoiding different phase-failure composition, isotropic MAE is **0.190771 → 0.106690 mm**, anisotropic **0.211795 → 0.092083 mm**. Cross-sectional performs better in aggregate on both grids here. In-plane spacing also changes between these grids, so this does not isolate the causal effect of anisotropy. Anisotropic phase failures remove some diameters from expanded aggregates and must not be treated as accurate zero-error runs.

### 4–5. Curved and tapered agreement

Core geometry mean MAEs across both spacings and trims:

| Geometry | EDT MAE (mm) | Cross-sectional MAE (mm) |
|---|---:|---:|
| Straight | 0.209241 | 0.093269 |
| Oblique | 0.248489 | 0.091836 |
| Curved | 0.145833 | 0.098646 |
| Tapered | 0.193752 | 0.127116 |

Curved and tapered aggregate agreement improves. An important exception is the **2 mm isotropic curved phantom**, trimmed: EDT MAE **0.016425 mm**, cross-sectional **0.181633 mm**. EDT's agreement in that particular voxel arrangement is better; it is not evidence that either method is always preferable in curvature.

Taper correlation averages **0.9410 EDT** versus **0.9484 cross-sectional**, but MAE/RMSE remain the agreement measures. Area-equivalent cross-sectional diameter and nearest-boundary sphere diameter are different estimands in a taper. Mesh voxelization, local tangent choice and end-cap effects remain relevant.

### 6. Does it introduce invalid samples?

There are **0 invalid cross-sections out of 3,649 attempted sample slots** across both trim variants (1,977 untrimmed slots; the trimmed slots reuse retained positions). EDT is also valid on every recovered position. This is **conditional** on a recovered centerline.

There are **12 shared failed geometries out of 58 (20.7%)**, corresponding to 24 unscorable pairs. All ten half-voxel-shifted straight cylinders fail Lee skeletonization, as do the 1 mm and 5 mm quarter-shifted anisotropic cylinders. Their masks are nonempty. They are recorded rather than modified to ensure success. This exposes a dependency shared by both estimators; a new area estimator alone cannot solve it.

Controlled tests separately verify open, disconnected, empty, outside-origin, degenerate-tangent and mesh-failure behavior. Zero section failures in this simple matrix does not imply reliability on branches or noisy surfaces.

### Endpoint trimming

Expanded scorable results:

| Estimator | MAE, trim 0 (mm) | MAE, trim 2 (mm) | Aggregate RMSE, trim 0 → 2 (mm) |
|---|---:|---:|---:|
| EDT | 0.231125 | 0.228785 | 0.306555 → 0.300227 |
| Cross-sectional | 0.080225 | 0.080310 | 0.102834 → 0.103113 |

Among 46 geometries with centerlines, trimming improves/worsens/leaves unchanged the EDT MAE in **16/16/14** cases and cross-sectional MAE in **12/13/21**. It slightly helps EDT in aggregate and very slightly worsens cross-sectional. There is no basis for a universal trimming benefit. Trimming refers to recovered path endpoints and preserves the production short-path fallback.

### 7–8. Superiority and default recommendation

Cross-sectional performs better in most tested matched experiments, especially oblique cases, but **strengths remain geometry- and phase-dependent**. EDT wins the 2 mm isotropic curved phantom, 2 mm centered anisotropic straight cylinder, 1 mm quarter-phase isotropic straight cylinder, and 4 mm quarter-phase isotropic straight cylinder (both trim settings for each).

**Keep EDT as the current NeuroVasc production default for now.** This task provides evidence to advance cross-sectional as an experimental candidate, not to change existing case measurements. Before a default change, test phase/orientation combinations, realistic partial volume and segmentation perturbations, branching/multiple contours, taper/bend scales, and tangent-window sensitivity. Address the common skeleton failure mode as a separate priority. A later comparison could use an independently validated centerline or an analytic centerline only as a clearly labeled diagnostic experiment; it must not silently replace the pipeline under evaluation.

## Tests, artifacts and limitations

The suite has **59 passing tests**, including 35 existing tests, physical tangents and bases, a known analytic surface area, moderate binary-cylinder area error, oblique/anisotropic coordinates, affine rotation equivariance, curved geometry, explicit failures, matched-support logic, and all 52 frozen v1 numerical results. For the 4 mm isotropic binary cylinder, area is **12.125 mm²** versus analytic **12.566371 mm²** and equivalent diameter **3.929126 mm**. Regression limits of **0.6 mm² area error** and **0.12 mm diameter error** provide modest headroom over those observed errors; they are not clinical acceptance thresholds. An analytic 256-sided cylinder isolates contour integration and is required to match true area within 0.1%.

Pytest reports **55 upstream scikit-image/NumPy deprecation warnings** about direct array-shape assignment. They originate in the existing marching-cubes dependency, not failed assertions. They are not suppressed or patched in production code. The frontend retains its pre-existing bundle-size warning.

Under `validation/outputs/caliber_v2/`:

- `results.csv`: all 232 estimator rows, including failures.
- `matched_comparison.csv`: all 116 pairs, support counts, deltas and winners.
- `profiles.csv`: exact sample locations, raw estimates, areas, physical tangent vectors, analytic truth, signed errors and status (no rows are invented for missing upstream paths).
- `summary.json`: main/core/grouped results, coverage, shared failures, parameters, dependency versions and source hashes.
- `geometries.json`: phase/orientation, spacings, affines, shapes and analytic metadata.
- `masks/`: all 58 binary NIfTI phantoms, including those with skeleton failures.
- `figures/`: 13 separate Matplotlib PNGs, using the default color cycle. The required comparison figures are supplemented by both spacings for taper/curve profiles and an orientation-sensitivity plot.
- `v1_reproduction/`: isolated v1 command verification artifacts, leaving archived v1 untouched.
- `manifest.txt`: exact generated-file list.

The phase plot explicitly retains the half-voxel category with missing accuracy points and coverage counts; missing values are not plotted as zero. Its successful-run aggregates can contain different diameter subsets, so interpret it with the failure counts and per-diameter CSV rows.

Limitations: a binary 0.5 mesh is not a measured continuous lumen boundary; no partial-volume/noise experiment was added. Only one curvature radius, taper length and tangent window are tested. Cross-sectional area of noncircular vessels need not agree with inscribed-sphere diameter. Contours with multiple components/holes are rejected rather than clinically interpreted. Closed contours are not proof of correct branch selection. There is no clinical reference standard, no disease/stenosis inference, no patient validation, no confidence-interval claim, and no production method switch.
