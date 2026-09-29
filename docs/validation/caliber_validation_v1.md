# Caliber Validation v1

## Purpose and scope

This is geometric/synthetic validation of a segmentation-derived caliber estimator. It is **not clinical validation**, a diagnostic accuracy study, or evidence that a caliber reduction represents disease. No patient/case data are used in these experiments. Production analysis code, frontend, backend and Case 001 outputs are unchanged.

## Reproduce

From the repository root, using its existing Python 3.14 environment:

```bash
source .venv/bin/activate
# Only needed when setting up the test dependency:
python3 -m pip install '.[validation]'
python3 scripts/run_caliber_validation.py
python3 -m pytest -q
```

The script works when launched by absolute path from another working directory as well. Outputs default to `validation/outputs/caliber_v1/`; `--output` may select another directory strictly inside `validation/outputs/`. Re-running replaces only validation artifacts. No randomness or timestamp changes are introduced. Package versions and source SHA-256 hashes are recorded in `summary.json`. The run takes a few seconds after the first Matplotlib font-cache initialization on this development Mac.

## Production method tested

The adapter calls the unchanged production functions:

1. `extract_centerline`: Lee skeletonization of a binary volume.
2. `build_vascular_graph` and `extract_main_path`: ordered longest terminal path in the largest component by physical edge length.
3. `build_diameter_profile`: twice the physical EDT radius at centerline voxels.
4. `analyze_profile`: existing Gaussian smoothing (sigma 1.25 **samples**) for secondary comparison/plotting.

SciPy EDT measures distance to the nearest background **voxel center**, with `sampling` supplying physical spacing. It is not an exact distance to the continuous segmentation boundary. See the [SciPy EDT documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.distance_transform_edt.html). The production skeletonization is a voxel-grid thinning algorithm; see [scikit-image's Lee skeletonization description](https://scikit-image.org/docs/stable/auto_examples/edges/plot_skeleton.html). Passing physical spacing to EDT does not make skeletonization itself spacing-aware.

The primary quantitative columns evaluate the **raw** production diameter. `smooth_mae_mm` and saved smoothed profiles permit secondary inspection. Smoothing is not substituted for the production estimator to improve headline accuracy.

## Analytic phantoms and experiment matrix

All lengths are 30 mm, with positive background padding, an odd-sized grid, and physical origin on a voxel center. Inclusion is evaluated at voxel centers, with the analytic boundary included. Affines encode physical mm coordinates. There is no resampling or partial-volume model. Flat caps are perpendicular to the analytic axis/tangent. NIfTI masks and metadata are saved for every phantom.

| Geometry | Diameters (mm) | Definition | Geometries across both spacings |
|---|---|---|---:|
| Straight | 1, 2, 3, 4, 5 | Cylinder along z | 10 |
| Oblique | 2, 3, 4 | Cylinder along normalized (1,1,1) | 6 |
| Curved | 2, 3, 4 | Circular arc in x-z, bend radius 24 mm, arc angle 1.25 radians | 6 |
| Tapered | 4→2, 3→1.5 | Linear diameter change along z | 4 |

Both `(0.5, 0.5, 0.5)` and `(0.46875, 0.46875, 0.8)` mm are tested. Each of 26 phantoms is evaluated with `trim_mm=0` and `trim_mm=2`: **52 paired runs**, not 52 independent statistical replicates. All 52 succeeded. No bifurcation is included in v1.

For taper truth, `D(u)=D_proximal + u*(D_distal-D_proximal)`. Normalized position is projected from each sample's **world coordinates** onto the known physical axis (arc angle for curved phantoms). It is not inferred by dividing the recovered path distance by its length. Thus path reversal, skeleton endpoint shortening, and trimming do not shift the truth curve. `profiles.csv` preserves both recovered path distance and analytic normalized position, including signed sample errors.

## Metrics and aggregation

Each results row reports mean/median/min/max raw diameter; MAE; RMSE; maximum absolute sample error; signed bias; absolute mean bias (`absolute_error_mm`); mean absolute percentage error (`percent_error`); signed percentage bias; and sample count. Constant-geometry `true_diameter_mm` is blank for tapers, which instead have proximal/distal diameters. Correlation is omitted for constant truth or constant predictions; it is computed for tapered profiles when defined. Correlation is not an agreement measure.

MAE is `mean(abs(estimate - truth))`; RMSE is `sqrt(mean((estimate-truth)^2))`. `absolute_error_mm=abs(mean(estimate-truth))` can be smaller than MAE when errors cancel. Percent error is `100*mean(abs(estimate-truth)/truth)` at corresponding samples.

Overall MAE is the equal-run mean of run MAEs. Overall RMSE is the square root of the equal-run mean of run MSEs, **not the mean of run RMSEs**. Both trim settings receive equal weight. No confidence intervals are inferred from deterministic phantoms. Failures, if any, are retained in results and excluded explicitly from aggregates; the runner exits unsuccessfully after writing diagnostics. `mae_delta_vs_straight_mm` compares each constant geometry with its same-diameter, same-spacing, same-trim straight baseline.

## Results from the recorded v1 run

- Overall MAE: **0.201283 mm**; overall RMSE: **0.252097 mm**.
- Isotropic MAE: **0.190771 mm**; Case-001-like anisotropic MAE: **0.211795 mm**.
- Equal-run signed bias: **+0.093299 mm**. There are 35 positive-bias runs and 17 negative-bias runs.
- Largest run MAE and percentage error: **1 mm straight cylinder, isotropic, trim=2 mm**. Mean estimate **1.414214 mm**, MAE **0.414214 mm**, percentage error **41.4214%**.
- Largest individual absolute sample error is a separate result: **1.0 mm**, in the **4 mm isotropic straight cylinder, untrimmed**.
- The moderate 4 mm isotropic straight cylinder with trimming has mean **4.123106 mm**, MAE **0.123106 mm**, or **3.0776%**.
- Taper raw-versus-truth correlations range from **0.9082 to 0.9657**; stair-step errors remain despite the strong trend correlation.

Mean MAE / signed bias by geometry and spacing (both trim settings):

| Geometry | Isotropic MAE (mm) | Isotropic bias (mm) | Anisotropic MAE (mm) | Anisotropic bias (mm) |
|---|---:|---:|---:|---:|
| Straight | 0.20972 | +0.19548 | 0.20871 | +0.20871 |
| Oblique | 0.22527 | +0.22527 | 0.27171 | -0.17116 |
| Curved | 0.11046 | -0.10344 | 0.18121 | -0.06850 |
| Tapered | 0.21200 | +0.20367 | 0.17550 | +0.17550 |

There is a modest overall tendency to overestimate in this matrix, **not a universal systematic sign**. Curved phantoms generally underestimate, and the oblique group changes mean bias sign with spacing. Curvature does not uniformly worsen MAE: isotropic curved MAE is lower than its straight baseline for 2 and 3 mm, and higher for 4 mm. Anisotropic curved behavior also varies by diameter. The baseline-difference CSV column supports these comparisons without averaging away diameter effects.

### Voxel spacing

The Case-001-like spacing increases aggregate MAE by **0.021024 mm (~11.0%)**. Oblique and curved groups worsen, while the tapered group improves and the straight group is nearly unchanged. Both in-plane and through-plane spacing differ between the two grids, so this is **not an isolated causal estimate of anisotropy**. A single grid alignment/orientation cannot establish general performance across acquisitions.

### Endpoint trimming

| Requested trim | Mean MAE (mm) | Aggregate RMSE (mm) | Signed bias (mm) |
|---|---:|---:|---:|
| 0 mm | 0.202009 | 0.255936 | +0.090574 |
| 2 mm | 0.200556 | 0.248198 | +0.096023 |

Trimming improves MAE in **8** phantoms, worsens it in **10**, and leaves it unchanged in **8**. Aggregate MAE decreases by only **0.001453 mm (~0.72%)**. All selected phantoms have enough samples for trimming to apply. The framework records actual sample removal and `trim_applied`; the production fallback for short paths remains unchanged.

The production trim removes 2 mm from the **recovered main-path endpoints**, not the analytic caps, and rebases retained profile distance to zero. Skeletonization may already shorten the path. Removing low endpoint estimates can expose a positive interior bias and increase MAE, as in the 1 mm cylinder. This prevents interpreting trimming as a universal accuracy correction.

## Artifacts

Under `validation/outputs/caliber_v1/`:

- `results.csv`: 52 run rows with trim status, metrics and straight-baseline differences.
- `profiles.csv`: sample coordinates, analytic truth, raw/smoothed estimates and signed errors.
- `summary.json`: counts, macro metrics, worst run, grouped results, versions and source hashes.
- `geometries.json`: analytic dimensions, axis, bend radius, affine-related spacing and voxelization conventions.
- `masks/`: 26 NIfTI binary phantoms.
- `figures/`: 9 separate PNG figures (no subplots; Matplotlib default color cycle).
- `manifest.txt`: exact generated-file inventory.

The first three figures show estimated versus true diameter, MAE versus true diameter, and percentage error versus true diameter for constant-caliber geometries at trim=2 mm. There are four tapered-profile plots and two 3 mm curved-profile examples, covering both spacings. Untrimmed data remain available in both CSVs.

## Automated regression tests

**35 tests passed** in the initial implementation. Tests cover physical bounds, metadata, NIfTI round-trip, background padding, analytic truth under path reversal, oblique/arc projections, every matrix geometry with both trim settings, trim sample membership, error formulas, invalid inputs, and a moderate-cylinder accuracy guard.

The 4 mm isotropic trimmed cylinder has observed MAE 0.123106 mm. The initial regression limit is **0.25 mm MAE and absolute mean error**, half the voxel spacing (6.25% of diameter). It allows headroom for dependency-level numerical changes while catching meaningful regressions such as radius/diameter or spacing mistakes. It is specific to this phantom and **is not a general accuracy requirement or clinical acceptance threshold**. A sample-count and path-coverage check prevents a trivially shortened path from passing on a few favorable samples.

## Limitations and next comparison

- Discrete distance to background centers differs from continuous distance to a surface. Thin 1 mm structures have only about two in-plane voxels across and are especially sensitive to voxelization.
- One centered grid phase, one diagonal, one curvature radius, and one taper length are tested. These can favor particular digital cross-sections. More phases and orientations are needed before broad bias conclusions.
- Skeleton thinning is not guaranteed to give the physical medial axis in anisotropic voxels. Off-axis samples may underestimate caliber even when EDT spacing is correct.
- EDT samples a nearest-boundary radius, not cross-sectional area. In a taper, even a continuous inscribed-sphere diameter differs slightly from the perpendicular cross-sectional diameter. Noncircular and branching vessels will introduce additional estimand differences.
- No segmentation errors, partial volume, imaging blur, noise, vessel-wall model, motion, bifurcations, or reference-image uncertainty are represented.
- Smoothing uses sample-index sigma, not a fixed physical distance; its physical effect changes with path sampling and spacing.
- Main-path selection and trimming can discard endpoints or branches. Recorded normalized coverage describes only the recovered samples; missing segments are not given zero error.

Recommended next step: compare a **physical-space cross-sectional equivalent diameter** adapter against the unchanged EDT baseline, alongside a phase/orientation sweep and controlled partial-volume voxelization. Do not apply a blanket EDT correction from this small matrix. `evaluate(..., estimator=..., estimator_name=...)` accepts an alternate adapter returning the production `DiameterProfile` shape, preserving common coordinate-based truth and metrics. No VMTK dependency or alternate estimator is implemented in v1.
