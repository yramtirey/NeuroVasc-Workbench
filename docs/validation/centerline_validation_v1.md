# Centerline Validation v1

## Scope and conclusion

This is geometric/synthetic validation, **not clinical validation**. The unchanged Lee pipeline succeeds on 48/60 geometries (80%); the experimental clearance-weighted geodesic tree succeeds on 60/60 (100%). It recovers **all 12 known caliber-v2 upstream failures**, including all ten half-voxel cylinders. However, better extraction and coverage do not imply a more accurate centerline: on the 48 shared successful geometries, Lee has lower mean positional and tangent error. **Keep Lee as the production default for now.** The alternative is a candidate for further study, not an automatic fallback deployed to Case 001.

This implementation modifies no production call site, existing algorithm, frontend, API, Case 001 measurement, source data or archived caliber-v1/v2 artifact. The mesh cross-sectional caliber estimator remains experimental. No downstream diameter measurements are recomputed in this milestone.

## Reproduce

From the repository root in the existing environment:

```bash
source .venv/bin/activate
python3 scripts/run_centerline_validation.py
python3 -m pytest -q
```

No new dependency is needed. The runner also works by absolute path. `--output` must resolve to `validation/outputs/centerline_v1/` or a descendant; it rejects the caliber archives and paths outside this tree. Expected estimator failures are data and do not produce an unsuccessful process exit. Unexpected evaluation, serialization or plotting exceptions propagate as pipeline errors.

The previous commands were also executed with isolated destinations:

```bash
python3 scripts/run_caliber_validation.py --output validation/outputs/centerline_v1/regressions/caliber_v1
python3 scripts/run_caliber_validation_v2.py --output validation/outputs/centerline_v1/regressions/caliber_v2
```

Both reproduced archived `results.csv` and `profiles.csv` byte-for-byte; v2 also reproduced `matched_comparison.csv` exactly. Source metadata differs because new experimental modules exist. SHA-256 verification covered 250 pre-existing files captured before this task. Both archives, source data, existing analysis/backend code and Case 001 outputs are unchanged. The final check found 249 files unchanged and an external concurrent change to `frontend/src/App.tsx`; this task never wrote that file and left the concurrent edit intact. Regression evidence is retained in `regressions/checks.json`.

## Methods and endpoint inference

### Baseline

The adapter calls the existing `extract_centerline()`, `build_vascular_graph()` and `extract_main_path()` unchanged. Lee thinning provides skeleton voxels. The production graph connects all 26 neighbors with physical edge lengths. The main-path method selects the component with greatest total edge length, then the longest weighted shortest path between endpoints, with its existing endpoint-free fallback. Validation separately records full graph topology and rejects a disconnected recovered graph as an extraction failure; it does not alter production component selection.

### Experimental clearance-weighted geodesic tree

`neurovasc/geometry/geodesic_centerline.py` constructs a sparse graph containing **all foreground voxel centers** and their 26-neighbor edges. This is a geodesic/clearance method, not another call to skeletonization and not a subvoxel medial-axis reconstruction.

1. Threshold the input at 0.5 and require a nonempty connected 3D foreground. Add one background layer for EDT calculation, without modifying the input mask.
2. Calculate physical EDT radii `r_i`. EDT measures distance to background voxel centers with the supplied sampling, rather than continuous distance to the surface. See [SciPy EDT documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.distance_transform_edt.html).
3. Set physical neighbor distance `l_ij = norm(spacing * offset)` and clearance cost `c_ij = l_ij * (r_max / sqrt(r_i*r_j))**2`. Symmetric positive costs favor interior points while retaining path-length cost. Dijkstra is run on this sparse graph; its multi-source mode supplies distances to the current tree. See [SciPy Dijkstra documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.dijkstra.html).
4. Seed endpoint inference at a maximum-EDT voxel. Use physical-distance geodesic sweeps to find opposite terminal regions. Each terminal region is the connected patch in the farthest shell, whose width is one maximum voxel spacing. Recenter that patch using an EDT-weighted physical centroid, then snap to its nearest patch voxel. Numerical ties select lexicographic voxel order. Neither analytic endpoints nor phantom geometry labels enter the estimator.
5. Trace the clearance-weighted shortest path between the inferred endpoints. This path supplies the initial tree.
6. Compute physical distance from every foreground voxel to the tree. If the largest residual exceeds `max(3 mm, 2.5 * r_max)`, infer a residual terminal patch and graft its clearance-weighted path to the tree. Repeat, allowing at most eight grafts; otherwise fail explicitly. Only chosen path edges enter the output, preventing incidental 26-neighbor cross-connections from creating cycles.
7. Return a `VascularGraph`; reuse production `extract_main_path()` for an ordered principal path. No resampling, path smoothing, mask alteration or analytic-truth fitting is performed.

Graph weights use physical spacing, and coordinates use the complete affine, including rotation and translation. Orthogonal affine axes must agree with declared spacing; shear and inconsistent affines are rejected. Rotated-affine tests verify identical voxel nodes/edges and correctly transformed world coordinates. This tests coordinate handling, not invariance to revoxelizing a vessel at another orientation.

All parameters are recorded in `summary.json`. A 200,000-foreground-voxel limit prevents unbounded graph construction. Fixed ordering makes runs deterministic in the recorded environment. Equal-cost Dijkstra choices may change across dependency versions; versions and source hashes are recorded. The method forces a **tree**, so it cannot preserve a true anatomical loop. It is not ready to replace a Circle-of-Willis network extractor.

## Matrix and truth

There are **60 geometries and 120 untrimmed estimator runs**: the exact 58 caliber-v2 masks/affines plus two Y phantoms. Geometry IDs match v2 for direct joins with its failures. There are no paired trim variants in this experiment.

The reused matrix contains straight diameters 1–5 mm; oblique and curved diameters 2, 3, 4 mm; tapers 4→2 and 3→1.5 mm; spacings `(0.5,0.5,0.5)` and `(0.46875,0.46875,0.8)` mm; centered/quarter/half voxel phases; and oblique directions `(1,1,1)`, `(1,2,1)`, `(2,1,3)`. Length is 30 mm. The curve is the same radius-24 mm arc in x-z. Phase and orientation are not fully crossed.

The Y phantom is the union of three flat-capped binary tubes. Its analytic junction is the origin, parent length 22 mm and diameter 4 mm, and two daughter lengths 18 mm and diameters 3 mm. Daughters diverge ±35° from positive z in x-z. Diameters, lengths, angle and spacing are configurable. The truth network length is 58 mm and either parent-to-daughter principal path is 40 mm. It is generated once at each spacing, centered on the grid, and analyzed separately from the nonbranching cohort.

For straight/oblique/taper samples, truth is the closest point on the **finite physical axis segment**; beyond a cap, longitudinal endpoint error therefore also counts. Curved truth is obtained by projecting world coordinates onto the finite analytic arc with clamped arc angle. No matching by voxel or sample index is used.

For a Y main path, endpoint positions identify the two analytic branches being traversed. Samples project to that branch pair, excluding the unused third branch even near the junction. Branch-local evaluation first assigns recovered endpoints one-to-one to analytic endpoints by minimum total physical distance; paths from the recovered logical junction to those endpoints are then compared with their assigned finite analytic branch. Analytic truth is used only after extraction.

## Metrics and success definitions

`extraction_success` requires a connected finite graph and a nonrepeating ordered path with at least two points, matching graph coordinates, actual graph adjacency, positive physical edge lengths and correct cumulative distances. Failure rows remain in `results.csv` with zero recovered main-path samples and explicit reasons; inapplicable accuracy values are blank/null, never zero. Graph diagnostics are retained when available. `path_found` is distinct from valid extraction.

Position error is Euclidean distance from each recovered sample to its projected truth point. Rows report mean, median, maximum and RMS error. Length is the unsmoothed voxel-center polyline length, compared with full analytic length; absolute, signed and percentage errors are retained. Voxel stair steps can inflate length even when coverage is incomplete.

Coverage is the span of projected analytic positions divided by full truth length, accompanied by normalized start/end positions and monotonic-order checks. Reversing a correctly ordered path is allowed. A simple path passes `correct_principal_path` only with the expected topology, monotonic projected ordering, and endpoints within the outer 10% at each end. Thus a short favorable segment is not treated as a fully recovered vessel. These coverage criteria are descriptive validation thresholds, not clinical requirements. `failure_inclusive_mean_coverage` also includes failed extractions as zero **coverage**, while their accuracy remains undefined.

A Y principal route is parameterized from the first matched branch tip through the analytic junction to the second tip. It must connect parent to a daughter, have correct topology, monotonic ordering and at least 80% coverage. Branch rows separately report `branch_path_found`; `branch_recovered` additionally requires at least 80% branch coverage and endpoint error at most two maximum voxel spacings. This distinction prevents correct topology from hiding inaccurate branch endpoints. A nonqualifying branch retains its measurements and reason.

Tangents reuse caliber-v2's physical 6 mm least-squares window on the full recovered path, with clipped one-sided endpoint fits. Angular error is `acos(abs(dot(unit_estimate, unit_truth)))`, so path reversal is not a 180° error. Branch and Y-route tangents within 3 mm of the analytic junction are excluded; valid tangent counts remain explicit. Even an exact curved path has finite-window endpoint angular bias.

Endpoints have degree one. Degree≥3 nodes form logical junctions by connected components of their induced subgraph. Raw junction-voxel counts are also saved. Independent cycles are `E - N + components`; they are not silently collapsed away. The experimental tree has zero cycles by construction, not because it has demonstrated correct loop anatomy.

All averages give equal weight to geometries, not to sample counts. Own-successful-support results and common-successful-geometry results are separate. Methods do not recover identical sample positions or extents, even on common geometries; this comparison is not an identical-coordinate paired error test. No missing interval is assigned invented positional error. Coverage and extraction success must accompany accuracy.

Matched wins use declared tie tolerances: 0.01 mm positional mean error, 0.01 coverage fraction, 0.1 mm length error and 1° tangent error. These are descriptive bins; the raw differences remain in CSV. They are not confidence intervals or statistical superiority tests.

## Recorded results

### Extraction and known failures

| Cohort | Lee | Experimental |
|---|---:|---:|
| All geometries | 48/60 (80%) | 60/60 (100%) |
| Exact nonbranching caliber-v2 cohort | 46/58 (79.31%) | 58/58 (100%) |
| Y phantoms | 2/2 | 2/2 |
| Known v2 failures | 0/12 | 12/12 |

All ten half-phase straight cylinders and both quarter-phase anisotropic cylinders (1 and 5 mm) still fail unchanged Lee thinning. Their masks are nonempty; all are recovered by the alternative. There are 48 both-successful pairs, 12 Lee-failed/alternative-successful pairs, no reverse failures and no both-failed pairs. Every known failure has a diagnostic showing binary foreground, absent Lee result, alternative and analytic truth.

47/60 Lee runs and 60/60 alternative runs meet the stricter principal-path criterion. The one otherwise-successful Lee exception is the 4 mm anisotropic `(1,1,1)` cylinder: coverage 81.20%, normalized endpoints 0.0818 and 0.8938; one end misses the outer-10% criterion. This is separate from its successful extraction.

### Position, coverage, length and tangent tradeoffs

On the **same 48 successful geometries**:

| Equal-geometry metric | Lee | Experimental |
|---|---:|---:|
| Mean position error (mm) | 0.127582 | 0.143830 |
| Mean coverage | 89.71% | 98.19% |
| Mean absolute length error (mm) | 2.262111 | 1.517388 |
| Mean tangent error (degrees) | 0.906520 | 2.553158 |

Lee has better average position and tangent accuracy; the alternative has better coverage and length accuracy. Positional wins are Lee 19 / alternative 13 / ties 16. Coverage wins are 0 / 45 / 3; length wins 15 / 30 / 3; tangent wins 22 / 1 / 25.

Across each method's own successful support (48 vs 60), alternative mean positional error is 0.184467 mm, coverage 98.26%, length error 1.273971 mm and tangent error 2.448944°. Its additional recovered phase cases are harder positional tests, so these means must not be confused with the common-cohort comparison. Including failed extractions only in coverage gives 71.77% Lee versus 98.26% alternative over all 60 geometries.

### Grid phase

For straight cylinders only, both methods succeed on all ten centered cases. Quarter phase gives 8/10 Lee versus 10/10 alternative; half phase gives 0/10 versus 10/10. Mean positional errors for centered/quarter/half phases are approximately `0 / 0.2256 / undefined` mm for Lee and `0.0387 / 0.2436 / 0.3619` mm for the alternative. Quarter-phase means have different successful diameter support.

The alternative fixes the empty-path failure but remains voxel-phase sensitive. On the half-phase 1 mm isotropic cylinder, the nearest foreground axes lie 0.353553 mm from the analytic axis; the recovered mean error equals that value. Its coverage is 29.5/30 = 98.33%, and its length error is 0.133975 mm. A voxel-center method cannot exactly recover a centerline lying between those axes.

### Spacing and orientation

Lee success is 25/30 isotropic (83.33%) and 23/30 anisotropic (76.67%); alternative success is 30/30 on each. On common successful geometries within each spacing:

| Spacing | Lee position / alternative (mm) | Lee coverage / alternative | Lee length error / alternative (mm) | Lee tangent / alternative (°) |
|---|---:|---:|---:|---:|
| Isotropic | 0.105707 / 0.153411 | 91.61% / 98.99% | 1.644621 / 1.446694 | 0.383753 / 2.895518 |
| Anisotropic | 0.151360 / 0.133415 | 87.65% / 97.31% | 2.933296 / 1.594230 | 1.474746 / 2.181028 |

The alternative improves mean position on the shared anisotropic subset, but worsens it on the isotropic subset. Both in-plane and through-plane spacing change, and the surviving diameter cohorts differ: this is not an isolated causal estimate of anisotropy.

All 18 oblique cases succeed with both methods. Means across diameters and spacings:

| Direction | Lee position / alternative (mm) | Lee length error / alternative (mm) | Lee tangent / alternative (°) |
|---|---:|---:|---:|
| (1,1,1) | 0.116061 / 0.139299 | 2.313753 / 1.413207 | 1.199473 / 3.167595 |
| (1,2,1) | 0.246627 / 0.212075 | 1.703385 / 3.870436 | 1.701389 / 2.612674 |
| (2,1,3) | 0.203156 / 0.231241 | 1.367762 / 2.418959 | 2.495530 / 4.981788 |

The alternative improves coverage in every orientation group (roughly 98% versus 87–88%), but not position, length or tangent universally. In particular, voxel polyline length inflation is substantial for direction (1,2,1).

### Curves and tapers

Both methods succeed on all six curves and all four tapers.

| Geometry | Lee position / alternative (mm) | Lee coverage / alternative | Lee length error / alternative (mm) | Lee tangent / alternative (°) |
|---|---:|---:|---:|---:|
| Curved | 0.134765 / 0.144535 | 90.34% / 98.23% | 1.883845 / 1.164591 | 1.660344 / 2.889945 |
| Tapered | approximately 0 / 0.017858 | 91.50% / 98.00% | 2.550000 / 0.834789 | 0 / 0.818948 |

Longer recovered extent helps coverage and absolute length error, while off-axis terminal choices and voxel steps can worsen local position and tangent estimates. This experiment does not separate terminal-region effects from interior error by trimming.

### Y topology and branches

Both methods recover three endpoints, one logical junction, zero cycles, and parent-to-daughter principal routes at both spacings. All six branch paths per method exist. This is topology recovery, not perfect geometric recovery.

Junction errors (isotropic, anisotropic) are **0, 0.8 mm Lee** and **0.707107, 0 mm alternative**. Principal lengths versus 40 mm truth are **38.4350, 36.4961 mm Lee**, and **42.0813, 40.7608 mm alternative**. Mean principal coverage is 91.87% versus 98.99%. Mean principal positional error is 0.057940 versus 0.096139 mm.

The additional branch coverage/endpoint criterion is met by **3/6 Lee branches and 5/6 alternative branches**. Lee misses the isotropic daughter endpoint tolerance (1.493 mm errors versus 1 mm allowed) and the anisotropic parent tolerance (2.8 mm versus 1.6 mm). The alternative misses the isotropic parent tolerance (1.118 mm versus 1 mm). Branch-local tangent errors, lengths, coverage and endpoint errors remain in `branch_results.csv`; no failed criterion is discarded. Only the alternative anisotropic Y passes all three branch criteria simultaneously.

## Failure modes, limitations and recommendation

The alternative has no extraction failures in this small matrix, but introduces meaningful limitations: off-axis cap endpoints; voxel stair-step length inflation; less stable endpoint tangents; memory cost of a full interior graph; explicit rejection of disconnected masks, shear and oversized graphs; and a radius-dependent branch threshold that can discard short/narrow branches near a wide parent. It constructs a tree and cannot preserve true vascular loops. The Y example does not test complex branching, junction angle/diameter sweeps, nearby vessel contacts or endpoint ambiguities.

The graph's 26-neighbor diagonal edges connect foreground voxel centers; this does not prove that an entire diagonal segment lies inside a reconstructed continuous lumen. EDT ridge preferences inherit voxel-center distance bias. A recentered terminal patch is still a heuristic, particularly for tapering, curved caps and asymmetric branches. There is no subvoxel centering or graph-parameter sensitivity experiment. Results are deterministic examples, not independent clinical replicates or estimates with statistical confidence bounds. No imaging blur, partial volume, segmentation perturbation, pathology or patient reference standard is tested.

**Retain the unchanged Lee production method and keep this alternative experimental.** Next work should test endpoint refinement and physically regularized paths, window and clearance/branch-threshold sensitivity, loops and richer branch phantoms, segmentation perturbations, and downstream caliber accuracy on both recovered paths. Do not infer that recovering all twelve failures automatically improves the downstream cross-sectional estimator. Any eventual fallback or default switch needs explicit review and further validation.

## Tests and artifacts

**89 tests pass**: the existing 59 plus 30 new tests. New tests cover all analytic centerlines, branch truth, physical projection and short-path coverage, tangent sign invariance, baseline equivalence, all known Lee failures, deterministic endpoint selection, anisotropic/rigid-affine coordinates, empty/disconnected/sheared/oversized inputs, malformed graph/path rejection, cycles/logical junctions, exact matrix mask/affine reuse, failure accounting and output-directory protection. Existing caliber regression tests remain intact. The 55 existing scikit-image/NumPy deprecation warnings remain unsuppressed.

Regression limits have physical meaning: the thin half-phase position error has a known voxel-center lower bound; its length error must remain below one 0.5 mm voxel; Y junction error must remain within one spatial voxel diagonal (observed worst 0.7071 mm isotropic); analytic curved tangents must stay within the finite-window endpoint bound. These are synthetic regression guards, not clinical acceptance thresholds.

Under `validation/outputs/centerline_v1/`:

- `results.csv`: 120 rows, with extraction, topology, coverage, position, length and tangent measurements, including 12 explicit Lee failures.
- `profiles.csv`: actual recovered main/branch samples in world mm, projected truth, analytic branch/position, point error and tangent error. Branch rows overlap main-path samples; do not pool them as independent replicates.
- `branch_results.csv`: 12 rows (three branches × two spacings × two methods), including branch-quality failures.
- `matched_comparison.csv`: 60 pairs, extraction outcomes, metric differences and declared-tolerance winners.
- `summary.json`: counts, own/common support aggregates, geometry/phase/orientation/spacing/cohort groups, all known failure pairs, parameters, dependencies and source hashes.
- `geometries.json` and `masks/`: complete metadata/affines and all 60 binary NIfTI masks.
- `figures/`: 28 separate Matplotlib figures: 16 summary/example plots and 12 known-failure diagnostics, using the default color cycle. Summary error plots show each estimator's successful support with failure/coverage context; the common-support table above is the fairer aggregate comparison.
- `manifest.txt`: exact inventory of the 95 primary artifacts, including itself. Nested regression artifacts have their own manifests; they are not attributed to the centerline runner.
- `regressions/`: isolated caliber reruns and preservation/comparison evidence.

New source files are `neurovasc/geometry/geodesic_centerline.py`; `validation/synthetic_centerline/{__init__,geometries,truth,evaluate,report,figures}.py`; `scripts/run_centerline_validation.py`; `tests/test_centerline_validation.py`; and this document. No pre-existing source file was modified by this implementation; the separate concurrent App.tsx edit is noted above. Nothing was pushed to GitHub.
