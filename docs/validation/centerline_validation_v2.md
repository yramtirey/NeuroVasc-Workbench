# Centerline Validation v2: experimental hybrid and physical refinement

## Scope and decision

This is synthetic geometric validation, **not clinical validation**. The experiment preserves the 60 Centerline-v1 geometries and adds four loop networks. The hybrid retains Lee whenever mask-only QC passes and uses refined geodesic fallback otherwise. On the v1 cohort it retains all 48 successful Lee paths **without altering their coordinates**, recovers all 12 known failures, and improves failure-inclusive coverage from **71.77% to 91.73%**. Refinement substantially improves raw geodesic position and smoothness in aggregate, but cap localization and tangent errors remain important exceptions.

**Keep production unchanged and do more validation before deployment.** An experimental fallback now exists, but it is not installed as a production fallback. Its tree topology is unsafe for general Circle-of-Willis recovery. Case 001, frontend/backend, default `extract_centerline()`, EDT caliber and the experimental status of cross-sectional caliber are unchanged. No pre-existing file was edited by this implementation.

## Reproduce and preservation

From the repository root:

```bash
source .venv/bin/activate
python3 scripts/run_centerline_validation_v2.py
python3 -m pytest -q
```

No new dependency is required. The runner also accepts an absolute script path. Output defaults to `validation/outputs/centerline_v2/`; `--output` must resolve to this directory or a descendant. Caliber-v1/v2 and Centerline-v1 directories are rejected. Expected extraction/refinement failures remain data; unexpected pipeline errors produce an unsuccessful exit. No timestamps or randomness are introduced. Source SHA-256 hashes, dependency versions, QC thresholds and every sweep configuration are recorded.

The previous caliber commands were executed with `--output validation/outputs/centerline_v2/regressions/caliber_v1` and `.../caliber_v2`. The unchanged Centerline-v1 runner only permits destinations inside its own archive tree, so it was executed in an isolated temporary mirror of the source tree using the same Python environment, then its results copied into `centerline_v2/regressions/centerline_v1/`. This exercises the original runner without changing its output guard or writing into its archive. The original v1 command remains `python3 scripts/run_centerline_validation.py`.

All nine comparison CSVs reproduce their archives **byte-for-byte**: results/profiles for each caliber version, caliber-v2 matched comparisons, and Centerline-v1 results/profiles/matched comparisons/branch results. Source metadata differs because new modules exist. Preservation hashes cover 476 pre-existing files, including production, source data, Case 001 and all three archives; evidence is saved under `regressions/checks.json`.

## Architecture and selection

`Recovery` separates the full underlying NetworkX graph from an ordered principal path. Node identifiers are opaque; refined coordinates need not be integer voxel centers. A Lee graph retains all existing edges, including cycles. The principal path is never substituted for the full network when evaluating loops.

The experimental hybrid performs:

1. Run unchanged Lee skeletonization, production graph construction and main-path extraction.
2. Apply mask-only QC to that graph/path.
3. If critical QC passes, return the original Lee recovery unchanged.
4. Otherwise run the unchanged experimental v1 clearance-weighted geodesic tree, refine it, and apply the same QC.
5. Return selected method, Lee extraction/QC results, fallback reason, geodesic/refinement success, final QC and diagnostics. Failed refinement/final QC is explicit; there is no silent switch back to raw geodesic.

Selection functions accept a `Volume` and engineering parameters, never a synthetic-vessel object, analytic endpoints or truth errors. Tests monkeypatch analytic projection to raise during selection and check that none of these modules imports validation truth. Optional cached recoveries in the runner avoid recomputation, without changing the decision rule.

`supports_cycles=False` is attached to raw/refined fallback. Hybrid diagnostics always expose `fallback_supports_cycles=False` and `production_eligible=False`. The selected Lee graph can preserve cycles; this does **not** mean the fallback can. No production call site imports the hybrid.

## Explicit QC criteria

The headline experiment uses **network mode on every geometry**, rather than selecting modes from truth labels. Simple mode is also available and tested. Every executed check records criterion, observed value, threshold, pass/fail and whether it is critical in `qc_results.csv`. Checks requiring coordinates are omitted when extraction fails. `summary.json` includes numeric thresholds. Inclusive numerical comparisons allow `1e-8` roundoff tolerance, preventing a computed 60.00000000000001° from incorrectly failing a 60° limit.

| Criterion | Threshold | Engineering rationale |
|---|---|---|
| Extraction exists | Nonempty graph and ordered path with at least 2 points | Establish a recovery to inspect |
| Finite world coordinates | Every graph/path coordinate finite | Prevent unusable downstream geometry |
| Connected components | Exactly 1 | This experiment expects one connected vessel/network |
| Minimum path samples | At least 5 | Avoid treating a minimal fragment as a supported local path |
| Repeated path samples | None | Reject degenerate ordering and zero-step fits |
| Edge lengths | Finite and >1e-8 mm | Reject collapsed/invalid edges |
| Ordered graph path | Consecutive coordinates correspond to actual graph edges | Check path/network consistency |
| Physical edge consistency | Euclidean world length, relative tolerance 1e-5 and absolute 1e-8 mm | Verify physical units and edge bookkeeping |
| Extent coverage proxy | At least 0.80 | Detect substantially truncated recovery |
| Terminal extent gap | At most 0.15 of mask extent | Detect one-sided shortening even with a plausible length |
| Simple-mode cycles | Zero; network mode permits cycles | Avoid imposing tree topology on network-mode inputs |
| Mean local turn | At most 60° | Flag severe voxel-scale zigzagging |
| Maximum tangent step | At most 45° | Flag abrupt fitted direction changes |
| Path length / mask extent | At most 1.8 simple, 4.0 network | Broad guard against path inflation, allowing curved network paths |
| Outside path samples | Zero, using interpolated mask ≥0.5 | Reject sampled excursions from the declared lumen model |
| Cycle capability | Reported, noncritical | A capability flag, not proof of anatomical topology |

The extent proxy projects **all foreground centers and all recovered graph nodes** onto the first physical PCA axis of the mask. Coverage is recovered projection span divided by foreground span. The terminal-gap proxy is the larger missing projection extent at either end, normalized by foreground span. These proxies avoid analytic caps and remain meaningful for branches, but they cannot identify missed loops or guarantee accurate cap centers. The tangent-step check uses the existing 6 mm physical tangent fit. Thresholds were fixed before the final comparison; they are transparent heuristics, not learned cutoffs or clinical acceptance criteria.

All successful Lee results pass this QC, including the v1 4 mm anisotropic oblique case whose analytic coverage is only 81.2%. QC is deliberately a broad structural screen, not an oracle for geometric accuracy. Raw geodesic trees on loops also pass these numerical checks despite losing the true cycle; the separate topology comparison exposes that limitation.

## Physical refinement algorithm

The raw graph is split into maximal degree-two chains between terminals/junctions. Each chain is refined separately and reconnected with shared junction identities. Original junction coordinates and the 3 mm path neighborhood around each junction are protected. The algorithm verifies tree topology after rebuilding; it explicitly rejects a cyclic graph rather than cutting it into a tree silently.

The declared headline configuration is `center1_step0.5_smooth1_end1`, fixed before evaluating the sweep:

- **Physical resampling:** piecewise-linear interpolation in cumulative world-space arc length, with `ceil(length / requested_spacing)` equal intervals. Both endpoints are retained. Requested spacing is 0.5 mm, near the smaller input voxel dimensions. A second resampling occurs after refinement; straight-line chord intervals on curves can be slightly shorter than the parameter interval.
- **Local recentering:** estimate a 6 mm-window physical tangent. Find nearby foreground centers within `2 * local_EDT_radius + max_spacing`. Compute a centroid weighted by physical EDT radius and an axial Gaussian with sigma 0.5 mm. Apply only the displacement transverse to the tangent. This is a local clearance-weighted centroid, **not** a true mesh cross-sectional area centroid or exact medial axis.
- **Regularization:** fit a local quadratic polynomial to xyz as a function of physical arc length. The Gaussian weighting sigma is 1 mm with support ±2 mm; evaluate the fitted intercept at each sample. Endpoint windows are one-sided. This uses physical distances, not sample-index smoothing.
- **Motion constraint:** each recentering or regularization proposal is bounded by two maximum voxel spacings. Backtrack by factors 1, 1/2, 1/4, 1/8, 0 toward its reference point until it satisfies the lumen constraint. Regularization is bounded relative to the initial resampled path; endpoint extension below is a separate bounded operation.
- **Endpoint experiment:** after centering/smoothing, estimate the local tangent and extend each true graph terminal outward in 0.1 mm increments, at most two maximum voxel spacings. Stop at the first interpolated-mask value below 0.5. Junctions are never extended. This is a segmentation-boundary heuristic, not projection onto an analytic cap.
- **Lumen checks:** transform world coordinates through the inverse affine and sample the binary mask by trilinear interpolation (`order=1`, zero outside the array). The interpolated field defines the ≥0.5 constraint; see [SciPy map_coordinates](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.map_coordinates.html). Every final chain is additionally sampled at one quarter of the minimum voxel spacing. Failure is explicit if any sampled segment leaves this model.

Dense containment checks are still sampled checks, not a continuous geometric proof, and this interpolated binary-mask model is not a measured lumen boundary. Positive edges, finite coordinates and duplicate-free paths are also checked. Refinement has no access to analytic truth.

Displacement statistics compare refined points with the raw chain at the same **normalized arc fraction**. They include reparameterization and endpoint-extension effects; they are not pure nearest-point motion. Junction protection may preserve a junction that was already misplaced in the raw tree.

## Experiment matrix and sweep

There are **64 geometries**, **256 headline rows** (four methods), and **960 refinement-sweep rows** (15 configurations × 64). The unchanged 60 v1 IDs, masks and affines remain directly joinable. Headline accuracy tables use those 60 geometries; loop-network metrics are separate because their estimand differs from principal-path metrics.

The four additional geometries are a radius-8 mm circular ring and that ring with two 8 mm radial branches, each at isotropic and anisotropic spacing. Tube diameter is 3 mm, the ring lies in x-z, and its center is the physical origin. Analytic network lengths are 50.265482 mm and 66.265482 mm. Both have one true cycle; endpoint/junction counts are respectively 0/0 and 2/2. Binary voxel-center inclusion and background padding remain deterministic.

The controlled sweep contains:

- Recenter off/on × resample spacing 0.5/1.0 mm × smoothing sigma 0/1/2 mm, with endpoint extension enabled: 12 combinations.
- The headline candidate with endpoint extension disabled: one matched control.
- Resampling alone at 0.5 and 1.0 mm, with recentering, smoothing and extension disabled: two controls.

All 15 variants are retained, including failures if any. All 960 attempted refinements succeeded in this matrix. No per-geometry winner is used for the hybrid. The candidate was not replaced by whichever variant minimized a reported error.

## Evaluation rules

The v1 physical axis/arc/branch projections, positional summaries, finite-cap treatment, coverage, length errors and sign-invariant tangent errors are reused. Y branch matching and quality thresholds are unchanged. Endpoint error uses one-to-one matching of recovered graph terminals to analytic terminals; unmatched terminal counts remain visible through topology. Interior positional error excludes 3 mm from each recovered main-path endpoint. This region can change when endpoints move, so the endpoint-control interior comparison is not a perfectly identical-location test.

Smoothness diagnostics include cumulative, mean and maximum local turning angle; turning radians divided by adjacent mean segment length as a curvature proxy; path-length/endpoint-chord ratio; and maximum step in the fitted tangent. These are geometric diagnostics, **not clinical tortuosity measures**. Mean turn depends on resampling density, so the curvature proxy and fixed-spacing controls also matter.

For loops, position and length use the **whole graph**, not a selected half-ring path. Analytic ring/branch truth is sampled approximately every 0.1 mm and recovered graph edges every ≤0.2 mm; coverage is the fraction of truth samples within one maximum voxel spacing of the recovered network. This tolerance-based coverage is explicitly different from the v1 projected-extent fraction and cannot certify topology. Network cycle rank is `E−N+components`, with endpoint/logical-junction counts checked separately. Logical junction centers are assigned one-to-one to analytic attachment points for mean/max junction errors; unmatched and spurious junction counts remain explicit. An unbranched ring has no applicable junction-position error. Loop tangent/smoothness metrics refer to the diagnostic principal path, with junction neighborhoods excluded from tangent error. A tree with high proximity coverage still fails cycle preservation.

Rows retain extraction and QC as separate outcomes. Error averages use successful geometries with equal geometry weight. Failure-inclusive coverage assigns failed extraction zero coverage, never zero positional error. Common-Lee-success results have identical geometry support, though raw/refined methods generally sample different positions. The all-64 aggregates are also saved for completeness but combine principal-path and network estimands; they should not replace the separated tables below. No confidence intervals or clinical superiority claims are inferred.

## Main results on the 60 v1 geometries

| Metric | Lee | Raw geodesic | Refined geodesic | Hybrid |
|---|---:|---:|---:|---:|
| Extraction success | 48/60 | 60/60 | 60/60 | 60/60 |
| QC pass rate | 80% | 100% | 100% | 100% |
| Mean positional error, mm | 0.127582 | 0.184467 | 0.049762 | 0.112664 |
| Mean tangent error, ° | 0.906520 | 2.448944 | 1.558906 | 1.015778 |
| Mean coverage on successful runs | 89.71% | 98.26% | 99.51% | 91.73% |
| Failure-inclusive coverage | 71.77% | 98.26% | 99.51% | 91.73% |
| Mean absolute length error, mm | 2.262111 | 1.273971 | 0.547076 | 1.906035 |
| Mean endpoint error, mm | 1.567705 | 0.642736 | 0.534740 | 1.354351 |
| Interior positional error, mm | 0.129044 | 0.150842 | 0.015816 | 0.106797 |
| Mean local turn, ° | 16.799747 | 14.538102 | 0.990398 | 13.555892 |
| Mean curvature proxy, 1/mm | 0.401277 | 0.352741 | 0.034983 | 0.325108 |

Lee's error columns use 48 successes; the other columns use 60. Refinement improves positional error in 49 geometries, worsens it in four and ties in seven (1e-6 numeric tolerance). Tangent error improves in 44, worsens in seven and ties in nine. Absolute length error improves in 32 and worsens in 28, even though its mean decreases substantially. Aggregate improvements are not universal corrections.

The hybrid selects **Lee 48 times and refined fallback 12 times (20%)**. All selected fallbacks pass final QC. On all 64 geometries including loops, Lee extraction/QC is 52/64 = 81.25%; each other method extracts and passes numerical QC on 64/64. Hybrid selects Lee 52 times and fallback 12 times (18.75%). This numerical QC success must be read alongside the loop failures of raw/refined fallback.

### Common-success comparison

On the same 48 Lee-success geometries:

| Metric | Lee and hybrid, exactly equal | Raw geodesic | Refined geodesic |
|---|---:|---:|---:|
| Position error, mm | 0.127582 | 0.143830 | 0.048955 |
| Tangent error, ° | 0.906520 | 2.553158 | 1.585431 |
| Coverage | 89.71% | 98.19% | 99.44% |
| Absolute length error, mm | 2.262111 | 1.517388 | 0.563413 |

The hybrid preserves Lee's successful geometry exactly, including its limitations. It does **not** improve coverage on retained Lee paths. Refining every path has lower average positional error but still worse tangents than Lee on this cohort; this is not evidence to switch the default.

### Known failures, phases, spacing and orientation

All ten half-voxel cylinders and the 1 mm/5 mm quarter-phase anisotropic cylinders trigger fallback and are recovered. On these 12 cases, raw→refined position error is **0.347017→0.052990 mm**, tangent error **2.032087→1.452806°**, and coverage **98.53%→99.79%**. Absolute length error increases **0.300302→0.481730 mm**. Hybrid results equal refined fallback for each case. Every failure has its own diagnostic figure.

On half-phase cylinders alone, position error decreases **0.361917→0.040046 mm**. The hybrid retains all eight Lee-success quarter-phase geometries and recovers the other two. Centered cases retain Lee.

On isotropic/anisotropic v1 cohorts respectively, refined positional errors are **0.065702/0.033822 mm** versus raw **0.193239/0.175696 mm**; refined tangent errors are **2.407649/0.710163°** versus raw **2.917961/1.979927°**. Hybrid succeeds on all 30 cases per spacing, compared with Lee's 25/30 and 23/30. Different in-plane as well as through-plane spacings and failure composition prevent interpreting this as a causal anisotropy estimate.

For directions (1,1,1)/(1,2,1)/(2,1,3), raw positional error means are **0.139299/0.212075/0.231241 mm**, refined **0.042638/0.025361/0.044489 mm**. Tangent means decrease **3.167595/2.612674/4.981788° → 1.607290/0.852113/1.675876°**. Hybrid retains Lee on every oblique case. Finite voxel-step mean turns exactly at 60° pass the inclusive QC threshold with roundoff tolerance.

Curved positional/tangent means improve **0.144535 mm / 2.889945° → 0.050222 mm / 1.991720°**. Taper values change **0.017858 mm / 0.818948° → 0.016932 mm / 0.671942°**. The hybrid retains Lee for all curves and tapers. Only one bend radius and two taper patterns are represented.

## Refinement sensitivity and endpoint experiment

| Configuration | Position, mm | Tangent, ° | Absolute length error, mm | Coverage |
|---|---:|---:|---:|---:|
| Raw | 0.184467 | 2.448944 | 1.273971 | 98.26% |
| Resample only, 0.5 mm | 0.171789 | 2.276585 | 0.689835 | 98.26% |
| Resample only, 1.0 mm | 0.172211 | 2.309302 | 0.464660 | 98.26% |
| Centered, 0.5 mm, no smoothing, extension | 0.050181 | 1.549883 | 0.575930 | 99.52% |
| Declared candidate: centered, 0.5 mm, sigma 1 mm, extension | 0.049762 | 1.558906 | 0.547076 | 99.51% |
| Centered, 0.5 mm, sigma 2 mm, extension | 0.050675 | 1.551239 | 0.495142 | 99.48% |
| Candidate without extension | 0.041851 | 1.317210 | 0.594333 | 98.25% |

Resampling alone reduces stair-step polyline length inflation, but also replaces small corners with chords; it is not merely relabeling coordinates. Coarser sampling can reduce reported length error without improving endpoint localization or physical truth everywhere.

The six centered, extension-enabled configurations are stable in extraction (all 60/60) and coverage (99.45–99.52%). Their position errors range **0.049762–0.059636 mm** and tangent errors **1.549883–1.762326°**. The 0.5 mm variants are more stable geometrically in this matrix than 1 mm variants; smoothing sigma 0–2 mm has relatively small aggregate effects once recentering is enabled. No-centering variants retain much larger positional error (0.165–0.181 mm). The evidence favors recentering, but does not establish a unique optimal smoothing scale.

Endpoint extension reduces mean endpoint error **0.585821→0.534740 mm**, increases coverage **98.25%→99.51%**, and slightly improves mean length error **0.594333→0.547076 mm**. It improves endpoint error in 34 geometries, worsens it in 24, and ties in two. Length error improves/worsens in 29/29 cases (two ties). Interior error changes only **0.015782→0.015816 mm**, but whole-path position and tangent errors worsen. Tangent error increases **1.317210→1.558906°**, worsening in 46 cases. Therefore endpoint extension is **not** a universal improvement and should remain a separately configurable experiment.

The worst refined tangent example is the centered 5 mm isotropic cylinder: mean tangent error **9.8233°**, positional error **0.251315 mm**, endpoint error **2.443668 mm**. Large off-axis raw cap regions can yield a biased local tangent; following that tangent toward the interpolated boundary does not identify the analytic cap center. Passing QC and staying within the mask do not resolve this error. Do not hide this exception behind good interior centering.

## Branches and loops

Both Y geometries retain three endpoints, one junction and zero cycles with every method. Refinement preserves raw junction coordinates exactly: errors remain 0.707107 mm isotropic and 0 mm anisotropic. Lee/hybrid junction errors remain 0 and 0.8 mm. Refined principal positional/tangent means improve **0.096139 mm / 2.773172° → 0.051892 mm / 1.646790°**, with coverage **98.99%→99.88%**. The stricter v1 branch criterion passes **3/6 Lee, 5/6 raw, 5/6 refined and 3/6 hybrid**. Hybrid inherits Lee's Y endpoints because QC passes. Geometric refinement does not automatically improve every branch's endpoint criterion.

All four loop networks retain exactly one cycle with Lee and hybrid. Raw and refined geodesic outputs have **zero cycles**, four endpoints and two junctions, even on the unbranched ring. They remain explicitly wrong networks despite high proximity coverage:

| Geometry | Lee/hybrid coverage | Raw coverage | Refined coverage | Lee/hybrid cycles | Fallback cycles |
|---|---:|---:|---:|---:|---:|
| Isotropic ring | 100% | 88.27% | 93.64% | 1 | 0 |
| Isotropic ring + branches | 100% | 88.27% | 92.33% | 1 | 0 |
| Anisotropic ring | 100% | 90.85% | 96.62% | 1 | 0 |
| Anisotropic ring + branches | 100% | 91.73% | 95.34% | 1 | 0 |

The hybrid succeeds on these loops **because Lee passes and is retained**, not because fallback is loop-aware. A future failed Lee loop could still receive a tree with passing numerical QC. The exposed capability flags and `production_eligible=False` prevent any claim of readiness, but they are not an implemented mask-derived loop detector. Full loop-aware fallback and topology QC remain necessary before Circle-of-Willis use.

## Direct answers and recommendation

Refinement reduces geodesic positional and tangent error on average, preserves 100% extraction on v1, and preserves high coverage. Physical resampling reduces aggregate voxel-path length inflation. Recentered 0.5 mm configurations are the strongest stable region in this small sweep; this is descriptive sensitivity evidence, not external parameter validation. Endpoint extension has mixed benefits and clearly worsens tangent error in many cases.

The hybrid chooses Lee on all QC-passing recoveries, falls back on all 12 known failures, improves failure-inclusive coverage, and preserves **exactly** Lee's position/tangent results on the common-success cohort. It retains Y topology and Lee's four tested loops. It does not make tree fallback safe for Circle-of-Willis networks, and there is insufficient evidence for a production change.

Continue experimental fallback work with independent/held-out phases, orientations and acquisition perturbations, improved endpoint centering, conservative topology detection, loop-aware extraction, wider Y angle/diameter/branch-length sweeps and downstream caliber comparisons. No patient validation, disease inference, statistical confidence bound or automatic promotion is provided.

## Tests and artifacts

**107 tests pass**: 89 existing tests plus 18 new tests, with 55 existing upstream scikit-image/NumPy deprecation warnings. Tests cover mask-only selection, correct fallback on failed/short Lee, numeric QC boundaries, physical resampling, straight/curved preservation, dense lumen constraints, anisotropic fallback, endpoint controls, fixed Y junctions, loop limitation, determinism, finite/unique samples, exact v1 mask/ID reuse and output protection. Existing regression tolerances are untouched. Declared regression bounds use fractions of voxel spacing and observed physical behavior: half-phase 1 mm cylinders must remain below 0.05 mm positional error and 1° tangent error; the 3 mm curve must stay below a quarter-voxel mean positional error while retaining its analytic-scale bend.

Under `validation/outputs/centerline_v2/`:

- `results.csv`: all 256 headline rows, including extraction/QC, smoothness, endpoints, refinement and hybrid decisions.
- `matched_comparison.csv`: 192 Lee-versus-comparator pairs; missing comparisons stay blank.
- `profiles.csv`: actual recovered main/branch/network positions and truth errors; branch/main samples can overlap and must not be pooled as independent observations.
- `qc_results.csv`: every executed criterion with observed value, threshold, result and criticality.
- `refinement_sweep.csv`: all 960 runs, including every parameter variant and control.
- `branch_results.csv`: 24 Y branch rows; `loop_results.csv`: 16 whole-network diagnostics.
- `summary.json`, `geometries.json`: separate cohorts, parameters, source hashes, versions, analytic network metadata and affines.
- `masks/`: 64 binary NIfTI phantoms.
- `figures/`: 35 separate Matplotlib PNGs, including 12 known failures, two Y examples, four loop diagnostics, required summary plots and two parameter-sweep plots. Network examples are explicitly x-z projections; no projection is used for quantitative errors.
- `manifest.txt`: all 109 primary artifacts. Nested isolated regressions have separate manifests and `checks.json`.

Created source files: `neurovasc/geometry/{centerline_qc,centerline_refinement,hybrid_centerline}.py`; `validation/synthetic_centerline/{experiments_v2,comparison_v2,figures_v2}.py`; `scripts/run_centerline_validation_v2.py`; `tests/test_centerline_validation_v2.py`; and this document. Existing files and production behavior remain unchanged. Nothing was pushed.
