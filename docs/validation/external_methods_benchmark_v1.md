# External Methods Benchmark v1: NeuroVasc and actual VMTK

## Decision and scope

**Actual VMTK executed successfully in a separate environment. Keep all NeuroVasc production defaults unchanged.** This is synthetic geometric validation against analytic phantom truth. VMTK is an established external comparator, not the truth or a clinical reference standard. No patient data, disease inference, clinical acceptance threshold or deployment claim is involved.

The benchmark includes **68 existing geometries**: 58 direct single-tube comparisons, two partial Y-bifurcation comparisons, and eight cyclic method-assumption mismatches. Of 120 attempted VMTK jobs, **114 produced valid paths**: 58/60 oracle-seeded and 56/60 automatically seeded. The eight loops generate 16 explicitly unattempted/mismatch rows, excluded from success denominators. Every actual job was repeated; all **120/120** reproduced status, array shapes/identities, coordinates, radii and branch outputs within absolute tolerance **1e-10** without registration or reordering.

The main finding is a tradeoff, not a winner declaration. VMTK paths have low directed positional error, but substantial truncation under the declared closed-surface/vertex-seed/no-endpoint-append setup. Refined geodesic has better coverage and length accuracy in this matrix. NeuroVasc CE has lower diameter error than MIS on common valid support here, but these methods have different estimands. Actual SlicerVMTK CE was **not executed**; those values remain blank.

## Literature and implementation basis

The existing [bibliography](../literature/references.bib) already contains all required references; it was not modified. The following sources motivate the comparator, not the numerical conclusions of this experiment:

- **`vmtkCenterlines`**: [official centerline tutorial](https://vmtk.github.io/tutorials/Centerlines.html). VMTK traces paths through a surface-derived Voronoi representation and returns maximal-inscribed-sphere radii. The documented seed-to-pole behavior explains why paths need not reach surface seeds when endpoint appending is disabled.
- **`antiga2004decomposition`**: Antiga and Steinman, *Robust and objective decomposition and mapping of bifurcating vessels*, IEEE TMI 23(6):704–713, 2004, [DOI 10.1109/TMI.2004.826946](https://pubmed.ncbi.nlm.nih.gov/15191145/). This is the branch-decomposition methodology, not evidence that any arbitrary voxel mask has correct anatomy.
- **`piccinelli2009geometry`**: Piccinelli et al., *A framework for geometric analysis of vascular structures: application to cerebral aneurysms*, IEEE TMI 28(8):1141–1155, 2009, [DOI 10.1109/TMI.2009.2021652](https://doi.org/10.1109/TMI.2009.2021652). The [official geometric-analysis tutorial](https://vmtk.github.io/tutorials/GeometricAnalysis.html), **`vmtkGeometry`**, explicitly uses tree-like branching assumptions and describes a bifurcation origin weighted by MIS sphere surface information.
- **`slicerExtract`**: [SlicerVMTK Extract Centerline](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/ExtractCenterline.md) distinguishes approximate network extraction from seeded accurate tree extraction. Its documented circular-network workflow interrupts the circle with nearby outlet points, leaving a gap. That is not whole-network cycle preservation.
- **`slicerCrossSection`**: [SlicerVMTK Cross-section analysis](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/CrossSectionAnalysis.md) reports MIS diameter separately from surface cross-sectional area and area-derived circular-equivalent diameter, with spatial coordinates for the samples.

These official pages are mutable and were checked on 2026-09-28. Actual executed provenance is pinned by the installed wheel hashes, package versions and compiled-module/wrapper hashes in the metadata. Documentation on a moving branch is not represented as a release commit. The report draws conclusions only from the declared implementation and experiments below.

## Installation, runtime and reproducibility

| Environment | Python | VTK | Other relevant versions |
|---|---|---|---|
| Existing NeuroVasc `.venv` | 3.14.6 | 9.7.0 | Full dependency versions in `environment.json` |
| `external/vmtk/.venv` | 3.14.6 | 9.6.2 | VMTK 1.5.1, NumPy 2.5.3 |
| Slicer / SlicerVMTK CE | unavailable | unavailable | No version or results invented |

The [official installation route](https://vmtk.github.io/download/) is now PyPI. A compatible macOS arm64 wheel was installed into the separate venv, leaving the production environment intact. Conda was available but unnecessary. `external/vmtk/install-report.json` retains actual wheel URLs/hashes; the platform-specific `requirements.lock` pins all 13 installed wheels. No package was installed into NeuroVasc's existing `.venv`.

The wheel's umbrella `vtkvmtk` import fails on a segmentation-library dependency, `libITKLabelMap-5.4.1.dylib`. The direct `vtkvmtkComputationalGeometryPython` module imports and executes. The adapter calls its real compiled `vtkvmtkPolyDataCenterlines`, branch extractor, branch geometry and bifurcation reference-system filters. The loader error is preserved; no binary patch, library alias, substitute algorithm, or claim that the general PypeS CLI works is made.

From the repository root, the existing environments can rerun the full benchmark with:

```bash
.venv/bin/python scripts/run_external_methods_benchmark.py
.venv/bin/python scripts/verify_external_methods_repeatability.py
.venv/bin/python -m pytest -q
```

Or use the explicit interchange stages:

```bash
.venv/bin/python scripts/run_external_methods_benchmark.py --stage export
.venv/bin/python scripts/run_external_methods_benchmark.py --stage run
.venv/bin/python scripts/run_external_methods_benchmark.py --stage report
```

The first stage verifies archived masks/affines, evaluates unchanged NeuroVasc methods and exports inputs. The second launches only VMTK in its isolated interpreter. The third validates request/surface hashes before importing external measurements. `--reuse` additionally requires matching external-module and adapter-source hashes. Missing runtime or results remain unavailable; they are not fabricated or silently replaced with local estimates. Output paths are restricted to `validation/outputs/external_methods_v1/` and descendants. See the [runtime README](../../external/vmtk/README.md) for initial installation, wrapper command, optional integration tests and manual Slicer instructions.

Repeatability uses fresh worker processes, fixed VTK RNG seed 0 and one VTK SMP thread, with four independent jobs scheduled concurrently. `repeatability.json` contains all 120 primary/repeat comparisons and hashes; raw repeats live under `repeatability/`. This establishes repeatability for this environment and input matrix, not across platforms/VTK releases. It does not claim byte-identical XML compression or rerun every older milestone.

## Experiment matrix and unchanged baselines

The 58 tubes reuse the caliber-v2/centerline-v1 matrix: 26 core geometries, 20 extra grid phases and 12 extra oblique orientations. Straight diameters are 1–5 mm; oblique and curved diameters are 2, 3 and 4 mm; tapers are 4→2 and 3→1.5 mm. Length is 30 mm. Curvature is the existing 24 mm-radius arc. Phase shifts are 0, 0.25 and 0.5 voxels; oblique directions are normalized (1,1,1), (1,2,1), (2,1,3). The matrix is not fully crossed over phase and orientation.

Both `(0.5,0.5,0.5)` and `(0.46875,0.46875,0.8)` mm grids are retained. Two existing centered Y geometries have a 22 mm parent, 18 mm daughters, 4/3 mm diameters and ±35° daughter directions. Eight centered cyclic masks reuse simple ring, ring with branches, figure-eight and communicating bridge at each spacing. All 68 voxel arrays are checked exactly against existing NIfTI archives, and affines checked for numeric agreement; their SHA-256 hashes are in `matrix.json`. No phantom definition was changed.

| Scope | Geometries | External rows | Interpretation |
|---|---:|---:|---|
| `SUPPORTED_DIRECT_COMPARISON` | 58 | 116 | Single-tube centerline/radius/caliber, with support and seed caveats |
| `PARTIAL_COMPARISON` | 2 | 4 | Y branches/junctions; branch boundary definitions differ |
| `METHOD_ASSUMPTION_MISMATCH` | 8 | 16 | Whole-cycle network comparison outside the declared tree workflow |

Lee, raw geodesic, refined geodesic and hybrid run unchanged on all 60 noncyclic geometries. The raw standalone topology extractor is added for the two Y cases. That gives **242 NeuroVasc centerline rows + 136 VMTK rows = 378 centerline rows**. Seed modes are related paired diagnostics, not independent statistical replicates. All four local methods succeeding on a case does not validate clinical performance.

EDT diameter and experimental NeuroVasc CE reuse the unchanged Lee/main-path pipeline, physical EDT and experimental cross-sectional implementation. CE uses its existing physical tangent window and polygon-area estimator. No Gaussian caliber smoothing, phase correction or estimator fallback is introduced. Shared Lee failures remain shared missing caliber results. The different evaluation station grid below means these aggregate values need not equal previous milestones' per-recovered-sample averages.

## Surface preparation and seeds

All methods start with the same binary segmentation. The external surface uses the existing `mask_to_mesh` 0.5 marching-cubes surface and full affine to world mm, followed by zero-tolerance VTK cleaning and triangulation. No smoothing, decimation, resampling of the mask, axis flip, end clipping or new capping is applied. The phantom's binary end caps remain; no open profiles are created. Per-surface bounds/point/triangle counts and all settings are saved.

**Oracle diagnostic:** analytic terminal coordinates determine seeds, explicitly `oracle_seed=true`. This isolates some endpoint-detection effects, but does not guarantee optimal VMTK seed poles. **Automatic benchmark:** endpoints come from the unchanged raw geodesic mask graph without analytic truth, sorted `(z,x,y)` with first source and remaining targets. This is a declared NeuroVasc-derived seeding strategy, not VMTK's or Slicer's automatic endpoint detector. Seeds in both modes snap to nearest surface vertices; requested/snapped coordinates and vertex IDs are recorded.

Centerline cost is `1/R`; endpoint appending and centerline resampling are off, Voronoi simplification is off, normal flip is off, and Delaunay tolerance is 0.001. The unused resampling step is recorded as 0.5 mm. Centerline coordinates, `MaximumInscribedSphereRadius`, Voronoi surface, branch tracts and reference systems are actual filter outputs.

**This specific capped, surface-vertex-seed setup has a major endpoint limitation.** Several surfaces produce very short paths or poorly covered interiors even with analytic seeds. No post-hoc endpoint extension or optimized open-profile preparation was selected to improve the score. Accordingly, the measured coverage/length ranking must not be generalized to a carefully prepared open-profile VMTK workflow. That is a clear next controlled experiment.

## Evaluation, common support and failure handling

Centerline position is directed recovered-point-to-analytic-centerline distance after 0.1 mm physical sampling, with duplicate sampled points removed. Path length sums unique physical graph edges; shared identical source-target segments are not counted twice. Endpoints receive one-to-one assignments, with missing/extra counts alongside error; this endpoint-error assignment has no distance gate. Tangents use the same 6 mm physical fit for every method and exclude a 3 mm Y-junction neighborhood. These are sampled geometric measures, not exact continuous integrals.

Coverage asks which of 601 analytic tube stations (401 per Y branch) are within one maximum voxel spacing of recovered points. A successful, short centerline can have low coverage. Failure-inclusive coverage assigns failed **attempted** runs zero coverage; unavailable and mismatch cases are reported separately, not treated as accurate zeros. Position and length means condition on successful extraction and therefore have different support sets.

Centerline matched comparisons additionally evaluate **truth-to-path distance on common covered analytic stations**. This is a different direction/metric from recovered-to-truth positional error. It exposes the support intersection and cannot erase the missing vessel outside that intersection. Unscorable pairs remain explicit.

For tube caliber, 105 fixed analytic stations span 2–28 mm at 0.25 mm increments. Recovered world coordinates are projected onto the analytic axis/arc; truth for tapers is evaluated at those physical positions. Interpolation requires monotone projected paths, proximity within one maximum spacing, and a bracket gap ≤2 maximum spacings. No extrapolation or NaN bridging is allowed. These engineering gates intentionally reject sparsely sampled/remote support and are not clinical tolerances. All original VMTK samples remain available even when a station is rejected. Caliber success is distinct from extraction success.

MAE, RMSE, signed bias and mean absolute percentage error use valid corresponding samples; radius is half diameter. Common-support comparisons use the intersection of valid stations, preserve coverage and failed pairs, and never compare unlabelled own-support means as if they were paired agreement. Aggregate MAE/bias are equal-run means; aggregate RMSE is the square root of mean run MSE. No confidence intervals are inferred from these deterministic phantoms.

## Centerline results: 58 direct geometries

Position and length columns below use each method's successful support. Coverage columns distinguish success-only from all attempted runs.

| Method | Success | Position error (mm) | Coverage, success (%) | Coverage, all (%) | Length error (mm) |
|---|---:|---:|---:|---:|---:|
| Lee | 46/58 | 0.123357 | 93.4312 | 74.1006 | 2.250271 |
| Raw geodesic | 58/58 | 0.173537 | 97.7767 | 97.7767 | 1.268899 |
| Refined geodesic | 58/58 | 0.043563 | 99.2828 | 99.2828 | 0.532111 |
| Hybrid | 58/58 | 0.107507 | 94.6239 | 94.6239 | 1.884366 |
| VMTK oracle | 56/58 | 0.038957 | 70.4154 | 67.9873 | 9.903927 |
| VMTK automatic | 54/58 | 0.059168 | 78.4618 | 73.0507 | 6.840763 |

On own recovered support, oracle VMTK has slightly lower position error than refined geodesic and both VMTK modes improve on Lee. This is not enough to rank entire recovered vessels: refined geodesic covers much more truth with lower length error. VMTK endpoint errors are 5.0440 mm oracle and 3.7566 mm automatic versus 1.5687 mm Lee and 0.5363 mm refined geodesic.

| Common covered truth support | Scorable pairs / 58 | Mean common coverage among those pairs (%) | First truth-to-path error (mm) | VMTK truth-to-path error (mm) |
|---|---:|---:|---:|---:|
| Lee / VMTK oracle | 45 | 72.6271 | 0.142567 | 0.074381 |
| Lee / VMTK automatic | 45 | 76.7018 | 0.139348 | 0.074224 |
| Refined / VMTK oracle | 56 | 70.3352 | 0.038216 | 0.090368 |
| Refined / VMTK automatic | 54 | 78.4618 | 0.033646 | 0.081571 |

The common-support directional comparison also favors VMTK over Lee but refined geodesic over VMTK. These numbers should not be substituted into the preceding table: direction, conditioning and physical support differ explicitly.

## MIS radius and diameter

The surface MIS and voxel EDT measure inscribed-sphere-related quantities on different geometric representations. Area-equivalent CE is another estimand, especially at tapers/noncircular sections. Synthetic circular truth enables a useful comparison but does not make the estimands identical on arbitrary vessels.

| Diameter comparison on common stations | Scorable pairs / 58 | Common fraction over all 58 (%) | First MAE (mm) | Second MAE (mm) |
|---|---:|---:|---:|---:|
| EDT / NeuroVasc CE | 46 | 78.2102 | 0.206308 | 0.077008 |
| EDT / VMTK oracle MIS | 41 | 57.3563 | 0.227378 | 0.260642 |
| EDT / VMTK automatic MIS | 43 | 60.7553 | 0.214121 | 0.260886 |
| NeuroVasc CE / oracle MIS | 41 | 57.3563 | 0.070754 | 0.260642 |
| NeuroVasc CE / automatic MIS | 43 | 60.7553 | 0.072154 | 0.260886 |
| Oracle MIS / automatic MIS | 49 | 67.1429 | 0.244614 | 0.245245 |

Thus matched **radius** MAE is 0.113689 EDT versus 0.130321 oracle MIS, and 0.107061 EDT versus 0.130443 automatic MIS (exactly half the corresponding diameter values). In this protocol, MIS does not improve average matched error over EDT; CE has lower matched diameter MAE than MIS. This is not evidence for a universal estimator ordering.

Own-support radius summaries, retained to expose bias and RMSE, are:

| Radius method | Scorable / 58 | MAE (mm) | Aggregate RMSE (mm) | Bias (mm) | Absolute percentage error |
|---|---:|---:|---:|---:|---:|
| EDT | 46 | 0.103154 | 0.137948 | -0.018140 | 8.2260% |
| Oracle MIS | 50 | 0.122396 | 0.133746 | -0.100461 | 8.8247% |
| Automatic MIS | 52 | 0.121749 | 0.135005 | -0.103388 | 8.5982% |

Own-support NeuroVasc CE diameter MAE/RMSE are 0.077008/0.099550 mm across 46 runs. The differing scorable counts are important: six successfully extracted oracle paths and two successfully extracted automatic paths have no valid gated caliber stations. They are recorded as `no_valid_gated_analytic_station_support`, not assigned zero error.

### The 1 mm vessel

Each table cell is **mean diameter / MAE / absolute percentage error**, on valid stations. `—` means no valid caliber support; it can be an extraction failure or a successfully extracted path with no accepted analytic stations, as distinguished in the CSVs.

| Grid / phase | EDT | NeuroVasc CE | Oracle MIS | Automatic MIS |
|---|---|---|---|---|
| Isotropic / centered | 1.414214 / 0.414214 / 41.4214% | 1.196827 / 0.196827 / 19.6827% | 1.224745 / 0.224745 / 22.4745% | same as oracle |
| Anisotropic / centered | 1.325825 / 0.325825 / 32.5825% | 1.122025 / 0.122025 / 12.2025% | 1.318572 / 0.318572 / 31.8572% | same as oracle |
| Isotropic / quarter | 1.000000 / 0 / 0% | 0.892062 / 0.107938 / 10.7938% | — | — |
| Isotropic / half | — | — | 1.224745 / 0.224745 / 22.4745% | — |
| Anisotropic / quarter | — | — | — | — |
| Anisotropic / half | — | — | — | — |

Centered EDT/CE have all 105 stations; MIS has only **30/105 isotropic and 22/105 anisotropic**, and half-phase isotropic oracle MIS has only **9/105**. Missing stations are not evidence of accuracy. The centered MIS estimates reduce EDT's error appreciably only on the isotropic grid here; both retain large relative error and CE is closer on the centered cases. Quarter-phase isotropic EDT's zero error is one favorable voxel arrangement, not general submillimeter reliability. Phase changes alter both accuracy and extraction/support availability. External Slicer CE remains unavailable for all six cases.

### Spacing and orientation

These are descriptive own-support diameter MAEs, with successful caliber counts. Spacing changes both in-plane and through-plane dimensions, so it does not isolate a causal anisotropy effect.

| Grid | EDT | NeuroVasc CE | Oracle MIS | Automatic MIS |
|---|---:|---:|---:|---:|
| Isotropic | 0.186226 (24) | 0.078730 (24) | 0.252925 (23) | 0.249498 (25) |
| Anisotropic | 0.228217 (22) | 0.075129 (22) | 0.237863 (27) | 0.237942 (27) |

For centerline position, isotropic→anisotropic errors are Lee 0.104975→0.143411 mm; refined geodesic 0.056558→0.030568; oracle VMTK 0.031754→0.045664; automatic VMTK 0.070811→0.047524. Successful-support VMTK coverage also changes substantially: oracle 58.2363→81.7545%, automatic 67.5356→89.3881%. Successful sets differ. These figures do not show invariance to spacing or establish an anisotropic correction.

Across oblique 2/3/4 mm diameters and both spacings:

| Direction | Lee position (mm) | Refined position (mm) | Oracle VMTK position (mm) | Automatic VMTK position (mm) | EDT / CE / oracle MIS / automatic MIS diameter MAE (mm) |
|---|---:|---:|---:|---:|---|
| (1,1,1) | 0.104215 | 0.036066 | 0.029426 | 0.049384 | 0.222730 / 0.089618 / 0.270589 / 0.272226 |
| (1,2,1) | 0.234639 | 0.021882 | 0.088746 | 0.098525 | 0.280506 / 0.054337 / 0.346007 / 0.346740 |
| (2,1,3) | 0.187283 | 0.038029 | 0.055069 | 0.080627 | 0.196186 / 0.044214 / 0.266854 / 0.258811 |

All six centerline cases succeed per direction/method. MIS caliber has only five scorable (1,1,1) cases; other entries have six. VMTK retains orientation-dependent position and MIS errors, especially for (1,2,1). Its oracle coverage for these directions is 79.70%, 78.51%, 85.08%; automatic is 80.45%, 79.03%, 80.45%. Small direction sets and different valid support prevent a broad rotational-robustness claim.

## Y bifurcation: actual branches and reference origins

The local logical-graph matching uses existing analytic branches and the existing 2.5 mm physical gate. For VMTK, **actual nonblanked GroupIds** define the three branches. Group terminal positions are assigned one-to-one to analytic terminals under the same gate; branch attachment is read through neighboring blanked junction groups in each ordered CenterlineIds/TractIds sequence. The reference origin must match the analytic junction under the gate. Reducing nearly duplicated source-target paths is not treated as VMTK's branch definition: those graph counts are explicitly secondary `duplicate_path_embedding_*` diagnostics.

Every actual VMTK Y run returns three matched branch groups attached to one blanked junction group and one bifurcation reference system: **4/4 full Y connectivity matches**, two modes × two spacings. This verifies these two Y configurations only. VMTK's sphere-weighted reference origin is not identical in definition to a local graph junction, so localization is a partial comparison.

| Method | Grid | Branches | Branch F1 | Junction error (mm) | Coverage (%) | Mean matched branch length error (mm) | Mean absolute branch mean-radius error (mm) |
|---|---|---:|---:|---:|---:|---:|---:|
| Lee | iso | 3 | 1.000 | 0.000 | 95.345 | 0.710 | 0.049 |
| Lee | aniso | 3 | 0.667 | 0.800 | 95.096 | 0.096 | 0.012 |
| Refined geodesic | iso | 3 | 1.000 | 0.707 | 99.584 | 0.748 | 0.077 |
| Refined geodesic | aniso | 3 | 1.000 | 0.000 | 100.000 | 0.195 | 0.083 |
| Raw topology | iso | 3 | 0.667 | 0.000 | 91.106 | 3.000 | 0.071 |
| Raw topology | aniso | 3 | 1.000 | 0.000 | 96.426 | 4.608 | 0.109 |
| VMTK oracle | iso | 3 | 1.000 | 0.832 | 97.589 | 1.322 | 0.200 |
| VMTK oracle | aniso | 3 | 1.000 | 0.469 | 97.506 | 2.360 | 0.094 |
| VMTK automatic | iso | 3 | 1.000 | 0.882 | 95.511 | 1.478 | 0.201 |
| VMTK automatic | aniso | 3 | 1.000 | 0.469 | 97.672 | 2.304 | 0.094 |

A count of three does not guarantee F1=1: endpoint shortening/matching can miss a branch. The local radius values are EDT interpolated at recovered branch points; VMTK uses the actual MIS array. Both exclude points within 3 mm of the analytic junction and 2 mm of the caps, but remain **descriptive own-branch-support** measurements, not matched cross-method radius agreement. VMTK blanking cuts branches away from the mathematical junction, so raw group lengths and graph-to-junction lengths do not have identical boundaries. Errors of unmatched branches are blank, not zero. Detailed branch assignments, lengths, radii and raw geometry arrays are retained in `branch_results.csv`.

For context, oracle VMTK mean radii for parent/left/right are **1.8345/1.2822/1.2822 mm isotropic**, **1.9315/1.3924/1.3934 mm anisotropic**, against truth **2/1.5/1.5 mm**. VMTK total unique-path length errors are 0.4266/1.5974 mm oracle and 0.4987/1.4278 mm automatic; these differ from branch-length errors because junction regions/cut boundaries differ. Y overlays show all methods against analytic truth in x-z projection; quantitative calculations use all three coordinates.

### Tortuosity definitions

The installed `vmtkbranchgeometry.py` output description specifies tortuosity as **`L/chord - 1`**; the [upstream C++ implementation](https://github.com/vmtk/vmtk/blob/master/vtkVmtk/ComputationalGeometry/vtkvmtkCenterlineBranchGeometry.cxx) agrees. Its input description is less precise, so the output definition and actual straight-cylinder integration test are used. The executed binary/wrapper hashes identify the release evidence; a mutable source URL is not asserted to be its exact commit. NeuroVasc's Distance Metric and the Slicer Extract Centerline documentation use **`L/chord`**.

Raw VMTK values and explicitly converted `+1` values are both retained. For isotropic oracle Y branches, converted parent/left/right values are **1.0189/1.0388/1.0472**, versus Lee **1.0000/1.0562/1.0562** and refined geodesic **1.0286/1.0050/1.0166**. Ideal analytic branches are straight (Distance Metric 1). Anisotropic oracle VMTK values are **1.0004/1.0110/1.0111**. Different blanked-group boundaries and group weighting prevent claiming point-for-point metric equivalence simply by adding one. Actual Slicer tortuosity was not executed.

## Seeds, loops and failures

Automatic versus oracle VMTK changes direct extraction success **54 versus 56**, position error **0.0592 versus 0.0390 mm**, successful coverage **78.46 versus 70.42%**, and length error **6.84 versus 9.90 mm**. Automatic seeds improve coverage here while worsening average directed position. Neither is uniformly preferable, and oracle anatomical endpoints are not necessarily optimal surface poles. On 49 common caliber pairs, the MIS MAE difference is only **+0.000631 mm** automatic minus oracle, despite different overall support. Do not mix modes in a single headline result.

All six external extraction failures contain nonempty but effectively zero-length lines, rejected at ≤1e-7 mm:

- Oracle: isotropic centered 2 mm curved tube; isotropic half-phase 3 mm straight tube.
- Automatic: isotropic 1 mm quarter/half-phase tubes and anisotropic 1 mm quarter/half-phase tubes.

The parser checks finite coordinates/radii and valid polylines instead of accepting nonempty output as successful extraction. Centerline logs and JSON reasons remain available. Caliber failures additionally include gated-support absence; 12 Lee extraction failures cause shared EDT/CE absence. These are separate from unavailable Slicer CE and assumption mismatches in `failures.csv`.

The loop investigation is **source-based scope analysis with exported test surfaces, not an executed whole-cycle recovery experiment**. The official branch/geometric workflow assumes trees; Slicer's documented cycle handling introduces a cut. Ring and figure-eight have no anatomical terminals, and communicating networks have multiple routes not certified by a seeded shortest-path tree. We therefore do not invent seeds/cuts and score the resulting tree against closed-cycle truth. All eight masks retain scope labels, seed-failure reasons and unattempted status. This does not claim every VMTK algorithm is incapable of cycle-related processing; approximate network extraction and a separately declared cut-network experiment remain untested here.

## SlicerVMTK CE and fair-comparison limits

No verified Slicer application/extension runtime was found; no headless CE calculation was claimed. `cross_sections.py` defines required coordinate/area/diameter/version fields and reports unavailable. The runtime README provides a manual procedure using the exported unmodified surfaces and real VMTK centerline VTPs, explicit RAS/LPS handling, raw CSV retention, source hashes and version recording. The future comparison must validate actual Slicer output and common physical support. The current report deliberately contains no Slicer CE measurement values.

Fair direct questions are reproducibility, extraction/support coverage, positional/length errors under the declared preparation, and raw/MIS/CE error against analytic truth with support exposed. Y branch identity, branch length, radius and reference-origin questions are partial because decomposition boundaries/weighting differ. A whole-cycle comparison against this tree workflow is not methodologically fair. None of these results imply clinical stenosis accuracy or transfer to real CoW segmentation.

## Tests, artifacts and preservation

**236 tests passed**, comprising all **192 existing tests**, **42 new infrastructure unit tests**, and **2 real VMTK subprocess integration tests**. There are 58 upstream scikit-image/NumPy deprecation warnings: the original 55 plus three surface-export calls in new tests. Warnings were not suppressed, and existing tolerances were not weakened. Unit tests do not require VMTK; the separate integration file explicitly skips if its external probe is unavailable. Coverage includes environment detection, no fabricated fallback, VTP interchange, affine coordinates, seeds, parser failures, actual Y groups, oracle/scope flags, common-support missingness, failure accounting, output protection and repeatability comparison behavior.

`preservation.json` audits **5,158 pre-existing files**: **4,921 match their baseline SHA-256 hashes**, and **237** files that timed out during the initial baseline read retain unchanged size/mtime. Those 237 are not claimed byte-verified. No compared metadata changed. The audit includes production code, Case 001/data, frontend/backend, documentation and all six previous validation archives and regression copies. Dependency directories, caches, build output and Finder metadata are excluded. Older archive runners were not executed or rewritten; their existing regression tests passed. Only new files were added for this milestone; nothing was pushed.

Artifacts in `validation/outputs/external_methods_v1/`:

- Required `results.csv`, `centerline_results.csv`, `caliber_results.csv`, `branch_results.csv`, `matched_comparison.csv`, `failures.csv`; additional `caliber_profiles.csv` stores all analytic-station samples including missing estimates.
- `summary.json`: scope/counts, successful/failure-inclusive aggregates, common-support results, 1 mm, spacing/orientation and Y details, versions/source hashes and recommendation.
- `environment.json`, `surface_parameters.json`, `matrix.json`, `requests.json`: actual runtime/install evidence, masks/truth, seeds and preparation.
- `inputs/`: 68 exported surfaces and both seed-mode request files. Archived NIfTI masks are reused, not duplicated or overwritten.
- `vmtk_outputs/`: real centerlines, Voronoi, branch/geometry/reference VTPs, full numerical JSON and process logs; mismatch/failure records retain their explicit statuses.
- `figures/`: **14 separate Matplotlib PNGs** covering the requested A–M topics, including two Y overlays. External CE has an unavailable notice, not a fabricated curve. Scatter/own-support figures state their support limits; common-support tables carry the paired inference.
- `repeatability.json`, `repeatability/`: fresh execution and comparisons for all 120 attempted jobs; `pytest.log`, `preservation.json` provide test and preservation evidence.
- `manifest.txt`: generated evidence inventory. It excludes caches, hidden metadata and early `smoke/`/`development/` probes; those local diagnostics are not part of the declared experiment. `validation_checks.json` records table counts, source-hash checks and explicit missing-result checks.

New source/setup/test/documentation files:

```text
external/vmtk/__init__.py
external/vmtk/.gitignore
external/vmtk/README.md
external/vmtk/requirements.txt
external/vmtk/requirements.lock
external/vmtk/install-report.json
external/vmtk/run_vmtk_benchmark.sh
external/vmtk/adapter.py
external/vmtk/surface_preparation.py
external/vmtk/centerlines.py
external/vmtk/branch_metrics.py
external/vmtk/cross_sections.py
validation/external_methods/__init__.py
validation/external_methods/evaluate.py
validation/external_methods/matching.py
validation/external_methods/report.py
validation/external_methods/figures.py
scripts/run_external_methods_benchmark.py
scripts/verify_external_methods_repeatability.py
tests/test_external_methods.py
tests/test_external_methods_vmtk_integration.py
docs/validation/external_methods_benchmark_v1.md
```

The ignored local `external/vmtk/.venv/` is the only new dependency environment. Exact generated artifact names are in the manifest.

## Limitations and next decision

The matrix has ideal binary phantoms, few phases/directions, two Y geometries, no acquisition blur, partial volume, noise, segmentation perturbation, patient anatomy or independent holdout. Surface voxelization, MIS versus area estimands, interpolation gates, seed-pole placement and different valid support all influence results. Analytic truth is used only for oracle seeds and evaluation; automatic endpoint selection still depends on NeuroVasc. Short path acceptance is explicitly separated from coverage and caliber support. Loop recovery and actual Slicer CE remain unexecuted, not negative accuracy evidence.

The benchmark adds an independently implemented external method and reproducible evidence, but does not justify a NeuroVasc default change. Keep production Lee and EDT, all experimental statuses, frontend/API and Case 001 unchanged. Next, compare declared open-profile/capping/endpoint-append variants without changing the primary archived protocol, inspect sparse/degenerate VMTK paths, validate true Slicer CE on common coordinates, and expand branching/phase/segmentation experiments. Parameter tuning and evaluation should then be separated with a holdout matrix before making a production recommendation.
