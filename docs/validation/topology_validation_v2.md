# Topology Validation v2: junction-aware geometry with strict cycle preservation

## Scope and decision

This is **synthetic geometric/topological validation, not clinical validation**. Production Lee, the current hybrid's experimental deployment status, EDT caliber, the experimental cross-sectional estimator, frontend, FastAPI and Case 001 are unchanged. All earlier validation archives are preserved.

V2 retains **76/76 correct cycle ranks and 58/76 strict connectivity successes** on the exact v1 cohort. It improves average physical length and position error and modestly improves junction localization, but **does not fix figure-eight structure or increase branch/adjacency F1**. Coverage regresses from **97.6076% to 90.8395%**. The additional junction cohort exposes further inherited structural failures. **Keep the topology method standalone and experimental; do more validation before integrating a hybrid fallback or changing production.**

## Reproduce

From the repository root:

```bash
source .venv/bin/activate
python3 scripts/run_topology_validation_v2.py
python3 -m pytest -q
# Rerun all five earlier milestones without overwriting their archives:
python3 scripts/verify_topology_v2_regressions.py
```

The v2 runner also works by absolute path and permits `--output` only at `validation/outputs/topology_v2/` or a descendant. Artifacts are deterministic for this source/dependency environment, with no random seed or timestamps. The summary records parameters, dependency versions and source SHA-256 hashes. Reruns replace this milestone's artifacts only.

The verification helper executes these original commands in a temporary source mirror, using the same Python environment:

```bash
python3 scripts/run_topology_validation.py
python3 scripts/run_centerline_validation.py
python3 scripts/run_centerline_validation_v2.py
python3 scripts/run_caliber_validation.py
python3 scripts/run_caliber_validation_v2.py
```

The helper orders caliber before centerline because the latter reads caliber-v2 results. It copies verification artifacts to `topology_v2/regressions/`. All five commands complete successfully. All 16 caliber/centerline CSVs reproduce byte-for-byte. Topology-v1 metrics reproduce as well, with a documented existing qualification: opaque recovered branch/node IDs and row order can vary across Python processes; floating summation order can change a metric by approximately 1e-16. The helper compares topology tables after ignoring opaque recovered IDs, sorting membership sets/rows and rounding numeric fields to 1e-10. Truth assignments, statuses, counts and quantitative metrics must remain equal. It does not modify v1 or weaken existing test tolerances. The final run has 20/23 CSVs byte-identical, with the three opaque-ID detail tables semantically equal. Exact per-file byte and semantic comparisons are saved in `regressions/checks.json`. All 1,703 archived files are unchanged. The broader before/after audit found 1,849 existing files unchanged; only two Finder `.DS_Store` metadata files changed, with no source/data changes (`regressions/preservation.json`).

## Topology and geometry are separate

The unchanged v1 extractor supplies a `TopologyNetwork`. Its logical MultiGraph holds node identities, branch identities, endpoint pairs, adjacency and fundamental cycle memberships. Each branch separately holds physical polyline coordinates, length and refined tangents. V2 never runs a topology reduction on smoothed coordinates to decide the accepted logical graph.

Geometry is embedded using explicit node/branch IDs. Interior samples receive branch-specific identities, including parallel routes and self-loop branches. Coincident or nearby coordinates do not silently add graph edges. Junction-attached paths end exactly at their shared logical node coordinates.

The default contract is **geometry-only**: no branch deletion, creation, junction merge or split is authorized. Consequently, it cannot repair incorrect combinatorial topology inherited from v1. This is a deliberate limitation, not a claim that component count and cycle rank alone certify the correct anatomy.

### Invariants and rejection

`topology_invariants.py` checks before/after:

- Raw and logical component counts; cycle rank; branch, endpoint and junction counts.
- Exact branch and node identity sets, node kinds, keyed logical edges, branch endpoint pairs and cycle memberships.
- Exact branch attachment to shared logical coordinates, finite geometry and agreement between embedded and logical cycle rank.
- Actual embedded edge/node identities and coordinates against a fresh embedding of the declared branch paths. This catches broken/reassigned branch edges even if global counts happen to agree.

A proposal that fails an invariant or geometric check is returned with `accepted=False`, an explicit reason, the proposal checks and the original network. The original is not mutated. The runner distinguishes `refined`, `baseline` and `original_retained_after_rejection`, and reports acceptance separately from extraction success. It never presents a retained original as successful refinement.

Geometric checks densely sample every branch at **minimum voxel spacing / 8**, requiring interpolated mask occupancy ≥0.5. They reject contact between distinct logical components at that sampling scale even if node identities remain distinct. These are finite-resolution checks, not a mathematical proof against all possible segment intersections. They do not infer new connectivity from geometry.

## Junction localization and disambiguation

For each existing logical junction, nearby foreground voxel centers within the configured physical radius are averaged with squared EDT-clearance weights and Gaussian distance weights. Movement is limited to the smaller of the neighborhood radius and maximum voxel spacing, with occupancy-constrained backtracking. All incident branch endpoints are reattached to that one coordinate. The default radius is **1 mm**.

The disambiguation stage inspects local branch incidence, logical degree, cycle memberships and bridge status. A short connection between two junctions that carries no cycle membership can be either a split digital junction or a real communicating connector. It is recorded as `ambiguous_split_or_true_connector`; a connection participating in a cycle is classified as `cycle_passage`. High-degree junctions and ordinary branch junctions are recorded separately. **Neither ambiguity nor spatial closeness triggers a merge.** Every decision is in `junction_decisions.csv`.

This conservative classification prevents a close double-junction pair from being deliberately merged, but leaves the figure-eight's split junction unresolved. A future declared junction-merge experiment needs additional validated evidence and its own identity mapping/invariant contract. V2 does not pretend that moving two distinct node coordinates is equivalent to fixing their adjacency.

## Spur classification

Every terminal branch is evaluated for four signals: length below two maximum voxel spacings, no cycle membership, endpoint inside the adjacent junction's EDT clearance region, and endpoint clearance below 75% of that region's clearance. At least three signals mark a candidate. Records include score, reasons, removal decision and the topology effect of hypothetical removal.

**Conservative mode marks only; it removes nothing.** Deleting even a cycle-free terminal branch would violate branch/endpoint identity preservation. Cycle-carrying communicating branches have no terminal endpoint and are protected from this classifier. Switching spur handling off changes diagnostic labels, not geometry.

Across the 106 full-refinement runs, 7 of 148 terminal branches are marked. Five are unmatched/spurious under analytic evaluation; **two are real matched short terminal branches** in anisotropic quarter/half-phase `near_junction_branch` cases. Removing all marked branches would therefore damage real anatomy even in this synthetic matrix. No marked or unmarked branch is deleted. This is evidence against automatic pruning using the current heuristic, not evidence that spur removal improves connectivity.

## Physical branch refinement and endpoints

The full stage uses the existing experimental physical chain-refinement routine without invoking its tree-only network method:

1. Attach the branch ends to refined shared nodes.
2. Resample at ≤0.5 mm physical arc-length intervals.
3. Fit physical tangents and recenter interior points using EDT-weighted foreground centroids in local tangent-normal slabs.
4. Apply a local quadratic fit in physical arc distance with a 1 mm smoothing scale.
5. Preserve both branch-end identities and coordinates; resample again and reconnect exactly.
6. Check dense containment, topology and component separation.

Movement proposals use the existing two-maximum-voxel displacement bound and occupancy backtracking. Closed branches retain their closing node identity; endpoint-constrained fits are not periodic splines, so a loop-anchor seam remains a limitation. Raw paths remain available for comparison. Junction-only refinement performs the shared-coordinate step without the branch recentering/smoothing stage.

Calculations run in an orthonormal voxel-axis frame measured in mm, then transform back to world coordinates. This does not resample the mask or change distances. Coordinates are rounded at 1e-12 mm in that frame to stabilize floating boundary/tie behavior under rigid affine relabeling. A strict affine-rotation test exposed the issue before this fix and now passes without relaxing its tolerance. Physical rotation of the vessel relative to the grid is a separate experiment.

Endpoints have four configurable experiments:

- `raw`: retain the current endpoint; default for the declared main candidate.
- `recentered`: transverse EDT-weighted recentering using a local physical tangent.
- `tangent`: propose an outward extension by two maximum voxel spacings, with occupancy backtracking and dense path checks.
- `boundary`: advance in 0.05 mm steps, stopping at the first occupancy boundary or the same extension limit.

These modes cannot restore an already deleted branch. They do not extend through a later foreground island. Tangent stability, endpoint error, coverage and branch/network lengths are recorded in the endpoint sweep. Extension is not selected separately for each phantom.

## Experiment matrix and parameters

The core cohort reuses **all 76 v1 geometry IDs, masks, dimensions and physical truth**. The runner checks voxel-array equality and affine agreement with every archived NIfTI before evaluation. No v1 file is rewritten.

Five new analytic networks are each evaluated at both existing spacings and all three phases, adding **30 geometries**:

| New geometry | Analytic definition | Components / endpoints / junctions / branches / rank |
|---|---|---|
| Four-way junction | Four 8 mm radial branches, diameter 2 mm | 1 / 4 / 1 / 4 / 0 |
| Double junction | Junctions at x=±1.5 mm, 3 mm connector; four outward branches to (±8,0,±6), diameter 2 mm | 1 / 4 / 2 / 5 / 0 |
| Asymmetric figure-eight | Two shared-junction diamond loops at scales 1 and 0.65, diameter 2 mm | 1 / 0 / 1 / 2 / 2 |
| Short communicator | Two long outside paths and a 1.5 mm central connector; major diameter 1.5 mm, connector diameter 1 mm | 1 / 0 / 2 / 3 / 2 |
| Branch near junction | Main junctions 3 mm apart; additional 2 mm terminal of diameter 1 mm beside 2 mm major branches | 1 / 4 / 2 / 5 / 0 |

New masks use exact finite-segment distances at voxel centers, with round caps and physical affine metadata. Truth branch polylines sample at ≤0.08 mm. Spacings are `(0.5,0.5,0.5)` and `(0.46875,0.46875,0.8)` mm; phases are 0, 0.25 and 0.5 voxels in each axis. Selected rotated cases are inherited unchanged from v1.

There are **106 geometries × four methods = 424 main runs**, comparing Lee, raw topology v1, v2 junction-only and full v2 refinement. There are additionally **210 refinement-sweep runs** and **72 endpoint-sweep runs**. These related deterministic observations are not independent statistical replicates.

The controlled one-factor-at-a-time sweep has seven configurations: the declared candidate; radius 0.5 or 1.5 mm instead of 1; spacing 1 mm instead of 0.5; smoothing 0 or 2 mm instead of 1; and spur handling off instead of conservative. It uses 30 geometries: all spacing/phase variants of ring+branches, figure-eight, 1 mm small bridge and near-gap vessels, the four rotated cases, and the two centered double-junction cases. It is not a full factorial interaction study.

There is no uniformly best configuration. The central candidate remains a declared experimental compromise, not an optimized winner: 1 mm radius, 0.5 mm resampling, 1 mm smoothing, conservative mark-only spurs and raw endpoints. The sweep records the tradeoffs below; it does not justify deployment or per-geometry selection.

## Evaluation conventions

V1's one-to-one physical node and structural branch matching are reused. The **2.5 mm gate** is unchanged. Nodes are matched by kind; parallel branches require mapped endpoint pairs and geometric assignment. Unmatched structures are explicit. Physical metrics reuse v1's analytic sampling and coverage tolerance of one maximum voxel spacing. Topology metrics use the retained logical graph, independently of any geometric embedding reduction used internally for physical measurements.

Branch precision/recall/F1 use matched branches, unmatched recovered branches and missed true branches. Adjacency TP/FP/FN compare multiplicity-aware endpoint-pair counts after node assignment; parallel connections count separately. Junction tables include matched/missed/spurious status, physical error and true/recovered degree. Errors average matched structures only and must be read with match counts. Strict connectivity requires every branch matched and no extras; strict topology additionally requires all counts/ranks correct.

Cycles are evaluated by mapped analytic branch membership, closure and perimeter; rank alone is insufficient. The full raw embedding is used for total physical network length. `length_results.csv` retains raw voxel-edge length, resampled branch-path length and returned network length, each with signed bias, absolute and percentage error. Resampling a polyline can shortcut a voxel corner; it is not an exact arc-length-preserving operation after reconnecting its samples. Also, v1 logical reduction contracts junction-tree edges, so resampled branch totals and raw voxel totals are different representations. The raw total is never hidden.

Aggregates are equal-run means, with extraction/acceptance counts shown. No confidence intervals or clinical acceptance thresholds are inferred. All main runs in this recorded matrix extracted and all proposed refinements were accepted; controlled tests separately verify rejection/rollback.

## Core results: the exact 76 v1 geometries

| Metric | Lee | Topology v1 | V2 junction | V2 full |
|---|---:|---:|---:|---:|
| Extracted | 76/76 | 76/76 | 76/76 | 76/76 |
| Component count correct | 70/76 | 76/76 | 76/76 | 76/76 |
| Cycle rank correct | 59/76 | 76/76 | 76/76 | 76/76 |
| Endpoint count correct | 66/76 | 65/76 | 65/76 | 65/76 |
| Junction count correct | 64/76 | 62/76 | 62/76 | 62/76 |
| Branch count correct | 53/76 | 62/76 | 62/76 | 62/76 |
| Strict connectivity correct | 50/76 | 58/76 | 58/76 | 58/76 |
| Branch precision | 0.820468 | 0.831579 | 0.831579 | 0.831579 |
| Branch recall | 0.841009 | 0.830044 | 0.830044 | 0.830044 |
| Branch F1 | 0.820269 | 0.827820 | 0.827820 | 0.827820 |
| Adjacency precision | 0.825731 | 0.831579 | 0.831579 | 0.831579 |
| Adjacency recall | 0.854167 | 0.830044 | 0.830044 | 0.830044 |
| Adjacency F1 | 0.827788 | 0.827820 | 0.827820 | 0.827820 |
| Junction position error (mm) | 0.562848 | 0.349624 | 0.319685 | 0.319685 |
| Endpoint position error (mm) | 0.761569 | 0.708623 | 0.708623 | 0.708623 |
| Network coverage (%) | 90.8276 | 97.6076 | 97.6076 | 90.8395 |
| Network position error (mm) | 0.248008 | 0.211671 | 0.210815 | 0.152562 |
| Network length error (mm) | 5.074722 | 7.198325 | 7.157193 | 2.770747 |

Full refinement reduces mean absolute length error by about **61.5%**, but coverage drops by **6.7681 percentage points**. Interior smoothing/recentering can shortcut a bend while staying inside the binary lumen and keeping graph connectivity. A lower directed position error does not guarantee better analytic network coverage. Junction localization improves modestly; endpoint positions remain fixed by the main configuration. Branch and adjacency F1 do not improve on this cohort.

### Length comparison

| Core length representation | Mean absolute error (mm) | Signed bias (mm) | Mean absolute percent error |
|---|---:|---:|---:|
| Raw voxel graph | 7.198325 | +5.789443 | 12.5824% |
| 0.5 mm resampled branch paths | 5.983632 | +4.174938 | 10.5449% |
| Full refined embedding | 2.770747 | -1.779515 | 5.0304% |

Resampling reduces some stair-step inflation. Full refinement changes the average bias to underestimation; it must not be interpreted as recovering the exact analytic length. The rotated core cases improve from **18.345384 to 1.690040 mm** absolute length error; canonical core cases improve from **6.579044 to 2.830787 mm**. Only four rotated examples were tested.

### Junction, cycle and bridge questions

- **Figure-eight:** all six retain rank two, but correct central junction structure and strict connectivity remain **0/6**. No merge is performed. Both independent memberships are not consistently certified just because rank is correct.
- **Double junction:** correct junction count is **4/6**, strict connectivity **2/6**, unchanged from its raw v1 extraction. The centered test retains two junctions. Refinement never merges a recovered pair, but cannot restore a pair already misrepresented upstream. Mean length error worsens from **10.364329 to 12.222293 mm**.
- **Four-way junction:** one logical junction is retained in all six, but complete connectivity is only **3/6** because other inherited structure/endpoint issues remain.
- **Asymmetric figure-eight:** rank two in 6/6, strict connectivity 0/6. Geometry improves without solving combinatorial structure.
- **Very short communicator:** full connectivity and cycle rank correct in **6/6**; the 1.5 mm-long bridge is retained.
- **Branch near junction:** strict connectivity **3/6**, unchanged; two real short branches are spur candidates but neither is pruned.
- **Existing small bridges:** each of 1.0, 1.5 and 2.0 mm diameter remains recovered and cyclic in **6/6** variants, for 18/18 total.
- **Near-touching:** all 12 retain two components and zero false connections; dense geometric-contact checks pass. These tests do not establish behavior for arbitrary near-crossings or diagonal-only true connections.
- **Incomplete loops:** retain rank zero; refinement does not close the missing connection.

Across all **106 geometries**, v1 and v2 retain 106/106 correct component counts and cycle ranks. Full v2 strict connectivity is **72/106** (58/76 core plus 14/30 added); junction counts are correct in **81/106**. Branch/adjacency F1 is **0.792655**, with precision **0.804403** and recall **0.790881**. Expanded network position error is **0.135679 mm**, length error **3.489001 mm**, coverage **90.3188%**, junction error **0.347122 mm** and endpoint error **0.686646 mm**. Keep these aggregates separate from core comparisons.

## Phase, spacing and endpoint experiments

Core canonical groups exclude the four extra rotations, providing balanced phase comparisons:

| Phase | Strict connectivity | Branch F1 | V2 length error (mm) | V2 junction error (mm) | V2 coverage (%) |
|---|---:|---:|---:|---:|---:|
| Centered | 20/24 | 0.861111 | 1.943556 | 0.140570 | 90.5934 |
| Quarter | 19/24 | 0.848016 | 2.351984 | 0.340654 | 95.5354 |
| Half | 15/24 | 0.745635 | 4.196820 | 0.436738 | 85.5265 |

**Half-phase structural performance does not improve.** Both core spacing groups retain 29/38 strict successes and 38/38 correct ranks. The summary includes full spacing-specific geometry/junction metrics. Changing both in-plane and through-plane spacing does not isolate a causal anisotropy effect.

The endpoint sweep contains Y, incomplete-loop and near-junction-terminal geometries across both spacings/all phases: **18 geometries × four modes**.

| Endpoint mode | Matched endpoint error (mm) | Coverage (%) | Length error (mm) | Strict connectivity |
|---|---:|---:|---:|---:|
| Raw | 0.974222 | 87.9994 | 5.938890 | 9/18 |
| Recentered | 0.887392 | 88.3630 | 5.998270 | 9/18 |
| Tangent extension | 0.924496 | 90.5639 | 4.973199 | 13/18 |
| Boundary constrained | 0.979673 | 90.5639 | 5.035853 | 13/18 |

Extension improves matched connectivity by bringing endpoints within the fixed gate, without changing graph adjacency. It is not a repair of a missing branch. Boundary extension does not improve endpoint error over raw, and tangent extension worsens directed network position error in this subset. Recentered endpoints improve localization but slightly worsen length error. The main candidate retains raw endpoints; extension remains a configurable experiment, not an automatic per-case choice.

## Controlled sweep and failure interpretation

On the 30 sweep geometries, mean absolute length errors for the candidate, radius 0.5, radius 1.5, resampling 1, smoothing 0, smoothing 2 and spur-off are approximately **2.6092, 2.5908, 2.6242, 2.7211, 2.3696, 2.9616 and 2.6092 mm**. Smoothing 0 improves length error and coverage over the central candidate but slightly worsens directed position error. Resampling at 1 mm gives better coverage/position error here but worse length error and coarser sampling. A larger junction radius improves mean junction localization on this cohort but can blend nearby structure. No setting fixes the frozen logical structure. The off/conservative spur settings produce identical paths because neither prunes.

All **212 main refinement proposals** preserve every checked topology invariant and pass containment/contact checks. There are no accepted topology regressions. This perfect before/after agreement is a constraint on refinement, not proof of correct extraction: incorrect branch/junction topology remains incorrect. The default candidate has 34 strict topology failures across the full matrix; they remain in the tables and receive overlays. Important geometric regressions include coverage loss, double-junction length worsening, and inherited endpoint/branch loss that smoothing cannot repair.

## Tests, artifacts and files

**192 tests passed**, including all 139 pre-existing tests and 53 new tests. The same 55 upstream scikit-image/NumPy deprecation warnings remain visible. Existing tolerances were unchanged. New tests cover all 76 frozen raw count results and refined cycle ranks, all 12 near-touching cases, all 18 narrow-bridge cases, invariants/identities, explicit rollback, geometric component contact, true short-branch protection, spur marking, double/four-way junctions, F1 and degree matching, dense containment, deterministic physical resampling, affine rotation, endpoint modes and output protection.

`validation/outputs/topology_v2/manifest.txt` inventories the main outputs; regression outputs have a separate manifest. The primary manifest contains **628 files**, including **81 separate figures**, 106 masks and 424 network JSON records. Main artifacts include:

- `results.csv`: extraction, refinement acceptance, topology/geometry metrics and explicit state.
- `junction_results.csv`: matched/unmatched junction and endpoint positions/degrees.
- `branch_results.csv`, `cycle_results.csv`, `connectivity_results.csv`: structural assignment, lengths, membership and adjacency TP/FP/FN.
- `spur_results.csv`, `junction_decisions.csv`: every diagnostic candidate and conservative decision.
- `length_results.csv`, `matched_comparison.csv`: raw/resampled/refined estimates and paired method deltas.
- `refinement_sweep.csv`, `endpoint_sweep.csv`: all controlled configurations and outcomes.
- `qc_results.csv`, `invariant_results.csv`: criteria, thresholds, before/after values and acceptance.
- `summary.json`, `geometries.json`: grouped results, analytic truth, parameters, versions and hashes.
- `masks/`, `networks/`, `figures/`: all 106 masks, 424 graph/geometry JSON records, summary plots and before/after failure/example overlays. Plots use Matplotlib only and separate figures; projected overlays are labelled x-z while metrics remain fully 3D.
- `regressions/`: five isolated previous runs, comparison evidence and preservation checks.

Only new source/documentation files were created:

```text
neurovasc/graph/topology_invariants.py
neurovasc/graph/junction_refinement.py
neurovasc/geometry/topology_refinement.py
validation/synthetic_topology/experiments_v2.py
validation/synthetic_topology/evaluate_v2.py
validation/synthetic_topology/figures_v2.py
scripts/run_topology_validation_v2.py
scripts/verify_topology_v2_regressions.py
tests/test_topology_validation_v2.py
docs/validation/topology_validation_v2.md
```

Generated artifacts are confined to `validation/outputs/topology_v2/`. No earlier source file or archive was edited, and nothing was pushed.

## Limitations and next decision

The unchanged six-connected extractor remains dependent on mask topology and voxel phase. Geometry-only invariance deliberately cannot solve a wrong logical junction, missing terminal, or spurious branch. Spur candidates include real terminals. Local weighted recentering can blend neighboring structures; smoothing can cut analytic corners inside a wide junction, improving directed error and length while reducing coverage. Fixed endpoint and loop-anchor fits remain scale-dependent. Finite dense containment/contact checks do not prove continuous embedding validity or correct same-component branch interpretation.

The matrix has no real patient anatomy, segmentation noise, imaging partial volume, broad vessel-angle/diameter combinations or clinical reference standard. Parameter selection and evaluation use synthetic cases from the same family; there is no independent holdout or clinical performance claim.

The results support further standalone experiments. Before an experimental hybrid integration, develop and independently validate a declared junction-rewrite policy that distinguishes a true short connector from a split junction, protect real terminal branches, and add a geometric fidelity/coverage criterion to complement containment and topology. **Do not modify the current hybrid or any production default on this evidence.**
