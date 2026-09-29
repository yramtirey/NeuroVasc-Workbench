# Topology Validation v1: loop-aware vascular networks

## Scope and decision

This is **synthetic geometric/topological validation, not clinical validation**. Centerline v2's geodesic fallback produces a tree, so its recovery success does not establish preservation of Circle-of-Willis connections. This milestone evaluates full networks, including zero-endpoint rings, parallel routes and independent cycles.

The experimental method preserves component count and cycle rank in all **76** tested geometries, but only **58/76** have fully correct matched branch connectivity. In particular, all six figure-eight cases retain rank two while failing the expected junction/branch structure. Length inflation is substantial. **Keep production Lee, the EDT caliber default, and the hybrid's experimental status unchanged. Do more validation and geometry/junction work before adding a general CoW fallback.** No frontend, API, Case 001, source data, previous analysis algorithms or validation archives were modified.

## Reproduce and preservation

From the repository root:

```bash
source .venv/bin/activate
python3 scripts/run_topology_validation.py
python3 -m pytest -q
```

The runner also works by absolute path. `--output` accepts only `validation/outputs/topology_v1/` or a descendant, protecting the four previous archives. It writes deterministic artifacts without timestamps or random seeds. Dependency versions, source SHA-256 hashes and parameters are in `summary.json`. Re-running replaces this milestone's artifacts; it does not alter production outputs.

All four original runners were executed successfully. Caliber runners used their supported isolated `--output` destinations under `topology_v1/regressions/`. Centerline runners were executed unchanged in a temporary source mirror because their output protection targets their own archive directories; their outputs were copied into the new regression directory. All **16 archived CSVs** reproduced byte-for-byte. SHA-256 comparison found **807 existing files unchanged**, including prior archives, production code and Case 001/data files. Evidence is in `regressions/checks.json`. Regression reproductions have their own manifests and are excluded from the primary artifact manifest.

Original commands remain:

```bash
python3 scripts/run_caliber_validation.py
python3 scripts/run_caliber_validation_v2.py
python3 scripts/run_centerline_validation.py
python3 scripts/run_centerline_validation_v2.py
```

Running those commands directly normally rewrites their respective outputs; the verification above used isolation to preserve the archives.

## Methods and representation

Five methods receive the same binary mask: unchanged Lee full voxel graph, unchanged raw geodesic tree, unchanged refined geodesic tree, current experimental hybrid, and new `extract_topology_network`. Lee's full graph is evaluated even for disconnected inputs; the existing hybrid retains its own single-component QC. Thus success rates can differ without changing the baseline algorithms. Raw/refined trees cannot recover cycles; hybrid can retain them only when it selects Lee.

The new method is **physical-clearance-ordered sequential simple-point thinning**, followed by a cycle-preserving multigraph reduction. It calls neither Lee nor geodesic recovery and receives no analytic truth. Physical EDT sets deletion priority, favoring high-clearance interior voxels. Boundary voxels are processed in increasing clearance with deterministic voxel-index ties. Every deletion rechecks its current local neighborhood; affected neighbors are requeued. Potential endpoints are protected when only one 26-neighbor remains.

The digital topology convention is **6-connected foreground / 26-connected background**. A removable point requires one foreground component touching its six face neighbors within the local 18-neighborhood and one background 26-component in the surrounding neighborhood. This follows the local simple-point approach described by [Bertrand and Malandain](https://www.sciencedirect.com/science/article/pii/0167865594900469); the exact implementation is experimental, not a claim of validated behavior for arbitrary masks. Component invariance is checked at runtime; tests check component and Euler invariance on the declared phantom cases. The surviving voxel graph uses face-neighbor edges with affine-derived physical coordinates and lengths. Orthogonal affines consistent with spacing are required; inputs above 150,000 foreground voxels are explicitly rejected to bound this experimental implementation.

Six-connectivity prevents diagonal-only connections in the near-touching masks, but also cannot represent a genuine connection present only diagonally. This is an explicit digital adjacency choice, not an anatomical classifier or a way to repair a broken segmentation. All matrix masks have the correct six-connected component count before extraction.

`TopologyNetwork` contains the full raw graph and a logical **MultiGraph**, permitting self-loop branches and parallel branches. Raw nodes record physical/voxel coordinates, node kind, clearance where available and branch memberships. Logical nodes infer endpoints from degree one and junction regions from connected clusters of degree-three-or-higher nodes. Degree-two chains become branches; a pure closed component receives a regular anchor, never an invented endpoint.

Inside a junction cluster, only a spanning tree is contracted. Remaining internal edges become explicit cycle branches, so voxel adjacency cycles are retained rather than silently discarded. The reduction checks that `E - N + C` is unchanged. Branches retain polylines, physical lengths, adjacency and cycle memberships. A deterministic **fundamental cycle basis**, not a minimum basis, supports multigraph self-loops and parallel routes. No cycle is pruned by length or any other rule. This exposes baseline voxel-adjacency artifacts and new-method junction splitting rather than hiding either.

Total network length sums raw physical edges exactly once; summing branch-cycle polylines could duplicate junction-cluster paths and is not used for this metric.

## Analytic phantoms and matrix

Voxel centers are classified against analytic physical tubes, with background padding and recorded affines. Phases shift the grid by 0, 0.25 or 0.5 voxels along all axes. Most tubes have round caps; the existing Y's flat-cap inclusion and dimensions are reused. All truth is constructed separately and consulted only during evaluation.

| Geometry | Definition | Components / endpoints / junctions / branches / rank |
|---|---|---|
| Simple ring | Radius 8 mm, diameter 3 mm | 1 / 0 / 0 / 1 / 1 |
| Ring + two branches | Ring with opposite 6 mm radial terminals | 1 / 2 / 2 / 4 / 1 |
| Figure-eight | Two diamond loops sharing a central junction; vertices at origin, (±6,0,±6), (±12,0,0); diameter 3 mm | 1 / 0 / 1 / 2 / 2 |
| Communicating bridge | Theta network: rails at x=±3 mm, z=±8 mm; two 22 mm outside paths and a 6 mm central bridge | 1 / 0 / 2 / 3 / 2 |
| Incomplete loop | Radius-8 ring with 0.7 radian gap, diameter 3 mm | 1 / 2 / 0 / 1 / 0 |
| Y | Existing 22 mm parent, two 18 mm daughters at ±35 degrees; diameters 4/3 mm | 1 / 3 / 1 / 3 / 0 |
| Near-touching | Two 24 mm parallel tubes, diameter 2 mm; gaps 0.4 and 1.2 mm along the transverse diagonal | 2 / 4 / 0 / 2 / 0 |
| Small bridge | Theta network with bridge diameter 1, 1.5 or 2 mm; major routes diameter 3 mm | 1 / 0 / 2 / 3 / 2 |
| Ring + multiple branches | Ring with four 6 mm radial terminals | 1 / 4 / 4 / 8 / 1 |

The figure-eight uses angular loops, not two tangent circles whose thick tubes overlap across a large region. Each ideal loop has perimeter `4*sqrt(72)` mm. This remains a challenging junction phantom, not a smooth model of anatomy.

There are 12 base configurations × two spacings × three phases = 72 geometries. Spacings are `(0.5,0.5,0.5)` and `(0.46875,0.46875,0.8)` mm. Four additional centered cases rotate the ring and communicating-bridge geometry by 0.43 rad about y then 0.31 rad about z, relative to each grid. **76 geometries × five methods = 380 runs**, with 304 paired baseline/new-method comparisons. These deterministic, related cases are not independent statistical replicates.

`geometries.json` stores each analytic node, branch, adjacency, branch length, cycle membership, total length, physical dimensions, phase and orientation. All 76 masks are saved.

## Matching, metrics and QC

Endpoint, junction and regular-anchor matches use separate one-to-one physical assignments under a **2.5 mm gate**. Dummy costs prioritize maximum valid cardinality, then minimum distance. Pure-loop anchors are compared by incident-branch centroid so an arbitrary loop cut does not determine success. Branches must agree on mapped endpoint/junction pairs and pass a symmetric mean nearest-sample distance gate of 2.5 mm; parallel routes are assigned one-to-one by geometry. Unmatched nodes and branches are retained explicitly. These engineering gates are not clinical tolerances.

Count accuracy and connectivity accuracy use all requested runs as denominator, including extraction failures. Full connectivity requires every analytic branch matched with no extra recovered branches. Strict topology additionally requires all component, endpoint, junction, branch and cycle counts correct. Merely recovering rank two does not certify a figure-eight's structure.

Analytic cycle membership is recovered only when all its branches match and form one connected degree-two closed subgraph. Each true cycle records coverage, perimeter and perimeter error when matched. Recovered basis cycles are also recorded. Different valid basis choices are not automatically spurious; `number_of_spurious_cycles` is positive rank excess, while uncertified recovered cycles have a separate count. Missing cycles have no fabricated perimeter error.

Network position error is the mean distance from densely sampled recovered edges to analytic branch samples (maximum analytic step 0.08 mm). Coverage is the fraction of analytic samples within **one maximum voxel spacing** of recovered edges sampled at at most 0.15 mm. Nearest-sample distances approximate continuous distances; their discretization scales are approximately 0.04 and 0.075 mm respectively. Duplicate edge endpoints mean these are sample averages, not exact arc-length integrals. Endpoint/junction errors average matched structures only and are accompanied by matched/unmatched counts. Network length error is absolute raw-edge length minus exact analytic length. Mean/minimum clearance sample the physical EDT. Failures retain reasons and missed truth structures, without invented zero geometry errors.

Aggregate geometry means are equal-run means over successful extractions. Failure-inclusive coverage is also reported, assigning failed runs zero coverage; positional/length comparisons must account for the different success sets. `matched_comparison.csv` gives pairwise differences where both methods succeeded. No confidence intervals are inferred.

For near-touching vessels, recovered nodes are attributed to the nearest analytic component. Actual edges crossing those labels and connected recovered components covering both true vessels identify false connections. An entirely absent vessel is recorded separately; fewer recovered components alone is not evidence of a false bridge.

Truth-independent QC checks finite coordinates, positive edge lengths, isolated nodes, raw self-loops, duplicate edges, component existence, simplification rank invariance and cycle capability. Multiple components are permitted. Short branches below one maximum voxel spacing, degree above six, and cycle perimeter below four maximum spacings are diagnostic engineering heuristics. They do not trigger pruning. Legitimate logical self-loop branches are distinct from malformed raw self-loop edges. `qc_results.csv` records criterion, value, threshold, pass/fail and criticality. A QC pass is not proof of correct anatomy.

## Recorded results

Counts below use denominator 76 for every method. Connectivity also equals strict topology success in this run.

| Method | Extracted | Components correct | Endpoints correct | Junctions correct | Branches correct | Rank correct | Connectivity correct |
|---|---:|---:|---:|---:|---:|---:|---:|
| Lee | 76 | 70 | 66 | 64 | 53 | 59 | 50 |
| Raw geodesic | 68 | 64 | 13 | 39 | 12 | 16 | 12 |
| Refined geodesic | 66 | 63 | 12 | 38 | 12 | 15 | 12 |
| Current hybrid | 69 | 64 | 61 | 62 | 47 | 53 | 46 |
| New topology | 76 | 76 | 65 | 62 | 62 | 76 | 58 |

| Method | Coverage, successful (%) | Coverage, all runs (%) | Mean position error (mm) | Mean absolute network length error (mm) |
|---|---:|---:|---:|---:|
| Lee | 90.8276 | 90.8276 | 0.248008 | 5.074722 |
| Raw geodesic | 88.6254 | 79.2964 | 0.230595 | 4.691962 |
| Refined geodesic | 93.0156 | 80.7767 | 0.175770 | 3.461176 |
| Current hybrid | 94.1512 | 85.4794 | 0.226730 | 3.750573 |
| New topology | 97.6076 | 97.6076 | 0.211671 | 7.198325 |

The new method improves cycle-rank accuracy from Lee's **77.63% to 100%** and connectivity from **65.79% to 76.32%**. Endpoint and junction count accuracy are slightly worse than Lee. Positional error improves on average, but total length error increases by **2.123603 mm (41.8%)**. Six-neighbor stair-step paths explain an important geometric cost. The four rotated cases have correct connectivity but mean length error **18.345384 mm**, versus 6.579044 mm across canonical cases. This is not a claim of topology improvement without material geometry degradation.

### Loop and branch design review

| New-method geometry | Rank correct | Complete connectivity correct | Interpretation |
|---|---:|---:|---|
| Simple ring | 8/8 | 6/8 | Cycle survives; two cases have terminal spurs |
| Ring + two branches | 6/6 | 4/6 | Half-phase terminal/junction loss in two cases |
| Figure-eight | 6/6 | 0/6 | Both independent ranks survive; central junction splits and/or spurs appear; no full-network success |
| Communicating bridge | 8/8 | 8/8 | Parallel routes and central bridge match, including selected rotations |
| Incomplete loop | 6/6 | 5/6 | Never creates a false cycle; one endpoint fails the matching gate |
| Y | 6/6 | 1/6 | Endpoint erosion/gating and two half-phase branch losses; geodesic trees match all six Y cases |
| Ring + four branches | 6/6 | 4/6 | Half-phase cases lose a terminal branch/junction |
| Small bridges | 18/18 | 18/18 | All tested bridge variants preserve the theta network |
| Near-touching | 12/12 | 12/12 | Both components and complete branch connectivity retained |

The new method outperforms geodesic trees on true loops: the latter have no loop recovery by construction. Rank preservation alone is insufficient for the figure-eight; both analytic cycle memberships are not consistently certified when the central junction is split.

### Narrow communicating bridges

At each diameter **1.0, 1.5 and 2.0 mm**, the new method recovers the bridge and correct rank in **6/6** spacing/phase variants, as do Lee and the current hybrid. Raw/refined trees have **0/6** complete bridge matches and correct cycle ranks at each diameter; this does not mean no fragment of the physical connector is present in their trees. No cycles or bridges were pruned by the new method. Its mean full-network positional error is 0.205029 mm at each tested diameter; branch-specific positional errors are in `branch_results.csv`.

**1.0 mm is the smallest tested diameter with complete observed recovery**, not an established general reliability threshold. No diameter below 1 mm, noise, partial-volume effect or broad branch-angle sweep was tested; the point at which bridge recovery becomes unreliable remains unknown.

### Near-touching vessels

| Gap (mm) | Method | Extracted / 6 | False-connection runs | Missing-vessel runs among extractions |
|---|---|---:|---:|---:|
| 0.4 | Lee | 6 | 4 | 1 |
| 0.4 | Raw geodesic | 4 | 4 | 0 |
| 0.4 | Refined geodesic | 3 | 3 | 0 |
| 0.4 | Hybrid | 4 | 3 | 1 |
| 0.4 | New topology | 6 | 0 | 0 |
| 1.2 | Lee | 6 | 0 | 1 |
| 1.2 | Raw geodesic | 0 | — | — |
| 1.2 | Refined geodesic | 0 | — | — |
| 1.2 | Hybrid | 1 | 0 | 1 |
| 1.2 | New topology | 6 | 0 | 0 |

All masks have two six-connected components. A full foreground 26-neighbor graph instead connects the vessels in **4/6** near-gap masks and **0/6** clear-gap masks. The new method makes no false connections in 12 cases. Lee's missing-vessel cases must not be counted as false bridges. The existing geodesic input restriction rejects disconnected masks, explaining the clear-gap failures; zero successful extractions cannot establish a zero false-connection rate. Structured rates use successful extractions as denominator, return null when none succeeded, and retain failure counts. Per-run cross-component edge counts and recovered component counts are in `results.csv`.

### Phase and spacing sensitivity

For the balanced canonical matrix, new-method connectivity is **20/24 (83.33%)** centered, **19/24 (79.17%)** quarter-phase and **15/24 (62.50%)** half-phase. Rank remains correct in all groups. Mean position error across the recorded phase groups increases from **0.144587** to **0.206988** to **0.294620 mm**; the first group includes four additional rotated centered cases, so it is not perfectly balanced. Summary/figures preserve the group counts rather than implying identical cohorts.

Connectivity is **29/38 (76.32%)** for each spacing, with rank correct for all. Position error increases from **0.199564 to 0.223779 mm** on the anisotropic grid; coverage is **97.6052% versus 97.6100%**, and length error **7.215841 versus 7.180809 mm**. In-plane spacing changes as well, so this does not isolate a causal anisotropy effect. Four selected orientation variants are insufficient to establish rotational robustness.

## Failures, tests and artifacts

All **18 new-method strict topology failures** have diagnostic overlays. They include spurious terminal spurs on rings, junction splitting on figure-eights, shortened endpoints outside the match gate, and lost terminal branches at half phase. Other methods have 25 extraction failures in total: 22 `disconnected_mask`, two `refined_segment_outside_lumen`, one `resampled_point_outside_lumen`. Every failed extraction and topology mismatch remains in the result tables with reasons; count accuracy does not exclude them.

**139 tests passed**, including all 107 existing tests and 32 new cases. The 55 existing scikit-image/NumPy deprecation warnings remain visible. Tests cover all phantom truths, both spacings, bridge diameters, near-touching 6/26 components, loop/parallel-edge ranks, junction contraction, branch matching and node reordering, false-cycle prevention, incomplete loops, physical affine rotation, deterministic extraction, absence of truth/Lee access in extraction, explicit failures, and output-directory protection. Existing tolerances were not weakened. The four isolated complete reruns provide additional numerical regression evidence.

Primary manifest contains **689 files**, including **248 separate Matplotlib figures** (12 summary figures plus 236 examples/failure diagnostics), 76 masks and 355 recovered-network JSON files. Examples overlay analytic and recovered networks, distinguishing matched, missed and spurious branches. Diagnostic projections are explicitly x-z; all metrics use 3D coordinates. A failed extraction still gets an analytic-only diagnostic.

- `results.csv` and `network_metrics.csv`: complete per-run metrics and failure diagnostics; the latter intentionally provides the same full metric view.
- `branch_results.csv`: every matched, missed or spurious branch, lengths and positional errors when defined.
- `node_matches.csv`: matched and unmatched endpoint/junction/regular nodes.
- `cycle_results.csv`: analytic cycle membership, coverage, perimeter and recovered basis diagnostics.
- `matched_comparison.csv`: 304 paired comparisons, retaining extraction failures.
- `qc_results.csv`: truth-independent QC values, thresholds and criticality.
- `summary.json`, `geometries.json`: aggregates, complete truth, parameters, versions and hashes.
- `networks/`, `masks/`, `figures/`, `manifest.txt`: full reproducible evidence.
- `regressions/`: isolated previous-run artifacts and preservation checks, separately inventoried.

## New files and limitations

Only new files were added:

- `neurovasc/geometry/topology_centerline.py`: independent loop-aware thinning extractor.
- `neurovasc/graph/topology_network.py`: full graph/multigraph, conservative reduction and cycle basis.
- `neurovasc/graph/topology_qc.py`: truth-independent diagnostics.
- `validation/synthetic_topology/__init__.py`, `geometries.py`, `evaluate.py`, `figures.py`: phantom library, matching, metrics and plotting.
- `scripts/run_topology_validation.py`: protected deterministic experiment runner.
- `tests/test_topology_validation.py`: new regression tests.
- `docs/validation/topology_validation_v1.md`: this report.
- `validation/outputs/topology_v1/**`: exact generated-file names are in the primary and regression manifests.

The method is not ready as a general Circle-of-Willis fallback. Its good cycle counts are conditional on suitable digital mask topology. Six-connectivity, index-order tie breaking, endpoint protection, local junction clustering and stair-step geometry introduce orientation/phase dependence. It does not infer anatomical branch identity or repair segmentation gaps. The matrix excludes patient anatomy, realistic CoW asymmetry, vessel segmentation errors, noise, partial volume, close crossings in many orientations, stenoses and large volumes. No disease inference or clinical validation is supported.

Next work should address endpoint/junction stability and physical path length while retaining verified cycles, then test diagonal-only true connections, broader phase/orientation combinations, curved/nonplanar communicating networks and segmentation perturbations. Keep the new extractor available as an **experimental standalone adapter**, without modifying or deploying the current hybrid. There is insufficient evidence to change any production default.
