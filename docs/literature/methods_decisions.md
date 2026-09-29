# Method decision records

These records describe the inspected implementation as of 2026-09-28. They do not change code or approve deployment. See the [review](literature_review.md), [matrix](evidence_matrix.csv), [claim map](claim_map.md), [bibliography](references.bib) and [internal evidence index](README.md#internal-evidence-index). Production is a deployment status, not a clinical quality designation.

## NV-METHOD-001 — Network rather than a universal tree

- **Method:** retain loops, parallel branches and disconnected components in the general CoW representation.
- **Current NeuroVasc implementation:** experimental [TopologyNetwork](../../neurovasc/graph/topology_network.py) is a MultiGraph. The production per-label main path and legacy branch graph are narrower representations.
- **External evidence:** anatomical variation and topology-aware segmentation motivate retaining connections; they do not prescribe this implementation. [@ayre2022cow](https://doi.org/10.1111/joa.13616), [@yang2026topcow](https://doi.org/10.1056/AIdbp2500994).
- **Internal validation evidence:** [Topology v1](../validation/topology_validation_v1.md): 76/76 correct ranks but 58/76 complete connectivity.
- **Known limitations:** rank and component count cannot certify vessel identity or junction arrangement.
- **Current status:** experimental network adapter; validation-only structural benchmarks.
- **Decision:** keep a network-capable representation as the target; never silently force a CoW loop into a tree.
- **What evidence would change the decision:** validated anatomy mapping could support integration; any intentionally tree-only use must declare its restricted region and discarded routes.

## NV-METHOD-002 — Retain anatomical labels

- **Method:** keep anatomical vessel labels separate from graph branch IDs.
- **Current NeuroVasc implementation:** [labels.py](../../neurovasc/segmentation/labels.py) maps 1–12 and optional 15; Case 001 has 12 labelled vessels. Background is zero.
- **External evidence:** multiclass TopCoW conventions preserve distinctions lost in a merged binary mask. [@yang2023topcow](https://arxiv.org/abs/2312.17670v1), [@yang2026topcow](https://doi.org/10.1056/AIdbp2500994).
- **Internal validation evidence:** current case/report code preserves labels; topology phantoms test structural truth, not patient label correctness.
- **Known limitations:** dataset provenance and label-version mapping require documentation; absent output is not proof of absent anatomy.
- **Current status:** production label handling; anatomical variant classification is future work.
- **Decision:** preserve identity/laterality and unknown/absent states; do not infer diagnoses from labels.
- **What evidence would change the decision:** a documented dataset ontology migration, expert annotation and per-label/variant evaluation would justify a new mapping or classifier.

## NV-METHOD-003 — Retain the Lee production baseline

- **Method:** Lee voxel skeletonization followed by ordered path extraction.
- **Current NeuroVasc implementation:** [centerline.py](../../neurovasc/geometry/centerline.py) calls scikit-image Lee; world coordinates are obtained afterward.
- **External evidence:** published thinning method and official implementation documentation. [@lee1994thinning](https://doi.org/10.1006/cgip.1994.1042), [@skimageSkeletonize](https://scikit-image.org/docs/stable/api/skimage.morphology.html#skimage.morphology.skeletonize).
- **Internal validation evidence:** [Centerline v1](../validation/centerline_validation_v1.md): 48/60 extraction success; lower position/tangent errors than geodesic on shared successes.
- **Known limitations:** empty skeletons, phase dependence, spacing-insensitive thinning and endpoint shortening; main-path selection loses other routes.
- **Current status:** production baseline, not a clinical accuracy claim.
- **Decision:** retain for continuity while recording failures; literature is not evidence of universal superiority.
- **What evidence would change the decision:** independent comparison demonstrating improved recovery, geometry, topology and downstream caliber without hidden failure exclusions.

## NV-METHOD-004 — Keep geodesic/hybrid recovery experimental

- **Method:** clearance-weighted voxel paths, physical refinement and QC-gated hybrid selection.
- **Current NeuroVasc implementation:** [geodesic](../../neurovasc/geometry/geodesic_centerline.py), [refinement](../../neurovasc/geometry/centerline_refinement.py), [hybrid](../../neurovasc/geometry/hybrid_centerline.py); no production fallback installed.
- **External evidence:** radius-weighted path methods have precedent, but VMTK operates on a different representation. [@vmtkCenterlines](https://vmtk.github.io/tutorials/Centerlines.html), [@antiga2003geometry](https://doi.org/10.1109/TMI.2003.812261).
- **Internal validation evidence:** [Centerline v2](../validation/centerline_validation_v2.md) preserves 48 Lee paths and recovers 12 known core failures.
- **Known limitations:** fallback is a tree; QC passing is not evidence of loop retention. Cost/seed/branch rules are local engineering choices.
- **Current status:** experimental; sweeps validation-only.
- **Decision:** no automatic deployment to general CoW anatomy.
- **What evidence would change the decision:** independent, failure-inclusive geometry/caliber tests plus explicit topology capability and supported-input boundaries.

## NV-METHOD-005 — Keep loop-aware extraction standalone

- **Method:** simple-point thinning and cycle-preserving multigraph reduction.
- **Current NeuroVasc implementation:** [topology_centerline.py](../../neurovasc/geometry/topology_centerline.py), with 6/26 adjacency; v2 refines geometry under exact logical identities.
- **External evidence:** simple-point theory is a methodological basis; topology-aware segmentation motivates structural checks. [@bertrand1994simple](https://doi.org/10.1016/0167-8655(94)90046-9), [@shit2021cldice](https://doi.org/10.1109/CVPR46437.2021.01629).
- **Internal validation evidence:** [Topology v1/v2](README.md#internal-evidence-index): rank preserved; figure-eight structure unresolved and refinement coverage regresses.
- **Known limitations:** adjacency convention, phase dependence, endpoint loss and junction splitting. Local theory does not prove all implementation behavior.
- **Current status:** standalone experimental; not part of the hybrid.
- **Decision:** retain failures and exact graph invariants; no automatic spur removal or junction merge.
- **What evidence would change the decision:** independently tested rewrite policy with identity mappings, real-terminal protection and geometry/coverage criteria, including diagonal connections and perturbations.

## NV-METHOD-006 — Retain EDT caliber as production baseline

- **Method:** twice physical EDT sampled along the ordered centerline.
- **Current NeuroVasc implementation:** [diameter_profile.py](../../neurovasc/geometry/diameter_profile.py); raw diameter is distinct from downstream smoothed statistics.
- **External evidence:** official EDT definition supports the computation, not clinical lumen accuracy. [@scipyEDT](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.distance_transform_edt.html).
- **Internal validation evidence:** [Caliber v1](../validation/caliber_validation_v1.md) records thin-vessel/grid errors; [v2](../validation/caliber_validation_v2.md) retains EDT wins as well as losses.
- **Known limitations:** background voxel-center distance, non-medial samples, sample-index smoothing and endpoint effects. No universal bias correction follows.
- **Current status:** production.
- **Decision:** retain named baseline pending external comparison; label output segmentation-derived caliber.
- **What evidence would change the decision:** held-out, common-support and failure-inclusive comparison against a declared geometric reference, with reproducibility and acquisition/segmentation sensitivity quantified.

## NV-METHOD-007 — Lead with cross-sectional equivalent diameter in experiments

- **Method:** physical tangent-normal mesh section area and `2 sqrt(A/pi)`.
- **Current NeuroVasc implementation:** [cross_sectional_caliber.py](../../neurovasc/geometry/cross_sectional_caliber.py), single closed contour required; invalid sections remain invalid.
- **External evidence:** SlicerVMTK provides the corresponding area-derived diameter concept. [@slicerCrossSection](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/CrossSectionAnalysis.md).
- **Internal validation evidence:** [Caliber v2](../validation/caliber_validation_v2.md): 78 wins/8 losses/6 ties on scorable pairs; 24 upstream failures remain.
- **Known limitations:** binary surface bias, tangent window, noncircular estimand differences and untested complex contours; no patient validation.
- **Current status:** experimental; leading candidate only within current comparisons.
- **Decision:** prioritize controlled external benchmarking, without replacing production EDT.
- **What evidence would change the decision:** independent perturbation/holdout results and area/diameter agreement against analytic truth and a pinned comparator, including contour failures and downstream centerline sensitivity.

## NV-METHOD-008 — Name Distance Metric tortuosity precisely

- **Method:** `L/C`, path length divided by endpoint chord, for open paths with positive chord.
- **Current NeuroVasc implementation:** [branch_graph.py](../../neurovasc/graph/branch_graph.py) and [metrics.py](../../neurovasc/graph/metrics.py); zero-chord fallback returns 1.0.
- **External evidence:** review distinguishes metric families; software documentation defines this ratio. [@bernaus2025tortuosity](https://doi.org/10.1016/j.compbiomed.2025.109990), [@slicerExtract](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/ExtractCenterline.md).
- **Internal validation evidence:** centerline-v2 smoothness diagnostics exist, but no dedicated clinical tortuosity validation or comprehensive degenerate-path benchmark.
- **Known limitations:** endpoint/sampling dependence; ratio undefined for a zero chord. Legacy 1.0 is not anatomical evidence.
- **Current status:** implemented analysis utility (production code); standalone metric validation and explicit degenerate handling are future work.
- **Decision:** lock terminology and disclose formula; no standardized-score or biomarker claim.
- **What evidence would change the decision:** dedicated analytic/degenerate tests and comparator agreement; any clinical interpretation additionally needs a population/reference study.

## NV-METHOD-009 — Candidate reduction is not stenosis

- **Method:** reduction relative to the smoothed profile's 75th-percentile reference.
- **Current NeuroVasc implementation:** [profile_analysis.py](../../neurovasc/geometry/profile_analysis.py), exposed as `candidate_reduction_percent` by [case_report.py](../../neurovasc/report/case_report.py).
- **External evidence:** the Samuels protocol uses an anatomically selected normal reference. [@samuels2000stenosis](https://pubmed.ncbi.nlm.nih.gov/10782772/).
- **Internal validation evidence:** caliber studies test synthetic geometric diameter, not lesions, normal reference selection or diagnostic performance.
- **Known limitations:** taper/voxelization/path selection can affect the statistic; similar algebra does not make reference definitions equivalent.
- **Current status:** production geometric descriptor; clinical stenosis interpretation unsupported.
- **Decision:** use candidate/relative caliber reduction and disclose reference/smoothing; no disease-severity thresholds.
- **What evidence would change the decision:** a separately specified clinical measurement protocol, expert reference selection and appropriate agreement/diagnostic validation; this percentile statistic must not silently inherit another method's name.

## NV-METHOD-010 — Attach resolution and provenance

- **Method:** distinguish acquisition, reconstruction and segmentation grids in measurement provenance.
- **Current NeuroVasc implementation:** [analyze_case.py](../../scripts/analyze_case.py) writes segmentation spacing into the case report; source NIfTI affine supports physical coordinates. Complete original acquisition metadata is not established.
- **External evidence:** TOF studies demonstrate resolution-dependent depiction and acquisition tradeoffs. [@cosottini2024tof](https://doi.org/10.1186/s41747-024-00463-z), [@lakhani2023tof](https://doi.org/10.1177/19714009221129576).
- **Internal validation evidence:** caliber-v1/v2 phase/spacing effects; binary masks omit imaging physics.
- **Known limitations:** voxel spacing is not effective resolution or an uncertainty interval; protocol and segmentation history may be missing.
- **Current status:** production spacing metadata; richer provenance is future work.
- **Decision:** future exports should carry acquisition/reconstruction/resampling/segmentation provenance and explicitly unknown fields; no fabricated metadata for Case 001.
- **What evidence would change the decision:** verified source records and repeat/perturbation measurements could support richer uncertainty estimates, not just more decimal places.

## NV-METHOD-011 — Benchmark junction definitions externally

- **Method:** compare branch-aware/radius-aware junction definitions before structural rewriting.
- **Current NeuroVasc implementation:** [junction_refinement.py](../../neurovasc/graph/junction_refinement.py) moves existing nodes with EDT-weighted centroids; does not merge/split graph identities.
- **External evidence:** published decomposition and VMTK reference-system definitions provide comparator choices. [@antiga2004decomposition](https://doi.org/10.1109/TMI.2004.826946), [@piccinelli2009geometry](https://doi.org/10.1109/TMI.2009.2021652), [@vmtkGeometry](https://vmtk.github.io/tutorials/GeometricAnalysis.html).
- **Internal validation evidence:** [Topology v2](../validation/topology_validation_v2.md): junction localization improves modestly; connectivity/F1 stay unchanged on the core; marked spurs include two real terminals.
- **Known limitations:** shared geometric coordinates and correct adjacency are different properties; tree assumptions need scoping.
- **Current status:** refinement experimental; two Y geometries compared against actual VMTK reference systems in the external benchmark.
- **Decision:** no automatic merge or prune; test true short connectors against split junctions.
- **What evidence would change the decision:** independently successful structural assignments with an explicit permitted-change contract and failure/rollback rules.

## NV-METHOD-012 — Choose VMTK/SlicerVMTK as comparators

- **Method:** external surface/centerline/branch/cross-section comparison.
- **Current NeuroVasc implementation:** [actual VMTK benchmark](../validation/external_methods_benchmark_v1.md) with an isolated pinned runtime; see the [contract](literature_review.md#external-benchmark-contract-and-remaining-extensions).
- **External evidence:** official capabilities and implementation, supported by methodological papers. [@vmtkCenterlines](https://vmtk.github.io/tutorials/Centerlines.html), [@slicerExtractSource](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/ExtractCenterline/ExtractCenterline.py), [@slicerCrossSection](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/CrossSectionAnalysis.md).
- **Internal validation evidence:** 68 geometries; 58 direct, 2 partial and 8 assumption mismatches. Actual VMTK results use separate seeds and common support; Slicer CE remains unexecuted.
- **Known limitations:** seed/cap/surface preparation, coordinate systems and tree assumptions can change the task. Software agreement is not ground truth.
- **Current status:** standalone external synthetic benchmark; Slicer CE and open-profile variants remain future work.
- **Decision:** retain all production defaults; current capped-surface seed results do not establish general VMTK performance.
- **What evidence would change the decision:** inability to run a compatible reproducible configuration, or evidence that a different independently maintained comparator better matches the declared estimand/topology.
