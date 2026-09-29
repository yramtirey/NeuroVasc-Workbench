# Literature review: method questions

Read with the [evidence index](README.md), [source matrix](evidence_matrix.csv), [decision records](methods_decisions.md), [claim map](claim_map.md) and [bibliography](references.bib). External evidence and our synthetic results are identified separately. This focused review establishes precedent and comparison targets, not novelty or clinical validity.

## 1. Circle-of-Willis anatomy and network topology

**External evidence.** Ayre et al.'s systematic review catalogues 82 configurations, including absent, hypoplastic and accessory segments. A fixed branching tree is not an adequate universal anatomical model. TopCoW makes anatomical identity and connections relevant to evaluation. [@ayre2022cow](https://doi.org/10.1111/joa.13616), [@yang2026topcow](https://doi.org/10.1056/AIdbp2500994).

**NeuroVasc / our evidence.** The experimental [TopologyNetwork](../../neurovasc/graph/topology_network.py) uses a MultiGraph with parallel edges and self-loop branches. [Topology v1](../validation/topology_validation_v1.md) preserves cycle rank yet fails figure-eight structure; rank alone is not anatomy. Production's per-label main path is not a complete CoW network representation.

**Alignment: partially aligned; network extraction experimental.** A missing segmented connector may reflect anatomy, image visibility or segmentation error. Neither a label nor an unlabelled cycle identifies a clinically correct variant. No automated variant classifier is validated here. Decision NV-METHOD-001.

## 2. Segmentation and topology-aware evaluation

**External evidence.** The initial TopCoW benchmark and later journal report support multiclass evaluation of CTA/MRA anatomy. clDice complements overlap with skeleton-mask agreement; its theoretical topology result has assumptions and is not a guarantee for arbitrary extracted graphs. [@yang2023topcow](https://arxiv.org/abs/2312.17670v1), [@yang2026topcow](https://doi.org/10.1056/AIdbp2500994), [@shit2021cldice](https://doi.org/10.1109/CVPR46437.2021.01629).

**NeuroVasc / our evidence.** [labels.py](../../neurovasc/segmentation/labels.py) retains laterality and vessel identity: labels 1–12 plus optional 15 (`3rd-A2`); Case 001 contains 12. The pipeline consumes a segmentation; these studies do not train or benchmark a segmentation model. Synthetic topology evaluation measures components, rank, structural matching and coverage, not clDice.

**Alignment: partially aligned.** clDice and per-class segmentation benchmarking are **future work**, not existing features. Dice, clDice, label correctness and matched branch connectivity answer different questions. Preserve missing/extra structures in denominators. NV-METHOD-002/005.

## 3. Centerline extraction

**External evidence.** Lee provides a voxel-thinning method, used by scikit-image. VMTK instead traces surface-derived Voronoi paths with a radius metric. Antiga et al. provide computational-geometry precedent for vascular reconstruction. These are distinct methods, not interchangeable implementations. [@lee1994thinning](https://doi.org/10.1006/cgip.1994.1042), [@skimageSkeletonize](https://scikit-image.org/docs/stable/api/skimage.morphology.html#skimage.morphology.skeletonize), [@vmtkCenterlines](https://vmtk.github.io/tutorials/Centerlines.html), [@antiga2003geometry](https://doi.org/10.1109/TMI.2003.812261).

**NeuroVasc / our evidence.** Production [centerline.py](../../neurovasc/geometry/centerline.py) calls `skeletonize(method="lee")`, then applies the affine. [main_path.py](../../neurovasc/geometry/main_path.py) selects an ordered terminal path. The [geodesic](../../neurovasc/geometry/geodesic_centerline.py) candidate uses a foreground-voxel graph and clearance-weighted costs; it does not compute VMTK's surface Voronoi diagram. [Centerline v1/v2 evidence](README.md#internal-evidence-index) separates recovered failures from positional accuracy.

**Alignment: aligned for the named Lee baseline; alternatives experimental.** Physical EDT spacing does not make thinning spacing-aware. Main-path selection discards other routes. A tree fallback cannot represent cycles; literature precedent does not validate its custom costs, seed inference or thresholds. Retain Lee pending external comparisons. NV-METHOD-003/004.

## 4. Branch decomposition

**External evidence.** Antiga–Steinman and Piccinelli describe vessel decomposition and geometric analysis; the official Slicer implementation invokes both `vtkvmtkCenterlineBranchExtractor` and `vtkvmtkMergeCenterlines`. These grouped centerline tracts are not automatically anatomical vessel labels. [@antiga2004decomposition](https://doi.org/10.1109/TMI.2004.826946), [@piccinelli2009geometry](https://doi.org/10.1109/TMI.2009.2021652), [@slicerExtractSource](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/ExtractCenterline/ExtractCenterline.py).

**NeuroVasc / our evidence.** Legacy [branch_graph.py](../../neurovasc/graph/branch_graph.py) groups junction voxels and degree-two chains into a simple graph. Experimental `topology_network.py` preserves parallel branches, closed branches and cycle memberships. [Topology v1](../validation/topology_validation_v1.md) exposes junction splitting and missing terminals rather than pruning them away.

**Alignment: partially aligned; cycle-aware decomposition experimental.** Do not claim VMTK-equivalent branch definitions or a bijection between labels and graph edges. Benchmark branch identity, multiplicity and attachment, not only counts. NV-METHOD-001/011/012.

## 5. Bifurcation and junction representation

**External evidence.** VMTK's geometric-analysis workflow defines a bifurcation reference system using branch reference points, weighting its origin by inscribed-sphere surface areas. It also provides parallel-transport frames, avoiding some Frenet-frame degeneracies. The documented workflow assumes tree-like topology. [@vmtkGeometry](https://vmtk.github.io/tutorials/GeometricAnalysis.html), [@antiga2004decomposition](https://doi.org/10.1109/TMI.2004.826946), [@piccinelli2009geometry](https://doi.org/10.1109/TMI.2009.2021652).

**NeuroVasc / our evidence.** [junction_refinement.py](../../neurovasc/graph/junction_refinement.py) uses local EDT/distance-weighted voxel centroids and mark-only ambiguity/spur rules. This is not VMTK's reference-point construction. [Topology v2](../validation/topology_validation_v2.md) improves localization but keeps figure-eight connectivity at 0/6. Two real terminals meet its spur-candidate heuristic.

**Alignment: experimental.** Moving junction coordinates cannot repair frozen adjacency. No cited source validates the local 1 mm radius, merge policy or automatic pruning. Compare definitions on bifurcations and close junction pairs before any declared topology rewrite. NV-METHOD-011.

## 6. Vessel radius and caliber

**External evidence.** SciPy's EDT finds distance to background elements with specified axis sampling. VMTK associates a surface-Voronoi centerline with maximum-inscribed-sphere radius. Those quantities do not share the same boundary representation. [@scipyEDT](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.distance_transform_edt.html), [@vmtkCenterlines](https://vmtk.github.io/tutorials/Centerlines.html).

**NeuroVasc / our evidence.** [caliber.py](../../neurovasc/geometry/caliber.py) and [diameter_profile.py](../../neurovasc/geometry/diameter_profile.py) sample twice the spacing-aware binary EDT on the ordered path. [Caliber v1](../validation/caliber_validation_v1.md) quantifies grid-dependent error. It evaluates raw diameters; production [profile_analysis.py](../../neurovasc/geometry/profile_analysis.py) additionally smooths with sigma 1.25 **samples**. Default trimming removes 2 mm from recovered endpoints only if enough samples remain.

**Alignment: aligned as segmentation-derived EDT caliber; clinical interpretation unsupported.** Distance to background voxel centers is not exact distance to a continuous lumen surface. A thinned point is not guaranteed to be a physical medial point. Do not rename EDT output `MaximumInscribedSphereRadius` or infer uncertainty from decimal precision. Literature does not establish these smoothing/trimming settings as optimal. NV-METHOD-006/010.

## 7. Cross-sectional equivalent diameter

**External evidence.** SlicerVMTK distinguishes cross-sectional area-derived circular-equivalent diameter from inscribed-sphere diameter and plots both against centerline distance. This is precedent for a useful geometric estimand, not independent evidence that NeuroVasc computes it accurately. [@slicerCrossSection](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/CrossSectionAnalysis.md).

**NeuroVasc / our evidence.** [cross_sectional_caliber.py](../../neurovasc/geometry/cross_sectional_caliber.py) cuts the affine-transformed binary 0.5 mesh with physical tangent-normal planes, integrates a closed polygon area and returns `D_eq = 2 sqrt(A/pi)`. Its 6 mm tangent window and conservative contour rejection are local choices. [Caliber v2](../validation/caliber_validation_v2.md) supplies the strongest direct evidence for this implementation: better common-support MAE in 78/92 scorable comparisons, with eight EDT wins, six ties and 24 failed pairs upstream. Expanded scorable mean run MAE is 0.229955 mm EDT versus 0.080268 mm cross-sectional.

**Alignment: experimental, leading caliber candidate within the tested matrix.** Area-equivalent diameter is neither a maximal transverse diameter nor necessarily an inscribed-sphere diameter, especially for noncircular sections and tapers. Binary mesh bias, tangent errors, branches/multiple contours and missing centerlines remain. No imaging partial-volume or patient reference study was performed. NV-METHOD-007.

## 8. Tortuosity

**External evidence.** Bernaus et al.'s systematic review distinguishes multiple metric families; a common simple Distance Metric does not constitute a universal cerebrovascular standard. The Slicer documentation exposes a length-to-endpoint-distance ratio. [@bernaus2025tortuosity](https://doi.org/10.1016/j.compbiomed.2025.109990), [@slicerExtract](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/ExtractCenterline.md).

**NeuroVasc / our evidence.** `branch_graph.py` calculates branch length divided by endpoint chord; [metrics.py](../../neurovasc/graph/metrics.py) exposes it. This is an analysis utility, not evidence of a validated patient risk score or a field in every case report. The legacy zero-chord fallback is **1.0**. A closed loop's length/chord ratio is undefined; that fallback must not be interpreted as a straight vessel. Centerline-v2 smoothness diagnostics also depend on sampling and endpoints.

**Alignment: partially aligned.** Use **Distance Metric tortuosity (L/C)** for valid open paths; disclose endpoint/path/sampling choices. Dedicated metric and degenerate-path validation remains future work. Do not equate curvature, torsion, L/C and L/C−1. NV-METHOD-008.

## 9. Intracranial stenosis terminology

**External evidence.** Samuels et al. specify `100(1 − D_stenosis/D_normal)`, including a reference-location protocol: proximal normal artery, with distal or feeding-artery contingencies. Reference selection is part of the measurement, not just the algebra. [@samuels2000stenosis](https://pubmed.ncbi.nlm.nih.gov/10782772/).

**NeuroVasc / our evidence.** `profile_analysis.py` uses the smoothed profile's **75th percentile**, then evaluates reduction at the smoothed minimum; [case_report.py](../../neurovasc/report/case_report.py) exposes `candidate_reduction_percent`. No clinician-selected normal segment, lesion adjudication or agreement against a clinical reference is present. Synthetic caliber accuracy studies do not test stenosis diagnosis.

**Alignment: unsupported as clinical stenosis; aligned with descriptive geometric reporting only.** Name it candidate/relative caliber reduction with its reference definition. Shared ratio algebra does not establish equivalence to the cited protocol. NV-METHOD-009.

## 10. MRA spatial-resolution limitations

**External evidence.** Cosottini et al. and Lakhani et al. illustrate small-vessel depiction using high-resolution 7 T TOF-MRA, with resolution, signal and acquisition tradeoffs. Their limited imaging studies do not calibrate a universal segmentation-derived diameter error bound. [@cosottini2024tof](https://doi.org/10.1186/s41747-024-00463-z), [@lakhani2023tof](https://doi.org/10.1177/19714009221129576).

**NeuroVasc / our evidence.** Case reports retain segmentation voxel spacing; the two caliber grids show measurement sensitivity. Their binary voxel-center masks have no partial-volume, flow-signal, motion or reconstruction model. Changing in-plane and through-plane spacing together does not isolate anisotropy causally.

**Alignment: partially aligned.** Acquired resolution, reconstruction/resampling spacing and segmentation-grid spacing must be distinguished. Original sequence/field strength and segmentation provenance should accompany measurements when available; mark unknown metadata explicitly. The code does not recover these facts from a NIfTI filename. Visibility is not a caliber accuracy certificate. NV-METHOD-010.

## 11. Existing vascular analysis software

**External evidence.** VMTK documents surface Voronoi centerlines, radius metrics and geometric analysis. SlicerVMTK documents approximate network extraction/endpoint discovery, a separate accurate centerline workflow, and cross-section tables. The source confirms branch extraction/merging. [@vmtkCenterlines](https://vmtk.github.io/tutorials/Centerlines.html), [@vmtkGeometry](https://vmtk.github.io/tutorials/GeometricAnalysis.html), [@slicerExtract](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/ExtractCenterline.md), [@slicerExtractSource](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/ExtractCenterline/ExtractCenterline.py), [@slicerCrossSection](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/CrossSectionAnalysis.md).

**NeuroVasc / our evidence.** The [external benchmark](../validation/external_methods_benchmark_v1.md) executes actual VMTK 1.5.1 filters in an isolated environment. Local PyVista/VTK slicing remains a distinct implementation; actual SlicerVMTK CE is unavailable. **Alignment: executed synthetic comparator with seed/support limitations, not clinical equivalence.** Select VMTK/SlicerVMTK as methodological comparators, not ground truth or proof of equivalence. NV-METHOD-012.

### External benchmark contract and remaining extensions

| Comparator | Suitable comparison | Controls and limitations |
|---|---|---|
| VMTK surface centerlines and `MaximumInscribedSphereRadius` | Open tubes and tree-like Y networks; physical path position, endpoints, coverage and radius/diameter | Use identical physical surfaces; record source/target seeds, caps, cleaning, simplification, smoothing and failures. Compare sphere and area estimands separately. |
| SlicerVMTK network extraction | Endpoint discovery and network structure | Evaluate its network output separately from its accurate centerline-tree output; do not assume equivalent graphs. |
| Branch extraction/merging and VMTK bifurcation reference systems | Branch attachments, multiplicity, junction position and geometry | Map identities explicitly; pin inlet/outlet policy. Do not interpret branch GroupIds as anatomical labels. |
| SlicerVMTK Cross-section analysis | Area, equivalent diameter and inscribed-sphere diameter versus distance | First use identical sample locations and surface to isolate section calculation; then compare complete pipelines. Record tangent/section choices, failed contours and RAS/LPS transforms. |
| Documented tree-based workflows on loops | A deliberately limited diagnostic comparison | Slicer's accurate workflow documents interrupting circular topology with nearby outlets. Such a cut changes the problem; report it separately and never score it as preservation of the original CoW network. |

Pin application, extension, VMTK and dependency versions/commits before execution; current documentation links are mutable. Keep an unmodified analytic truth cohort, common-support errors and failure-inclusive coverage. Record own-support metrics as secondary. Include phase/orientation, junctions, noncircular sections and segmentation perturbations in a separate held-out cohort; do not tune on every test shape. An external package agreeing with NeuroVasc is not independent anatomical truth.

## 12. Implications for NeuroVasc

**External evidence.** The cited methods establish substantial prior art; neither centerline profiles nor equivalent diameter nor vascular branch analysis should be presented as a novel invention here. Digital simple-point theory provides a basis for topology-preserving thinning under an explicit adjacency convention, not proof of this application's correctness. [@bertrand1994simple](https://doi.org/10.1016/0167-8655(94)90046-9), [@antiga2003geometry](https://doi.org/10.1109/TMI.2003.812261), [@piccinelli2009geometry](https://doi.org/10.1109/TMI.2009.2021652).

**NeuroVasc / our evidence.** [topology_centerline.py](../../neurovasc/geometry/topology_centerline.py) uses 6-connected foreground/26-connected background and sequential deletion. [topology_refinement.py](../../neurovasc/geometry/topology_refinement.py) holds graph identities fixed. The six [internal reports](README.md#internal-evidence-index) demonstrate tradeoffs, not a single dominant method. Our interpretation is to retain production Lee/EDT and keep alternatives standalone.

**Alignment: partially aligned overall; clinical extension unsupported.** Default retention is a local risk/continuity decision supported by observed regressions, not a literature finding that Lee/EDT is universally best. Put estimand, units, support, failures and method status beside every future benchmark result. Use the [twelve decisions](methods_decisions.md).

## 13. Remaining evidence gaps

**External evidence.** The anatomy and tortuosity reviews describe heterogeneity; neither validates local thresholds. The clinical reference-diameter protocol concerns a different measurement from NeuroVasc's percentile statistic. [@ayre2022cow](https://doi.org/10.1111/joa.13616), [@bernaus2025tortuosity](https://doi.org/10.1016/j.compbiomed.2025.109990), [@samuels2000stenosis](https://pubmed.ncbi.nlm.nih.gov/10782772/).

**NeuroVasc / our evidence.** An [executed VMTK benchmark](../validation/external_methods_benchmark_v1.md) now provides external synthetic comparisons with explicit support and seed limitations. There is still no independent holdout, patient reference standard, reader/retest study or validated uncertainty model. Dataset provenance is documented in [Data Sources](../../DATA_SOURCES.md). Binary phantoms do not establish robustness to partial volume/noise/segmentation variation. Specific local choices remain engineering hypotheses: sigma 1.25 samples, 75th-percentile reference, trim 2 mm, tangent window 6 mm, structural matching gate 2.5 mm, junction radius 1 mm, spur rules and geometric containment/coverage thresholds. The topology convention can reject diagonal-only real connections. Exact graph invariance does not certify correct anatomy.

**Alignment: experimental for these choices; clinical accuracy and risk claims unsupported.** Next implement a **separate, version-pinned VMTK/SlicerVMTK comparison adapter and provenance manifest**, initially on open tubes/Y junctions, with no production switch. Follow with held-out perturbation/loop-aware studies and a declared junction-rewrite validation protocol. Clinical study design would be a separate undertaking. Citation-access and version gaps are recorded in the [source audit](README.md#source-verification-and-provenance).
