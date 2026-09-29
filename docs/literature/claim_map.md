# Claim map and terminology lock

Use this when drafting public descriptions. See the [method review](literature_review.md), [decisions](methods_decisions.md), [evidence matrix](evidence_matrix.csv), [bibliography](references.bib) and [validation index](README.md#internal-evidence-index). External precedent does not validate our implementation. The evidence currently supports engineering/synthetic claims, not clinical claims.

## Claim map

| Claim | Evidence supporting claim (external) | Internal validation supporting claim | Allowed wording | Disallowed wording | Citation(s) |
|---|---|---|---|---|---|
| Caliber profiles are computed | EDT definition provides computational meaning | Production source; caliber v1 documents its errors | NeuroVasc computes segmentation-derived caliber profiles using spacing-aware EDT. | NeuroVasc measures the true lumen diameter precisely. | [@scipyEDT](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.distance_transform_edt.html) |
| Candidate caliber reduction is available | Clinical reference protocol illustrates a different estimand | Code uses smoothed 75th-percentile reference; no diagnostic study | Candidate caliber reduction is an internal geometric descriptor. | NeuroVasc accurately diagnoses stenosis or grades disease severity. | [@samuels2000stenosis](https://pubmed.ncbi.nlm.nih.gov/10782772/) |
| Cross-sectional method is promising | Existing area-derived diameter software precedent | Caliber v2: 78 wins, 8 losses, 6 ties; 24 unscorable pairs | Experimental cross-sectional equivalent diameter outperformed EDT on most tested scorable synthetic comparisons. | Cross-sectional diameter is clinically accurate or always superior. | [@slicerCrossSection](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/CrossSectionAnalysis.md) |
| Network topology is analyzed | Anatomy variation and topology-aware segmentation motivate connections | Topology v1/v2 expose rank versus connectivity differences | NeuroVasc experimentally analyzes vascular network topology. | NeuroVasc identifies clinically correct CoW variants. | [@ayre2022cow](https://doi.org/10.1111/joa.13616), [@yang2026topcow](https://doi.org/10.1056/AIdbp2500994) |
| Cycle preservation was observed | Digital topology supplies a method family, not a software certificate | Core topology: rank correct in 76/76; connectivity 58/76 | The experimental extractor preserved cycle rank in the tested cohort. | It guarantees correct anatomical connectivity for every vessel network. | [@bertrand1994simple](https://doi.org/10.1016/0167-8655(94)90046-9) |
| Geometry refinement helps some metrics | Branch geometry has established methodological precedent | Topology v2 length/position improve; coverage worsens | Refinement reduced core mean length error while reducing analytic coverage. | Refinement fixed CoW topology or improved every geometric measure. | [@piccinelli2009geometry](https://doi.org/10.1109/TMI.2009.2021652) |
| Anatomical labels are retained | TopCoW multiclass convention | Local label map and Case 001 report, not label-accuracy validation | NeuroVasc retains input anatomical vessel labels and laterality. | Twelve labels describe every normal CoW, or label presence proves variant diagnosis. | [@yang2023topcow](https://arxiv.org/abs/2312.17670v1) |
| Centerline recovery was improved | Alternative path constructions exist | Centerline v2 recovered 12 known failures; tree limitation remains | Experimental hybrid recovery succeeded on the 60-case core cohort. | A clinically validated, loop-preserving fallback is deployed. | [@vmtkCenterlines](https://vmtk.github.io/tutorials/Centerlines.html) |
| Distance Metric tortuosity is calculated | Multiple metric families; no universal standard | Legacy L/C implementation; degenerate chord handling unresolved | The analysis utility calculates Distance Metric tortuosity for valid open branches. | NeuroVasc provides the standardized cerebrovascular tortuosity score or a validated risk biomarker. | [@bernaus2025tortuosity](https://doi.org/10.1016/j.compbiomed.2025.109990) |
| Resolution matters | TOF imaging studies show acquisition-dependent depiction | Caliber studies vary grid spacing/phase, not imaging physics | Results depend on segmentation-grid resolution and require acquisition context. | Submillimeter sampling proves submillimeter accuracy. | [@cosottini2024tof](https://doi.org/10.1186/s41747-024-00463-z), [@lakhani2023tof](https://doi.org/10.1177/19714009221129576) |
| Topology-aware segmentation metrics are relevant | clDice provides skeleton-aware comparison | Not implemented in the current validation suite | clDice is a proposed additional segmentation comparator. | NeuroVasc already achieves clDice-guaranteed topology. | [@shit2021cldice](https://doi.org/10.1109/CVPR46437.2021.01629) |
| External software is a planned comparator | Official documentation and branch extraction/merging source | No benchmark against VMTK/SlicerVMTK has run | VMTK/SlicerVMTK are selected methodological comparison targets. | NeuroVasc matches VMTK accuracy or is clinically validated by Slicer. | [@slicerExtractSource](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/ExtractCenterline/ExtractCenterline.py) |
| Methods have prior art | Published centerline, decomposition and geometry methods | This project integrates and evaluates its own implementations | NeuroVasc combines established geometric ideas with documented synthetic comparisons. | NeuroVasc invented centerline profiles, equivalent diameter or vascular branch decomposition. | [@antiga2003geometry](https://doi.org/10.1109/TMI.2003.812261), [@antiga2004decomposition](https://doi.org/10.1109/TMI.2004.826946) |

All numerical internal claims above refer to the archived reports linked in the [evidence index](README.md#internal-evidence-index). Describe the cohort and failure denominator when reusing a number. The paired trim settings do not create independent statistical replicates.

## Terminology lock

| Preferred term | Required meaning / qualification |
|---|---|
| **Segmentation-derived caliber** | Geometric estimate from a segmentation; state EDT versus surface-area method, units and smoothing. |
| **Cross-sectional equivalent diameter** | `2 sqrt(A/pi)` from a defined plane and section area; not the maximum chord or necessarily the inscribed-sphere diameter. “Circular-equivalent diameter” is the corresponding Slicer documentation term. |
| **Maximum-inscribed-sphere radius/diameter** | Use for an actual inscribed-sphere construction and state radius versus twice radius. Do not relabel the current voxel EDT as VMTK's surface-derived measurement. |
| **Relative caliber reduction** | Descriptive ratio with its reference definition, smoothing and support stated. |
| **Candidate caliber reduction** | Preferred public description of `candidate_reduction_percent`; reference is the smoothed profile's 75th percentile. |
| **Distance Metric tortuosity** | Here `L/C`, for positive endpoint chord; identify path/endpoints. Zero chord is undefined, despite the legacy 1.0 fallback. |
| **Network topology** | Connections, components and cycles under a declared graph/digital adjacency model. |
| **Cycle rank** | `E − N + C`, independent-cycle count; not full anatomical network identity or proof of correct cycle membership. |
| **Branch connectivity** | Matched branch attachments and multiplicity, including missing and extra structures; more specific than rank. |
| **Anatomical vessel label** | Input anatomical identity, separate from a recovered branch ID or inferred variant. |

Avoid **stenosis**, **disease severity**, **diagnostic finding**, **clinical-grade**, **clinical validation**, **patient-specific risk**, and **validated biomarker** as claims about current NeuroVasc outputs. These words can appear when discussing an external study or explicitly stating the absence of such validation; they must not be inherited from a citation. Do not apply a clinical threshold to the percentile-based statistic.

## Evidence required for stronger language

“Synthetic validation” requires a stated phantom matrix, truth, metrics and failures. “Externally benchmarked” requires an executed, version-pinned comparison, not a software citation. “Clinically validated” would require a separately designed study with suitable patient/reference data and intended-use endpoints. None of the current six synthetic milestones supplies that clinical evidence. The [decision records](methods_decisions.md) specify what would justify revisiting each method.
