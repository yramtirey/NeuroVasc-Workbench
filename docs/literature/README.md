# NeuroVasc methods and evidence

Evidence review date: **2026-09-28**. This is a focused, source-verified methods review, not a systematic literature search or clinical validation. It records the current implementation; it does not authorize a production method change.

## Navigate

- [Literature review](literature_review.md): thirteen method questions, implementation alignment and comparator protocol.
- [Evidence matrix](evidence_matrix.csv): source type, findings, limitations, implementation status and decision links.
- [Method decisions](methods_decisions.md): NV-METHOD-001 through NV-METHOD-012, with evidence required to revisit each decision.
- [Claim map and terminology lock](claim_map.md): permitted language and limits.
- [BibTeX library](references.bib): stable citation keys used throughout this folder. Linked `@keys` in Markdown identify these entries without requiring a citation renderer.

**External evidence** establishes methodological precedent or limitations. **Our synthetic validation** evaluates this implementation under declared conditions. Neither automatically validates patient measurements. `Production` means used by the current case pipeline, not clinically approved. `Experimental` means implemented but not a production default; `validation-only` refers to experiment/evaluation machinery; `future work` means proposed and not implemented. Alignment labels describe correspondence with a source, independently of deployment status.

## Internal evidence index

These are pointers to the recorded reports, not new experiments. Denominators include the documented failures; paired trims are not independent replicates.

| Our synthetic validation | Evidence relevant to decisions | Boundary | Machine-readable record |
|---|---|---|---|
| [Caliber v1](../validation/caliber_validation_v1.md) | Raw EDT mean run MAE 0.201283 mm across 52 runs; centered 1 mm isotropic trimmed cylinder has 41.4214% error | Binary voxelization, two spacings; no universal correction or trimming benefit | [summary](../../validation/public_results/caliber_v1/summary.json) |
| [Caliber v2](../validation/caliber_validation_v2.md) | Cross-sectional wins 78, EDT wins 8, ties 6 among 92 scorable pairs; 24 upstream-failed pairs retained | Common recovered positions; 0.01 mm tie rule; no patient accuracy inference | [summary](../../validation/public_results/caliber_v2/summary.json), [pairs](../../validation/public_results/caliber_v2/matched_comparison.csv) |
| [Centerline v1](../validation/centerline_validation_v1.md) | Lee extracts 48/60, geodesic 60/60; Lee has lower position/tangent error on shared successes | Recovery is not accuracy; tree cannot preserve cycles | [summary](../../validation/public_results/centerline_v1/summary.json) |
| [Centerline v2](../validation/centerline_validation_v2.md) | Hybrid retains 48 Lee paths and recovers 12 failures in the 60-case core; failure-inclusive coverage 71.77% → 91.73% | Four added loop tests do not establish a cycle-capable fallback | [summary](../../validation/public_results/centerline_v2/summary.json) |
| [Topology v1](../validation/topology_validation_v1.md) | New method has correct rank 76/76 but complete connectivity 58/76; figure-eight connectivity 0/6 | Rank is insufficient; junction/branch and length failures remain | [summary](../../validation/public_results/topology_v1/summary.json) |
| [Topology v2](../validation/topology_validation_v2.md) | Core length error 7.198325 → 2.770747 mm, coverage 97.6076% → 90.8395%; connectivity stays 58/76 | Geometry-only refinement preserves wrong logical structure too | [summary](../../validation/public_results/topology_v2/summary.json) |

## Source verification and provenance

Titles, authors and dates were checked against publisher/conference records, PubMed, author-institution records and official software documentation. Source access was sometimes limited to an abstract or indexed publisher text; this review does not claim full-text appraisal of every paper. No secondary commercial summary is used as methodological evidence. Reviews are identified as reviews, and software capabilities are not clinical evidence.

| Citation keys | Verification record and access scope |
|---|---|
| `yang2026topcow` | [PubMed 42631290](https://pubmed.ncbi.nlm.nih.gov/42631290/) and [author institution](https://researchinformation.umcutrecht.nl/en/publications/the-topcow-challenge-topology-aware-circle-of-willis-segmentation/): NEJM AI 3(8), 2026, DOI verified; publisher full text was not directly accessible. |
| `yang2023topcow` | [arXiv v1](https://arxiv.org/abs/2312.17670v1), initial 2023 preprint, also linked by the [official 2023 challenge](https://topcow23.grand-challenge.org/results/). Distinct version from the expanded 2026 paper, not independent corroboration. |
| `shit2021cldice` | [CVF paper](https://openaccess.thecvf.com/content/CVPR2021/html/Shit_clDice_-_A_Novel_Topology-Preserving_Loss_Function_for_Tubular_Structure_CVPR_2021_paper.html) for authors/title; [IEEE](https://ieeexplore.ieee.org/document/9578225/footnotes) for DOI and 2021 proceedings. Pagination omitted rather than reconciling differing index page ranges. |
| `ayre2022cow`, `bernaus2025tortuosity`, `samuels2000stenosis` | PubMed [34936097](https://pubmed.ncbi.nlm.nih.gov/34936097/), [40117796](https://pubmed.ncbi.nlm.nih.gov/40117796/), [10782772](https://pubmed.ncbi.nlm.nih.gov/10782772/): metadata and abstracts. Anatomy review is **Ayre et al., 2022**, online 2021, not Jones et al. No DOI was verified for Samuels; PMID/URL used. |
| `antiga2003geometry`, `antiga2004decomposition`, `piccinelli2009geometry` | PubMed [12846436](https://pubmed.ncbi.nlm.nih.gov/12846436/), [15191145](https://pubmed.ncbi.nlm.nih.gov/15191145/), [Piccinelli author-hosted paper](https://user.engineering.uiowa.edu/~raghavan/images/chung/04915792.pdf), [Antiga publication list](https://lantiga.github.io/publications.html) and VMTK references. |
| `cosottini2024tof`, `lakhani2023tof` | [Springer paper](https://link.springer.com/article/10.1186/s41747-024-00463-z), PubMed [38844683](https://pubmed.ncbi.nlm.nih.gov/38844683/) and [36173305](https://pubmed.ncbi.nlm.nih.gov/36173305/). Lakhani journal issue year is 2023; online date is 2022. |
| `lee1994thinning`, `bertrand1994simple` | Publisher records for [Lee](https://www.sciencedirect.com/science/article/abs/pii/S104996528471042X) and [Bertrand/Malandain](https://www.sciencedirect.com/science/article/pii/0167865594900469); metadata/abstracts, supplemented by scikit-image's Lee reference. |
| `vmtkCenterlines`, `vmtkGeometry`, `slicerExtract`, `slicerExtractSource`, `slicerCrossSection`, `scipyEDT`, `skimageSkeletonize` | Official documentation/source URLs in [references.bib](references.bib), accessed on the review date. Slicer source was inspected for branch extraction and merging. Mutable `master`/`stable` pages are **not pinned benchmark versions**; no publication year or release number is invented. |

The long TopCoW author lists are deliberately abbreviated in BibTeX using `and others`; the source records retain full authorship. The remaining entries list verified authors or a corporate software author. The original methods review cites 12 peer-reviewed articles/proceedings papers (including two systematic reviews), one preprint and seven official documentation/source records. The release adds the IXI website and official TopCoW external-dataset record. The reviews summarize external work; they are not new NeuroVasc validation.

Release update: [Data Sources](../../DATA_SOURCES.md) now documents TopCoW/IXI provenance, licenses and citations. Real-example identification is distinguished from synthetic validation. Original archives remain local; curated numerical evidence is published below.

## Original documentation-milestone verification (historical)

The creation check passed for all six files, 20 unique BibTeX keys, 20 evidence rows with the required 17 columns, all 12 decision IDs, all 13 review sections, cited-key coverage and repository-relative links. The only root-relative example link is explicitly part of the proposed README text below. Citation metadata/access qualifications are in the source audit above; missing software publication dates and the unverified Samuels DOI are intentionally omitted.

Only this new documentation directory was written. Readable-file SHA-256 comparisons found 4,911 existing files unchanged, and size/mtime comparisons found 5,151 unchanged records. An earlier partial baseline also matched. Reads of 241 existing regression PNGs timed out; those files received metadata checks only, so this is **not a complete byte-level archive audit**. No archive was regenerated, no production/case/frontend/backend file was edited and nothing was pushed. Dependency/build/cache directories and Finder metadata were excluded from the audit. This checkout has no Git metadata, so preservation checks did not use `git diff`. Application builds and algorithm tests were not rerun for documentation-only additions.

## Original suggested README wording

> NeuroVasc Workbench explores anatomically labelled cerebral vessel segmentations using 3D visualization and segmentation-derived caliber profiles. Anatomical labels follow the TopCoW convention; retaining network connectivity is important for this anatomy. [TopCoW](https://doi.org/10.1056/AIdbp2500994), [Ayre et al.](https://doi.org/10.1111/joa.13616).
>
> The current case pipeline uses Lee skeletonization and spacing-aware EDT caliber. Separate synthetic studies compare experimental cross-sectional equivalent diameter, centerline recovery and loop-aware topology methods. Their evidence and limitations are recorded in the [methods and evidence layer](README.md). Cross-sectional measurements have established software precedent in [SlicerVMTK](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/CrossSectionAnalysis.md); that precedent does not establish NeuroVasc's clinical accuracy.
>
> This is research software. Candidate caliber reduction is an internal geometric descriptor, not a clinical stenosis measurement or diagnosis. Synthetic validation does not establish clinical validation.

The block above is retained as original suggested wording. The public root README now contains the implemented methods, data and validation sections; within this folder use the navigation above. Bibliographic keys for its sources are `yang2026topcow`, `ayre2022cow` and `slicerCrossSection`.

## Public release audit

The original preservation notes above describe the documentation-only milestone, before Git initialization. For current release changes, validation and publication status, see [release audit](../release/release_audit.md). Dataset citations were rechecked against the [official release](https://zenodo.org/records/15692630), [IXI website](https://brain-development.org/ixi-dataset/) and [authors’ updated publication record](https://arxiv.org/abs/2312.17670v5). The standalone VMTK benchmark is now implemented; Slicer-based cross-sectional comparison remains unavailable.
