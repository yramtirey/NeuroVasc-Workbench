# Data Sources

## TopCoW / IXI MRA external dataset

NeuroVasc development uses selected anatomical Circle-of-Willis annotations released by the TopCoW challenge. The `MRA_IXI_HH` external testset contains 20 MRA cases selected from the IXI dataset's HH (Hammersmith Hospital) cohort. TopCoW supplies the derived anatomical annotations; the underlying MRA originates from IXI. This relationship is documented in the [official TopCoW release](https://zenodo.org/records/15692630).

The local `cow_case_001_seg.nii.gz`, analyzed as `case_001`, is an **example annotated MRA case from the TopCoW/IXI external test dataset**, not an acquisition by this project. Its compressed file bytes exactly match this member of the local official-release archive:

```text
MRA_IXI_HH/cow_seg_labelsTr/IXI012-HH-1211-MRA_LPS_DFmasked_ROI-braincase_full_LPS_Mask.nii.gz
SHA-256: 7ee1bb10c9c08d5a653b91e8fd18790726e0bf9240ef96b3c578bd61ef9d00b6
```

This verifies the working copy against a published dataset identifier; no personal identity or clinical history is inferred. The segmentation is a TopCoW-derived annotation rather than an annotation created by NeuroVasc.

Raw scans, NIfTI segmentations and case-derived VTP/JSON/CSV outputs are not distributed in this source repository. Obtain data from the official sources and keep it locally:

- [TopCoW Training Data and External Testsets, version 1](https://zenodo.org/records/15692630), TopCoW Challenge Organizers, 2025. Dataset DOI: [10.5281/zenodo.15692630](https://doi.org/10.5281/zenodo.15692630). Select `MRA_IXI_HH.zip`; follow its accompanying terms and documentation.
- [IXI Dataset — Brain Development](https://brain-development.org/ixi-dataset/), the original acquisition source.
- [TopCoW challenge](https://topcow24.grand-challenge.org/), for challenge context and label conventions.

## Licensing

**Original IXI license: Creative Commons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0).** The [IXI source](https://brain-development.org/ixi-dataset/) and TopCoW's external-testset description explicitly identify this license. See the [license](https://creativecommons.org/licenses/by-sa/3.0/).

The repository's MIT license covers NeuroVasc source code and original documentation; it does not relicense TopCoW/IXI data or third-party dependencies. Retain the original data terms and any accompanying TopCoW annotation conditions. Do not apply the separately described TopCoW training-set license to every external dataset. Dataset-derived screenshots have separate attribution and licensing in [assets/README.md](assets/README.md).

## Required citation

When using this derived release, cite **both the TopCoW challenge publication and the IXI dataset website**. The current verified publication is:

Yang, K., Musio, F., Ma, Y., et al. **The TopCoW Challenge — Topology-Aware Circle of Willis Segmentation for CT and MR Angiography.** *NEJM AI*. 2026;3(8). DOI: [10.1056/AIdbp2500994](https://doi.org/10.1056/AIdbp2500994).

Title, lead authors, journal reference and DOI are corroborated by the authors' [arXiv record, version 5](https://arxiv.org/abs/2312.17670v5). The 2025 dataset record cites the earlier challenge preprint; this is the same evolving work, not independent evidence. The bibliography retains the explicitly versioned original preprint as well.

IXI attribution: **“IXI data were obtained from the IXI Dataset, Brain Development: https://brain-development.org/ixi-dataset/.”** This is a website/data-source citation, not an invented journal article or DOI. The source itself asks users to acknowledge its website.

Machine-readable entries are in [references.bib](docs/literature/references.bib): `yang2026topcow`, `yang2023topcow`, `ixiDataset`, `topcowExternalRelease`.

## Synthetic validation data

All caliber, centerline and topology validation phantoms are generated analytically by NeuroVasc's validation code. They are **not TopCoW/IXI scans** and contain no patient data. Public summaries and selected figures are in [validation/public_results](validation/public_results). Synthetic geometric validation does not establish clinical validation.

## Repository policy

Raw medical images, NIfTI segmentations, downloaded dataset archives, complete local validation outputs and real-case derivatives are excluded from Git. No dataset is downloaded automatically. Public evidence is limited to source, documentation, curated synthetic results and explicitly attributed genuine screenshots. The original complete local archives are preserved.

The release audit verifies dataset provenance using official sources and exact local-archive member hashing; no personal identity is inferred. Source pages were checked on 2026-09-28.
