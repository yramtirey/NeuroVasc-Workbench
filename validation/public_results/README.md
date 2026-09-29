# Public validation evidence

Curated snapshots from seven completed synthetic milestones: caliber v1/v2, centerline v1/v2, topology v1/v2, and external methods v1. Each folder includes its archived numerical summary and results table, plus 2–3 representative figures. `curation.json` records source/public hashes and transformations. CSVs/figures are byte-identical; JSON numeric values are unchanged, with only local path strings generalized. `topology_v1/mask_fingerprints.json` preserves exact uint8 mask-array hashes, shapes and archived affines without redistributing generated NIfTI files.

No real-case imaging or case outputs are included. See [validation reports](../../docs/validation) for full methods, failures, denominators and limitations. These are synthetic geometric experiments, not clinical validation.

## Regenerate complete archives

From the repository root, in the Python 3.14 environment with `.[validation]` installed, run in this dependency order:

```bash
python scripts/run_caliber_validation.py
python scripts/run_caliber_validation_v2.py
python scripts/run_centerline_validation.py
python scripts/run_centerline_validation_v2.py
python scripts/run_topology_validation.py
python scripts/run_topology_validation_v2.py
# Optional external environment: see external/vmtk/README.md first.
python scripts/run_external_methods_benchmark.py
```

These commands create the ignored `validation/outputs/` archives, including masks needed by later runners. Run in a fresh checkout to avoid rewriting an existing local archive. The existing scripts retain their output-directory protections. Tests use curated frozen CSVs when full archives are absent; the mask regression still checks exact archived array hashes and affine tolerances. Optional VMTK tests skip explicitly if its separate runtime is absent. Public evidence is not regenerated automatically during tests or installation.
