# Actual VMTK comparator: isolated runtime and interchange

This directory is benchmark infrastructure, not a production dependency. The experiment report is [External Methods Benchmark v1](../../docs/validation/external_methods_benchmark_v1.md). No clinical claims or production-method switch follow from these synthetic experiments.

## Runtime verified on this machine

Actual VMTK 1.5.1 compiled computational-geometry filters run under macOS arm64, CPython 3.14.6, VTK 9.6.2 and NumPy 2.5.3. The existing NeuroVasc `.venv` retains its own VTK 9.7.0. JSON and VTP files cross the process boundary; the two VTK installations are never imported into the same process.

The [official installation page](https://vmtk.github.io/download/) recommends PyPI. Conda is present on this machine, but was unnecessary. `requirements.txt` specifies the comparator versions; `requirements.lock` pins all 13 wheels actually installed, including SHA-256 hashes for **macOS arm64 / CPython 3.14**. This is a platform-specific lock, not a promise that these wheel hashes work on Linux or another Python ABI. `install-report.json` records wheel URLs and hashes from the successful pip installation. Runtime probes additionally hash the installed computational-geometry binary, library and relevant VMTK wrapper sources.

There is a real packaging limitation: the umbrella import `from vmtk import vtkvmtk` fails because the segmentation library requests `libITKLabelMap-5.4.1.dylib`. The adapter directly imports `vtkvmtkComputationalGeometryPython`, which loads and executes the actual VMTK centerline, branch, geometry and bifurcation-reference filters. It does not implement replacements, edit installed libraries, or pretend the full VMTK application works. `environment.json` preserves the full loader error. The general CLI/PypeS workflow has not been certified by this benchmark.

## Setup and run

From the repository root, when creating the optional environment on this same platform:

```bash
.venv/bin/python -m venv external/vmtk/.venv
external/vmtk/.venv/bin/python -m pip install --only-binary=:all: --require-hashes -r external/vmtk/requirements.lock
external/vmtk/.venv/bin/python -m external.vmtk.adapter --probe
```

Do not install VMTK into the production `.venv`. Keep the venv executable path intact; resolving its symlink to the base interpreter bypasses the environment. On another platform, create a separate environment and install `requirements.txt`, record that platform's actual wheel provenance, and rerun the probe and integration tests before making comparisons. The initial successful install used `pip install --only-binary=:all: --report /tmp/neurovasc-vmtk-install.json vmtk==1.5.1`; the report was copied here without changing it.

The staged workflow uses the production interpreter as an orchestrator:

```bash
# 1. Read archived phantoms, verify exact masks/affines, export surfaces and inputs.
.venv/bin/python scripts/run_external_methods_benchmark.py --stage export
# 2. Execute actual VMTK in its isolated interpreter.
.venv/bin/python scripts/run_external_methods_benchmark.py --stage run
# 3. Import provenance-bound outputs, evaluate and write figures/report data.
.venv/bin/python scripts/run_external_methods_benchmark.py --stage report
# Independent repeat of all 120 external jobs; primary runs remain intact.
.venv/bin/python scripts/verify_external_methods_repeatability.py
# Tests without optional VMTK, then optional real subprocess integration:
.venv/bin/python -m pytest -q tests/test_external_methods.py
.venv/bin/python -m pytest -q tests/test_external_methods_vmtk_integration.py
# Complete regression suite:
.venv/bin/python -m pytest -q
```

The equivalent single command is `.venv/bin/python scripts/run_external_methods_benchmark.py`, or `sh external/vmtk/run_vmtk_benchmark.sh`. `--external-python /absolute/path/to/python` selects another isolated runtime. `--reuse` skips existing jobs only when request bytes, surface bytes, VMTK modules and adapter sources have matching hashes. Repeatability verification uses fresh subprocesses, fixed VTK RNG seed 0 and one VTK SMP thread; the runner can schedule four independent jobs concurrently. A 120-second per-job timeout is recorded as failure.

Only `validation/outputs/external_methods_v1/` and descendants are accepted destinations. A development subset must use a descendant, e.g. `--output validation/outputs/external_methods_v1/trial --limit 1`. `--stage report` does not execute or synthesize external measurements. Missing runtime/results produce explicit unavailable rows and blank metrics. Integration tests skip explicitly if the runtime probe fails. Re-running replaces this milestone's generated files, never earlier archives.

## Protocol and data contracts

- Surface: unchanged NeuroVasc binary 0.5 marching cubes with the volume affine, then zero-tolerance VTK cleaning and triangulation. Existing binary caps remain. No decimation, smoothing, axis flip or clipping. `surface_parameters.json` records per-surface counts/bounds.
- Seeds: oracle analytic terminal points versus automatic unchanged raw-geodesic graph endpoints, sorted by `(z,x,y)`; first point is source, all others targets. Both are snapped to nearest surface vertices. Requested/snapped coordinates and vertex IDs are retained. Automatic mode is not Slicer endpoint auto-detection and is not independent of NeuroVasc's endpoint strategy.
- Centerlines: actual `vtkvmtkPolyDataCenterlines`, cost `1/R`, no endpoint append, no resampling, no Voronoi simplification, no normal flip, Delaunay tolerance 0.001. Voronoi and MIS arrays are saved. Nonfinite/missing arrays or paths of length ≤1e-7 mm are invalid.
- Branches: actual compiled branch extractor, branch geometry and bifurcation reference systems; raw GroupIds, CenterlineIds, TractIds, Blanking and geometry arrays remain available. Partial branch comparisons use nonblanked GroupIds and attachments through blanked junction groups, not a reduction of duplicated source-target polylines.
- Tortuosity: the installed VMTK branch-geometry output reports `L/chord - 1`. Raw values and an explicitly named `+1` conversion are saved beside NeuroVasc's `L/chord` Distance Metric. Slicer Extract Centerline documents `L/chord`; do not silently equate these outputs. Branch cut boundaries and weighting also differ.
- Cycles: the declared tree workflow has no whole-network loop contract. Eight cyclic masks have surfaces/manifests but `METHOD_ASSUMPTION_MISMATCH` results, not scored failures or invented paths. See the report's source-based scope investigation.

## SlicerVMTK CE: manual protocol, **not executed**

No verified 3D Slicer/SlicerVMTK runtime is available here. `cross_sections.py` deliberately reports unavailable, with a table schema; no Slicer CE values or version are invented. These instructions are a proposed manual execution protocol, not a claim of tested Slicer automation.

1. On a machine with 3D Slicer and the VMTK extension, record exact Slicer version/revision, extension Git revision/package version, OS and bundled VTK/VMTK versions. Save the scene, logs, parameters and untouched exported CSV under a new `external_methods_v1/slicer_manual/<run>/` directory.
2. Choose a **single-tube** successful job first, e.g. `04_straight_4_centered_z_z0.5/oracle`. Input files are `inputs/<geometry_id>/surface.vtp` and `vmtk_outputs/<geometry_id>/<mode>/centerlines.vtp`. Record SHA-256 hashes, geometry ID and oracle flag. Do not replace them with a smoothed/decimated surface or rerun Extract Centerline silently.
3. Import the surface and VMTK centerline as model nodes. The files contain NIfTI-affine world mm with no LPS flip. To avoid file-loader coordinate conversion, the following **Slicer Python console** sketch assigns the raw VTK coordinates directly to RAS model nodes (replace `path`/`name` for each file). Check several coordinates and all bounds against the saved JSON before analysis:

   ```python
   import vtk, slicer
   reader = vtk.vtkXMLPolyDataReader()
   reader.SetFileName(path)
   reader.Update()
   model = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', name)
   model.SetAndObservePolyData(reader.GetOutput())
   model.CreateDefaultDisplayNodes()
   print(model.GetPolyData().GetBounds())
   ```

4. In **Cross-section analysis**, select the centerline model and the surface model; select/create an output table and run the module's calculation. Request split RAS coordinate columns and keep the original centerline point order/origin. Export the actual output table with distance, MIS diameter, area, CE diameter and coordinates. Do not calculate local NeuroVasc CE and relabel it as Slicer output. Record invalid/missing sections with reasons; never fill them with zero or MIS.
5. Normalize column names to the schema below in a separate copy, retaining the raw CSV. Require finite coordinates/units, positive finite area for valid rows, `CE=2*sqrt(area/pi)` as an internal consistency check, and agreement of MIS with twice the supplied radius array where present. Coordinates must match original sample locations within a declared numeric tolerance, not a fitted alignment. If exporting LPS, explicitly convert x/y signs back once and record it. Do not compare distances alone when origins or path direction differ.
6. For a later evaluated CE comparison, evaluate NeuroVasc's cross-sectional method on these same physical positions/tangents, or declare and quantify differences before matching common analytic stations. The current benchmark does **not** automatically import unverified manual files. A later adapter must validate provenance/schema and rerun the common-support analysis. Y junctions and multiple contours need separate scope/QC; do not average them into the single-tube CE result.

Expected normalized CSV columns (blank numeric values for invalid sections):

```text
geometry_id,oracle_seed,sample_index,x_mm,y_mm,z_mm,distance_mm,area_mm2,ce_diameter_mm,mis_diameter_mm,slicer_version,extension_revision,coordinate_system
```

Add `status`, `failure_reason`, surface/centerline SHA-256 and parameter-file reference. `coordinate_system` must be explicit (`RAS_mm` after verified conversion). Store the schema without example measurement rows to avoid fabricated evidence.

The official [Cross-section analysis documentation](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/CrossSectionAnalysis.md) distinguishes MIS diameter from surface area and CE diameter. [Extract Centerline](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/Docs/ExtractCenterline.md) also separates approximate network extraction from accurate seeded tree extraction; its default preprocessing is not used in this benchmark.
