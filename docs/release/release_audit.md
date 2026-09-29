# NeuroVasc Workbench 0.1.0 release preparation

Local checks completed 2026-09-28. **Prepared locally; GitHub publication and portfolio deployment are blocked, not completed.** The verified GitHub identity is `yramtirey`. No NeuroVasc remote/repository exists, GitHub CLI is unavailable, local Git has no stored HTTPS credential, and the connected GitHub tools cannot create a repository. No push or force-push was attempted.

## Changes and preservation

- Added the public README, MIT code license, Data Sources, version-0.1.0 changelog, exclusions, screenshot attribution and these release records.
- Curated seven synthetic milestones (including the executed external-method benchmark): summaries, concise result tables and 15 representative figures. All 30 copied evidence files preserve original numerical values; CSVs/figures are byte-identical. Complete archives remain local.
- Updated literature links to public evidence, added verified IXI/dataset citations and corrected stale statements about the now-executed external benchmark.
- Three regression-test files now read public frozen evidence. The 76-mask test verifies generated array hashes, shapes and affines without requiring distributed NIfTI files. Existing scientific assertions/tolerances remain intact.
- Fixed a real frontend integration regression: `App.tsx` imported the chart and stored the selected point but still rendered a placeholder and omitted `markerPoint`. It now mounts the existing chart and passes the point to the existing viewer. No VTK architecture or scientific algorithm was refactored.
- Set frontend package metadata to 0.1.0, documented frontend startup and added a NeuroVasc favicon to remove a browser 404.
- Added genuine dashboard/linked-marker screenshots. No fabricated image or hosted-demo URL is used.

All 22 geometry/graph source hashes recorded by the external benchmark still match. No scientific module, real-case measurement, NIfTI, mesh or archived validation output was intentionally rewritten. Ignored package metadata/caches and frontend dependency/build directories changed during installation and testing. A damaged prior `node_modules` was moved to an ignored backup before reinstalling from the lockfile. The existing Python environment was retained; tests used a fresh isolated Python 3.14 environment after imports in the original environment timed out.

## Validation results

| Check | Result |
|---|---|
| Working-tree pytest | **234 passed, 2 skipped**, 56 upstream deprecation warnings |
| Public source-only checkout pytest | **234 passed, 2 skipped**, same warnings; no full archives/raw data present |
| Optional compiled VMTK tests | Skipped because the separate runtime was unavailable to the probe; historical recorded benchmark results are not represented as a fresh run |
| Analysis dependencies and FastAPI import/start | Passed in isolated release environment; API version 0.1.0 |
| Frontend clean install | `npm ci` passed; initial `npm install` also completed |
| Frontend TypeScript/build/lint | `npm run build` and `npm run lint` passed |
| Local API | Cases/report plus all 12 geometry and 12 profile routes returned HTTP 200 |
| Browser | 12 vessels, R-ICA/ACom/PCA/PCom selection, profile clicks, minimum selection, rotation/zoom and reset checked; no uncaught errors or console errors in final run |
| Viewer lifecycle | 12 geometry requests before and after selection changes; one canvas retained |
| Public candidate | No flagged private paths, credential patterns, excluded data or files above 5 MB |
| Portfolio | Lint, TypeScript and static export passed with `npm run build -- --webpack` |
| Portfolio browser | Personal Projects card/detail route and screenshot work; no uncaught errors; résumé/CV HTTP 200; no horizontal overflow at 390 px |

The normal Next/Turbopack build was attempted but failed under the local execution environment's process/port restrictions. The supported webpack build produced the static export without source or deployment-config changes. GitHub Actions and deployment have not been run for these changes. This is not a claim that the unchanged CI command has passed remotely.

Warnings retained: Vite's approximately 1.81 MB JavaScript chunk warning; 56 upstream scikit-image/NumPy deprecations; six low-severity npm audit entries in the polyfill dependency chain (elliptic advisory GHSA-848j-6mx2-7j84, no supplied fix). No moderate/high/critical npm advisory was reported. Dependency upgrades or viewer refactors were not introduced merely to suppress these warnings.

## Public contents

Included: analysis/API/UI source, scripts, tests, configuration, complete literature/validation reports, curated synthetic evidence, actual attributed screenshots and reproducibility instructions. Excluded: virtual environments, node_modules, distributions/caches, package egg-info, real-case `outputs/`, raw/sample imaging, dataset archives, full `validation/outputs/`, regression mirrors and private release scratch files. No Git LFS was introduced.

The [README](../../README.md) documents installation, separate backend/frontend startup, local case preparation, architecture, features, methods, limitations and citations. [Data Sources](../../DATA_SOURCES.md) and the [provenance audit](provenance_audit.md) identify the exact archive member, official sources, both required citations and the code/data license boundary.

## Portfolio preparation

The verified existing repository is `yramtirey/yramtirey.github.io`; existing site URL: https://yramtirey.github.io. Work was prepared in an isolated local clone so the original checkout and its untracked deployment file remain untouched.

Files changed/created in that clone:

- `app/page.tsx`: Personal Projects section using the existing project card.
- `app/projects/neurovasc-workbench/page.tsx`: detail page, route `/projects/neurovasc-workbench/`.
- `public/projects/neurovasc/neurovasc-dashboard.png` and `neurovasc-linked-view.png`: actual captures.
- `public/projects/neurovasc/ATTRIBUTION.md`: dataset/image citation and license.

Card copy: “Interactive cerebrovascular imaging and quantitative geometry for Circle-of-Willis analysis. Connects labeled MRA segmentations with anatomy-aware 3D visualization, caliber profiles and spatially linked measurements, supported by synthetic validation.”

Provenance copy identifies TopCoW annotations derived from IXI MRA, states that raw imaging is not redistributed, and distinguishes independently generated synthetic phantoms. The page states research/engineering prototype, not diagnostic software, and no clinical validation.

**No GitHub source URL has been inserted:** `repositoryUrl` remains null until a public source repository exists and its push is verified. Source/Data Sources/methods links are intentionally pending, so this portfolio draft must not be deployed yet. There is no fabricated live-demo link. Existing deployment configuration and site visual identity are preserved.

## Remaining publication steps

1. Authenticate GitHub CLI for the verified account (or create the empty public `NeuroVasc-Workbench` repository through GitHub and provide its confirmed URL).
2. Create/push the audited source repository on `main`, with commit message `Release NeuroVasc Workbench v0.1.0`. Verify files, rendered docs, images and absence of raw data.
3. Set the portfolio's `repositoryUrl` to the verified public URL. Rebuild and verify its source/Data Sources/methods links.
4. Push the portfolio change (`Add NeuroVasc Workbench project`) to its existing `main`, after checking upstream changes. Verify the Pages workflow and deployed route. Never force-push.

The local handoff records actual commit hashes and precise workstation commands. No GitHub repository URL, pushed branch or deployment success is claimed before these steps complete.
