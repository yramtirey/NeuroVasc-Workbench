# NeuroVasc frontend

React, TypeScript, Vite, VTK.js and Recharts dashboard. See the [project README](../README.md) for data preparation, scientific scope and backend setup.

Using Node.js 24:

```bash
npm ci
npm run build
npm run dev -- --host localhost
```

Open http://localhost:5173 with FastAPI at http://127.0.0.1:8000. The UI loads local `case_001`; obtain data separately as described in [Data Sources](../DATA_SOURCES.md). No real-case data are included in the source repository.

Selecting a vessel updates metrics, highlighting and its caliber profile. Clicking a profile sample or “Show minimum caliber” places the existing 3D marker at that sample's physical coordinates. Reset view frames the complete anatomy.
