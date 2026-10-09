import {
  useEffect,
  useState,
} from "react";
import CaliberProfileChart from "./CaliberProfileChart";
import {
  Activity,
  Brain,
} from "lucide-react";

import { getCase, getCases } from "./api";

import type {
  CaseReport,
  VesselMetrics,
  CaliberSample,
} from "./types";

import "./App.css";
import VesselViewer from "./VesselViewer";

function MetricCard({
  label,
  value,
  subtext,
}: {
  label: string;
  value: string;
  subtext?: string;
}) {
  return (
    <div className="metric-card">
      <span className="metric-label">
        {label}
      </span>

      <strong className="metric-value">
        {value}
      </strong>

      {subtext && (
        <span className="metric-subtext">
          {subtext}
        </span>
      )}
    </div>
  );
}


function App() {
  const [
    report,
    setReport,
  ] = useState<CaseReport | null>(null);

  const [
    selected,
    setSelected,
  ] = useState<VesselMetrics | null>(null);

  const [
    error,
    setError,
  ] = useState<string | null>(null);

  const [
    selectedPoint,
    setSelectedPoint,
  ] = useState<CaliberSample | null>(
    null,
  );
  const [cases, setCases] = useState<{case_id: string}[]>([]);
  const [caseId, setCaseId] = useState<string | null>(new URLSearchParams(window.location.search).get("case"));
  useEffect(() => {
    const controller = new AbortController();
    getCases(controller.signal).then((items) => {
      if (controller.signal.aborted) return;
      if (!items.length && !new URLSearchParams(window.location.search).has("case")) {
        setError("No prepared cases. Run python scripts/prepare_synthetic_demo.py first.");
        return;
      }
      setCases(items);
      setCaseId((requested) => requested ?? (items.find((item) => item.case_id === "synthetic_demo")?.case_id ?? items[0]?.case_id ?? "case_001"));
    }).catch((err) => { if (!controller.signal.aborted) setError(String(err)); });
    return () => controller.abort();
  }, []);
  useEffect(() => {
    if (caseId === null) return;
    const controller = new AbortController();
    getCase(caseId, controller.signal)
      .then((data) => {
        if (controller.signal.aborted) return;
        if (!data.vessels.length) throw new Error("Case has no vessel metrics");
        setReport(data);
        setSelectedPoint(null);

        const initial =
          data.vessels.find(
            (vessel) =>
              vessel.abbreviation === "R-ICA",
          ) ?? data.vessels[0];

        setSelected(initial);
      })
      .catch((err) => {
        if (controller.signal.aborted) return;
        setError(String(err));
      });
    return () => controller.abort();
  }, [caseId]);


  if (error) {
    return (
      <div className="status-screen">
        Failed to load NeuroVasc:
        {error}
      </div>
    );
  }


  if (!report || report.case_id !== caseId || !selected) {
    return (
      <div className="status-screen">
        Loading NeuroVasc Workbench...
      </div>
    );
  }


  return (
    <div className="app">

      <header className="topbar">

        <div className="brand">
          <div className="brand-icon">
            <Brain size={23} />
          </div>

          <div>
            <h1>
              NeuroVasc Workbench
            </h1>

            <span>
              Cerebrovascular Imaging Analysis
            </span>
          </div>
        </div>

        <div className="case-summary">
          <label>
            Case <select aria-label="Case" value={report.case_id} onChange={(event) => {
              const next = event.target.value;
              const url = new URL(window.location.href);
              url.searchParams.set("case", next);
              window.history.replaceState(null, "", url);
              setSelectedPoint(null);
              setError(null);
              setCaseId(next);
            }}>
              {cases.map((item) => <option key={item.case_id} value={item.case_id}>{item.case_id}</option>)}
            </select>
          </label>

          <small>
            {report.vessel_count} vessels
          </small>
        </div>

      </header>


      {report.case_id === "synthetic_demo" && <p className="synthetic-banner" role="note">Synthetic demonstration — stylized generated tubes, not real anatomical imaging. No anatomical or clinical accuracy is claimed.</p>}
      <div className="workspace">

        <aside className="sidebar">

          <div className="sidebar-title">
            <Activity size={17} />

            <span>
              {report.case_id === "synthetic_demo" ? "Synthetic tubes" : "Anatomy"}
            </span>
          </div>

          <div className="vessel-list">

            {report.vessels.map(
              (vessel) => (

                <button
                  key={vessel.label}
                  className={
                    selected.label === vessel.label
                      ? "vessel-button selected"
                      : "vessel-button"
                  }
                  onClick={() => {
                    setSelected(vessel);
                    setSelectedPoint(null);
                  }}
                >

                  <div>
                    {vessel.abbreviation}
                  </div>

                  <small>
                    {vessel.name}
                  </small>

                </button>

              ),
            )}

          </div>

        </aside>


        <main className="content">

          <section className="vessel-heading">

            <div>
              <span className="eyebrow">
                Selected vessel
              </span>

              <h2>
                {selected.name}
              </h2>

              <span className="code-badge">
                {selected.abbreviation}
              </span>
            </div>

            <div className="coordinates">
              Minimum-caliber location

              <strong>
                {selected.minimum_location_mm.toFixed(1)}
                {" mm"}
              </strong>
            </div>

          </section>


          <section className="metrics">

            <MetricCard
              label="Main path length"
              value={
                `${selected.main_path_length_mm.toFixed(1)} mm`
              }
            />

            <MetricCard
              label="Mean caliber"
              value={
                `${selected.mean_diameter_mm.toFixed(2)} mm`
              }
            />

            <MetricCard
              label="Minimum caliber"
              value={
                `${selected.minimum_diameter_mm.toFixed(2)} mm`
              }
            />

            <MetricCard
              label="Geometric reduction"
              value={
                `${selected.candidate_reduction_percent.toFixed(1)}%`
              }
              subtext="Segmentation-derived"
            />

          </section>


          <div className="analysis-grid">

            <section className="panel viewer-panel">

              <div className="panel-header">
                <span>
                  {report.case_id === "synthetic_demo" ? "3D stylized tubes" : "3D Anatomy"}
                </span>

                <small>
                  {selected.abbreviation}
                </small>
              </div>

              <div className="viewer-container">

                <VesselViewer
                  key={report.case_id}
                  caseId={report.case_id}
                  selectedLabel={selected.label}
                  labels={report.labels_present}
                  markerPoint={selectedPoint}
                />

              </div>

            </section>


            <section className="panel details-panel">

              <div className="panel-header">
                Vessel Details
              </div>

              <div className="detail-row">
                <span>
                  Reference caliber
                </span>

                <strong>
                  {selected.reference_diameter_mm.toFixed(2)}
                  {" mm"}
                </strong>
              </div>

              <div className="detail-row">
                <span>
                  Maximum caliber
                </span>

                <strong>
                  {selected.maximum_diameter_mm.toFixed(2)}
                  {" mm"}
                </strong>
              </div>

              <div className="detail-row">
                <span>
                  Surface area
                </span>

                <strong>
                  {selected.surface_area_mm2.toFixed(1)}
                  {" mm²"}
                </strong>
              </div>

              <div className="detail-row">
                <span>
                  Segmented volume
                </span>

                <strong>
                  {selected.voxel_volume_mm3.toFixed(1)}
                  {" mm³"}
                </strong>
              </div>

              <div className="detail-row">
                <span>
                  Components
                </span>

                <strong>
                  {selected.centerline_components}
                </strong>
              </div>

            </section>

          </div>


          <section className="panel profile-panel">

            <div className="panel-header">
              Caliber Profile
            </div>

            <CaliberProfileChart
              caseId={report.case_id}
              label={selected.label}
              onPointSelect={setSelectedPoint}
            />

          </section>


          <footer className="measurement-note">
            {report.measurement_note}
          </footer>

        </main>

      </div>

    </div>
  );
}


export default App;
