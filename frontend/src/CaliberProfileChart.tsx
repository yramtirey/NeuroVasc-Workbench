import { useEffect, useState } from "react";
import {
  CartesianGrid, Line, LineChart, ReferenceDot, ReferenceLine,
  ResponsiveContainer, Tooltip, XAxis, YAxis,
} from "recharts";
import type { MouseHandlerDataParam } from "recharts";
import { getVesselProfile } from "./api";
import type { CaliberSample, VesselProfile } from "./types";

interface CaliberProfileChartProps {
  caseId: string;
  label: number;
  onPointSelect?: (sample: CaliberSample) => void;
}

function CaliberProfileChart({ caseId, label, onPointSelect }: CaliberProfileChartProps) {
  const requestKey = JSON.stringify([caseId, label]);
  const [result, setResult] = useState<{
    key: string;
    profile?: VesselProfile;
    error?: string;
  } | null>(null);
  const [selection, setSelection] = useState<{ key: string; sample: CaliberSample } | null>(null);

  useEffect(() => {
    const controller = new AbortController();
    getVesselProfile(caseId, label, controller.signal)
      .then((profile) => {
        if (controller.signal.aborted) return;
        if (profile.case_id !== caseId || profile.label !== label) {
          throw new Error("The returned profile does not match the selected vessel.");
        }
        setSelection(null);
        setResult({ key: requestKey, profile });
      })
      .catch((error: unknown) => {
        if (controller.signal.aborted) return;
        setResult({ key: requestKey, error: error instanceof Error ? error.message : String(error) });
      });
    return () => controller.abort();
  }, [caseId, label, requestKey]);

  // Hide the previous vessel immediately, including before the fetching effect runs.
  const current = result?.key === requestKey ? result : null;
  const profile = current?.profile;
  const selectedSample = selection?.key === requestKey ? selection.sample : null;

  function selectSample(sample: CaliberSample) {
    setSelection({ key: requestKey, sample });
    onPointSelect?.(sample);
  }

  function showMinimum() {
    if (!profile?.samples.length) return;
    const sample = profile.samples.reduce((nearest, candidate) =>
      Math.abs(candidate.distance_mm - profile.minimum_location_mm)
        < Math.abs(nearest.distance_mm - profile.minimum_location_mm) ? candidate : nearest,
    );
    selectSample(sample);
  }

  function handleChartClick(state: MouseHandlerDataParam) {
    if (!state.isTooltipActive || state.activeTooltipIndex == null) return;
    const index = Number(state.activeTooltipIndex);
    if (!Number.isInteger(index)) return;
    const sample = profile?.samples[index];
    if (!sample) return;
    selectSample(sample);
  }

  if (current?.error) {
    return <div className="profile-status error" role="alert">{current.error}</div>;
  }
  if (!profile) {
    return <div className="profile-status" role="status">Loading caliber profile...</div>;
  }
  if (!profile.samples.length) {
    return <div className="profile-status" role="status">No caliber samples available for {profile.name}.</div>;
  }

  return (
    <div className="profile-chart-wrapper">
      <div className="profile-chart-title">
        <div>
          <strong>{profile.name}</strong>
          <span>Diameter along ordered centerline</span>
        </div>
        <div className="profile-minimum">
          <span>Minimum at {profile.minimum_location_mm.toFixed(1)} mm</span>
          <strong>{profile.minimum_diameter_mm.toFixed(2)} mm</strong>
          <button type="button" className="profile-minimum-button" onClick={showMinimum}>
            Show minimum caliber
          </button>
        </div>
      </div>
      <div className="profile-legend" aria-label="Chart legend">
        <span><i className="profile-key raw" />Raw diameter</span>
        <span><i className="profile-key smooth" />Smoothed diameter</span>
        <span><i className="profile-key reference" />Reference caliber</span>
        <span><i className="profile-key minimum" />Minimum caliber</span>
      </div>
      <ResponsiveContainer width="100%" height={280} minWidth={0}>
        <LineChart
          key={requestKey}
          data={profile.samples}
          onClick={handleChartClick}
          accessibilityLayer
          margin={{ top: 24, right: 30, bottom: 22, left: 12 }}
        >
          <CartesianGrid stroke="#e4ebeb" strokeDasharray="3 3" />
          <XAxis
            dataKey="distance_mm" type="number" domain={["dataMin", "dataMax"]}
            tickFormatter={(value: number) => value.toFixed(1)}
            label={{ value: "Distance along centerline (mm)", position: "insideBottom", offset: -16 }}
          />
          <YAxis
            domain={[0, "auto"]}
            tickFormatter={(value: number) => value.toFixed(1)}
            label={{ value: "Diameter (mm)", angle: -90, position: "insideLeft" }}
          />
          <Tooltip
            formatter={(value, name) => [`${Number(value).toFixed(2)} mm`, String(name)]}
            labelFormatter={(value) => `${Number(value).toFixed(1)} mm along vessel`}
          />
          <ReferenceLine
            y={profile.reference_diameter_mm} stroke="#71888b" strokeDasharray="6 4"
            ifOverflow="extendDomain"
          />
          <Line
            type="linear" dataKey="raw_diameter_mm" name="Raw diameter"
            stroke="#9eb2b4" strokeWidth={1.5} dot={false} activeDot={{ r: 4 }}
            isAnimationActive={false}
          />
          <Line
            type="linear" dataKey="smooth_diameter_mm" name="Smoothed diameter"
            stroke="#27666a" strokeWidth={2.5} dot={false} activeDot={{ r: 5 }}
            isAnimationActive={false}
          />
          <ReferenceDot
            x={profile.minimum_location_mm} y={profile.minimum_diameter_mm}
            r={6} fill="#a84c3d" stroke="#ffffff" strokeWidth={2}
            ifOverflow="extendDomain" label={{ value: "Minimum", position: "top" }}
          />
          {selectedSample && (
            <ReferenceDot
              x={selectedSample.distance_mm} y={selectedSample.smooth_diameter_mm}
              r={7} fill="none" stroke="#27666a" strokeWidth={2}
            />
          )}
        </LineChart>
      </ResponsiveContainer>
      <div className="profile-selection" aria-live="polite">
        {selectedSample ? (
          <>
            <strong>Selected location</strong>
            <span>{selectedSample.distance_mm.toFixed(2)} mm along vessel</span>
            <span>{selectedSample.smooth_diameter_mm.toFixed(2)} mm diameter (smoothed)</span>
            <span className="profile-raw-readout">Raw {selectedSample.raw_diameter_mm.toFixed(2)} mm</span>
          </>
        ) : "Hover to inspect caliber; click to select a sample."}
      </div>
      <div className="profile-note">{profile.note}</div>
    </div>
  );
}

export default CaliberProfileChart;
