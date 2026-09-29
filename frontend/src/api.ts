import type {
  CaseReport,
  VesselProfile,
} from "./types";


const API_BASE = "http://127.0.0.1:8000";

export function getGeometryUrl(caseId: string, label: number): string {
  return `${API_BASE}/api/cases/${encodeURIComponent(caseId)}/geometry/${label}`;
}

export async function getCase(
  caseId: string,
  signal?: AbortSignal,
): Promise<CaseReport> {
  const response = await fetch(
    `${API_BASE}/api/cases/${encodeURIComponent(caseId)}`,
    { signal },
  );

  if (!response.ok) {
    throw new Error(
      `Unable to load case ${caseId}`,
    );
  }

  return response.json();
}

export async function getVesselProfile(
  caseId: string,
  label: number,
  signal?: AbortSignal,
): Promise<VesselProfile> {

  const response = await fetch(
    `${API_BASE}/api/cases/${encodeURIComponent(caseId)}/profiles/${label}`,
    { signal },
  );

  if (!response.ok) {

    throw new Error(
      `Unable to load vessel profile ${label} (HTTP ${response.status})`,
    );
  }

  return response.json();
}
