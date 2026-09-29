export interface VesselMetrics {
  label: number;
  abbreviation: string;
  name: string;

  voxel_count: number;
  voxel_volume_mm3: number;
  mesh_volume_mm3: number;
  surface_area_mm2: number;

  centerline_components: number;
  main_path_length_mm: number;

  mean_diameter_mm: number;
  median_diameter_mm: number;
  minimum_diameter_mm: number;
  maximum_diameter_mm: number;

  reference_diameter_mm: number;
  candidate_reduction_percent: number;

  minimum_location_mm: number;

  minimum_x_mm: number;
  minimum_y_mm: number;
  minimum_z_mm: number;
}

export interface CaseReport {
  case_id: string;
  source_file: string;

  shape: number[];
  spacing_mm: number[];

  labels_present: number[];
  vessel_count: number;

  vessels: VesselMetrics[];

  warnings: string[];

  measurement_note: string;
}

export interface CaliberSample {
  distance_mm: number;
  raw_diameter_mm: number;
  smooth_diameter_mm: number;

  x_mm: number;
  y_mm: number;
  z_mm: number;
}


export interface VesselProfile {
  case_id: string;

  label: number;

  abbreviation: string;
  name: string;

  reference_diameter_mm: number;

  minimum_diameter_mm: number;
  minimum_location_mm: number;

  candidate_reduction_percent: number;

  samples: CaliberSample[];

  note: string;
}