"""Contract for real SlicerVMTK output; never substitutes local CE values."""
REQUIRED_COLUMNS = ('geometry_id', 'oracle_seed', 'sample_index', 'x_mm', 'y_mm', 'z_mm',
                    'distance_mm', 'area_mm2', 'ce_diameter_mm', 'mis_diameter_mm',
                    'slicer_version', 'extension_revision', 'coordinate_system')


def status():
    return {'status': 'unavailable', 'executed': False, 'version': None,
            'reason': 'No verified 3D Slicer + SlicerVMTK CrossSectionAnalysis runtime; no external CE executed',
            'expected_output_columns': list(REQUIRED_COLUMNS)}
