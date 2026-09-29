"""Separate Matplotlib figures; failure coverage is shown alongside accuracy."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def create_figures(output: Path, rows, matched, samples):
    directory = output / 'figures'
    directory.mkdir(exist_ok=True)
    data, pairs, profiles = pd.DataFrame(rows), pd.DataFrame(matched), pd.DataFrame(samples)
    valid = data[data.valid_sample_count > 0]
    core = valid[(valid.cohort == 'core') & (valid.trim_mm == 2)]
    plt.rcParams.update({'font.size': 10, 'figure.figsize': (7, 4.8), 'savefig.dpi': 180})

    def save(name, xlabel, ylabel, title):
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(alpha=0.2)
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(directory / name)
        plt.close()

    for estimator in ('edt', 'cross_sectional'):
        plt.figure()
        for (geometry, sz), group in core[(core.estimator_name == estimator) & (core.geometry != 'tapered')].groupby(['geometry', 'spacing_z_mm']):
            group = group.sort_values('true_diameter_mm')
            plt.plot(group.true_diameter_mm, group.estimated_mean_mm, marker='o' if sz == 0.5 else 's', linestyle='-' if sz == 0.5 else '--', label=f'{geometry}, z={sz} mm')
        plt.plot([1, 5], [1, 5], ':', label='Identity y=x')
        save(f'{estimator}_vs_true.png', 'True diameter (mm)', 'Estimated mean diameter (mm)', f'{estimator}: core phantoms, trim=2 mm')

    scored = pairs[pairs.winner != 'unscorable']
    plt.figure()
    for geometry, group in scored.groupby('geometry'):
        plt.scatter(group.edt_mae, group.cross_sectional_mae, label=geometry, alpha=0.7)
    limit = max(scored.edt_mae.max(), scored.cross_sectional_mae.max()) * 1.05
    plt.plot([0, limit], [0, limit], ':', label='Equal MAE')
    save('mae_head_to_head.png', 'EDT MAE (mm)', 'Cross-sectional MAE (mm)', 'Matched valid samples · all scorable experiments')

    plt.figure()
    for cohort, group in scored.groupby('cohort'):
        means = group.groupby('geometry').cross_sectional_mae_minus_edt_mae.mean()
        plt.plot(means.index, means.values, 'o', label=cohort)
    plt.axhline(0, linestyle=':', label='Equal MAE')
    save('mae_delta_by_geometry.png', 'Geometry', 'Cross-sectional MAE minus EDT MAE (mm)', 'Negative values favor cross-sectional')

    plt.figure()
    for estimator, group in core.groupby('estimator_name'):
        means = group.groupby('geometry').percent_error.mean()
        plt.plot(means.index, means.values, 'o-', label=estimator)
    save('percent_error.png', 'Geometry', 'Mean absolute percentage error (%)', 'Core phantoms · trim=2 mm')

    plt.figure()
    for estimator, group in core.groupby('estimator_name'):
        means = group.groupby('spacing_z_mm').mae_mm.mean()
        plt.plot(['Isotropic', 'Case-like anisotropic'], means.values, 'o-', label=estimator)
    save('spacing_comparison.png', 'Voxel spacing', 'Equal-run mean MAE (mm)', 'Core phantoms · trim=2 mm')

    for geometry in ('tapered', 'curved'):
        for sz in (0.5, 0.8):
            diameter = 4 if geometry == 'tapered' else 3
            selected = core[(core.geometry == geometry) & (core.spacing_z_mm == sz) & (core.proximal_diameter_mm == diameter)]
            experiment_id = selected.iloc[0].experiment_id
            plt.figure()
            for estimator in ('edt', 'cross_sectional'):
                group = profiles[(profiles.experiment_id == experiment_id) & (profiles.trim_mm == 2) & (profiles.estimator_name == estimator)].sort_values('normalized_position')
                plt.plot(group.normalized_position, group.diameter_mm, label=estimator)
            plt.plot(group.normalized_position, group.true_diameter_mm, '--', label='Analytic truth')
            save(f'{geometry}_profile_z{sz}.png', 'Analytic normalized position', 'Diameter (mm)', f'{geometry.title()} · z={sz} mm · trim=2 mm')

    plt.figure()
    straight = data[(data.geometry == 'straight') & (data.trim_mm == 2)]
    labels = ['centered', 'shift_0.25', 'shift_0.5']
    for (estimator, sz), group in straight.groupby(['estimator_name', 'spacing_z_mm']):
        means, descriptions = [], []
        for phase in labels:
            phase_group = group[group.phase_id == phase]
            good = phase_group[phase_group.valid_sample_count > 0]
            means.append(good.mae_mm.mean() if len(good) else np.nan)
            descriptions.append(f'{len(good)}/{len(phase_group)}')
        plt.plot(np.arange(len(labels)), means, 'o-', label=f'{estimator}, z={sz}; valid runs ' + ', '.join(descriptions))
    plt.xticks(np.arange(len(labels)), labels)
    plt.xlim(-0.1, 2.1)
    plt.text(0.76, 0.44, 'No recovered centerline\nat half-voxel phase\n(not zero error)', transform=plt.gca().transAxes, ha='center', fontsize=8)
    save('grid_phase_sensitivity.png', 'Grid offset in voxel units (all axes)', 'Mean MAE over successful runs (mm)', 'Straight phantoms · trim=2 mm · coverage in legend')

    plt.figure()
    sections = data[(data.estimator_name == 'cross_sectional') & (data.sample_count > 0)]
    rates = sections.groupby('geometry')[['valid_sample_count', 'sample_count']].sum()
    plt.bar(rates.index, rates.valid_sample_count / rates.sample_count, label='Valid sections / attempted sections')
    plt.ylim(0, 1.08)
    plt.text(0.02, 0.98, 'Shared centerline failures excluded; see summary.json', transform=plt.gca().transAxes, va='top', fontsize=8)
    save('valid_fraction.png', 'Geometry', 'Valid fraction', 'Cross-sectional coverage conditional on recovered centerline')

    plt.figure()
    oblique = valid[(valid.geometry == 'oblique') & (valid.trim_mm == 2)]
    for (estimator, sz), group in oblique.groupby(['estimator_name', 'spacing_z_mm']):
        means = group.groupby('orientation_id').mae_mm.mean()
        plt.plot(means.index, means.values, 'o-', label=f'{estimator}, z={sz}')
    save('orientation_sensitivity.png', 'Physical direction (normalized)', 'Mean MAE across diameters (mm)', 'Oblique cylinders · trim=2 mm')
