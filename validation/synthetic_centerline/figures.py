"""Separate Matplotlib figures; failures retained as missing values, not zeros."""
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import nibabel as nib
import numpy as np
from .truth import sample_truth
from .report import aggregate


def create_figures(output, rows, profiles, examples):
    folder = output / 'figures'
    folder.mkdir(exist_ok=True)
    def save(name):
        plt.tight_layout()
        plt.savefig(folder / f'{name}.png', dpi=150)
        plt.close()
    methods = ('lee', 'alternative')
    for name, metric, ylabel in [('A_success', 'success_rate', 'Extraction success fraction'),
                                 ('D_position', 'mean_position_error_mm', 'Mean position error (mm)'),
                                 ('E_length', 'mean_length_error_mm', 'Mean absolute length error (mm)'),
                                 ('F_tangent', 'mean_tangent_error_deg', 'Mean tangent error (degrees)'),
                                 ('G_coverage', 'mean_coverage', 'Mean path coverage fraction')]:
        plt.figure(figsize=(6, 4))
        plt.bar(methods, [aggregate([r for r in rows if r['estimator_name'] == m])[metric] for m in methods])
        plt.ylabel(ylabel)
        plt.title('All 60 geometries' if metric == 'success_rate' else 'Successful runs only; coverage and failures reported separately')
        save(name)
    for name, field, metric, ylabel in [('B_phase_success', 'phase_id', 'success_rate', 'Success fraction'),
                                       ('C_spacing_success', 'spacing_z_mm', 'success_rate', 'Success fraction'),
                                       ('H_diameter_position', 'diameter_mm', 'mean_position_error_mm', 'Position error (mm)'),
                                       ('I_phase_position', 'phase_id', 'mean_position_error_mm', 'Position error (mm)'),
                                       ('J_orientation_position', 'orientation_id', 'mean_position_error_mm', 'Position error (mm)')]:
        selected = [r for r in rows if r['geometry'] != 'bifurcation']
        if field == 'diameter_mm':
            selected = [r for r in selected if r['geometry'] == 'straight']
        if field == 'phase_id':
            selected = [r for r in selected if r['geometry'] == 'straight']
        if field == 'orientation_id':
            selected = [r for r in selected if r['geometry'] == 'oblique']
        categories = sorted(set(r[field] for r in selected))
        plt.figure(figsize=(7, 4.5))
        for m in methods:
            groups = [[r for r in selected if r[field] == c and r['estimator_name'] == m] for c in categories]
            values = [aggregate(g)[metric] for g in groups]
            plt.plot(range(len(categories)), [np.nan if v is None else v for v in values], 'o-', label=m)
            for i, (group, value) in enumerate(zip(groups, values)):
                if value is not None:
                    plt.annotate(f"{sum(r['extraction_success'] for r in group)}/{len(group)}", (i, value), xytext=(0, 8 if m == 'lee' else -14), textcoords='offset points', fontsize=8)
        plt.xticks(range(len(categories)), categories)
        plt.xlabel(field)
        plt.ylabel(ylabel)
        plt.title('Labels: successful/attempted geometries; missing errors stay blank')
        plt.legend()
        save(name)
    plt.figure(figsize=(8, 4.5))
    for eid, item in examples.items():
        vessel = item['vessel']
        if vessel.geometry == 'curved' and vessel.proximal_diameter_mm == 3:
            for method in methods:
                p = [s for s in profiles if s['geometry_id'] == eid and s['estimator_name'] == method and s['path_kind'] == 'main']
                p.sort(key=lambda s: s['analytic_position'])
                plt.plot([s['analytic_position'] * vessel.length_mm for s in p], [s['position_error_mm'] for s in p], label=f'{method}, z={vessel.volume.spacing[2]} mm')
    plt.xlabel('Projected analytic arc distance (mm)')
    plt.ylabel('Position error (mm)')
    plt.legend()
    save('K_curved_error')
    for eid, item in examples.items():
        vessel = item['vessel']
        diagnostic = item['known_failure']
        example = vessel.geometry == 'bifurcation' or (vessel.geometry in ('straight', 'oblique', 'curved') and vessel.proximal_diameter_mm == 3 and item['cohort'] == 'core' and vessel.volume.spacing[2] == .5)
        if not diagnostic and not example:
            continue
        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(projection='3d')
        if diagnostic:
            points = nib.affines.apply_affine(vessel.volume.affine, np.argwhere(vessel.volume.data > .5))
            ax.scatter(*points.T, s=8, color='gray', alpha=.3, label='binary foreground voxels')
        names = list(vessel.branches) if vessel.geometry == 'bifurcation' else [None]
        for i, name in enumerate(names):
            truth = sample_truth(vessel, name)
            ax.plot(*truth.T, '--', color='C2', linewidth=2, label='analytic truth' if i == 0 else None)
        for index, method in enumerate(methods):
            graph = item[method][0]
            if graph is None:
                ax.plot([], [], [], color=f'C{index}', label=f'{method}: extraction failed')
            else:
                for i, (a, b) in enumerate(graph.edges):
                    segment = np.array([graph.nodes[a]['world'], graph.nodes[b]['world']])
                    ax.plot(*segment.T, color=f'C{index}', label=method if i == 0 else None)
        ax.set(xlabel='x (mm)', ylabel='y (mm)', zlabel='z (mm)')
        ax.set_title(eid, fontsize=10)
        # Equal physical scale; nonzero lateral width remains visible.
        bounds = points if diagnostic else np.concatenate([sample_truth(vessel, n) for n in names])
        ranges = np.maximum(np.ptp(bounds, axis=0) * 1.1, 8)
        center = (bounds.min(axis=0) + bounds.max(axis=0)) / 2
        ax.set_xlim(center[0] - ranges[0] / 2, center[0] + ranges[0] / 2)
        ax.set_ylim(center[1] - ranges[1] / 2, center[1] + ranges[1] / 2)
        ax.set_zlim(center[2] - ranges[2] / 2, center[2] + ranges[2] / 2)
        ax.set_box_aspect(ranges)
        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.zaxis.set_major_locator(MaxNLocator(5))
        ax.legend(fontsize=8)
        save(('failure_' if diagnostic else 'M_' if vessel.geometry == 'bifurcation' else 'L_') + eid)
