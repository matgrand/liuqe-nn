# Figure description: These figures show traceable examples of the final LiuqeNN dataset constructed from user-selected raw TCV shots and time indices. The four upper panels follow cell 12 of ``prepare_dataset.ipynb`` and display the 500 spatial targets retained for one training example: poloidal flux $\psi$, plasma-current distribution $I_y$, and the radial and vertical magnetic fields $B_r$ and $B_z$. The spatial coordinates are drawn with the same flux-based importance-sampling rule used during dataset preparation, concentrating samples near the approximate LCFS while retaining points throughout the vessel. The dashed LCFS overlay is the $\psi$ contour whose value is estimated by interpolating the raw flux at the boundary of the nonzero-current support and applying a positive correction equal to 5% of the full flux-map span. Flux and current are interpolated from their raw LIUQE grids; the magnetic fields are obtained by differentiating flux on a seven-times-finer grid and then interpolating to the sampled coordinates. The lower panel shows the 134 physical inputs after normalization with the training-dataset statistics, grouped into magnetic probes $B_m$, flux loops $F_f$, poloidal-field-coil currents $I_a$, flux-loop voltages $U_f$, and $rB_t$. The first selection is written to ``final_ds_example.pdf`` for the main narrative, while subsequent selections are written to numbered PDFs for appendix use.
from pathlib import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ._common import *
from ._importance_sampling import sample_flux_importance
from ._raw_data import load_raw_sample


RAW_DIR = Path('raw_ds')
NORMALIZATION_PATH = Path('dss/eval_ds.npz')
SELECTIONS = ((86168, 152), (83825, 184), (80265, 195))  # First selection is the main-text figure.
N_POINTS, IMPORTANCE_SIGMA, RANDOM_SEED = 500, 0.3, 19
FIGURE_HEIGHT = 4.7  # [in]
MAP_NAMES = (r'$\psi$ (Fx)', r'$j$ (Iy)', r'$B_r$', r'$B_z$')
MAP_SCALES = (1.0, 1e-3, 1.0, 1.0)
MAP_UNITS = ('Wb', 'kA', 'T', 'T')


def interpolate_targets(flux, current, points, k_finer=7):
    """Reproduce the notebook interpolation and fine-grid magnetic-field reconstruction."""
    rf = np.linspace(RRD[0,0], RRD[0,-1], (RRD.shape[1]-1)*k_finer+1)
    zf = np.linspace(ZZD[0,0], ZZD[-1,0], (RRD.shape[0]-1)*k_finer+1)
    rrf, zzf = np.meshgrid(rf, zf); fine_points = np.stack((rrf, zzf), axis=-1)
    fine_flux = interp_pts(flux, fine_points.reshape(-1, 2)).reshape(rrf.shape)
    fine_br, fine_bz = [value.reshape(rrf.shape) for value in meqBrBz(fine_flux, rr=rrf, zz=zzf)]
    psi = interp_pts(flux, points)
    iy = interp_pts(current, points)
    br = interp_pts(fine_br, points, gr=rf, gz=zf, method='linear')
    bz = interp_pts(fine_bz, points, gr=rf, gz=zf, method='linear')
    return psi, iy, br, bz


def prepare_final_sample(raw_dir, shot, time_index, normalization, n_points=N_POINTS, sigma=IMPORTANCE_SIGMA, seed=RANDOM_SEED):
    """Load one raw slice and construct its sampled, normalized final-dataset representation."""
    sample = load_raw_sample(Path(raw_dir)/f'{shot}.mat', time_index)
    flux, current = sample.fg[:2]; rq, zq = find_approx_lcfs(current)
    points = sample_flux_importance(n_points, flux, rq, zq, sigma=sigma, rng=np.random.default_rng(seed))
    targets = interpolate_targets(flux, current, points)
    mean, std = normalization
    if mean.shape != sample.x.shape or std.shape != sample.x.shape: raise ValueError(f'Normalization shape {mean.shape, std.shape} does not match inputs {sample.x.shape}.')
    if np.any(std <= 0): raise ValueError('All input standard deviations must be positive.')
    return sample, points, targets, (sample.x-mean)/std


def _colorbar(image, axis):
    color_axis = make_axes_locatable(axis).append_axes('right', size='12%', pad='3%')
    locator = color_axis.get_axes_locator()
    def colorbar_locator(ax, renderer):
        position = locator(ax, renderer)
        return position.shrunk(1, 0.9).anchored('C', position)
    colorbar_locator.get_subplotspec = locator.get_subplotspec
    color_axis.set_axes_locator(colorbar_locator)
    plt.colorbar(image, cax=color_axis); plt.sca(axis)


def _plot_map(slot, panel, points, values, name, scale, unit, flux, current):
    plt.subplot(slot); axis = plt.gca()
    cmap = plt.rcParams['image.cmap'] + ('_r' if np.mean(values) < 0 else '')
    image = plt.scatter(points[:,0], points[:,1], c=values*scale, s=4, marker='o', linewidths=0, cmap=cmap, rasterized=True)
    plot_approx_lcfs(flux, current, axis)
    prep_plot_tcv(axis, lw=0.7, alpha=0.8, c='k')
    plt.title(rf'({chr(96+panel)}) {name} [{unit}]'); _colorbar(image, axis)


def _plot_inputs(slot, values):
    plt.subplot(slot); channels = np.arange(1, len(values)+1)
    sizes = np.array([DS_SIZES[name][0] for name in INPUT_NAMES]); edges = np.concatenate(([0], np.cumsum(sizes)))
    centers = (edges[:-1]+edges[1:]+1)/2
    plt.bar(channels, values, color=C2)
    plt.xticks(centers, ('Bm', 'Ff', 'Ia', 'Uf', 'rBt')); plt.xlim(0.2, len(values)+0.8)
    plt.title(r'(e) Normalized physical inputs [-]')


def _save_name(number): return 'final_ds_example.pdf' if number == 1 else f'final_ds_example_{number}.pdf'


def generate(selections=SELECTIONS, raw_dir=RAW_DIR, normalization_path=NORMALIZATION_PATH):
    """Generate one final-dataset example per selected raw shot and time index."""
    if not selections: raise ValueError('Select at least one (shot, time index) pair.')
    with np.load(normalization_path) as data:
        if 'x_mean_std' not in data: raise KeyError(f'Missing x_mean_std in {normalization_path}.')
        normalization = np.asarray(data['x_mean_std'])
    for number, (shot, time_index) in enumerate(selections, start=1):
        sample, points, targets, normalized_inputs = prepare_final_sample(raw_dir, shot, time_index, normalization, seed=RANDOM_SEED+number-1)
        plt.figure(figsize=(THESIS_W, FIGURE_HEIGHT)); grid = plt.GridSpec(2, 1, height_ratios=(1.3, 0.7))
        map_grid = grid[0].subgridspec(1, 4); input_grid = grid[1].subgridspec(1, 3, width_ratios=(0.08, 1, 0.08))
        for panel, (values, name, scale, unit) in enumerate(zip(targets, MAP_NAMES, MAP_SCALES, MAP_UNITS), start=1): _plot_map(map_grid[0,panel-1], panel, points, values, name, scale, unit, sample.fg[0], sample.fg[1])
        _plot_inputs(input_grid[0,1], normalized_inputs)
        plt.suptitle(rf'TCV \#{sample.shot}, $t={sample.time:.3f}\,$s (index {sample.time_index})', y=0.965)
        save_figure(_save_name(number), rect=(0, 0, 1, 0.985))


if __name__ == '__main__': generate(); compile_inspection_pdf()
