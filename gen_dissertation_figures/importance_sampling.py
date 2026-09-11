# Figure description: This eight-panel figure isolates the effect of the flux-based spatial sampling used to construct the LiuqeNN dataset. Each row shows one traceable equilibrium from the raw LIUQE dataset, with its poloidal-flux map and approximate last closed flux surface held fixed across the row. The left column draws points uniformly inside the TCV first wall. The remaining columns apply progressively narrower flux-importance sampling with $\sigma_\psi=1.0$, the dataset's production value $\sigma_\psi=0.3$, and $\sigma_\psi=0.1$. Smaller widths concentrate samples more strongly near the LCFS according to their normalized flux distance from it, while larger widths spread samples more broadly through the vessel. Every panel contains 800 samples, and fixed random seeds make the comparison reproducible. Black foreground markers identify the retained spatial samples; the dashed white contour identifies the approximate LCFS estimated from the boundary of the nonzero-current support.
from pathlib import Path
from matplotlib.path import Path as MplPath

from ._common import *
from ._importance_sampling import sample_flux_importance, sample_uniform_vessel
from ._raw_data import load_raw_samples


RAW_DIR = Path('raw_ds')
SELECTIONS = ((86168, 152), (83825, 184))
N_POINTS, IMPORTANCE_SIGMAS, RANDOM_SEED = 800, (1, 0.3, 0.1), 31
FIGURE_HEIGHT = 5.0  # [in]
MAP_LEVELS, MAP_FILL_ALPHA, MAP_LINE_ALPHA = 16, 0.75, 0.9
MARKER_SIZE = 8.8


def _lcfs_level(flux, current):
    rq, zq = find_approx_lcfs(current)
    if not rq.size: raise ValueError('Cannot estimate the LCFS from an empty current support.')
    points = np.stack((rq, zq), axis=-1)
    return float(np.mean(interp_pts(flux, points, method='linear'))), rq, zq


def _sample_columns(flux, current, seed):
    lcfs_flux, rq, zq = _lcfs_level(flux, current)
    uniform = sample_uniform_vessel(N_POINTS, rng=np.random.default_rng(seed))
    importance = [sample_flux_importance(N_POINTS, flux, rq, zq, sigma=sigma, rng=np.random.default_rng(seed+i+1)) for i, sigma in enumerate(IMPORTANCE_SIGMAS)]
    return (uniform, *importance), lcfs_flux


def _plot_lcfs(flux, current, lcfs_flux):
    contours = plt.contour(RRD, ZZD, flux, levels=[lcfs_flux], linewidths=0)
    segments = contours.allsegs[0]; contours.remove()
    weights = np.abs(current); center = (np.average(RRD, weights=weights), np.average(ZZD, weights=weights))
    enclosing = [segment for segment in segments if len(segment) > 2 and MplPath(segment, closed=True).contains_point(center)]
    if not enclosing: raise ValueError('No LCFS contour encloses the plasma-current centroid.')
    lcfs = max(enclosing, key=lambda segment: len(segment))
    plt.plot(lcfs[:,0], lcfs[:,1], color='white', linestyle='--', linewidth=1.0, zorder=4)


def _plot_panel(slot, panel, sample, points, lcfs_flux, title):
    plt.subplot(slot); flux = sample.fg[0]
    field = np.ma.masked_where(~IV, flux)
    cmap = plt.rcParams['image.cmap'] + ('_r' if np.mean(flux[IV]) < 0 else '')
    contours = plt.contourf(RRD, ZZD, field, levels=MAP_LEVELS, cmap=cmap, alpha=MAP_FILL_ALPHA)
    plt.contour(RRD, ZZD, field, levels=contours.levels, cmap=cmap, linewidths=0.4, alpha=MAP_LINE_ALPHA)
    _plot_lcfs(flux, sample.fg[1], lcfs_flux)
    plt.scatter(points[:,0], points[:,1], color='black', s=MARKER_SIZE, marker='.', linewidths=0, rasterized=False, zorder=5)
    prep_plot_tcv(plt.gca(), lw=0.7, alpha=0.8, c='k')
    plt.title(rf'({chr(97+panel)}) {title}')
    if panel % 4 == 0: plt.ylabel(rf'TCV \#{sample.shot}' + '\n' + rf'$t={sample.time:.3f}\,$s')


def generate(selections=SELECTIONS, raw_dir=RAW_DIR):
    """Generate the uniform and flux-importance sampling comparison."""
    if len(selections) != 2: raise ValueError(f'Exactly two equilibria are required, got {len(selections)}.')
    samples = load_raw_samples(raw_dir, selections)
    titles = ('Uniform', *(rf'$\sigma_\psi={sigma:g}$' for sigma in IMPORTANCE_SIGMAS))
    plt.figure(figsize=(THESIS_W, FIGURE_HEIGHT))
    for row, sample in enumerate(samples):
        points_by_column, lcfs_flux = _sample_columns(sample.fg[0], sample.fg[1], RANDOM_SEED+10*row)
        for col, (points, title) in enumerate(zip(points_by_column, titles)):
            panel = 4*row+col; _plot_panel(2*100+4*10+panel+1, panel, sample, points, lcfs_flux, title)
    save_figure('importance_sampling.pdf')


if __name__ == '__main__': generate(); compile_inspection_pdf()
