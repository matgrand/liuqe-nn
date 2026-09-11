# Figure description: These figures show traceable snapshots of the raw LIUQE dataset before spatial point sampling and neural-network training. Each figure corresponds to one user-selected TCV discharge and the middle time slice that passes the dataset filters, unless an explicit time-index override is supplied; the shot and time provenance appear in the figure suptitle. The four-panel upper row follows cell 7 of ``prepare_dataset.ipynb``: it shows the stored poloidal-flux map $\psi$ and plasma-current distribution $I_y$, together with the radial and vertical magnetic fields $B_r$ and $B_z$ reconstructed from the flux map by finite differences. Each field is rendered with translucent filled contours and matching contour lines inside the vessel, suppressing the visual grid while retaining the raw gridded values and their spatial structure. The dashed LCFS overlay is the $\psi$ contour whose value is estimated by interpolating the flux at the boundary of the nonzero-current support and applying a positive correction equal to 5% of the full flux-map span. The three wider panels in the lower row show the stored magnetic-probe measurements $B_m$, flux-loop measurements $F_f$, and fitted poloidal-field-coil currents $I_a$. The first selected shot is written to ``raw_ds_example.pdf`` for the main narrative, while every subsequent selection is written to a numbered PDF for appendix use.
from pathlib import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ._common import *
from ._raw_data import _raw_array, raw_valid_mask


RAW_DIR = Path('raw_ds')
SHOTS = (86168, 83825, 80265)  # First shot is the main-text figure; the others are appendix figures.
TIME_INDICES = {}  # Optional overrides, for example: {86168: 120}.
FIGURE_HEIGHT = 4.7  # [in]
MAP_LEVELS, MAP_FILL_ALPHA, MAP_LINE_ALPHA = 16, 0.9, 1
MAP_NAMES = (r'$\psi$'+' (Fx)', r'$j$'+' (Iy)', r'$B_r$', r'$B_z$')
MAP_SCALES = (1.0, 1e-3, 1.0, 1.0)
MAP_UNITS = ('Wb', 'kA', 'T', 'T')


def _load_snapshot(raw_dir, shot, time_indices):
    path = Path(raw_dir)/f'{shot}.mat'
    if not path.is_file(): raise FileNotFoundError(f'Raw shot file not found: {path}')
    data = loadmat(path, variable_names=['t', IP, BM, FF, IA, FX, IY])
    nt = np.asarray(data['t']).size; valid = np.flatnonzero(raw_valid_mask(data))
    if not len(valid): raise ValueError(f'Shot {shot} has no time slices that pass the dataset filters.')
    index = time_indices.get(shot, valid[len(valid)//2])
    if index < 0: index += nt
    if index not in valid: raise ValueError(f'Shot {shot}, time index {index} does not pass the dataset filters.')
    values = {name:_raw_array(data, name, nt)[...,index] for name in (BM, FF, IA, FX, IY)}
    iy = np.zeros_like(values[FX]); iy[1:-1,1:-1] = values[IY]
    br, bz = (value.squeeze(0) for value in meqBrBz(values[FX]))
    return int(index), float(np.asarray(data['t']).reshape(-1)[index]), values, (values[FX], iy, br, bz)


def _plot_map(slot, panel, values, name, scale, unit, flux, current):
    plt.subplot(slot)
    axis = plt.gca()
    cmap = plt.rcParams['image.cmap'] + ('_r' if np.mean(values[IV]) < 0 else '')
    field = np.ma.masked_where(~IV, values)*scale
    image = plt.contourf(RRD, ZZD, field, levels=MAP_LEVELS, cmap=cmap, alpha=MAP_FILL_ALPHA)
    plt.contour(RRD, ZZD, field, levels=image.levels, cmap=cmap, linewidths=0.45, alpha=MAP_LINE_ALPHA)
    plot_approx_lcfs(flux, current, axis, colors='white', linestyles='-')
    # prep_plot_tcv(plt.gca(), lw=0.7, alpha=0.4, c='gray')
    # prep_plot_tcv(plt.gca(), lw=0.7, alpha=0.9, c='gray')
    prep_plot_tcv(plt.gca(), lw=0.7, alpha=0.8, c='k')
    plt.title(rf'({chr(96+panel)}) {name} [{unit}]'); 
    # plt.colorbar(image, fraction=0.07, pad=0.02)
    # plt.colorbar(image, fraction=0.06)
    color_axis = make_axes_locatable(axis).append_axes('right', size='12%', pad='3%')
    locator = color_axis.get_axes_locator()
    def colorbar_locator(ax, renderer):
        position = locator(ax, renderer)
        return position.shrunk(1, 0.9).anchored('C', position)
    colorbar_locator.get_subplotspec = locator.get_subplotspec
    color_axis.set_axes_locator(colorbar_locator)
    colorbar = plt.colorbar(image, cax=color_axis); colorbar.solids.set_edgecolor('face'); colorbar.ax.grid(False); plt.sca(axis)


def _plot_vector(slot, panel, values, name, scale, unit):
    plt.subplot(slot)
    channels = np.arange(1, values.size+1)
    plt.bar(channels, values.reshape(-1)*scale, color=C2)
    plt.xlabel('Channel index'); plt.title(f'({chr(96+panel)}) {name.replace("_", "")} [{unit}]')
    plt.xlim(0.2, values.size+0.8)


def _save_name(number): return 'raw_ds_example.pdf' if number == 1 else f'raw_ds_example_{number}.pdf'


def generate(shots=SHOTS, raw_dir=RAW_DIR, time_indices=TIME_INDICES):
    """Generate one raw-dataset snapshot per shot; the first is reserved for the main text."""
    if not shots: raise ValueError('Select at least one shot.')
    for number, shot in enumerate(shots, start=1):
        index, time, inputs, maps = _load_snapshot(raw_dir, shot, time_indices)
        # plt.figure(figsize=(THESIS_W, FIGURE_HEIGHT)); grid = plt.GridSpec(2, 1, height_ratios=(1.2, 0.8))
        plt.figure(figsize=(THESIS_W, FIGURE_HEIGHT)); grid = plt.GridSpec(2, 1, height_ratios=(1.3, 0.7))
        map_grid = grid[0].subgridspec(1, 4)
        vector_grid = grid[1].subgridspec(1, 4, width_ratios=(0.1, 1, 1, 1), wspace=0.3)  # Offset compensates the left tick-label overhang.
        map_slots = [map_grid[0, i] for i in range(4)]
        vector_slots = [vector_grid[0, i] for i in range(1, 4)]
        for panel, (slot, values, name, scale, unit) in enumerate(zip(map_slots, maps, MAP_NAMES, MAP_SCALES, MAP_UNITS), start=1): _plot_map(slot, panel, values, name, scale, unit, maps[0], maps[1])
        _plot_vector(vector_slots[0], 5, inputs[BM], 'B_m', 1.0, 'T')
        _plot_vector(vector_slots[1], 6, inputs[FF], 'F_f', 1.0, 'Wb')
        _plot_vector(vector_slots[2], 7, inputs[IA], 'I_a', 1e-3, 'kA')
        # plt.suptitle(rf'TCV \#{shot}, $t={time:.3f}\,$s (index {index})', y=0.965)
        # save_figure(_save_name(number), rect=(0, 0, 1, 0.985))
        plt.suptitle(rf'TCV \#{shot}, $t={time:.3f}\,$s (index {index})', y=0.97)
        save_figure(_save_name(number), rect=(0, 0, 1, 1))


if __name__ == '__main__': generate(); compile_inspection_pdf()
