# Figure description: This twelve-panel figure recreates the TCV plasma-configuration sequence in fig_useful_material/tcv_shape_flexibility_draft.png as three rows of four panels at dissertation width. The panels show a limited plasma, upper and lower single-null divertors, a double-null divertor, high squareness, high plasma current, high fully ECCD-driven current, a pear shape, a doublet, high elongation, negative triangularity, and a snowflake. The doublet uses shot 87751 at 0.687 s from ft_ds in place of the draft's shot 10159. Each panel identifies its shot and instant, outlines the TCV vessel, draws pale exterior poloidal-flux contours and nested filled contours inside the last closed flux surface, and marks the magnetic axis. The colors use the dissertation's inferno palette. The eleven older cases come from MATLAB exports of LIUQE results; the raw doublet's boundary flux is estimated from the edge of its nonzero-current support. Unavailable cases remain visibly labelled so the other panels can still be inspected.
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.ndimage import binary_dilation, binary_erosion, label

from . import _common as common
from ._common import THESIS_W, RRD, ZZD, VESSO, inpoly, plot_vessel, compile_inspection_pdf

DATA_DIR = Path('data/tcv_shape_flexibility')
FT_DATA_DIR = Path('ft_ds')
OUTPUT_PATH = Path('dissertation_figures/tcv_shape_flexibility.pdf')
SHOW_PLOTS = False
DOUBLET, DOUBLET_TIME = 87751, 0.687 
# DOUBLET, DOUBLET_TIME = 87809, 0.495
# DOUBLET, DOUBLET_TIME = 87869, 0.880

FONTSZ = 7

CASES = [
    (9849, 0.5, r'Limited' + '\n' + r'$\delta=0$, $I_p=230\,\mathrm{kA}$'),
    (6010, 0.8, r'SN upper' + '\n' + r'$I_p=330\,\mathrm{kA}$'),
    (5650, 0.8, r'SN lower' + '\n' + r'$I_p=335\,\mathrm{kA}$'),
    (8856, 0.29, r'DN' + '\n' + r'$I_p=325\,\mathrm{kA}$'),
    (11962, 1.0, r'Highest' + '\n' + r'squareness' + '\n' + r'$\lambda=0.5$'),
    (11368, 0.65, r'Highest current' +  '\n' + r'$I_p=1.06\,\mathrm{MA}$'),
    (18548, 1.5, r'Max purely' + '\n' +r'driven ECCD current' + '\n' + r'$I_p=210\,\mathrm{kA}$'),
    (6442, 0.5, r'Pear shape' + '\n' + r'$I_p=360\,\mathrm{kA}$'),
    (DOUBLET, DOUBLET_TIME, r'Doublet shape' + '\n' + r'$I_p=368\,\mathrm{kA}$'),
    (19373, 0.42, r'Highest elongation' + '\n' + r'$\kappa=2.80$'),
    (11928, 0.67, r'Lowest' + '\n' + r'triangularity' + '\n' + r'$\delta=-0.77$'),
    (36151, 0.457, r'Snowflake' + '\n' + r'$I_p=230\,\mathrm{kA}$'),
]

doublet_in_ft_ds = [(87809, 0.495), (87812, 0.591), (87869, 0.880), (87751, 0.687)] # dont delte this


def _load_equilibrium(path, shot, time):
    if shot == DOUBLET:
        d = loadmat(path, variable_names=['t', 'Fx', 'Iy'])
        times = d['t'].ravel(); k = np.argmin(abs(times-time))
        if abs(times[k]-time) > 0.0005: raise ValueError(f'No time slice near {time} s in {path}')
        psi = d['Fx'][..., k]
        current = np.zeros_like(psi); current[1:-1, 1:-1] = d['Iy'][..., k]
        support = current != 0; border = support & ~binary_erosion(support)
        if not np.any(border): raise ValueError(f'No plasma-current support in {path}')
        psi_boundary = float(np.mean(psi[border]) + 0.026*np.ptp(psi))
        psi_axis = float(np.min(psi[support]))
        return RRD[0], ZZD[:, 0], psi, psi_axis, psi_boundary
    d = loadmat(path)
    return d['r'].ravel(), d['z'].ravel(), np.squeeze(d['psi']), float(d['psi_axis'].item()), float(d['psi_boundary'].item())


def _plot_equilibrium(path, shot, time):
    r, z, psi, psi_axis, psi_boundary = _load_equilibrium(path, shot, time)
    if psi.shape != (len(z), len(r)): raise ValueError(f'Unexpected flux shape {psi.shape} in {path}')
    if not np.all(np.isfinite(psi)) or psi_axis == psi_boundary: raise ValueError(f'Invalid flux values in {path}')
    rr, zz = np.meshgrid(r, z)
    inside = inpoly(np.stack((rr, zz), axis=-1))
    normalized = (psi - psi_boundary) / (psi_axis - psi_boundary)
    flux = np.ma.masked_where(~inside, normalized)
    regions, nregions = label((normalized > 0) & inside)
    if not nregions: raise ValueError(f'No positive-flux plasma region in {path}')
    areas = np.bincount(regions.ravel()); areas[0] = 0
    main = np.argmax(areas)
    plasma = regions == main
    if shot == DOUBLET: plasma |= regions == np.argsort(areas)[-2]  # The doublet has two plasma lobes.
    plasma_flux = np.ma.masked_where(~binary_dilation(plasma) | ~inside, normalized)
    cmap = plt.get_cmap(plt.rcParams['image.cmap'])
    plt.contour(r, z, flux, levels=np.arange(-1.0, 0.0, 0.1), colors='0.83', linewidths=0.35, zorder=1)
    plt.contourf(r, z, plasma_flux, levels=np.linspace(0, 1, 12), cmap=cmap, zorder=2)
    plt.contour(r, z, flux, levels=[0], colors='black', linewidths=0.8, zorder=3)
    k = np.argmax(np.where(regions == main, normalized, -np.inf))
    plt.plot(rr.flat[k], zz.flat[k], '+', color=common.C2, markersize=5, markeredgewidth=0.7, zorder=5)


def generate(data_dir=DATA_DIR):
    """Plot the requested snapshots; leave labelled panels for missing exports."""
    plt.figure(figsize=(THESIS_W, 7))
    missing = []
    for i, (shot, time, caption) in enumerate(CASES, start=1):
        plt.subplot(3, 4, i)
        path = FT_DATA_DIR/f'{shot}.mat' if shot == DOUBLET else Path(data_dir)/f'shot_{shot}_t_{time:0.3f}.mat'
        try: _plot_equilibrium(path, shot, time)
        except (OSError, KeyError, ValueError) as err:
            missing.append((shot, time)); print(f'Unavailable {path}: {err}')
            plt.text(0.5, 0.5, 'No equilibrium\navailable', ha='center', va='center', transform=plt.gca().transAxes)
        plot_vessel(c='0.3', lw=0.6, alpha=0.85)
        plt.xlim(VESSO[:, 0].min()-0.02, VESSO[:, 0].max()+0.02)
        plt.ylim(VESSO[:, 1].min()-0.02, VESSO[:, 1].max()+0.02)
        plt.gca().set_aspect('equal', adjustable='box'); plt.axis('off')
        # plt.title(rf'\#{shot}, $t={time:g}\,\mathrm{{s}}$', fontsize=FONTSZ)
        titl = rf'\#{shot}, $t={time:g}\,\mathrm{{s}}$' 
        plt.text(0.5, -0.04, rf'({chr(96+i)}) ' + titl + '\n' + caption, ha='center', va='top', transform=plt.gca().transAxes, fontsize=FONTSZ)
    OUTPUT_PATH.parent.mkdir(exist_ok=True)
    plt.tight_layout()
    plt.savefig(OUTPUT_PATH); print(f'Saved {OUTPUT_PATH} ({len(CASES)-len(missing)}/{len(CASES)} equilibria)')
    if SHOW_PLOTS or common.SHOW_PLOTS: plt.show(block=common.WAIT_FOR_PLOTS)
    else: plt.close()
    return missing


if __name__ == '__main__': generate(); compile_inspection_pdf()
