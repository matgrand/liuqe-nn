# Figure description: This single-panel figure introduces the TCV magnetic-equilibrium and shape-control geometry. The vessel encloses an example poloidal-flux reconstruction, with the last closed flux surface and representative shape-control evaluation points highlighted. Around the vessel, the figure locates the 19 independently driven poloidal-field circuits: the OH1 and OH2 transformer coils, the eight E and eight F shaping coils, and the in-vessel G coil used for fast vertical stabilization. It also shows the positions and orientations of the poloidal magnetic-field pick-up coils Bm and the positions of the poloidal-flux loops Ff used by LIUQE. The composition follows the information hierarchy of Fig. 1 in Mele et al., “First experimental demonstration of plasma shape control in a tokamak through Model Predictive Control,” but uses the dissertation plotting style and the locally saved LIUQE example rather than copying the source figure's colors or typography.
from ._common import *
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import ScaledTranslation


EQUILIBRIUM_DATA_PATH = 'data/pid_shape_net.mat'
HARDWARE_DATA_PATH = 'tcv_params/tcv_geometry_hardware.mat'
PROBE_HALF_LENGTH = 0.018  # [m]
COIL_LABEL_OFFSET = 1.0  # [pt], corrects the optical center for unused font descender space.


def _connected_components(r, z, w, h):
    adjacent = (np.abs(r[:,None]-r) <= (w[:,None]+w)/2+0.005) & (np.abs(z[:,None]-z) <= (h[:,None]+h)/2+0.005)
    unseen, components = set(range(len(r))), []
    while unseen:
        stack, component = [unseen.pop()], []
        while stack:
            i = stack.pop(); component.append(i)
            neighbors = set(np.flatnonzero(adjacent[i])) & unseen
            unseen -= neighbors; stack.extend(neighbors)
        components.append(np.asarray(component))
    return components


def _plot_coils(hardware):
    rw, zw = hardware['rw'].reshape(-1), hardware['zw'].reshape(-1)
    ww, hw, twa = hardware['ww'].reshape(-1), hardware['hw'].reshape(-1), hardware['Twa']
    names = [str(name[0]) for name in hardware['dima'].reshape(-1)]
    for j, raw_name in enumerate(names):
        selected = np.abs(twa[:,j]) > 0; r, z, w, h = rw[selected], zw[selected], ww[selected], hw[selected]
        name = raw_name.replace('_00', '')
        finite = np.isfinite(w) & np.isfinite(h)
        if not finite.any():
            for upper in (False, True):
                side = z > 0 if upper else z < 0
                plt.scatter(r[side], z[side], s=18, marker='o', color=C1)
                plt.text(r[side].mean()-0.035, z[side].mean(), r'$G$', ha='right', va='center')
            continue
        r, z, w, h = r[finite], z[finite], w[finite], h[finite]
        for component in _connected_components(r, z, w, h):
            left = np.min(r[component]-w[component]/2); right = np.max(r[component]+w[component]/2)
            bottom = np.min(z[component]-h[component]/2); top = np.max(z[component]+h[component]/2)
            plt.fill([left, right, right, left], [bottom, bottom, top, top], color=C1, alpha=0.28)
            plt.plot([left, right, right, left, left], [bottom, bottom, top, top, bottom], color=C1)
            rotation = 90 if top-bottom > right-left else 0
            angle = np.deg2rad(rotation); offset = ScaledTranslation(COIL_LABEL_OFFSET*np.sin(angle)/72, -COIL_LABEL_OFFSET*np.cos(angle)/72, plt.gcf().dpi_scale_trans)
            plt.text((left+right)/2, (bottom+top)/2, name, rotation=rotation, ha='center', va='center', transform=plt.gca().transData+offset)


def _plot_legend():
    handles = [
        Line2D([], [], color=C2, linewidth=0.8, label=r'Flux contours $psi$'),
        Line2D([], [], color='k', linestyle='--', linewidth=1.5, label='LCFS'),
        Patch(facecolor=C1, edgecolor=C1, alpha=0.28, label='PF coils'),
        Line2D([], [], color=C2, linewidth=1.6, label=r'\texttt{Bm} probes'),
        Line2D([], [], color=C3, marker='x', linestyle='none', markersize=6, label=r'\texttt{Ff} loops'),
        Line2D([], [], color=C4, marker='o', linestyle='none', markersize=5, label='Shape-control points'),
    ]
    plt.legend(handles=handles, loc='center left', bbox_to_anchor=(1.03, 0.5))


def generate():
    data, hardware = loadmat(EQUILIBRIUM_DATA_PATH), loadmat(HARDWARE_DATA_PATH)
    grid_r = data['fge_grid_r'].reshape(-1); grid_z = data['fge_grid_z'].reshape(-1)
    rr, zz = np.meshgrid(grid_r, grid_z); flux = np.squeeze(data['fge_final_flux'])
    boundary_flux = float(data['fge_final_boundary_flux'].reshape(-1)[0])
    i_sc = np.argmin(np.abs(data['sc_t'].reshape(-1)-data['t'].reshape(-1)[-1]))
    rc, zc = data['sc_r'][:,i_sc], data['sc_z'][:,i_sc]
    valid = (rc > 0) & data['active_shape_mask'].reshape(-1).astype(bool)
    rm, zm, am = (hardware[name].reshape(-1) for name in ('rm', 'zm', 'am'))
    rf, zf = (hardware[name].reshape(-1) for name in ('rf', 'zf'))
    inside = inpoly(np.stack((rr, zz), axis=-1)); masked_flux = np.ma.masked_where(~inside, flux)
    probe_r = np.column_stack((rm-PROBE_HALF_LENGTH*np.cos(am), rm+PROBE_HALF_LENGTH*np.cos(am), np.full_like(rm, np.nan))).ravel()
    probe_z = np.column_stack((zm-PROBE_HALF_LENGTH*np.sin(am), zm+PROBE_HALF_LENGTH*np.sin(am), np.full_like(zm, np.nan))).ravel()

    plt.figure(figsize=(THESIS_W, 6.76))
    plt.contour(grid_r, grid_z, masked_flux, levels=16, cmap=plt.rcParams['image.cmap'], linewidths=0.6, alpha=0.55)
    plt.contour(grid_r, grid_z, flux, levels=[boundary_flux], linestyles='--', linewidths=1.5)
    plot_vessel(); _plot_coils(hardware)
    plt.plot(probe_r, probe_z, color=C2, linewidth=1.6)
    plt.scatter(rf, zf, s=20, marker='x', color=C3)
    plt.scatter(rc[valid], zc[valid], s=24, marker='o', color=C4, zorder=6)
    plt.xlim(0.35, 1.42); plt.ylim(-1.25, 1.25); plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel(r'$r$ [m]'); plt.ylabel(r'$z$ [m]'); _plot_legend()
    save_figure('mag_meas.pdf')


if __name__ == '__main__': generate(); compile_inspection_pdf()
