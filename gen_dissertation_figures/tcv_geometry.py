# Figure description: This three-panel figure introduces the TCV geometry and the spatial interface used by the equilibrium surrogate. Panel (a) shows the native 28-by-65 FGE reconstruction grid inside the first-wall vessel outline. Panel (b) shows the final surrogate equilibrium from the baseline PID shape-control simulation as a poloidal-flux map, with the corresponding last closed flux surface identified by the saved boundary-flux value. Panel (c) overlays that LCFS with the active PID shape-control points from the same final simulation instant, showing where the controller queries the surrogate relative to the target plasma boundary.
from ._common import *


PID_DATA_PATH = 'data/pid_shape_net.mat'


def _format_panel():
    plt.xlim(min(RRD.min(), VESSO[:, 0].min()), max(RRD.max(), VESSO[:, 0].max()))
    plt.ylim(min(ZZD.min(), VESSO[:, 1].min()), max(ZZD.max(), VESSO[:, 1].max()))
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel(r'$r$ [m]'); plt.ylabel(r'$z$ [m]')
    for spine in plt.gca().spines.values(): spine.set_visible(False)


def generate():
    data = loadmat(PID_DATA_PATH)
    grid_r = data['fge_grid_r'].reshape(-1); grid_z = data['fge_grid_z'].reshape(-1)
    rr, zz = np.meshgrid(grid_r, grid_z); flux = np.squeeze(data['fge_final_flux'])
    boundary_flux = float(data['fge_final_boundary_flux'].reshape(-1)[0])
    i_sc = np.argmin(np.abs(data['sc_t'].reshape(-1)-data['t'].reshape(-1)[-1]))
    rc, zc = data['sc_r'][:,i_sc], data['sc_z'][:,i_sc]
    valid = (rc > 0) & data['active_shape_mask'].reshape(-1).astype(bool)
    assert flux.shape == rr.shape == zz.shape == (65, 28), 'Full-grid arrays must use (z, r) ordering'

    plt.figure(figsize=(THESIS_W, 5.0))

    plt.subplot(1, 3, 1)
    plot_vessel(); plt.scatter(rr, zz, s=2, marker='.', rasterized=True)
    plt.title(r'(a) Native grid' + '\n' + r'$28\times65$ points')
    _format_panel()

    plt.subplot(1, 3, 2)
    flux_map = plt.pcolormesh(rr, zz, flux, shading='nearest', cmap=plt.rcParams['image.cmap'] + '_r', rasterized=True)
    plt.contour(grid_r, grid_z, flux, levels=[boundary_flux], linestyles='--', linewidths=1.2)
    plot_vessel(); colorbar_axis = plt.gca().inset_axes([1.04, 0.0, 0.08, 1.0])
    plt.colorbar(flux_map, cax=colorbar_axis, label=r'Poloidal flux $\psi$ [Wb]')
    plt.title(r'(b) PID equilibrium' + '\n' + r'last closed flux surface')
    _format_panel()

    plt.subplot(1, 3, 3)
    plot_vessel(); plt.contour(grid_r, grid_z, flux, levels=[boundary_flux], linestyles='--', linewidths=1.2, alpha=0.5)
    plt.scatter(rc[valid], zc[valid], s=34, marker='x', color=C2)
    plt.title(r'(c) Shape control' + '\n' + rf'${valid.sum()}$ active points')
    _format_panel()

    plt.figtext(0.5, 0.035, 'Active points at the final PID simulation instant.\n' + r'The inference interface accepts dynamic $(r,z)$ coordinates.', ha='center')
    save_figure('tcv_geometry.pdf', rect=(0, 0.075, 1, 1))

if __name__ == '__main__': generate(); compile_inspection_pdf()