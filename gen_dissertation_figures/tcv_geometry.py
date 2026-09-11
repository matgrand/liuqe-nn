# Figure description: This three-panel figure introduces the spatial quantities connecting LIUQE to the control-oriented surrogate. Panel (a) shows the native 28-by-65 equilibrium-reconstruction grid relative to the TCV vessel. Panel (b) gives an example poloidal-flux map on that grid using translucent filled contours and matching contour lines, and identifies the last closed flux surface (LCFS). Panel (c) overlays example shape-control evaluation points on the LCFS, illustrating the arbitrary spatial queries supplied to the surrogate. Together, the panels connect LIUQE's grid-based equilibrium reconstruction to the arbitrary spatial evaluations required for shape control. The equilibrium and control-point arrays remain the existing saved example and are presented generically here rather than as a comparison of controllers.
from ._common import *


PID_DATA_PATH = 'data/pid_shape_net.mat'
MAP_LEVELS, MAP_FILL_ALPHA, MAP_LINE_ALPHA = 16, 0.85, 1


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
    points = np.stack((rr, zz), axis=-1); inside = inpoly(points)
    assert flux.shape == rr.shape == zz.shape == (65, 28), 'Full-grid arrays must use (z, r) ordering'

    plt.figure(figsize=(THESIS_W, 5.0))

    plt.subplot(1, 3, 1)
    plot_vessel(c='k', lw=0.7, alpha=0.8); plt.scatter(rr, zz, s=2, marker='.', rasterized=False)
    plt.title(r'(a) LIUQE grid')
    _format_panel()

    plt.subplot(1, 3, 2)
    field = np.ma.masked_where(~inside, flux); cmap = plt.rcParams['image.cmap'] + '_r'
    contours = plt.contourf(grid_r, grid_z, field, levels=MAP_LEVELS, cmap=cmap, alpha=MAP_FILL_ALPHA)
    plt.contour(grid_r, grid_z, field, levels=contours.levels, cmap=cmap, linewidths=0.45, alpha=MAP_LINE_ALPHA)
    plt.contour(grid_r, grid_z, flux, levels=[boundary_flux], linestyles='--', linewidths=1.2)
    plot_vessel(c='k', lw=0.7, alpha=0.8); plt.title(r'(b) Poloidal flux $\psi$' + '\n' + r'and LCFS')
    _format_panel()

    plt.subplot(1, 3, 3)
    plot_vessel(c='k', lw=0.7, alpha=0.8); plt.contour(grid_r, grid_z, flux, levels=[boundary_flux], linestyles='--', linewidths=1.2, alpha=0.5)
    plt.scatter(rc[valid], zc[valid], s=34, marker='x', color=C2)
    plt.title(r'(c) Shape-control' + '\n' + r'example points')
    _format_panel()

    save_figure('tcv_geometry.pdf')

if __name__ == '__main__': generate(); compile_inspection_pdf()
