from ._common import *


EQUILIBRIUM_INDEX = 297  # representative single-null evaluation example
QUERY_R = np.array([1.0864, 1.0313, 0.9357, 0.8253, 0.7297, 0.6746, 0.6746, 0.7297, 0.8253, 0.9357, 1.0313, 1.0864, 1.0137, 0.9747, 0.8805, 0.7863, 0.7473, 0.7863, 0.8805, 0.9747, 0.9182, 0.9182, 0.8428, 0.8428, 0.8805])
QUERY_Z = np.array([0.1574, 0.4299, 0.5873, 0.5873, 0.4299, 0.1574, -0.1574, -0.4299, -0.5873, -0.5873, -0.4299, -0.1574, 0.0, 0.2687, 0.38, 0.2687, 0.0, -0.2687, -0.38, -0.2687, -0.1075, 0.1075, 0.1075, -0.1075, 0.0])


def _format_panel():
    plt.xlim(min(RRD.min(), VESSO[:, 0].min()), max(RRD.max(), VESSO[:, 0].max()))
    plt.ylim(min(ZZD.min(), VESSO[:, 1].min()), max(ZZD.max(), VESSO[:, 1].max()))
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel(r'$r$ [m]'); plt.ylabel(r'$z$ [m]')


def generate():
    fx, iy = ds.fg[EQUILIBRIUM_INDEX, :2]
    assert fx.shape == iy.shape == RRD.shape == (65, 28), 'Full-grid arrays must use (z, r) ordering'
    assert QUERY_R.shape == QUERY_Z.shape == (25,), 'The real-time analysis set must contain 25 points'

    plt.figure(figsize=(THESIS_W, 5.0))

    plt.subplot(1, 3, 1)
    plot_vessel(); plt.scatter(RRD, ZZD, s=2, marker='.', rasterized=True)
    plt.title(r'(a) Native grid' + '\n' + r'$28\times65$ points')
    _format_panel()

    plt.subplot(1, 3, 2)
    flux_map = plt.pcolormesh(RRD, ZZD, fx, shading='nearest', cmap=plt.rcParams['image.cmap'] + '_r', rasterized=True)
    plt.contour(RRD, ZZD, np.abs(iy) > 0, levels=[0.5], linestyles='--', linewidths=1.2)
    plot_vessel(); colorbar_axis = plt.gca().inset_axes([1.04, 0.0, 0.08, 1.0])
    plt.colorbar(flux_map, cax=colorbar_axis, label=r'Poloidal flux $\psi$ [Wb]')
    plt.title(fr'(b) Equilibrium {EQUILIBRIUM_INDEX}' + '\n' + r'current-support boundary')
    _format_panel()

    plt.subplot(1, 3, 3)
    plot_vessel(); plt.scatter(QUERY_R, QUERY_Z, s=14)
    plt.title(r'(c) Analysis set' + '\n' + r'$25$ query points')
    _format_panel()

    plt.figtext(0.5, 0.035, 'Fixed point set for this analysis.\n' + r'The inference interface accepts dynamic $(r,z)$ coordinates.', ha='center')
    save_figure('tcv_geometry.pdf', rect=(0, 0.075, 1, 1))
