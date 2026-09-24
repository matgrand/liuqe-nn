# Figure description: Three TCV-section schematics explain the local vertical-force response of an elongated plasma. The same small upward displacement is shown in all panels. Curved external poloidal-field lines imply positive, zero, or negative radial-field gradient at the plasma, producing downward restoring force, zero first-order force, or upward destabilizing force for positive toroidal plasma current. The vessel is drawn with the repository's TCV geometry helper, and the equilibrium outline is dashed. This is a local force schematic, not a quantitative TCV equilibrium or an eddy-current calculation; passive and active stabilization are explained in the surrounding text.
from ._common import *
from matplotlib.patches import Ellipse


Z_SHIFT = 0.16
FIELD_CURVATURE = 0.12
_r0s = list(np.linspace(0.62, 1.14, 4))
ell_r = 0.88

def _panel(index, gradient, heading, ell_za=0.95):
    plt.subplot(1, 3, index)
    plot_vessel(c='k', lw=0.7, alpha=0.35)
    z = np.linspace(-0.97, 0.97, 180)
    # for r0 in (0.62, ell_r, 1.14):
    for r0 in _r0s:
        # B_Z < 0: dR/dZ = B_R/B_Z, so leftward curvature gives dB_R/dZ > 0.
        r = r0 - gradient*FIELD_CURVATURE*z*z
        plt.plot(r, z, color=C1, lw=1.0, zorder=2)
        plt.annotate('', xy=(r[31], z[31]), xytext=(r[45], z[45]), arrowprops=dict(arrowstyle='-|>', color=C1, lw=1.0))
        _off = 110
        plt.annotate('', xy=(r[31+_off], z[31+_off]), xytext=(r[45+_off], z[45+_off]), arrowprops=dict(arrowstyle='-|>', color=C1, lw=1.0))
    axis = plt.gca()
    axis.add_patch(Ellipse((ell_r, 0), 0.33, ell_za, fill=False, ls='--', lw=0.8, edgecolor='0.45', zorder=3))
    axis.add_patch(Ellipse((ell_r, Z_SHIFT), 0.33, ell_za, facecolor=C2, edgecolor=C2, alpha=0.23, lw=1.2, zorder=4))
    axis.add_patch(Ellipse((ell_r, Z_SHIFT), 0.33, ell_za, fill=False, edgecolor=C2, lw=1.2, zorder=5))
    plt.plot(ell_r, Z_SHIFT, marker='x', color='k', ms=6, mew=1.2, zorder=6)
    if gradient:
        br = 1 if gradient > 0 else -1
        plt.annotate('', xy=(ell_r+0.15*br, Z_SHIFT), xytext=(ell_r, Z_SHIFT), arrowprops=dict(arrowstyle='-|>', color=C3, lw=1.4), zorder=7)
        # plt.annotate('', xy=(1.08, Z_SHIFT-0.22*br), xytext=(1.08, Z_SHIFT), arrowprops=dict(arrowstyle='-|>', color='k', lw=1.4), zorder=7)        
        plt.annotate('', xy=(ell_r, Z_SHIFT-0.22*br), xytext=(ell_r, Z_SHIFT), arrowprops=dict(arrowstyle='-|>', color='k', lw=1.4), zorder=7)
        plt.text(ell_r, Z_SHIFT-0.27*br, r'$F_Z$', ha='center', va='center')
    else:
        # plt.text(ell_r, -0.87, r'$B_R=0,\ F_Z=0$', ha='center', va='center')
        plt.text(ell_r, 1.02, r'$B_R=0,\ F_Z=0$', ha='center', va='center')
    plt.title(heading)
    plt.xlim(0.48, 1.28); plt.ylim(-1.08, 1.08)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xticks([]); plt.yticks([])
    for spine in plt.gca().spines.values(): spine.set_visible(False)


def generate():
    plt.figure(figsize=(THESIS_W, 5))
    _panel(1, 1, r'(a) $\partial_Z B_R>0$'+'\n'+r'restoring', ell_za=0.65)
    _panel(2, 0, r'(b) $\partial_Z B_R=0$'+'\n'+r'neutral', ell_za=0.75)
    _panel(3, -1, r'(c) $\partial_Z B_R<0$'+'\n'+r'destabilizing', ell_za=0.85)
    save_figure('magnetic_vertical_stability.pdf')


if __name__ == '__main__': generate(); compile_inspection_pdf()
