from matplotlib.lines import Line2D

from ._common import *
from ._common import _cmap


def generate():
    """Compare MPC regulation, its applied first move, and QP iteration counts."""
    sim_net = loadmat('data/mpc_shape_net.mat')
    sim_std = loadmat('data/mpc_shape_std.mat')
    t = sim_net['t'].reshape(-1)
    ctrl_start, opt_start, dec, n_metric, n_move = 0.6, 0.62, 5, 10, 16

    se_net = sim_net['shape_errs'][:n_metric].T * 1e3
    se_std = sim_std['shape_errs'][:n_metric].T * 1e3
    du_net = sim_net['duopt'][:n_move].T
    du_std = sim_std['duopt'][:n_move].T
    post = t >= opt_start
    score = np.sqrt(np.mean(du_net[post]**2, axis=0) + np.mean(du_std[post]**2, axis=0))
    shown = np.sort(np.argsort(score)[-6:])
    rms_t_net = np.sqrt(np.mean(se_net**2, axis=1)); rms_t_std = np.sqrt(np.mean(se_std**2, axis=1))
    rms_net = np.sqrt(np.mean(se_net[post]**2)); rms_std = np.sqrt(np.mean(se_std[post]**2))
    updates = np.r_[True, np.any(np.diff(du_net, axis=0) != 0, axis=1)] & post

    plt.figure(figsize=(THESIS_W, 5.6))
    plt.subplot(2, 2, 1)
    colors = _cmap(np.linspace(0, 1, n_metric))
    for i, color in enumerate(colors):
        plt.plot(t[::dec], se_net[::dec, i], color=color, label='Surrogate' if i == 0 else None)
        plt.plot(t[::dec], se_std[::dec, i], '--', color=color, label='FGE/qintmex' if i == 0 else None)
    plt.axvline(ctrl_start, color='gray', linestyle=':'); plt.axvline(opt_start, color='gray', linestyle='--')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'Shape error $e_i$ [mm]'); plt.legend()

    plt.subplot(2, 2, 2)
    colors = _cmap(np.linspace(0, 1, len(shown)))
    for i, mode in enumerate(shown):
        plt.plot(t[::dec], du_net[::dec, mode], color=colors[i], label=rf'$j={mode + 1}$')
        plt.plot(t[::dec], du_std[::dec, mode], '--', color=colors[i])
    style_handles = [Line2D([0], [0], color='black', label='Surrogate'), Line2D([0], [0], color='black', linestyle='--', label='FGE/qintmex')]
    mode_legend = plt.legend(title='Largest saved components', ncol=2, loc='upper right'); plt.gca().add_artist(mode_legend)
    plt.legend(handles=style_handles, loc='lower right')
    plt.axvline(opt_start, color='gray', linestyle='--'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'First QP move $\Delta u_{k|k,j}$ [stored units]')

    plt.subplot(2, 2, 3)
    plt.plot(t[::dec], rms_t_net[::dec], color=C1, label=rf'Surrogate ({rms_net:.2f} mm)')
    plt.plot(t[::dec], rms_t_std[::dec], '--', color=C3, label=rf'FGE/qintmex ({rms_std:.2f} mm)')
    plt.axvline(ctrl_start, color='gray', linestyle=':'); plt.axvline(opt_start, color='gray', linestyle='--')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'Instantaneous RMS [mm]'); plt.legend(title=r'RMS for $t\geq0.62$ s')

    plt.subplot(2, 2, 4)
    plt.step(t[updates], sim_net['niter'].reshape(-1)[updates], where='post', color=C1, label='Surrogate')
    plt.step(t[updates], sim_std['niter'].reshape(-1)[updates], where='post', linestyle='--', color=C3, label='FGE/qintmex')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel('QP solver iterations'); plt.ylim(3.8, 5.2); plt.yticks([4, 5]); plt.legend(title='2 ms updates')

    save_figure('mpc_shape_control_simulation.pdf')
