from ._common import *
from ._common import _cmap


def generate():
    """Compare the PID closed loops using the ten channels in the reported RMS metric."""
    sim_net = loadmat('data/pid_shape_net.mat')
    sim_std = loadmat('data/pid_shape_std.mat')
    t = sim_net['t'].reshape(-1)
    ctrl_start, dec, n_metric = 0.6, 5, 10

    se_net = sim_net['shape_errs'][:n_metric].T * 1e3
    se_std = sim_std['shape_errs'][:n_metric].T * 1e3
    active_proj = np.flatnonzero(np.max(np.abs(sim_net['proj_errs']), axis=1) > 0)
    pe_net = sim_net['proj_errs'][active_proj].T * 1e3
    pe_std = sim_std['proj_errs'][active_proj].T * 1e3
    post = t >= ctrl_start
    rms_t_net = np.sqrt(np.mean(se_net**2, axis=1)); rms_t_std = np.sqrt(np.mean(se_std**2, axis=1))
    rms_ch_net = np.sqrt(np.mean(se_net[post]**2, axis=0)); rms_ch_std = np.sqrt(np.mean(se_std[post]**2, axis=0))
    rms_net = np.sqrt(np.mean(se_net[post]**2)); rms_std = np.sqrt(np.mean(se_std[post]**2))

    plt.figure(figsize=(THESIS_W, 5.6))
    plt.subplot(2, 2, 1)
    colors = _cmap(np.linspace(0, 1, n_metric))
    for i, color in enumerate(colors):
        plt.plot(t[::dec], se_net[::dec, i], color=color, label='Surrogate' if i == 0 else None)
        plt.plot(t[::dec], se_std[::dec, i], '--', color=color, label='FGE/qintmex' if i == 0 else None)
    plt.axvline(ctrl_start, color='gray', linestyle=':'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'Shape error $e_i$ [mm]'); plt.legend()

    plt.subplot(2, 2, 2)
    colors = _cmap(np.linspace(0, 1, len(active_proj)))
    for i, color in enumerate(colors):
        plt.plot(t[::dec], pe_net[::dec, i], color=color, label=rf'Mode {active_proj[i] + 1}')
        plt.plot(t[::dec], pe_std[::dec, i], '--', color=color)
    plt.axvline(ctrl_start, color='gray', linestyle=':'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'Projected error [mm]'); plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(t[::dec], rms_t_net[::dec], color=C1, label=rf'Surrogate ({rms_net:.2f} mm)')
    plt.plot(t[::dec], rms_t_std[::dec], '--', color=C3, label=rf'FGE/qintmex ({rms_std:.2f} mm)')
    plt.axvline(ctrl_start, color='gray', linestyle=':'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'Instantaneous RMS [mm]'); plt.legend(title='Post-activation RMS')

    plt.subplot(2, 2, 4)
    channels = np.arange(1, n_metric + 1)
    plt.plot(channels, rms_ch_net, 'o-', color=C1, label='Surrogate')
    plt.plot(channels, rms_ch_std, 's--', color=C3, label='FGE/qintmex')
    plt.xticks(channels); plt.xlabel('Saved shape-error channel'); plt.ylabel(r'Post-activation RMS [mm]'); plt.legend()

    save_figure('pid_shape_control_simulation.pdf')
