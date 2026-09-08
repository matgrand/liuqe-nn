from ._common import *
from ._common import _cmap


def generate():
    # %% ## FMAG control simulation
    ## FMAG control simulation
    CTRL_SIM_SMALL = 'data/fmag_net.mat'
    CTRL_SIM_STD = 'data/fmag_std.mat'
    # expected keys: 't', 'rY', 'zY', 'Ip', 'Va', 'mixed_errors', 'observables', 'xpt_flux', 'xpt_flux_lih'
    ctrl_start = 0.6  # controller start time
    vlinec = 'gray'
    vlines = '--'
    std_style = '--'

    sim_net = loadmat(CTRL_SIM_SMALL)
    sim_std = loadmat(CTRL_SIM_STD)

    def _vec(x): return np.asarray(x).reshape(-1)

    def _ch_time(x):
        """Return array shaped (n_channels, n_time)."""
        a = np.asarray(x)
        if a.ndim == 1:
            return a.reshape(1, -1)
        if a.shape[0] > a.shape[1]:
            # likely time x channels -> transpose
            return a.T
        return a

    tN = _vec(sim_net['t'])
    tL = _vec(sim_std['t'])
    IpN, IpL = _vec(sim_net['Ip']), _vec(sim_std['Ip'])
    mixN = sim_net['mixed_errors'].T 
    mixL = sim_std['mixed_errors'].T
    err_dr_net = _vec(sim_net['xpt_err_delta_r'])
    err_dr_std = _vec(sim_std['xpt_err_delta_r'])
    err_dz_net = _vec(sim_net['xpt_err_delta_z'])
    err_dz_std = _vec(sim_std['xpt_err_delta_z'])

    # Figure 4: Mixed errors (first 6 channels)
    n_show = 6
    plt.figure(figsize=(THESIS_W, 2.8))

    plt.subplot(1,2,1)
    plt.plot(tN, 100*err_dr_net / IpN, color=C1, label=r'$\Delta r$, surrogate')
    plt.plot(tL, 100*err_dr_std / IpL, '--', color=C1, label=r'$\Delta r$, LIH')
    plt.plot(tN, -100*err_dz_net / IpN, color=C2, label=r'$\Delta z$, surrogate')
    plt.plot(tL, -100*err_dz_std / IpL, '--', color=C2, label=r'$\Delta z$, LIH')
    plt.axhline(y=0, color='k', linestyle=':')
    plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
    plt.xlabel(r'Time $t$ [s]'), plt.ylabel(r'X-point error [cm]')
    plt.legend(loc='best')

    plt.subplot(1,2,2)
    colors = _cmap(np.linspace(0, 1, n_show))
    print(f'max values: {np.max(mixN[:, :n_show], axis=0)}')
    scalings = [2e3, 1e3, 1e1, 1e1, 1, 1] # scale first channel (shape error) to be visible
    for i in range(n_show):
        plt.plot(tN, mixN[:, i]/scalings[i], color=colors[i], label='Surrogate' if i==0 else None)
        plt.plot(tL, mixL[:, i]/scalings[i], std_style, color=colors[i], label='LIH' if i==0 else None)
    # plt.plot(tN, mixN[:, :n_show], lw=lw)
    plt.xlabel(r'Time $t$ [s]'), plt.ylabel(r'Mixed error [a.u.]'), plt.legend(loc='best')

    save_figure('fmag_control_simulation.pdf')
