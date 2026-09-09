from ._common import *


def generate():
    """Compare FMag X-point regulation and its plant-level consequences."""
    sim_net = loadmat('data/fmag_net.mat')
    sim_std = loadmat('data/fmag_std.mat')

    def vec(data, key): return np.asarray(data[key]).reshape(-1)

    t = vec(sim_net, 't'); dec = 5
    ip_net, ip_std = vec(sim_net, 'Ip'), vec(sim_std, 'Ip')
    dr_net = 1e3 * vec(sim_net, 'xpt_err_delta_r') / ip_net
    dr_std = 1e3 * vec(sim_std, 'xpt_err_delta_r') / ip_std
    dz_net = 1e3 * vec(sim_net, 'xpt_err_delta_z') / ip_net
    dz_std = 1e3 * vec(sim_std, 'xpt_err_delta_z') / ip_std
    enorm_net = np.hypot(dr_net, dz_net); enorm_std = np.hypot(dr_std, dz_std)
    rms_net = np.sqrt(np.mean(enorm_net**2)); rms_std = np.sqrt(np.mean(enorm_std**2))
    obs_dr_net = 1e3 * vec(sim_net, 'xpt_delta_r'); obs_dr_lih = 1e3 * vec(sim_net, 'xpt_delta_r_lih')
    obs_dz_net = 1e3 * vec(sim_net, 'xpt_delta_z'); obs_dz_lih = 1e3 * vec(sim_net, 'xpt_delta_z_lih')
    dr_y = 1e3 * (vec(sim_net, 'rY') - vec(sim_std, 'rY'))
    dz_y = 1e3 * (vec(sim_net, 'zY') - vec(sim_std, 'zY'))
    dip = 1e-3 * (ip_net - ip_std)
    va_net = np.asarray(sim_net['Va']); va_std = np.asarray(sim_std['Va'])
    va_norm_net = np.linalg.norm(va_net, axis=0); va_norm_std = np.linalg.norm(va_std, axis=0)
    va_diff = np.linalg.norm(va_net - va_std, axis=0)

    plt.figure(figsize=(THESIS_W, 7.2))
    plt.subplot(3, 2, 1)
    plt.plot(t[::dec], dr_net[::dec], color=C1, label=r'$\Delta r$, surrogate')
    plt.plot(t[::dec], dr_std[::dec], '--', color=C1, label=r'$\Delta r$, LIH')
    plt.plot(t[::dec], dz_net[::dec], color=C3, label=r'$\Delta z$, surrogate')
    plt.plot(t[::dec], dz_std[::dec], '--', color=C3, label=r'$\Delta z$, LIH')
    plt.axhline(0, color='gray', linestyle=':'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('X-point control error [mm]'); plt.legend(ncol=2)

    plt.subplot(3, 2, 2)
    plt.plot(t[::dec], enorm_net[::dec], color=C1, label=rf'Surrogate ({rms_net:.2f} mm)')
    plt.plot(t[::dec], enorm_std[::dec], '--', color=C3, label=rf'LIH ({rms_std:.2f} mm)')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'Euclidean error $\|e_X\|_2$ [mm]'); plt.legend(title='Full-window RMS')

    plt.subplot(3, 2, 3)
    plt.plot(t[::dec], obs_dr_net[::dec], color=C1, label=r'$\delta r$, surrogate')
    plt.plot(t[::dec], obs_dr_lih[::dec], '--', color=C1, label=r'$\delta r$, LIH shadow')
    plt.plot(t[::dec], obs_dz_net[::dec], color=C3, label=r'$\delta z$, surrogate')
    plt.plot(t[::dec], obs_dz_lih[::dec], '--', color=C3, label=r'$\delta z$, LIH shadow')
    plt.axhline(0, color='gray', linestyle=':'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Estimated X-point displacement [mm]'); plt.legend(ncol=2)

    plt.subplot(3, 2, 4)
    plt.plot(t[::dec], dr_y[::dec], color=C1, label=r'$r_Y$')
    plt.plot(t[::dec], dz_y[::dec], color=C3, label=r'$z_Y$')
    plt.axhline(0, color='gray', linestyle=':'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Surrogate-run minus LIH-run [mm]'); plt.legend(title='Current-centroid position')

    plt.subplot(3, 2, 5)
    plt.plot(t[::dec], dip[::dec], color=C2)
    plt.axhline(0, color='gray', linestyle=':'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'$I_p^{\mathrm{sur}}-I_p^{\mathrm{LIH}}$ [kA]')

    plt.subplot(3, 2, 6)
    plt.plot(t[::dec], va_norm_net[::dec], color=C1, label=r'$\|V_a^{\mathrm{sur}}\|_2$')
    plt.plot(t[::dec], va_norm_std[::dec], '--', color=C3, label=r'$\|V_a^{\mathrm{LIH}}\|_2$')
    plt.plot(t[::dec], va_diff[::dec], ':', color=C2, label=r'$\|V_a^{\mathrm{sur}}-V_a^{\mathrm{LIH}}\|_2$')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Coil-voltage vector norm [V]'); plt.legend()

    save_figure('fmag_control_simulation.pdf')
