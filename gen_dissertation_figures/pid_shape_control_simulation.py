from matplotlib.lines import Line2D

from ._common import *
from ._common import _cmap


def _required(data, *keys):
    missing = [key for key in keys if key not in data]
    if missing: raise KeyError(f'Missing PID fields: {missing}. Run pid_shape_liuqenn_comparison.m and copy both MAT files into data/.')


def generate():
    """Plot the extended PID dataset produced by pid_shape_liuqenn_comparison.m."""
    net = loadmat('data/pid_shape_net.mat')
    std = loadmat('data/pid_shape_std.mat')
    geometry_keys = ('fge_grid_r', 'fge_grid_z', 'fge_final_flux', 'fge_final_boundary_flux', 'limiter_r', 'limiter_z')
    _required(net, 't', 'controller_on', 'shape_errs', 'proj_errs', 'active_shape_mask', 'active_projected_mask', 'sc_t', 'sc_r', 'sc_z', 'coil_voltages', *geometry_keys)
    _required(std, 't', 'shape_errs', 'proj_errs', 'coil_voltages', *geometry_keys)

    t = net['t'].reshape(-1); t_std = std['t'].reshape(-1)
    if not np.allclose(t, t_std): raise ValueError('PID runs do not have matching time vectors.')
    on = float(np.asarray(net['controller_on']).reshape(-1)[0]); dec, n_metric = 5, 10
    ish = np.flatnonzero(net['active_shape_mask'].reshape(-1).astype(bool))
    ipr = np.flatnonzero(net['active_projected_mask'].reshape(-1).astype(bool))
    e_net = net['shape_errs'][ish].T * 1e3; e_std = std['shape_errs'][ish].T * 1e3
    p_net = net['proj_errs'][ipr].T * 1e3; p_std = std['proj_errs'][ipr].T * 1e3
    metric = np.arange(min(n_metric, net['shape_errs'].shape[0])); post = t >= on
    rms_net_t = 1e3*np.sqrt(np.mean(net['shape_errs'][metric].T**2, axis=1))
    rms_std_t = 1e3*np.sqrt(np.mean(std['shape_errs'][metric].T**2, axis=1))
    rms_net = np.sqrt(np.mean(rms_net_t[post]**2)); rms_std = np.sqrt(np.mean(rms_std_t[post]**2))
    rms_net_ch = 1e3*np.sqrt(np.mean(net['shape_errs'][metric][:,post]**2, axis=1))
    rms_std_ch = 1e3*np.sqrt(np.mean(std['shape_errs'][metric][:,post]**2, axis=1))
    va_net = net['coil_voltages']; va_std = std['coil_voltages']

    plt.figure(figsize=(THESIS_W, 8.0))
    plt.subplot(3, 2, 1)
    colors = _cmap(np.linspace(0, 1, len(ish)))
    for i, color in enumerate(colors):
        plt.plot(t[::dec], e_net[::dec,i], color=color, label='Surrogate' if i == 0 else None)
        plt.plot(t[::dec], e_std[::dec,i], '--', color=color, label='FGE/qintmex' if i == 0 else None)
    plt.axvline(on, color='gray', linestyle=':'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Active shape errors [mm]'); plt.legend()

    plt.subplot(3, 2, 2)
    colors = _cmap(np.linspace(0, 1, len(ipr)))
    for i, color in enumerate(colors):
        plt.plot(t[::dec], p_net[::dec,i], color=color, label=rf'$j={ipr[i]+1}$')
        plt.plot(t[::dec], p_std[::dec,i], '--', color=color)
    plt.axvline(on, color='gray', linestyle=':'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Projected errors [mm]'); plt.legend(title='Active PID directions')

    plt.subplot(3, 2, 3)
    plt.plot(t[::dec], rms_net_t[::dec], color=C1, label=rf'Surrogate ({rms_net:.2f} mm)')
    plt.plot(t[::dec], rms_std_t[::dec], '--', color=C3, label=rf'FGE/qintmex ({rms_std:.2f} mm)')
    plt.axvline(on, color='gray', linestyle=':'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Ten-channel RMS [mm]'); plt.legend(title='Post-activation RMS')

    plt.subplot(3, 2, 4)
    channels = metric + 1
    plt.plot(channels, rms_net_ch, 'o-', color=C1, label='Surrogate')
    plt.plot(channels, rms_std_ch, 's--', color=C3, label='FGE/qintmex')
    plt.xticks(channels); plt.xlabel('Saved shape-error channel'); plt.ylabel('Post-activation RMS [mm]'); plt.legend()

    plt.subplot(3, 2, 5)
    grid_r = net['fge_grid_r'].reshape(-1); grid_z = net['fge_grid_z'].reshape(-1)
    flux_net = np.squeeze(net['fge_final_flux']); flux_std = np.squeeze(std['fge_final_flux'])
    boundary_flux_net = float(net['fge_final_boundary_flux'].reshape(-1)[0]); boundary_flux_std = float(std['fge_final_boundary_flux'].reshape(-1)[0])
    i_sc = np.argmin(np.abs(net['sc_t'].reshape(-1)-t[-1])); rc, zc = net['sc_r'][:,i_sc], net['sc_z'][:,i_sc]
    valid = (rc > 0) & net['active_shape_mask'].reshape(-1).astype(bool)
    plt.plot(net['limiter_r'].reshape(-1), net['limiter_z'].reshape(-1), color='gray', linewidth=0.8, label='Limiter')
    plt.contour(grid_r, grid_z, flux_net, levels=[boundary_flux_net], colors=[C1], linewidths=1.4)
    plt.contour(grid_r, grid_z, flux_std, levels=[boundary_flux_std], colors=[C3], linewidths=1.4, linestyles='--')
    plt.plot(rc[valid], zc[valid], 'x', color=C2, label='Control points')
    boundary_handles = [Line2D([0],[0],color=C1,label='Surrogate LCFS'), Line2D([0],[0],color=C3,linestyle='--',label='FGE/qintmex LCFS')]
    handles, labels = plt.gca().get_legend_handles_labels(); plt.xlabel(r'Radial coordinate $R$ [m]'); plt.ylabel(r'Vertical coordinate $Z$ [m]'); plt.axis('equal')
    geometry_axis = plt.gca(); limiter_r = net['limiter_r'].reshape(-1); limiter_z = net['limiter_z'].reshape(-1)
    vertical_span = np.ptp(limiter_z); horizontal_span = max(1.25*vertical_span, 1.10*np.ptp(limiter_r))
    shifted_left = np.min(limiter_r) - 0.04*horizontal_span
    geometry_axis.set_xlim(shifted_left, shifted_left+horizontal_span); geometry_axis.set_ylim(np.min(limiter_z)-0.04*vertical_span, np.max(limiter_z)+0.04*vertical_span)
    geometry_axis.set_aspect('equal', adjustable='box')
    geometry_axis.yaxis.tick_left(); geometry_axis.yaxis.set_label_position('left')
    plt.legend(boundary_handles + handles, [h.get_label() for h in boundary_handles] + labels, loc='center right')

    plt.subplot(3, 2, 6)
    plt.plot(t[::dec], np.linalg.norm(va_net,axis=0)[::dec], color=C1, label=r'$\|V_a^{\mathrm{sur}}\|_2$')
    plt.plot(t[::dec], np.linalg.norm(va_std,axis=0)[::dec], '--', color=C3, label=r'$\|V_a^{\mathrm{base}}\|_2$')
    plt.plot(t[::dec], np.linalg.norm(va_net-va_std,axis=0)[::dec], ':', color=C2, label=r'$\|\Delta V_a\|_2$')
    plt.axvline(on, color='gray', linestyle=':'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Coil-voltage vector norm [V]'); plt.legend()

    save_figure('pid_shape_control_simulation.pdf')
