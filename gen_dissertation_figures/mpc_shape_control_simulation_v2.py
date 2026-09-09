from matplotlib.lines import Line2D

from ._common import *
from ._common import _cmap


def _required(data, *keys):
    missing = [key for key in keys if key not in data]
    if missing: raise KeyError(f'Missing MPC v2 fields: {missing}. Run mpc_shape_liuqenn_comparison_v2.m and copy both MAT files into data/.')


def generate():
    """Plot regulation, applied MPC moves, solver work, constraints, and geometry."""
    net = loadmat('data/mpc_shape_net_v2.mat')
    std = loadmat('data/mpc_shape_std_v2.mat')
    keys = ('t','controller_on','optimizer_on','shape_errs','active_shape_mask','qp_first_move','qp_iterations','mpc_update_mask','ief_upper_margin','ief_lower_margin','fge_boundary_r','fge_boundary_z','sc_t','sc_r','sc_z')
    _required(net, *keys); _required(std, *keys)
    t = net['t'].reshape(-1); t_std = std['t'].reshape(-1)
    if not np.allclose(t, t_std): raise ValueError('MPC v2 runs do not have matching time vectors.')
    ctrl_on = float(net['controller_on'].reshape(-1)[0]); opt_on = float(net['optimizer_on'].reshape(-1)[0]); dec, n_metric = 5, 10
    active = np.flatnonzero(net['active_shape_mask'].reshape(-1).astype(bool))
    metric = np.arange(min(n_metric, net['shape_errs'].shape[0])); post = t >= opt_on
    e_net = net['shape_errs'][active].T*1e3; e_std = std['shape_errs'][active].T*1e3
    e_metric_net = net['shape_errs'][metric].T*1e3; e_metric_std = std['shape_errs'][metric].T*1e3
    rms_net_t = np.sqrt(np.mean(e_metric_net**2,axis=1)); rms_std_t = np.sqrt(np.mean(e_metric_std**2,axis=1))
    rms_net = np.sqrt(np.mean(rms_net_t[post]**2)); rms_std = np.sqrt(np.mean(rms_std_t[post]**2))
    du_net = net['qp_first_move'].T; du_std = std['qp_first_move'].T
    score = np.sqrt(np.mean(du_net[post]**2,axis=0)+np.mean(du_std[post]**2,axis=0)); shown = np.sort(np.argsort(score)[-6:])
    updates = net['mpc_update_mask'].reshape(-1).astype(bool)
    margin_net = np.minimum(net['ief_upper_margin'],net['ief_lower_margin'])
    margin_std = np.minimum(std['ief_upper_margin'],std['ief_lower_margin'])

    plt.figure(figsize=(THESIS_W, 8.0))
    plt.subplot(3, 2, 1)
    colors = _cmap(np.linspace(0,1,len(active)))
    for i,color in enumerate(colors):
        plt.plot(t[::dec],e_net[::dec,i],color=color,label='Surrogate' if i==0 else None)
        plt.plot(t[::dec],e_std[::dec,i],'--',color=color,label='FGE/qintmex' if i==0 else None)
    plt.axvline(ctrl_on,color='gray',linestyle=':'); plt.axvline(opt_on,color='gray',linestyle='--'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Active shape errors [mm]'); plt.legend()

    plt.subplot(3, 2, 2)
    colors = _cmap(np.linspace(0,1,len(shown)))
    for i,j in enumerate(shown):
        plt.plot(t[::dec],du_net[::dec,j],color=colors[i],label=rf'$j={j+1}$')
        plt.plot(t[::dec],du_std[::dec,j],'--',color=colors[i])
    modes = plt.legend(title='Largest applied components',ncol=2,loc='upper right'); plt.gca().add_artist(modes)
    plt.legend(handles=[Line2D([0],[0],color='black',label='Surrogate'),Line2D([0],[0],color='black',linestyle='--',label='FGE/qintmex')],loc='lower right')
    plt.axvline(opt_on,color='gray',linestyle='--'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'First QP move $\Delta u_{k|k,j}$ [stored units]')

    plt.subplot(3, 2, 3)
    plt.plot(t[::dec],rms_net_t[::dec],color=C1,label=rf'Surrogate ({rms_net:.2f} mm)')
    plt.plot(t[::dec],rms_std_t[::dec],'--',color=C3,label=rf'FGE/qintmex ({rms_std:.2f} mm)')
    plt.axvline(opt_on,color='gray',linestyle='--'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Ten-channel RMS [mm]'); plt.legend(title=rf'RMS for $t\geq{opt_on:.2f}$ s')

    plt.subplot(3, 2, 4)
    plt.step(t[updates],net['qp_iterations'].reshape(-1)[updates],where='post',color=C1,label='Surrogate')
    plt.step(t[updates],std['qp_iterations'].reshape(-1)[updates],where='post',linestyle='--',color=C3,label='FGE/qintmex')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel('QP solver iterations'); plt.legend(title='Explicit 2 ms update mask')

    plt.subplot(3, 2, 5)
    plt.plot(t[::dec],np.min(margin_net,axis=0)[::dec]/1e3,color=C1,label='Surrogate')
    plt.plot(t[::dec],np.min(margin_std,axis=0)[::dec]/1e3,'--',color=C3,label='FGE/qintmex')
    plt.axhline(0,color='gray',linestyle=':'); plt.axvline(opt_on,color='gray',linestyle='--'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Minimum E/F-current margin [kA]'); plt.legend()

    plt.subplot(3, 2, 6)
    rb_net,zb_net = net['fge_boundary_r'][:,-1],net['fge_boundary_z'][:,-1]
    rb_std,zb_std = std['fge_boundary_r'][:,-1],std['fge_boundary_z'][:,-1]
    i_sc = np.argmin(np.abs(net['sc_t'].reshape(-1)-t[-1])); rc,zc = net['sc_r'][:,i_sc],net['sc_z'][:,i_sc]; valid = (rc>0) & net['active_shape_mask'].reshape(-1).astype(bool)
    plt.plot(rb_net,zb_net,color=C1,label='Surrogate-run boundary'); plt.plot(rb_std,zb_std,'--',color=C3,label='Baseline-run boundary')
    plt.plot(rc[valid],zc[valid],'x',color=C2,label='Active control points'); plt.xlabel(r'Radial coordinate $R$ [m]'); plt.ylabel(r'Vertical coordinate $Z$ [m]'); plt.axis('equal'); plt.legend()

    save_figure('mpc_shape_control_simulation_v2.pdf')
