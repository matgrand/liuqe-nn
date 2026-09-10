from matplotlib.lines import Line2D

from ._common import *
from ._common import _cmap


def _required(data, *keys):
    missing = [key for key in keys if key not in data]
    if missing: raise KeyError(f'Missing MPC fields: {missing}. Run mpc_shape_liuqenn_comparison.m and copy both MAT files into data/.')


def generate():
    """Compare surrogate and FGE/qintmex MPC regulation, moves, solver work, constraints, and geometry; use C1, C2, and C3 in series order where at most three colors are needed."""
    net = loadmat('data/mpc_shape_net.mat')
    std = loadmat('data/mpc_shape_std.mat')
    keys = ('t','controller_on','optimizer_on','shape_errs','active_shape_mask','qp_first_move','qp_iterations','mpc_update_mask','ief_upper_margin','ief_lower_margin','fge_grid_r','fge_grid_z','fge_final_flux','fge_final_boundary_flux','limiter_r','limiter_z','sc_t','sc_r','sc_z')
    _required(net, *keys); _required(std, *keys)
    t = net['t'].reshape(-1); t_std = std['t'].reshape(-1)
    if not np.allclose(t, t_std): raise ValueError('MPC runs do not have matching time vectors.')
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
    plt.plot(t[::dec],rms_std_t[::dec],'--',color=C2,label=rf'FGE/qintmex ({rms_std:.2f} mm)')
    plt.axvline(opt_on,color='gray',linestyle='--'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Ten-channel RMS [mm]'); plt.legend(title=rf'RMS for $t\geq{opt_on:.2f}$ s')

    plt.subplot(3, 2, 4)
    plt.step(t[updates],net['qp_iterations'].reshape(-1)[updates],where='post',color=C1,label='Surrogate')
    plt.step(t[updates],std['qp_iterations'].reshape(-1)[updates],where='post',linestyle='--',color=C2,label='FGE/qintmex')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel('QP solver iterations'); plt.legend(title='Explicit 2 ms update mask')

    plt.subplot(3, 2, 5)
    plt.plot(t[::dec],np.min(margin_net,axis=0)[::dec]/1e3,color=C1,label='Surrogate')
    plt.plot(t[::dec],np.min(margin_std,axis=0)[::dec]/1e3,'--',color=C2,label='FGE/qintmex')
    plt.axhline(0,color='gray',linestyle=':'); plt.axvline(opt_on,color='gray',linestyle='--'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Minimum E/F-current margin [kA]'); plt.legend()

    plt.subplot(3, 2, 6)
    grid_r=net['fge_grid_r'].reshape(-1); grid_z=net['fge_grid_z'].reshape(-1)
    flux_net=np.squeeze(net['fge_final_flux']); flux_std=np.squeeze(std['fge_final_flux'])
    boundary_flux_net=float(net['fge_final_boundary_flux'].reshape(-1)[0]); boundary_flux_std=float(std['fge_final_boundary_flux'].reshape(-1)[0])
    i_sc = np.argmin(np.abs(net['sc_t'].reshape(-1)-t[-1])); rc,zc = net['sc_r'][:,i_sc],net['sc_z'][:,i_sc]; valid = (rc>0) & net['active_shape_mask'].reshape(-1).astype(bool)
    plt.plot(net['limiter_r'].reshape(-1),net['limiter_z'].reshape(-1),color='gray',linewidth=0.8,label='Limiter')
    plt.contour(grid_r,grid_z,flux_net,levels=[boundary_flux_net],colors=[C1],linewidths=1.4)
    plt.contour(grid_r,grid_z,flux_std,levels=[boundary_flux_std],colors=[C2],linewidths=1.4,linestyles='--')
    plt.plot(rc[valid],zc[valid],'x',color=C3,label='Control points')
    boundary_handles=[Line2D([0],[0],color=C1,label='Surrogate LCFS'),Line2D([0],[0],color=C2,linestyle='--',label='FGE/qintmex LCFS')]
    handles,labels=plt.gca().get_legend_handles_labels(); plt.xlabel(r'Radial coordinate $R$ [m]'); plt.ylabel(r'Vertical coordinate $Z$ [m]'); plt.axis('equal')
    geometry_axis=plt.gca(); limiter_r=net['limiter_r'].reshape(-1); limiter_z=net['limiter_z'].reshape(-1)
    vertical_span=np.ptp(limiter_z); horizontal_span=max(1.25*vertical_span,1.10*np.ptp(limiter_r))
    shifted_left=np.min(limiter_r)-0.04*horizontal_span
    geometry_axis.set_xlim(shifted_left,shifted_left+horizontal_span); geometry_axis.set_ylim(np.min(limiter_z)-0.04*vertical_span,np.max(limiter_z)+0.04*vertical_span)
    geometry_axis.set_aspect('equal',adjustable='box')
    geometry_axis.yaxis.tick_left(); geometry_axis.yaxis.set_label_position('left')
    plt.legend(boundary_handles+handles,[h.get_label() for h in boundary_handles]+labels,loc='center right')

    save_figure('mpc_shape_control_simulation.pdf')


if __name__ == '__main__': generate(); compile_inspection_pdf()