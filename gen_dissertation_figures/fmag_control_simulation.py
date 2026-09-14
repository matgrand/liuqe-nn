from ._common import *


def _required(data, *keys):
    missing = [key for key in keys if key not in data]
    if missing: raise KeyError(f'Missing FMag fields: {missing}. Run fmag_liuqenn_comparison.m and copy both MAT files into data/.')


def _validate_schema(data):
    version = int(np.asarray(data.get('schema_version', 0)).reshape(-1)[0])
    if version < 3: raise ValueError('FMag MAT data predates the LiuqeNN E-F-G-OH coil-order fix. Rerun fmag_liuqenn_comparison.m.')
    permutation = np.asarray(data['nn_coil_permutation']).reshape(-1).astype(int)
    if not np.array_equal(permutation, np.r_[np.arange(1,17),19,17,18]): raise ValueError(f'Unexpected LiuqeNN coil permutation: {permutation}.')


def _hessian_condition(data, suffix=''):
    drr = data[f'xpt_raw_d2Fdrr{suffix}'].reshape(-1)
    drz = data[f'xpt_raw_d2Fdrdz{suffix}'].reshape(-1)
    dzz = data[f'xpt_raw_d2Fdzz{suffix}'].reshape(-1)
    discriminant = np.sqrt((drr-dzz)**2 + 4*drz**2)
    eigenvalue_1 = 0.5*(drr+dzz+discriminant)
    eigenvalue_2 = 0.5*(drr+dzz-discriminant)
    largest = np.maximum(np.abs(eigenvalue_1),np.abs(eigenvalue_2))
    smallest = np.minimum(np.abs(eigenvalue_1),np.abs(eigenvalue_2))
    return np.divide(largest,smallest,out=np.full_like(largest,np.inf),where=smallest>0)


def generate():
    """Plot FMag control errors, X-point displacement, observer accuracy, and Hessian conditioning for the neural and LIH observer configurations."""
    net = loadmat('data/fmag_net.mat')
    std = loadmat('data/fmag_std.mat')
    _validate_schema(net); _validate_schema(std)
    keys = ('t','ip_scaling_reference','xpt_err_delta_r','xpt_err_delta_z','xpt_delta_r','xpt_delta_z','xpt_delta_r_lih','xpt_delta_z_lih','design_xpoint_r','design_xpoint_z','fge_active_xpoint_r','fge_active_xpoint_z','xpt_raw_d2Fdrr','xpt_raw_d2Fdrdz','xpt_raw_d2Fdzz','xpt_raw_d2Fdrr_lih','xpt_raw_d2Fdrdz_lih','xpt_raw_d2Fdzz_lih','coil_voltages')
    _required(net,*keys); _required(std,*keys)
    t = net['t'].reshape(-1); t_std = std['t'].reshape(-1)
    if not np.allclose(t,t_std): raise ValueError('FMag runs do not have matching time vectors.')
    dec = 5; ip_net = net['ip_scaling_reference'].reshape(-1); ip_std = std['ip_scaling_reference'].reshape(-1)
    dr_net = 1e3*net['xpt_err_delta_r'].reshape(-1)/ip_net; dr_std = 1e3*std['xpt_err_delta_r'].reshape(-1)/ip_std
    dz_net = 1e3*net['xpt_err_delta_z'].reshape(-1)/ip_net; dz_std = 1e3*std['xpt_err_delta_z'].reshape(-1)/ip_std
    ctrl_norm_net = np.hypot(dr_net,dz_net); ctrl_norm_std = np.hypot(dr_std,dz_std)
    r0 = float(net['design_xpoint_r'].reshape(-1)[0]); z0 = float(net['design_xpoint_z'].reshape(-1)[0])
    true_dr = 1e3*(net['fge_active_xpoint_r'][0]-r0); true_dz = 1e3*(net['fge_active_xpoint_z'][0]-z0)
    obs_dr = 1e3*net['xpt_delta_r'].reshape(-1); obs_dz = 1e3*net['xpt_delta_z'].reshape(-1)
    lih_dr = 1e3*net['xpt_delta_r_lih'].reshape(-1); lih_dz = 1e3*net['xpt_delta_z_lih'].reshape(-1)
    obs_err = np.hypot(obs_dr-true_dr,obs_dz-true_dz); lih_err = np.hypot(lih_dr-true_dr,lih_dz-true_dz)
    condition_net = _hessian_condition(net); condition_lih = _hessian_condition(net,'_lih')
    va_net = net['coil_voltages']; va_std = std['coil_voltages']

    plt.figure(figsize=(THESIS_W, 8.0))
    plt.subplot(3, 2, 1)
    plt.plot(t[::dec],dr_net[::dec],color=C1,label=r'$\Delta r$, surrogate'); plt.plot(t[::dec],dr_std[::dec],'--',color=C1,label=r'$\Delta r$, LIH')
    plt.plot(t[::dec],dz_net[::dec],color=C2,label=r'$\Delta z$, surrogate'); plt.plot(t[::dec],dz_std[::dec],'--',color=C2,label=r'$\Delta z$, LIH')
    plt.axhline(0,color='gray',linestyle=':'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel('X-point control error [mm]'); plt.legend(ncol=2)

    plt.subplot(3, 2, 2)
    plt.plot(t[::dec],ctrl_norm_net[::dec],color=C1,label=rf'Surrogate ({np.sqrt(np.mean(ctrl_norm_net**2)):.2f} mm)')
    plt.plot(t[::dec],ctrl_norm_std[::dec],'--',color=C2,label=rf'LIH ({np.sqrt(np.mean(ctrl_norm_std**2)):.2f} mm)')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'Euclidean control error $\|e_X\|_2$ [mm]'); plt.legend(title='Full-window RMS')

    plt.subplot(3, 2, 3)
    plt.plot(t[::dec],true_dr[::dec],color=C1,label='FGE truth'); plt.plot(t[::dec],obs_dr[::dec],color=C2,label='Surrogate')
    plt.plot(t[::dec],lih_dr[::dec],'--',color=C3,label='LIH shadow'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'X-point radial displacement $\delta r$ [mm]'); plt.legend()

    plt.subplot(3, 2, 4)
    plt.plot(t[::dec],true_dz[::dec],color=C1,label='FGE truth'); plt.plot(t[::dec],obs_dz[::dec],color=C2,label='Surrogate')
    plt.plot(t[::dec],lih_dz[::dec],'--',color=C3,label='LIH shadow'); plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'X-point vertical displacement $\delta z$ [mm]'); plt.legend()

    plt.subplot(3, 2, 5)
    plt.plot(t[::dec],obs_err[::dec],color=C1,label='Surrogate'); plt.plot(t[::dec],lih_err[::dec],'--',color=C2,label='LIH shadow')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Observer error relative to FGE [mm]'); plt.legend()

    plt.subplot(3, 2, 6)
    plt.plot(t[::dec],condition_net[::dec],color=C1,label='Surrogate')
    plt.plot(t[::dec],condition_lih[::dec],'--',color=C2,label='LIH shadow')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'Hessian condition number $\kappa_2(H_\psi)$'); plt.legend()

    save_figure('fmag_control_simulation.pdf')
    print(f'FMag voltage-limit samples: surrogate={np.count_nonzero(net.get("voltage_at_limit",0))}, LIH={np.count_nonzero(std.get("voltage_at_limit",0))}')
    print(f'FMag maximum voltage-vector difference: {np.max(np.linalg.norm(va_net-va_std,axis=0)):.3f} V')


if __name__ == '__main__': generate(); compile_inspection_pdf()
