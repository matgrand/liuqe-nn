# Figure description: This six-panel figure compares the same single-X-point FMag controller using either its standard LIH observer or NN, with the nonlinear FGE plant, references, timing, and all non-X-point observer paths held fixed. The top row shows the signed radial and vertical displacement errors perceived by each controller together with their Euclidean norm. Each central panel shows six traces: radial and vertical X-point displacements from LIH, NN, and FGE, evaluated on the NN-controlled trajectory at center-left and the LIH-controlled trajectory at center-right. The inactive observer is a passive shadow and does not affect that run's feedback. The bottom row compares the true FGE X-point distance from the target and the resulting coil-voltage demand in the two closed-loop runs. C1 denotes LIH, C2 NN, and C3 FGE throughout; solid and dashed lines denote radial and vertical quantities. The collector and plot validate the corrected E-F-G-OH neural coil order, the single fixed X-point query, the absence of boundary collocation points, and the inactive X-point-flux feedback channel.

from ._common import *


def _required(data, *keys):
    missing = [key for key in keys if key not in data]
    if missing: raise KeyError(f'Missing FMag fields: {missing}. Run fmag_liuqenn_comparison.m and copy both MAT files into data/.')


def _first(data, key): return np.asarray(data[key]).reshape(-1, data['t'].size)[0]


def _hessian_condition(data, suffix=''):
    new = (f'xpt_d2Fdr2{suffix}',f'xpt_d2Fdrdz{suffix}',f'xpt_d2Fdz2{suffix}')
    old = (f'xpt_raw_d2Fdrr{suffix}',f'xpt_raw_d2Fdrdz{suffix}',f'xpt_raw_d2Fdzz{suffix}')
    names = new if all(name in data for name in new) else old
    drr,drz,dzz = (np.asarray(data[name]).reshape(-1) for name in names)
    discriminant = np.sqrt((drr-dzz)**2+4*drz**2)
    eigenvalue_1 = 0.5*(drr+dzz+discriminant); eigenvalue_2 = 0.5*(drr+dzz-discriminant)
    largest = np.maximum(np.abs(eigenvalue_1),np.abs(eigenvalue_2)); smallest = np.minimum(np.abs(eigenvalue_1),np.abs(eigenvalue_2))
    return np.divide(largest,smallest,out=np.full_like(largest,np.inf),where=smallest>0)


def _validate(data, name):
    version = int(np.asarray(data.get('schema_version',0)).item())
    if version < 5: raise ValueError(f'{name} does not contain both LIH and NN shadows. Regenerate it with the current fmag_liuqenn_comparison.m.')
    permutation = np.asarray(data['nn_coil_permutation']).reshape(-1).astype(int)
    if not np.array_equal(permutation,np.r_[np.arange(1,17),19,17,18]): raise ValueError(f'{name} has unexpected NN coil order: {permutation}.')
    if np.count_nonzero(data['is_boundary_flux_observer']): raise ValueError(f'{name} contains boundary collocation points that do not affect this controller.')
    if 'observer_query_r' in data and np.asarray(data['observer_query_r']).size != 1: raise ValueError(f'{name} must contain exactly one neural query point.')
    flux_mask = np.asarray(data['nonlin_obs_type']).reshape(-1)==26
    if np.any(np.asarray(data['Tpre'])[:,flux_mask,:]): raise ValueError(f'{name} unexpectedly feeds X-point flux into the controller.')
    if np.any(data['controller_error_code']) or np.any(data['fmag_checks_error_code']): raise ValueError(f'{name} reports an FMag error.')
    if not np.all(data['fge_isconverged']): raise ValueError(f'{name} contains unconverged FGE samples.')


def generate():
    """Compare perceived control error, observer estimates, physical X-point regulation, and actuator demand."""
    net = loadmat('data/fmag_net.mat'); std = loadmat('data/fmag_std.mat')
    keys = ('schema_version','t','ip_scaling_reference','xpt_err_delta_r','xpt_err_delta_z','xpt_delta_r','xpt_delta_z','xpt_delta_r_lih','xpt_delta_z_lih','xpt_delta_r_nn','xpt_delta_z_nn','design_xpoint_r','design_xpoint_z','fge_active_xpoint_r','fge_active_xpoint_z','coil_voltages','nn_coil_permutation','is_boundary_flux_observer','nonlin_obs_type','Tpre','controller_error_code','fmag_checks_error_code','fge_isconverged')
    _required(net,*keys); _required(std,*keys); _validate(net,'fmag_net.mat'); _validate(std,'fmag_std.mat')
    t = np.asarray(net['t']).reshape(-1); t_std = np.asarray(std['t']).reshape(-1)
    if not np.allclose(t,t_std): raise ValueError('FMag runs do not have matching time vectors.')

    dec = 5; ip_net = np.asarray(net['ip_scaling_reference']).reshape(-1); ip_std = np.asarray(std['ip_scaling_reference']).reshape(-1)
    er_net = 1e3*np.asarray(net['xpt_err_delta_r']).reshape(-1)/ip_net; er_std = 1e3*np.asarray(std['xpt_err_delta_r']).reshape(-1)/ip_std
    ez_net = 1e3*np.asarray(net['xpt_err_delta_z']).reshape(-1)/ip_net; ez_std = 1e3*np.asarray(std['xpt_err_delta_z']).reshape(-1)/ip_std
    perceived_net = np.hypot(er_net,ez_net); perceived_std = np.hypot(er_std,ez_std)
    r0 = float(np.asarray(net['design_xpoint_r']).item()); z0 = float(np.asarray(net['design_xpoint_z']).item())
    fge_r_net = 1e3*(_first(net,'fge_active_xpoint_r')-r0); fge_z_net = 1e3*(_first(net,'fge_active_xpoint_z')-z0)
    fge_r_std = 1e3*(_first(std,'fge_active_xpoint_r')-r0); fge_z_std = 1e3*(_first(std,'fge_active_xpoint_z')-z0)
    nn_r = 1e3*np.asarray(net['xpt_delta_r_nn']).reshape(-1); nn_z = 1e3*np.asarray(net['xpt_delta_z_nn']).reshape(-1)
    lih_net_r = 1e3*np.asarray(net['xpt_delta_r_lih']).reshape(-1); lih_net_z = 1e3*np.asarray(net['xpt_delta_z_lih']).reshape(-1)
    lih_r = 1e3*np.asarray(std['xpt_delta_r_lih']).reshape(-1); lih_z = 1e3*np.asarray(std['xpt_delta_z_lih']).reshape(-1)
    nn_std_r = 1e3*np.asarray(std['xpt_delta_r_nn']).reshape(-1); nn_std_z = 1e3*np.asarray(std['xpt_delta_z_nn']).reshape(-1)
    physical_net = np.hypot(fge_r_net,fge_z_net); physical_std = np.hypot(fge_r_std,fge_z_std)
    observer_error_nn = np.hypot(nn_r-fge_r_net,nn_z-fge_z_net); observer_error_lih = np.hypot(lih_r-fge_r_std,lih_z-fge_z_std)
    va_net = np.linalg.norm(net['coil_voltages'],axis=0); va_std = np.linalg.norm(std['coil_voltages'],axis=0)
    central_values = np.concatenate((lih_net_r,lih_net_z,nn_r,nn_z,fge_r_net,fge_z_net,lih_r,lih_z,nn_std_r,nn_std_z,fge_r_std,fge_z_std))
    central_min,central_max = np.min(central_values),np.max(central_values); central_pad = 0.05*(central_max-central_min); central_ylim = (central_min-central_pad,central_max+central_pad)

    plt.figure(figsize=(THESIS_W,8.0))
    plt.subplot(3,2,1)
    plt.plot(t[::dec],er_std[::dec],color=C1,label=r'LIH, $e_R$'); plt.plot(t[::dec],ez_std[::dec],'--',color=C1,label=r'LIH, $e_Z$')
    plt.plot(t[::dec],er_net[::dec],color=C2,label=r'NN, $e_R$'); plt.plot(t[::dec],ez_net[::dec],'--',color=C2,label=r'NN, $e_Z$')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel('Perceived control error [mm]'); plt.legend(ncol=2)

    plt.subplot(3,2,2)
    plt.plot(t[::dec],perceived_std[::dec],color=C1,label=rf'LIH ({np.sqrt(np.mean(perceived_std**2)):.2f} mm)')
    plt.plot(t[::dec],perceived_net[::dec],color=C2,label=rf'NN ({np.sqrt(np.mean(perceived_net**2)):.2f} mm)')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'Perceived error $\|\bm e_X\|_2$ [mm]'); plt.legend(title='Control RMS')

    plt.subplot(3,2,3)
    plt.plot(t[::dec],lih_net_r[::dec],color=C1,label=r'LIH, $\delta R_X$'); plt.plot(t[::dec],lih_net_z[::dec],'--',color=C1,label=r'LIH, $\delta Z_X$')
    plt.plot(t[::dec],nn_r[::dec],color=C2,label=r'NN, $\delta R_X$'); plt.plot(t[::dec],nn_z[::dec],'--',color=C2,label=r'NN, $\delta Z_X$')
    plt.plot(t[::dec],fge_r_net[::dec],color=C3,label=r'FGE, $\delta R_X$'); plt.plot(t[::dec],fge_z_net[::dec],'--',color=C3,label=r'FGE, $\delta Z_X$')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'X-point displacement [mm]'); plt.ylim(central_ylim); plt.legend(ncol=2,title='NN control')

    plt.subplot(3,2,4)
    plt.plot(t[::dec],lih_r[::dec],color=C1,label=r'LIH, $\delta R_X$'); plt.plot(t[::dec],lih_z[::dec],'--',color=C1,label=r'LIH, $\delta Z_X$')
    plt.plot(t[::dec],nn_std_r[::dec],color=C2,label=r'NN, $\delta R_X$'); plt.plot(t[::dec],nn_std_z[::dec],'--',color=C2,label=r'NN, $\delta Z_X$')
    plt.plot(t[::dec],fge_r_std[::dec],color=C3,label=r'FGE, $\delta R_X$'); plt.plot(t[::dec],fge_z_std[::dec],'--',color=C3,label=r'FGE, $\delta Z_X$')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'X-point displacement [mm]'); plt.ylim(central_ylim); plt.legend(ncol=2,title='LIH control')

    plt.subplot(3,2,5)
    plt.plot(t[::dec],physical_std[::dec],color=C1,label=rf'LIH control ({np.sqrt(np.mean(physical_std**2)):.2f} mm)')
    plt.plot(t[::dec],physical_net[::dec],color=C2,label=rf'NN control ({np.sqrt(np.mean(physical_net**2)):.2f} mm)')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'Physical error $\|\bm x_X-\bm x_X^{\mathrm{ref}}\|_2$ [mm]'); plt.legend(title='Physical RMS')

    plt.subplot(3,2,6)
    plt.plot(t[::dec],va_std[::dec],color=C1,label='LIH'); plt.plot(t[::dec],va_net[::dec],color=C2,label='NN')
    plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'Coil-voltage vector norm [V]'); plt.legend()

    save_figure('fmag_control_simulation.pdf')
    condition_nn = _hessian_condition(net); condition_lih = _hessian_condition(std)
    print(f'FMag active-observer error RMS: NN={np.sqrt(np.mean(observer_error_nn**2)):.2f} mm, LIH={np.sqrt(np.mean(observer_error_lih**2)):.2f} mm')
    print(f'FMag Hessian condition range: NN={np.min(condition_nn):.3f}-{np.max(condition_nn):.3f}, LIH={np.min(condition_lih):.3f}-{np.max(condition_lih):.3f}')
    print(f'FMag voltage-limit samples: NN={np.count_nonzero(net.get("voltage_at_limit",0))}, LIH={np.count_nonzero(std.get("voltage_at_limit",0))}')


if __name__ == '__main__': generate(); compile_inspection_pdf()
