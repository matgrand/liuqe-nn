from ._common import *
from ._common import _cmap
from . import _common


def generate():
    # %% ## Comparison with LIUQE LIUQE-RT and LIH
    ## Comparison with LIUQE LIUQE-RT and LIH
    # realtime_file_path = 'data/realtime_analysis_shot88903.mat'
    # realtime_file_path = 'data/realtime_analysis_shot88566.mat'
    realtime_file_path = 'realtime_analysis/rt_shot88564.mat'

    # 't_or', 'phys_or', 'Fx_or', 'Br_or', 'Bz_or', 'Fx_online', 'Br_online', 'Bz_online', 'Fx_lrt', 'Br_lrt', 'Bz_lrt', 'Fx_lih', 'Br_lih', 'Bz_lih', 'Fx_liu', 'Br_liu', 'Bz_liu', 'Fx_small', 'Br_small', 'Bz_small', 'Fx_small2', 'Br_small2', 'Bz_small2', 'Fx_turbo', 'Br_turbo', 'Bz_turbo', 'Fx_turbo2', 'Br_turbo2', 'Bz_turbo2', 't', 'phys', 'Fx_map', 'Iy_map', 't_map', 'Ip_map', 'r_map', 'z_map');

    m = loadmat(realtime_file_path)

    online = (m['Fx_online'], np.sqrt(m['Br_online']**2 + m['Bz_online']**2))
    lrt = (m['Fx_lrt'], np.sqrt(m['Br_lrt']**2 + m['Bz_lrt']**2))
    lih = (m['Fx_lih'], np.sqrt(m['Br_lih']**2 + m['Bz_lih']**2))
    liu = (m['Fx_liu'], np.sqrt(m['Br_liu']**2 + m['Bz_liu']**2))

    t, phys = m['t'].reshape(-1), m['phys']

    # t_or = m['t_or'].reshape(-1)
    # phys_or = m['phys_or']
    # orig = (m['Fx_or'], m['Br_or'], m['Bz_or'])

    # print(f'online: {online[0].shape}|{online[1].shape}, lrt: {lrt[0].shape}|{lrt[1].shape}, lih: {lih[0].shape}|{lih[1].shape}, liu: {liu[0].shape}|{liu[1].shape}')

    eonline = list(map(lambda a, b: np.abs(a-b), online, liu))
    elrt = list(map(lambda a, b: np.abs(a-b), lrt, liu))
    elih = list(map(lambda a, b: np.abs(a-b), lih, liu))

    plt.figure(figsize=(THESIS_W, 4.6))
    cols = [C2, C1, C3]
    titles = [r'Poloidal-flux error $\lvert\Delta\psi\rvert$ [Wb]', r'Magnetic-field error $\lvert\Delta B\rvert$ [T]']
    for k in range(len(online)):
        plt.subplot(len(online), 1, k+1)
        for i, (errs, name) in enumerate(zip([eonline, elrt, elih], ['Surrogate', 'LIUQE-RT', 'LIH'])):
            min_err, max_err, avg_err = np.min(errs[k], axis=1), np.max(errs[k], axis=1), np.mean(errs[k], axis=1)
            plt.plot(t, avg_err, color=cols[i], label=name)
            plt.fill_between(t, avg_err, max_err, color=cols[i], alpha=0.3)
        plt.ylim(2e-4, 1e-1)
        plt.yscale('log')
        plt.ylabel(titles[k])
        plt.legend(loc='upper right')
        if k==len(online)-1: plt.xlabel(r'Time $t$ [s]')
    # plt.suptitle('Error wrt Offline LIUQE')
    save_figure('comparison1.pdf')

    if _common.SHOW_PLOTS: plt.show()

if __name__ == '__main__': generate(); compile_inspection_pdf()