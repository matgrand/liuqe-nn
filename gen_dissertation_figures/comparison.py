r"""Compare surrogate, warm- and cold-start LIUQE-RT, and LIH errors against offline LIUQE during TCV discharge 88564.

The three vertically stacked panels show the absolute poloidal-flux, radial-field, and vertical-field errors at 25 reference control points from $0.2$ to $1.6\,\mathrm{s}$. The LiuqeNN surrogate trace was evaluated in real time during the discharge, whereas warm-start LIUQE-RT, cold-start LIUQE-RT, and LIH were evaluated in post-processing; offline LIUQE supplies the common reference. Warm-start LIUQE-RT emulates the deployed observer parameters and advances from preceding equilibria, while cold-start LIUQE-RT independently initializes each reconstruction from offline-LIUQE inputs with one real-time iteration. Color identifies the observer, each solid line shows the moving-averaged spatial mean error, and its translucent band extends to the equivalently filtered maximum error across the 25 points. Logarithmic axes expose differences among the low-error surrogate and LIUQE-RT variants alongside the larger, comparatively steady LIH errors without merging the magnetic-field components.
"""
from matplotlib.lines import Line2D

from ._common import *


REALTIME_DATA_PATH = 'realtime_analysis/rt_shot88564.mat'
TIME_LIMITS = (0.2, 1.6)  # [s]
MOVING_AVERAGE_WINDOW = 11  # [samples], with 1 ms between samples in this trace.
FIELDS = ('Fx', 'Br', 'Bz')
YLABELS = (r'$\lvert\Delta\psi\rvert$ [Wb]', r'$\lvert\Delta B_r\rvert$ [T]', r'$\lvert\Delta B_z\rvert$ [T]')
OBSERVERS = (('online', 'Surrogate', C2), ('lrtW', 'LIUQE-RT warm', C1), ('lrtC', 'LIUQE-RT cold', C4), ('lih', 'LIH', C3))


def moving_average(values, window):
    """Return a centered moving average with edge-value padding."""
    if window < 1: raise ValueError(f'Moving-average window must be positive, got {window}.')
    if window == 1: return values
    left = window//2; right = window-1-left
    return np.convolve(np.pad(values, (left, right), mode='edge'), np.ones(window)/window, mode='valid')


def generate(moving_average_window=MOVING_AVERAGE_WINDOW):
    """Generate the three-panel comparison against offline LIUQE."""
    # %% ## Comparison with LIUQE, warm- and cold-start LIUQE-RT, and LIH
    ## Comparison with LIUQE, warm- and cold-start LIUQE-RT, and LIH
    # realtime_file_path = 'data/realtime_analysis_shot88903.mat'
    # realtime_file_path = 'data/realtime_analysis_shot88566.mat'
    # 't_or', 'phys_or', 'Fx_or', 'Br_or', 'Bz_or', 'Fx_online', 'Br_online', 'Bz_online', 'Fx_lrtW', 'Br_lrtW', 'Bz_lrtW', 'Fx_lrtC', 'Br_lrtC', 'Bz_lrtC', 'Fx_lih', 'Br_lih', 'Bz_lih', 'Fx_liu', 'Br_liu', 'Bz_liu', 'Fx_small', 'Br_small', 'Bz_small', 'Fx_small2', 'Br_small2', 'Bz_small2', 'Fx_turbo', 'Br_turbo', 'Bz_turbo', 'Fx_turbo2', 'Br_turbo2', 'Bz_turbo2', 't', 'phys', 'Fx_map', 'Iy_map', 't_map', 'Ip_map', 'r_map', 'z_map');
    data = loadmat(REALTIME_DATA_PATH); time = data['t'].reshape(-1)
    # t_or = m['t_or'].reshape(-1)
    # phys_or = m['phys_or']
    # orig = (m['Fx_or'], m['Br_or'], m['Bz_or'])
    # print(f'online: {online[0].shape}|{online[1].shape}, lrtW: {lrtW[0].shape}|{lrtW[1].shape}, lrtC: {lrtC[0].shape}|{lrtC[1].shape}, lih: {lih[0].shape}|{lih[1].shape}, liu: {liu[0].shape}|{liu[1].shape}')
    plt.figure(figsize=(THESIS_W, 5.5))
    for panel, (field, ylabel) in enumerate(zip(FIELDS, YLABELS), start=1):
        plt.subplot(3, 1, panel)
        reference = data[f'{field}_liu']
        for source, _, color in OBSERVERS:
            error = np.abs(data[f'{field}_{source}']-reference)
            mean_error, max_error = np.mean(error, axis=1), np.max(error, axis=1)
            mean_error = moving_average(mean_error, moving_average_window); max_error = moving_average(max_error, moving_average_window)
            plt.plot(time, mean_error, color=color)
            plt.fill_between(time, mean_error, max_error, color=color, alpha=0.3)
        plt.xlim(*TIME_LIMITS); plt.ylim(2e-4, 1e-1); plt.yscale('log'); plt.ylabel(ylabel)
        if panel == len(FIELDS): plt.xlabel(r'Time $t$ [s]')
        else: plt.tick_params(labelbottom=False)
    handles = [Line2D([], [], color=color, label=name) for _, name, color in OBSERVERS]
    plt.figlegend(handles=handles, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 0.985))
    # plt.suptitle('Error wrt Offline LIUQE')
    save_figure('comparison1.pdf', rect=(0, 0, 1, 0.955))


if __name__ == '__main__': generate(); compile_inspection_pdf()
