from ._common import *
from ._common import _cmap


def generate():
    # %% ## Control simulation (old version)
    ## Control simulation
    CTRL_SIM_SMALL = 'data/sim_res_small.mat'
    CTRL_SIM_TURBO = 'data/sim_res_turbo.mat'
    CTRL_SIM_STD = 'data/sim_res_std.mat'
    # 'Br', 'Br2', 'Bz', 'Bz2', 'Fx', 'Fx2', 'currs', 'proj_errs', 'rz', 'shape_errs'

    ctrl_start = 0.6 # controller start time
    vlinec = 'gray'
    vlines = '--'

    sim_net = loadmat(CTRL_SIM_SMALL)
    # sim_net = loadmat(CTRL_SIM_TURBO)
    sim_std = loadmat(CTRL_SIM_STD)

    FxN, BrN, BzN, IaN, peN, seN = sim_net['Fx'], sim_net['Br'], sim_net['Bz'], sim_net['currs'], sim_net['proj_errs'], sim_net['shape_errs']
    FxL, BrL, BzL, IaL, peL, seL = sim_std['Fx'], sim_std['Br'], sim_std['Bz'], sim_std['currs'], sim_std['proj_errs'], sim_std['shape_errs']
    FxNL, BrNL, BzNL, IaNL, peNL, seNL = sim_net['Fx2'], sim_net['Br2'], sim_net['Bz2'], sim_net['currs'], sim_net['proj_errs2'], sim_net['shape_errs2']
    t = loadmat(CTRL_SIM_SMALL)['t'].reshape(-1)

    # print shapes of all the variables
    print(f'FxN:{FxN.shape}, BrN:{BrN.shape}, BzN:{BzN.shape}, IaN:{IaN.shape}, peN:{peN.shape}, seN:{seN.shape}')
    print(f'FxL:{FxL.shape}, BrL:{BrL.shape}, BzL:{BzL.shape}, IaL:{IaL.shape}, peL:{peL.shape}, seL:{seL.shape}')
    print(f'FxNL:{FxNL.shape}, BrNL:{BrNL.shape}, BzNL:{BzNL.shape}, t:{t.shape}')

    # analyze Ia, pe, se, to remove the 0 raws
    μIaN, μpeN, μseN = np.mean(IaN, axis=1), np.mean(peN, axis=1), np.mean(seN, axis=1)
    μIaL, μpeL, μseL = np.mean(IaL, axis=1), np.mean(peL, axis=1), np.mean(seL, axis=1)
    idxIaN = np.where(np.abs(μIaN) > 0)[0]
    idxpeN = np.where(np.abs(μpeN) > 0)[0]
    idxseN = np.where(np.abs(μseN) > 0)[0]
    idxIaL = np.where(np.abs(μIaL) > 0)[0]
    idxpeL = np.where(np.abs(μpeL) > 0)[0]
    idxseL = np.where(np.abs(μseL) > 0)[0]
    assert np.array_equal(idxIaN, idxIaL), f'idxIaN: {idxIaN}, idxIaL: {idxIaL}'
    assert np.array_equal(idxpeN, idxpeL), f'idxpeN: {idxpeN}, idxpeL: {idxpeL}'
    assert np.array_equal(idxseN, idxseL), f'idxseN: {idxseN}, idxseL: {idxseL}'
    print(f'Non-zero rows indices: Ia: {idxIaN}, pe: {idxpeN}, se: {idxseN}')
    IaN, peN, seN, IaL, peL, seL = IaN[idxIaN], peN[idxpeN], seN[idxseN], IaL[idxIaL], peL[idxpeL], seL[idxseL]

    # assert False

    dec = 5
    n_ctrl_pts = 10
    t = t[::dec]
    FxN, BrN, BzN, FxL, BrL, BzL, FxNL, BrNL, BzNL = map(lambda x: x[:n_ctrl_pts,::dec].T, [FxN, BrN, BzN, FxL, BrL, BzL, FxNL, BrNL, BzNL])
    IaN, peN, seN, IaL, peL, seL, IaNL, peNL, seNL = map(lambda x: x[:,::dec].T, [IaN, peN, seN, IaL, peL, seL, IaNL, peNL, seNL])

    # # error wrt liuqe on its sim
    # eFx = np.abs(FxN - FxL)
    # eBr = np.abs(BrN - BrL)
    # eBz = np.abs(BzN - BzL)
    # error wrt liuqe on same sim
    eFx = np.abs(FxN - FxNL)
    eBr = np.abs(BrN - BrNL)
    eBz = np.abs(BzN - BzNL)

    color_cycle = _cmap(np.linspace(0, 1, n_ctrl_pts))
    α = 0.7
    # spl, fs = (2, 3), (8, 6)
    # spl, fs = (3, 2), (6, 8)
    spl, fs = (3, 2), (THESIS_W, 7.2)

    plt.figure(figsize=fs)
    std_style = '--'
    plt.subplot(*spl, 1) # top: Fx Br Bz. bot: se, pe, Ia
    for cp in reversed(range(n_ctrl_pts)):
        zero_line = 0.0025 * cp * np.ones_like(t)
        er = eFx[:, cp] + zero_line
        plt.fill_between(t, zero_line, er, alpha=α, color=color_cycle[cp])
    plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
    plt.ylabel(r'Stacked $\lvert\Delta\psi\rvert$ [Wb]')
    plt.subplot(*spl, 2)
    for cp in reversed(range(n_ctrl_pts)):
        zero_line = 0.0025 * cp * np.ones_like(t)
        er = eBr[:, cp] + zero_line
        plt.fill_between(t, zero_line, er, alpha=α, color=color_cycle[cp])
    plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
    plt.ylabel(r'Stacked $\lvert\Delta B_r\rvert$ [T]')
    plt.subplot(*spl, 3)
    for cp in reversed(range(n_ctrl_pts)):
        zero_line = 0.0025 * cp * np.ones_like(t)
        er = eBz[:, cp] + zero_line
        plt.fill_between(t, zero_line, er, alpha=α, color=color_cycle[cp])
    plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
    plt.ylabel(r'Stacked $\lvert\Delta B_z\rvert$ [T]')

    plt.subplot(*spl, 4)
    colors = _cmap(np.linspace(0, 1, seN.shape[1]))
    for i in range(seN.shape[1]):
        plt.plot(t, seN[:, i], color=colors[i], label='Surrogate' if i==0 else None)
        plt.plot(t, seL[:, i], std_style, color=colors[i], label='LIUQE-RT' if i==0 else None)
    plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
    plt.ylabel(r'Shape error [m]'), plt.legend()

    plt.subplot(*spl, 5)
    colors = _cmap(np.linspace(0, 1, peN.shape[1]))
    for i in range(peN.shape[1]):
        plt.plot(t, peN[:, i], color=colors[i])
        plt.plot(t, peL[:, i], std_style, color=colors[i])
    plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
    plt.xlabel(r'Time $t$ [s]'), plt.ylabel(r'Projection error [m]')

    plt.subplot(*spl, 6)
    colors = _cmap(np.linspace(0, 1, IaN.shape[1]))
    for i in range(IaN.shape[1]):
        plt.plot(t, IaN[:, i], color=colors[i])
        plt.plot(t, IaL[:, i], std_style, color=colors[i])
    plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
    plt.xlabel(r'Time $t$ [s]'), plt.ylabel(r'Coil current [kA]')

    save_figure('control_simulation.pdf')
