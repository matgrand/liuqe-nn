# %% ## Script for the dissertation figures
## Script for the dissertation figures
from utils import *
FIGS_DIR = 'dissertation_figures'
THESIS_W = 6.30  # [in]
SHOW_PLOTS = True
WAIT_FOR_PLOTS = False
os.makedirs(FIGS_DIR, exist_ok=True)

# plt.style.use('dissertation-thesis')
plt.style.use('dissertation-mono')

def save_figure(filename, rect=None):
    path = os.path.join(FIGS_DIR, filename)
    plt.tight_layout(rect=rect); plt.savefig(path); print(f'Saved {path}')
    if SHOW_PLOTS: plt.show(block=WAIT_FOR_PLOTS)
    else: plt.close()


from matplotlib.colors import LinearSegmentedColormap
# _cmap = LinearSegmentedColormap.from_list('my_custom_colormap', plt.cm.inferno(np.linspace(0.1, 0.85, 256)))
_base_cmap = plt.get_cmap(plt.rcParams['image.cmap'])
_cmap = LinearSegmentedColormap.from_list('my_custom_colormap', _base_cmap(np.linspace(0.05, 0.9, 256)))
# _cmap = LinearSegmentedColormap.from_list('my_custom_colormap', plt.cm.plasma(np.linspace(0.1, 0.85, 256)))

# suppress SyntaxWarning
import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)

# create defualt colors
C1, C2, C3, C4 = _cmap(32/255), _cmap(140/255), _cmap(200/255), _cmap(254/255)

# print rgb triplet for each color as integers (0-255)
print(f"C1: {[int(C1[i]*255) for i in range(3)]}, C2: {[int(C2[i]*255) for i in range(3)]}, C3: {[int(C3[i]*255) for i in range(3)]}, C4: {[int(C4[i]*255) for i in range(3)]}") 

# plot 4 big sqares with these colors
if SHOW_PLOTS:
    plt.figure(figsize=(THESIS_W, 2))
    plt.subplot(241), plt.imshow([[C1]], aspect='auto'), plt.axis('off'), plt.title('C1')
    plt.subplot(242), plt.imshow([[C2]], aspect='auto'), plt.axis('off'), plt.title('C2')
    plt.subplot(243), plt.imshow([[C3]], aspect='auto'), plt.axis('off'), plt.title('C3')
    plt.subplot(244), plt.imshow([[C4]], aspect='auto'), plt.axis('off'), plt.title('C4')
    # full-width colormap (manual axes instead of broken span)
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    plt.axes([0.125, 0.1, 0.775, 0.35]), plt.imshow(gradient, aspect='auto', cmap=_cmap), plt.yticks([])
    plt.show(block=WAIT_FOR_PLOTS); plt.close()


# list all colormaps
# print(plt.colormaps())

# %% ## config
## config
# DATA_DIR_TEST = 'data/3452622/' # best with MSE loss (small)
# DATA_DIR_TEST = 'data/3452764/' # best with L1 loss (small)
# DATA_DIR_TEST = 'data/3456052/' # finetuned over 3452764 (small),
# DATA_DIR_TEST = 'data/3456791/' # best l1 loss turbosmall [change utils.py to use turbosmall config]
DATA_DIR_TEST = 'data/3535303/' # best with L1 loss (small) (3.5 mil ds) <-
# DATA_DIR_TEST = 'data/3884413/' # recent only mag

# %% ## load data
## load data
ds = LiuqeDataset(EVAL_DS_PATH)
lnn = load_net_from_data(data_dir=DATA_DIR_TEST)
N = len(ds.fg)

# %% ## Plot the losses [FX, IY, B, BD, IYD]
## Plot the losses [FX, IY, B, BD, IYD]
train_loss = np.load(f"{DATA_DIR_TEST}train_losses.npy")
eval_loss = np.load(f"{DATA_DIR_TEST}eval_losses.npy")
plt.figure(figsize=(THESIS_W, 2.8))
plt.subplot(1,2,1)
# plt.plot(train_loss[:, 0], label='train $\psi$', color=C1, lw=lw)
plt.plot(train_loss[:, 0], color=C1, label=r'Training $\bm{\psi}$')
plt.plot(eval_loss[:, 0], '--', color=C2, label=r'Validation $\bm{\psi}$')
plt.yscale('log'), plt.xlabel('Epoch'), plt.ylabel(r'Poloidal-flux loss'), plt.legend(), plt.ylim(4e-4, 1e-2)

plt.subplot(1,2,2)
plt.plot(train_loss[:, 3], color=C1, label=r'Training $\bm{B}_d$')
plt.plot(eval_loss[:, 3], '--', color=C2, label=r'Validation $\bm{B}_d$')
plt.yscale('log'), plt.xlabel('Epoch'), plt.ylabel(r'Magnetic-field loss'), plt.legend(), plt.ylim(4e-4, 1e-2)#plt.ylim(1e-3, 1e-2)

save_figure('losses.pdf')

# %% ## Average relative accuracy
## Average relative accuracy
def get_pred_on_grid(ds, lnn):
    fg = ds.fg
    fgp = np.zeros_like(fg)
    n = len(fg)
    p = torch.stack([torch.tensor(RRD.reshape(-1), dtype=torch.float32), torch.tensor(ZZD.reshape(-1), dtype=torch.float32)], dim=1)  # shape (n, 2)
    for i in range(n):
        x = ds[i][0]
        r = p.clone()[:,0]; p.requires_grad = True # need r for PINN
        fxp, _, iyp = lnn(x.unsqueeze(0), p.unsqueeze(0))
        dfx_drz = torch.autograd.grad(outputs=fxp.sum(), inputs=p, create_graph=True)[0].unsqueeze(0) # compute gradients wrt p
        dfx_dr, dfx_dz = dfx_drz[:,:,0], dfx_drz[:,:,1]  # shape (1, n)
        brp = -dfx_dz / (2*π*r) # Br from derivative
        bzp = dfx_dr / (2*π*r)  # Bz from derivative
        d2fx_dr2 = torch.autograd.grad(outputs=dfx_dr.sum(), inputs=p, create_graph=True)[0][:,0].unsqueeze(0)
        d2fx_dz2 = torch.autograd.grad(outputs=dfx_dz.sum(), inputs=p, create_graph=True)[0][:,1].unsqueeze(0)
        iyp = - (1/(2*π*r) * d2fx_dz2 + 1/(2*π*r) * d2fx_dr2 - 1/(2*π*r**2) * dfx_dr)
        fgp[i] = np.array([fxp.detach().numpy(), bzp.detach().numpy(), brp.detach().numpy(), iyp.detach().numpy()]).reshape(fg[i].shape)
    return fgp

fgp = get_pred_on_grid(ds, lnn)

def print_errors(fg, fgp):
    n = len(fg)
    fx, fxp = fg[:,0].reshape(n,-1), fgp[:,0].reshape(n,-1)
    br, brp = fg[:,2].reshape(n,-1), fgp[:,2].reshape(n,-1)
    bz, bzp = fg[:,3].reshape(n,-1), fgp[:,1].reshape(n,-1)
    iy, iyp = fg[:,1].reshape(n,-1), fgp[:,3].reshape(n,-1)
    fx_abs_err = np.abs(fx - fxp)
    bz_abs_err = np.abs(bz - bzp)
    br_abs_err = np.abs(br - brp)
    iy_abs_err = np.abs(iy - iyp)
    print(f'Average abs errors -> Fx: {fx_abs_err.mean():.2e}, Br: {br_abs_err.mean():.2e}, Bz: {bz_abs_err.mean():.2e}, Iy: {iy_abs_err.mean():.2e}')
    fx_spans = (np.max(fx, axis=1) - np.min(fx, axis=1)).reshape(-1, 1)
    bz_spans = (np.max(bz, axis=1) - np.min(bz, axis=1)).reshape(-1, 1)
    br_spans = (np.max(br, axis=1) - np.min(br, axis=1)).reshape(-1, 1)
    iy_spans = (np.max(iy, axis=1) - np.min(iy, axis=1)).reshape(-1, 1)
    fx_perc_err = 100 * fx_abs_err / fx_spans
    bz_perc_err = 100 * bz_abs_err / bz_spans
    br_perc_err = 100 * br_abs_err / br_spans
    iy_perc_err = 100 * iy_abs_err / iy_spans
    print(f'Average % errors   -> Fx: {fx_perc_err.mean():.2f}%, Br: {br_perc_err.mean():.2f}%, Bz: {bz_perc_err.mean():.2f}%, Iy: {iy_perc_err.mean():.2f}%')
    all_abs_err = np.concatenate([fx_abs_err.reshape(n,1,-1), br_abs_err.reshape(n,1,-1), bz_abs_err.reshape(n,1,-1), iy_abs_err.reshape(n,1,-1)], axis=1)
    all_perc_err = np.concatenate([fx_perc_err.reshape(n,1,-1), br_perc_err.reshape(n,1,-1), bz_perc_err.reshape(n,1,-1), iy_perc_err.reshape(n,1,-1)], axis=1)
    return all_abs_err, all_perc_err
abs_err, perc_err = print_errors(ds.fg, fgp)
print(abs_err.shape, perc_err.shape)

# %% ## plot errors by position
## plot errors by position
use_abs_err = True
use_abs_err = False
err = abs_err if use_abs_err else perc_err
title = 'absolute error' if use_abs_err else r'relative error'
error_unit = '' if use_abs_err else r' [\%]'
efx = np.mean(err[:,0].reshape(N, -1), axis=0)
ebr = np.mean(err[:,1].reshape(N, -1), axis=0)
ebz = np.mean(err[:,2].reshape(N, -1), axis=0)
eiy = np.mean(err[:,3].reshape(N, -1), axis=0)
plt.figure(figsize=(THESIS_W, 5.0))
plt.subplot(2, 2, 1)
plt.scatter(RRD.reshape(-1), ZZD.reshape(-1), c=efx, s=20, marker='s', rasterized=True)
prep_plot_tcv(plt.gca(), c='black', alpha=.5), plt.colorbar(label=fr'Mean error{error_unit}')
plt.xlabel(r'$R$ [m]'), plt.ylabel(r'$Z$ [m]'), plt.title(fr'$\psi$: {title}')
plt.subplot(2, 2, 2)
plt.scatter(RRD.reshape(-1), ZZD.reshape(-1), c=ebr, s=20, marker='s', rasterized=True)
prep_plot_tcv(plt.gca(), c='black', alpha=.5), plt.colorbar(label=fr'Mean error{error_unit}')
plt.xlabel(r'$R$ [m]'), plt.ylabel(r'$Z$ [m]'), plt.title(fr'$B_r$: {title}')
plt.subplot(2, 2, 3)
plt.scatter(RRD.reshape(-1), ZZD.reshape(-1), c=ebz, s=20, marker='s', rasterized=True)
prep_plot_tcv(plt.gca(), c='black', alpha=.5), plt.colorbar(label=fr'Mean error{error_unit}')
plt.xlabel(r'$R$ [m]'), plt.ylabel(r'$Z$ [m]'), plt.title(fr'$B_z$: {title}')
plt.subplot(2, 2, 4)
plt.scatter(RRD.reshape(-1), ZZD.reshape(-1), c=eiy, s=20, marker='s', rasterized=True)
prep_plot_tcv(plt.gca(), c='black', alpha=.5), plt.colorbar(label=fr'Mean error{error_unit}')
plt.xlabel(r'$R$ [m]'), plt.ylabel(r'$Z$ [m]'), plt.title(fr'$I_y$: {title}')
save_figure('error_heatmaps.pdf')

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

# %% ## PID shape control simulation 
## PID shape control simulation
CTRL_SIM_SMALL = 'data/pid_shape_net.mat'
CTRL_SIM_STD = 'data/pid_shape_std.mat'
# 'Br', 'Br2', 'Bz', 'Bz2', 'Fx', 'Fx2', 'currs', 'proj_errs', 'rz', 'shape_errs'

ctrl_start = 0.6 # controller start time
vlinec = 'gray'
vlines = '--'

sim_net = loadmat(CTRL_SIM_SMALL)
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

dec = 5
n_ctrl_pts = 10
t = t[::dec]
FxN, BrN, BzN, FxL, BrL, BzL, FxNL, BrNL, BzNL = map(lambda x: x[:n_ctrl_pts,::dec].T, [FxN, BrN, BzN, FxL, BrL, BzL, FxNL, BrNL, BzNL])
IaN, peN, seN, IaL, peL, seL, IaNL, peNL, seNL = map(lambda x: x[:,::dec].T, [IaN, peN, seN, IaL, peL, seL, IaNL, peNL, seNL])

color_cycle = _cmap(np.linspace(0, 1, n_ctrl_pts))
α = 0.7
spl, fs = (1,2), (THESIS_W, 2.8)

plt.figure(figsize=fs)
std_style = '--'

plt.subplot(*spl, 1)
colors = _cmap(np.linspace(0, 1, seN.shape[1]))
for i in range(seN.shape[1]):
    plt.plot(t, seN[:, i], color=colors[i], label='Surrogate' if i==0 else None)
    plt.plot(t, seL[:, i], std_style, color=colors[i], label='LIUQE-RT' if i==0 else None)
plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
plt.xlabel(r'Time $t$ [s]'), plt.ylabel(r'Shape error [m]'), plt.legend()

plt.subplot(*spl, 2)
colors = _cmap(np.linspace(0, 1, peN.shape[1]))
for i in range(peN.shape[1]):
    plt.plot(t, peN[:, i], color=colors[i], label='Surrogate' if i==0 else None)
    plt.plot(t, peL[:, i], std_style, color=colors[i], label='LIUQE-RT' if i==0 else None)
plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
plt.xlabel(r'Time $t$ [s]'), plt.ylabel(r'Projection error [m]'), plt.legend()

# plt.subplot(*spl, 3)
# colors = _cmap(np.linspace(0, 1, IaN.shape[1]))
# for i in range(IaN.shape[1]):
#     plt.plot(t, IaN[:, i], color=colors[i], label='Surrogate' if i==0 else None, lw=lw)
#     plt.plot(t, IaL[:, i], std_style, color=colors[i], label='Liuqe-RT' if i==0 else None, lw=lw)
# plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
# plt.ylabel(r'Currents $[\mathrm{kA}]$'), plt.legend()

save_figure('pid_shape_control_simulation.pdf')

# %% ## MPC shape control simulation
## MPC shape control simulation
CTRL_SIM_SMALL = 'data/mpc_shape_net.mat'
CTRL_SIM_STD = 'data/mpc_shape_std.mat'
# 'shape_errs', 'shape_errs2', 'Fx', 'Fx2', 'Br', 'Bz', 'Br2', 'Bz2', 't', 'currs', 'rz', 'niter', 'mpc_cost', 'duopt'

ctrl_start = 0.6 # controller start time
vlinec = 'gray'
vlines = '--'
n_svd = 6 

sim_net = loadmat(CTRL_SIM_SMALL)
sim_std = loadmat(CTRL_SIM_STD)

FxN, BrN, BzN, IaN, peN, seN = sim_net['Fx'], sim_net['Br'], sim_net['Bz'], sim_net['currs'], sim_net['duopt'][-n_svd:,:], sim_net['shape_errs']
FxL, BrL, BzL, IaL, peL, seL = sim_std['Fx'], sim_std['Br'], sim_std['Bz'], sim_std['currs'], sim_std['duopt'][-n_svd:,:], sim_std['shape_errs']
FxNL, BrNL, BzNL, IaNL, peNL, seNL = sim_net['Fx2'], sim_net['Br2'], sim_net['Bz2'], sim_net['currs'], sim_net['duopt'][-n_svd:,:], sim_net['shape_errs2']
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

dec = 5
n_ctrl_pts = 10
t = t[::dec]
FxN, BrN, BzN, FxL, BrL, BzL, FxNL, BrNL, BzNL = map(lambda x: x[:n_ctrl_pts,::dec].T, [FxN, BrN, BzN, FxL, BrL, BzL, FxNL, BrNL, BzNL])
IaN, peN, seN, IaL, peL, seL, IaNL, peNL, seNL = map(lambda x: x[:,::dec].T, [IaN, peN, seN, IaL, peL, seL, IaNL, peNL, seNL])

color_cycle = _cmap(np.linspace(0, 1, n_ctrl_pts))
α = 0.7
spl, fs = (1,2), (THESIS_W, 2.8)

plt.figure(figsize=fs)
std_style = '--'

plt.subplot(*spl, 1)
colors = _cmap(np.linspace(0, 1, seN.shape[1]))
for i in range(seN.shape[1]):
    plt.plot(t, seN[:, i], color=colors[i], label='Surrogate' if i==0 else None)
    plt.plot(t, seL[:, i], std_style, color=colors[i], label='LIUQE-RT' if i==0 else None)
plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
plt.xlabel(r'Time $t$ [s]'), plt.ylabel(r'Shape error [m]'), plt.legend()

plt.subplot(*spl, 2)
colors = _cmap(np.linspace(0, 1, peN.shape[1]))
for i in range(peN.shape[1]):
    plt.plot(t, peN[:, i], color=colors[i], label='Surrogate' if i==0 else None)
    plt.plot(t, peL[:, i], std_style, color=colors[i], label='LIUQE-RT' if i==0 else None)
plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
plt.xlabel(r'Time $t$ [s]'), plt.ylabel(r'Optimized $\Delta U$ [kA/s]'), plt.legend(title='6 SVD modes')

# plt.subplot(*spl, 3)
# colors = _cmap(np.linspace(0, 1, IaN.shape[1]))
# for i in range(IaN.shape[1]):
#     plt.plot(t, IaN[:, i], color=colors[i], label='Surrogate' if i==0 else None, lw=lw)
#     plt.plot(t, IaL[:, i], std_style, color=colors[i], label='Liuqe-RT' if i==0 else None, lw=lw)
# plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
# plt.ylabel(r'Currents $[\mathrm{kA}]$'), plt.legend()

save_figure('mpc_shape_control_simulation.pdf')

# %% ## FMAG control simulation
## FMAG control simulation
CTRL_SIM_SMALL = 'data/fmag_net.mat'
CTRL_SIM_STD = 'data/fmag_std.mat'
# expected keys: 't', 'rY', 'zY', 'Ip', 'Va', 'mixed_errors', 'observables', 'xpt_flux', 'xpt_flux_lih'
ctrl_start = 0.6  # controller start time
vlinec = 'gray'
vlines = '--'

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

# %% ## Plot network outputs for the dissertation
## Plot network outputs for the dissertation
def plot_network_outputs_dissertation(ds:LiuqeDataset, model:liuqeML, idxs=[0], ylab=[''], save_name='net_example.pdf'):
    assert ds.fg is not None, "Dataset does not have full grid data (fg)"
    assert len(idxs) == len(ylab), "Length of idxs and ylab must be the same"
    model.to(CPU)
    n_plt = len(idxs)

    def prep_plot_tcv(ax, cols=['black','black','gray'], lw=1.5, alpha=0.5):
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal', adjustable='box')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.plot(VESS[:,0], VESS[:,1], color=cols[2], lw=lw, alpha=alpha) # most inner
        ax.plot(VESSI[:,0], VESSI[:,1], color=cols[1], lw=lw, alpha=alpha) # inner
        ax.plot(VESSO[:,0], VESSO[:,1], color=cols[0], lw=lw, alpha=alpha) # outer
        fill_between_polygons(ax, VESSI, VESSO, color='gray', alpha=alpha*0.8, lw=0)
        fill_between_polygons(ax, VESS, VESSI, color='gray', alpha=alpha*0.3, lw=0)
        return ax

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    plt.figure(figsize=(THESIS_W, 4*n_plt))
    gs = plt.GridSpec(n_plt, 1)
    for pi, i in enumerate(idxs):
        x, p, fx, iy, br, bz = ds[i]
        x, p, fx, iy, br, bz = map(lambda x: x.to(CPU), [x, p, fx, iy, br, bz])
        μ0 = 4e-7 * np.pi  # Permeability of free space
        p = torch.stack([torch.tensor(RRD.reshape(-1), dtype=torch.float32), torch.tensor(ZZD.reshape(-1), dtype=torch.float32)], dim=1)  # shape (n, 2)
        fx, iy, br, bz = [to_tensor(v.reshape(-1)) for v in ds.fg[i]] # reshape to (n,)
        iy = iy*μ0/GRA # workaround, on grid points dont have the conversion applied
        r = p.clone()[:,0]; p.requires_grad = True # need r for PINN
        reverse_cmap = iy.mean() < 0  # reverse cmap if current density is negative
        fxp, _, _ = model(x.unsqueeze(0), p.unsqueeze(0))
        dfx_drz = torch.autograd.grad(outputs=fxp.sum(), inputs=p, create_graph=True)[0].unsqueeze(0) # compute gradients wrt p
        π = torch.pi
        dfx_dr, dfx_dz = dfx_drz[:,:,0], dfx_drz[:,:,1]  # shape (1, n)
        brp = -dfx_dz / (2*π*r) # Br from derivative
        bzp = dfx_dr / (2*π*r)  # Bz from derivative
        d2fx_dr2 = torch.autograd.grad(outputs=dfx_dr.sum(), inputs=p, create_graph=True)[0][:,0].unsqueeze(0)
        d2fx_dz2 = torch.autograd.grad(outputs=dfx_dz.sum(), inputs=p, create_graph=True)[0][:,1].unsqueeze(0)
        iyp = - (1/(2*π*r) * d2fx_dz2 + 1/(2*π*r) * d2fx_dr2 - 1/(2*π*r**2) * dfx_dr)

        fxp, iyp, brp, bzp = map(lambda t: t.squeeze(0), [fxp, iyp, brp, bzp]) # remove batch dimension
        p, fx, iy, br, bz, fxp, iyp, brp, bzp = map(lambda t: t.detach().cpu().numpy(), [p, fx, iy, br, bz, fxp, iyp, brp, bzp]) # move to CPU and numpy
        iv = inpoly(p) # inside vessel indices
        # mfx, Mfx = np.min([fx[iv], fxp[iv]]), np.max([fx[iv], fxp[iv]]) # both maps
        # miy, Miy = np.min([iy[iv], iyp[iv]]), np.max([iy[iv], iyp[iv]])
        # mbr, Mbr = np.min([br[iv], brp[iv]]), np.max(lw=1.5, alpha=0.3[br[iv], brp[iv]])
        # mbz, Mbz = np.min([bz[iv], bzp[iv]]), np.max([bz[iv], bzp[iv]])
        mfx, Mfx = np.min(fx[iv]), np.max(fx[iv]) # only actual
        miy, Miy = np.min(iy[iv]), np.max(iy[iv])
        mbr, Mbr = np.min(br[iv]), np.max(br[iv])
        mbz, Mbz = np.min(bz[iv]), np.max(bz[iv])
        # mfx, Mfx = np.min(fxp[iv]), np.max(fxp[iv]) # only predicted
        # miy, Miy = np.min(iyp[iv]), np.max(iyp[iv])
        # mbr, Mbr = np.min(brp[iv]), np.max(brp[iv])
        # mbz, Mbz = np.min(bzp[iv]), np.max(bzp[iv])
        ms, cm = 4, plt.rcParams['image.cmap'] + ('_r' if reverse_cmap else '') # line width, color, marker size, colormap
        limits = [(mfx, Mfx), (miy, Miy), (mbr, Mbr), (mbz, Mbz)] # variable limits
        vals = [(fx, fxp), (iy, iyp), (br, brp), (bz, bzp)]
        perc_max = [1.0, 5, 2.0, 2.0] # max percentage error for colorbar scaling, set higher for Fx and Iy which have larger ranges
        names = [r'$\psi(r,z)$', r'$I_y(r,z)$', r'$B_r(r,z)$', r'$B_z(r,z)$'] # with (r,z)
        # names = [r'$\psi$', r'$B_r$', r'$B_z$', r'$I_y$'] # without (r,z)
        # lvls = [15, 12, 12, 8] # <-
        lvls = [20, 12, 16, 16]
        lw, shrink = 0.75, 0.9
        group = gs[pi].subgridspec(2, 4)
        plt.subplot(gs[pi], frameon=False); plt.axis('off'); plt.title(ylab[pi].capitalize())
        for j, ((v, vp), l, n) in enumerate(zip(vals, limits, names)):
            row, col = j//2, 2*(j % 2)
            plt.subplot(group[row, col])
            plt.contour(RRD, ZZD, vp.reshape(65,28), levels=lvls[j], linestyles='solid', cmap=cm, linewidths=lw)
            plt.contour(RRD, ZZD, v.reshape(65,28), levels=lvls[j], linestyles='dashed', cmap=cm, linewidths=lw)
            prep_plot_tcv(plt.gca(), cols=['k','k','gray'])
            sm = ScalarMappable(cmap=cm, norm=Normalize(vmin=l[0], vmax=l[1]))
            sm.set_array([])
            plt.colorbar(sm, ax=plt.gca(), shrink=shrink, fraction=0.08, pad=0.03)
            plt.ylabel(n)
            plt.subplot(group[row, col+1])
            ep = 100*np.abs(v[iv]-vp[iv])/(l[1]-l[0]) # error in %
            plt.scatter(p[iv,0], p[iv,1], c=ep, s=ms, vmin=0, vmax=perc_max[j], marker='s', rasterized=True) # fixed bar but not for Iy
            prep_plot_tcv(plt.gca(), cols=['k','k','gray']); plt.colorbar(shrink=shrink, fraction=0.08, pad=0.03); plt.ylabel(n + r' $\%$ err')

    save_figure(save_name)
    return


# interesting_idxs = [297, 250, 317, 28, 347, 12, 65]
# interesting_idxs = [297, 250, 123, 28, 347, 12, 157]
# interesting_idxs = [297, 250, 123, 28, 347, 157]
# interesting_idxs, labels = [297, 123, 28, 347, 157], ["single null", "double null", "limiter", 'negative\ntriangularity', 'high vertical\nposition']
# interesting_idxs, labels = [297, 123, 347, 157], ['single null', 'double null', 'negative triangularity', 'high vertical position']
# interesting_idxs, labels = [297, 123, 347], ['single null', 'double null', 'negative triangularity']

interesting_idxs1, labels1 = [297, 123], ["single null", "double null", ]
interesting_idxs2, labels2 = [28, 347], ["limiter", 'negative triangularity']
interesting_idxs3, labels3 = [157], ['high vertical position']

plot_network_outputs_dissertation(ds, lnn, idxs=interesting_idxs1, ylab=labels1, save_name='net_example1.pdf')
plot_network_outputs_dissertation(ds, lnn, idxs=interesting_idxs2, ylab=labels2, save_name='net_example2.pdf')
plot_network_outputs_dissertation(ds, lnn, idxs=interesting_idxs3, ylab=labels3, save_name='net_example3.pdf')


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

if SHOW_PLOTS: plt.show()
