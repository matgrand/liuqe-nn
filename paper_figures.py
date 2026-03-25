# %%
## Script for the figures in the paper
from utils import *
FIGS_DIR = 'paper_figures'
os.makedirs(FIGS_DIR, exist_ok=True)
plt.style.use('default')
plt.rcParams['figure.figsize'] = (8,3)
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.7
plt.rcParams['font.family'] = 'monospace'
plt.rcParams['image.cmap'] = 'inferno'
# plt.rcParams['image.cmap'] = 'plasma'
# plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.colormaps['plasma'](np.array([0.2,0.4,0.6,0.8,1.0])))
# plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.colormaps['plasma'](np.array([.1,.5,.2,.6,.3,.7,.4,.8,.9,1.0]))) # <-
plt.rcParams['text.usetex'] = True
plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
plt.rcParams["text.latex.preamble"] = r"\usepackage{bm}"
plt.rcParams['font.size'] = 12

from matplotlib.colors import LinearSegmentedColormap
_cmap = LinearSegmentedColormap.from_list('my_custom_colormap', plt.cm.inferno(np.linspace(0.1, 0.85, 256)))
# _cmap = LinearSegmentedColormap.from_list('my_custom_colormap', plt.cm.plasma(np.linspace(0.1, 0.85, 256)))

# create defualt colors
C1, C2, C3, C4 = _cmap(0.0), _cmap(0.5), _cmap(0.7), _cmap(0.999)

# plot 4 big sqares with these colors
plt.figure(figsize=(8,2))
plt.subplot(1,4,1), plt.imshow([[C1]], aspect='auto'), plt.axis('off'), plt.title('C1')
plt.subplot(1,4,2), plt.imshow([[C2]], aspect='auto'), plt.axis('off'), plt.title('C2')
plt.subplot(1,4,3), plt.imshow([[C3]], aspect='auto'), plt.axis('off'), plt.title('C3')
plt.subplot(1,4,4), plt.imshow([[C4]], aspect='auto'), plt.axis('off'), plt.title('C4')
plt.tight_layout()
# plt.savefig(f"{FIGS_DIR}/colors.pdf")
plt.show()

# suppress SyntaxWarning
import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)

# list all colormaps
# print(plt.colormaps())

# %%
## config
# DATA_DIR_TEST = 'data/3452622/' # best with MSE loss (small)
# DATA_DIR_TEST = 'data/3452764/' # best with L1 loss (small)
# DATA_DIR_TEST = 'data/3456052/' # finetuned over 3452764 (small),
# DATA_DIR_TEST = 'data/3456791/' # best l1 loss turbosmall [change utils.py to use turbosmall config]
DATA_DIR_TEST = 'data/3535303/' # best with L1 loss (small) (3.5 mil ds) <-
# DATA_DIR_TEST = 'data/3884413/' # recent only mag

# %%
## load data
ds = LiuqeDataset(EVAL_DS_PATH)
lnn = load_net_from_data(data_dir=DATA_DIR_TEST)
N = len(ds.fg)

# %%
## Plot the losses [FX, IY, B, BD, IYD]
train_loss = np.load(f"{DATA_DIR_TEST}train_losses.npy")
eval_loss = np.load(f"{DATA_DIR_TEST}eval_losses.npy")

plt.figure(figsize=(8, 4))
plt.subplot(2,1,1)
# plt.plot(train_loss[:, 0], label='train $\psi$', color='tab:blue')
plt.plot(train_loss[:, 0], label=r'train $\bm\psi$', color='tab:blue')
plt.plot(eval_loss[:, 0], label='eval $\psi$', color='tab:red')
plt.yscale('log'), plt.xlabel('Epoch'), plt.ylabel('$\psi$ Loss'), plt.legend(), plt.ylim(4e-4, 1e-2)

plt.subplot(2,1,2)
plt.plot(train_loss[:, 3], label='train $\mathbf{B_d}$', color='tab:blue')
plt.plot(eval_loss[:, 3], label='eval $\mathbf{B_d}$', color='tab:red')
plt.yscale('log'), plt.xlabel('Epoch'), plt.ylabel('$\mathbf{B_d}$ Loss'), plt.legend(), plt.ylim(1e-3, 1e-2)

# plt.subplot(1, 2, 3)
# plt.plot(train_loss[:, 4], label='train_Iyd')
# plt.plot(eval_loss[:, 4], label='eval_Iyd')
# plt.yscale('log'), plt.xlabel('Epoch'), plt.ylabel('Loss'), plt.legend()

plt.tight_layout()
plt.savefig(f"{FIGS_DIR}/losses.pdf")
plt.show()

# %%
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

# %%
## plot errors by position
use_abs_err = True
use_abs_err = False
err = abs_err if use_abs_err else perc_err
title = 'Abs Error' if use_abs_err else '$\%$ error'
efx = np.mean(err[:,0].reshape(N, -1), axis=0)
ebr = np.mean(err[:,1].reshape(N, -1), axis=0)
ebz = np.mean(err[:,2].reshape(N, -1), axis=0)
eiy = np.mean(err[:,3].reshape(N, -1), axis=0)
plt.figure(figsize=(8,3))
plt.subplot(1, 4, 1)
plt.scatter(RRD.reshape(-1), ZZD.reshape(-1), c=efx, s=20, marker='s')
prep_plot_tcv(plt.gca(), c='white', alpha=.5), plt.colorbar()
plt.xlabel('R'), plt.ylabel('Z'), plt.title(fr'$\psi$ {title}')
plt.subplot(1, 4, 2)
plt.scatter(RRD.reshape(-1), ZZD.reshape(-1), c=ebr, s=20, marker='s')
prep_plot_tcv(plt.gca(), c='white', alpha=.5), plt.colorbar()
plt.xlabel('R'), plt.ylabel('Z'), plt.title(fr'$B_r$ {title}')
plt.subplot(1, 4, 3)
plt.scatter(RRD.reshape(-1), ZZD.reshape(-1), c=ebz, s=20, marker='s')
prep_plot_tcv(plt.gca(), c='white', alpha=.5), plt.colorbar()
plt.xlabel('R'), plt.ylabel('Z'), plt.title(fr'$B_z$ {title}')
plt.subplot(1, 4, 4)
plt.scatter(RRD.reshape(-1), ZZD.reshape(-1), c=eiy, s=20, marker='s')
prep_plot_tcv(plt.gca(), c='white', alpha=.5), plt.colorbar()
plt.xlabel('R'), plt.ylabel('Z'), plt.title(fr'$I_y$ {title}')
plt.tight_layout()
plt.savefig(f"{FIGS_DIR}/error_heatmaps.pdf")
plt.show()

# %%
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
spl, fs = (3, 2), (8, 8)

plt.figure(figsize=fs)
std_style = '--'
plt.subplot(*spl, 1) # top: Fx Br Bz. bot: se, pe, Ia
for cp in reversed(range(n_ctrl_pts)):
    zero_line = 0.0025 * cp * np.ones_like(t)
    er = eFx[:, cp] + zero_line
    plt.fill_between(t, zero_line, er, alpha=α, color=color_cycle[cp])
plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
plt.ylabel(f'Error $\psi$')
plt.subplot(*spl, 2)
for cp in reversed(range(n_ctrl_pts)):
    zero_line = 0.0025 * cp * np.ones_like(t)
    er = eBr[:, cp] + zero_line
    plt.fill_between(t, zero_line, er, alpha=α, color=color_cycle[cp])
plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
plt.ylabel(f'Error $B_r$')
plt.subplot(*spl, 3)
for cp in reversed(range(n_ctrl_pts)):
    zero_line = 0.0025 * cp * np.ones_like(t)
    er = eBz[:, cp] + zero_line
    plt.fill_between(t, zero_line, er, alpha=α, color=color_cycle[cp])
plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
plt.ylabel(f'Error $B_z$')

plt.subplot(*spl, 4)
colors = _cmap(np.linspace(0, 1, seN.shape[1]))
for i in range(seN.shape[1]):
    plt.plot(t, seN[:, i], color=colors[i])
    plt.plot(t, seL[:, i], std_style, color=colors[i])
plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
plt.ylabel(f'Shape Err')

plt.subplot(*spl, 5)
colors = _cmap(np.linspace(0, 1, peN.shape[1]))
for i in range(peN.shape[1]):
    plt.plot(t, peN[:, i], color=colors[i])
    plt.plot(t, peL[:, i], std_style, color=colors[i])
plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
plt.ylabel(f'Proj Err')

plt.subplot(*spl, 6)
colors = _cmap(np.linspace(0, 1, IaN.shape[1]))
for i in range(IaN.shape[1]):
    plt.plot(t, IaN[:, i], color=colors[i])
    plt.plot(t, IaL[:, i], std_style, color=colors[i])
plt.axvline(x=ctrl_start, color=vlinec, linestyle=vlines)
plt.ylabel(f'Currents')

plt.tight_layout()
plt.savefig(f"{FIGS_DIR}/control_simulation.pdf")
plt.show()

# %%
## plot network outputs for the paper
def plot_network_outputs_paper(ds:LiuqeDataset, model:liuqeML, idxs=[0], ylab=['']):
    assert ds.fg is not None, "Dataset does not have full grid data (fg)"
    assert len(idxs) == len(ylab), "Length of idxs and ylab must be the same"
    model.to(CPU)
    n_plt = len(idxs)

    def prep_plot_tcv(ax, cols=['black','black','gray'], lw=1.5, alpha=0.5):
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
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

    plt.figure(figsize=(14, 3*n_plt))
    for pi, i in enumerate(idxs):
        x, p, fx, iy, br, bz = ds[i]
        x, p, fx, iy, br, bz = map(lambda x: x.to(CPU), [x, p, fx, iy, br, bz])
        μ0 = 4e-7 * np.pi  # Permeability of free space
        p = torch.stack([torch.tensor(RRD.reshape(-1), dtype=torch.float32), torch.tensor(ZZD.reshape(-1), dtype=torch.float32)], dim=1)  # shape (n, 2)
        fx, iy, br, bz = [to_tensor(v.reshape(-1)) for v in ds.fg[i]] # reshape to (n,)
        iy = iy*μ0/GRA # workaround, on grid points dont have the conversion applied
        r = p.clone()[:,0]; p.requires_grad = True # need r for PINN
        reverse_cmap = iy.mean() < 0  # reverse cmap if current density is negative
        fxp, brzp, iyp = model(x.unsqueeze(0), p.unsqueeze(0))
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
        limits = [(mfx, Mfx), (mbr, Mbr), (mbz, Mbz), (miy, Miy)] # variable limits
        vals = [(fx, fxp), (br, brp), (bz, bzp), (iy, iyp)]
        perc_max = [1.0, 2.0, 2.0, 5] # max percentage error for colorbar scaling, set higher for Fx and Iy which have larger ranges
        names = ['$\psi(r,z)$', '$B_r(r,z)$', '$B_z(r,z)$', '$I_y(r,z)$']
        lvls = [15, 12, 12, 8]
        lw = 1
        for j, ((v, vp), l, n) in enumerate(zip(vals, limits, names)):
            plt.subplot(n_plt, 8, pi*8+2*j+1),
            if pi==0: plt.title(f"{n}")
            plt.contour(RRD, ZZD, vp.reshape(65,28), levels=lvls[j], linestyles='solid', cmap=cm, linewidths=lw)
            plt.contour(RRD, ZZD, v.reshape(65,28), levels=lvls[j], linestyles='dashed', cmap=cm, linewidths=lw)
            prep_plot_tcv(plt.gca(), cols=['k','k','gray']); cbar = plt.colorbar()
            cbar.ax.tick_params(labelsize=8)  # smaller tick labels

            if j==0: plt.ylabel(ylab[pi])

            plt.subplot(n_plt, 8, pi*8+2*j+2)
            if pi==0: plt.title(f"{n} $\%$ err")
            ep = 100*np.abs(v[iv]-vp[iv])/(l[1]-l[0]) # error in %
            plt.scatter(p[iv,0], p[iv,1], c=ep, s=ms, vmin=0, vmax=perc_max[j], marker='s') # fixed bar but not for Iy
            prep_plot_tcv(plt.gca(), cols=['k','k','gray']); cbar = plt.colorbar()
            cbar.ax.tick_params(labelsize=8)  # smaller tick labels

    plt.tight_layout()
    plt.savefig(f"{FIGS_DIR}/net_example.pdf")
    if LOCAL: plt.show()
    plt.close()
    return


idxs = list(np.random.choice(len(ds), size=10, replace=False))
# interesting_idxs = [297, 250, 317, 28, 347, 12, 65]
# interesting_idxs = [297, 250, 123, 28, 347, 12, 157]
# interesting_idxs = [297, 250, 123, 28, 347, 157]
interesting_idxs = [297, 123, 28, 347, 157]
labels = ["single null", "double null", "limiter", "negative triangularity", "high vertical position"]
plot_network_outputs_paper(ds, lnn, idxs=interesting_idxs, ylab=labels)
# plot_network_outputs_paper(ds, lnn, idxs=idxs, ylab=labels)

# %%
## Comparison with LIUQE LIUQE-RT and LIH
realtime_file_path = 'data/realtime_analysis_shot88903.mat'
realtime_file_path = 'data/realtime_analysis_shot88566.mat'

# 't_or', 'phys_or', 'Fx_or', 'Br_or', 'Bz_or', 'Fx_net', 'Br_net', 'Bz_net', 'Fx_lrt', 'Br_lrt', 'Bz_lrt', 'Fx_lih', 'Br_lih', 'Bz_lih', 'Fx_liu', 'Br_liu', 'Bz_liu', 'Fx_rec', 'Br_rec', 'Bz_rec', 'Fx_rec2', 'Br_rec2', 'Bz_rec2', 't', 'phys');

m = loadmat(realtime_file_path)
net = (m['Fx_net'], m['Br_net'], m['Bz_net'])
lrt = (m['Fx_lrt'], m['Br_lrt'], m['Bz_lrt'])
lih = (m['Fx_lih'], m['Br_lih'], m['Bz_lih'])
liu = (m['Fx_liu'], m['Br_liu'], m['Bz_liu'])
# rec = (m['Fx_rec'], m['Br_rec'], m['Bz_rec'])
# rec2 = (m['Fx_rec2'], m['Br_rec2'], m['Bz_rec2'])

t, phys = m['t'].reshape(-1), m['phys']

# t_or = m['t_or'].reshape(-1)
# phys_or = m['phys_or']
# orig = (m['Fx_or'], m['Br_or'], m['Bz_or'])

# print(f'net: {net[0].shape}|{net[1].shape}|{net[2].shape}, lrt: {lrt[0].shape}|{lrt[1].shape}|{lrt[2].shape}, lih: {lih[0].shape}|{lih[1].shape}|{lih[2].shape}, liu: {liu[0].shape}|{liu[1].shape}|{liu[2].shape}')

enet = list(map(lambda a, b: np.abs(a-b), net, liu))
elrt = list(map(lambda a, b: np.abs(a-b), lrt, liu))
elih = list(map(lambda a, b: np.abs(a-b), lih, liu))

plt.figure(figsize=(10, 5))
cols = ['tab:red', 'tab:blue', 'tab:orange']
titles = ['$\psi$', '$B_r$', '$B_z$']
for k in range(3):
    plt.subplot(3, 1, k+1)
    for i, (errs, name) in enumerate(zip([enet, elrt, elih], ['Surrogate', 'Liuqe-RT', 'LIH'])):
        min_err, max_err, avg_err = np.min(errs[k], axis=1), np.max(errs[k], axis=1), np.mean(errs[k], axis=1)
        plt.plot(t, avg_err, label=name, color=cols[i], linewidth=1.5)
        plt.fill_between(t, avg_err, max_err, color=cols[i], alpha=0.3)
    plt.ylim(2e-4, 1e-1)
    plt.yscale('log')
    plt.ylabel(titles[k])
    plt.legend(loc='upper right')
# plt.suptitle('Error wrt Offline LIUQE')
plt.tight_layout()
plt.savefig(f"{FIGS_DIR}/comparison1.pdf")
if LOCAL: plt.show()
