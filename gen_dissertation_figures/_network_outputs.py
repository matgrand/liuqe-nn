from ._common import *


def predict_network_output_fields(x, model):
    """Evaluate flux and its field/current derivatives on the full LIUQE grid."""
    model.to(CPU)
    x = torch.as_tensor(x, dtype=torch.float32, device=CPU)
    p = torch.stack([torch.tensor(RRD.reshape(-1), dtype=torch.float32), torch.tensor(ZZD.reshape(-1), dtype=torch.float32)], dim=1)
    r = p.clone()[:,0]; p.requires_grad = True
    fx, _, _ = model(x.unsqueeze(0), p.unsqueeze(0))
    dfx_drz = torch.autograd.grad(outputs=fx.sum(), inputs=p, create_graph=True)[0].unsqueeze(0)
    dfx_dr, dfx_dz = dfx_drz[:,:,0], dfx_drz[:,:,1]
    br = -dfx_dz/(2*torch.pi*r); bz = dfx_dr/(2*torch.pi*r)
    d2fx_dr2 = torch.autograd.grad(outputs=dfx_dr.sum(), inputs=p, create_graph=True)[0][:,0].unsqueeze(0)
    d2fx_dz2 = torch.autograd.grad(outputs=dfx_dz.sum(), inputs=p, create_graph=True)[0][:,1].unsqueeze(0)
    iy = -(d2fx_dz2/(2*torch.pi*r) + d2fx_dr2/(2*torch.pi*r) - dfx_dr/(2*torch.pi*r**2))
    return tuple(value.squeeze(0).detach().cpu().numpy().reshape(RRD.shape) for value in (fx, iy, br, bz))


# %% ## Plot network outputs for the dissertation
## Plot network outputs for the dissertation
def plot_network_output_samples_dissertation(samples, model:liuqeML, ylab=None, save_name='net_example.pdf'):
    """Plot preprocessed ``(x, full_grid)`` samples as reference, prediction, and percentage-error panels."""
    ylab = ['']*len(samples) if ylab is None else ylab
    assert len(samples) == len(ylab), "Length of samples and ylab must be the same"
    model.to(CPU)
    n_plt = len(samples)

    def prep_plot_tcv(ax, cols=['black','black','gray'], lw=plt.rcParams['lines.linewidth'], alpha=0.5):
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
    for pi, (x, fg) in enumerate(samples):
        μ0 = 4e-7 * np.pi  # Permeability of free space
        fx, iy, br, bz = [np.asarray(v).reshape(-1) for v in fg]
        iy = iy*μ0/GRA # workaround, on grid points dont have the conversion applied
        reverse_cmap = iy.mean() < 0  # reverse cmap if current density is negative
        fxp, iyp, brp, bzp = [value.reshape(-1) for value in predict_network_output_fields(x, model)]
        p = np.stack((RRD, ZZD), axis=-1).reshape(-1, 2)
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
        plt.subplot(gs[pi], frameon=False); plt.axis('off'); plt.title(ylab[pi])
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


def plot_network_outputs_dissertation(ds:LiuqeDataset, model:liuqeML, idxs=[0], ylab=[''], save_name='net_example.pdf'):
    """Compatibility wrapper for plotting samples from a ``LiuqeDataset``."""
    assert ds.fg is not None, "Dataset does not have full grid data (fg)"
    samples = [(ds[i][0], ds.fg[i]) for i in idxs]
    return plot_network_output_samples_dissertation(samples, model, ylab=ylab, save_name=save_name)


# interesting_idxs = [297, 250, 317, 28, 347, 12, 65]
# interesting_idxs = [297, 250, 123, 28, 347, 12, 157]
# interesting_idxs = [297, 250, 123, 28, 347, 157]
# interesting_idxs, labels = [297, 123, 28, 347, 157], ["single null", "double null", "limiter", 'negative\ntriangularity', 'high vertical\nposition']
# interesting_idxs, labels = [297, 123, 347, 157], ['single null', 'double null', 'negative triangularity', 'high vertical position']
# interesting_idxs, labels = [297, 123, 347], ['single null', 'double null', 'negative triangularity']
