from ._common import *
from ._common import _cmap


def generate():
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
