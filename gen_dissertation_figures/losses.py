from ._common import *
from ._common import _cmap


def generate():
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
