# Figure description: This four-panel figure connects the prescribed optimization schedules to the observed convergence of LiuqeNN. The upper row shows the exponentially decaying Adam learning rate and the loss weights used to balance the poloidal-flux, current-density, and physics-informed magnetic-field objectives. The lower row compares training and validation losses for the directly predicted poloidal flux and the magnetic field recovered from its spatial derivatives, demonstrating stable convergence without a widening generalization gap.
from ._common import *
from ._common import _cmap


LR_START, LR_END = 5e-3, 5e-5
LAMBDA_PSI, LAMBDA_B = 1.0, 0.5
LAMBDA_J_START, LAMBDA_J_END = 1e-2, 1e-6


def generate():
    ## Plot the losses [FX, IY, B, BD, IYD]
    train_loss = np.load(f"{DATA_DIR_TEST}train_losses.npy")
    eval_loss = np.load(f"{DATA_DIR_TEST}eval_losses.npy")
    epochs = np.arange(1, len(train_loss)+1)
    lr = np.geomspace(LR_START, LR_END, len(epochs))
    lambda_psi = LAMBDA_PSI*np.ones_like(epochs); lambda_b = LAMBDA_B*np.ones_like(epochs)
    lambda_j = np.geomspace(LAMBDA_J_START, LAMBDA_J_END, len(epochs))

    plt.figure(figsize=(THESIS_W, 4))
    plt.subplot(2, 2, 1)
    plt.plot(epochs, lr, color=C1)
    plt.yscale('log'); plt.xlabel('Epoch'); plt.ylabel(r'Learning rate $\eta$'); plt.title(r'(a) Learning-rate schedule')

    plt.subplot(2, 2, 2)
    plt.plot(epochs, lambda_psi, color=C1, label=r'$\lambda_{\psi}$')
    plt.plot(epochs, lambda_b, color=C2, label=r'$\lambda_{B}$')
    plt.plot(epochs, lambda_j, color=C3, label=r'$\lambda_{j}$')
    plt.yscale('log'); plt.xlabel('Epoch'); plt.ylabel(r'Loss weight $\lambda$'); plt.title(r'(b) Loss-weight schedules'); plt.legend()

    plt.subplot(2, 2, 3)
    # plt.plot(train_loss[:, 0], label='train $\psi$', color=C1, lw=lw)
    plt.plot(epochs, train_loss[:, 0], color=C1, label='Training')
    plt.plot(epochs, eval_loss[:, 0], '--', color=C2, label='Validation')
    plt.yscale('log'); plt.xlabel('Epoch'); plt.ylabel(r'Poloidal-flux loss'); plt.title(r'(c) Poloidal flux $\psi$'); plt.legend(); plt.ylim(4e-4, 1e-2)

    plt.subplot(2, 2, 4)
    plt.plot(epochs, train_loss[:, 3], color=C1, label='Training')
    plt.plot(epochs, eval_loss[:, 3], '--', color=C2, label='Validation')
    plt.yscale('log'); plt.xlabel('Epoch'); plt.ylabel(r'Magnetic-field loss'); plt.title(r'(d) Derived field $\bm{B}$'); plt.legend(); plt.ylim(4e-4, 1e-2)#plt.ylim(1e-3, 1e-2)

    save_figure('sched_losses.pdf')

if __name__ == '__main__': generate()
