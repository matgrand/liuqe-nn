# Figure description: This two-panel figure illustrates the trainable-slope Swish activation used by the coordinate-conditioned neural operator. The left panel compares the activation profiles $\operatorname{Swish}_\beta(x)=x\sigma(\beta x)$ for several slope parameters, showing their transition from a smooth, nearly linear response to a sharper ReLU-like response. The right panel shows the corresponding analytical derivatives used for coordinate differentiation, with each curve retaining the same color as its activation profile.
from ._common import *
from ._common import _cmap


N_BETAS = 7
BETA_RANGE = (0.2, 2.0)
X_LIM = (-5.0, 5.0)


def _swish(x, beta):
    sigmoid = 1 / (1 + np.exp(-beta*x))
    return x*sigmoid, sigmoid + beta*x*sigmoid*(1-sigmoid)


def generate(n_betas=N_BETAS, beta_range=BETA_RANGE):
    assert 5 <= n_betas <= 10, 'Use between 5 and 10 slope values for a readable figure'
    x = np.linspace(*X_LIM, 500); betas = np.geomspace(*beta_range, n_betas)
    colors = _cmap(np.linspace(0.08, 0.92, n_betas))

    plt.figure(figsize=(THESIS_W, 2.8))
    plt.subplot(1, 2, 1)
    for beta, color in zip(betas, colors):
        swish, _ = _swish(x, beta)
        plt.plot(x, swish, color=color, label=rf'$\beta={beta:.3f}$')
    plt.xlabel(r'Input $x$'); plt.ylabel(r'$\operatorname{Swish}_\beta(x)$'); plt.title(r'(a) Swish activation'); plt.legend(ncol=2)

    plt.subplot(1, 2, 2)
    for beta, color in zip(betas, colors):
        _, derivative = _swish(x, beta)
        plt.plot(x, derivative, color=color, label=rf'$\beta={beta:.3f}$')
    plt.xlabel(r'Input $x$'); plt.ylabel(r'$\mathrm{d}\operatorname{Swish}_\beta/\mathrm{d}x$'); plt.title(r'(b) Analytical derivative'); plt.ylim(-0.15, 1.15); #plt.legend(ncol=2)

    save_figure('swish_act.pdf')
