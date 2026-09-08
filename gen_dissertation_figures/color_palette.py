from ._common import *
from ._common import _cmap
from . import _common


def generate():
    # plot 4 big sqares with these colors
    if _common.SHOW_PLOTS:
        plt.figure(figsize=(THESIS_W, 2))
        plt.subplot(241), plt.imshow([[C1]], aspect='auto'), plt.axis('off'), plt.title('C1')
        plt.subplot(242), plt.imshow([[C2]], aspect='auto'), plt.axis('off'), plt.title('C2')
        plt.subplot(243), plt.imshow([[C3]], aspect='auto'), plt.axis('off'), plt.title('C3')
        plt.subplot(244), plt.imshow([[C4]], aspect='auto'), plt.axis('off'), plt.title('C4')
        # full-width colormap (manual axes instead of broken span)
        gradient = np.linspace(0, 1, 256).reshape(1, -1)
        plt.axes([0.125, 0.1, 0.775, 0.35]), plt.imshow(gradient, aspect='auto', cmap=_cmap), plt.yticks([])
        plt.tight_layout()
        plt.show(block=_common.WAIT_FOR_PLOTS); plt.close()
