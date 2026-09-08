from ._common import *
from ._common import _cmap
from ._network_outputs import plot_network_outputs_dissertation


def generate():
    interesting_idxs2, labels2 = [28, 347], ["limiter", 'negative triangularity']
    plot_network_outputs_dissertation(ds, lnn, idxs=interesting_idxs2, ylab=labels2, save_name='net_example2.pdf')
