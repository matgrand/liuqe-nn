from ._common import *
from ._common import _cmap
from ._network_outputs import plot_network_outputs_dissertation


def generate():
    interesting_idxs3, labels3 = [157, 250], ['high vertical position', 'tbd']
    plot_network_outputs_dissertation(ds, lnn, idxs=interesting_idxs3, ylab=labels3, save_name='net_example3.pdf')
