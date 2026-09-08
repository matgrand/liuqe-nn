from ._common import *
from ._common import _cmap
from ._network_outputs import plot_network_outputs_dissertation


def generate():
    interesting_idxs1, labels1 = [297, 123], ["single null", "double null"]
    plot_network_outputs_dissertation(ds, lnn, idxs=interesting_idxs1, ylab=labels1, save_name='net_example1.pdf')
