from ._common import *
from ._common import _cmap
from ._network_outputs import plot_network_outputs_dissertation


def generate():
    interesting_idxs1, labels1 = [297, 123], ["single null", "double null"]
    plot_network_outputs_dissertation(ds, lnn, idxs=interesting_idxs1, ylab=labels1, save_name='net_example1.pdf')

    interesting_idxs2, labels2 = [28, 347], ["limiter", 'negative triangularity']
    plot_network_outputs_dissertation(ds, lnn, idxs=interesting_idxs2, ylab=labels2, save_name='net_example2.pdf')

    interesting_idxs3, labels3 = [157, 250], ['high vertical position', 'tbd']
    plot_network_outputs_dissertation(ds, lnn, idxs=interesting_idxs3, ylab=labels3, save_name='net_example3.pdf')