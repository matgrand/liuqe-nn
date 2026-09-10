"""Network-output examples loaded directly from raw LIUQE shot files.

The figure reproduces the six evaluation examples and panel styling of the original dataset-backed implementation while bypassing ``LiuqeDataset``. Each raw time slice is filtered and preprocessed as in ``prepare_dataset.ipynb``: diagnostic vectors are concatenated in training order, current density is padded onto the full grid, and magnetic-field targets are reconstructed through fine-grid flux differentiation. The row titles identify the source shot, making each example traceable to the original measurements.
"""
from pathlib import Path

from ._common import load_dissertation_net
from ._network_outputs import plot_network_output_samples_dissertation
from ._raw_data import load_raw_samples

RAW_DIR = Path('raw_ds')
SELECTION_GROUPS = [
    [(86168, 120), (80265, 120)],
    [(83825, 325), (80265, 313)],
    # [(87037, 116), (85803, 130)],
    # [(82223, 212), (85803, 114)],
]


def generate(raw_dir=RAW_DIR, model=None):
    """Generate the three raw-shot dissertation network-example figures."""
    if model is None: model = load_dissertation_net()
    for i, selections in enumerate(SELECTION_GROUPS, start=1):
        samples = load_raw_samples(raw_dir, selections)
        plot_samples = [(sample.x, sample.fg) for sample in samples]
        plot_network_output_samples_dissertation(plot_samples, model, ylab=[sample.shot for sample in samples], save_name=f'net_example{i}.pdf')


if __name__ == '__main__': generate()
