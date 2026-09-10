"""Interactive browser for finding interesting plasma configurations in raw shots.

Each random shot is summarized by ten equispaced valid time slices. The window compares reference and network-predicted poloidal-flux and radial-field contours and displays their pointwise relative-error maps, allowing useful configurations to be identified rapidly before they are promoted into dissertation figures.
"""
from pathlib import Path

from ._common import *
from ._network_outputs import predict_network_output_fields
from ._raw_data import RAW_NAMES, preprocess_raw_sample, raw_valid_mask

RAW_DIR = Path('raw_ds')
N_TIMES = 10
ERROR_MAX = {FX:1.0, BR:2.0}


def relative_error(actual, predicted):
    """Return pointwise error as a percentage of the in-vessel reference range."""
    span = np.ptp(actual[IV])
    return 100*np.abs(actual-predicted)/span if span > 0 else np.zeros_like(actual)


class ShapeBrowser:
    """Keyboard-driven random-shot browser using ``n`` for next and ``q`` for quit."""
    def __init__(self, raw_dir=RAW_DIR, model=lnn, n_times=N_TIMES, seed=None):
        self.raw_dir, self.model, self.n_times = Path(raw_dir), model, n_times
        self.files = sorted(self.raw_dir.glob('*.mat')); self.rng = np.random.default_rng(seed)
        if not self.files: raise FileNotFoundError(f'No .mat files found in {self.raw_dir}.')
        self.current_path = None
        self.figure = plt.figure(figsize=(18, 8)); self.figure.canvas.mpl_connect('key_press_event', self.on_key)
        try: self.figure.canvas.manager.set_window_title('Raw-shot configuration browser — n: next, q: quit')
        except AttributeError: pass

    def load_random_shot(self):
        """Load a random shot containing enough valid slices and preprocess equispaced selections."""
        for path in self.rng.permutation(self.files):
            if len(self.files) > 1 and path == self.current_path: continue
            data = loadmat(path, variable_names=['t', *RAW_NAMES]); valid = np.flatnonzero(raw_valid_mask(data))
            if len(valid) < self.n_times: continue
            positions = np.rint(np.linspace(0, len(valid)-1, self.n_times)).astype(int)
            self.current_path = path
            return path, [preprocess_raw_sample(data, path.stem, int(index), require_valid=False) for index in valid[positions]]
        raise RuntimeError(f'No raw shot contains {self.n_times} valid time slices.')

    def draw_random_shot(self):
        """Replace the current canvas with a new random-shot overview."""
        path, samples = self.load_random_shot(); self.figure.clear()
        for col, sample in enumerate(samples):
            fx, _, br, _ = sample.fg; fxp, _, brp, _ = predict_network_output_fields(sample.x, self.model)
            for row, (actual, predicted, name) in enumerate(((fx, fxp, FX), (br, brp, BR))):
                plt.subplot(4, self.n_times, row*2*self.n_times+col+1)
                plt.contour(RRD, ZZD, predicted, levels=12, linestyles='solid', linewidths=0.6)
                plt.contour(RRD, ZZD, actual, levels=12, linestyles='dashed', linewidths=0.6)
                prep_plot_tcv(plt.gca(), lw=0.7, alpha=0.45, c='gray')
                if col == 0: plt.ylabel(r'$\psi$' if name == FX else r'$B_r$')
                if row == 0: plt.title(rf'{sample.shot}:{sample.time_index}' + '\n' + rf'$t={sample.time:.3f}\,$s')
                plt.subplot(4, self.n_times, (row*2+1)*self.n_times+col+1)
                error = np.where(IV, relative_error(actual, predicted), np.nan)
                plt.pcolormesh(RRD, ZZD, error, shading='auto', vmin=0, vmax=ERROR_MAX[name], rasterized=True)
                prep_plot_tcv(plt.gca(), lw=0.7, alpha=0.45, c='gray')
                if col == 0: plt.ylabel((r'$\psi$' if name == FX else r'$B_r$') + r' err. [\%]')
        self.figure.suptitle(f'Shot {path.stem} — solid: network, dashed: reference — n: next, q: quit')
        plt.tight_layout(rect=(0, 0, 1, 0.96)); self.figure.canvas.draw_idle()
        print(f'Showing shot {path.stem}: n = next random shot, q = quit')

    def on_key(self, event):
        """Handle the two browser controls."""
        if event.key and event.key.lower() == 'n': self.draw_random_shot()
        elif event.key and event.key.lower() == 'q': plt.close(self.figure)

    def show(self):
        """Draw the first shot and start the GUI event loop."""
        self.draw_random_shot(); plt.show(block=True)


def main(): ShapeBrowser().show()


if __name__ == '__main__': main()
