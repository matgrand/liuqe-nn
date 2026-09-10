"""Interactive browser for finding interesting plasma configurations in raw shots.

Each random shot is summarized by equispaced valid time slices. Full mode compares reference and network-predicted poloidal flux and radial field and displays their pointwise relative-error maps; reference-only mode skips model loading and inference for rapid configuration browsing.
"""
from pathlib import Path

from ._common import *
from ._network_outputs import predict_network_output_fields
from ._raw_data import RAW_NAMES, preprocess_raw_sample, raw_valid_mask

RAW_DIR = Path('raw_ds')
N_TIMES = 8
FAST = True
SKIP_INFERENCE = False
FAST_FIELD_K_FINER = 1
FAST_FIELD_INTERP_METHOD = 'linear'
ERROR_MAX = {FX:1.0, BR:2.0}


def relative_error(actual, predicted):
    """Return pointwise error as a percentage of the in-vessel reference range."""
    span = np.ptp(actual[IV])
    return 100*np.abs(actual-predicted)/span if span > 0 else np.zeros_like(actual)


class ShapeBrowser:
    """Keyboard-driven random-shot browser using ``n`` for next and ``q`` for quit."""
    def __init__(self, raw_dir=RAW_DIR, model=None, n_times=N_TIMES, fast=FAST, skip_inference=SKIP_INFERENCE, seed=None):
        self.raw_dir, self.n_times, self.fast, self.skip_inference = Path(raw_dir), n_times, fast, skip_inference
        self.model = None if skip_inference else (load_dissertation_net() if model is None else model)
        self.files = sorted(self.raw_dir.glob('*.mat')); self.rng = np.random.default_rng(seed)
        if not self.files: raise FileNotFoundError(f'No .mat files found in {self.raw_dir}.')
        self.current_path = None
        self.figure = plt.figure(figsize=(18, 4 if skip_inference else 8)); self.figure.canvas.mpl_connect('key_press_event', self.on_key)
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
            field_k_finer = FAST_FIELD_K_FINER if self.fast else 7
            field_interp_method = FAST_FIELD_INTERP_METHOD if self.fast else INTERP_METHOD
            return path, [preprocess_raw_sample(data, path.stem, int(index), require_valid=False, field_k_finer=field_k_finer, field_interp_method=field_interp_method) for index in valid[positions]]
        raise RuntimeError(f'No raw shot contains {self.n_times} valid time slices.')

    def draw_random_shot(self):
        """Replace the current canvas with a new random-shot overview."""
        path, samples = self.load_random_shot(); self.figure.clear()
        nrows = 2 if self.skip_inference else 4
        for col, sample in enumerate(samples):
            fx, iy, br, _ = sample.fg
            if self.skip_inference: predictions = (None, None)
            else:
                fxp, _, brp, _ = predict_network_output_fields(sample.x, self.model); predictions = (fxp, brp)
            rq, zq = find_approx_lcfs(iy); lcfs_flux = np.mean(interp_pts(fx, np.stack((rq, zq), axis=1), method='linear'))
            for row, (actual, predicted, name) in enumerate(zip((fx, br), predictions, (FX, BR))):
                reference_row = row if self.skip_inference else 2*row
                plt.subplot(nrows, self.n_times, reference_row*self.n_times+col+1)
                if self.fast: plt.pcolormesh(RRD, ZZD, np.where(IV, actual, np.nan), shading='auto', rasterized=True)
                else:
                    if predicted is not None: plt.contour(RRD, ZZD, predicted, levels=12, linestyles='solid', linewidths=0.6)
                    plt.contour(RRD, ZZD, actual, levels=12, linestyles='solid' if predicted is None else 'dashed', linewidths=0.6)
                plt.contour(RRD, ZZD, fx, levels=[lcfs_flux], colors='cyan', linewidths=1.0)
                prep_plot_tcv(plt.gca(), lw=0.7, alpha=0.45, c='gray')
                if col == 0: plt.ylabel(r'$\psi$' if name == FX else r'$B_r$')
                if row == 0: plt.title(rf'{sample.shot}:{sample.time_index}' + '\n' + rf'$t={sample.time:.3f}\,$s')
                if not self.skip_inference:
                    plt.subplot(nrows, self.n_times, (2*row+1)*self.n_times+col+1)
                    error = np.where(IV, relative_error(actual, predicted), np.nan)
                    plt.pcolormesh(RRD, ZZD, error, shading='auto', vmin=0, vmax=ERROR_MAX[name], rasterized=True)
                    plt.contour(RRD, ZZD, fx, levels=[lcfs_flux], colors='cyan', linewidths=1.0)
                    prep_plot_tcv(plt.gca(), lw=0.7, alpha=0.45, c='gray')
                    if col == 0: plt.ylabel((r'$\psi$' if name == FX else r'$B_r$') + r' err. [\%]')
        if self.skip_inference: plot_description = 'reference maps'
        elif self.fast: plot_description = 'reference and relative-error maps'
        else: plot_description = 'solid: network, dashed: reference'
        self.figure.suptitle(f'Shot {path.stem} — {plot_description}, cyan: LCFS — n: next, q: quit')
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
