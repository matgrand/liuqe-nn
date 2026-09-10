"""Interactively inspect raw LIUQE flux and current-density maps.

The browser loads one random MATLAB shot at a time and displays ten equispaced valid time slices of the stored poloidal flux and plasma current density. It intentionally performs no interpolation, magnetic-field reconstruction, dataset construction, or neural-network loading, providing a minimal and responsive view of the raw equilibria.
"""
from pathlib import Path

from ._common import *
from ._raw_data import raw_valid_mask

RAW_DIR = Path('raw_ds')
N_TIMES = 7


class RawDatasetBrowser:
    """Keyboard-driven raw-shot browser using ``n`` for next and ``q`` for quit."""
    def __init__(self, raw_dir=RAW_DIR, n_times=N_TIMES, seed=None):
        self.raw_dir, self.n_times = Path(raw_dir), n_times
        self.files = sorted(self.raw_dir.glob('*.mat')); self.rng = np.random.default_rng(seed)
        if not self.files: raise FileNotFoundError(f'No .mat files found in {self.raw_dir}.')
        self.current_path = None
        self.figure = plt.figure(figsize=(18, 4)); self.figure.canvas.mpl_connect('key_press_event', self.on_key)
        try: self.figure.canvas.manager.set_window_title('Raw LIUQE browser — n: next, q: quit')
        except AttributeError: pass

    def load_random_shot(self):
        """Load a random shot with enough valid samples and select equispaced time indices."""
        for path in self.rng.permutation(self.files):
            if len(self.files) > 1 and path == self.current_path: continue
            data = loadmat(path, variable_names=['t', IP, FX, IY]); valid = np.flatnonzero(raw_valid_mask(data))
            if len(valid) < self.n_times: continue
            positions = np.rint(np.linspace(0, len(valid)-1, self.n_times)).astype(int)
            self.current_path = path
            return path, data, valid[positions]
        raise RuntimeError(f'No raw shot contains {self.n_times} valid time slices.')

    def draw_random_shot(self):
        r"""Replace the current canvas with raw $\psi$ and $I_y$ maps from a new shot."""
        path, data, indices = self.load_random_shot(); self.figure.clear()
        nt = np.asarray(data['t']).size
        times = np.asarray(data['t']).reshape(-1)
        fx = np.asarray(data[FX]).reshape(*DS_SIZES[FX], nt)
        iy = np.asarray(data[IY]).reshape(*DS_SIZES[IY], nt)
        for col, index in enumerate(indices):
            iy_full = np.zeros(RRD.shape, dtype=iy.dtype); iy_full[1:-1,1:-1] = iy[...,index]
            rq, zq = find_approx_lcfs(iy_full)
            plt.subplot(2, self.n_times, col+1)
            plt.pcolormesh(RRD, ZZD, fx[...,index], shading='auto', rasterized=True)
            plt.scatter(rq, zq, s=3, c='cyan', linewidths=0, rasterized=True)
            prep_plot_tcv(plt.gca(), lw=0.7, alpha=0.45, c='gray')
            if col == 0: plt.ylabel(r'$\psi$')
            plt.title(rf'{path.stem}:{index}' + '\n' + rf'$t={times[index]:.3f}\,$s')
            plt.subplot(2, self.n_times, self.n_times+col+1)
            plt.pcolormesh(RRD2, ZZD2, iy[...,index], shading='auto', rasterized=True)
            plt.scatter(rq, zq, s=3, c='cyan', linewidths=0, rasterized=True)
            prep_plot_tcv(plt.gca(), lw=0.7, alpha=0.45, c='gray')
            if col == 0: plt.ylabel(r'$I_y$')
        self.figure.suptitle(f'Shot {path.stem} — raw maps, cyan: approximate LCFS — n: next, q: quit')
        plt.tight_layout(rect=(0, 0, 1, 0.95)); self.figure.canvas.draw_idle()
        print(f'Showing shot {path.stem}: n = next random shot, q = quit')

    def on_key(self, event):
        """Handle the browser controls."""
        if event.key and event.key.lower() == 'n': self.draw_random_shot()
        elif event.key and event.key.lower() == 'q': plt.close(self.figure)

    def show(self):
        """Draw the first shot and start the GUI event loop."""
        self.draw_random_shot(); plt.show(block=True)


def main(): RawDatasetBrowser().show()


if __name__ == '__main__': main()
