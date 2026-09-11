# %% ## Script for the dissertation figures
## Script for the dissertation figures
from utils import *
from pathlib import Path
import shutil
import subprocess

FIGS_DIR = 'dissertation_figures'
THESIS_W = 6.30  # [in]
# SHOW_PLOTS = True
SHOW_PLOTS = False
WAIT_FOR_PLOTS = False
SORT_INSPECTION_FIGURES_ALPHABETICALLY = True
os.makedirs(FIGS_DIR, exist_ok=True)

# plt.style.use('dissertation-thesis')
plt.style.use('dissertation-mono')

def save_figure(filename, rect=None):
    path = os.path.join(FIGS_DIR, filename)
    plt.tight_layout(rect=rect); plt.savefig(path); print(f'Saved {path}')
    if SHOW_PLOTS: plt.show(block=WAIT_FOR_PLOTS)
    else: plt.close()


def plot_approx_lcfs(flux, current, axis=None, flux_offset=0.005, **kwargs):
    r"""Estimate the LCFS flux, offset it by a fraction of the map span, and plot its $\psi$ contour."""
    if flux.shape != RRD.shape or current.shape != RRD.shape: raise ValueError(f'Flux and current must both have shape {RRD.shape}, got {flux.shape} and {current.shape}.')
    rq, zq = find_approx_lcfs(current)
    if not rq.size: raise ValueError('Cannot estimate the LCFS from an empty current support.')
    boundary = np.stack((rq, zq), axis=-1)
    boundary_flux = np.mean(interp_pts(flux, boundary, method='linear')) + flux_offset*np.ptp(flux)
    axis = plt.gca() if axis is None else axis
    defaults = {'colors':'white', 'linestyles':':', 'linewidths':1.0, 'zorder':4}
    # return axis.contour(RRD, ZZD, flux, levels=[boundary_flux], **(defaults | kwargs))
    contour_levels = [boundary_flux - flux_offset*np.ptp(flux), boundary_flux, boundary_flux + flux_offset*np.ptp(flux)]
    return axis.contour(RRD, ZZD, flux, levels=contour_levels, **(defaults | kwargs))


def compile_inspection_pdf():
    """Compile a PDF that places every generated dissertation figure on its own page."""
    script_dir = Path(__file__).resolve().parent
    figures_dir = script_dir.parent / FIGS_DIR
    build_dir = script_dir / 'build'
    build_dir.mkdir(exist_ok=True)
    sort_key = (lambda pdf: pdf.name.lower()) if SORT_INSPECTION_FIGURES_ALPHABETICALLY else (lambda pdf: pdf.stat().st_mtime)
    figure_pdfs = sorted(figures_dir.glob('*.pdf'), key=sort_key)
    figure_entries = '\n'.join(
        rf'\inspectfigure{{{pdf.name}}}{{\texttt{{\detokenize{{{pdf.stem}}}}}}}{{fig:inspection-{index:02d}}}'
        for index, pdf in enumerate(figure_pdfs, start=1)
    )
    figs_tex = build_dir / 'figs.tex'
    figs_tex.write_text(
        rf'''\documentclass[a4paper]{{article}}
\usepackage[margin=20mm]{{geometry}}
\usepackage{{graphicx}}

\newcommand{{\inspectfigure}}[3]{{%
  \begin{{figure}}[p]
    \centering
    \includegraphics[width=\textwidth]{{../dissertation_figures/#1}}
    \caption{{#2}}
    \label{{#3}}
  \end{{figure}}
  \clearpage
}}

\begin{{document}}
{figure_entries}
\end{{document}}
''',
        encoding='utf-8',
    )
    if shutil.which('latexmk') is None:
        print('Figure inspection PDF: FAIL')
        return
    result = subprocess.run(
        ['latexmk', '-pdf', '-file-line-error', '-f', '-outdir=build', 'build/figs.tex'],
        cwd=script_dir,
        capture_output=True,
        text=True,
    )
    print(f'Figure inspection PDF: {"FAIL" if result.returncode else "SUCCESS"}')


from matplotlib.colors import LinearSegmentedColormap
# _cmap = LinearSegmentedColormap.from_list('my_custom_colormap', plt.cm.inferno(np.linspace(0.1, 0.85, 256)))
_base_cmap = plt.get_cmap(plt.rcParams['image.cmap'])
_cmap = LinearSegmentedColormap.from_list('my_custom_colormap', _base_cmap(np.linspace(0.05, 0.9, 256)))
# _cmap = LinearSegmentedColormap.from_list('my_custom_colormap', plt.cm.plasma(np.linspace(0.1, 0.85, 256)))

# suppress SyntaxWarning
import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)

# create defualt colors
C1, C2, C3, C4 = _cmap(32/255), _cmap(140/255), _cmap(200/255), _cmap(254/255)

# print rgb triplet for each color as integers (0-255)
print(f"C1: {[int(C1[i]*255) for i in range(3)]}, C2: {[int(C2[i]*255) for i in range(3)]}, C3: {[int(C3[i]*255) for i in range(3)]}, C4: {[int(C4[i]*255) for i in range(3)]}") 

# list all colormaps
# print(plt.colormaps())

# %% ## config
## config
# DATA_DIR_TEST = 'data/3452622/' # best with MSE loss (small)
# DATA_DIR_TEST = 'data/3452764/' # best with L1 loss (small)
# DATA_DIR_TEST = 'data/3456052/' # finetuned over 3452764 (small),
# DATA_DIR_TEST = 'data/3456791/' # best l1 loss turbosmall [change utils.py to use turbosmall config]
DATA_DIR_TEST = 'data/3535303/' # best with L1 loss (small) (3.5 mil ds) <-
# DATA_DIR_TEST = 'data/3884413/' # recent only mag

# %% ## load model
## load model
def load_dissertation_net(): return load_net_from_data(data_dir=DATA_DIR_TEST)
