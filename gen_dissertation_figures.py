# %% Generate dissertation figures
from pathlib import Path
import shutil
import subprocess

import gen_dissertation_figures as figures

# figures.configure(SHOW_PLOTS=True) # uncomment to show plots
# figures.configure(WAIT_FOR_PLOTS=True) # uncomment to not wait for plots to be closed before continuing

# Comment out completed figures while working on the remaining ones.
# figures.color_palette.generate()
# figures.losses.generate()
# figures.error_heatmaps.generate()
# figures.control_simulation.generate()
# figures.pid_shape_control_simulation.generate()
# figures.mpc_shape_control_simulation.generate()
# figures.fmag_control_simulation.generate()
# figures.tcv_geometry.generate()
# figures.net_example1.generate()
# figures.net_example2.generate()
# figures.net_example3.generate()
# figures.comparison.generate()

figures.pid_shape_control_simulation.generate()
figures.mpc_shape_control_simulation.generate()
figures.fmag_control_simulation.generate()

# figures.tcv_geometry.generate()

script_dir = Path(__file__).resolve().parent / 'gen_dissertation_figures'
figures_dir = script_dir.parent / 'dissertation_figures'
build_dir = script_dir / 'build'
build_dir.mkdir(exist_ok=True)
figure_pdfs = sorted(figures_dir.glob('*.pdf'))
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
    print('Skipping figure inspection PDF compilation: latexmk is not installed.')
else:
    result = subprocess.run(
        ['latexmk', '-pdf', '-file-line-error', '-f', '-outdir=build', 'build/figs.tex'],
        cwd=script_dir,
    )
    if result.returncode:
        print(f'Figure inspection PDF compilation failed (exit code {result.returncode}).')
    else:
        print(f'Saved {script_dir / "build" / "figs.pdf"}')

import matplotlib.pyplot as plt
plt.show()
