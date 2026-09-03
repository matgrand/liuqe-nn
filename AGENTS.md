# AGENTS.md

# General instructions
- for python writing, follow instructions at /home/mg/llm_instructions/python_code_style.md




# DISSERTATION PLOT INSTRUCTIONS
## Dissertation Plotting
These rules apply to figures intended for the dissertation. They override the generic plotting rules above when there is a conflict.

### Style

Use the dissertation Matplotlib style for all dissertation figures:

```python
plt.style.use('dissertation-mono') #or
plt.style.use('dissertation-thesis')
```

The default dissertation style is `dissertation-mono`:

* LaTeX text rendering;
* monospace text;
* standard LaTeX math;
* 10 pt text throughout;
* grid enabled by the style;
* `inferno` default colormap;
* to have matching colors across the figures, do the same as in paper_figures.py, see this for example:
    ```python
    from matplotlib.colors import LinearSegmentedColormap
    _cmap = LinearSegmentedColormap.from_list('my_custom_colormap', plt.cm.inferno(np.linspace(0.05, 0.9, 256)))
    # create defualt colors
    C1, C2, C3, C4 = _cmap(32/255), _cmap(140/255), _cmap(200/255), _cmap(254/255)
    ```
* PDF as the preferred output format;
* 300 DPI for rasterized content embedded in figures.

Do not manually duplicate or override style settings inside individual plotting scripts unless required by a specific figure.

If exact serif matching with the dissertation is explicitly requested, use:

```python
plt.style.use('dissertation-thesis')
```

This uses TeX Gyre Pagella, matching the dissertation's pdfLaTeX fallback. Do not introduce PGF/XeLaTeX/EB Garamond plotting workflows unless explicitly requested.

### Physical figure size

The dissertation text width is:

```python
THESIS_W = 6.30  # [in]
```

The fundamental rule is:

> Generate figures at approximately their final physical size in the dissertation. Do not generate oversized figures and shrink them later in LaTeX.

The default full-width figure is:

```python
plt.figure()  # 6.30 x 3.50 in from dissertation style
```

The width should normally remain `6.30` in. Change the height as required by the content:

```python
plt.figure(figsize=(6.30, 4.2))
```

For a genuinely half-width figure:

```python
plt.figure(figsize=(3.09, 2.4))
```

Prefer full-width figures unless the plot is simple enough to remain readable at half width.

Do not reduce font sizes to make a crowded figure fit. Instead:

* increase the figure height;
* simplify the plot;
* change the subplot layout;
* move to a full-width figure.

### Text size

Use **10 pt consistently throughout dissertation figures**.

Do not set per-figure `fontsize` values unless explicitly requested. In particular, keep approximately the same physical size for:

* axis labels;
* tick labels;
* legends;
* titles, when used;
* annotations.

Changing figure height does not require changing text size.

Avoid substantially resizing exported figures in LaTeX, since resizing also changes the effective text, line, marker, and spacing sizes.

### LaTeX notation

Use LaTeX consistently for mathematical quantities and symbols:

```python
plt.xlabel(r'Time $t$ [s]')
plt.ylabel(r'Plasma current $I_p$ [MA]')
plt.legend(title=r'$I_p = 250\,\mathrm{kA}$')
```

Conventions:

* use dissertation notation exactly when known;
* write physical quantities as mathematical symbols, e.g. `$I_p$`, `$q_{95}$`, `$\beta_p$`;
* put units in square brackets;
* keep notation consistent across related figures;
* use raw Python strings for LaTeX-heavy labels.

Do not redefine the LaTeX preamble inside individual scripts. The dissertation style already loads the required `amsmath` and `bm` packages.

### Plot construction

Continue following the global repository plotting conventions:

* use `matplotlib.pyplot` directly;
* use `plt.figure()` to start figures;
* use `plt.subplot()` for multi-panel figures;
* prefer `plt.*` over explicit `Axes` objects;
* do not add explicit colors unless requested;
* do not configure grids inside plotting code;
* guard interactive display with `SHOW_PLOTS`;
* keep plotting code dense and readable.

Example:

```python
plt.figure()
plt.plot(t, ip, label=r'$I_p$')
plt.plot(t, ip_ref, '--', label=r'$I_p^{\mathrm{ref}}$')
plt.xlabel(r'Time $t$ [s]'); plt.ylabel(r'Plasma current $I_p$ [MA]'); plt.legend()
plt.tight_layout(); plt.savefig('data/plasma_current.pdf'); print('Saved data/plasma_current.pdf')
if SHOW_PLOTS: plt.show()
```

For related multi-panel plots, prefer one physical Matplotlib figure:

```python
plt.figure()
plt.subplot(1, 2, 1); ...
plt.subplot(1, 2, 2); ...
plt.tight_layout()
```

Prefer this over generating independent oversized panels and scaling them individually in LaTeX.

### Titles, legends, and annotations

For dissertation figures:

* avoid figure titles when the dissertation caption already explains the figure;
* use legends only when they improve interpretation;
* avoid redundant labels and annotations;
* keep annotations concise;
* keep the same semantic meaning for colors, line styles, and markers across related figures whenever possible.

Scientific readability takes priority over filling unused space.

### Saving figures

Prefer vector PDF:

```python
plt.tight_layout(); plt.savefig('data/figure_name.pdf'); print('Saved data/figure_name.pdf')
if SHOW_PLOTS: plt.show()
```

Always call:

```python
plt.tight_layout()
```

before saving or showing.

Save **before** `plt.show()`.

Do not use:

```python
bbox_inches='tight'
```

and do not globally set:

```python
plt.rcParams['savefig.bbox'] = 'tight'
```

for dissertation figures.

`bbox_inches='tight'` changes the physical dimensions of the exported canvas and breaks the simple relationship between the Matplotlib figure size and its intended dissertation size.

`plt.tight_layout()` is preferred because it adjusts the content inside the existing physical canvas without changing its nominal size.

Raster-like plots such as `imshow` or `pcolormesh` may still be saved as PDF; rasterized content can remain raster while axes, labels, and other elements remain vector graphics.

### LaTeX inclusion

A normal full-width figure should be generated at approximately `6.30` in and included as:

```latex
\begin{figure}[t]
    \centering
    \includegraphics[width=\textwidth]{my_figure.pdf}
    \caption{...}
    \label{fig:my_figure}
\end{figure}
```

The exported physical width and final LaTeX width should remain approximately equal.

Avoid using LaTeX as a generic resizing mechanism for plots.

If:

```text
W_gen   = generated figure width
W_final = final dissertation width
```

then all visual elements are effectively scaled by:

```text
scale = W_final / W_gen
```

Therefore generating a 12 in figure and inserting it at 6.30 in is not acceptable: fonts, markers, lines, and spacing would all be reduced to approximately 52.5% of their intended size.

### Agent constraints

When editing or generating dissertation plotting code:

* do not change the dissertation `.mplstyle` files unless explicitly instructed;
* do not change the standard 10 pt typography on a figure-by-figure basis;
* do not change the normal 6.30 in full-width convention without a reason;
* figure height may be changed freely to suit the content;
* do not add plotting dependencies beyond Matplotlib unless strictly necessary;
* do not introduce plotting wrapper libraries;
* do not use `bbox_inches='tight'`;
* do not replace the existing `plt.*` coding style with the object-oriented Matplotlib API;
* do not modify unrelated plotting/configuration code while working on a figure.

When a figure is crowded, fix the composition rather than shrinking its typography.

The invariant for dissertation figures is:

> **Create the figure at its final physical width, keep text sizes fixed in points, use consistent LaTeX notation, save as PDF, and avoid substantial post-export scaling.**
