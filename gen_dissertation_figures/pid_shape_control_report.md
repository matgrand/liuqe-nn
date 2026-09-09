# PID shape-control comparison

## What the controller is and what this simulation tests

The generalized TCV shape controller is an outer loop around the core hybrid magnetic controller. The target boundary is represented by control points. Flux and poloidal-field values evaluated at those points are transformed by `STP.SC2FB` into the saved shape-feedback vector; subtracting this vector from `SC.ref` gives `shapeErrs`. The PID path then projects the errors with `STP.Tpre`, filters the projected errors, applies diagonal PID laws, maps the result back through `STP.Tpost`, and adds the resulting corrections to hybrid references. The SVD-based pre/post maps retain the directions that the closed-loop plant can control most effectively. This follows the controller architecture and the motivation for truncating weak directions described by Anand et al., while the current implementation includes later changes in the `hybrid` repository.

The data-generation script is `/home/mg/repos/hybrid/pid_shape_liuqenn_comparison.m`. It uses TCV shot 79742, a 100 microsecond FGE and hybrid time step, a 2 ms shape-controller period, controller activation at 0.6 s, a fixed shape reference, no allocator, and a linearized FGE plant (`Lfge.P.lin = true`). The simulation ends at 0.8999 s even though the controller tunables specify a later switch-off. The surrogate run calls `liuqenn_inference` from synthetic magnetic measurements. The comparison run does not execute LIUQE-RT: `shapectrl.m` evaluates the simulated FGE flux map directly with `qintmex`. It is therefore a comparison with an FGE-map interpolation reference, not a software-in-the-loop execution of the complete LIUQE-RT reconstruction algorithm. The previous plot and dissertation text label this baseline as LIUQE-RT, which should be corrected.

## What is saved and why

- `shape_errs` and `shape_errs2` store the active-observer and direct-FGE shape errors. They are the regulated outputs and the main closed-loop comparison.
- `proj_errs` and `proj_errs2` store `STP.Tpre * shapeErrs`. These are the errors seen by the diagonal PID channels and show how physical boundary errors are combined into controllable directions.
- `Fx`, `Br`, and `Bz` store the active observer at the 25 allocated query slots; the `*2` variants store direct interpolation of the FGE map in the same run. These permit observer comparisons on an identical plant trajectory. In the surrogate files, unused query slots must not be interpreted: the code replaces nonpositive radial coordinates by 1 mm before inference, and several unused outputs have nonphysical values.
- `currs` combines selected hybrid-reference corrections and shape-current outputs, while `rz` stores two selected hybrid signals. Their channel names, transformation matrices, and units are not saved, so the arrays are unsuitable for a thesis plot with physical labels.
- `t` is the controller debug time base. The file contains 25 allocated shape-error rows, of which 11 are nonzero in this case, and 16 projected rows, of which four are nonzero. The reported dissertation metric uses the first ten shape-error rows. Over $t\geq0.6$ s, their pooled RMS values are 2.75 mm for the surrogate and 3.60 mm for direct FGE interpolation.

## Level 1: improvement implemented from the existing files

`pid_shape_control_simulation.py` now produces a 6.30 inch by 5.60 inch four-panel figure. It keeps the first ten error trajectories so the visual result and scalar RMS use exactly the same channels, converts meters to millimeters, shows the four active projected modes, adds the instantaneous ten-channel RMS, and adds per-channel post-activation RMS values. Solid and dashed lines retain the observer distinction, and the baseline is now labelled `FGE/qintmex`. The figure no longer plots unlabeled current arrays. The generated PDF is `dissertation_figures/pid_shape_control_simulation.pdf`.

The choice to retain ten channels is conservative: it reproduces the existing table exactly. The eleventh nonzero saved shape-error row is not added because the data file does not contain the `SC.dims` and point-type masks needed to identify it, and silently changing the metric would break comparison with the stated 2.75 mm and 3.60 mm results.

## Level 2: changes recommended for the MATLAB collection script

Save the simulation and controller metadata with each result: the exact baseline name, `Lfge.P.lin`, shot, time interval, `FP.hybTs`, `FP.shapeTs`, controller activation, and observer model identifier. Save `SC.dims`, all point-type masks, control-point coordinates, references, `STP.Tpre`, `STP.Tpost`, retained singular values, and PID gains so every raw and projected channel can be named. Save the final voltage command, active-coil currents with `L.G.dima`, hybrid-reference corrections with `cpar.refs.dim`, FGE convergence flags, the actual boundary coordinates, and geometric point-to-boundary errors. These additions would support a physical actuator panel, a boundary cross-section at selected times, and an explicit account of the eleventh active error channel.

For validation beyond this deterministic case, collect matched runs for nonlinear FGE, several equilibria and topologies, reference steps, measurement noise, profile perturbations, and actuator saturation. The present linear fixed-reference case is useful for detecting observer-induced changes, but it cannot establish robustness over the controller operating domain.

## Level 3: recommended dissertation narrative

Introduce the signal chain before the figure: flux and field at target points are converted into shape errors, projected into four active controller directions, passed through the PID laws, and mapped into corrections to the inner hybrid loop. State that FGE supplies the plant and synthetic diagnostics and that the baseline is direct interpolation of the FGE flux map. Describe the test as a controlled comparison against the simulation state, not as execution of LIUQE-RT.

Use the expanded figure in this order: raw ten-channel errors, projected controller errors, their instantaneous aggregate, and the per-channel post-activation summary. Define the pooled RMS next to the figure and state why only ten channels enter it. The supported result is that both closed loops settle to millimeter-scale residuals in this one linear scenario and that the surrogate does not degrade the reported aggregate. Do not infer superiority from the lower surrogate RMS without repeated scenarios and disturbances.

## Sources read

- Dissertation section `Controller-level validation in FGE` in `/home/mg/repos/phd_thesis/2_eqrec.tex`.
- H. Anand et al., *A Novel Plasma Position and Shape Controller for Advanced Configuration Development on the TCV Tokamak*, local file `/home/mg/repos/phd_thesis/material/tcv_ctrl_infra/P27_Anand2017.pdf`.
- `/home/mg/repos/phd_thesis/guides/guides.md`, `hybrid_guide.md`, and `fge_guide.md`.
- `/home/mg/repos/hybrid/pid_shape_liuqenn_comparison.m`, `ctrls/shape/shapectrl.m`, `ctrls/shape/shapectrl_step.m`, and the associated controller-building code.

