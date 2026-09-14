# PID shape control

## Purpose and outputs

The data collector is `/home/mg/repos/hybrid/pid_shape_liuqenn_comparison.m`, the figure script is `pid_shape_control_simulation.py`, and the generated figure is `dissertation_figures/pid_shape_control_simulation.pdf`. The collector performs matched closed-loop simulations using a neural magnetic observer and direct interpolation of the FGE flux map. It saves controller configuration and timing, point definitions and masks, observer signals, shape errors, projected PID coordinates, actuator signals, and FGE geometry and convergence diagnostics.

The figure contains six panels: all active shape-error trajectories, the active projected controller errors, the first-ten-channel RMS history, the post-activation RMS of each metric channel, the final FGE last closed flux surface (LCFS) and control points, and the norm of the coil-voltage vector. Solid curves denote the neural observer and dashed curves denote direct FGE interpolation with `qintmex`.

## Controller architecture

The generalized TCV shape controller is an outer loop around the fast hybrid magnetic controller. The plasma boundary target is represented by control points. Flux and poloidal-field values at those points are transformed by `STP.SC2FB` into the shape-feedback vector, and subtracting this feedback from `SC.ref` gives `shapeErrs`. `STP.Tpre` projects the physical errors into controllable directions, diagonal PID laws act on the filtered projected errors, and `STP.Tpost` maps the controller outputs into corrections to the hybrid references. The inner hybrid controller converts those corrected references into coil actuation.

The pre/post transformation is constructed from the controlled plant response and retains four active PID directions in this configuration. Eleven of the 25 allocated shape-error rows are active. The metric uses the first ten saved rows, matching the declared dissertation statistic, while the trajectory panel shows every active row identified by `active_shape_mask`.

## Simulation configuration

The simulation uses TCV shot 79742. FGE and hybrid execute with a $100\,\mu\mathrm{s}$ step, the shape controller executes every 2 ms, and shape control activates at 0.6 s. The simulated interval spans the 50 ms preceding activation through 0.8999 s. The FGE plant is linearized (`Lfge.P.lin = true`), the shape reference is fixed, the allocator is disabled, and no reference step is applied.

In the neural run, `shapectrl` evaluates `liuqenn_inference` using synthetic magnetic measurements. In the direct-observer run, `shapectrl` evaluates the simulated FGE flux map at the control points with `qintmex`. This path is an observer based directly on the FGE simulation state; it is not a software-in-the-loop execution of LIUQE-RT. Both runs use the same controller design, initial state, reference, and plant configuration.

## Saved data and interpretation

- `shape_errs` contains the regulated shape errors produced by the active observer. `shape_errs_direct` contains direct-FGE values evaluated within the same run, and the paired MAT files record the closed-loop trajectory obtained under each active observer.
- `proj_errs` contains `STP.Tpre * shapeErrs`, the errors presented to the diagonal PID channels. `active_projected_mask` identifies the four active rows.
- `Fx`, `Br`, and `Bz` contain observer values at the 25 allocated query positions; the `_direct` arrays contain direct FGE interpolation. Only positions identified by the saved shape masks and positive radial coordinates are physically meaningful. Unused neural-query rows can contain nonphysical values because the inference interface substitutes 1 mm for nonpositive radial coordinates.
- `shape_names_char`, `projected_names_char`, `sc_r`, `sc_z`, `sc_ref`, and the `sc_ls*` masks define the controlled quantities, their locations, and their roles. `Tpre`, `Tpost`, `SC2FB`, `Kp`, `Ki`, and `Kd` record the controller transformation and gains.
- `shape_hybrid`, `shape_current`, `shape_voltage`, and the corresponding `delta_*` arrays record the corrections produced by the outer loop. `coil_currents` and `coil_voltages` record realized plant actuation.
- `fge_final_flux` is the final flux map on `fge_grid_r` and `fge_grid_z`; the LCFS is the contour $\psi(R,Z)=\psi_B$ using `fge_final_boundary_flux`. `limiter_r` and `limiter_z` give the vessel limiter. `fge_boundary_r` and `fge_boundary_z` are the boundary-determining X-point or limiter point for each plasma domain, not a sampled LCFS curve.
- `fge_isconverged`, `fge_residual`, and `fge_iterations` record the numerical state of the FGE solution.

## Figure and quantitative result

The first panel shows all 11 active shape-error channels and reveals their transient and residual structure. The second shows the four projected coordinates acted on by the PID laws. The third reports the instantaneous RMS across the first ten saved error channels, and its legend gives the pooled post-activation RMS. The fourth separates that statistic by channel. The fifth plots the final LCFS as the boundary-flux contour together with the limiter and active control points. The sixth shows the Euclidean voltage-vector magnitude for both runs and the norm of their difference.

For $t\geq0.6$ s, the pooled RMS over the first ten shape-error rows is 2.75 mm with the neural observer and 3.60 mm with direct FGE interpolation. Both closed loops settle to millimetre-scale residuals in this deterministic linear scenario. The difference between these two values is a result for this initial condition, reference, and plant model; it is not by itself evidence of general observer or controller superiority.

## Scope and limitations

The experiment covers one linearized equilibrium with a fixed reference, no diagnostic noise, no profile perturbation, and no deliberately active actuator saturation. It does not determine robustness across equilibria, topologies, reference transients, or nonlinear plant evolution. The saved schema does not contain the design-time singular-value spectrum, geometric normal errors at each boundary point, exact component-wise controller saturation state, or measured controller execution time.

## Sources

- Dissertation section `Controller-level validation in FGE` in `/home/mg/repos/phd_thesis/2_eqrec.tex`.
- H. Anand et al., *A Novel Plasma Position and Shape Controller for Advanced Configuration Development on the TCV Tokamak*, local file `/home/mg/repos/phd_thesis/material/tcv_ctrl_infra/P27_Anand2017.pdf`.
- `/home/mg/repos/phd_thesis/guides/guides.md`, `hybrid_guide.md`, `fge_guide.md`, and `meq_infra_report.md`.
- `/home/mg/repos/hybrid/pid_shape_liuqenn_comparison.m`, `ctrls/shape/shapectrl.m`, `ctrls/shape/shapectrl_step.m`, `utils/plot/shape_ctrl_plot.m`, and the associated controller-building code.
- `/home/mg/repos/meq/meqpdom.m` and `meqpost.m`, used to establish the meaning of `rB,zB`.
