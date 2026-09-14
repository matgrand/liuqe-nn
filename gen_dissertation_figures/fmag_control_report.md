# FMag X-point control

## Purpose and outputs

The data collector is `/home/mg/repos/hybrid/fmag_liuqenn_comparison.m`, the figure script is `fmag_control_simulation.py`, and the generated figure is `dissertation_figures/fmag_control_simulation.pdf`. The experiment runs two matched closed-loop simulations so that the nonlinear X-point observer can be evaluated with either neural local-flux quantities or the standard LIH quantities while the controller structure, references, plant initialization, and simulation timing remain fixed.

The MAT files contain the controller errors and state, observer outputs and metadata, the six local flux quantities used by the nonlinear X-point calculation, FGE geometry and convergence signals, controller transformations, references, coil currents, current derivatives, voltages, voltage-limit indicators, and plant quantities. The figure presents the radial and vertical X-point control errors, their Euclidean norm, the inferred radial and vertical displacement alongside the FGE X-point displacement, the Euclidean observer error relative to FGE, and the spectral condition number of the local flux Hessian.

## Controller architecture

FMag executes at 10 kHz. It constructs observable signals from magnetic measurements, subtracts those signals from their references, mixes the resulting errors through `Tpre`, applies diagonal PID state-space blocks, and distributes the PID outputs through `Tpost` into active-coil current-derivative requests. A coil resistance and inductance model converts those requests to voltage commands. The mixed control coordinates cover slow radial and vertical current-centroid position, X-point control directions, coil-current directions orthogonal to the position and X-point actuation, fast vertical stabilization through the G coil, plasma-current control through the OH coils, and a differential-OH direction.

The X-point observer uses a quadratic local model of the poloidal flux. Given the flux gradient $\nabla\psi$ and Hessian $H_\psi$ at the design X-point, it solves $\nabla\psi+H_\psi\Delta x=0$ for the displacement $\Delta x=(\Delta r,\Delta z)$. The neural observer supplies $\psi$, $\partial_r\psi$, $\partial_z\psi$, $\partial_{rr}\psi$, $\partial_{rz}\psi$, and $\partial_{zz}\psi$ from synthetic magnetic measurements. The LIH observer supplies the corresponding local quantities from the linearized plasma response. In the neural configuration, the boundary-flux value subtracted from the X-point flux remains LIH-based; the neural inference therefore replaces the local X-point quantities but not the boundary-flux observer. The linear coil-current, plasma-current, and current-centroid observers are common to both configurations.

## Simulation configuration

The simulation uses TCV shot 78071 from 0.6 s to 1.2 s with a $100\,\mu\mathrm{s}$ step. The nonlinear FGE plant uses the `OhmTor_0D` current-density evolution model, 100 passive-structure modes, the E, F, OH, and G active-coil groups, and one active separatrix X-point selected from the initial equilibrium. No measurement noise is added. Controller references are produced by the FMag controller builder from the initial FGE state and are shared by the two simulations.

The LIH simulation sets `useliunn = false`. The neural simulation sets `useliunn = true` and provides the selected design X-point coordinates to the inference path. Each simulation starts from the same FGE model and equilibrium trajectory definition. This matched construction isolates the closed-loop consequences of the nonlinear X-point observer choice.

## Saved data and interpretation

- `xpt_err_flux`, `xpt_err_delta_r`, and `xpt_err_delta_z` are the controller errors in the three nonlinear X-point output channels. FMag scales each displacement error by the plasma-current reference before error mixing, so `xpt_err_delta_r / ip_scaling_reference` and `xpt_err_delta_z / ip_scaling_reference` recover the radial and vertical errors in metres.
- `xpt_delta_r`, `xpt_delta_z`, and `xpt_flux` contain the active observer outputs. Their `_lih` counterparts contain an LIH shadow estimate evaluated on the same plant trajectory. The shadow signals diagnose the local observer calculation without constituting an independent FGE truth signal.
- `xpt_raw_F`, `xpt_raw_dFdr`, `xpt_raw_dFdz`, `xpt_raw_d2Fdrr`, `xpt_raw_d2Fdrdz`, and `xpt_raw_d2Fdzz`, together with their `_lih` counterparts, contain the local Taylor-model inputs. `xpt_hessian_det` and `xpt_hessian_singular` provide direct degeneracy diagnostics.
- `fge_active_xpoint_r`, `fge_active_xpoint_z`, and `fge_active_xpoint_flux` contain the active X-point extracted from FGE. Subtracting `design_xpoint_r` and `design_xpoint_z` gives the geometric displacement used as external truth in the figure.
- `references`, `feedforwards`, `coil_current_references`, `Tpre`, `Tpost`, the PID state-space matrices, and `controller_state` record the complete control transformation and dynamic state needed to interpret the controller signals.
- `coil_voltages`, `coil_current_derivative`, `coil_currents`, and `voltage_at_limit` describe actuator demand and voltage-limit contact. `current_centroid_r`, `current_centroid_z`, and `plasma_current` record other regulated plant quantities.
- `fge_isconverged`, `fge_residual`, and `fge_iterations` expose the numerical status of the nonlinear plant solution.

The allocated `mixed_errors` array combines coordinates with different dimensions and scalings; it should only be interpreted together with `mixed_names_char` and the controller transformations. A norm of the coil-voltage vector is a useful aggregate actuation measure, but it does not identify the coils responsible for a change.

## Figure and quantitative result

The six-panel figure is read from controller residuals to geometric validation. The first panel shows the signed radial and vertical controller errors. The second shows $\|e_X\|_2$ and reports its full-window RMS. The third and fourth panels compare the observer displacement estimates with the displacement of the FGE X-point from its design location. The fifth gives the Euclidean displacement error relative to FGE. The sixth gives $\kappa_2(H_\psi)$, which directly measures the sensitivity of the quadratic X-point solve to perturbations in the gradient.

Over the complete saved interval, the RMS Euclidean controller error is 2.06 mm with the neural observer and 2.11 mm with LIH. These controller residuals include the activation transient and quantify the error perceived by each observer. They are distinct from the geometric error relative to the FGE X-point. In the neural-observer trajectory, the maximum Euclidean discrepancy from FGE is approximately 29.1 mm for the neural local estimate and 11.2 mm for the LIH shadow estimate. The Hessian condition number ranges from approximately 1.000 to 1.816 for the neural estimate and from 1.007 to 1.057 for the LIH shadow; both Hessians remain indefinite, as expected at a saddle point, and neither has an exactly zero determinant.

The experiment therefore shows that both observer configurations drive their perceived X-point residuals to similar millimetre-scale levels in this nonlinear simulation, while their local displacement estimates differ materially from the FGE X-point trajectory. Because the neural path retains an LIH boundary-flux estimate, the result characterizes replacement of the six local X-point quantities rather than a fully LIH-independent observer.

## Scope and limitations

This is a deterministic single-shot test with a fixed X-point target and one active X-point. It does not establish behavior across equilibria, plasma topologies, reference steps, diagnostic noise, profile perturbations, or actuator saturation. Exact internal anti-windup scaling, feedback/feedforward current-derivative components, per-cycle execution time, and general observer-validity flags are not exposed by the controller wrapper. Multiple-X-point studies also require an explicit rule for tracking X-point identity over time.

## Sources

- Dissertation section `Controller-level validation in FGE` in `/home/mg/repos/phd_thesis/2_eqrec.tex`.
- `/home/mg/repos/phd_thesis/guides/guides.md`, `hybrid_guide.md`, `fge_guide.md`, and `meq_infra_report.md`.
- `/home/mg/repos/hybrid/fmag_liuqenn_comparison.m`, `ctrls/SCDalgofmag/FMag_step.m`, `FMag_control.m`, `FMag_observer_linearization.m`, and `Examples/build_RZIp_Xpoint_controller.m`.
