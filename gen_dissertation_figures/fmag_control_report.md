# FMag controller comparison

## Current implementation and result

The extended collector is `/home/mg/repos/hybrid/fmag_liuqenn_comparison.m`, and the corresponding figure script is `fmag_control_simulation.py`. The intended comparison, the copied MAT files, and generated PDF use shot 79742 over 0.6--1.1999 s.

The collected files save the six raw local X-point quantities from the active observer and the LIH shadow, the exact plasma-current reference used to scale displacement errors, design and FGE X-point coordinates, controller matrices and state, observer masks and names, actuator signals, Hessian diagnostics, boundary-flux provenance, and FGE convergence diagnostics. The latest six-panel figure shows radial and vertical controller errors, their Euclidean norm, surrogate and LIH-shadow local displacements against the FGE X-point displacement, Euclidean observer error relative to FGE, and the Hessian condition number $\kappa_2(H_\psi)$ on a linear axis. The output is `dissertation_figures/fmag_control_simulation.pdf`.

## What the controller is and what this simulation tests

FMag is a 10 kHz magnetic controller that constructs observables from magnetic measurements, subtracts them from references, mixes the errors through `Tpre`, applies diagonal PID state-space blocks, distributes their outputs through `Tpost` into active-coil current-derivative requests, and converts those requests to voltage commands using the coil inductance and resistance model. In the tested controller, the mixed channels include slow radial and vertical current-centroid position, X-point control directions, coil-current directions orthogonal to position and X-point actuation, fast vertical stabilization through the G coil, plasma current through the OH coils, and a differential-OH direction.

The X-point observer is nonlinear. LIH supplies local flux and its first and second derivatives at the design X-point. FMag uses a quadratic Taylor expansion and solves $\nabla\psi+H_\psi\Delta x=0$ for the estimated displacement $\Delta x=(\Delta r,\Delta z)$. The neural path overwrites the six local quantities $\psi$, $\partial_r\psi$, $\partial_z\psi$, $\partial_{rr}\psi$, $\partial_{rz}\psi$, and $\partial_{zz}\psi$, then repeats the quadratic solve. It does not replace every LIH contribution: the boundary-flux value subtracted from the neural X-point flux remains LIH-based in the current `FMag_step.m`. The linear current, plasma-current, and position observers are identical in both runs.

The original data-generation script is `/home/mg/repos/hybrid/fmag_liuqenn_comparison.m`. It uses shot 79742 from 0.6 to 1.2 s with a 100 microsecond step, nonlinear FGE with `OhmTor_0D`, 100 passive-structure modes, one active separatrix X-point, no added noise, and fixed controller references taken from the first FGE snapshot. Unlike the PID and MPC cases, the standard observer is genuinely LIH.

## What the original files save and why

- `xpt_err_flux`, `xpt_err_delta_r`, and `xpt_err_delta_z` are the controller errors for the X-point output channels. FMag multiplies displacement-channel errors by the plasma-current reference before mixing, so the saved radial and vertical errors must both be divided by $I_p$ to recover meters. The previous Python plot divided both channels but introduced an additional minus sign only for the vertical channel; that sign change has no basis in `FMag_step.m`.
- `xpt_delta_r`, `xpt_delta_z`, and `xpt_flux` store the active observer outputs. The `_lih` variants store the LIH estimates evaluated on the same plant trajectory, allowing a shadow-observer comparison within the surrogate run. A shadow estimate is not the true FGE X-point position.
- `rY`, `zY`, and `Ip` store plant current-centroid position and plasma current interpolated onto the controller debug time. They show whether changing only the X-point observer perturbs the other regulated plant quantities.
- `Va` stores the 19 coil-voltage commands. It permits comparison of total actuation and run-to-run command differences, although the file does not save coil names, saturation flags, feedforward/feedback decomposition, or current-derivative requests.
- `mixed_errors` stores 32 allocated mixed-channel rows. Their names and scalings are not saved, and the channels combine incompatible physical quantities. The former hand-selected scaling factors therefore supported visualization only and could not support quantitative comparisons.
- `observables` and `observables_lih` store all 256 allocated observer rows. Their type masks and channel names are absent, so the extracted X-point fields are the safely interpretable subset.

Over the full saved window, the RMS Euclidean control error is 2.06 mm with the surrogate and 2.11 mm with LIH, matching the dissertation table. These values include the initial transient at 0.6 s. The later steady errors are smaller, so the time window must be stated whenever the metric is quoted.

## Level 1: improvement implemented from the existing files

`fmag_control_simulation.py` now converts both displacement-error channels with the same $e/I_p$ rule and removes the unexplained vertical sign reversal. The new six-panel figure shows radial and vertical X-point control errors, their Euclidean norm and full-window RMS, the surrogate and LIH shadow displacement estimates on the same surrogate-driven trajectory, run-to-run differences in current-centroid position, run-to-run plasma-current difference, and norms of the voltage command vectors and their difference. This replaces the arbitrary-unit mixed-error panel with quantities that have defined units or a defined vector norm. The generated PDF is `dissertation_figures/fmag_control_simulation.pdf`.

The shadow-observer panel is included to expose what is actually changed inside FMag, but it must not be read as observer error against ground truth. The saved data do not contain the FGE X-point trajectory needed for that comparison. The voltage-vector norm also summarizes actuation magnitude; it does not identify which coils account for the difference.

## Level 2: implemented collection and remaining changes

The collector implements the FGE and design X-point coordinates, fluxes, active-X-point mask, six raw local flux quantities for both observers, Hessian determinant and exact-zero flags, LIH boundary flux used by the nonlinear path, observer metadata, controller matrices and state, channel names, current derivatives, currents, voltages, references, feedforwards, and FGE convergence diagnostics proposed above. The Python script computes the spectral condition number from the three independent entries of the saved symmetric Hessian.

Exact internal saturation flags, anti-windup scaling, feedback and feedforward current-derivative components, per-cycle execution time, and broader observer-validity criteria are not exposed by the current wrapper and remain absent. Repeat the experiment with X-point reference steps, multiple equilibria and topologies, noise, and perturbations. A matched 2 cm X-point step would be informative if its reference and purpose are recorded explicitly.

In the current shot-79742 surrogate-driven trajectory, $\kappa_2(H_\psi)$ ranges from approximately 1.000 to 1.816 for the surrogate Hessian and from 1.007 to 1.057 for the LIH shadow. Both Hessians remain indefinite, as required for a saddle point, and no exact-zero determinant occurs. This linear, dimensionless plot replaces the earlier logarithmic $|\det H_\psi|$ panel, whose narrow range produced nonuniform-looking grid spacing and whose dimensional magnitude was not itself a condition number.

The FGE-truth panels change the interpretation of the controller result. The observer-defined controller residuals approach zero, but the FGE X-point moves relative to the fixed design point. The Euclidean discrepancy reaches approximately 29.1 mm for the surrogate local estimate and 11.2 mm for the LIH shadow in the current surrogate-driven run. This comparison is unambiguous here because the saved arrays contain one X-point; multi-X-point cases require explicit tracking of X-point identity.

The neural path should also be tested with an independently obtained boundary-flux estimate. At present the LIH boundary flux is retained, so the experiment isolates local X-point derivative replacement but does not demonstrate a fully LIH-independent X-point observer.

## Level 3: recommended dissertation narrative

Describe FMag's common and changed paths before presenting trajectories. The common path contains linear current and position observers and the complete PID/actuator chain; only six local quantities entering the nonlinear X-point solve change, while the boundary-flux subtraction remains LIH-based. This makes the isolation claim precise.

Read the figure from controller perception to external geometric verification: component feedback errors, their Euclidean aggregate, local radial and vertical estimates against the FGE X-point, Euclidean observer error, and Hessian conditioning. State that the 2.06 mm and 2.11 mm values include the activation transient and are observer-defined controller errors, not geometric errors against FGE. The supported conclusion is that the two controllers reduce their respective perceived X-point residuals with similar full-window RMS in this nonlinear simulation. The current data do not show equivalent tracking of the actual FGE X-point and do not establish observer accuracy over other configurations or active TCV operation.

## Sources read

- Dissertation section `Controller-level validation in FGE` in `/home/mg/repos/phd_thesis/2_eqrec.tex`.
- `/home/mg/repos/phd_thesis/guides/guides.md`, `hybrid_guide.md`, `fge_guide.md`, and `meq_infra_report.md`.
- `/home/mg/repos/hybrid/fmag_liuqenn_comparison.m`, `fmag_liuqenn_comparison.m`, `ctrls/SCDalgofmag/FMag_step.m`, `FMag_control.m`, `FMag_observer_linearization.m`, and `Examples/build_RZIp_Xpoint_controller.m`.
- The FGE references in `ref.bib` were checked. Their `localfile` fields are `-`, so no local FGE paper or thesis was available to read.
