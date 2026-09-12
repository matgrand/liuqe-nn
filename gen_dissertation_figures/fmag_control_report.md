# FMag controller comparison

## Controller and simulation

FMag is a $10\,\mathrm{kHz}$ magnetic controller that constructs observables from magnetic measurements, subtracts references, mixes the resulting errors through `Tpre`, applies diagonal PID state-space blocks, distributes their outputs through `Tpost` into active-coil current-derivative requests, and converts those requests to voltage commands using the coil inductance and resistance model. Its controlled quantities include plasma current, slow radial and vertical current-centroid position, X-point directions, coil-current directions, and fast vertical motion through the G coil.

The nonlinear X-point observer uses a quadratic Taylor expansion of the flux around a design point and solves $\nabla\psi+H_\psi\Delta x=0$ for the estimated displacement $\Delta x=(\Delta r,\Delta z)$. The PINO replaces the six local quantities entering this calculation: $\psi$, $\partial_r\psi$, $\partial_z\psi$, $\partial_{rr}\psi$, $\partial_{rz}\psi$, and $\partial_{zz}\psi$. The linear current and position observers, controller states, transformations, and actuator path are common to both simulations. The boundary-flux value subtracted from the PINO X-point flux remains supplied by LIH.

The simulation uses shot 79742 from $0.6$ to $1.2\,\mathrm{s}$ with a $100\,\mu\mathrm{s}$ step, nonlinear FGE with `OhmTor_0D`, 100 passive-structure modes, one active separatrix X-point, and fixed controller references from the first FGE snapshot. The comparison observer is LIH. The dissertation figure gives full-window RMS Euclidean controller errors of $2.78\,\mathrm{mm}$ with the PINO and $2.84\,\mathrm{mm}$ with LIH.

## Figure

The collector is `/home/mg/repos/hybrid/fmag_liuqenn_comparison.m`, and the figure script is `fmag_control_simulation.py`. The figure shows radial and vertical X-point controller errors, their Euclidean norm, the PINO and LIH local displacement estimates against the FGE X-point displacement, Euclidean observer error relative to FGE, and the Hessian condition number $\kappa_2(H_\psi)$. On the PINO-driven trajectory, the maximum geometrical discrepancy relative to the FGE X-point is approximately $29.1\,\mathrm{mm}$ for the PINO local estimate and $11.2\,\mathrm{mm}$ for the LIH estimate. The output is `dissertation_figures/fmag_control_simulation.pdf`.

## Dissertation narrative

Introduce FMag and its controlled quantities before presenting the Taylor expansion. Explain how the local flux gradient and Hessian determine the X-point displacement, then state exactly which six quantities the PINO replaces and that LIH continues to supply the boundary flux. Read the figure from the radial and vertical feedback errors to their Euclidean aggregate, the local estimates against the FGE X-point, the geometrical discrepancy, and Hessian conditioning. Both controller configurations reduce their X-point feedback residuals with similar full-window RMS values, while the direct comparison with the simulated FGE geometry shows a larger displacement discrepancy for the PINO estimate in this case.

## Sources

- `/home/mg/repos/hybrid/fmag_liuqenn_comparison.m`, `ctrls/SCDalgofmag/FMag_step.m`, `FMag_control.m`, `FMag_observer_linearization.m`, and `Examples/build_RZIp_Xpoint_controller.m`.
- `/home/mg/repos/meq/meqpdom.m` and `meqpost.m` for the FGE equilibrium geometry.
