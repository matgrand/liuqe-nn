# PID shape-control comparison

## Controller and simulation

The generalized TCV shape controller is an outer loop around the core hybrid magnetic controller. The target plasma boundary is represented by control points. Flux and poloidal-field values at those points are transformed by `STP.SC2FB` into the shape-feedback vector, which is subtracted from `SC.ref` to obtain the shape errors. `STP.Tpre` projects these errors into the directions handled by diagonal PID blocks, and `STP.Tpost` maps the controller outputs into corrections to the hybrid references. The SVD-based pre/post maps retain the directions that the closed-loop plant can control most effectively. This follows Anand et al., while the current implementation includes later changes in the `hybrid` repository.

The simulation uses TCV shot 79742, a linearized FGE plant, a $100\,\mu\mathrm{s}$ FGE and hybrid step, a $2\,\mathrm{ms}$ shape-controller period, a fixed shape reference, and controller activation at $0.6\,\mathrm{s}$. FGE supplies the simulated plasma, circuit evolution, magnetic measurements, and equilibrium. Direct `qintmex` interpolation of the FGE flux map provides the reference equilibrium quantities that LIUQE-RT would estimate in an experiment. The PINO run instead reconstructs the quantities required by the shape controller from the synthetic magnetic measurements.

Eleven of the 25 allocated shape-error rows and four of the 16 projected PID directions are active. The dissertation aggregate uses the first ten shape-error channels. For $t\geq0.6\,\mathrm{s}$, their pooled RMS is $2.75\,\mathrm{mm}$ with the PINO and $3.60\,\mathrm{mm}$ with direct FGE interpolation.

## Figure

The collector is `/home/mg/repos/hybrid/pid_shape_liuqenn_comparison.m`, and the figure script is `pid_shape_control_simulation.py`. The figure shows the active shape errors, four projected PID errors, their aggregate and per-channel RMS values, the final FGE LCFS with limiter and control points, and the coil-voltage vector norms. The LCFS is the contour $\psi(R,Z)=\psi_B$ obtained from the saved final flux grid. The output is `dissertation_figures/pid_shape_control_simulation.pdf`.

## Dissertation narrative

Introduce the complete signal chain before discussing the figure: FGE evolves the simulated equilibrium and produces magnetic measurements; either direct FGE interpolation or the PINO supplies flux and poloidal-field values at the boundary points; the shape transformation constructs the physical errors; the pre-transformation projects them into four active control directions; the PID blocks act on those directions; and the post-transformation sends corrections to the inner hybrid controller. Then read the figure from physical boundary errors to projected controller errors, aggregate regulation, final plasma geometry, and actuator commands. Both simulations regulate the selected shape errors to millimetre-scale residuals, with pooled RMS values of $2.75\,\mathrm{mm}$ for the PINO and $3.60\,\mathrm{mm}$ for direct FGE interpolation in this simulation.

## Sources

- H. Anand et al., *A Novel Plasma Position and Shape Controller for Advanced Configuration Development on the TCV Tokamak*, `/home/mg/repos/phd_thesis/material/tcv_ctrl_infra/P27_Anand2017.pdf`.
- `/home/mg/repos/hybrid/pid_shape_liuqenn_comparison.m`, `ctrls/shape/shapectrl.m`, `ctrls/shape/shapectrl_step.m`, `utils/plot/shape_ctrl_plot.m`, and associated controller-building code.
- `/home/mg/repos/meq/meqpdom.m` and `meqpost.m` for the FGE equilibrium geometry.
