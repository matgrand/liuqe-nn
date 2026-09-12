# MPC shape-control comparison

## Controller and simulation

The isoflux MPC is an outer shape-control loop around TCV's inner hybrid controller. Its prediction model couples a linearized FGE plasma-and-circuit model to a state-space representation of hybrid and uses a Kalman filter to estimate a 30-state reduced model. In velocity form, the quadratic program chooses increments to the outer-loop control vector over a prediction horizon of 15 samples and a control horizon of three samples. With the $2\,\mathrm{ms}$ shape period used here, these correspond to $30\,\mathrm{ms}$ and $6\,\mathrm{ms}$. The objective penalizes predicted shape error and control increments and can include E/F-coil-current constraints. Only the first optimized move is applied before the horizon is recomputed.

The simulation uses shot 79742, a fixed isoflux reference derived from FBT, a linearized FGE plant, hybrid execution at $100\,\mu\mathrm{s}$, and MPC updates at $2\,\mathrm{ms}$. The shape loop switches on at $0.60\,\mathrm{s}$, the Kalman filter runs for ten shape steps, and the quadratic program starts at $0.62\,\mathrm{s}$. FGE supplies the simulated plasma, circuit evolution, magnetic measurements, and equilibrium. Direct `qintmex` interpolation of the FGE map supplies the reference equilibrium quantities that LIUQE-RT would estimate in an experiment, while the PINO reconstructs the quantities required by the controller from synthetic magnetic measurements.

For $t\geq0.62\,\mathrm{s}$, the pooled RMS of the first ten shape errors is $2.57\,\mathrm{mm}$ with the PINO and $3.03\,\mathrm{mm}$ with direct FGE interpolation. The quadratic-program solver uses four or five iterations in both simulations, and the minimum realized E/F-coil-current margin is approximately $0.495\,\mathrm{kA}$.

## Figure

The collector is `/home/mg/repos/hybrid/mpc_shape_liuqenn_comparison.m`, and the figure script is `mpc_shape_control_simulation.py`. The figure shows the active shape-error trajectories, the six largest components of the applied first QP move using neutral controller-component labels, the first-ten-channel RMS, solver iterations at controller updates, minimum realized E/F-current margin, and the final FGE LCFS with limiter and control points. The output is `dissertation_figures/mpc_shape_control_simulation.pdf`.

## Dissertation narrative

Introduce the nested architecture before discussing the figure: MPC does not command coil voltages directly; it modifies the references supplied to hybrid, which closes the fast magnetic loop. Explain the state estimate, $2\,\mathrm{ms}$ update period, 15-step prediction horizon, three-step control horizon, first-move receding-horizon implementation, and the distinction between shape-loop activation and optimizer start. Then read the figure from regulated shape errors to applied control increments, aggregate regulation, solver behavior, realized current margin, and final plasma geometry. Replacing direct FGE interpolation with the PINO produces similar regulation and the same four-to-five QP iteration range, with pooled RMS values of $2.57\,\mathrm{mm}$ for the PINO and $3.03\,\mathrm{mm}$ for direct FGE interpolation.

## Sources

- A. Mele et al., *First Experimental Demonstration of Plasma Shape Control in a Tokamak through Model Predictive Control*, `/home/mg/repos/phd_thesis/material/mele_mpc.pdf`.
- `/home/mg/repos/hybrid/mpc_shape_liuqenn_comparison.m`, `ctrls/shape/shapectrl.m`, `ctrls/shape/MPC/shapempc_step.m`, `build_mpc_tunable_params.m`, and `utils/plot/shape_ctrl_plot.m`.
- `/home/mg/repos/meq/meqpdom.m` and `meqpost.m` for the FGE equilibrium geometry.
