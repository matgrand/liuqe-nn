# MPC shape control

## Purpose and outputs

The data collector is `/home/mg/repos/hybrid/mpc_shape_liuqenn_comparison.m`, the figure script is `mpc_shape_control_simulation.py`, and the generated figure is `dissertation_figures/mpc_shape_control_simulation.pdf`. The collector performs matched closed-loop simulations using a neural magnetic observer and direct interpolation of the FGE flux map. It records the controller and estimator models, shape errors, complete QP decision vector, applied first move, optimizer timing, realized current margins, actuator signals, and final FGE geometry.

The six-panel figure shows all active shape errors, the six largest components of the applied first QP move, the first-ten-channel RMS, QP iteration counts at controller updates, the minimum realized E/F-coil-current margin, and the final LCFS with the limiter and control points. Solid curves denote the neural observer and dashed curves denote direct FGE interpolation with `qintmex`.

## Controller architecture

The MPC is an outer shape-control loop around TCV's inner hybrid magnetic controller. Its prediction model couples a linearized FGE plasma-and-circuit model to a state-space representation of hybrid and reduces the combined state to 30 variables. A Kalman filter estimates that state from the available magnetic and actuator signals.

The controller uses a velocity-form quadratic program. At each shape-control update it chooses increments to the outer-loop control vector over a prediction horizon of 15 samples and a control horizon of three samples. With a 2 ms shape period, these horizons correspond to 30 ms and 6 ms. The objective penalizes predicted shape error and control increments and supports E/F-coil-current constraints. Only the first move $\Delta u_{k|k}$ is applied to the hybrid references; the complete horizon is recomputed at the next update.

This configuration sets `svdtol = 0`, so `build_mpc_tunable_params.m` uses identity pre/post mappings rather than an SVD-truncated input basis. The controller has 16 input coordinates, and `duopt` contains 48 rows because it stores 16 increments for each of the three control-horizon moves. Rows 1--16 form the applied first move, while rows 17--48 are planned future moves.

## Simulation configuration

The simulation uses TCV shot 79742. FGE and hybrid execute with a $100\,\mu\mathrm{s}$ step, and the shape loop executes every 2 ms. Shape control switches on at 0.6 s. The Kalman filter then runs for ten shape samples before QP actuation begins at 0.62 s. The plant is a linearized FGE model, the isoflux reference is fixed and derived from FBT, and the allocator is disabled.

In the neural run, the controller evaluates `liuqenn_inference` from synthetic magnetic measurements. In the direct-observer run, the controller evaluates the simulated FGE flux map with `qintmex`. This direct path uses the simulation state and is not a software-in-the-loop execution of LIUQE-RT. Both runs share the controller matrices, Kalman-filter design, QP parameters, plant initialization, and shape reference.

## Saved data and interpretation

- `shape_errs` contains the regulated shape errors produced by the active observer, and `active_shape_mask` identifies the 14 active rows. The reported scalar metric uses the first ten rows.
- `qp_decision` contains all three planned moves. `qp_first_move` contains the 16-component move actually applied at each update. The figure selects the six components with the largest combined post-activation RMS across the two runs and labels them by controller-coordinate index because no physical channel scaling is defined for these stored coordinates.
- `qp_iterations` and `qp_objective_reduced` record solver diagnostics. The reduced `quadprog` objective omits the constant term and can be negative, so it is not an absolute physical performance measure. Iteration count measures algorithmic work, not wall-clock execution time.
- `mpc_update_mask` identifies the 2 ms optimizer events on the $100\,\mu\mathrm{s}$ debug time base. Controller signals are otherwise held between shape-loop updates.
- `kf_output_estimate`, `kf_reference_output`, and `kf_state` record the Kalman-filter trajectory. `MPC_Akf`, `MPC_Bkf`, `MPC_Ckf`, `MPC_Dkf`, and `MPC_Lkf` record its state-space model, while `MPC_Psi`, `MPC_Phi`, `MPC_H`, `MPC_f`, `MPC_FY`, and `MPC_FYPhi` record the optimization model.
- `ief_currents`, `ief_upper_limit`, `ief_lower_limit`, `ief_upper_margin`, and `ief_lower_margin` describe the realized E/F-coil currents and their distance from the configured bounds. A positive minimum margin indicates that every saved realized current remains within its bound.
- `fge_final_flux` is the final flux map on `fge_grid_r` and `fge_grid_z`; the LCFS is the contour $\psi(R,Z)=\psi_B$ using `fge_final_boundary_flux`. `fge_boundary_r` and `fge_boundary_z` identify the boundary-determining X-point or limiter point and are not LCFS coordinate sequences.
- `coil_currents`, `coil_voltages`, `shape_hybrid`, `shape_current`, and `shape_voltage` record the nested controller's command and actuator signals. `fge_isconverged`, `fge_residual`, and `fge_iterations` record plant-solver diagnostics.

## Figure and quantitative result

The first panel displays all 14 active regulated errors, with vertical markers distinguishing shape-loop activation from optimizer activation. The second panel shows the dominant components of the applied first move rather than future planned moves. The third gives the instantaneous ten-channel RMS and its pooled value after optimizer activation. The fourth shows QP iterations only at the saved optimizer events. The fifth shows the smallest realized distance to an E/F-current limit, and the sixth reconstructs the final LCFS from the saved boundary-flux contour.

For $t\geq0.62$ s, the pooled RMS over the first ten shape-error rows is 2.57 mm with the neural observer and 3.03 mm with direct FGE interpolation. The QP uses four or five iterations in both runs, and the minimum saved realized current margin is approximately 0.495 kA. Thus both observer configurations regulate the fixed target to millimetre-scale residuals with similar solver iteration counts while remaining inside the realized E/F-current bounds in this deterministic linear simulation.

## Scope and limitations

The saved update mask begins at 0.622 s, while the recorded iteration signal changes at 0.620 s, so the first optimization cycle is absent from the event-only iteration panel. The controller debug interface does not expose QP exit flags, solver-reported constraint violation, active predicted constraints, predicted output trajectories, the objective's constant term, or solve time. The solver is the local interior-point `quadprog` implementation; iteration counts cannot be interpreted as real-time feasibility or equated with the fixed-budget solver used in experimental MPC work.

The experiment covers one linearized equilibrium with a fixed reference, no diagnostic noise, no profile disturbance, and no deliberately active current constraint. It does not establish behavior for nonlinear FGE evolution, reference transients, active constraints, other equilibria, or other plasma topologies.

## Sources

- Dissertation section `Controller-level validation in FGE` in `/home/mg/repos/phd_thesis/2_eqrec.tex`.
- A. Mele et al., *First Experimental Demonstration of Plasma Shape Control in a Tokamak through Model Predictive Control*, local file `/home/mg/repos/phd_thesis/material/mele_mpc.pdf`.
- `/home/mg/repos/phd_thesis/guides/guides.md`, `hybrid_guide.md`, `fge_guide.md`, and `meq_infra_report.md`.
- `/home/mg/repos/hybrid/mpc_shape_liuqenn_comparison.m`, `ctrls/shape/shapectrl.m`, `ctrls/shape/MPC/shapempc_step.m`, `build_mpc_tunable_params.m`, and `utils/plot/shape_ctrl_plot.m`.
- `/home/mg/repos/meq/meqpdom.m` and `meqpost.m`, used to establish the meaning of `rB,zB`.
