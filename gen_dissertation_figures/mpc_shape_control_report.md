# MPC shape control

## Background and control objective

In an axisymmetric tokamak, the plasma boundary is the last closed magnetic-flux surface (LCFS). Its position and shape affect plasma-wall clearance, divertor strike-point placement, vertical stability, and confinement, so discrepancies between a nominal equilibrium design and the realized discharge must be corrected in feedback. Shape feedback has two coupled tasks: an equilibrium observer must infer the poloidal flux $\psi$ and field $(B_R,B_Z)$ from magnetic diagnostics, and a controller must alter the coil actuation so that the reconstructed boundary approaches the requested one. The simulation studied here isolates the effect of replacing direct interpolation of the simulated flux map with a neural magnetic observer while holding the plant, target, and MPC design fixed.

The target is formulated through isoflux constraints. Control points are placed on the desired boundary; for a diverted configuration, the controller drives the flux differences $\psi_i-\psi_X$ between those points and the primary X-point toward zero, while $B_R$ and $B_Z$ at the requested X-point are driven toward zero. Satisfying these local flux and field conditions makes a common flux contour pass through the requested locations and places the magnetic null, thereby controlling the LCFS without explicitly tracking a sampled boundary curve. In the implementation, `STP.SC2FB` converts the observer queries into a shape-feedback vector $y_{\mathrm{sh}}$, and the regulated error is

$$e_{\mathrm{sh}}=r_{\mathrm{sh}}-y_{\mathrm{sh}}.$$

This makes the equilibrium observer part of the closed loop: an observer error changes the control move and can therefore change the subsequent plasma trajectory.

## Model-predictive control formulation

Model-predictive control (MPC) is a receding-horizon feedback method. Instead of applying a fixed gain to the present error, it uses a dynamic model to predict how candidate control moves will affect future outputs, chooses the sequence that best balances tracking and actuation objectives while satisfying explicit constraints, applies only the first move, and repeats the optimization at the next sample using new measurements. Re-optimizing makes the strategy feedback rather than open-loop trajectory planning; predicting over several samples lets it coordinate strongly coupled actuators and anticipate limits before they are reached.

The TCV MPC is an outer shape loop around the faster inner `hybrid` magnetic controller. Its prediction plant is a linearized FGE plasma-and-circuit model is interconnected with a state-space model of `hybrid`, so the MPC predicts the response from modifications of the hybrid references to both shape quantities and E/F-coil currents. This preserves the established inner loop for fast vertical stabilization, current and position control, and voltage generation, while the MPC coordinates the slower shape task. The continuous combined model is discretized at the shape-controller period and reduced by Hankel balanced truncation to 30 states to make estimation and optimization tractable in real time.

Because the reduced dynamic state is not directly measured, a Kalman filter combines the model, the applied hybrid-reference correction, and the available shape, current, and magnetic measurements to estimate it. The filter is initialized from outputs sampled when shape control activates and is allowed to run for several samples before optimization begins. The controller uses a velocity-form model in the increments $\Delta x_k=x_k-x_{k-1}$ and $\Delta u_k=u_k-u_{k-1}$. Stacking predictions over a horizon gives the affine relation

$$Y_k=\Psi\Delta\hat x_k+\Phi\Delta U_k+\Omega y_{k-1}=\widehat Y_k+\Phi\Delta U_k,$$

where $Y_k$ contains the predicted outputs, $\Delta U_k$ contains the planned hybrid-reference increments over the control horizon, and $\widehat Y_k$ is the no-new-move prediction. The incremental formulation embeds the previous measured output in the prediction and provides integral-like, offset-free tracking for constant references and disturbances under the usual steady-state assumptions.

At each update the controller solves a quadratic program of the form

$$\min_{\Delta U_k}\;(R_k-Y_k)^{\mathsf T}W_y(R_k-Y_k)+\Delta U_k^{\mathsf T}W_{\Delta u}\Delta U_k \quad \text{subject to}\quad \underline I_{EF}\leq I_{EF,j}\leq\overline I_{EF},$$

for every predicted sample $j$ in the horizon. The first term penalizes predicted shape-tracking error and the second discourages abrupt or excessive changes to the hybrid references; optional weights can also trade shape accuracy against coil-current magnitude. Substituting the prediction equation produces the standard convex QP stored through the Hessian and linear term in the MAT files. The E/F-coil-current limits are represented directly as constraints on the predicted outputs.

Only the first optimized increment is used,

$$u_k=u_{k-1}+\Delta u_{k|k},$$

and all later planned increments are discarded when the next measurement arrives. For the present controller, the prediction horizon is 15 samples and the control horizon is three samples. At a 2 ms update period these correspond to 30 ms of prediction and 6 ms of independently optimized moves; inputs are held constant after the third planned move. With 16 control coordinates, the QP therefore has 48 decision variables. The paper reports that this horizon was chosen as a practical compromise between anticipation and real-time computational cost.

Unlike the generalized PID design, this simulation sets `svdtol = 0`: `build_mpc_tunable_params.m` uses identity pre/post maps instead of truncating the control space to dominant singular modes. Consequently, the 16 entries in the applied first move are controller coordinates passed to hybrid, and `duopt` has 48 rows arranged as three consecutive 16-component increments. The optimization supports a hard iteration budget and current constraints, but the local simulation uses MATLAB's interior-point `quadprog`; its iteration count is therefore only a diagnostic of this offline solve.

## Purpose and outputs

The data collector is `/home/mg/repos/hybrid/mpc_shape_liuqenn_comparison.m`, the figure script is `mpc_shape_control_simulation.py`, and the generated figure is `dissertation_figures/mpc_shape_control_simulation.pdf`. The collector performs matched closed-loop simulations using a neural magnetic observer and direct interpolation of the FGE flux map. It records the controller and estimator models, shape errors, complete QP decision vector, applied first move, optimizer timing, realized current margins, actuator signals, and final FGE geometry.

The six-panel figure shows all active shape errors, the six largest components of the applied first QP move, the first-ten-channel RMS, QP iteration counts at controller updates, the minimum realized E/F-coil-current margin, and the final LCFS with the limiter and control points. Solid curves denote the neural observer and dashed curves denote direct FGE interpolation with `qintmex`.

## Controller dimensions in this simulation

The controller has 16 input coordinates, and `duopt` contains 48 rows because it stores 16 increments for each of the three control-horizon moves. Rows 1--16 form the applied first move, while rows 17--48 are planned future moves.

## Simulation configuration

The simulation uses TCV shot 79742. FGE and hybrid execute with a $100\,\mu\mathrm{s}$ step, and the shape loop executes every 2 ms. Shape control switches on at 0.6 s. The Kalman filter then runs for ten shape samples before QP actuation begins at 0.62 s. The plant is a linearized FGE model, the isoflux reference is fixed and derived from FBT, and the allocator is disabled.

In the neural run, the controller evaluates `liuqenn_inference` from synthetic magnetic measurements. In the direct-observer run, the controller evaluates the simulated FGE flux map with `qintmex`. Both runs share the controller matrices, Kalman-filter design, QP parameters, plant initialization, and shape reference.

## Saved data and interpretation

- `shape_errs` contains the regulated shape errors produced by the active observer, and `active_shape_mask` identifies the 14 active rows. The reported scalar metric uses the first ten rows.
- `qp_decision` contains all three planned moves. `qp_first_move` contains the 16-component move actually applied at each update. The figure selects the six components with the largest combined post-activation RMS across the two runs and labels them by controller-coordinate index because no physical channel scaling is defined for these stored coordinates.
- `qp_iterations` and `qp_objective_reduced` record solver diagnostics. 
- `mpc_update_mask` identifies the 2 ms optimizer events on the $100\,\mu\mathrm{s}$ debug time base. Controller signals are otherwise held between shape-loop updates.
- `kf_output_estimate`, `kf_reference_output`, and `kf_state` record the Kalman-filter trajectory. `MPC_Akf`, `MPC_Bkf`, `MPC_Ckf`, `MPC_Dkf`, and `MPC_Lkf` record its state-space model, while `MPC_Psi`, `MPC_Phi`, `MPC_H`, `MPC_f`, `MPC_FY`, and `MPC_FYPhi` record the optimization model.
- `ief_currents`, `ief_upper_limit`, `ief_lower_limit`, `ief_upper_margin`, and `ief_lower_margin` describe the realized E/F-coil currents and their distance from the configured bounds. A positive minimum margin indicates that every saved realized current remains within its bound.
- `fge_final_flux` is the final flux map on `fge_grid_r` and `fge_grid_z`; the LCFS is the contour $\psi(R,Z)=\psi_B$ using `fge_final_boundary_flux`. `fge_boundary_r` and `fge_boundary_z` identify the boundary-determining X-point or limiter point.
- `coil_currents`, `coil_voltages`, `shape_hybrid`, `shape_current`, and `shape_voltage` record the nested controller's command and actuator signals. `fge_isconverged`, `fge_residual`, and `fge_iterations` record plant-solver diagnostics.

## Figure and quantitative result

The first panel displays all 14 active regulated errors, with vertical markers distinguishing shape-loop activation from optimizer activation. The second panel shows the dominant components of the applied first move. The third gives the instantaneous ten-channel RMS and its pooled value after optimizer activation. The fourth shows QP iterations only at the saved optimizer events. The fifth shows the smallest realized distance to an E/F-current limit, and the sixth reconstructs the final LCFS from the saved boundary-flux contour.

For $t\geq0.62$ s, the pooled RMS over the first ten shape-error rows is 2.57 mm with the neural observer and 3.03 mm with direct FGE interpolation. The QP uses four or five iterations in both runs, and the minimum saved realized current margin is approximately 0.495 kA. Thus both observer configurations regulate the fixed target to millimetre-scale residuals with similar solver iteration counts while remaining inside the realized E/F-current bounds in this deterministic linear simulation.

## Sources

- Dissertation section `Controller-level validation in FGE` in `/home/mg/repos/phd_thesis/2_eqrec.tex`.
- A. Mele et al., *First Experimental Demonstration of Plasma Shape Control in a Tokamak through Model Predictive Control*, local file `/home/mg/repos/phd_thesis/material/mele_mpc.pdf`.
- `/home/mg/repos/phd_thesis/guides/guides.md`, `hybrid_guide.md`, `fge_guide.md`, and `meq_infra_report.md`.
- `/home/mg/repos/hybrid/mpc_shape_liuqenn_comparison.m`, `ctrls/shape/shapectrl.m`, `ctrls/shape/MPC/shapempc_step.m`, `build_mpc_tunable_params.m`, and `utils/plot/shape_ctrl_plot.m`.
- `/home/mg/repos/meq/meqpdom.m` and `meqpost.m`, used to establish the meaning of `rB,zB`.
