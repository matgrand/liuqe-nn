# MPC shape-control comparison

## What the controller is and what this simulation tests

The MPC is an outer shape-control loop around TCV's inner hybrid controller. Its prediction model couples a linearized FGE plasma-and-circuit model to a state-space representation of hybrid, reduces the state to 30 variables, and estimates the state with a Kalman filter. In velocity form, the QP chooses increments to the outer-loop control vector over a prediction horizon of 15 samples and a control horizon of three samples. With the 2 ms shape period used here, these correspond to 30 ms and 6 ms. The objective penalizes predicted shape error and control increments, and the formulation can impose E/F-coil-current constraints. Only the first optimized move is applied; the horizon is recomputed at the next shape update.

The data-generation script is `/home/mg/repos/hybrid/mpc_shape_liuqenn_comparison.m`. It uses shot 79742, a fixed isoflux reference derived from FBT, a linearized FGE plant, no allocator, hybrid execution at 100 microseconds, and MPC updates at 2 ms. The shape loop switches on at 0.6 s, after which the Kalman filter runs for ten shape steps; the QP therefore starts at 0.62 s. As in the PID script, the comparison path is direct `qintmex` interpolation of the FGE flux map and is not a run of LIUQE-RT.

This particular configuration sets `svdtol = 0`. In `build_mpc_tunable_params.m`, that selects the identity pre/post mapping, not an SVD-truncated input basis. The saved `duopt` name is generic. Consequently, the former figure caption describing six retained SVD components is incorrect for these files.

## What is saved and why

- `shape_errs` and `shape_errs2` store the regulated errors obtained with the active observer and direct FGE interpolation. Fourteen of the 25 allocated rows are nonzero, while the dissertation metric uses the first ten.
- `Fx`, `Br`, `Bz` and their `*2` counterparts save active-observer and direct-FGE values at the query slots. Unused surrogate query rows contain nonphysical values and must be masked using controller metadata that is currently absent.
- `duopt` has shape $48\times3500$ because it stores $16$ controller-input increments for each of three control-horizon moves. Rows 1--16 are $\Delta u_{k|k}$, the first move applied by `shapempc_step`; rows 17--32 and 33--48 are future planned moves. The old Python figure used the final six rows, which belong to the third planned move, and labelled them as six SVD modes.
- `niter` stores QP solver iterations and `mpc_cost` stores the `quadprog` objective value. Both are held between 2 ms updates on the 100 microsecond debug time base. The saved cost omits the objective's constant term and can be negative, so it is not an absolute physical performance metric. Iteration count is not measured execution time.
- `currs` and `rz` omit channel names and transformation metadata; `rz` values also show that they are hybrid weighted signals, not directly plot-ready positions in meters.

For $t\geq0.62$ s, the pooled RMS of the first ten errors is 2.57 mm for the surrogate and 3.03 mm for direct FGE interpolation, reproducing the dissertation table. There are 140 distinct MPC updates in the saved post-delay interval. The solver uses four or five iterations in both runs.

## Level 1: improvement implemented from the existing files

`mpc_shape_control_simulation.py` now marks both the 0.60 s shape-loop switch-on and the 0.62 s optimizer start. It shows the first ten shape errors in millimeters, their instantaneous RMS with the reported post-0.62 s aggregate, the six largest saved components selected from the applied first move only, and the QP iteration count at distinct 2 ms updates. The decision-vector axis says `stored units` because the files do not contain the channel map and scaling required for a physical unit. The baseline is labelled `FGE/qintmex`. The generated PDF is `dissertation_figures/mpc_shape_control_simulation.pdf`.

The first-move panel deliberately does not call its entries SVD modes. It uses component indices and explains the selection criterion. This is less physically descriptive than named control directions, but it is correct for the available files and exposes the part of the solution that actually enters the controller.

## Level 2: changes recommended for the MATLAB collection script

Save the full MPC definition needed to interpret `duopt`: `FP.Ny`, `FP.Nu`, `FP.nShapeCtrl`, `svdtol`, `MTP.Tpre`, `MTP.Tpost`, channel names, weights, reference scaling, and the map from decision coordinates to hybrid references. Save the applied first move separately, together with the resulting hybrid-reference vector and coil-voltage command. Save QP exit flags, constraint violation, active constraints, predicted E/F currents, current limits, predicted shape trajectories, Kalman innovations, and measured solve time at each 2 ms update. This would permit figures of constraint margin, prediction versus realization, and actual timing.

Save only one optimizer record per MPC update or save an explicit update mask. Also save the complete objective, including the constant tracking-error term, if objective values are to be compared across observer runs. The data should record the exact solver and options because the local MATLAB path currently uses interior-point `quadprog`, while the cited real-time work discusses a different fixed-budget implementation.

For a stronger observer comparison, repeat the matched test with nonlinear FGE, reference transients, active current constraints, diagnostic noise, profile disturbances, and multiple equilibria. The current fixed-reference, unconstrained-looking trajectory mainly shows regulation and optimizer sensitivity in one linear case.

## Level 3: recommended dissertation narrative

Explain the nested architecture first: MPC does not command coil voltages directly; it modifies references supplied to hybrid, which closes the fast magnetic loop. Then state the 2 ms update period, 15-step prediction horizon, three-step control horizon, Kalman-filter warm-up, and the distinction between the switch-on time and optimizer start.

The figure discussion should proceed from regulated errors to the applied QP move and then to solver behavior. State that the decision components are unnamed in the saved data and that only their indices can be shown. Use the RMS panel to connect the trajectories to the table. The supported conclusion is that replacing direct FGE interpolation with the surrogate leads to similar regulation and a QP requiring the same four-to-five iteration range in this simulation; the data do not provide real-time solver latency or a test of constraint handling.

## Sources read

- Dissertation section `Controller-level validation in FGE` in `/home/mg/repos/phd_thesis/2_eqrec.tex`.
- A. Mele et al., *First Experimental Demonstration of Plasma Shape Control in a Tokamak through Model Predictive Control*, local file `/home/mg/repos/phd_thesis/material/mele_mpc.pdf`.
- `/home/mg/repos/phd_thesis/guides/guides.md`, `hybrid_guide.md`, and `fge_guide.md`.
- `/home/mg/repos/hybrid/mpc_shape_liuqenn_comparison.m`, `ctrls/shape/shapectrl.m`, `ctrls/shape/MPC/shapempc_step.m`, `build_mpc_tunable_params.m`, and the associated fixed-parameter definitions.

