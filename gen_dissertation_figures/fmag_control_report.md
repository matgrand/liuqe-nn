# FMag X-point control

## Purpose

The data collector is `/home/mg/repos/hybrid/fmag_liuqenn_comparison.m`, the figure script is `fmag_control_simulation.py`, and the generated figure is `dissertation_figures/fmag_control_simulation.pdf`. The experiment compares the closed-loop behavior of the same FMag controller and nonlinear FGE plant when its X-point observer is supplied either by LiuqeNN or by the standard LIH linear observer. The observer comparison is intentionally limited to the quantities used by this FMag configuration: it is not a replacement of every equilibrium or controller signal.

This report has two goals. First, it explains the FMag control law from its theoretical observable-feedback structure through the implementation that produces coil voltages. Second, it identifies exactly where an observer enters this configuration, which quantities LiuqeNN replaces, which paths remain common, and which saved signals establish whether the replacement remains accurate in closed loop.

## FMag theory and implementation

FMag runs with a controller interval of $100\,\mu\mathrm{s}$, corresponding to 10 kHz. Its measurement vector contains flux-loop signals, magnetic-probe signals, active-coil currents, toroidal-current input, flux-loop voltages, and selected time derivatives of magnetic-probe signals. A static observer matrix first maps this measurement vector to linear observable channels. Nonlinear processing then converts selected local flux quantities into geometrically meaningful X-point or shape displacements.

For an X-point with design position $x_0=(R_0,Z_0)$, FMag uses the second-order local model $\psi(x_0+\delta x)\simeq\psi_0+g^{\mathsf T}\delta x+\tfrac{1}{2}\delta x^{\mathsf T}H\delta x$, where $g=(\partial_R\psi,\partial_Z\psi)$ and $H$ contains $\partial_{RR}\psi$, $\partial_{RZ}\psi$, and $\partial_{ZZ}\psi$. Setting the model gradient to zero gives the Newton displacement $\delta x=-H^{-1}g$. Substitution into the quadratic model gives the estimated absolute flux at the displaced saddle, $\hat\psi_X=\psi_0+\tfrac{1}{2}g^{\mathsf T}\delta x$.

The boundary flux is estimated at fixed candidate positions sampled from the design last closed flux surface. For the plasma-current sign $s_I=\operatorname{sign}(I_p)$, FMag selects $\hat\psi_B=s_I\max_i(s_I\hat\psi_i)$. The X-point flux observable is the relative value $\hat\psi_X-\hat\psi_B$. The candidate coordinates are fixed before either simulation and are identical in the LIH and LiuqeNN runs, so the comparison changes the estimator rather than the geometry supplied to it.

FMag subtracts its observables from their references, scales geometric displacement errors by the plasma-current reference, and mixes the errors through `Tpre`. Independent diagonal PID state-space blocks act on the mixed coordinates. `Tpost` maps their outputs into active-coil current-derivative requests, after which the coil resistance and inductance model computes voltage commands. The implementation also provides controller switching, bumpless transfer, voltage saturation, and anti-windup.

The controller built by `build_RZIp_Xpoint_controller` contains slow $rI_p$ and $zI_p$ regulation, X-point directions, coil-current directions orthogonal to the position and X-point actuation, fast vertical stabilization through the G coil, plasma-current control through the OH coils, and a differential-OH direction. In the selected shot there is one active X-point. The controller construction removes the common X-point-flux direction, leaving two active X-point displacement modes; correspondingly, the absolute-flux channel has zero weight in `Tpre`. A consistent boundary estimate is still required to define and validate the complete X-point observable, but changing it is not expected to change coil commands in this single-X-point case. A causal test of boundary-flux feedback would require a separate configuration with multiple controlled X-points or an explicitly retained absolute-flux mode.

## Observer substitution

In the LIH run, the local X-point flux, its two first derivatives, its three independent second derivatives, and all boundary-candidate fluxes are linear functions of the simulated diagnostic vector constructed from LIH observer matrices. The nonlinear X-point calculation then uses these LIH quantities without further access to the FGE equilibrium solution.

In the LiuqeNN run, one vectorized neural inference evaluates the design X-point and the fixed boundary candidates. LiuqeNN supplies $\psi$, $\partial_R\psi$, $\partial_Z\psi$, $\partial_{RR}\psi$, $\partial_{RZ}\psi$, and $\partial_{ZZ}\psi$ at the design X-point and supplies $\psi$ at every boundary candidate. FMag forms both $\hat\psi_X$ and $\hat\psi_B$ from those neural outputs. Both runs retain the established LIH-based bumpless initialization, after which the LiuqeNN branch replaces the active X-point and boundary-flux quantities from the first normal controller step; using neural quantities inside the initialization iteration produced a large controller-state offset and an unstable FGE trajectory, so it is intentionally outside this comparison.

The remaining paths are intentionally unchanged. Active-coil current channels are direct measurements, while $I_p$, $rI_p$, $zI_p$, and the fast vertical observable retain their existing linear observers. The current-density distribution is not an input to this FMag control law. The experiment therefore tests substitution of the complete local X-point and boundary-flux observer used by this controller, not substitution of every state estimate in FMag or FGE.

## Simulation configuration

The default comparison uses TCV shot 79742 from 0.6 s to 1.2 s with a $100\,\mu\mathrm{s}$ simulation and controller step; shot 78071 is retained as the commented alternative in the collector. The nonlinear FGE plant uses the `OhmTor_0D` current-density evolution model, 100 passive-structure modes, the E, F, OH, and G active-coil groups, and the active separatrix X-point selected from the initial equilibrium. No artificial diagnostic noise is added. The plant initialization, controller matrices, references, feedforwards, and boundary-candidate coordinates are shared between the two runs; the active observer source is the intended difference.

FGE supplies the simulated plant evolution and the diagnostic signals seen by the controller. Its reconstructed X-point coordinates and boundary flux are saved only as external truth for analysis; they are not passed to either online observer during the simulation.

After each closed-loop simulation has finished, the collector also runs full nonlinear LIUQE as a shadow observer on every `LIUQE_SAMPLE_STRIDE`-th FGE sample, with a default stride of 100. The shadow receives that run's simulated flux loops, magnetic probes, active-coil currents, flux-loop voltages, plasma current, and toroidal-field signal. If FGE supplies a simulated diamagnetic-flux signal it is also substituted; otherwise the LIUQE diamagnetic constraint is disabled instead of retaining unrelated experimental data. This sparse calculation never feeds FMag and therefore cannot change either closed-loop trajectory. Its purpose is to distinguish disagreement with FGE from disagreement with the full-LIUQE reconstruction on which LiuqeNN was trained.

## Saved data

Schema version 4 saves the controller errors and state, active and LIH-shadow observable vectors, the six local X-point Taylor quantities, the absolute and boundary-subtracted X-point fluxes, boundary-candidate coordinates and flux values, the active boundary estimate, FGE geometric and flux truth, controller transformations and references, coil currents and current derivatives, voltages and voltage-limit indicators, FGE convergence diagnostics, and the sparse full-LIUQE shadow reconstruction.

- `xpt_err_delta_r` and `xpt_err_delta_z` are the scaled errors perceived by the active controller. Dividing by `ip_scaling_reference` recovers displacement errors in metres.
- `xpt_delta_r` and `xpt_delta_z` are the active observer displacement estimates. Their comparison with the corresponding run's `fge_active_xpoint_r` and `fge_active_xpoint_z` measures geometric observer accuracy.
- `boundary_query_r` and `boundary_query_z` contain the fixed design-LCFS candidates. `boundary_flux_candidates` contains their active observer estimates, while `boundary_flux_candidates_lih` retains the LIH shadow on the same plant trajectory.
- `boundary_flux_used` is the signed extremum actually subtracted by the active observer. `boundary_flux_lih_shadow` records the LIH counterpart.
- `xpt_flux_absolute` contains the Taylor-estimated saddle flux before subtraction, and `xpt_flux` contains $\hat\psi_X-\hat\psi_B$. The `_lih` counterparts retain the LIH calculation.
- `fge_active_xpoint_flux-fge_boundary_flux` is the external FGE reference for the relative X-point flux.
- `Tpre` confirms which observed channels influence the retained control coordinates. In this single-X-point experiment its X-point-flux column is zero, while the radial and vertical displacement columns are active.
- `liuqe_sample_stride`, `liuqe_sample_indices`, and `liuqe_t` identify the sparse shadow samples. `liuqe_active_xpoint_r`, `liuqe_active_xpoint_z`, and `liuqe_active_xpoint_flux` are the full-LIUQE X-point nearest the design X-point, while `liuqe_boundary_flux` is its reconstructed boundary flux. `liuqe_uses_dml` records whether a simulated diamagnetic measurement was available.

## Figure interpretation

The figure proceeds from closed-loop performance to observer validation and actuator consequences. The first two panels show the signed radial and vertical residuals perceived by the controller and their Euclidean norm. The central four panels separate the LiuqeNN-controlled and LIH-controlled trajectories and, for each radial and vertical coordinate, show the FGE plant state, the online observer estimate, and the sparsely sampled full-LIUQE shadow estimate. This makes the apparent neural drift interpretable: agreement between LiuqeNN and full LIUQE would indicate chiefly a LIUQE-versus-FGE model difference, whereas disagreement between LiuqeNN and both full LIUQE and FGE would indicate surrogate error. The seventh panel applies the same comparison to boundary-subtracted X-point flux, and the final panel compares the coil-voltage vector norms and their difference.

The previous schema-2 and schema-3 MAT files must not be used for this figure. Schema 2 contained no active boundary candidates and stored a zero value under the misleading name `lih_boundary_flux_used`; schema 3 corrected the boundary observer but did not contain the full-LIUQE shadow. Quantitative claims for the new three-way comparison must be filled from newly generated schema-4 files after the TCV run.

## Scope

This is a deterministic single-shot, single-X-point test with fixed target geometry. It establishes whether the LiuqeNN X-point and boundary-flux observables can replace the corresponding LIH quantities without losing closed-loop displacement regulation in this scenario. It does not by itself test diagnostic noise, other equilibria or topologies, reference changes, profile perturbations, deliberately active actuator limits, execution time, or the causal influence of boundary-flux feedback.

## Sources

- `/home/mg/repos/phd_thesis/guides/guides.md`, `hybrid_guide.md`, `fge_guide.md`, `liuqe_guide.md`, and `meq_infra_report.md`.
- `/home/mg/repos/hybrid/fmag_liuqenn_comparison.m`, `ctrls/SCDalgofmag/FMag_step.m`, `FMag_control.m`, `FMag_init.m`, `FMag_observer_linearization.m`, `get_fmag_xpoint_observer.m`, and `Examples/build_RZIp_Xpoint_controller.m`.
