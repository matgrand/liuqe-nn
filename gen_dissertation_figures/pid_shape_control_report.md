# PID shape control

## Background and control objective

In an axisymmetric tokamak, the plasma boundary is the last closed magnetic-flux surface (LCFS). Its position and shape affect plasma-wall clearance, divertor strike-point placement, vertical stability, and confinement, so errors in the programmed poloidal-field-coil currents cannot in general be left to feedforward control alone. Shape feedback has two coupled tasks: an equilibrium observer must infer the poloidal flux $\psi$ and field $(B_R,B_Z)$ from magnetic diagnostics, and a controller must alter the coil actuation so that the reconstructed boundary approaches the requested one. The simulation studied here isolates the effect of replacing the direct FGE flux-map interpolation used as its baseline observer with a neural magnetic observer while leaving the plant, target, and controller design unchanged.

The target is expressed through the isoflux principle. A set of control points is placed on the desired boundary. For a diverted plasma, the primary X-point supplies the reference flux $\psi_X$, and boundary or divertor-leg points are controlled by driving $\psi_i-\psi_X$ to zero; the X-point itself is positioned by driving $B_R$ and $B_Z$ to zero at its requested location. For a limited plasma, flux differences between boundary points play the corresponding role. If all these conditions are satisfied locally, the same flux contour passes through the requested points and therefore approximates the desired boundary. The controller feedback vector $y_{\mathrm{sh}}$ is assembled from the observer queries as

$$y_{\mathrm{sh}} = T_{\mathrm{SC}\rightarrow\mathrm{FB}}\begin{bmatrix}\psi & B_R & B_Z\end{bmatrix}^{\mathsf T}, \qquad e_{\mathrm{sh}} = r_{\mathrm{sh}}-y_{\mathrm{sh}},$$

where $T_{\mathrm{SC}\rightarrow\mathrm{FB}}$ is stored as `STP.SC2FB`. This construction turns flux and field values with different geometrical roles into a single vector of shape errors. It also makes the observer an integral part of the feedback loop: an observer error changes the error seen by the controller and can consequently change the entire closed-loop plasma trajectory.

## Generalized multivariable PID formulation

A tokamak shape controller is intrinsically multiple-input, multiple-output: every E/F-coil combination changes the flux at many control points, and each shape error can require several actuators. Applying an independent PID directly to every physical error would therefore ignore strong cross-coupling and would attempt to control directions to which the plant may respond only weakly. Following Anand et al., the generalized controller uses a linearized response about a reference equilibrium to identify well-conditioned joint directions in error and actuation space before applying scalar feedback laws.

Let $G_0$ denote the static response from corrections to the inner hybrid-controller references to the shape-feedback vector, and let $W_y$ and $W_u$ weight the importance of shape quantities and the availability of actuation directions. The controller design forms the weighted sensitivity and takes its singular-value decomposition,

$$W_yG_0W_u = U\Sigma V^{\mathsf T}.$$

The left singular vectors describe orthogonal combinations of shape errors, the right singular vectors give the corresponding combinations of hybrid-reference corrections, and each singular value measures how effectively that input direction changes its paired output direction. Retaining the $n_r$ largest singular directions rejects weakly controllable combinations that would demand disproportionate actuator effort or amplify reconstruction and modelling errors. In the present implementation the transformations are, schematically,

$$e_m = T_{\mathrm{pre}}e_{\mathrm{sh}}, \qquad T_{\mathrm{pre}}=U_r^{\mathsf T}W_y, \qquad \Delta r_h=T_{\mathrm{post}}c, \qquad T_{\mathrm{post}}=W_uV_r,$$

with zero padding for discarded directions. Unlike a literal static pseudoinverse, the current implementation deliberately does not include $\Sigma_r^{-1}$ in these transformations; the modal loop gains are tuned from the dynamics of each transformed channel instead. For this simulation `svdtol = 4` means that four dominant singular directions are retained. The resulting projected errors are not individual geometrical displacements: each is a weighted linear combination of the physical isoflux and field errors.

Each retained modal error is first low-pass filtered and then acted on by a diagonal discrete-time PID law. For mode $j$ the implemented update can be written as

$$I_{j,k}=I_{j,k-1}+K_{i,j}\bar e_{j,k}, \qquad c_{j,k}=K_{p,j}\bar e_{j,k}+I_{j,k}+K_{d,j}(\bar e_{j,k}-\bar e_{j,k-1}),$$

where $\bar e_{j,k}$ is the filtered projected error; the stored discrete gains already contain the appropriate sampling-period factors. The proportional term reacts to the current error, the integral state accumulates persistent error and removes steady offsets, and the derivative term reacts to changes in the error and can add damping. The gains are obtained by tuning the diagonalized loop to the requested crossover frequency and phase margin. Although the Anand et al. paper presents a PI realization of the same generalized architecture, the repository implementation calls MATLAB's PID tuner and supports nonzero derivative gains, so “PID” is the accurate name for the controller simulated here.

The mapped correction $\Delta r_h=T_{\mathrm{post}}c$ is applied to the references of TCV's faster inner `hybrid` magnetic controller. The outer loop therefore does not directly command arbitrary coil voltages: `hybrid` continues to regulate plasma current and position and to generate actuator requests, while the shape controller biases its reference channels. This nesting preserves the established fast magnetic-control layer and separates it from the slower equilibrium-based shape loop. Low-pass filtering suppresses poorly modelled high-frequency dynamics, while reference blending and integral-state initialization provide bumpless activation. The controller code also supports saturation-aware anti-windup, although this particular call does not pass component-wise saturation bounds to `PID_step`.

## Purpose and outputs

The data collector is `/home/mg/repos/hybrid/pid_shape_liuqenn_comparison.m`, the figure script is `pid_shape_control_simulation.py`, and the generated figure is `dissertation_figures/pid_shape_control_simulation.pdf`. The collector performs matched closed-loop simulations using a neural magnetic observer and direct interpolation of the FGE flux map. It saves controller configuration and timing, point definitions and masks, observer signals, shape errors, projected PID coordinates, actuator signals, and FGE geometry and convergence diagnostics.

The figure contains six panels: all active shape-error trajectories, the active projected controller errors, the first-ten-channel RMS history, the post-activation RMS of each metric channel, the final FGE last closed flux surface (LCFS) and control points, and the norm of the coil-voltage vector. Solid curves denote the neural observer and dashed curves denote direct FGE interpolation with `qintmex`.

## Controller dimensions in this simulation

Eleven of the 25 allocated shape-error rows are active and the SVD retains four PID directions. The reported metric uses the first ten saved rows, matching the declared dissertation statistic, while the trajectory panel shows every active row identified by `active_shape_mask`.

## Simulation configuration

The simulation uses TCV shot 79742. FGE and hybrid execute with a $100\,\mu\mathrm{s}$ step, the shape controller executes every 2 ms, and shape control activates at 0.6 s. The simulated interval spans the 50 ms preceding activation through 0.8999 s. The FGE plant is linearized (`Lfge.P.lin = true`), the shape reference is fixed.

In the neural run, `shapectrl` evaluates `liuqenn_inference` using synthetic magnetic measurements. In the direct-observer run, `shapectrl` evaluates the simulated FGE flux map at the control points with `qintmex`. This path is an observer based directly on the FGE simulation state. Both runs use the same controller design, initial state, reference, and plant configuration.

## Saved data and interpretation

- `shape_errs` contains the regulated shape errors produced by the active observer. `shape_errs_direct` contains direct-FGE values evaluated within the same run, and the paired MAT files record the closed-loop trajectory obtained under each active observer.
- `proj_errs` contains `STP.Tpre * shapeErrs`, the errors presented to the diagonal PID channels. `active_projected_mask` identifies the four active rows.
- `Fx`, `Br`, and `Bz` contain observer values at the 25 allocated query positions; the `_direct` arrays contain direct FGE interpolation. Only positions identified by the saved shape masks and positive radial coordinates are physically meaningful. Unused neural-query rows can contain nonphysical values because the inference interface substitutes 1 mm for nonpositive radial coordinates.
- `shape_names_char`, `projected_names_char`, `sc_r`, `sc_z`, `sc_ref`, and the `sc_ls*` masks define the controlled quantities, their locations, and their roles. `Tpre`, `Tpost`, `SC2FB`, `Kp`, `Ki`, and `Kd` record the controller transformation and gains.
- `shape_hybrid`, `shape_current`, `shape_voltage`, and the corresponding `delta_*` arrays record the corrections produced by the outer loop. `coil_currents` and `coil_voltages` record realized plant actuation.
- `fge_final_flux` is the final flux map on `fge_grid_r` and `fge_grid_z`; the LCFS is the contour $\psi(R,Z)=\psi_B$ using `fge_final_boundary_flux`. `limiter_r` and `limiter_z` give the vessel limiter. `fge_boundary_r` and `fge_boundary_z` are the boundary-determining X-point or limiter point for each plasma domain.
- `fge_isconverged`, `fge_residual`, and `fge_iterations` record the numerical state of the FGE solution.

## Figure and quantitative result

The first panel shows all 11 active shape-error channels and reveals their transient and residual structure. The second shows the four projected coordinates acted on by the PID laws. The third reports the instantaneous RMS across the first ten saved error channels, and its legend gives the pooled post-activation RMS. The fourth separates that statistic by channel. The fifth plots the final LCFS as the boundary-flux contour together with the limiter and active control points. The sixth shows the Euclidean voltage-vector magnitude for both runs and the norm of their difference.

For $t\geq0.6$ s, the pooled RMS over the first ten shape-error rows is 2.75 mm with the neural observer and 3.60 mm with direct FGE interpolation. Both closed loops settle to millimetre-scale residuals in this deterministic linear scenario. The difference between these two values is a result for this initial condition, reference, and plant model.

## Sources

- Dissertation section `Controller-level validation in FGE` in `/home/mg/repos/phd_thesis/2_eqrec.tex`.
- H. Anand et al., *A Novel Plasma Position and Shape Controller for Advanced Configuration Development on the TCV Tokamak*, local file `/home/mg/repos/phd_thesis/material/tcv_ctrl_infra/P27_Anand2017.pdf`.
- `/home/mg/repos/phd_thesis/guides/guides.md`, `hybrid_guide.md`, `fge_guide.md`, and `meq_infra_report.md`.
- `/home/mg/repos/hybrid/pid_shape_liuqenn_comparison.m`, `ctrls/shape/shapectrl.m`, `ctrls/shape/shapectrl_step.m`, `utils/plot/shape_ctrl_plot.m`, and the associated controller-building code.
- `/home/mg/repos/meq/meqpdom.m` and `meqpost.m`, used to establish the meaning of `rB,zB`.
