# FMag X-point control

## Purpose

The data collector is `/home/mg/repos/hybrid/fmag_liuqenn_comparison.m`, the figure script is `fmag_control_simulation.py`, and the generated figure is `dissertation_figures/fmag_control_simulation.pdf`. The experiment compares the closed-loop capability of the same FMag controller and nonlinear FGE plant when the local X-point quantities are supplied either by LiuqeNN or by the standard LIH linear observer. The comparison is deliberately restricted to the observer outputs that influence this particular controller.

## Signals and closed-loop structure

FGE is the simulated plant. At each time step it evolves the plasma and coil state, then produces simulated magnetic measurements $M$. FMag does not receive the FGE equilibrium geometry directly. Its observer maps $M$ to estimated controlled quantities $U_C$, FMag compares $U_C$ with its references, and the resulting control law produces coil-voltage commands. FGE then applies those commands and produces the next plant state and measurements. The FGE equilibrium quantities $Y_C$, including the actual X-point coordinates, are saved only as external truth for evaluating the closed loop.

The LiuqeNN measurement input is ordered as magnetic probes, flux loops, active-coil currents, flux-loop voltages, and toroidal-field signal: `Bm,Ff,Ia(E,F,G,OH),Uf,rBt`. FMag stores its active-coil currents in `E-F-OH-G` order, whereas LiuqeNN was trained with `E-F-G-OH`. The neural path must therefore apply the permutation `[FP.iE; FP.iF; FP.iG; FP.iOH]` before inference. Feeding the unpermuted FMag current vector to LiuqeNN assigns the OH currents to the G input and the G current to an OH input, invalidating every neural flux derivative and the resulting control comparison. The collector records both the permutation and reordered coil names, and the plotting script refuses data that do not contain the corrected ordering.

The LIH and LiuqeNN simulations begin from the same plant initialization, controller matrices, references, feedforwards, and timing. Once their observer estimates cause different voltage commands, their FGE trajectories and later measurement vectors are no longer identical; this is the intended closed-loop consequence of changing the observer. Within each trajectory, both observer outputs are recorded from the same measurements while only the selected observer enters feedback, allowing the other observer to act as a passive same-$M$ shadow.

## FMag X-point observable

For the fixed design query point $x_0=(R_0,Z_0)$, the observer supplies the local flux Taylor model

$$
\psi(x_0+\delta x)\simeq\psi_0+g^{\mathsf T}\delta x+\frac{1}{2}\delta x^{\mathsf T}H\delta x,
$$

where $g=(\partial_R\psi,\partial_Z\psi)$ and $H$ contains $\partial_{RR}\psi$, $\partial_{RZ}\psi$, and $\partial_{ZZ}\psi$. Setting the model gradient to zero gives the estimated displacement from the design point,

$$
\delta x=-H^{-1}g.
$$

Consequently, LiuqeNN is queried at exactly one $(R,Z)$ position per controller step and returns six quantities: $\psi$, the two first derivatives, and the three independent second derivatives. The derivatives are all required to calculate the radial and vertical displacement estimates $\delta R_X$ and $\delta Z_X$.

FMag subtracts its observable vector from the reference vector, scales the geometric displacement channels by the plasma-current reference, and mixes the errors through `Tpre`. Diagonal PID state-space blocks act on the mixed coordinates, `Tpost` distributes their outputs into active-coil current-derivative requests, and the coil model converts those requests into voltages.

The selected configuration has one active X-point. Its two displacement outputs have nonzero `Tpre` weights and therefore affect the voltage command. Its X-point-flux output has zero `Tpre` weight, and the standard controller builder creates no boundary-flux observer rows. Boundary flux and boundary-subtracted X-point flux therefore cannot affect the actuator command in this experiment. No LCFS or boundary collocation points are generated or passed to LiuqeNN. The collector asserts these structural facts before starting either simulation.

The remaining feedback paths are unchanged. Active-coil currents are direct measurements, while $I_p$, $rI_p$, $zI_p$, and the fast vertical observable retain their existing LIH-based paths. Current density is not consumed by this FMag controller. Controller initialization also remains on the established LIH path; LiuqeNN replaces the X-point Taylor quantities during normal controller steps only.

## Data collection

The default collector configuration uses TCV shot 78071 from 0.6 s to 1.2 s with a $100\,\mu\mathrm{s}$ FGE and controller step; shot 79742 is retained as the commented alternative. FGE uses the nonlinear `OhmTor_0D` current-density evolution, 100 passive-structure modes, and the E, F, OH, and G coil groups. Artificial diagnostic noise is disabled.

Schema version 5 records the observer identity, the single query coordinate, corrected neural coil permutation and names, perceived X-point displacement errors, simultaneous LIH and LiuqeNN shadow estimates, local Hessian entries, controller matrices and outputs, coil currents and voltages, FGE X-point truth, and convergence diagnostics. Both observers are evaluated at every normal controller step, but only the selected active observer supplies feedback in each run.

The stored `xpt_delta_r` and `xpt_delta_z` are the active observer estimates. The corresponding `_lih` and `_nn` arrays contain the two observer outputs evaluated from the same measurements in either file, making both same-trajectory shadows available. The stored `xpt_err_delta_r` and `xpt_err_delta_z` are FMag's internally scaled feedback errors; division by `ip_scaling_reference` recovers displacement error in metres. The FGE truth displacement is obtained from `fge_active_xpoint_r-design_xpoint_r` and `fge_active_xpoint_z-design_xpoint_z`.

## Figure interpretation

The top-left panel shows the signed radial and vertical control errors perceived by the two active observers, while the top-right panel shows their Euclidean norm. These are reference-minus-estimate errors formed from the controlled quantity supplied to FMag, and therefore are the signals that the controller directly drives toward zero. In both runs the transient decays and the perceived error converges without subsequent drift. This demonstrates that both observer-controller combinations close a stable loop: LIH and LiuqeNN each provide a controlled quantity that FMag can regulate. It does not yet establish where the physical X-point converges, because the controller can make its observer estimate coincide with the reference even when that estimate is biased relative to FGE.

The central and bottom-left panels answer the separate question of where the system actually converges. Each central panel contains the radial and vertical displacements estimated by LIH and LiuqeNN together with the corresponding FGE ground truth. The center-left panel uses the LiuqeNN-controlled trajectory and the center-right panel uses the LIH-controlled trajectory; the inactive observer in each panel is a passive shadow evaluated from the same measurement history and cannot affect feedback. Solid lines denote radial displacement and dashed lines vertical displacement. Both panels use the same vertical scale so their accuracy can be compared directly. The traces settle rather than drift, but their steady values need not agree: the LiuqeNN estimate remains substantially closer to the FGE X-point than the LIH estimate. Thus both observers stabilize their perceived controlled quantity, while LiuqeNN is considerably more precise with respect to the simulated physical ground truth.

The bottom-left panel reduces the FGE radial and vertical offsets to the physical distance from the target. It is the direct closed-loop capability metric: it confirms that the absence of drift seen in the central panels corresponds to convergence, and that the LiuqeNN-controlled plant settles closer to the intended X-point than the LIH-controlled plant. The bottom-right panel is different in purpose: it compares the coil-voltage vector norms and therefore shows the actuator consequence of substituting the observer rather than the achieved X-point position.

The color meaning is invariant across the six panels: C1 denotes LIH, C2 denotes LiuqeNN, and C3 denotes the FGE physical reference. FGE appears explicitly only where plant truth and observer estimates are shown on the same trajectory; in the physical-error panel, each curve is already calculated from the FGE X-point of its corresponding closed-loop run.

## Scope

This is a deterministic single-shot, single-X-point test with fixed target geometry. It tests the causal effect of replacing the X-point displacement observer in the controller while retaining the existing plant and other observer paths. It does not test boundary-flux feedback, multiple X-points, different plasma topologies, measurement noise, reference changes, profile perturbations, or broad operating-domain robustness.

## Sources

- `/home/mg/repos/phd_thesis/guides/guides.md`, `hybrid_guide.md`, `fge_guide.md`, `liuqe_guide.md`, and `meq_infra_report.md`.
- `/home/mg/repos/hybrid/fmag_liuqenn_comparison.m`, `ctrls/SCDalgofmag/FMag_control.m`, `ctrls/SCDalgofmag/FMag_step.m`, `ctrls/SCDalgofmag/FMag_init.m`, `get_fmag_xpoint_observer.m`, and `Examples/build_RZIp_Xpoint_controller.m`.
