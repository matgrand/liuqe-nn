# FGE — Forward Evolutive Grad-Shafranov Solver: Technical Reference

**Codebase**: MEQ (MatlabEQuilibrium) toolbox, Swiss Plasma Center EPFL  
**Module**: `fge.m` and supporting files  
**Purpose**: Time-stepping free-boundary MHD equilibrium solver for TCV and other tokamaks  
**Status**: Production code, Apache 2.0 license, 2022-2025

---

## 1. Overview

FGE (Forward Grad-Shafranov Evolution) is the time-stepping solver in the MEQ toolbox. It propagates a plasma free-boundary equilibrium forward in time by simultaneously evolving:

1. **Coil currents** — driven by applied voltage waveforms through circuit equations
2. **Vessel eddy currents** — passively induced by flux changes
3. **Plasma current distribution** — either prescribed via constraints or evolved through a Current Diffusion Equation (CDE)
4. **Poloidal flux map** — from the Grad-Shafranov equation at each time step

The physics picture: given an initial equilibrium state and coil voltage waveforms `Va(t)`, FGE computes the self-consistent evolution of the coupled electromagnetic-MHD system. This is distinct from static solvers (FGS/LIUQE) which find a single equilibrium given boundary conditions.

**Key relationship to other MEQ solvers:**
- FGS: static Grad-Shafranov solver — FGE uses FGS internally for initialization
- LIUQE/FBT: reconstruction solvers — FGE uses their output as initial conditions
- RZIP: a linearized 1D plasma model — FGE is the full nonlinear version

---

## 2. The Full Equation System

### 2.1 State Variables

The FGE state vector `x` (packed into a nonlinear solver vector `xnl`) has three components, depending on the `algoNL` parameter:

- **`algoNL = 'all-nl'`**: `x = [Iy(:); ag(:); Ia; Iu] ./ xscal`
- **`algoNL = 'all-nl-Fx'`**: `x = [Fx(:); ag(:); Ia; Iu] ./ xscal`
- **`algoNL = 'Newton-GS'`**: `x = [Fx(:); ag(:); Ia; Iu] ./ xscal` (different residual formulation)

Where:
- `Iy(nzy, nry)` — plasma toroidal current density on the inner grid [A/m²·m² = A per cell]
- `Fx(nzx, nrx)` — poloidal flux on the full computational grid [Wb]
- `ag(ng, 1)` — basis-function coefficients for plasma current profiles (p', TT')
- `Ia(na, 1)` — active coil currents [A]
- `Iu(nu, 1)` — vessel eddy current components [A]

The combined external current vector is `Ie = [Ia; Iu]`, dimension `ne = na + nu`.

`ag` are the coefficients of the Grad-Shafranov RHS: the plasma current profile is `j_phi = R·p'(psi_N) + TT'(psi_N)/(mu0·R)` expanded in basis functions parameterized by `ag`.

### 2.2 Electromagnetic Circuit Equations

The core circuit equation governing coil and vessel current evolution is:

$$M_{ee} \frac{dI_e}{dt} + R_e I_e + M_{ey} \frac{dI_y}{dt} = V_e$$

or more precisely as a block system:

$$\begin{bmatrix} M_{aa} & M_{au} \\ M_{au}^T & M_{uu} \end{bmatrix} \frac{d}{dt} \begin{bmatrix} I_a \\ I_u \end{bmatrix} + \begin{bmatrix} R_a & 0 \\ 0 & R_u \end{bmatrix} \begin{bmatrix} I_a \\ I_u \end{bmatrix} = \begin{bmatrix} V_a \\ 0 \end{bmatrix}$$

and the back-coupling from plasma:

$$M_{ey} \frac{dI_y}{dt} \text{ contributes to all rows}$$

In code:
```matlab
L.Mee = [G.Maa, G.Mau; G.Mau.', G.Muu]  % ne x ne inductance matrix
L.Re  = [G.Ra; G.Ru]                      % ne x 1 resistance vector
L.Mey = P.mey * [Mya.'; Myu.']            % ne x ny plasma-to-external mutual inductance
```

`P.mey` is a flag (0 or 1) that enables/disables back-coupling between plasma and vessel currents (used for testing). When `mey=0`, vessel eddy currents do not see the plasma's changing flux.

The **implicit Euler discretization** of this equation (as applied at each step in `fgeF.m`) is:

```
residualcirc = (M_ee*(Ie - Iep) + M_ey*(Iy - Iyp)) / dt  +  R_e .* Ie  -  Ve  =  0
```

where subscript `p` denotes the previous time step value. This is computed in `fgeF.m`:

```matlab
Fedot = (L.Mee*(Ie-[LYp.Ia;LYp.Iu]) + L.Mey*(Iy(:)-LYp.Iy(:)))*idt;
residualcirc = (Fedot + L.Re.*Ie - Ve);
```

### 2.3 The Grad-Shafranov Equation

At each time step, the equilibrium must satisfy the Grad-Shafranov (GS) equation:

$$\Delta^* \psi = -\mu_0 R j_\phi(\psi)$$

where $\Delta^*$ is the Grad-Shafranov operator (not the Laplacian — it includes the R-weighting from toroidal symmetry). In MEQ notation, the flux $F_x$ satisfies:

$$L.\text{dlst} \cdot F_x(:) = L.\text{rhsf} \cdot I_y(:) \quad \text{(on inner grid)}$$
$$F_x(\text{boundary}) = M_{be} I_e + M_{by} I_y \quad \text{(on outer grid)}$$

where:
- `L.dlst` is the finite-element discretization of `Delta*`
- `L.rhsf` is the RHS scaling factor (includes `mu0`, grid area element)
- `Mbe` is the Green's function matrix from external currents to boundary flux
- `Mby` is the Green's function matrix from interior plasma currents to boundary flux
- `L.Mxe = [Mxa, Mxu]` is the full grid Green's function from external currents

The GS equation is inverted to give `Iy` from `Fx` (or vice versa) via:

```matlab
Iy = -(L.dlst * Fx(:)) ./ L.rhsf - Iyie
```

where `Iyie = L.Tye * Ie` is the vacuum current density (the part of `Iy` due to external currents as seen by the GS operator).

### 2.4 The Constraint Equations (Profile Closure)

The `ag` coefficients are not free — they are determined by `nC` constraint equations of the form:

```
C(F0, F1, rA, ag, Fx, ...) = Co(t)
```

where `Co(t)` are prescribed target quantities (from LX, i.e., from measurements or reference trajectories). Typical constraints:
- `Ip` — total plasma current
- `bp` — poloidal beta
- `qA` — safety factor at the magnetic axis
- `li` — internal inductance
- `bt` — toroidal beta
- `ag` — direct prescription of basis function coefficients

This is the `agcon` parameter in `fgep.m`/`fgeptcv.m`. The number of constraints `nC` must equal the number of free basis function coefficients `ng` when no CDE is used (so that `ag` is fully determined).

### 2.5 The Current Diffusion Equation (CDE)

When a CDE is active (`P.cde` is not empty), `np` additional residual equations replace `np` of the profile constraints. The CDE replaces prescribing `ag` directly and instead evolves it through a physics model. Available CDE types (from `meqcdefun.m`):

| CDE Name | Description |
|---|---|
| `'OhmTor_rigid_0D'` | Rigid plasma Ohm's law: `Lp*(dIp/dt) + Rp*(Ip-Ini) + Fext_dot*(Iy/Ip) = 0` |
| `'OhmTor_0D'` | Toroidal Ohm's law: `int(psi_dot * Iy/Ip) * (1/dt) + Rp*(Ip-Ini) = 0` |
| `'cde_ss_0D'` | Steady-state 0D CDE from RAPTOR formulation (Carpanese thesis Eq. G.17) |
| `'cde_1D'` | Full 1D resistive diffusion in flux coordinates |
| `'cde_OhmTor_1D'` | Toroidal Ohm's law projected onto test functions: `int(psi_dot/R * g)dRdZ + int(2pi/sigma * (j_phi - j_ni) * g)dRdZ = 0` |
| `'static_cde_1D'` | Static 1D CDE constraint (for initialization) |
| `'static_Iy_1D'` | Static constraint matching Iy to reference under test functions |

**The `OhmTor_rigid_0D` model** (most commonly used) implements the projected circuit equation for the plasma:

$$L_p \frac{dI_p}{dt} + R_p (I_p - I_{ni}) + \langle \dot{\psi}_{ext} \rangle \frac{I_y}{I_p} = 0$$

where:
- $L_p = \frac{I_y^T M_{yy} I_y}{I_p^2}$ is the plasma self-inductance (computed or provided externally via `LX.Lp`)
- $R_p$ is the plasma resistance (provided as `LX.Rp`, Ohm)
- $I_{ni}$ is the non-inductive current drive (provided as `LX.Ini`, A)
- $\langle \dot{\psi}_{ext} \rangle$ is the loop voltage induced by all external sources (coils + other plasma domains)

**The `cde_OhmTor_1D` model** implements the full toroidal Ohm's law projected under test functions:

$$\int \frac{\partial \psi/\partial t}{R} g_k \, dR \, dZ + \int \frac{2\pi}{\sigma} (j_\phi - j_{ni}) g_k \, dR \, dZ = 0$$

for each test function $g_k$. This is a 1D PDE in flux-coordinate space (discretized via the basis function expansion) and includes the spatial resistivity profile $\sigma(\psi)$.

The resistivity profile enters via `LX.signeo` — a neoclassical conductivity profile on the `pQ` grid (normalized flux coordinate squared, `pQ = sqrt(psi_N)`). The profile must be strictly positive everywhere.

**Boundary conditions for 1D CDE**: The axis BC is implicit (regularity), and the edge/separatrix BC is contained in the plasma-domain definition (current must vanish at `F1 = FB`).

---

## 3. The fgeF Function — Core Physics Operator

`fgeF.m` implements the complete residual `F(x) = 0` and (optionally) its Jacobian. It is the heart of FGE.

### 3.1 Residual Structure

The full residual vector `res` has `nN` components, partitioned as:

```
res = [residualGS(:);    % nGS = ny or nx components (GS equation)
       residualag(:);    % nC components (profile constraints)
       residualcde(:);   % np components (CDE, optional)
       residualcirc(:)]  % ne components (circuit equation, FGE only)
```

The `nN = nGS + nC + np + ne` total size is set up in `meqalgoNLc`.

### 3.2 GS Residual Variants

Three residual formulations for the GS part:

**`all-nl`** (default): state is `Iy`, residual measures self-consistency:
```
residualGS = resscalGS .* (Iy1(:) - Iy(:))
```
where `Iy1 = Tyg * ag` is the new current from the profile model, `Iy` is the state.

**`all-nl-Fx`**: state is `Fx`, residual is flux self-consistency:
```
Fx1 = meqFx(L, Iy1, Ie, ...)
residualGS = resscalGS .* (Fx1(:) - Fx(:))
```

**`Newton-GS`**: mixed formulation — inner grid uses GS equation directly, boundary uses flux matching:
```
residualGS(inner) = rhsf .* (Iy1 - Iy)
residualGS(boundary) = Fx(bdy) - meqfbp(Iy1) - Fbe
```

### 3.3 Jacobians

`fgeF` computes four Jacobian matrices when requested:

- **`Jx = dF/dx`** at `dt=0` (static part: `K` matrix in the linearization)
- **`Ju = dF/du`** where `u = [Va; Co; Ini; Rp]` (`U` matrix)
- **`Jxdot = dF/dxdot`** at `xdot=0` (`S` matrix — mass/inductance matrix)

The pre-computed vacuum Jacobian (independent of plasma state) is stored in `L.Jx`, `L.Ju`, `L.Jxdot` (from `fgeFvacJac.m`). The full nonlinear Jacobian is obtained by adding plasma-state-dependent corrections.

The Jacobian action (matrix-vector product) is implemented separately in `fgeFJac.m` for efficiency with JFNK solvers.

---

## 4. Time Integration (fget.m)

### 4.1 Integration Scheme

FGE uses **backward Euler (implicit Euler)** time integration. At each time step from `t_p` to `t`, the implicit system is:

```
F(x_t, u_t, x_p) = 0
```

where `F` includes the discretized circuit equation with `1/dt` factors. This is first-order accurate in time (consistent with linearized convergence tests in `fge_convergence_tests.m` showing slope ≈ 1 in log-log).

The discrete-time linearized state-space form (in `fgess.m`) is:

```
A = (S - Ts*K) \ S
B = (S - Ts*K) \ (Ts*U)
```

confirming backward Euler discretization (`Ts = dt`).

### 4.2 Time Loop Structure (`fget.m`)

```
Initialize: LYp = LX(t=0) [from fgeinit]
For each time step t = t_2, ..., t_N:
  1. Get LXt = LX at current time
  2. Apply power supply model (voltage delays)
  3. If P.lin == true:
       Step linear model via fgetkl
  Else:
       Step nonlinear model via fgetk_implicit
  4. Compute derivatives (meqLYdot)
  5. Check protection limits (meqlim)
  6. Store LY, update LYp
  7. Call controller step (if fbcontrol)
```

### 4.3 Nonlinear Step: `fgetk_implicit.m`

At each time step, solves `F(x) = 0` where `F = fgeF` or `F = fgeFpicard` depending on `algoNL`:

```matlab
[xnlt, solverinfo] = solveF(F, xnl0, P);
```

The nonlinear solver is selected via `P.algoF`:
- `'newton'`: full Newton method with analytic or JFNK Jacobian
- `'jfnk'`: Jacobian-Free Newton-Krylov (GMRES inner solver)
- `'broyden'`: Broyden quasi-Newton

**Initial guess strategies** (parameter `initguess`, evaluated and the one with lowest `||F||^2` chosen):

| Strategy | Description |
|---|---|
| `'vacuum'` | Linear circuit step without plasma (vacuum approximation) |
| `'extrap'` | Linear extrapolation: `x_0 = x_p + xdot*dt` |
| `'previous'` | Use previous solution unchanged |
| `'LX'` | Use provided LX data |

**Convergence**: Solution accepted when `solverinfo.res < P.tolF` (default `1e-8`). If `keepok=true`, acceptable solutions with `res < tolok` (default `1e-2`) are accepted.

### 4.4 Linear Step: `fgetkl.m`

When `P.lin = true`, the nonlinear solver is replaced by stepping the linearized state-space model:

```matlab
dxt = lin.sysA * (zxt - lin.x0L) + lin.sysB * (ut - lin.u0L) + lin.x0dotL
yt  = lin.sysC * dxt + lin.sysD * (ut - lin.u0L) + lin.y0
```

The state is `zxt = [Ia; Iu; ag_i]` (reduced linear state), and inputs are `ut = [Va; Co; Ini; Rp; dCodt]`.

**Multiple linearizations**: FGE supports switching between multiple linearization points during a simulation run. The `idx_lin = ifloor([L.lin.t], LYp.t)` logic always selects the most recent linearization in the past. At transition points, `fgel_update_lin.m` ensures measurement continuity by computing an offset `dy0` and resetting the `ag` state to match plasma current (`fgel_update_ag.m`).

### 4.5 Plasma Domain Update

At each step, `meqpdom` is called within `fgeF` to identify:
- `FA`, `rA`, `zA` — axis flux value and location
- `FB`, `rB`, `zB` — separatrix (LCFS) flux value and location
- `FX`, `rX`, `zX` — X-point location(s)
- `F0`, `F1` — limiting flux values used for profile normalization
- `Opy(nzy, nry)` — integer domain map (0=vacuum, 1=first domain, 2=second for doublet)
- `lB`, `lX` — number of boundaries and X-points

This is fully dynamic — the plasma boundary moves with each Newton iteration.

### 4.6 Preconditioner (`fgepre.m`)

FGE computes an ILU or user-supplied preconditioner for the GMRES inner solver. The preconditioner is based on the linearized Jacobian `Jx + (1/dt)*Jxdot`, computed via `fgeFlin`. Update frequency is controlled by `P.precupdate` ('once', 'time', or 'newton').

The Picard-approximation preconditioner (`fgeFJac_inverse_picard_approx.m`) provides a block-triangular inverse approximation assuming `dIy/dFx ≈ 0` and `dCo/dFx ≈ 0`. It decomposes the Jacobian as:

```
inv(J) ≈ block-triangular involving:
   iJgg = inv(Jgg)            [ag block]
   iJee = inv(Jee + Jxdot/dt) [circuit block]
   JGSe, JeGS                  [coupling terms]
```

---

## 5. The CDE in Detail

### 5.1 The 0D Rigid Plasma Model (`OhmTor_rigid_0D`)

This is the most commonly used model for TCV simulations. The plasma is treated as a single inductor-resistor circuit element:

```
res = Lp * (IpD - IpD_prev)/dt  +  Rp * (IpD - Ini)  +  <psi_dot_ext, IyD/IpD>  =  0
```

Parameters entering from `LX`:
- `LX.Rp` — plasma resistance [Ohm]. Could be Spitzer or any model.
- `LX.Ini` — non-inductive driven current [A]
- `LX.Lp` — plasma self-inductance [H] (if `P.iLpext = true`, provided externally; if `false`, computed as `Iy^T * Myy * Iy / Ip^2`)

This model assumes the plasma current distribution shape (`Iy/Ip`) is frozen (rigid); only the amplitude `Ip` evolves. The profile coefficients `ag` that control the shape are thus constrained by the remaining `nC - 1` profile constraints (typically `bp` and `qA`), and the `Ip` degree of freedom is replaced by the CDE.

### 5.2 The Toroidal Ohm's Law (`OhmTor_0D`)

A projection of the full 2D Ohm's law onto the current profile:

```
res = (1/dt) * int[(psi_new - psi_old) * Iy/Ip dRdZ]  +  Rp * (Ip - Ini)  =  0
```

This accounts for the actual flux change at each grid point, weighted by the current distribution. More accurate than the rigid model for non-rigid current evolution.

### 5.3 The 1D Formulation (`cde_OhmTor_1D`, `cde_1D`)

The 1D CDE projecting toroidal Ohm's law under test functions `g_k`:

```
int [ (dpsi/dt) / R * g_k ] dR dZ  +  int [ (2*pi/sigma) * (j_phi - j_ni) * g_k ] dR dZ  =  0
```

The conductivity `sigma` is provided as a profile on the normalized flux grid via `LX.signeo` (dimension `nQ x nD`). In `fge_CDE_1D_analytical_test.m`, this is verified against an analytical screw-pinch solution using Bessel functions.

For the full 1D CDE (`cde_1D`), the resistive diffusion is in flux coordinates; for `cde_OhmTor_1D`, it is the toroidal Ohm's law in R,Z space projected onto test functions (valid in the small diamagnetic approximation).

### 5.4 Steady-State 0D from RAPTOR (`cde_ss_0D`)

Based on Carpanese's thesis (Eq. G.17, G.20, G.21, G.24), integrating the RAPTOR CDE from `rho=0` to `rho=1` assuming `dPhi/dt = 0`:

```
res = Ip - Ini  +  (rBt / (2*pi)) * int[sigma/T^2] dPhi * d(FA+FB)/2 / dt  =  0
```

This model uses the toroidal flux profiles and poloidal current function `T(psi)`.

### 5.5 CDE Initialization

In `fgeinit.m`, when a 0D CDE is used, it is replaced by `Ip` constraints for initialization (via `replace_cde_constraint`). For 1D CDEs, a `static_Iy_1D` or `static_cde_1D` model is used instead. This ensures a well-posed FGS problem for the initial condition.

---

## 6. Linearization (fgel, fgelinall, fgess)

This section covers the most important capability of FGE for control design.

### 6.1 Conceptual Framework (`fgel.m` header)

The FGE equations are split into:

- **Static equations** $e(\mathbf{x}, \mathbf{u}) = 0$: the GS equation + profile constraints
- **Dynamic equations** $f(\mathbf{x}) \frac{d}{dt}[\mathbf{g}(\mathbf{x})] + \mathbf{h}(\mathbf{x}, \mathbf{u}) = 0$: circuit + CDE

Linearized around an equilibrium $(\mathbf{x}^*, \mathbf{u}^*)$:

$$\frac{\partial e}{\partial x} d\mathbf{x} + \frac{\partial e}{\partial u} d\mathbf{u} = 0 \quad \text{(static)}$$
$$S \frac{d(d\mathbf{x})}{dt} + K \, d\mathbf{x} + U \, d\mathbf{u} = 0 \quad \text{(dynamic)}$$

where:
- $K = \frac{\partial F}{\partial x}$ evaluated at $dt=0$ (Jacobian of the static+dynamic parts)
- $S = \frac{\partial F}{\partial \dot{\mathbf{x}}}$ (the mass/inductance matrix)
- $U = \frac{\partial F}{\partial u}$ at $dt=0$

The continuous-time state-space system is:
$$S \frac{d(d\mathbf{x})}{dt} = K \, d\mathbf{x} + U \, d\mathbf{u} \quad \Rightarrow \quad \frac{d(d\mathbf{x})}{dt} = A \, d\mathbf{x} + B \, d\mathbf{u}$$
$$\mathbf{y} = C \, d\mathbf{x} + D \, d\mathbf{u} + \mathbf{y}_0$$

where $A = S^{-1}K$, $B = S^{-1}U$.

### 6.2 State Elimination (`fgselim`)

The full state `x = [Iy/Fx; ag; Ia; Iu]` is too large for a practical LTI system. `fgselim` (called inside `fgel`) performs **static condensation**: the static equations `de/dx * dx = -de/du * du` are used to eliminate `Iy`/`Fx` (GS equation) and `ag` (constraint equations), leaving only the **dynamic states**:

- `Ia` — active coil currents (na components)
- `Iu` — vessel eddy current components (nu components)
- `ag_i` — the CDE-controlled basis function coefficients (np components, only when CDE is active)

The reduced state `xL = [Ia; Iu; ag_i]` lives in the linear subspace identified by `lin.maskL`. The mapping from full state to reduced state is `lin.dxdxL = dx/dxL`.

### 6.3 State-Space Matrices (`fgess.m`)

The standard Walker-et-al. (FST 2006) state-space form:

```
S * dxL/dt = K * xL + U * u
y = C * xL + D * u + yo
```

**Inputs** `u = [Va; Co; Ini; Rp; dCodt]`:
- `Va (na)` — active coil voltages
- `Co (nC)` — constraint quantities (Ip, bp, qA, etc.)
- `Ini (nD)` — non-inductive current drive (CDE only)
- `Rp (nD)` — plasma resistance (CDE only)
- `dCodt (nC)` — time derivatives of constraints

**State** `xL = [Ia; Iu; (ag_i if CDE)]`:

**Outputs** `y` (default 'all' list):
```matlab
measlist = {'Bm','Ff','Ia','Iu','Ip','Ft','rIp','zIp','Iy','Fx','FA','FB',
            'ag','rA','zA','bp','qA','li','Wk'}
```

Where:
- `Bm` — magnetic probe field measurements
- `Ff` — flux loop measurements
- `Ia`, `Iu` — active and vessel currents
- `rA`, `zA` — magnetic axis position
- `FA`, `FB` — axis and boundary flux values
- `bp`, `qA`, `li`, `Wk` — integral plasma shape parameters

**Output matrices** `C` and `D` are assembled from the linearizations of each quantity with respect to `xL` and `u` (computed in `fgel.m` by calling `meqdFdI` plus algebraic chain rules):

```matlab
lin.dBmdxL = dBmdx * dxdxL  % B-probe response to state change
lin.dFfdxL = dFfdx * dxdxL  % Flux-loop response to state change
lin.dIadxL = dIadx * dxdxL  % etc.
...
lin.drAdxL = drAdx * dxdxL  % Axis position response
```

The `yo` offset vector ensures the linearization is exact at the operating point: `yo = y0 - C*x0L - D*u0`.

### 6.4 Discrete-Time Realization

For time-step $T_s$, backward Euler gives:

$$A = (S - T_s K)^{-1} S$$
$$B = (S - T_s K)^{-1} (T_s U)$$

with $\dot{\mathbf{x}}_0 = -(S - T_s K)^{-1} (T_s \, \text{res}_{0L})$ accounting for non-zero residual at the initial point (non-stationary initial condition).

### 6.5 Multiple Linearizations (`fgelinall.m`)

For long simulations where the operating point drifts significantly, FGE supports an array `L.lin(1..N)` of linearization points at different times. During a linearized simulation (`P.lin=true`), `fget` selects `idx_lin = max(ifloor([L.lin.t], LYp.t), 1)` — always the most recent past linearization.

At a linearization switch, `fgel_update_lin.m` handles the transition:
1. Resets `ag_i` to match `Ip` continuity (via `fgel_update_ag`)
2. Computes measurement offset `dy0` to prevent discontinuities in controller inputs

### 6.6 Eigenvalue Analysis (`fgeeig.m`)

`fgeeig` extracts unstable eigenvalues of the linearized system. Key feature: it can operate directly from the Jacobians without building the full state-space (for speed). The most unstable mode (largest real part of eigenvalue) typically corresponds to the vertical instability (VDE).

---

## 7. Parameters (fgep.m)

Complete parameter catalog for FGE-specific settings:

### 7.1 FGE-Specific Parameters

| Parameter | Default | Description |
|---|---|---|
| `lin` | `false` | Run linearized simulation instead of NL |
| `mey` | `1` | Scale Mey (0=no plasma-vessel back-coupling, for testing) |
| `iLpext` | `true` | Lp source: true=from LX, false=computed as `Iy^T*Myy*Iy/Ip^2` |
| `ctrlfct` | `[]` | Feedback controller function handle (open-loop if empty) |
| `diagfct` | `[]` | Synthetic diagnostic model function handle |
| `psfct` | `@meqps` | Power supply model (handles voltage delays) |

### 7.2 Simulation Flow Parameters

| Parameter | Default | Description |
|---|---|---|
| `stoponalarm` | `true` | Stop if coil protection limits exceeded |
| `nnoc` | `0` | Max consecutive non-converging steps allowed |
| `keepok` | `false` | Accept solution if `res < tolok` even if not converged |
| `tolok` | `1e-2` | "Good enough" tolerance when `keepok=true` |
| `ssinit` | `true` | Initialize to steady-state (find consistent Ia, Va) |
| `wcompa` | `1e-9` | Regularization weight for compensation coil currents |
| `wioha` | `5e-11` | Regularization weight for OH current-driving currents |
| `compregion` | `'Oasx-sym'` | Region for compensation ('Iy', 'Oasx', 'Oasx-sym') |
| `initguess` | `{'vacuum'}` | Initial guess strategies list |

### 7.3 Inherited Solver Parameters (from fgsp.m)

| Parameter | Default | Description |
|---|---|---|
| `algoNL` | `'all-nl'` | NL formulation ('all-nl', 'all-nl-Fx', 'Newton-GS', 'Picard') |
| `algoF` | `'newton'` | Root finder ('newton', 'jfnk', 'broyden') |
| `tolF` | `1e-8` | Nonlinear residual tolerance (FGE default; FGS default is `1e-10`) |
| `kmax` | `50` | Max Newton iterations |
| `mkryl` | `200` | GMRES Krylov subspace dimension (FGE uses 200; FGS uses 150) |
| `algoGMRES` | `'aya'` | GMRES variant ('sim', 'giv', 'aya') |
| `usepreconditioner` | `true` | Enable ILU preconditioner |
| `precupdate` | varies | Preconditioner update frequency |
| `anajac` | `true` | Use analytical Jacobian (false = finite differences) |
| `xscal` | `true` | Scale state vector to O(1) |
| `cde` | `''` | CDE type string (empty = no CDE) |
| `Ipmin` | small | Below this Ip, assign Iy=LX.Iy (vacuum/zero-current mode) |

### 7.4 TCV-Specific Parameters (`fgeptcv.m`)

| Parameter | Default | Description |
|---|---|---|
| `insrc` | `'liu'` | Source for LX data: 'fbt' or 'liu' (LIUQE.M) |
| `selu` | `'v'` | Vessel current representation: 'v'=full, 'e'=eigenmode |
| `izgrid` | `true` | Include z-grid for flux interpolation |
| `eqfct` | `@meqetcv` | Equilibrium function |

TCV power supply delays (`fgegtcv.m`):
- E, F, OH coils: 0.5 ms delay
- G coils: 0.05 ms delay

---

## 8. Geometry (fgeg.m, fgegtcv.m)

`fgeg.m` calls `meqg` to assemble geometry and electromagnetic parameters. The key precomputed quantities:

### 8.1 Computational Grids

- `(rx, zx)`: full R,Z grid for flux (`nzx × nrx`)
- `(ry, zy)`: inner grid for plasma current (`nzy × nry = nzx-2 × nrx-2`)
- `lxy`: logical mask mapping inner grid indices to full grid indices
- `L.dlst`: sparse `Delta*` operator matrix
- `L.rhsf`: RHS factor vector (includes `mu0`, cell area)

### 8.2 Green's Function Matrices

- `G.Mxa (nx × na)`: flux on full grid from active coil unit currents
- `G.Mxu (nx × nu)`: flux on full grid from vessel unit currents
- `G.Maa (na × na)`: mutual inductance between active coils
- `G.Mau (na × nu)`: mutual inductance between active and vessel
- `G.Muu (nu × nu)`: mutual inductance between vessel segments
- `G.Ra (na × 1)`: active coil resistances [Ohm]
- `G.Ru (nu × 1)`: vessel segment resistances [Ohm]
- `G.Bmx (nm × nx)`: B-probe response to plasma current density
- `G.Bma (nm × na)`: B-probe response to active coil currents
- `L.Mbe (nb × ne)`: boundary flux from external currents (for GS BCs)
- `L.Mey (ne × ny)`: external current response to plasma current changes

The vessel can be represented as either the full set of segments (`selu='v'`, nu=up to 256 for TCV) or as a truncated eigenmode basis (`selu='e'`, nu selected by `P.nu`). The eigenmode representation uses the nu dominant modes of the vessel eddy current response, ordered by time constant. Using `selu='e'` with `nu=20-50` gives a good approximation at much reduced cost.

### 8.3 Relationship to LIUQE/FBT

FGE uses the same Green's function matrices as LIUQE and FBT (via `meqgtcv`). The equilibrium coordinate transforms (`meqetcv`) are shared. This ensures that an FGE simulation initialized from a LIUQE reconstruction gives consistent results.

---

## 9. Static Variant (fgess.m)

`fgess.m` has a dual role:
1. **As a state-space constructor**: called after `fgel` to produce MATLAB `ss` objects
2. The name is confusingly similar to `fgst.m` (the static FGS time-stepper) — these are different functions

The static solver used for FGE initialization is `fgst` (Forward Grad-Shafranov Time-stepper, no circuit dynamics), called inside `fgeinit.m`. The initial equilibrium is found by `fgst` with static parameters: the circuit equations are replaced by the current `Ia`, `Iu` as inputs (not states).

**Difference from FBT/LIUQE reconstruction**: FGS finds the equilibrium consistent with given external currents and profile constraints; FBT finds it consistent with magnetic measurements. FGE uses FGS-type initialization (constraints + currents).

---

## 10. Initialization (`fgeinit.m`, `fgex.m`)

### 10.1 fgex — Getting Initial and Time-Varying Data

`fgex` builds the `LX` structure (time-varying inputs) from an FGS/LIUQE solution:

1. Calls `fgsxtcv` (or equivalent) to get initial equilibrium `LX0` from the database
2. Calls `L.codeinit` = `fgeinit` to compute steady-state consistent initial conditions
3. Assembles time arrays for: constraint quantities (`Ip`, `bp`, etc.), `signeo`, `Lp`, `Rp`, `Ini`
4. Builds feedforward `Va` trajectory: `LX.Va = Va0 + Vadot0 * (t - t0)` to maintain the initial equilibrium over time

Key fields assembled in `LX`:
- `Va (na × nt)` — coil voltage waveforms
- `Ia (na × nt)` — reference active coil currents (for `initguess='LX'`)
- `signeo (nQ × nD × nt)` — conductivity profile (must be strictly positive)
- `Lp (nD × nt)` — plasma inductance (if `iLpext=true`)
- `Rp (nD × nt)` — plasma resistance
- `Ini (nD × nt)` — non-inductive current drive

### 10.2 fgeinit — Steady-State Initialization

`fgeinit` finds `Ia, Iu, Va` consistent with the initial equilibrium by:

1. Calling `meqinit` to find compensation coil currents that produce zero net flux change in the plasma region
2. Running FGS (`fgst`) to get a converged equilibrium starting from the LX data
3. Computing `Vadot0` and `Iadot0` — the rate of change needed to maintain the equilibrium (for a non-zero `Ip` ramp, coil currents must change)

The `fgeIaind.m` function computes the OH coil direction that drives loop voltage in the plasma without changing the field shape (used for `Ip` control at initialization). The `fgeIacomp.m` function finds compensation coil currents that null the field in a specified region.

---

## 11. Physics Accuracy and Tradeoffs

### 11.1 Approximations

FGE makes the following approximations relative to a full MHD code:

1. **Axisymmetry**: Full 3D effects (NTMs, RWMs, 3D field errors) not included
2. **Force-free equilibrium**: The GS equation assumes `j × B = ∇p`, i.e., inertia and flow effects are neglected
3. **Resistive diffusion**: Plasma evolution is governed by resistive diffusion rather than ideal MHD. Valid on resistive timescales (`tau_R >> tau_A`)
4. **Profile parametrization**: `j_phi` is represented by a finite number of basis functions (`ng` terms). The `ag` vector captures only the low-order content of the current profile
5. **Vessel linearization**: Eddy currents are modeled as linear circuits. Nonlinear ferromagnetic effects in the vessel walls are neglected

### 11.2 Validity Range

The resistive diffusion model is valid when:
- `dt >> tau_A` (Alfvén time, ~microseconds) — inertia negligible
- `dt <= tau_R` (resistive diffusion time, ~seconds for plasmas) — resistivity matters
- Plasma remains in quasi-static force balance between time steps

Typical FGE time steps: `dt = 50 µs - 1 ms` for TCV (well within the validity range).

### 11.3 Numerical Stability

- **Backward Euler** provides unconditional stability for the circuit equations (unlike forward Euler which would require `dt < tau_R_min`)
- **Vertical instability** (VDE): the most unstable mode of `A = S\K` has growth rate ~1000/s on TCV. FGE resolves this even without a controller, but `dt < 1/growth_rate ≈ 1 ms` is needed for accuracy
- **Non-convergence**: At large `dt` or large perturbations, Newton convergence may fail. The `nnoc` parameter controls how many consecutive non-converging steps to tolerate
- **Preconditioner**: The Picard-block preconditioner significantly reduces GMRES iterations. Without it (or with stale preconditioner), convergence degrades

---

## 12. Issues for Python Rewrite

### 12.1 Time Integration

The backward Euler scheme is:
```python
F(x_new, x_old, dt) = 0
```

**Recommendation**: Do NOT use `scipy.integrate.solve_ivp`. FGE is a DAE (differential-algebraic equation), not an ODE. The GS equation is an algebraic constraint at each time step. `solve_ivp` doesn't handle DAEs.

Instead, implement explicit backward Euler with a Newton solver per step:

```python
def fge_step(x_old, dt, LXt):
    def residual(x_new):
        return fgeF(x_new, x_old, dt, LXt)
    return scipy.optimize.root(residual, x_old, method='krylov')
```

For the linearized model, the discrete ABCD matrices can be directly applied without a solver.

### 12.2 Handling the GS Constraint

The GS equation is a large sparse linear system. In Python:
- Store `dlst` (the Delta* matrix) as `scipy.sparse.csr_matrix`
- Solve `dlst * Fx = rhsf * Iy` using `scipy.sparse.linalg.spsolve`
- The Green's function matrices `Mxa`, `Mxu` are dense for small grids but can be large (e.g., `1024×128` on TCV)

The GS solve is the computational bottleneck. In MATLAB MEQ this is done via the precomputed `dlst` matrix and dense/sparse matrix-vector products. Python `numpy`/`scipy` equivalents will work but may need careful memory management.

### 12.3 Linearization

The state-space matrices `A, B, C, D` can be extracted via `fgess` and saved as numpy arrays. For control design, use `scipy.signal.StateSpace` or `control.StateSpace` (python-control package):

```python
import scipy.signal as sig
sys = sig.StateSpace(A, B, C, D)
```

The step `S\K`, `S\U` matrix inversions in `fgess` should be done with `numpy.linalg.solve` or `scipy.linalg.solve` (not `inv`).

### 12.4 GMRES Solver

Replace MATLAB's GMRES with `scipy.sparse.linalg.gmres`. The Picard preconditioner involves:
- A block inverse of `Jgg` (ng×ng, small)
- A block inverse of `Jee` (ne×ne, up to ~256×256)
- Matrix-vector products

Implement the preconditioner as a `scipy.sparse.linalg.LinearOperator`.

### 12.5 Performance Bottleneck

The main cost per time step is:
1. `meqFx`: flux computation `Fx = Mxa*Ia + Mxu*Iu + Mxy*Iy` — dense matrix-vector products
2. `meqpdom`: plasma domain identification — involves grid search operations
3. `fgeFJac`: Jacobian assembly — complex chain-rule accumulation
4. GMRES: typically 5-20 iterations × Krylov dimension

On TCV grids (~64×128), step (1) takes ~1ms/call in MATLAB. In Python with NumPy, expect similar or slightly faster (BLAS-optimized). The bottleneck can be shifted to Python by vectorizing `meqpdom`.

### 12.6 Data Structures

The `L` structure maps to a Python dataclass or dict-of-arrays. Key sub-structures:
- `L.G` — geometry (all static, computed once)
- `L.P` — parameters
- `L.ind` — index arrays into the state vector
- `L.lin` — linearization data (after calling `fgel`)

The `LX`/`LY` time-series structures map to dicts of 2D numpy arrays `(quantity, time)`.

### 12.7 CDE Implementation

The CDE residuals (`meqcdefun`) are relatively straightforward to translate. The `signeo` profile interpolation:

```python
sigma = np.interp(FNy.ravel(), pQ**2, signeo_profile).reshape(FNy.shape)
```

The `OhmTor_rigid_0D` model needs `Lp` either provided or computed as `np.dot(IyD.ravel(), FyD) / IpD**2`.

---

## 13. File Index

| File | Role |
|---|---|
| `/home/grandin/repos/meq/fge.m` | Top-level entry point |
| `/home/grandin/repos/meq/fgep.m` | FGE parameter defaults |
| `/home/grandin/repos/meq/fgeg.m` | Geometry assembly |
| `/home/grandin/repos/meq/fgec.m` | Parameter consolidation, CDE setup, Jacobian pre-computation |
| `/home/grandin/repos/meq/fgex.m` | Initial and time-varying data preparation |
| `/home/grandin/repos/meq/fgeptcv.m` | TCV-specific parameters |
| `/home/grandin/repos/meq/fgegtcv.m` | TCV-specific geometry (power supply delays) |
| `/home/grandin/repos/meq/fgextcv.m` | TCV-specific data retrieval (wraps fgsxtcv) |
| `/home/grandin/repos/meq/fget.m` | Main time-stepping loop |
| `/home/grandin/repos/meq/fgetk_implicit.m` | Single implicit Euler time step |
| `/home/grandin/repos/meq/fgetkl.m` | Single linear model time step |
| `/home/grandin/repos/meq/fgeinit.m` | Initial condition computation |
| `/home/grandin/repos/meq/fgepre.m` | Preconditioner assembly |
| `/home/grandin/repos/meq/fgess.m` | State-space model constructor (also discrete-time) |
| `/home/grandin/repos/meq/fgeF.m` | Core physics residual F(x)=0 |
| `/home/grandin/repos/meq/fgeFJac.m` | Jacobian action (matrix-free J*v) |
| `/home/grandin/repos/meq/fgeFlin.m` | Jacobian evaluation at equilibrium |
| `/home/grandin/repos/meq/fgeFpicard.m` | Picard iteration operator |
| `/home/grandin/repos/meq/fgeFJac_inverse_picard_approx.m` | Block-triangular preconditioner |
| `/home/grandin/repos/meq/fgeFvacJac.m` | Vacuum (state-independent) Jacobian |
| `/home/grandin/repos/meq/fgeIacomp.m` | Compensation coil current computation |
| `/home/grandin/repos/meq/fgeIaind.m` | OH coil current direction for Ip drive |
| `/home/grandin/repos/meq/fgelinall.m` | Batch linearization over multiple time slices |
| `/home/grandin/repos/meq/fgel.m` | Linearization at a single equilibrium |
| `/home/grandin/repos/meq/fgel_update_ag.m` | ag reset at linearization transition |
| `/home/grandin/repos/meq/fgel_update_lin.m` | Full linearization transition handling |
| `/home/grandin/repos/meq/fgeeig.m` | Eigenvalue analysis of linearized system |
| `/home/grandin/repos/meq/meqcdefun.m` | CDE residual functions catalog |
| `/home/grandin/repos/meq/fgsp.m` | FGS/FGE base solver parameters |
| `/home/grandin/repos/meq/tests/fge_test.m` | Basic FGE tests |
| `/home/grandin/repos/meq/tests/fge_convergence_tests.m` | Time convergence + algoNL comparison |
| `/home/grandin/repos/meq/tests/fge_CDE_test.m` | CDE model tests |
| `/home/grandin/repos/meq/tests/fge_CDE_1D_analytical_test.m` | 1D CDE vs analytical screw-pinch |
| `/home/grandin/repos/meq/tests/fgelin_test.m` | Linear vs nonlinear simulation comparison |
| `/home/grandin/repos/meq/tests/fge_constraint_tests.m` | Profile constraint satisfaction tests |
| `/home/grandin/repos/meq/tests/fgetcv_test.m` | TCV-specific data loading test |
