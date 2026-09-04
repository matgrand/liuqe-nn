# LIUQE Technical Reference

Comprehensive technical reference for the LIUQE equilibrium reconstruction algorithm
as implemented in the MEQ (MatlabEQuilibrium) toolbox, Swiss Plasma Center EPFL.

Primary reference: Moret et al., Fusion Engineering and Design 91 (2015) 1–15.
Original FORTRAN LIUQE: Hofmann and Tonetti, Nuclear Fusion 28(10):1871, 1988.

---

## 1. Overview of Equilibrium Reconstruction

### 1.1 The Physical Problem

A tokamak plasma in MHD equilibrium satisfies the Grad-Shafranov (GS) equation, which
determines the poloidal magnetic flux ψ(R,Z) everywhere in the machine cross-section.
The equilibrium reconstruction problem is: given a set of magnetic diagnostic measurements
**y** (flux loop readings, magnetic probe readings, plasma current, etc.), find ψ(R,Z) and
the plasma current profiles p'(ψ), TT'(ψ) such that the forward model F(ψ, j_φ) matches
the observed measurements within their noise.

The measurement model is:

```
y = F(ψ, j_φ) + noise
```

where F is the operator that maps equilibrium fields to diagnostic signals (e.g., flux
loop readings are line integrals of ψ, magnetic probe readings are field components, etc.).

### 1.2 The Grad-Shafranov Equation

The GS equation in standard form is:

$$\Delta^* \psi = -\mu_0 R p'(\psi) - \frac{TT'(\psi)}{R}$$

where:
- $\psi(R,Z)$ is the poloidal magnetic flux divided by $2\pi$ (Wb/rad), defined on a 2D $(R,Z)$ grid in the machine cross-section.
- $\Delta^*$ is the Grad-Shafranov operator: $\Delta^* \psi = R \frac{\partial}{\partial R}\left(\frac{1}{R} \frac{\partial \psi}{\partial R}\right) + \frac{\partial^2 \psi}{\partial Z^2}$. In the code this is assembled as a sparse matrix `L.dlst` by `dlstc.m`.
- $p'(\psi) = dp/d\psi$ is the pressure gradient with respect to poloidal flux.
- $TT'(\psi) = T \, dT/d\psi$ where $T(\psi) = R \cdot B_\phi$ is the toroidal field function. Thus $TT'/R$ is the diamagnetic current contribution.
- $\mu_0 = 4\pi \times 10^{-7}$ H/m is the permeability of free space.
- $j_\phi = \mu_0 R p'(\psi) + TT'(\psi)/R$ is the toroidal plasma current density (right-hand side).

The code represents $p'(\psi)$ and $TT'(\psi)$ as a linear combination of basis functions:

$$p'(\psi) = \sum_g a_g \cdot g'_P(\psi)$$
$$TT'(\psi) = \sum_g a_g \cdot g'_T(\psi)$$

where $g'_P$ and $g'_T$ are the basis functions (see Section 3) and $a_g$ are the coefficients to be optimised. This makes $j_\phi$ linear in the coefficients $a_g$, and the GS equation becomes linear in $a_g$ for a fixed flux map $\psi$.

### 1.3 Nonlinearity and Iteration

The reconstruction is fundamentally nonlinear because:

1. **Boundary location depends on ψ**: The Last Closed Flux Surface (LCFS / separatrix /
   plasma boundary) is the flux contour that just touches a limiter or forms an X-point.
   Which grid cells are "inside the plasma" depends on ψ, but ψ depends on where the
   plasma current flows, which depends on the boundary. This is the main outer nonlinearity.

2. **Basis function evaluation depends on ψ normalization**: The profile functions
   p'(ψ_N) and TT'(ψ_N) are parameterized in terms of the normalized flux
   ψ_N = (ψ - ψ_axis)/(ψ_boundary - ψ_axis). Both ψ_axis and ψ_boundary depend on the
   current solution.

3. **The vacuum contribution is precomputed but the plasma contribution varies**: The
   total ψ = ψ_vacuum + ψ_plasma. The vacuum ψ is computed from known coil currents
   via precomputed mutual inductance matrices. The plasma contribution requires solving
   the GS problem.

LIUQE handles this through an outer fixed-point (Picard) iteration or Newton/SQP
iteration on the full nonlinear problem. At each outer iteration, the domain (which
cells are plasma) is updated, and the measurement Jacobian is relinearized.

---

## 2. LIUQE Algorithm: Calling Chain and Iteration Structure

### 2.1 Top-Level Calling Scheme

```
liuqe(shot, t, ...)
  → liu('TCV', shot, t, ...)
      → liuptcv(shot, ...)    % TCV-specific parameters → struct P
      → liugtcv(shot, P)      % TCV geometry → struct G
      → liup(P)               % validate/fill default params
      → liug(G, P)            % validate/fill geometry → G
      → liuc(P, G)            % assemble ancillary data → struct L
      → liuxtcv(shot, t, L)   % load measurements → struct LX
      → liux(L, LX)           % validate measurements
      → liut(L, LX)           % RECONSTRUCTION → struct LY
```

The structure `L` is the master parameter/geometry object passed throughout.
`L.P` = parameters, `L.G` = geometry, additional fields = precomputed matrices.

### 2.2 The Reconstruction Loop in liut.m

`liut.m` iterates over time slices. For each time slice `kt`:

**Step 1: Extract single-slice measurements.**
```matlab
LXt = meqxk(LX, kt);
```

**Step 2: Build the linearized optimization problem.**
```matlab
[s, lin, IDs, ws] = liulsqinit(L, LXt, s, liusolver);
```
This assembles the linearized measurement equations, constraint matrices,
and regularization terms for the current time slice.

**Step 3: Initial guess.**
```matlab
[s, L, LXt, LYp, chih, zfun, sinit, z0] = liuinit(...);
```
Two strategies controlled by `L.P.initguess`:
- `'LIH'` (default for first slice or after reset): A finite-element (FE) least-squares
  fit of the plasma current distribution is performed using `lihlsqsolve`. This gives an
  initial plasma current on the FE grid (`Iy`), which initializes the state vector.
- `'previous'`: The solution from the previous time slice is reused as the starting
  point, bypassing the FE fit. Used for time series to warm-start the solver.

**Step 4: Nonlinear solve.**
```matlab
[s, LYs, conv, nFeval, dzg] = meqopt(L, LXt, LYp, s, sinit, ...);
```
This is the core nonlinear solver. It calls `fgeF` (the GS residual function) repeatedly
inside a Newton/SQP/Picard loop until convergence.

**Step 5: Post-processing.**
After convergence, measurement fits are updated (including the dz stabilization term),
and the output structure `LYt` is assembled and appended to `LY`.

**Convergence criterion**: `L.P.psichco` (default 5×10⁻⁵ for TCV). This is the
relative change in poloidal flux ψ between iterations.

### 2.3 The Nonlinear Solver (meqopt / algoNL)

The outer nonlinear iteration strategy is controlled by `L.P.algoNL`:
- `'all-nl'` (default): The full nonlinear problem with plasma current density `Iy` as
  the state variable. This is the standard LIUQE mode.
- `'Newton-GS'` [GUESS]: A Newton-on-GS variant that uses the GS flux `Fx` as state.
- `'all-nl-Fx'` [GUESS]: Similar to all-nl but with Fx state.

The solver choice is controlled by `L.P.solver`:
- `'SQP'` (default): Sequential Quadratic Programming. A Picard warm-start phase
  (`L.P.wsSQP` iterations, default 3) is followed by SQP iterations. SQP provides
  better convergence near the solution and handles inequality constraints on the
  current profile coefficients.
- `'SQP-broyden'`: SQP with Broyden quasi-Newton updates for the Jacobian. Efficient
  for time series where solutions are close between time slices. Requires
  `initguess='previous'`.

In real-time mode (`L.P.itert > 0`), only `L.P.itert` iterations are performed per
time slice (typically 1), and the previous solution is always used as warm start.
In offline mode (`L.P.itert = 0`), iterations continue until `psichco` is satisfied
or `itera` iterations are reached (default 50 for TCV).

### 2.4 The Linearized Least-Squares Problem

At each outer (Picard/Newton) iteration, the plasma boundary and ψ-normalization are
fixed. The problem reduces to finding the coefficient vector **a** = [Ie; ag] where:
- `Ie` = external (coil + vessel) currents [ne elements]
- `ag` = basis function amplitude coefficients [ng elements, with ng = nq LIUQE base functions]

The linearized system being solved is a weighted least squares:

```
min_{a,Co,dz}  ||W_r (M_ry Iy(ag) + M_re Ie - Y_r)||²
              + ||W_e (Ie - Y_e)||²
              + ||W_i (diag terms - Y_i)||²
              + ||W_q (Q_qg ag    - X_q)||²
              + ||W_j (Co - Y_j)||²
```

subject to constraints on Co (equality constraints for exactly-constrained coils and
optional inequality constraints from `bfct` mode 7).

where:
- `W_r` = weights for flux loop + probe measurements
- `M_ry` = response matrix: measurement response to plasma current grid Iy [nr×ny]
- `M_re` = response matrix: measurement response to external currents [nr×ne]
- `Y_r` = measured flux loop / probe values
- `W_e`, `Y_e` = coil current constraint terms
- `W_i`, `Y_i` = plasma current Ip and diamagnetic flux terms
- `W_q` = regularization weight matrix (from `L.P.wreg`)
- `Q_qg`, `X_q` = regularization constraints from `bfct` mode 6
- `W_j`, `Co`, `Y_j` = constraints on base function coefficients ag

The `M_ry` matrix (called `L.Mry` in the code) maps plasma current distribution
`Iy` on the y-grid to measurement space. It is assembled in `liuc.m`:
```matlab
L.Mry = [G.Mfx(:,L.lxy) ; G.Bmx(:,L.lxy)];
```
where `Mfx` is the mutual inductance between flux loops and the x-grid and `Bmx`
between probes and the x-grid. `L.lxy` is the mask of inner (non-boundary) grid cells.

The `M_re` matrix maps external currents to measurements:
```matlab
L.Mre = [G.Mfa G.Mfu ; G.Bma G.Bmu];
```
where `Mfa`/`Mfu` = coil/vessel mutual inductances to flux loops,
`Bma`/`Bmu` = coil/vessel response at magnetic probe locations.

### 2.5 Block Solver Structure (iters)

When `L.P.iters > 0` (default 48), a block-iterative solver is used that alternates
between fitting:
1. The external currents `Ie` (coil + vessel) given current plasma Iy
2. The plasma current profile coefficients `ag` given current Ie

This block-Jacobi / block-Gauss-Seidel structure (implemented in `bslv`-related
routines) is important for the real-time version where the number of solver iterations
per time step is strictly limited.

### 2.6 State Scaling

When `L.P.xscal = true` (default), the state vector is scaled by characteristic
amplitudes stored in `L.xscal`. This is critical for numerical conditioning. The
characteristic scales are set in `meqalgoNLc.m` based on typical plasma current
amplitudes.

---

## 3. Base Functions for Current Profiles

### 3.1 The Base Function Concept

LIUQE represents p'(ψ) and TT'(ψ) as a linear superposition of basis functions:

```
p'(ψ_N)  = Σ_g  a_g · g'_P(ψ_N)
TT'(ψ_N) = Σ_g  a_g · g'_T(ψ_N)
```

where ψ_N = (ψ - ψ_axis)/(ψ_boundary - ψ_axis) ∈ [0,1] is the normalized poloidal
flux (0 on axis, 1 on boundary). Each basis function `g'` is assigned to either the
p' expansion or the TT' expansion via Boolean masks `fPg` (flag for p') and `fTg`
(flag for TT'). One basis function can contribute to both if needed [GUESS].

The total plasma toroidal current density is then:

```
j_φ(R,Z) = μ₀ R p'(ψ_N(R,Z)) + TT'(ψ_N(R,Z)) / R
          = Σ_g a_g [μ₀ R g'_P(ψ_N) + g'_T(ψ_N)/R]
```

This expresses j_φ as linear in the coefficients a_g. The plasma current distribution
on the y-grid is:

```
Iy = TYG * a_g
```

where `TYG` is the mode-1 basis function matrix [ny × ng] returned by the bfct function
call with mode=1. `TYG[k,g]` = value of basis function g (multiplied by R or 1/R) at
grid cell k.

### 3.2 The bfct Function Interface

All base function sets use the same calling interface (described in `bfhelp.m`):

```
[FPG, FTG, TDG]         = BF(0, PAR)                          % masks and domain matrix
[TYG, TPDG, ITPDG]      = BF(1, PAR, F, F0, F1, O, RY, IRY)  % basis on y-grid
[GN, IGN]               = BF(2, PAR, FN, F0, F1)              % basis on normalized flux grid
[APP, ATTP, AP, AHQT]   = BF(3, PAR, A, F0, F1, FP, FT, IDS) % profile coefficients
[QQG, XQ]               = BF(6, PAR, ~, F0, F1, R0, IR0, IDS) % regularization
[QCG, XC]               = BF(7, PAR, ~, F0, F1)               % inequality constraints
[DGY, ...]              = BF(11, PAR, F, F0, F1, ...)          % derivatives of TYG
[DGDF, ...]             = BF(12, PAR, F, F0, F1)               % derivatives of GN, IGN
```

Key arguments:
- `PAR`: basis function parameters (varies by set, e.g. `P.bfp`)
- `F, F0, F1`: flux values (F=flux map, F0=axis flux, F1=boundary flux)
- `O`: Boolean mask for plasma interior
- `RY, IRY`: R values and 1/R values on y-grid

The code handle `P.bfct` (e.g. `@bfabmex`, `@bf3pmex`, `@bf3imex`) is called via
`L.bfct` in all reconstruction routines.

### 3.3 bfabmex: Polynomial Basis Functions (Default)

`bfabmex` with `P.bfp = [nP, nT]` produces `nP + nT` polynomial basis functions.
For p' it uses the first `nP` of:
```
g'_1(ψ) = ψ - ψ_B                                      (linear, zero at boundary)
g'_2(ψ) = (ψ - ψ_B)(ψ - ψ_A)                           (quadratic, zero at both ends)
g'_3(ψ) = (ψ - ψ_B)(ψ - ψ_A)(2ψ - ψ_B - ψ_A)          (cubic)
```
For TT' it uses the first `nT` of the same polynomials. These are also called "BFAB"
functions (BF Algebraic/polynomial). All are zero at the boundary ψ = ψ_B; p' has a
free value at the axis.

Default TCV configuration: `bfp = [1, 2]`, giving 3 total basis functions:
- 1 for p': `g'_P1 = (ψ - ψ_B)`
- 2 for TT': `g'_T1 = (ψ - ψ_B)` and `g'_T2 = (ψ - ψ_B)(ψ - ψ_A)`

Note: the `bfabmex` file in the repository stubs out with an error; the actual
implementation is a compiled MEX file. The MATLAB source equivalent is in `mexm/bfab.m`.

### 3.4 bf3pmex: 3 Polynomial Basis Functions (Real-Time Default)

`bf3pmex` with `P.bfp = false` (no p') produces:
- 1 function for TT': `g'_T1 = ψ - ψ_B` (linear)
- 1 function for TT': `g'_T2 = (ψ - ψ_B)(ψ - ψ_A)` (quadratic)
- With `P.bfp = true` (default): additionally `g'_P1 = ψ - ψ_B` for p'

This is the default for Simulink/real-time: `D.bfct = @bf3pmex; D.bfp = false` when
`P.slx = 1`. The simplified profile parametrization reduces computation time.

### 3.5 bf3imex / bfli: Tabular Interpolation Basis Functions

`bf3imex` uses user-supplied tabulated basis functions. The parameter `P.bfp` must be
a struct with:
- `gNg(n, 1:3)`: 3 basis functions tabulated at `linspace(0,1,n)` in normalized flux
- `IgNg(n, 1:3)`: primitives of gNg (precomputed by `bfprmex`)
- `fPg(1:3)`: boolean mask which of the 3 are for p'
- `fTg(1:3)`: boolean mask which are for TT'

`bfli.m` is the pure-MATLAB equivalent of `bf3imex` (the tabular interpolation MEX).
The interpolation uses piecewise-linear interpolation of the tabulated profiles.
From `bfli_fN`:
```matlab
k = min(max(ifloor(FN, FQ, 1), 1), nN-1);
f = FQ*idFN - k + 1;   % fractional position within interval
gQ = P.gNg(k,:).*f1 + P.gNg(k+1,:).*f;
```

### 3.6 Physical Constraints on Profiles

Mode 6 (`BF(6,...)`) returns regularization constraints `Q_qg * a_g = X_q` that
penalize physically undesirable profile shapes (e.g., large oscillations). This
enters the weighted least-squares system with weight `L.P.wreg`.

Mode 7 (`BF(7,...)`) returns inequality constraints `X_c ≤ Q_cg * a_g` that can
enforce positivity of pressure gradient (p'≤0 everywhere [GUESS], since ψ increases
toward axis) or monotonicity constraints.

The number of basis functions `L.ng` equals `L.nq` (number of regularization
constraints) which is determined by the bfct choice.

---

## 4. Measurement Model

### 4.1 Measurement Types

For TCV, the measurements loaded by `liuxtcv.m` are:

| Field | Physical quantity | Count (TCV) | Unit |
|-------|------------------|-------------|------|
| `LX.Ff` | Poloidal flux at flux loop positions | 38 | Wb |
| `LX.Bm` | Poloidal magnetic field at probe positions | 38 | T |
| `LX.Ia` | Active coil currents | 19 | A |
| `LX.Iu` | Vessel segment currents (from voltage loops) | 38 | A |
| `LX.Ip` | Plasma current (Rogowski coil) | 1 | A |
| `LX.Ft` | Diamagnetic flux (corrected) | 1 | Wb |
| `LX.Xt` | Diamagnetic flux (before Iy correction) | 1 | Wb |
| `LX.rBt` | Vacuum toroidal field × R | 1 | m·T |
| `LX.Uf` | Loop voltages (used to compute Iu) | 38 | V |

For TCV, flux loops are read via `TCV_IDEALLOOP("FLUX")` (ideal loop conversion),
magnetic probes from `\MAGNETICS::BPOL_003`, active coil currents from
`\MAGNETICS::IPOL`, and the diamagnetic loop from `\RESULTS::DMLCOR`.

Plasma current is not directly measured as a diagnostic current but estimated from
the magnetic probe array using a boundary integration (Ampere's law):
```matlab
X.Ip = L.G.Ipm * X.Bm + L.Ipa * X.Ia + L.Ips * Is;
```
where `L.G.Ipm` is a vector of arc-length weights around the boundary.

### 4.2 Forward Model for Measurements

Given flux distribution ψ(R,Z) (encoded as `Fx` on the x-grid), the predicted
measurements are:

**Flux loops** (indexed `kdf`):
```
Ff_pred = Mfx * Fx  (mutual inductance × flux on grid)
```
plus contributions from external coils via `Mfa * Ia + Mfu * Iu`.

**Magnetic probes** (indexed `kdm`):
```
Bm_pred = Bmx * Fx
```
where `Bmx` encodes the field component orientation (cos(θ_m) B_R + sin(θ_m) B_Z at
probe angle θ_m). The response is obtained from ∂ψ/∂R and ∂ψ/∂Z at probe locations.

**Plasma current** (indexed `kdp`):
```
Ip_pred = L.Wip * Iy   (integral of current on y-grid)
```

**Diamagnetic flux** (indexed `kdt`):
If `idml=1`:
```
Ft_pred = Xt + Mty * Iy   (raw DML + plasma-current correction)
```
where `Mty` = DML plasma coupling matrix (loaded from MDS `\RESULTS::DML:BO###`).
The `Xt` field is the DML measurement before the plasma current correction.

The full predicted measurement vector is assembled as:
```
Y_pred = Mry * Iy + Mre * Ie
```
where `Iy` = plasma current on y-grid [ny × 1], `Ie` = [Ia; Iu] [ne × 1].

### 4.3 Measurement Weighting and Covariance

Measurements are weighted by their assumed noise standard deviations:
```matlab
Wf = wFf / Fferr;   % flux loop weight
Wm = wBm / Bmerr;   % probe weight
Wa = wIa / Iaerr;   % coil current weight
Wu = wIu / Iuerr;   % vessel current weight
```

TCV defaults (from `liuptcv_default`):
- `Fferr = 0.0012 * 2π ≈ 7.5e-3 Wb` (flux loop error)
- `Bmerr = 0.01 T` (probe error)
- `Iaerr = 200 A` (coil current error)
- `Iuerr = 200 A` (vessel segment current error)
- `Iperr = 24000 A` (plasma current error)

Per-channel weights `wFf`, `wBm`, `wIa`, `wIu` (each a vector of 38/38/19/38 values)
allow individual channels to be down-weighted or turned off. For TCV, flux loops 10–12
have `wFf = 0` by default (disabled due to calibration issues).

Setting a weight to `Inf` makes the corresponding current an exact (equality)
constraint: `ke = ~isinf(We)` identifies approximate constraints, while
`~ke` elements are exact constraints and excluded from the residual sum.

### 4.4 Measurement Jacobian

For the linearized problem, the Jacobian of measurements with respect to the
coefficient vector a_g is:

```
∂Y_r/∂a_g = Mry * (∂Iy/∂a_g) = Mry * TYG
```

where `TYG[ny × ng]` is the mode-1 output of the bfct function (basis functions
evaluated at all y-grid points, multiplied by the appropriate R or 1/R factors).
This matrix `Wrh = L.Wr .* Mry * TYG` is called `L.Wrh` (approximately; final
assembly in `liuc.m`) and forms the core of the measurement equation.

The Jacobian also includes derivatives with respect to Ie (direct measurement response
`Mre`) and the sensitivity of TYG to ψ normalization changes (handled in bfct mode 11).

---

## 5. Geometry Setup

### 5.1 The x-Grid and y-Grid

MEQ uses two overlapping grids:
- **x-grid** (nzx × nrx): the full computational grid including the boundary row/column.
  `rx`, `zx` are the R, Z coordinates. For TCV: nzx = nz = 32, nrx = nr = 27 by default.
- **y-grid** (nzy × nry = (nzx-2) × (nrx-2)): the inner grid, excluding boundary cells.
  `ry`, `zy` are coordinates. The GS equation is solved on the y-grid; boundary
  conditions are imposed via `Fx` on the x-boundary.

The mask `L.lxy` [nzx×nrx logical] identifies which x-grid cells correspond to y-grid
cells (i.e., the inner cells). `Iy[ny]` is the plasma current on the y-grid (vectorized).

Note: `nz` must be a power of 2 (checked in `liup.m`) for FFT-based solvers [GUESS].

### 5.2 Electromagnetic Matrices in G

The geometry structure `L.G` (originally assembled by `meqgtcv`/`liugtcv`) stores
precomputed electromagnetic coupling matrices between coils, vessel, plasma grid, and
diagnostics. All are computed once at setup time and reused every time slice:

| Matrix | Dimensions | Meaning |
|--------|-----------|---------|
| `G.Mfx` | [nf × nx] | Mutual inductance: flux loops ↔ x-grid cells |
| `G.Mfa` | [nf × na] | Mutual inductance: flux loops ↔ active coils |
| `G.Mfu` | [nf × nu] | Mutual inductance: flux loops ↔ vessel segments |
| `G.Maa` | [na × na] | Self/mutual inductance among active coils |
| `G.Muu` | [nu × nu] | Self/mutual inductance among vessel segments |
| `G.Mau` | [na × nu] | Cross inductance coils ↔ vessel |
| `G.Bmx` | [nm × nx] | Field response: probes ↔ x-grid cells |
| `G.Bma` | [nm × na] | Field response: probes ↔ active coils |
| `G.Bmu` | [nm × nu] | Field response: probes ↔ vessel |
| `G.Mtx` | [nx × nx] | DML plasma coupling matrix (loaded per shot) |
| `G.Mxu` | [nx × nu] | Mutual inductance: x-grid ↔ vessel |
| `G.Mxa` | [nx × na] | Mutual inductance: x-grid ↔ active coils |
| `G.Ra`  | [na]      | Active coil resistances |
| `G.Ru`  | [nu]      | Vessel segment resistances |
| `G.Brxa`| [nb × na] | B_R at boundary: active coils ↔ boundary |
| `G.Brxu`| [nb × nu] | B_R at boundary: vessel ↔ boundary |

The `Brxa`/`Brxu` matrices are used for the vertical stabilization dz correction:
they encode how a vertical shift of all coils/vessel changes the flux on the boundary.

### 5.3 TCV-Specific Geometry (liugtcv.m)

For TCV, `liugtcv.m`:
1. Calls `meqgtcv` to load the machine description (coil positions, probe locations,
   limiter wall, etc.) and compute all electromagnetic matrices from the static database
   `\STATIC::*`.
2. Loads the `dzMfx` matrix (derivative of flux-loop mutual inductance with respect to
   vertical shift) via `STATIC("BR_F_X")` entries. This is used in the vertical
   stabilization equations.
3. Similarly loads `dzBmx` (derivative of probe response to vertical shift).
4. Loads the DML coupling matrix from `\RESULTS::DML:BO###` for the FIX-PLA-LIN
   compensation term.
5. Computes the plasma current estimator weights `G.Ipm` (arc-length weights around
   the boundary for Ampere's law). For TCV: 38 magnetic probes are used.

### 5.4 Vessel Representation

TCV vessel currents are represented as segment currents (resistive shells). The
`P.selu` parameter controls which representation is used:
- `'sx'` (default TCV): Extended segment representation with x-grid coupling
- `'s'`: Standard segment representation
- `'e'`: Explicit modes (used in some tests)

Vessel segment currents `Iu` are computed from the loop voltages: `Iu = Uf / (-Ru)`.

### 5.5 L Structure: Key Fields

After `liuc.m`, the assembled `L` structure contains (in addition to `L.P`, `L.G`):

| Field | Meaning |
|-------|---------|
| `L.rx, L.zx` | x-grid coordinates [nzx], [nrx] |
| `L.ry, L.zy` | y-grid coordinates [nzy], [nry] |
| `L.ny` | Total y-grid cells (nzy × nry) |
| `L.nN` | State vector dimension |
| `L.nuN` | Input vector dimension (Ie + Co) |
| `L.ndz` | Dimension of vertical shift state (0 or nD) |
| `L.nD` | Number of plasma domains (1 for standard, 2 for doublets) |
| `L.dlst` | Sparse GS operator matrix [ny × nx], assembled by `dlstc.m` |
| `L.Tye` | Sparse mapping: external currents → y-grid flux |
| `L.Mry` | Response: plasma y-grid → measurements [nr × ny] |
| `L.Mre` | Response: external currents → measurements [nr × ne] |
| `L.Mdzry`| Response: vertical shift → measurements [nr × ny] |
| `L.Wr, L.We, L.Wi` | Diagonal measurement weight vectors |
| `L.Wq` | Regularization weight [nq × 1] |
| `L.Wj` | Constraint weight for Co and dz |
| `L.Mty` | DML plasma current coupling vector [1 × ny] |
| `L.Tyh` | FE basis functions: y-grid → FE elements [ny × nh] |
| `L.Wdh, L.Wde` | FE fit matrices (precomputed in `liuc.m`) |
| `L.Ahd, L.Aed` | FE direct solver matrices |
| `L.bfct` | Function handle to active basis function set |
| `L.bfp` | Parameters for bfct |
| `L.ind` | Index structure for state vector components |
| `L.xscal` | State vector scaling factors |
| `L.LX2x` | Function handle: LY → state vector x |
| `L.Jx, L.Ju` | Vacuum Jacobian (linear part of ∂F/∂x, ∂F/∂u) |

---

## 6. Parameters Reference

### 6.1 Key Algorithmic Parameters (liup.m / liuptcv.m)

**Grid and domain:**
| Parameter | TCV Default | Meaning |
|-----------|-------------|---------|
| `nr` | 27 | Number of R grid points |
| `nz` | 32 | Number of Z grid points (must be power of 2) |
| `nelem` | 10 | Number of FE elements for initial current guess |
| `fetype` | 'bil' | FE element type: 'bil' bilinear, 'pyr' pyramidal, 'tri' triangular |

**Iteration control:**
| Parameter | TCV Default | Meaning |
|-----------|-------------|---------|
| `itera` | 50 | Max outer iterations (offline mode) |
| `iterfrz` | 50 | Max iterations with frozen plasma domain |
| `itert` | 0 | Real-time: iterations per time slice (0 = offline) |
| `iters` | 48 | Block solver iterations |
| `iterh` | 30 | FE nonlinear solver iterations for initial guess |
| `tolh` | 1e-3 | FE solver tolerance |
| `psichco` | 5e-5 | Convergence criterion on ψ change |

**Measurement weights:**
| Parameter | TCV Default | Meaning |
|-----------|-------------|---------|
| `Fferr` | 7.54e-3 Wb | Flux loop noise σ |
| `Bmerr` | 0.01 T | Magnetic probe noise σ |
| `Iaerr` | 200 A | Coil current noise σ |
| `Iuerr` | 200 A | Vessel current noise σ |
| `Iperr` | 24000 A | Plasma current noise σ |
| `wFf` | ones(38,1), [10:12]=0 | Per-channel flux loop weights |
| `wBm` | ones(38,1) | Per-channel probe weights |
| `wIa` | ones(19,1) | Per-coil weight |
| `wIu` | ones(38,1) | Per-vessel-segment weight |
| `wIp` | 1 | Plasma current weight |

**Regularization:**
| Parameter | TCV Default | Meaning |
|-----------|-------------|---------|
| `wreg` | 2e-6 | Regularization weight on base function amplitudes |
| `wregadapt` | 15 | Decay rate of wreg with elongation |
| `elomin` | 1.1 | Min elongation below which wreg stays constant |

**Base functions:**
| Parameter | TCV Default | Meaning |
|-----------|-------------|---------|
| `bfct` | `@bfabmex` | Basis function set |
| `bfp` | [1;2] | nP=1 for p', nT=2 for TT' |
| `nFW` | 1 | Number of flux surfaces for profile output |

**Solver:**
| Parameter | Default | Meaning |
|-----------|---------|---------|
| `algoNL` | 'all-nl' | Nonlinear algorithm variant |
| `useSQP` | false | Use SQP vs Picard only |
| `solver` | 'SQP' | Solver type |
| `wsSQP` | 3 | Picard warm-start iterations before SQP |
| `initguess` | 'LIH' | Initial guess source |
| `xscal` | true | Scale state vector |
| `tolstep` | 1e-8 | Step size tolerance |

### 6.2 Regularization Weight Adaptation

The regularization weight is adapted based on plasma elongation κ:
```
wreg_eff = wreg * exp(-wregadapt * max(kappa - elomin, 0))
```
[GUESS based on parameter names; exact formula in meqc or liuc routines not shown].
The idea is that elongated plasmas have more constrained current profiles (fewer free
parameters needed), so regularization can be reduced.

### 6.3 Convergence and Acceptance

After reconstruction, the quality is assessed by `L.P.chimax` (default 3 for TCV):
a root-mean-square residual above this threshold causes the equilibrium to be flagged
as not converged. The fields `LY.chi`, `LY.chie`, `LY.werr`, `LY.uerr` track different
components of the fit residual.

---

## 7. Vertical Stabilization

### 7.1 Why Vertical Position Needs Special Treatment

Vertically elongated plasmas are unstable to rigid vertical displacement. Small
perturbations in Z are amplified by the surrounding conducting structures. The position
of the magnetic axis Z_A (and by extension the plasma boundary shape) is therefore
sensitive to the initial conditions and the iterative solution path.

In the reconstruction problem, an ambiguity arises: a plasma shifted vertically by δz
from its "true" position produces measurements that are almost identical to a plasma at
the correct position but with slightly different current profile. The reconstruction
algorithm can converge to an incorrect vertical position if the initial guess is poor.

LIUQE addresses this by introducing a **vertical shift dz** as an additional degree of
freedom in the cost function. The shift dz changes all flux loop and probe predictions
by:
```
δF_f = dzMfx * Iy * dz    (derivative of flux w.r.t. vertical shift)
δB_m = dzBmx * Iy * dz
```
This is assembled in `liuc.m` via:
```matlab
L.Mdzry = [G.dzMfx(:,L.lxy) ; G.dzBmx(:,L.lxy)];
```

### 7.2 Vertical Stabilization Functions (liudz, liuz0)

The stabilization logic is encapsulated in the function handle `P.zfct`:

**Default (`P.zfct = @(P,rst,z0,zA,dz) dz`)**: The fitted shift `dzg` from the
least-squares solve is used directly as `dz`. This is the "passive" mode with no
feedback.

**`liudz`** (proportional-integral feedback on dzg):
```matlab
dz = P.stab * dzg + P.stabint * zdzg;   % PI on fitted shift
zdzg = zdzg + dzg;                       % integrate dzg
```
Used when the fitted shift itself provides a useful feedback signal.

**`liuz0`** (PID feedback on Z_A - Z_0 error):
```matlab
S = (zA - Z0) + [0; S(2); -S(1)];  % PID states
dz = K * S;                          % PID output
Z0 = Z0 - dz * P.stabcint;          % update reference
```
where Z_0 is a reference vertical position and Z_A is the reconstructed magnetic axis
position. Parameters: `P.stab` (proportional), `P.stabint` (integral),
`P.stabder` (derivative), `P.stabcint` (integral of dz→Z_0 correction).

The condition `P.stabz ~= 0` must be true to include dz as a degree of freedom
(`L.ndz = nD` rather than 0). If `stabz = 0`, vertical stabilization is disabled and
`L.ndz = 0`.

Note: vertical stabilization (`stabz`) is only available when `algoNL='all-nl'`.

---

## 8. Real-Time Operation (liusimc)

### 8.1 Real-Time vs. Offline Mode

Setting `P.itert = k > 0` enables real-time mode where exactly `k` iterations are
performed per time slice. The Simulink real-time model uses `P.slx = 1`.

In real-time mode:
- `P.initguess` is forced to `'previous'` (no FE fit, warm-start from last slice)
- Only `bf3pmex` or `bf3imex` are supported as bfct (MEX files, fast)
- All coil/vessel weights must be finite (no exact equality constraints via Inf)
- `ilim = 1` (limited configuration assumed)
- The Picard block solver is used rather than full SQP

### 8.2 Simulink Model Setup (liusimc.m)

`liusimc` prepares the MATLAB workspace or Simulink data dictionaries for running the
Simulink models `liutslx.slx` (reconstruction) and `liuxslx.slx` (data source). Steps:
1. Sets variant subsystem switches (BFCTMODE, HFCTMODE, IPMODE, BSLVMODE, etc.) based
   on current `L.P` settings.
2. Creates Simulink Bus objects for `LX` and `LY` structures (with specific field lists
   and single-precision casting).
3. Removes sparse matrices from `L` that can't be represented in Simulink fixed types
   (`Tye`, `Txy`, `Txb`, `dlst`).
4. Casts `L` entries to single precision.
5. Converts function handles in `L.P.bfct`, `L.bfct`, `L.P.zfct` to strings for
   data dictionary storage.
6. Writes data dictionaries (`*.sldd`) linking `L`, `LX`, bus objects.

The LX Simulink bus includes: `{Ff, Bm, Ip, Ie, rBt, Xt}` — minimal set for RT.
Note that `Ie = [Ia; Iu]` is the combined external current vector.

### 8.3 Real-Time Performance Considerations

[GUESS based on code structure and Moret 2015 paper]:
- The reconstruction loop in Simulink is a discrete-time system with a fixed step
  equal to the magnetics acquisition time (100 µs for TCV at 10 kHz).
- Typical real-time LIUQE on TCV runs at ~1 ms update rate (10 magnetics samples per
  update), using `itert=1` (one Picard iteration per update).
- The block-solver (`iters`) with 48 iterations resolves the coil current / profile
  amplitude coupling without a full matrix factorization each step.
- Speed-limiting steps: basis function evaluation (MEX), boundary finding (`meqpdom`),
  and the measurement assembly and solve.

---

## 9. Output Structure (LY)

The output structure `LY` from `liut.m` contains reconstructed equilibrium quantities
at each time slice. Fields are documented in `meqt.m` (generic) and the LIUQE-specific
fields in `liut.m`.

### 9.1 Primary Equilibrium Fields

| Field | Size | Unit | Meaning |
|-------|------|------|---------|
| `LY.t` | [1,nt] | s | Time base |
| `LY.Fx` | [nzx,nrx,nt] | Wb/rad | Poloidal flux on x-grid ψ(R,Z) |
| `LY.Iy` | [nzy,nry,nt] | A/m² [GUESS] | Plasma toroidal current density on y-grid |
| `LY.Ia` | [na,nt] | A | Fitted active coil currents |
| `LY.Iu` | [nu,nt] | A | Fitted vessel segment currents |
| `LY.ag` | [ng,nt] | — | Basis function amplitude coefficients |
| `LY.Ip` | [1,nt] | A | Total plasma current |
| `LY.rIp` | [1,nt] | m | R-position of plasma current centroid |
| `LY.zIp` | [1,nt] | m | Z-position of plasma current centroid |

### 9.2 Plasma Boundary

| Field | Size | Unit | Meaning |
|-------|------|------|---------|
| `LY.rA` | [nD,nt] | m | R of magnetic axis (per domain) |
| `LY.zA` | [nD,nt] | m | Z of magnetic axis |
| `LY.FA` | [nD,nt] | Wb/rad | ψ at magnetic axis |
| `LY.rB` | [nD,nt] | m | R of LCFS reference point |
| `LY.zB` | [nD,nt] | m | Z of LCFS reference point |
| `LY.FB` | [nD,nt] | Wb/rad | ψ on LCFS (ψ_boundary) |
| `LY.rX` | [nX,nt] | m | R of X-point(s) |
| `LY.zX` | [nX,nt] | m | Z of X-point(s) |
| `LY.FX` | [nX,nt] | Wb/rad | ψ at X-point(s) |
| `LY.lB` | [1,nt] | — | Limiter/X-point flag (0=limited, 1=diverted) |
| `LY.FW` | [nrx,nt] [GUESS] | m | Flux at wall (first wall) |

### 9.3 Profile Quantities

| Field | Size | Unit | Meaning |
|-------|------|------|---------|
| `LY.aq` | [noq,nt] | various | q-profile on noq flux surfaces |
| `LY.aW` | [noq,nt] | various | Profile integrals on noq flux surfaces |
| `LY.PpQg` | [noq,ng] | various | p'(ψ_N) from basis functions |
| `LY.TTpQg`| [noq,ng] | various | TT'(ψ_N) from basis functions |
| `LY.TQ` | [nthg,nt] | T·m | Toroidal field T(ψ) = R·B_φ on ψ-grid |
| `LY.pQ` | [noq,nt] | Pa | Pressure profile |

### 9.4 Global Scalar Quantities

| Field | Size | Unit | Meaning |
|-------|------|------|---------|
| `LY.beta` | [1,nt] | — | Toroidal beta = 2μ₀<p>/B₀² |
| `LY.betap` | [1,nt] | — | Poloidal beta |
| `LY.li` | [1,nt] | — | Internal inductance per unit length |
| `LY.Vp` | [1,nt] | m³ | Plasma volume |
| `LY.Wp` | [1,nt] | J | Stored energy = (3/2) ∫p dV |
| `LY.qA` | [nD,nt] | — | Safety factor at magnetic axis |
| `LY.kappa`| [1,nt] | — | Elongation |
| `LY.delta`| [1,nt] | — | Triangularity |

### 9.5 Diagnostic / Quality Fields

| Field | Size | Unit | Meaning |
|-------|------|------|---------|
| `LY.chi` | [1,nt] | — | RMS fit residual (overall) |
| `LY.chie` | [1,nt] | — | RMS constraint residual |
| `LY.werr` | [nd,nt] | — | Weighted per-measurement residual |
| `LY.uerr` | [nd,nt] | — | Unweighted per-measurement residual |
| `LY.chih` | [1,nt] | — | FE fit residual (initial guess quality) |
| `LY.dz` | [nD,nt] | m | Vertical stabilization shift |
| `LY.dzg` | [nD,nt] | m | Fitted vertical shift |
| `LY.niter`| [1,nt] | — | Iterations to convergence |
| `LY.cycle`| [1,nt] | µs | Computation time |
| `LY.rst` | [1,nt] | — | Reset flag (1 = reconstruction restarted) |
| `LY.isconverged` | [1,nt] | — | Convergence flag |
| `LY.Ff` | [nf,nt] | Wb | Synthetic flux loop signals |
| `LY.Bm` | [nm,nt] | T | Synthetic probe signals |
| `LY.Xt` | [1,nt] | Wb | Synthetic diamagnetic flux |
| `LY.ID` | [nd] | — | Equation ID labels for cost terms |

The `isconverged` flag identifies slices where convergence was declared. By default
only converged slices are returned (controlled by `L.P.LYall`).

---

## 10. Physics Accuracy and Tradeoffs

### 10.1 What LIUQE Gets Right

- **Full GS equilibrium**: Unlike linear extrapolation codes, LIUQE solves the full
  2D GS equation at each iteration. The flux surface geometry (elongation, triangularity,
  X-point presence) is correctly captured.
- **Toroidal current distribution**: The polynomial/tabular basis function expansion
  gives a physically motivated parametrization of j_φ(ψ).
- **Simultaneous coil-current and profile fitting**: External coil currents (including
  vessel eddy currents) are fitted simultaneously with the profile, avoiding systematic
  biases from fixed-current assumptions.
- **Diamagnetic flux**: When `idml=1`, the DML measurement constrains the diamagnetic
  energy integral, connecting the pressure profile to a global measurement.

### 10.2 Approximations and Limitations

**Profile flexibility**: With default `bfp=[1,2]` (3 basis functions: 1 for p', 2 for
TT'), the profile parametrization is very rigid. This is appropriate for a
magnetics-only reconstruction where the profile information is limited, but cannot
capture complex profiles (hollow current, multiple peaks). More basis functions can be
added via `bfp=[2,3]` etc., but stability decreases.

**Toroidal symmetry**: The GS equation assumes perfect toroidal symmetry. In reality,
the TCV vacuum vessel is a 3D structure and coils have finite extent. The
`TCV_IDEALLOOP` correction converts discrete flux loops to idealized toroidal integrals,
partially mitigating this, but residual 3D effects are not captured.

**Current profile self-consistency**: The GS equation solution at each outer iteration
assumes the ψ_N normalization from the previous iteration. This fixed-point approach
converges but is not guaranteed to converge to the global minimum.

**Vessel current modeling**: Vessel eddy currents are modeled as lumped resistive
segments (38 for TCV). The actual continuous vessel has richer current distributions
not captured by this model.

**Pressure gradient only**: p'(ψ) is what enters the GS equation and is what LIUQE
fits. The actual pressure profile p(ψ) requires integration and has an undetermined
additive constant (determined by the boundary condition p(ψ_B) = 0, i.e., zero pressure
at the plasma boundary). This is a standard approximation.

**No kinetic / transport physics**: LIUQE is purely a magnetic equilibrium code.
It has no knowledge of transport, neoclassical physics, fast particle effects, or
non-equilibrium states. For time series, each time slice is reconstructed independently.

### 10.3 Speed vs. Accuracy Tradeoffs

| Parameter | Fast/RT | Accurate/Offline |
|-----------|---------|-----------------|
| `bfct` | `@bf3pmex` (2 functions) | `@bfabmex [1,2]` (3 functions) |
| `itera` | 1 | 50 |
| `nelem` | 6 | 10 |
| `iterh` | 0 | 30 |
| `initguess` | `'previous'` | `'LIH'` |
| `useSQP` | false | true |
| `nz×nr` | 32×27 | 32×27 (same for TCV) |

The real-time version accepts larger fit residuals in exchange for completing within
the ~1 ms cycle time. The offline version iterates to tight convergence and uses a
more flexible initial guess.

---

## 11. Issues for Python Rewrite

### 11.1 Core Computational Structure

The reconstruction at each time slice is:
1. Initial guess via FE least-squares (can be done with `scipy.linalg.lstsq`)
2. Outer iteration: boundary detection + ψ renormalization + basis function evaluation
3. Inner linear solve: weighted least squares (can use `numpy.linalg.lstsq` or
   `scipy.sparse.linalg.lsqr` for sparse systems)
4. Post-processing: flux surface tracing for q, beta, li, etc.

### 11.2 MEX File Replacements

The key MEX files that need Python equivalents:
- `bfabmex.c`: Polynomial basis functions. The formulas are given explicitly in
  `bfhelp.m`: `g'_k` are products of `(ψ - ψ_B)` and `(ψ - ψ_A)`. Pure NumPy.
- `bf3pmex.c`: Simplified polynomial set. Same approach.
- `bf3imex.c` / `bfli.m`: Piecewise-linear interpolation. `scipy.interpolate.interp1d`
  or direct NumPy indexing.
- `bfprmex.c`: Tabular integration of basis functions. Trapezoidal rule in NumPy.
- `dlstc.m`: Sparse GS operator. Already shows its construction clearly — use
  `scipy.sparse`.
- `bintc.m`: Bilinear interpolation. Already in MATLAB; direct NumPy equivalent.

### 11.3 Sparse vs. Dense Matrices

Most matrices in LIUQE are dense (mutual inductances, response matrices). The GS
operator `dlst` is sparse ([ny × nx] with 5 nonzeros per row). For a 32×27 grid,
`ny = 30×25 = 750` and `nx = 32×27 = 864`, so dense storage is feasible in Python.
The main advantage of sparse storage is in the GS solve: `scipy.sparse.linalg.spsolve`
or similar.

The measurement matrices `Mry` [nr × ny] with nr ≈ 76 (38 flux loops + 38 probes)
and ny ≈ 750 are about 57 kB dense — manageable.

### 11.4 The Iterative Solver

The Picard/SQP iteration in LIUQE maps to:
```python
for iteration in range(max_iter):
    Iy = TYG @ ag                          # Current from profile coefficients
    Fx = solve_GS(Iy, Ie)                  # GS solve for flux
    psi_A, psi_B = find_axis_boundary(Fx)  # Boundary detection
    TYG = bfct(Fx, psi_A, psi_B)          # Recompute basis functions
    # Build weighted LSQ system
    A_lsq = W @ M_ry @ TYG
    b_lsq = W @ (Y - M_re @ Ie)
    ag = np.linalg.lstsq(A_lsq, b_lsq)    # Solve for coefficients
    # Check convergence
    if converged(Fx, Fx_prev): break
```

The real SQP variant is more complex (handling inequality constraints from bfct mode 7),
but this captures the essential structure.

### 11.5 Boundary Detection

The `meqpdom` function (not read in detail but called from `bdet.m`) identifies the
plasma boundary (LCFS), X-points, and magnetic axis. This is a key computational step
involving:
- Finding the minimum/maximum of ψ (magnetic axis, for circular or elongated plasma)
- Tracing the LCFS as the outermost closed flux contour not crossing the limiter
- Detecting X-points as saddle points of ψ
- Computing ψ_N normalization

In Python this could be implemented using `scipy.interpolate.RectBivariateSpline` and
`scipy.optimize` for root finding and extremum detection.

### 11.6 Real-Time Constraints

The LIUQE RT constraints (1 iteration, fixed step, single precision) do not apply to a
Python offline reconstruction. A Python version would naturally:
- Use double precision throughout
- Run to convergence (many iterations)
- Not need Simulink bus structures or MEX compilation
- Use `scipy.optimize.minimize` or custom Newton for the nonlinear solve if SQP is needed

### 11.7 Data Access

Measurement loading (`liuxtcv.m`) uses MDSplus. In Python, the `mdsplus` Python
package provides direct access to TCV MDSplus trees. Key nodes needed:
- `\MAGNETICS::BPOL_003[dim,*]` — probes
- `TCV_IDEALLOOP("FLUX")[dim,*]` — flux loops
- `\MAGNETICS::IPOL[dim,*]` — coil currents
- `\RESULTS::DMLCOR` — diamagnetic flux
- `\MAGNETICS::RBPHI` — vacuum toroidal field

Geometry matrices should be loaded once from `\STATIC::*` via MDSplus and cached.

---

## 12. Mathematical Summary

### 12.1 The GS Operator

From `dlstc.m`, the GS operator Δ* is discretized as the sparse matrix `M` where
each row corresponds to an interior y-grid point (iz, ir):

```
(Δ*ψ)_{iz,ir} = [ -C * ψ_{iz,ir}
                 + ψ_{iz,ir+1}      (right)
                 + ψ_{iz,ir-1}      (left)
                 + A * ψ_{iz+1,ir}  (up)
                 + B * ψ_{iz-1,ir}  (down) ] / dz²
```

where:
```
A = (dz/dr)² * r / (r + dr/2)   (upward coefficient)
B = (dz/dr)² * r / (r - dr/2)   (downward coefficient)
C = 2 + A + B                   (diagonal)
```

This is a second-order finite difference discretization of the Grad-Shafranov operator.

### 12.2 Normalized Flux

The normalized flux $\psi_N$ used by basis functions is:
$$\psi_N = \frac{\psi - \psi_A}{\psi_B - \psi_A}$$

where $\psi_A = \psi$ at magnetic axis and $\psi_B = \psi$ at boundary. This maps axis to 0 and boundary to 1. Note: in many tokamak conventions $\psi_A > \psi_B$ ($\psi$ increases toward axis in TCV's ohmic phase), so $\psi_B - \psi_A < 0$. The basis functions are designed with $F_0=\psi_A$ and $F_1=\psi_B$ so that $g'_k(F_1) = 0$ (zero at boundary).

### 12.3 Picard Convergence

LIUQE's Picard iteration solves the sequence of linear problems:
1. Fix ψ^(n) → find domain Ω^(n) and ψ_N^(n) → evaluate TYG^(n)
2. Solve linear LSQ: a^(n+1) = argmin_a ||W (M_ry TYG^(n) a + M_re Ie - Y)||²
3. Solve GS: Δ*ψ^(n+1) = -μ₀ R p'(ψ_N^(n+1)) - TT'(ψ_N^(n+1))/R
4. Check ||ψ^(n+1) - ψ^(n)|| / ||ψ^(n)|| < psichco

This converges because the solution is typically close to linear in profile
coefficients for fixed boundary, and the boundary moves slowly between iterations.
