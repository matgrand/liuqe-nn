# FBT Static Free-Boundary Equilibrium Solver — Technical Reference

This document is a detailed technical reference for the FBT (Free Boundary Tokamak) solver as
implemented in the MEQ (MatlabEQuilibrium) toolbox at the Swiss Plasma Center, EPFL. The target
audience is someone who understands tokamak MHD equilibrium physics but is not yet familiar with
this particular code.

---

## 1. Overview of Static Free-Boundary Equilibrium

### 1.1 The Physical Problem

A tokamak plasma equilibrium is described by the Grad-Shafranov (GS) equation, which in axisymmetric
$(R, Z)$ cylindrical coordinates takes the form:

$$\Delta_F^* \psi = -\mu_0 R^2 p'(\psi) - \frac{\mu_0}{4\pi^2} TT'(\psi) = -\mu_0 R j_\phi(R, \psi)$$

where:
- $\psi(R, Z)$ is the poloidal flux function (Wb/rad) at grid point $(R, Z)$
- $\Delta_F^*$ is the Grad-Shafranov differential operator: $\Delta_F^* \psi = R \frac{\partial}{\partial R}\left(\frac{1}{R} \frac{\partial \psi}{\partial R}\right) + \frac{\partial^2 \psi}{\partial Z^2}$
- $p'(\psi)$ is the derivative of the kinetic pressure with respect to $\psi$
- $T(\psi) = R \cdot B_\phi$ is the poloidal current function (diamagnetic function), and $TT'(\psi)$ is its derivative
- $j_\phi$ is the toroidal current density

In a **free-boundary** problem the boundary of the plasma (the Last Closed Flux Surface, LCFS) is
not known in advance; it is part of the solution. The external conductors (coils, vessel) impose
boundary conditions on `ψ` at the computational grid edge through their Green's function
contributions.

**FBT solves the following forward problem:**

Given:
- Active coil currents `Ia` (or target coil currents with associated weights)
- Passive conductor currents `Iu` (prescribed or self-consistently evolved)
- Current profile shape functions parameterized by `bffbt` (or other basis function sets), and
  scalar constraints such as total plasma current `Ip`, safety factor on axis `qA`, poloidal
  beta `βp`
- Geometric constraints: set of control points where flux, field, or field gradients should take
  specified values

FBT finds:
- Self-consistent `ψ(R, Z)` on the computational grid
- Plasma current distribution `Iy(R, Z)` satisfying the GS equation
- Coil currents `Ia` that satisfy the cost function / constraints
- An offset flux `Fb` (boundary flux fitting parameter)

### 1.2 Why Free Boundary

"Free boundary" means the LCFS location — whether the plasma is limited by the wall or has an
X-point (separatrix) — is determined by the solution itself rather than being an imposed boundary.
At each iteration FBT re-identifies which grid points are inside the plasma domain (via `meqpdom`)
and updates the current distribution accordingly. The coil geometry and mutual inductances provide
the connection between external currents and the flux inside the vessel.

This contrasts with fixed-boundary codes (e.g., CHEASE), where the plasma boundary is given as
input.

### 1.3 Comparison with LIUQE (Inverse vs. Forward Problem)

LIUQE is the TCV equilibrium reconstruction code (inverse problem): given magnetic measurements
(flux loops, Rogowski coils, B probes) it finds the plasma current profile and `ψ(R, Z)`. LIUQE
"inverts" the GS equation using experimental data as constraints.

FBT is a **forward problem**: given coil currents and a current profile parametrization, it
computes the resulting equilibrium. FBT is used for **shot design** — computing what coil current
waveforms will produce a desired plasma shape, rather than reconstructing what happened.

Both codes share the same MEQ framework for the GS solver, base functions, and output structures.

---

## 2. FBT Algorithm — Top-Level Call Chain

### 2.1 Entry Points

**Primary entry point:** `fbt.m`
```matlab
[L, LX, LY] = fbt('TCV', shot, t, 'PAR', val, ...)
```

**TCV-specific entry point:** `fbte.m`
```matlab
[L, LX, LY] = fbte(shot, varargin)
```
which simply calls `fbt('TCV', shot, [], 'LYall', true, varargin{:})` and then writes results to
the TCV MDS+ tree.

### 2.2 Initialization Sequence in `fbt.m`

```matlab
fbtptok = str2func(['fbtp' tok]);    % e.g., fbtptcv
fbtgtok = str2func(['fbtg' tok]);    % e.g., fbtgtcv
fbtxtok = str2func(['fbtx' tok]);    % e.g., fbtxtcv

P = fbtptok(shot, varargin{:});      % Load configuration parameters
G = fbtgtok(shot, P);                % Load geometry / electromagnetic data
P = fbtp(P);                         % Fill in MEQ-level defaults
G = fbtg(G, P);                      % Fill in MEQ geometry defaults
L = fbtc(P, G);                      % Assemble ancillary data (L structure)

LX = fbtxtok(t, L);                  % Load time-dependent input data
LX = fbtx(L, LX);                   % Validate and fill LX defaults

LY = fbtt(L, LX);                   % THE CORE SOLVER
```

The three output structures are:
- `L`: precomputed static data (geometry, matrices, configuration) — reusable across runs
- `LX`: time-dependent input data (coil targets, shape constraints, profile parameters)
- `LY`: time-dependent output data (equilibrium solutions)

### 2.3 The Core Solver: `fbtt.m`

`fbtt.m` loops over time slices (or intervals for the circuit-equation case) and for each calls
the optimization engine. The key code structure is:

```matlab
for kstg = 1:nstage
  % Assemble cost function matrices
  [s, lin, IDs, ws, bbs] = fbtlsqinit(L, LXs, s, fbtsolver);
  % Compute initial guess
  s = fbtinit(L, LXs, s, bbs, fbtsolver);
  % Run optimization (Picard / SQP / Newton)
  [s, LYs, conv, nFeval] = meqopt(L, LXs, LYp, s, ...);
end
```

### 2.4 The Fixed-Point (Picard) Iteration Structure

FBT is fundamentally an iterative solver. The nonlinearity arises because the plasma current
`Iy(R, Z)` depends on `ψ(R, Z)` (through the plasma domain determination and profile functions),
but `ψ(R, Z)` itself depends on `Iy` (through the GS equation).

**Picard scheme (default for `algoNL = 'all-nl'`):**

1. Start with an initial current distribution `Iy⁰` (parabolic inside a bounding box)
2. Solve the linear least-squares problem: find `Ia` and `Fb` that satisfy the cost function
   given the current plasma contribution `Iy`
3. Compute the new flux: `ψ = MYY·Iy + MYE·Ie + Fb`
4. Find the plasma domain (LCFS, axis) from `ψ`
5. Update `Iy` from basis functions evaluated on the new `ψ`
6. Check convergence on the change in `ψ` or `Iy`; if not converged, go to step 2

**The SQP/Newton scheme (enabled with `useSQP = true`):**

The full non-linear problem is solved in one shot using sequential quadratic programming. The
Jacobian of the GS residual with respect to `(Iy, Ia, Fb)` is used to form Newton steps. The
parameter `wsSQP` (default 2) sets the number of Picard warm-start iterations before activating
the full Newton solver.

### 2.5 How the GS Equation is Solved (FGS/FGE Operator)

The discrete GS equation is represented through the operator `meqFx(L, Iy, Ie)`:

$$F_x = G \cdot I_y + M_{xa} \cdot I_a + M_{xu} \cdot I_u + F_b$$

where:
- $F_x$ is the flux vector on the interior computational grid ($n_{zy} \times n_{ry}$ points)
- $G$ is the Green's function matrix connecting plasma current to interior flux (mutual inductance)
- $M_{xa}$ (`G.Mxa`) connects active coil currents to interior flux
- $M_{xu}$ (`G.Mxu`) connects passive/vessel currents to interior flux
- $F_b$ is a uniform offset (the "boundary flux" fitting parameter)

The GS operator enforces that `Fx` satisfies the GS equation given `Iy`. The modes `gsxe=1,2,3`
select different inversion strategies for this operator (direct inversion, iterative, or hybrid).

### 2.6 Plasma Domain Identification

At each iteration, `meqpdom` identifies:
- `FA` = flux at the magnetic axis (maximum/minimum of `Fx`, depending on sign of `Ip`)
- `FB` = flux at the LCFS (either limited by the wall or at an X-point saddle)
- `Opy` = boolean mask of which grid points are inside the plasma
- `rX, zX, FX` = X-point location and flux (if diverted)

The normalized flux is then `FN = (ψ - FA)/(FB - FA)`, which ranges from 0 at the axis to 1
at the boundary.

The criterion for the LCFS is:
- **Limiter**: LCFS is where `ψ` first touches the first-wall contour `(rl, zl)`
- **Diverted (X-point)**: LCFS passes through a magnetic saddle point where
  `Br = -1/R ∂ψ/∂Z = 0` and `Bz = 1/R ∂ψ/∂R = 0`

### 2.7 Convergence Criteria

From `fbtp.m`:
- `tol = 1e-6`: tolerance on the RMS residual of the GS equation
- `tolstep = 1e-8`: tolerance on the step size in the Newton iteration
- `niter = 100`: maximum number of Picard iterations

Convergence is declared when the weighted cost function residual falls below `tol`. The output
`LY.chi` is the RMS cost function residual and `LY.isconverged` flags whether convergence was
achieved.

---

## 3. Current Profile Parametrization and Base Functions

### 3.1 The Role of Base Functions

The plasma current density at any grid point (R, Z) with normalized flux `FN` is:

```
j_φ(R, Z) = Σᵢ agᵢ · gᵢ'(FN) · (R · fPgᵢ + (1/R) · fTgᵢ)
```

where:
- `ag` = coefficient vector (ng × 1), the unknowns of the profile
- `gᵢ'(FN)` = the i-th basis function evaluated at `FN`
- `fPg(i)` = 1 if basis function i contributes to `p'` (pressure gradient), 0 otherwise
- `fTg(i)` = 1 if basis function i contributes to `TT'` (diamagnetic function gradient), 0 otherwise

This parametrization follows directly from the GS source term:
```
μ₀ j_φ = μ₀ R p'(ψ) + μ₀/(4π²R) TT'(ψ)
```

In matrix form over all interior grid points:
```
Iy = TYG · ag
```
where `TYG(ny, ng)` is assembled by the base function called in mode 1.

The scalar constraints on the profile (total current `Ip`, poloidal beta `βp`, on-axis q `qA`) are
then linear in `ag` and are assembled as rows of the equality constraint matrix in the optimization.

### 3.2 The `ag` Constraints (Profile Integrals)

The number of basis functions `ng` equals the number of profile constraints `nC`. These constraints
are assembled in `meqagconc` from the list `fbtagcon` (default: `{'Ip', 'Wk', 'qA'}` or for TCV
`{'Ip', 'qA', 'fbtlegacy'}`).

The constraints take the form of linear equations in `ag`:
```
agconc · ag = [Ip, Wk, qA, ...]ᵀ
```

where `agconc` is the matrix of integral operators (rows) for each constraint, and the right-hand
side values come from `LX`.

Common constraints:
- `Ip`: total plasma current = `(2π/μ₀) · ITPDG · ag` where `ITPDG` sums basis function integrals
  over the plasma domain
- `Wk`: kinetic energy = `(3/2) · μ₀ · p integral`
- `qA`: safety factor on axis, expressed via the GS equation at the axis
- `li` (or `liD`): normalized internal inductance

### 3.3 The `bffbt` Basis Function (TCV Default)

The default TCV basis function set is `bffbt`, defined in `bffbt.m`. It uses 3 basis functions
parameterized by a 6-element vector `P = [pf, pl, pe, tf, tl, te]`:

**For p' (2 functions):**
```
g₁(FN) = pf · (pl + (1-pl)·FN) · (1-FN)^pe
g₂'    = 0   (fPg = [1; 0; 0])
```

**For TT' (2 functions, contributing via fTg):**
```
g₂(FN) = tf · (1-FN)^te
g₃(FN) = tf · (tl·FN - 1 - tl) · FN · (1-FN)^te  [Note: -(tl·f+1)·F·x₂ where f=1-FN]
```

The masks are:
```
fPg = [1; 0; 0]   (g₁ contributes to p')
fTg = [0; 1; 1]   (g₂, g₃ contribute to TT')
```

So ng=3 with 1 p' basis function and 2 TT' basis functions.

The primitive integrals `IgQ` (for p, T²/2) are computed analytically:
```
IgQ₁ = ((-(pl+1/(pe+1))/(pe+2)) + ((pl-1)/(pe+2))·FN) · pf·(1-FN)^(pe+1)·(1-FN)
IgQ₂ = (-1/(te+1)) · tf·(1-FN)^(te+1)·(1-FN)
IgQ₃ = (...) · tf·(1-FN)^(te+1)·(1-FN)
```

The parameters have the following physical meaning:
- `pf`: amplitude of p' profile
- `pl`: peaking/hollowing parameter (pl=1 gives peaked at axis, pl<1 more hollow)
- `pe`: edge exponent for p' (steepness near LCFS)
- `tf`: amplitude of TT' profile (controls toroidal field variation and diamagnetic current)
- `tl`: controls the shape of the second TT' basis function
- `te`: edge exponent for TT'

The `bfpfbt.m` function computes the time-dependent bfp parameters from MGAMS control parameters:
```
bfp = [ppfac; ppal; ell1; ttfac1·rh²; omega; emm1]
```
where `rh` is the radial center of the bounding box of plasma boundary control points.

### 3.4 Full Catalog of Base Function Types

#### 3.4.1 `bfabmex` — Polynomial Basis (MEX)

Parameters: `P = [nP, nT]` with `nP ≤ 3`, `nT ≤ 3`.

Produces `nP + nT` polynomial basis functions. The polynomials are:
```
For p' (nP functions):
  g₁ = F - FB
  g₂ = (F - FB)(F - FA)
  g₃ = (F - FB)(F - FA)(2F - FB - FA)

For TT' (nT functions): same polynomials
```

All these vanish at the LCFS (`F = FB`) so boundary conditions are automatically satisfied.
Physical meaning: simple low-order polynomial profiles that can represent a range of shapes from
peaked (dominant g₁) to hollow or edge-peaked (g₂, g₃).

MEX implementation in `libmeq/bfab.c`.

#### 3.4.2 `bf3pmex` — 3-Function Polynomial (MEX, special case)

A particular case of `bfabmex` with:
- `p' = 0` (zero pressure, vacuum-like)
- 1 function for TT': `F - FB`
- OR: the standard 3-function polynomial (1 for p', 2 for TT')

Controlled by boolean parameter `bfp`:
- `bfp = false`: standard 3-function
- `bfp = true`: only TT', no p'

Physical meaning: minimal parametrization for testing or forcefree plasmas.

#### 3.4.3 `bf3imex` / `bfli` — Tabular Interpolation

Parameters: `P.bfp = struct('gNg', GN, 'IgNg', IGN, 'fPg', FP, 'fTg', FT)` where `GN(n, ng)` 
contains the basis function values tabulated on `linspace(0, 1, n)` and `IGN` their integrals.

`bf3imex` is a MEX implementation; `bfli` is the pure MATLAB equivalent. Both perform linear
interpolation between tabulated values.

Physical meaning: completely arbitrary user-defined profiles, provided as tables. Useful when
profiles come from external computations or experimental fits.

The integrals `IGN` can be precomputed from `GN` using `bfprmex(GN)`.

#### 3.4.4 `bfsp` — B-Spline Basis

Parameters: `P` structure with `P.tauP`, `P.MP` for p' basis and `P.tauT`, `P.MT` for TT',
obtained from `bfpsp`.

Uses order-5 B-splines (degree 4) for the integrated profiles, with their exact derivatives for
p' and TT'. The knot distribution is the same for p' and TT' by default.

Physical meaning: smooth, flexible profiles with many degrees of freedom. The number of basis
functions is determined by the number of knot intervals. More flexible than polynomials but
requires more profile constraints.

Note: no regularization implemented yet (as of 2025).

#### 3.4.5 `bfgenD` — Multi-Domain Generic Wrapper

Parameters: `P` = cell array, each row is `{bfct_handle, bfp, domain_assignment_vector}`.

Example:
```matlab
bfp = {@bf3pmex, false, [1;0;0];   % active only in domain 1
       @bf3pmex, false, [0;1;0];   % active only in domain 2  
       @bfabmex, [1,2], [1;1;1]};  % active in all domains
```

`bfgenD` assembles multiple single-domain basis function sets into one multi-domain set, handling
the mapping of normalized flux between domains. The `TDg` matrix (ng × nD) indicates which domains
each basis function operates in.

This is the primary tool for multi-domain (doublet, droplet, complex topology) plasmas.

#### 3.4.6 `bfdoublet` — Doublet-Optimized 3-Domain Wrapper

A specialization of `bfgenD` for the specific case of a doublet (two plasma lobes connected by a
mantle). It handles 3 domains:
- Domain 1: upper lobe
- Domain 2: lower lobe
- Domain 3: mantle (the common region around the internal X-point)

The mantle basis functions are "trivially continued" into the lobe regions: the derivative (p', TT')
is extended as a constant equal to the mantle axis value, and the integral (p, T²/2) is extended
as a linear ramp. This ensures physical continuity of pressure and toroidal field across the
domain boundaries.

See Section 7 for more details on doublet physics.

### 3.5 The `bfprmex` Function

`bfprmex(GN)` numerically integrates tabulated basis functions to obtain primitives. It is
specifically designed for use with `bf3imex`/`bfli`: given `GN(n, ng)` (the basis functions
evaluated on `linspace(0, 1, n)`), it returns `IGN(n, ng)` satisfying `IGN = 0` at `FN = 1`
(LCFS boundary condition).

Internally it uses trapezoidal integration with the constraint that the integral is zero on the LCFS.

### 3.6 Mode System for Base Function Calls

All base function sets share a common calling convention via `bfct` (the dispatcher). The `mode`
argument selects what is computed:

| Mode | Returns | Used for |
|------|---------|----------|
| 0 | fPg, fTg, TDg | Structural info (which functions contribute to p' vs TT', domain assignments) |
| 1 | TYG, TPDG, ITPDG | Full plasma current matrix on grid (used in GS iteration) |
| 2 | GN, IGN | Normalized basis function values at prescribed FN points |
| 3 | APP, ATTP, AP, AHQT | Profile coefficients (p', TT', p, T²/2) from ag |
| 4 | TDG, TGY | Chord integrals for diagnostic |
| 5 | G0G, IG0G | Values at magnetic axis |
| 6 | QQG, XQ | Regularization (smoothness) constraints |
| 7 | QCG, XC | Inequality constraints on profiles |
| 8 | BTY | Toroidal field profile |
| 11 | dTYG/d(F,F0,F1) | Jacobians for Newton solver |
| 12 | dGN/d(FN,F0,F1) | Jacobians of normalized basis functions |
| 91 | G, IG | Physical (unnormalized) basis functions |
| 92 | dG/dF, dIG/dF | Physical basis function derivatives |

---

## 4. Geometry Setup

### 4.1 The L, G, P Structure Hierarchy

The code uses a three-level hierarchy:
- `P` (parameters): configuration choices, algorithm settings, coil selections
- `G` (geometry): electromagnetic matrices, coil positions, grid coordinates
- `L` (library): combined working structure with all precomputed data

### 4.2 TCV Geometry (`fbtgtcv.m`, `meqgtcv`)

`fbtgtcv` calls `meqgtcv` with the selected active coils (`P.sela`), passive structures
(`P.selu`), and grid (`P.selx`). The TCV selections are:

```matlab
P.sela  = {'E', 'F'};   % 8 E-coils + 8 F-coils = 16 active coils
P.selu  = 'n';           % No vessel passive structures (for static FBT)
P.selx  = 'XF';         % Historical FBTE computational grid
```

TCV-specific constraint: `dohzero = true` enforces `IOH1 = IOH2` (if OH coils are included) via
the `G.Coa` matrix.

### 4.3 Electromagnetic Matrices (from `meqg`)

The key matrices precomputed at initialization are:

| Symbol | Shape | Meaning |
|--------|-------|---------|
| `G.Mxa` | (ny, na) | Green's function: flux at interior grid from active coils |
| `G.Mxu` | (ny, nu) | Green's function: flux at interior grid from passive conductors |
| `G.Maa` | (na, na) | Mutual inductance between active coils |
| `G.Mau` | (na, nu) | Mutual inductance between active and passive |
| `G.Muu` | (nu, nu) | Mutual inductance between passive conductors |
| `G.Ra`  | (na,)   | Active coil resistances |
| `G.Ru`  | (nu,)   | Passive conductor resistances |
| `G.Dda` | (na-1, na) | Finite difference matrix for dipole constraints |
| `G.Coa` | (no, na) | Linear combination matrix for coil constraints |

The mutual inductance matrices are computed using Green's functions for axisymmetric conductors.

### 4.4 Computational Grid

The computational grid uses `rx(nrx)` and `zx(nzx)` as the discrete R, Z coordinates. Interior
grid points are `ry = rx(2:end-1)`, `zy = zx(2:end-1)` with `ny = (nrx-2)*(nzx-2)` total
interior points.

The grid choice `P.selx = 'XF'` selects the historical TCV FBTE grid. Alternative `P.selx = 'X'`
selects a slightly different resolution.

Within `L` (from `meqg`/`fbtc`):
- `L.ry`, `L.iry` = R values of interior grid points and their reciprocals
- `L.rry`, `L.zzy` = matrices of R, Z for all interior grid points
- `L.Mxy` = interpolated Green's function matrix for control points
- `L.Tye` = external current contribution to interior flux (combined active + passive)

### 4.5 Coil Limit Matrices (from `fbtptcv`)

TCV coil current limits from the MDS+ tree:
```matlab
P.limc = [eye(na); Ceq]     % limc: constraint matrix (identity for individual limits + Ceq for protection)
P.limu = [INOM_A; IHIG]     % upper limits [A]
P.liml = [-INOM_A; ILOW]    % lower limits [A]
```
where `INOM_A` is the nominal coil current, and `Ceq`, `IHIG`, `ILOW` are from the coil protection
system.

---

## 5. The FBT Cost Function and Optimization

### 5.1 Problem Formulation

FBT solves a constrained optimization problem at each (nonlinear) iteration. The optimization
variables are:
- `Ia`: active coil currents (na × 1), in static case
- `Va`: active coil voltages (na × 1), in circuit-equation case
- `Fb`: boundary flux offset scalar (1 × 1)
- (In evolutive case also `Iu` and `IniD`, the passive conductor currents and non-inductive
  plasma current drive)

The state `x` contains the plasma current distribution `Iy` (or flux `Fx` depending on `algoNL`).

**Cost function (least squares):**
```
min W = (1/2) || R_a·x + R_a·Ie + R_a·Fb - Ra0 ||²_w
```
subject to:
```
Re·x + Re·Ie + Re·Fb = Re0     (equality constraints)
Ri·x + Ri·Ie + Ri·Fb ≤ Ri0    (inequality constraints: coil limits)
```

The subscript `a` = "approximate" (soft constraints in the cost function), `e` = "equality", 
`i` = "inequality". The subscripts `Iy`/`Ie`/`Fb` label which variable each matrix multiplies.

### 5.2 Cost Function Terms (`fbtcost.m`)

From `fbthelp.m`, the full cost function (before splitting into hard vs. soft) is:

```
1/(gpid·gpie) · [ Ia                    = gpia                           ]  ID=1  Active coil currents
1/(gpod·gpoe) · [ Coa·Ia                = gpoa                           ]  ID=1  Coil combinations
1/(gpdd·gpde) · [ Dda·(gpdw⊙Ia)        = gpda                           ]  ID=2  Dipoles (current differences)
1/(gpud·gpue) · [ Iu                    = gpua                           ]  ID=3  Passive currents
1/(gpfd·gpfe) · [ MSe·Ie - gpfb·Fb      = gpfa - gpfy·MSy·Iy             ]  ID=4  Flux at control points
1/(gpbd·gpbe) · [ BrSe·Ie               = gpbr - gpby·BrSy·Iy            ]  ID=5  Br at control points
1/(gpbd·gpbe) · [ BzSe·Ie               = gpbz - gpby·BzSy·Iy            ]  ID=5  Bz at control points
1/(gpcd·gpce) · [ MrrSe·Ie              = gpcr - gpcy·MrrSy·Iy           ]  ID=6  ∂²ψ/∂R² (flux Hessian)
1/(gpcd·gpce) · [ MrzSe·Ie              = gpcz - gpcy·MrzSy·Iy           ]  ID=6  ∂²ψ/∂R∂Z
1/(gpvd·gpve) · [ BrrSe·Ie              = gpvrr - gpvy·BrrSy·Iy          ]  ID=7  Vacuum ∂Br/∂R
1/(gpvd·gpve) · [ BrzSe·Ie              = gpvrz                          ]  ID=7  Vacuum ∂Br/∂Z
1/(gpvd·gpve) · [ BzzSe·Ie              = gpvzz                          ]  ID=7  Vacuum ∂Bz/∂Z
1/(gpad·gpae) · [ Va                    = gpaa                           ]  ID=8  Active coil voltage (circuit case)
```

Where:
- `MSy·Iy` is the flux at the control point due to the plasma current (`M = Green's function`)
- `MSe·Ie` is the flux due to external currents
- `gpf*`, `gpb*`, `gpc*`, `gpv*`, `gpi*` etc. are the control point constraint parameters
- `gp?e`: relative weight/error for each equation (Inf means disable, 0 means equality constraint)
- `gp?d`: global scale for the error

### 5.3 How Constraints are Split

In `fbtcost.m`, lines in the combined cost matrix are classified:
```matlab
ke = ~isfinite(err) | isnan(Ra0);  % Remove: infinite weight or NaN value
kC = ~ke & (err == 0);             % Convert to equality: zero error
```

Lines with `err = 0` (i.e., `gp?e = 0` at that row) become hard equality constraints. Lines
with finite non-zero `err` enter the weighted least-squares cost function. Lines with `err = Inf`
or `NaN` are removed.

### 5.4 Inequality Constraints (Coil Current Limits)

From `fbtcost.m`:
```matlab
limu = P.limm .* P.limu;    % upper limits with margin factor
liml = P.limm .* P.liml;    % lower limits
% Encode as: limc·Ia ≤ limu  and  -limc·Ia ≤ -liml
RiIa = [P.limc(~j,:); -P.limc(~j,:)];
Ri0  = [limu(~j);     -liml(~j)];
```

Coil pairs with `limu == liml` are promoted to equality constraints.

### 5.5 Optimization Algorithms

The optimization is handled by `meqopt` (called from `fbtt.m`). The available algorithms
(controlled by `P.useSQP` and `P.solver`):

**Picard (default for TCV FBT, `useSQP = false`):**
- Iterates solving the linear least-squares problem with fixed plasma state
- One Newton step at initialization to get reasonable coil currents
- Then Picard iterations updating the plasma state

**SQP (Sequential Quadratic Programming, `useSQP = true`):**
- Full nonlinear solver treating the combined (plasma + coil) problem
- Uses `P.wsSQP` Picard warm-start iterations before Newton
- Supports `solver = 'SQP'` or `solver = 'SQP-broyden'`
- Uses Block-GMRES or direct inversion for the linear subproblem

**Newton-GS / all-nl-Fx:**
- Alternative linearization where the state is `Fx` (flux) instead of `Iy` (current)
- Enabled by `P.algoNL = 'all-nl-Fx'` or `'Newton-GS'`

### 5.6 The Initial Guess (`fbtinit.m`)

The initial guess generation:
1. If `LX.Iy` is provided, use it directly as the initial plasma current distribution
2. Otherwise, use `meqinitIy`: generate a parabolic (with peaking `capaj`) current distribution
   inside the bounding box `bb` of control points with `gpb = 1`
3. Perform one Newton step (linear solve) to get initial coil currents matching the plasma

The bounding box `bb` is computed from boundary control points (`gpb = 1`) and any X-point
constraints, after eliminating strike points from the initial guess (in `fbtptcv.m`, function
`gpb`).

---

## 6. Input/Output Data Structures

### 6.1 LX: Input Data to the Solver

The `LX` structure contains all time-dependent inputs. Each field has shape `(nrows, nt)` where
`nt = numel(LX.t)`.

**Scalar/vector quantities (per time slice):**
- `LX.t`: time array (1 × nt)
- `LX.Ip`: total plasma current [A] (1 × nt) [or `LX.IpD` for multi-domain (nD × nt)]
- `LX.rBt`: R·Bt product at reference radius [T·m] (1 × nt)
- `LX.qA`: on-axis safety factor (nD × nt)
- `LX.bp` or `LX.bpD`: poloidal beta (1 × nt or nD × nt) — optional if `Wk` or `WkD` provided
- `LX.niter`: max iterations (1 × nt)
- `LX.tol`: convergence tolerance (1 × nt)
- `LX.capaj`: current peaking parameter for initial guess (1 × nt)
- `LX.idoublet`: flag for doublet topology (1 × nt)
- `LX.isaddl`: flag for additional (multi-domain) plasma (1 × nt)

**Control point specifications (npts × nt):**
- `LX.gpr`, `LX.gpz`: R, Z coordinates of control points
- `LX.gpb`: whether point is on LCFS (used for initial guess bounding box)
- `LX.gpfa`, `LX.gpfb`, `LX.gpfe`: flux constraint: value, boundary-relative flag, weight
- `LX.gpbr`, `LX.gpbz`, `LX.gpba`, `LX.gpbe`: Br, Bz, field angle, weight
- `LX.gpcr`, `LX.gpcz`, `LX.gpca`, `LX.gpce`: flux Hessian ∂²ψ/∂R², ∂²ψ/∂R∂Z, angle, weight
- `LX.gpvrr`, `LX.gpvrz`, `LX.gpvzz`, `LX.gpve`: vacuum field gradients and weight
- `LX.gpfd`, `LX.gpbd`, `LX.gpcd`, `LX.gpvd`: global scale factors for constraint groups

**Current targets (na × nt or nu × nt):**
- `LX.gpia`, `LX.gpie`, `LX.gpid`: target active coil current, weight, scale
- `LX.gpua`, `LX.gpue`, `LX.gpud`: target passive current, weight, scale
- `LX.gpda`, `LX.gpde`, `LX.gpdd`, `LX.gpdw`: dipole target, weight, scale, individual weights
- `LX.gpoa`, `LX.gpoe`, `LX.gpod`: coil combination constraint target/weight/scale

**Profile basis function parameters:**
- `LX.bfp`: time-varying basis function parameters (6 × nt for `bffbt`)

**Time-derivative constraints (for circuit case):**
- `LX.g1r`, `LX.g1z`, `LX.g1b`, `LX.g1f*`: first-order time derivative constraints
- `LX.g2r`, `LX.g2z`, `LX.g2b`, `LX.g2f*`: second-order time derivative constraints

### 6.2 LY: Solver Output

Standard MEQ outputs (from `meqt`):
- `LY.t`: time array
- `LY.aq`: profile coefficients (ag) vector (ng × nt)
- `LY.aW`: profile coefficients for the weighted (p', TT') representation
- `LY.Ia`: active coil currents found (na × nt)
- `LY.Iu`: passive conductor currents (nu × nt)
- `LY.Iy`: plasma current distribution (nzy × nry × nt) [only if `argout` requested]
- `LY.Fx` or `LY.Fy`: poloidal flux on grid (nzy × nry × nt)
- `LY.rA`, `LY.zA`: magnetic axis location (nD × nt)
- `LY.FA`, `LY.FB`: axis and LCFS flux (nD × nt)
- `LY.rX`, `LY.zX`, `LY.FX`: X-point location and flux
- `LY.rB`, `LY.zB`: LCFS touchdown (limiter/X-point)
- `LY.lX`, `LY.lB`: flags for X-point/limiter presence
- `LY.Opy`: plasma domain mask (nzy × nry × nt)
- `LY.Ip`, `LY.IpD`: total plasma current(s)
- `LY.bp`, `LY.bpD`: poloidal beta
- `LY.qA`: on-axis safety factor
- `LY.li`: normalized internal inductance
- `LY.isconverged`: convergence flag

FBT-specific outputs:
- `LY.Sp`: sum of plasma-containing cell areas [m²]
- `LY.Ep`: plasma elongation (height/width of bounding box)
- `LY.Fb`: boundary flux fitting parameter
- `LY.chi`: RMS cost function residual
- `LY.chie`: RMS constraint residual
- `LY.werr`: weighted cost function error
- `LY.uerr`: unweighted cost function error
- `LY.ID`: equation ID for each cost function term

### 6.3 `fbtxcat` and `fbtxinterp`: Time Series Handling

**`fbtxcat`**: concatenates multiple LX structures in time. After calling `meqcat`, it fixes NaN
values in weight fields (`gp?e`, `gp?d`) to `Inf` (meaning "inactive") rather than NaN.

**`fbtxinterp`**: interpolates LX in time using methods `'linear'`, `'nearest'`, `'previous'`,
`'next'`. Linear interpolation of control points is only allowed when the `gpr`/`gpz` patterns
and weight structures are compatible between adjacent time slices. For interpolated time points
without original control points, constraints are disabled (coordinates set to NaN, weights set
to Inf).

### 6.4 `fbt2pcs` — Converting FBT Results to MGAMS Parameters

`fbt2pcs(L, LX)` converts the FBT input structure `LX` into a MGAMS (TCV plasma control system)
parameter object `MG`. This is used for shot preparation: the engineer designs an equilibrium in
FBT, then writes the resulting parameters to the PCS tree for experimental use.

The conversion renames fields appropriately:
- `LX.tol` → `MG.testa`
- `LX.niter` → `MG.itamax`
- `LX.rBt` → `MG.bzero = rBt/0.88` (R₀·B₀ at R₀=0.88m)
- `LX.Ip` → `MG.placu1 = Ip/1e3` (kA)
- Basis function parameters serialized to `debug_p1` string

### 6.5 `fbt2SC` — Extracting Shape Control Descriptors

`fbt2SC(L, LX)` classifies control points by type for use in the TCV shape control system:
- `lcC = 1`: general contour points (LCFS boundary shape)
- `lcL = 1`: limiter touch point (plasma is limited)
- `lcS = 1`: strike points (on divertor target)
- `lcX1 = 1`: primary X-point (`Br = Bz = 0`, on LCFS)
- `lcX2 = 1`: secondary X-point (off-LCFS saddle)
- `lcD = 1`: separatrix clearance point (outermost contour point when secondary X exists)
- `lcI = 1`: ignored/inactive points

Classification logic:
- Limiter: closest control point to the wall, when no primary X-point found
- Primary X-point: points with `gpbr = gpbz = 0` and `gpb = 1`
- Strike points: points within 2 cm of the wall that are not on the LCFS

---

## 7. Doublet Plasma Support

### 7.1 What is a Doublet Plasma

A doublet (or double-droplet) plasma consists of two separate plasma lobes connected via an
internal X-point. The configuration has three topological domains:
- Domain 1: upper plasma lobe
- Domain 2: lower plasma lobe  
- Domain 3: the "mantle" region surrounding the internal X-point, between the two lobes

In TCV, doublets are enabled by setting `P.mvloop = 8` or `P.idoublet = true`.

### 7.2 Multi-Domain Infrastructure

The multi-domain support is built into the MEQ base through the `isaddl` flag (additional domain)
and `nD` (number of domains). For doublets, `nD = 3` (or 2 if no mantle current).

The GS equation must be solved simultaneously in all domains. Flux values `FA(1:nD)` and
`FB(1:nD)` are determined per domain. The normalized flux `FN` is computed separately for each
domain.

### 7.3 The `bfdoublet` Basis Function

`bfdoublet` handles the physically-motivated profile continuation:

**Mantle basis functions (domain 3):** normal basis function evaluation in the mantle region
`[F0(3), F1(3)]`.

**Continuation into lobes (domains 1, 2):** For lobe regions where `ψ < F0(3)` (inside the
internal separatrix):
- Derivative functions (p', TT') are extended as a **constant** equal to their mantle axis value
  `GA3 = g'(F0(3))`
- Integral functions (p, T²/2) are extended as a **linear ramp**: `IG(F) = IG(F0(3)) + (F - F0(3))·G(F0(3))`

This ensures:
1. The current density is constant inside the lobes (from the mantle profile perspective)
2. The pressure integral and toroidal field are continuous across the internal X-point
3. The constraint equations for integrated quantities remain valid

From `bfdoublet.m`, mode 1 (grid current computation):
```matlab
Tyg(:, igm) = Tyg(:, igm) + (fac1 + fac2) .* GA3;  % Add mantle contribution to lobes
```
where `fac1`, `fac2` are the R (for p') and 1/R (for TT') weighted sums over lobe grid points.

### 7.4 Input Specifications for Doublets

In LX, doublet-specific fields:
- `LX.IpD`: current per domain (3 × nt for 3-domain doublet)
- `LX.bpD`: poloidal beta per domain (3 × nt)
- `LX.qA`: on-axis q per lobe (2 × nt, since domain 3/mantle has no axis)

The constraint `fbtagcon` for multi-domain can be a cell array of cell arrays, one per domain:
```matlab
P.fbtagcon = {{'Ip','Wk','qA'}, {'Ip','Wk','qA'}, {'Ip','Wk'}};
```

In `fbtptcv.m`:
```matlab
if P.idoublet && ~any(strcmp(varargin(1:2:end),'isaddl'))
  P.isaddl = 1;
end
```

---

## 8. Time Derivative Computation (`fbttimeder.m`)

### 8.1 Purpose

`fbttimeder` computes finite-difference stencil matrices `d1dt` and `d2dt` for first and
second-order time derivatives on an arbitrary (non-uniform) time grid.

### 8.2 Stencil Types

For a time sequence `t(1:nt)`, the stencils satisfy:
```
dF/dt|ₜ ≈ sum_j d1dt(i,j) · F(t(j))
d²F/dt²|ₜ ≈ sum_j d2dt(i,j) · F(t(j))
```

First derivative options:
- `'backward'` (default): `d1dt(i,i-1) = -1/Δt`, `d1dt(i,i) = 1/Δt`
  (zero for the first time step, matching a steady-state initial condition)
- `'centered'`: second-order accurate, uses both previous and next points
- `'forward'`: `d1dt(i,i) = -1/Δt`, `d1dt(i,i+1) = 1/Δt` (zero for last step)

Second derivative (always centered, second-order):
```
d2dt(i,i-1) = 2·(1/Δt₋)·(1/Δt₊)·(1/Δt₋)/(1/Δt₋+1/Δt₊)
d2dt(i,i  ) = -2·(1/Δt₋)·(1/Δt₊)
d2dt(i,i+1) = 2·(1/Δt₋)·(1/Δt₊)·(1/Δt₊)/(1/Δt₋+1/Δt₊)
```

### 8.3 Where Time Derivatives Are Used

In the circuit-equation mode (`P.circuit = true`), FBT solves the GSPD (Grad-Shafranov with
Passive Dissipation) problem:

```
Mee·dIe/dt + Re·Ie + Mey·dIy/dt = [Va; 0]
```

This couples the plasma (GS equation) to the external circuit (Faraday's law). The stencils
`d1dt`, `d2dt` are used to discretize the time derivatives over the interval `kstg`.

In the static case (`P.circuit = false`), the stencils are used only for the optional time-
derivative constraints `g1*` and `g2*` in LX (e.g., constraining dψ/dt at a specific location
to match a target loop voltage).

### 8.4 Use in `fbtlsqinit.m`

```matlab
[d1dt, d2dt] = fbttimeder(LXs.t, L.P.timeder);
```
These are stored in `s.d1dt`, `s.d2dt` and used in `meqlsqprep` to assemble the derivative
contributions into the optimization problem.

---

## 9. Physics Accuracy and Tradeoffs

### 9.1 Approximations in FBT

1. **Axisymmetry**: FBT is 2D axisymmetric; all 3D effects (ELMs, kinks, etc.) are ignored.

2. **Equilibrium**: purely static equilibrium — no flows, no toroidal rotation, no kinetic effects.

3. **Resistive wall (static case)**: in the default static mode (`P.circuit = false`), vessel
   currents are prescribed as known inputs (`LX.gpua`), not self-consistently evolved. The
   actual induced wall currents from `dψ/dt` are ignored unless `P.circuit = true`.

4. **Profile parametrization**: the profiles p'(ψ) and TT'(ψ) are parametrized by a small number
   of basis functions (typically 3 with `bffbt`). This limits profile flexibility.

5. **Discrete GS operator**: the GS equation is solved on a finite rectangular grid with
   boundary conditions from external currents via Green's functions. Grid resolution sets a limit
   on the accuracy of the current distribution near the LCFS and X-points.

6. **Green's function approach**: FBT uses analytical (or tabulated) Green's functions for the
   mutual inductance between filament conductors and grid points. This is exact in the limit of
   infinitely thin conductors but uses finite-area coil models in practice.

### 9.2 Speed/Accuracy Tradeoffs

| Setting | Effect |
|---------|--------|
| `niter = 100`, `tol = 1e-6` | Standard accuracy |
| `tol = 1e-4` | Faster, less accurate |
| `gsxe = 1` | Direct GS inversion (accurate, slower for large grids) |
| `gsxe = 2` | Iterative GS inversion (faster for large grids) |
| `gsxe = 3` | Approximate GS (fastest, less accurate) |
| `useSQP = true` | Newton convergence (fewer iterations, each more expensive) |
| `useSQP = false` | Picard (more iterations, each cheaper) |
| `niter = 3` + `itert` | Warm-start from previous time step (fast for slow time evolution) |

For TCV shot design, `useSQP = false` (Picard) is default and adequate. For the `SQP-broyden`
solver (for research), Newton convergence is available.

### 9.3 The Role of Regularization

Profile regularization (preventing unphysical oscillations in p'(ψ), TT'(ψ)) enters through:

1. **Mode 6 of base functions** (`QQG, XQ`): regularization least-squares equations of the form
   `||QQG·ag||² = Xq` added to the cost function

2. **Dipole terms** (`gpdd`): penalize large differences between adjacent coil currents
   `sum |gpdw(i)·Ia(i) - gpdw(j)·Ia(j)|² / gpdd²`
   This prevents unrealistic coil current patterns.

3. **Direct profile constraints** (via `fbtagcon`): hard constraints on `Ip`, `qA`, `βp` limit
   the profile degrees of freedom.

4. **Number of basis functions**: with only 3 basis functions (`bffbt`), the profile is
   intrinsically smooth.

---

## 10. Notes for Python/NumPy Rewrite

### 10.1 Iterative Solver

The Picard iteration loop (currently in `fbtt.m` → `meqopt.m`) can be implemented with:
- A Python function `solve_equilibrium(Ia, ag, Fb)` → `(ψ, Opy, Ip, qA)`
- An outer loop updating `Ia, ag, Fb` using scipy's `lsq_linear` or `quadprog` for the QP step
- Convergence checked on `||Ia - Ia_prev|| / ||Ia||` or on `||ψ - ψ_prev|| / ||ψ_prev||`

For the Newton/SQP variant, the Jacobian structure is critical: the Jacobian of the GS operator
with respect to `Iy` is dense (size `ny × ny`) but the response matrix `dψ/dIa` (ny × na) is
moderate and precomputable.

### 10.2 Dense vs. Sparse Matrices

Current MEX implementation uses dense matrices throughout:
- Green's function matrices `Mxa` (ny × na), `Mxu` (ny × nu): moderate size (typically
  ~3000 × 16 for TCV), dense
- GS operator on the grid: the interior grid is ~50 × 100 = 5000 points; the GS operator is
  sparse (5-point stencil), but the MEX code inverts it as a dense system via LU factorization

For Python:
- `Mxa`, `Mxu`: numpy arrays, multiply with `@` operator
- GS operator: use `scipy.sparse.linalg.spsolve` with the 5-point stencil as a sparse matrix
- For large grids: consider iterative solvers (GMRES, BiCGSTAB)

### 10.3 Profile Parametrization in Python

The `bffbt` profiles translate directly to NumPy:

```python
def bffbt_profile(P, FN):
    pf, pl, pe, tf, tl, te = P
    f = 1 - FN
    # zero outside boundary
    f = np.where(f > 0, f, 0.0)
    x1 = pf * f**pe
    x2 = tf * f**te
    gQ = np.column_stack([
        (pl + (1 - pl) * FN) * x1,
        x2,
        -(tl * f + 1) * FN * x2
    ])
    return gQ

def bffbt_integral(P, FN):
    pf, pl, pe, tf, tl, te = P
    f = 1 - FN
    f = np.where(f > 0, f, 0.0)
    x1 = pf * f**pe * f
    x2 = tf * f**te * f
    IgQ = np.column_stack([
        (-(pl + 1/(pe+1))/(pe+2) + (pl-1)/(pe+2) * FN) * x1,
        (-1/(te+1)) * x2,
        (tl/((te+2)*(te+3)) + 1/((te+1)*(te+2)) + (tl*(te+1)/((te+2)*(te+3)) + 1/(te+2) - (tl/(te+3))*FN)*FN) * x2
    ])
    return IgQ
```

For `bfabmex` polynomials:
```python
def bfab_profile(P, F, FA, FB):
    nP, nT = int(P[0]), int(P[1])
    polys = [F - FB, (F-FB)*(F-FA), (F-FB)*(F-FA)*(2*F-FB-FA)]
    gQ_p = polys[:nP]
    gQ_t = polys[:nT]
    return np.column_stack(gQ_p + gQ_t)
```

### 10.4 Control Point / Cost Function Assembly

The matrices `RaIy`, `RaIe`, `RaFb`, `Ra0` in `fbtcost.m` are assembled per-constraint-category
using Green's functions evaluated at the control points. In Python, replace `greenem` calls with:
- `scipy.special.ellipk/ellipe` for the analytical Green's function of a circular current loop
- Or a precomputed lookup table interpolated via `scipy.interpolate.RegularGridInterpolator`

### 10.5 QP Solver

The linear constrained least-squares problem at each Picard step is:
```
min_u  ||Ra·u - Ra0 - RaIy·Iy||²_w
s.t.   Re·u = Re0 + ReIy·Iy
       Ri·u ≤ Ri0
```

In Python:
- `scipy.optimize.lsq_linear` handles simple bounds
- For general linear constraints, `scipy.optimize.minimize` with `method='SLSQP'`
- Or `qpsolvers` library (`osqp`, `ecos`, `quadprog` backends)
- The FBT `ipm-optimized` QP solver corresponds to an Interior Point Method

### 10.6 Key Data Types

| MATLAB | Python/NumPy |
|--------|-------------|
| `int8` mask `Opy` | `numpy.bool_` or `numpy.uint8` |
| Complex function handles `bfct` | Python callable/closure with same signature |
| `varargin{:}` | `*args` |
| `struct` | `dataclasses.dataclass` or `SimpleNamespace` |
| Cell arrays of cell arrays | `list of lists` |

### 10.7 Critical Dependencies to Implement

1. `greenem` — Green's function for mutual inductance between points and filaments
2. `meqpdom` — plasma domain identification (finding magnetic axis, LCFS, X-points)
3. `meqFx` — GS equation solver (inversion of the Grad-Shafranov operator)
4. `meqinitIy` — initial plasma current distribution
5. `meqagconc` — assembly of profile constraint equations
6. `meqopt` — optimization loop (at minimum: Picard iteration + linear QP)

The Green's function `greenem` is the most critical external dependency — it provides the mutual
inductance between any point in space and a set of filamentary conductors using the complete
elliptic integrals formula.

---

## Appendix: Key Symbol Reference

| Symbol | Meaning | Units |
|--------|---------|-------|
| ψ (psi, Fx) | Poloidal flux function | Wb/rad |
| FA, FB | Flux at magnetic axis, LCFS | Wb/rad |
| FN | Normalized flux (0 at axis, 1 at LCFS) | — |
| Iy | Plasma toroidal current density on grid | A |
| Ia | Active coil currents | A |
| Iu | Passive conductor currents | A |
| Va | Active coil voltages | V |
| Fb | Boundary flux fitting parameter | Wb/rad |
| ag | Profile coefficient vector | varies |
| ng | Number of profile basis functions | — |
| na | Number of active coils | — |
| nu | Number of passive conductors | — |
| nD | Number of plasma domains | — |
| ny | Number of interior grid points | — |
| gp* | Control point constraint parameters | varies |
| bfct | Basis function handle | — |
| bfp | Basis function parameters | — |
| TYG | Plasma current matrix (ny × ng) | A/Wb |
| Mxa | Coil-to-grid mutual inductance (ny × na) | Wb/(A·rad) |
| Ip | Total plasma current | A |
| rBt | R·Bφ at reference radius | T·m |
| qA | On-axis safety factor | — |
| βp | Poloidal beta | — |
| li | Normalized internal inductance | — |
