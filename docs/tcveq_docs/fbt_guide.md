# FBT Practical Usage Guide

This guide explains how to use FBT (Free Boundary Tokamak) for TCV equilibrium design and
analysis. It assumes you have the MEQ toolbox in your MATLAB path and access to the TCV MDS+
database (for actual shot data) or the ANA test tokamak (for standalone testing).

---

## 1. Minimal Example: Compute an Equilibrium for the ANA Test Tokamak

The ANA ("analytical") tokamak is a simple test machine built into MEQ. It requires no database
connection and is ideal for testing, learning, and unit testing.

```matlab
% Compute a single-null diverted equilibrium on the ANA test tokamak
% Shot 2 = single null diverted configuration (see fbtxana.m for all options)
[L, LX, LY] = fbt('ana', 2, []);

% L  = pre-computed data (geometry, matrices)
% LX = input data (constraints, target coil currents, profile parameters)
% LY = output equilibrium data

% Check convergence
fprintf('Converged: %d\n', LY.isconverged)
fprintf('Ip = %.1f kA\n', LY.Ip/1e3)
fprintf('qA = %.2f\n', LY.qA)
fprintf('chi (residual) = %.2e\n', LY.chi)

% Plot the flux surfaces (if you have a plotting function)
% meqplot(L, LY)
```

### 1.1 Pre-defined ANA Shots

From `fbtxana.m`:

| Shot | Description |
|------|-------------|
| 0 | No plasma |
| 1 | Limited, circular |
| 2 | Single null diverted |
| 3 | Single null + secondary X-point |
| 4 | Double null |
| 5 | Squashed (vertically stable) |
| 6 | Negative triangularity, limited |
| 11 | High elongation |
| 21 | Circular with X-point outside limiter |
| 81 | Two droplets (no X-point) |
| 82 | Doublet (two lobes + internal X-point) |
| 83 | Diverted droplet + second in private flux |
| 84 | Diverted droplet + second in public flux |
| 85 | Two diverted droplets |
| 88 | Doublet with mantle current |
| 91 | Vacuum, quadrupole field |

Time-dependent shots:

| Shot | Description |
|------|-------------|
| 101 | Current ramp-up to limiter and back |
| 201 | Same with time-derivative constraints |
| 301 | Transition from circular to elongated |

```matlab
% Time-dependent example
[L, LX, LY] = fbt('ana', 101, []);
% LY now contains equilibria at all time slices
fprintf('Found %d converged slices out of %d\n', sum(LY.isconverged), numel(LY.t))
```

---

## 2. How to Specify the Current Profile

The current profile is specified via:
1. The basis function type (`P.bfct`)
2. The basis function parameters (`P.bfp` or `LX.bfp`)
3. The profile constraints in `LX` (`Ip`, `qA`, `bp`, etc.)

### 2.1 The Default TCV Profile: `bffbt`

For TCV FBT (the standard shot design tool), the default is `bffbt` with 6 parameters:
```
P = [ppfac, ppal, ell1, ttfac1*rh², omega, emm1]
  = [pf,    pl,   pe,   tf,        tl,    te  ]
```

Physical interpretation:
- `ppfac (pf)`: amplitude of pressure gradient p'
- `ppal (pl)`: pressure peaking factor (1 = peaked at axis, 0 = flat)
- `ell1 (pe)`: pressure edge exponent (higher = steeper near boundary)
- `ttfac1*rh² (tf)`: amplitude of diamagnetic current TT' (∝ R₀² for dimensional consistency)
- `omega (tl)`: shape of second TT' basis function
- `emm1 (te)`: TT' edge exponent

Typical TCV defaults (from `bfpfbt`):
```matlab
D.ppfac  = 1;     % Normalized pressure amplitude
D.ppal   = 1;     % Peaked profile
D.ell1   = 1;     % Linear p' edge
D.ttfac1 = 1;     % TT' amplitude (scaled by rh²)
D.omega  = 0;     % Second TT' function weight
D.emm1   = 1.8;   % TT' edge exponent
```

These parameters, combined with the constraints `Ip`, `qA`, `bp`, fully determine the profile.

### 2.2 Specifying Profile Constraints

The profile constraints are set in `LX`:

```matlab
LX.Ip  = 200e3;  % Total plasma current [A]
LX.qA  = 1.0;    % On-axis safety factor
LX.bp  = 0.5;    % Poloidal beta (optional — used if 'Wk' is in fbtagcon)
```

The three constraints (`Ip`, `qA`, `bp`) are used by default to determine the 3 basis function
coefficients `ag`. This is exact: with 3 constraints and 3 unknowns, the solution is unique.

### 2.3 Choosing a Different Basis Function

For more flexibility, use `bfabmex` (polynomials):

```matlab
% 3 polynomial p' functions + 2 TT' functions (5 basis functions total)
% Needs 5 constraints: Ip, qA, bp, li, [one more]
[L, LX, LY] = fbt('ana', 1, [], 'bfct', @bfabmex, 'bfp', [3, 2]);
```

For tabular (arbitrary shape) profiles:
```matlab
% Create tabular profile
FN = linspace(0, 1, 50)';
gNg = [(1-FN).^2, (1-FN).^3];   % Two p' functions (peaked)
IgNg = bfprmex(gNg);              % Compute integrals (requires MEX)
bfp = struct('gNg', gNg, 'IgNg', IgNg, 'fPg', [1;1], 'fTg', [0;0]);

[L, LX, LY] = fbt('ana', 1, [], 'bfct', @bfli, 'bfp', bfp);
```

For spline profiles:
```matlab
% Requires bfpsp to set up knot distribution
P = bfpsp('knots', [0, 0.25, 0.5, 0.75, 1.0], 'order', 5);
[L, LX, LY] = fbt('ana', 1, [], 'bfct', @bfsp, 'bfp', P);
```

### 2.4 Time-Varying Profile Parameters

For time-varying `bffbt` profiles, set `LX.bfp` as a (6 × nt) matrix:

```matlab
LX.bfp = [ppfac_vec;   % 1 × nt row
           ppal_vec;    % 1 × nt row
           ell1_vec;    % 1 × nt row
           ttfac_vec;   % 1 × nt row
           omega_vec;   % 1 × nt row
           emm1_vec];   % 1 × nt row
```

This is enabled by `P.LXbfp = true` in `fbtptcv` (default for TCV).

---

## 3. How to Set Coil Currents

### 3.1 Soft Coil Current Targets (Default)

By default, FBT minimizes the deviation of coil currents from targets `gpia` (usually zero)
with a specified uncertainty `gpid`. This allows the optimizer to freely choose coil currents
to satisfy the shape constraints.

```matlab
% Set coil current targets and their uncertainties
LX.gpia = zeros(na, 1);       % Target = 0 A for all coils
LX.gpid = 1e4 * ones(1, nt);  % Allow ±10 kA deviation
LX.gpie = ones(na, nt);       % Equal weight for all coils
```

A small `gpid` tightly constrains the coil currents to their targets. A large `gpid` gives the
optimizer freedom.

### 3.2 Hard Coil Current Constraints

To fix a specific coil current to a prescribed value:
```matlab
LX.gpia(k, :) = Ia_target;    % Target value
LX.gpie(k, :) = 0;            % Zero weight = equality constraint
```

### 3.3 Coil Current Limits (Inequality Constraints)

Coil limits are set in `P` (configuration parameters), not `LX`:
```matlab
P.limc = eye(na);              % Individual coil limits
P.limu = 20e3 * ones(na, 1);  % Upper limit: +20 kA
P.liml = -20e3 * ones(na, 1); % Lower limit: -20 kA
P.limm = 0.95;                 % 5% safety margin
```

For TCV, these are loaded automatically from the MDS+ tree in `fbtptcv`.

### 3.4 Dipole Regularization

To regularize coil current differences (prevent large adjacent coil current swings):
```matlab
LX.gpdd = 1e3;               % Global dipole weight [A]
LX.gpdw = ones(na, 1);       % Equal weights
LX.gpda = zeros(na-1, 1);    % Target current differences = 0
```

Large `gpdd` strongly penalizes current differences. This is important for elongated plasmas
where adjacent coils often have opposite-sign currents; the dipole term prevents solutions where
non-neighboring coils dominate.

---

## 4. Specifying Shape Constraints (Control Points)

Shape constraints are the core of FBT — they define what equilibrium shape you want to achieve.
Use `fbtgp` to add control points to an LX (or P) structure.

### 4.1 Syntax

```matlab
LX = fbtgp(LX, r, z, b, fa, fb, fe, br, bz, ba, be, cr, cz, ca, ce, vrr, vrz, vzz, ve)
```

Arguments:
- `r, z`: R, Z coordinates of the control point(s) [m]
- `b`: 1 if point is on the plasma boundary (used for initial guess bounding box), 0 otherwise
- `fa`: target flux value
- `fb`: 0 = absolute flux value; 1 = relative to boundary flux (common offset `Fb`)
- `fe`: weight of flux constraint (Inf = inactive, 0 = equality)
- `br, bz`: target Br, Bz values (NaN = don't constrain)
- `ba`: target field angle (used when br/bz are NaN)
- `be`: weight of field constraint (Inf = inactive, 0 = equality)
- `cr, cz`: target ∂²ψ/∂R², ∂²ψ/∂R∂Z (NaN = don't constrain)
- `vrr, vrz, vzz`: vacuum field gradient targets (NaN = don't constrain)

### 4.2 Common Patterns

**Circular limiter plasma:**
```matlab
thp = linspace(0, 2*pi, 9)'; thp = thp(1:end-1);
r0 = 1.05; z0 = 0; a0 = 0.4;
ra = r0 + a0*cos(thp);
za = z0 + a0*sin(thp);
%                r, z, b, fa, fb, fe, br, bz, ba, be, ...
LX = fbtgp(LX, ra, za, 1,  0,  1,  1, [],[], [],[], [], [],[], [], [],[], [],  []);
% Make the first point an exact boundary constraint:
LX.gpfe(1) = 0;
```

**Single-null diverted plasma with X-point:**
```matlab
% Boundary contour points (with Fb=1 for relative flux)
LX = fbtgp(LX, ra, za, 1, 0, 1, 1, [], [], [], [], [], [], [], [], [], [], [], []);

% X-point: Br = Bz = 0, on boundary (gpb=1), exact constraints
rX = 0.85; zX = -0.4;
LX = fbtgp(LX, rX, zX, 1, 0, 1, 0, 0, 0, [], 0, [], [], [], [], [], [], [], []);
```

The key for X-points:
- `gpbr = 0, gpbe = 0`: Br = 0 exactly
- `gpbz = 0, gpbe = 0`: Bz = 0 exactly
- `gpfe = 0`: flux = Fb exactly (X-point is on the LCFS)

**Flux loop at a specific point:**
```matlab
% Constrain flux at (r=0.9, z=0) to be within 5 mWb of Fb
LX = fbtgp(LX, 0.9, 0, 0, 0, 1, 1, [], [], [], [], [], [], [], [], [], [], [], []);
% fe=1 means this enters cost function with weight gpfd
```

**B-probe constraint:**
```matlab
% Constrain Br = 0 at a specific point (saddle point below divertor)
LX = fbtgp(LX, r_probe, z_probe, 0, [], [], [], 0, [], [], 0, [], [], [], [], [], [], [], []);
% gpbr=0, gpbe=0 means exact constraint on Br
```

**Snowflake (second-order null):**
```matlab
% Point where Br = Bz = 0 AND psirr = psirz = 0 (gradient of B also zero)
LX = fbtgp(LX, r_sf, z_sf, 1, 0, 1, 0, 0, 0, [], 0, 0, 0, [], 0, [], [], [], []);
```

### 4.3 The `gpfb` Field: Absolute vs. Relative Flux

`gpfb = 0`: the constraint is `ψ(r,z) = gpfa` (absolute flux)
`gpfb = 1`: the constraint is `ψ(r,z) - Fb = gpfa` where `Fb` is the free boundary flux

In practice, the LCFS is defined by `ψ = Fb` (with `gpfa = 0` and `gpfb = 1`). The fitting
parameter `Fb` is adjusted by the optimizer to satisfy the shape constraints. Points near the
LCFS should have `gpfa = 0`, `gpfb = 1`.

Strike points on the divertor target typically have `gpfe = 0` (exact) and `gpfb = 1` (on LCFS),
with `gpbr = 0` or `gpbz = 0` to specify the field angle.

---

## 5. Understanding the Output

### 5.1 Convergence Diagnostics

```matlab
% After running FBT:
LY.isconverged   % 1 if converged, 0 if not
LY.chi           % RMS of cost function residual (should be ~tol after convergence)
LY.chie          % RMS of constraint residual (should be near machine epsilon)
LY.werr          % Weighted cost function error
LY.uerr          % Unweighted (absolute) error
```

If `LY.isconverged = 0`, common causes:
- Conflicting constraints (X-point too close to boundary points)
- Coil current limits too tight for the requested shape
- Profile constraints inconsistent (e.g., `qA` too high for given `Ip`, `bp`)
- Grid too coarse to resolve the plasma domain

### 5.2 Equilibrium Properties

```matlab
LY.Ia           % Active coil currents found [A] (na × nt)
LY.Ip           % Total plasma current [A]
LY.qA           % On-axis safety factor
LY.bp           % Poloidal beta
LY.li           % Normalized internal inductance
LY.rA, LY.zA   % Magnetic axis location [m]
LY.rX, LY.zX   % X-point location [m] (NaN if limited)
LY.lX           % 1 if X-point (diverted), 0 if limited
LY.lB           % 1 if limited, 0 if diverted
LY.FA           % Flux at magnetic axis [Wb/rad]
LY.FB           % Flux at LCFS [Wb/rad]
LY.Fb           % Boundary flux fitting parameter [Wb/rad]
LY.Sp           % Plasma cross-sectional area [m²]
LY.Ep           % Plasma elongation (height/width of bounding box)
```

### 5.3 Interpreting the Cost Function Residual

The cost function IDs (`LY.ID`) label each term:
- ID=1 `Coils`: active coil current term
- ID=2 `Dipoles`: adjacent coil current difference term
- ID=3 `Passives`: passive conductor current term
- ID=4 `Flux`: flux constraint at control points
- ID=5 `Fields`: Br/Bz field constraint
- ID=6 `Hessian`: flux curvature constraint
- ID=7 `Vacuum`: vacuum field gradient constraint
- ID=8 `Voltage`: active coil voltage (circuit case)

If the total `chi` is large but `chie` is small, the soft constraints (shape) cannot be fully
satisfied while respecting the hard constraints (coil limits, circuit equations). This indicates
a physically incompatible design.

### 5.4 Reconstructing the Flux Map

```matlab
% Get the full flux map (requires requesting it explicitly)
% [GUESS] By default LY.Fy is stored; check LY fields
% Alternatively reconstruct from Ia:
Fx_vac = L.G.Mxa * LY.Ia + LY.Fb;  % Vacuum contribution at interior grid
% Plasma contribution requires Iy (may not be stored by default)
```

---

## 6. TCV-Specific Usage

### 6.1 Running FBT for an Existing TCV Shot

```matlab
% Reconstruct equilibria for TCV shot 70000 at times 0.5, 0.6, 0.7 s
[L, LX, LY] = fbt('TCV', 70000, [0.5, 0.6, 0.7]);

% Check which slices converged
fprintf('%d of %d slices converged\n', sum(LY.isconverged), numel(LY.t))
```

### 6.2 Using `fbte` (FBTE: FBT + Write to Tree)

```matlab
% Run FBT and write results to PCS tree (for shot design)
% Shot -1 is the "MODEL" shot used for preparation
[L, LX, LY] = fbte(-1);
```

`fbte` calls `fbt('TCV', shot, [], 'LYall', true, ...)` and then writes results to the
TCV MDS+ tree using `fbtwtcv`, `meqwtcv`, and `psiwtcv`.

### 6.3 Loading Configuration from MGAMS

The MGAMS (Magnetic parameters for Gams) system stores shot design parameters in the PCS tree.
`fbtptcv` reads these via `mgp(shot)`:

```matlab
% Check what MGAMS loaded for a specific shot
P = fbtptcv(70000);

% Key parameters from MGAMS:
P.gpr    % Control point R coordinates
P.gpz    % Control point Z coordinates
P.gpb    % Boundary indicator
P.gpfa   % Flux targets
P.gpbr   % Br targets
P.gpia   % Coil current targets
P.Ip     % Plasma current target (= P.placu1)
P.rBt    % R*Bt (= P.bzero * 0.88)
P.qA     % On-axis q target (= P.qzero)
P.bfct   % Basis function type (usually @bffbt)
P.ppfac  % bffbt p' amplitude
P.ppal   % bffbt p' peaking
P.ell1   % bffbt p' edge exponent
P.ttfac1 % bffbt TT' amplitude factor
P.omega  % bffbt TT' shape
P.emm1   % bffbt TT' edge exponent
```

### 6.4 Converting FBT Design to MGAMS Parameters

After designing an equilibrium with FBT, convert it to MGAMS format for shot preparation:

```matlab
[L, LX, LY] = fbt('TCV', -1, []);

% Pick a single time slice
LXt = meqxk(LX, 1);

% Convert to MGAMS parameter object
MG = fbt2pcs(L, LXt);

% Inspect the parameters
MG.gpr    % Control point coordinates
MG.bzero  % B0 = rBt / 0.88
MG.placu1 % Ip in kA
% ... write to PCS tree with standard MGAMS tools
```

### 6.5 Classifying Control Points for Shape Control

```matlab
% Get shape control point classification
SC = fbt2SC(L, LX);

SC.rc(:, 1)    % R coordinates for slice 1
SC.zc(:, 1)    % Z coordinates for slice 1
SC.lcX1(:, 1)  % Primary X-point indicator
SC.lcL(:, 1)   % Limiter point indicator
SC.lcS(:, 1)   % Strike point indicator
SC.lcC(:, 1)   % General contour point indicator
SC.dimc(:, 1)  % Labels: 'Surf01', 'Xpri01', 'Limi01', etc.
```

---

## 7. Design Workflow with FBT

### 7.1 The Shot Design Process

The typical TCV shot design workflow using FBT:

1. **Define target shape**: specify boundary points using `fbtgp` (or the MGAMS GUI)
2. **Set plasma parameters**: `Ip`, `rBt`, `qA`, `bp` or `Wk`
3. **Run FBT**: compute the equilibrium and coil currents `Ia`
4. **Check feasibility**: verify coil currents are within limits, shape matches targets
5. **Iterate**: adjust shape points, coil weights, or plasma parameters as needed
6. **Write to PCS tree**: use `fbt2pcs` and MGAMS to prepare the shot

### 7.2 Using `fbtgp` for Shape Design

`fbtgp` is the main interface for building the LX structure:

```matlab
% Start from scratch
LX = struct();
LX.t = 0;

% Define plasma current and B field
LX.Ip = 200e3;    % 200 kA
LX.rBt = 0.92;    % R*Bt [T*m], for TCV R0=0.88m, Bt=1.04T
LX.qA = 1.0;
LX.bp = 0.5;

% Add boundary points using parametric shape generator bouD:
%   bouD(n, r0, z0, a0, kappa_top, kappa_bot, delta_top, delta_bot, ...)
[ra, za] = bouD(9, 0.88, 0, 0.25, 1.5, 1.5, 0.3, 0.3, 0, 0, 0, 0);
LX = fbtgp(LX, ra, za, 1, 0, 1, 1, [], [], [], [], [], [], [], [], [], [], [], []);

% Make one point an exact boundary constraint
LX.gpfe(1) = 0;

% Add X-point
rX = 0.65; zX = -0.55;
LX = fbtgp(LX, rX, zX, 1, 0, 1, 0, 0, 0, [], 0, [], [], [], [], [], [], [], []);
```

### 7.3 Time-Dependent Sequences

For shot design, FBT is typically run at multiple time points:

```matlab
% Build an LX with multiple time slices manually
LX1 = struct(); LX1.t = 0.2;
LX1.Ip = 100e3; LX1.rBt = 0.88; LX1.qA = 1.2;
% ... add control points for time 0.2 s

LX2 = struct(); LX2.t = 0.5;
LX2.Ip = 200e3; LX2.rBt = 0.88; LX2.qA = 1.0;
% ... add control points for time 0.5 s

% Concatenate
LX = fbtxcat(LX1, LX2);

% Interpolate to finer time grid if needed
t_dense = linspace(0.2, 0.5, 10);
LX_dense = fbtxinterp(LX, t_dense, 'linear');

% Run FBT on full sequence
LY = fbtt(L, LX_dense);
```

### 7.4 Typical Parameter Values for TCV

From practice and code inspection:

| Parameter | Typical TCV value | Meaning |
|-----------|-------------------|---------|
| `Ip` | 100–500 kA | Plasma current |
| `rBt` | 0.88 × (0.8–1.4) T·m | For BT = 0.9–1.6 T at R0=0.88m |
| `qA` | 0.9–2.0 | On-axis safety factor |
| `bp` | 0.1–3.0 | Poloidal beta |
| `gpfd` | 5e-3 Wb/rad | Flux constraint tolerance |
| `gpbd` | 5e-3 Wb/rad | Field constraint tolerance |
| `gpid` | ~2e5 A | Coil current tolerance |
| `ppfac` | 0.5–2.0 | p' amplitude |
| `ppal` | 0.5–1.0 | p' peaking (1=peaked, 0=flat) |
| `ell1` | 1–3 | p' edge steepness |
| `ttfac1` | 0.5–2.0 | TT' amplitude |
| `emm1` | 1–3 | TT' edge steepness |

### 7.5 Diagnosing Convergence Problems

**Problem: `chi` large, equilibrium not converging**

Checklist:
1. Check control point consistency: do boundary points form a closed, reasonable shape?
2. Check for conflicting constraints: e.g., two points with `gpb=1` far from each other and
   inconsistent flux targets
3. Relax tight constraints: set `gpfe > 0` for some points
4. Increase `niter`: `fbt('TCV', shot, t, 'niter', 200)`
5. Adjust initial guess: set `capaj` (default 0.5, higher = more peaked initial current)

**Problem: Coil current limits violated**

Solutions:
1. Use `limm = 0.9` for a 10% safety margin: `fbt(..., 'limm', 0.9)`
2. Relax shape requirements (move control points, reduce weights)
3. Reduce plasma current
4. Change the OH coil configuration

**Problem: Wrong sign of Ip or Bt**

In `fbtptcv`, signs are handled by `iohfb` (sign of plasma current relative to MGAMS convention)
and `if36fb` (sign of toroidal field). Both affect control point targets.

Internally, FBT works with consistent signs: the assertion
`all(~any(LX.qA .* LX.rBt .* IpD < 0))` checks that the product of signs is positive.

### 7.6 Multi-Domain (Doublet) Design

For TCV doublet plasmas:
```matlab
% Enable doublet mode
LX.idoublet = true;
LX.isaddl = 1;  % Additional domain enabled

% Specify current per lobe (3 domains: upper, lower, mantle)
LX.IpD = [100e3; 100e3; 0];  % kA per domain
LX.qA = [1.0; 1.0];           % q on axis per lobe (no axis in mantle)
LX.bpD = [0.5; 0.5; 0];       % bp per domain

% Upper lobe boundary
[ra, za] = bouD(11, 0.88, +0.4, 0.25, 1.2, 1.2, 0, 0, 0, 0, 0, 0);
LX = fbtgp(LX, ra, za, 1, 0, 1, 1, [], [], [], [], [], [], [], [], [], [], [], []);

% Lower lobe boundary  
[ra, za] = bouD(11, 0.88, -0.4, 0.25, 1.2, 1.2, 0, 0, 0, 0, 0, 0);
LX = fbtgp(LX, ra, za, 1, 0, 1, 1, [], [], [], [], [], [], [], [], [], [], [], []);

% Internal X-point
LX = fbtgp(LX, 0.88, 0, 0, [], [], [], 0, 0, [], 0, [], [], [], [], [], [], [], []);

% For bffbt in doublet mode, use per-domain constraints
P.fbtagcon = {'Ip', 'qA', 'fbtlegacy'};  % or per-domain cell arrays
```

### 7.7 Circuit Equation Mode (Time-Dependent with Vessel)

For time-dependent simulations where vessel eddy currents matter:
```matlab
[L, LX, LY] = fbt('TCV', shot, t, 'circuit', true, 'selu', 'all');
% 'selu', 'all' includes all vessel passive structures
% LX.gpue: vessel current equation weights (Inf for first slice = steady-state)
% LX.Va: active coil voltages (optimization variable instead of Ia)
```

In this mode, `Va` (coil voltages) replaces `Ia` as the optimization variable, and the system
solves:
```
Mee · dIe/dt + Re · Ie + Mey · dIy/dt = [Va; 0]
```

---

## 8. Key Gotchas and Common Mistakes

1. **Sign consistency**: `qA`, `rBt`, and `Ip` must satisfy `sign(qA) = sign(rBt) * sign(Ip)`.
   FBT asserts this and will error if violated.

2. **`gpfb` vs. absolute flux**: boundary points should use `gpfb = 1` (relative to the
   boundary flux `Fb`). Using `gpfb = 0` with `gpfa = 0` constrains the absolute flux to zero,
   which is almost never what you want.

3. **Conflicting X-point and boundary constraints**: if an X-point location `(rX, zX)` is also
   a boundary point, use `mergex = true` (default) to merge the duplicate constraints
   automatically.

4. **`gpfy`, `gpby` plasma contribution flags**: these (default 1) control whether the plasma
   contribution `MSy·Iy` appears in the constraint equation. For vacuum field curvature
   constraints (`gpv*`), `gpvy = 0` by default (vacuum only). Mixing these up leads to
   physically incorrect equilibria.

5. **Doublet with `SQP-broyden`**: the `SQP-broyden` solver is incompatible with doublet
   configurations. Use `useSQP = false` or `solver = 'SQP'`.

6. **`LXbfp` vs. fixed `bfp`**: when `P.LXbfp = true` (default for TCV), the basis function
   parameters are time-dependent and come from `LX.bfp`. If `LX.bfp` is missing, FBT will use
   the static `P.bfp` from `bfpfbt`.

7. **Grid mismatch**: if you pass custom control points at exactly the same (R,Z) as a grid
   point and `MHyinterp = false`, FBT throws an error. Use `MHyinterp = true` (default) to
   avoid this.
