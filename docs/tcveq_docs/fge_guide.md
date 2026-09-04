# FGE User Guide — TCV Practical Reference

**Codebase**: MEQ toolbox, Swiss Plasma Center EPFL  
**Target audience**: Users of FGE at TCV; Python surrogate developers

---

## 1. Minimal Working Example

### 1.1 Run FGE from a TCV Shot Database

```matlab
%% Basic FGE run: 100ms simulation starting at t=0.2s
tok   = 'tcv';
shot  = 65668;
t0    = 0.2;       % start time [s]
dt    = 0.5e-3;    % time step [s]
tend  = 0.3;       % end time [s]
t     = t0:dt:tend;

%% Get L structure (geometry + parameters) and LX (time-varying inputs)
[L, LX, LY] = fge(tok, shot, t);
```

The three output arguments trigger:
- `nargout >= 1`: compute `L` (parameters + geometry; calls `fgeptcv`, `fgegtcv`, `fgep`, `fgeg`, `fgec`)
- `nargout >= 2`: compute `LX` (input data; calls `fgex` → `fgsxtcv`)
- `nargout >= 3`: run the simulation (calls `fget`) and return `LY`

To decouple the steps (allows modifying LX before running):

```matlab
[L, LX] = fge(tok, shot, t);    % parameters + data, no simulation
LY = fget(L, LX);                % run simulation
```

### 1.2 Run from an Analytical Test Case (ANAMAK, no database)

```matlab
%% ANAMAK circular equilibrium test
tok  = 'ana';
shot = 1;
t    = (0:1e-4:5e-3);  % 5ms, 100µs steps

[L, LX, LY] = fge(tok, shot, t);
```

### 1.3 Run from an FBT/LIUQE Initial Condition (TCV, detailed)

```matlab
%% Step 1: Get L structure only (no time data yet)
L = fge(tok, shot, []);  % empty t → only L returned

%% Step 2: Get initial condition from LIUQE reconstruction
LX0 = fgsxtcv(t0, L);  % FGS initial condition from LIUQE at t0

%% Step 3: Build LX for the full time window
LX = fgex(tok, t, L, LX0);  % pass LX0 explicitly

%% Step 4: (Optional) modify inputs
LX.Va(2, :) = LX.Va(2, :) + 20;  % add 20V step on coil #2

%% Step 5: Run
LY = fget(L, LX);
```

---

## 2. How to Specify Coil Voltage Waveforms

### 2.1 Structure of LX.Va

`LX.Va` is a `(na × nt)` matrix where:
- `na` = number of active coils (e.g., 19 for TCV: OH + E + F + G coils)
- `nt` = number of time steps

Time is `LX.t` (1 × nt vector), uniformly spaced.

The coil names are in `L.G.dima` (cell array of strings, `na` elements).

```matlab
% Find index of a specific coil
icoil = find(strcmp(L.G.dima, 'E1'));  % coil named 'E1'
fprintf('E1 is coil index %d\n', icoil);
```

### 2.2 Open-Loop Voltage Command

By default (no controller), `LX.Va` is the complete voltage waveform applied to the coils. The initial LX.Va is set to the feedforward trajectory that maintains the initial equilibrium:

```matlab
LX.Va = LX.Va(:, 1) + LX.Vadot0(:, 1) * (LX.t - LX.t(1));
```

To override with a custom waveform:

```matlab
%% Constant voltage step on all coils
LX.Va = repmat(LX.Va(:, 1), 1, nt);  % hold initial voltage

%% Step on coil 'OH'
iOH = find(startsWith(L.G.dima, 'OH'));
t_step = 0.21;  % step start time
LX.Va(iOH, LX.t >= t_step) = LX.Va(iOH, LX.t >= t_step) + 50;  % +50V
```

### 2.3 Power Supply Delays

TCV power supply delays are stored in `L.G.Vadelay` (na × 1 vector, in seconds):
- E, F, OH coils: 0.5 ms
- G coils: 0.05 ms

The power supply model (`L.P.psfct = @meqps`) applies these delays automatically. The commanded voltage `Vacmd` is filtered to produce the actual voltage `Vact` applied to the circuit equation.

To disable delays (e.g., for linearization comparison):
```matlab
L.G.Vadelay(:) = 0;
```

### 2.4 Feedback Controller

To use a feedback controller:

```matlab
L.P.ctrlfct = @my_controller;   % function handle
L.P.ctrlpar.dt = 1e-3;          % controller time step [s]
% mycontroller must implement:
%   [cstate, Vactr] = my_controller('init', [], ctrlpar, L, LY0)
%   [cstate, Vactr, ~, debug] = my_controller('step', cstate, ctrlpar, L, LYt)
```

The controller runs at `ctrlpar.dt` intervals (must be >= simulation `dt`). It receives the previous LY and outputs a voltage command vector `Vactr (na × 1)`.

---

## 3. How to Set Up the CDE Model

### 3.1 Available CDE Types

```matlab
% No CDE (default): profile constraints fully determine ag
L.P.cde = '';

% 0D rigid plasma model (most common for TCV control)
L.P.cde = 'OhmTor_rigid_0D';

% 0D toroidal Ohm's law
L.P.cde = 'OhmTor_0D';

% 0D steady-state RAPTOR-derived
L.P.cde = 'cde_ss_0D';

% 1D full current diffusion in flux coordinates
L.P.cde = 'cde_1D';

% 1D toroidal Ohm's law under test functions
L.P.cde = 'cde_OhmTor_1D';
```

When a CDE is active, one of the `agcon` constraints must be removed (since the CDE replaces it). Typically remove `'Ip'` and use `cde='OhmTor_rigid_0D'`:

```matlab
[L, LX] = fge(tok, shot, t, ...
    'cde',   'OhmTor_rigid_0D', ...
    'agcon', {'bp', 'qA'});     % no Ip constraint — CDE provides it
```

### 3.2 CDE Input Fields in LX

These fields must be provided in `LX` (time-varying):

| Field | Dimension | Description |
|---|---|---|
| `LX.Rp` | `nD × nt` | Plasma resistance [Ohm] |
| `LX.Ini` | `nD × nt` | Non-inductive current drive [A] |
| `LX.Lp` | `nD × nt` | Plasma inductance [H] (if `iLpext=true`) |
| `LX.signeo` | `nQ × nD × nt` | Neoclassical conductivity profile [S/m] (must be > 0 everywhere) |

For TCV, these come from `fgsxtcv`. For synthetic tests:

```matlab
%% Synthetic 0D CDE inputs
Rp_val  = 1e-6;    % [Ohm] - typical TCV plasma resistance
Ini_val = 0;       % [A]   - no non-inductive drive
Lp_val  = 5e-7;    % [H]   - rough TCV plasma inductance

nt = numel(LX.t);
LX.Rp  = repmat(Rp_val,  1, nt);
LX.Ini = repmat(Ini_val, 1, nt);
LX.Lp  = repmat(Lp_val,  1, nt);
```

### 3.3 Computed vs Provided Lp

The `iLpext` parameter (default `true`) controls how `Lp` is used:

```matlab
% Option 1: Provide Lp externally (from database or model)
L.P.iLpext = true;
LX.Lp = ...;  % must fill this

% Option 2: Compute Lp from plasma current distribution
L.P.iLpext = false;
% Lp = (Iy^T * Myy * Iy) / Ip^2  — computed inside fgeF automatically
% LX.Lp is ignored
```

### 3.4 Setting Up for 1D CDE

For 1D CDE, the `agcon` list must be empty (the CDE replaces all profile constraints):

```matlab
[L, LX] = fge(tok, shot, t, ...
    'cde',   'cde_OhmTor_1D', ...
    'agcon', {});      % CDE provides all ng degrees of freedom
```

The `signeo` profile must cover all normalized flux coordinates `pQ^2 ∈ [0, 1]`.

### 3.5 Initialization with CDE

The initialization (`fgeinit`) automatically handles the CDE:
- For 0D CDEs: replaces CDE constraint with `Ip` for the static initial solve
- For 1D CDEs: uses `static_Iy_1D` to constrain the profile to match the input Iy

Set `ssinit = false` to skip the steady-state initialization step:
```matlab
[L, LX] = fge(tok, shot, t, 'cde', 'OhmTor_rigid_0D', 'agcon', {'bp','qA'}, ...
    'ssinit', false);
```

---

## 4. Extracting the Linearized Model for Control

### 4.1 Compute Linearization at a Single Point

```matlab
%% Get L and LX
[L, LX] = fge(tok, shot, t);
LX0 = meqxk(LX, 1);  % single time slice

%% Compute linearization
L = fgel(L, LX0);

%% Extract continuous-time state-space model
sys = fgess(L);         % continuous time (Ts=0)
sysd = fgess(L, dt);    % discrete time, backward Euler
```

### 4.2 The L.lin Structure After fgel

After `fgel`, `L.lin` contains:

```matlab
L.lin.S          % Mass/inductance matrix (S * dxdt = K*x + U*u)
L.lin.K          % "Stiffness" matrix
L.lin.U          % Input matrix
L.lin.Udot       % dCodt input matrix
L.lin.x0L        % Linearization state (reduced xL)
L.lin.u0L        % Linearization inputs [Va0; Co0; Ini0; Rp0; 0]
L.lin.res0L      % Residual at linearization point (non-zero if transient init)
L.lin.maskL      % Boolean mask: which x components are dynamic states
L.lin.dxdxL      % Full-state-to-reduced-state mapping dx = dxdxL * dxL
L.lin.Ip         % Plasma current at linearization point
L.lin.ag         % Basis function coefficients at linearization point
% ... many d(output)/dxL and d(output)/du fields
```

### 4.3 The State-Space sys Object

```matlab
sys = fgess(L);

% Access matrices
A = sys.A;  % [nS × nS] — typically ns = na + nu + np
B = sys.B;  % [nS × nU] — nU = na + nC + nD + nD + nC
C = sys.C;  % [nM × nS] — nM = total measurement size
D = sys.D;  % [nM × nU]

% Named groups
sys.InputGroup.Va     % indices of Va inputs in u
sys.InputGroup.Co     % indices of constraint inputs
sys.InputGroup.Ini    % indices of Ini inputs (CDE only)
sys.InputGroup.Rp     % indices of Rp inputs (CDE only)
sys.InputGroup.dCodt  % indices of dCo/dt inputs

sys.OutputGroup.Ia    % indices of Ia in y
sys.OutputGroup.Iu    % indices of Iu in y
sys.OutputGroup.rA    % index of rA (magnetic axis R)
sys.OutputGroup.zA    % index of zA (magnetic axis Z)
sys.OutputGroup.Bm    % indices of B-probe measurements
sys.OutputGroup.Ff    % indices of flux loop measurements
sys.OutputGroup.Ip    % index of plasma current
% ... etc.
```

### 4.4 Extract RZIP-Style Plasma-Position Model

The linearized response of axis position to coil voltages:

```matlab
% Step response of zA to Va perturbation
sys_zA_Va = sys('zA', 'Va');  % MATLAB ss subsystem selection

% Or manual extraction:
i_zA = sys.OutputGroup.zA;
i_Va = sys.InputGroup.Va;
A_zA_Va = sys.A;
B_zA_Va = sys.B(:, i_Va);
C_zA_Va = sys.C(i_zA, :);
D_zA_Va = sys.D(i_zA, i_Va);
```

### 4.5 Compute Eigenvalues (Growth Rates)

```matlab
%% Find unstable modes (vertical instability etc.)
[Ei, Vi, LYi] = fgeeig(L);

fprintf('Most unstable growth rate: %.1f /s (%.1f rad/s)\n', ...
    real(Ei(1)), imag(Ei(1)));

%% Plot the eigenmode current pattern
meqplott(L, LYi(1));
```

### 4.6 Batch Linearization Over Time

For a gain-scheduled controller, compute linearizations at multiple equilibria:

```matlab
%% First run a NL simulation
[L, LX, LY] = fge(tok, shot, t);

%% Compute linearizations at every 10th time step
Llin = L;  % use same geometry/parameters
t_lin_idx = 1:10:numel(LY.t);
LY_subset = meqxk(LY, t_lin_idx);

[Llin, LYlin, success_mask] = fgelinall(Llin, L, LY_subset, 'dt', dt, 'robust', true);

%% Access linearization at time index k
fprintf('%d linearizations computed\n', numel(Llin.lin));
sys_k = fgess(Llin);  % uses Llin.lin (last one)
% or access specific one:
% Llin.lin(k).sysA, Llin.lin(k).sysB, ...
```

### 4.7 Run Linearized FGE Simulation

To use the linearized model instead of the full nonlinear solver:

```matlab
%% Compute linearization
L = fgel(L, meqxk(LX, 1));

%% Run linear simulation
L.P.lin = true;
LY_lin = fget(L, LX);

%% Compare
figure;
subplot(2,1,1); plot(LY_lin.t, LY_lin.zA, LY.t, LY.zA); legend('lin','NL');
subplot(2,1,2); plot(LY_lin.t, LY_lin.Ip, LY.t, LY.Ip); legend('lin','NL');
```

---

## 5. Understanding the State and Output

### 5.1 LY — Output Structure

`LY = fget(L, LX)` returns a structure-of-arrays where each field has size `(dim × nt)`:

| Field | Size | Description |
|---|---|---|
| `LY.t` | `1 × nt` | Time vector [s] |
| `LY.Fx` | `nzx × nrx × nt` | Poloidal flux [Wb] |
| `LY.Iy` | `nzy × nry × nt` | Plasma current density [A/cell] |
| `LY.Ia` | `na × nt` | Active coil currents [A] |
| `LY.Iu` | `nu × nt` | Vessel eddy currents [A] |
| `LY.Ip` | `1 × nt` | Total plasma current [A] |
| `LY.rA` | `nD × nt` | Magnetic axis R position [m] |
| `LY.zA` | `nD × nt` | Magnetic axis Z position [m] |
| `LY.FA` | `nD × nt` | Axis flux value [Wb] |
| `LY.FB` | `nD × nt` | Boundary (LCFS) flux value [Wb] |
| `LY.rB` | `nB × nt` | Boundary R positions [m] |
| `LY.zB` | `nB × nt` | Boundary Z positions [m] |
| `LY.rX` | `nX × nt` | X-point R positions [m] |
| `LY.zX` | `nX × nt` | X-point Z positions [m] |
| `LY.bp` | `1 × nt` | Poloidal beta |
| `LY.qA` | `1 × nt` | Safety factor at axis |
| `LY.li` | `1 × nt` | Internal inductance |
| `LY.Wk` | `1 × nt` | Plasma kinetic energy [J] |
| `LY.ag` | `ng × nt` | Profile basis coefficients |
| `LY.Bm` | `nm × nt` | Magnetic probe measurements [T] |
| `LY.Ff` | `nf × nt` | Flux loop measurements [Wb] |
| `LY.Va` | `na × nt` | Applied coil voltages [V] |
| `LY.isconverged` | `1 × nt` | Convergence flag per time step |
| `LY.res` | `1 × nt` | Final residual per time step |
| `LY.niter` | `1 × nt` | Newton iterations per time step |

### 5.2 Accessing a Single Time Slice

```matlab
LYt = meqxk(LY, it);  % it = time index
% LYt has same fields but scalar in time
```

### 5.3 Plotting

```matlab
%% Time traces
figure;
subplot(3,1,1); plot(LY.t, LY.Ip/1e3); ylabel('Ip [kA]');
subplot(3,1,2); plot(LY.t, LY.zA*100); ylabel('zA [cm]');
subplot(3,1,3); plot(LY.t, LY.Ia); ylabel('Ia [A]'); xlabel('t [s]');

%% Equilibrium at one time slice
meqplott(L, meqxk(LY, end));  % plot last time step

%% Evolution animation (or multi-panel)
meqplotevo(L, LY);
```

### 5.4 L Structure Key Fields

```matlab
L.G.na    % number of active coils
L.G.nu    % number of vessel current components
L.G.dima  % coil names {na × 1 cell}
L.G.dimu  % vessel segment names {nu × 1 cell}
L.G.rx    % R coordinates of full grid (nzx × nrx)
L.G.zx    % Z coordinates of full grid (nzx × nrx)
L.G.rl    % R coordinates of limiter wall
L.G.zl    % Z coordinates of limiter wall
L.P.r0    % major radius [m]
L.ng      % number of profile basis functions
L.nC      % number of profile constraints
L.np      % number of CDE equations
L.nD      % number of plasma domains (1=normal, 2=doublet)
L.Ip0     % reference plasma current [A]
L.Ia0     % reference coil currents (na × 1) [A]
L.Iu0     % reference vessel currents (nu × 1) [A]
L.Fx0     % reference flux scale [Wb]
```

---

## 6. Common Failure Modes

### 6.1 "FGS did not converge for FGE initialization"

Cause: The initial equilibrium from FBT/LIUQE cannot be reproduced by FGS.

**Fixes**:
- Set `ssinit = false` to skip steady-state initialization
- Check that `signeo` is strictly positive
- Reduce the number of profile constraints or change `agcon`
- Try `initguess = {'LX', 'vacuum'}` for more robust initialization
- Check that the chosen `cde` type is compatible with the profile constraints

```matlab
[L, LX] = fge(tok, shot, t, 'ssinit', false);
```

### 6.2 Non-convergence During Simulation

Symptom: `LY.isconverged(it) = false` for some steps.

**Diagnostics**:
```matlab
L.P.debug = 2;  % Verbose output per time step
LY = fget(L, LX);
% Look for 'Newton solver did not converge' messages
```

**Fixes**:
- Reduce `dt` (improve time resolution)
- Increase `kmax` (more Newton iterations): `L.P.kmax = 100`
- Increase `mkryl` (larger Krylov space): `L.P.mkryl = 300`
- Allow non-converged solutions: `L.P.keepok = true; L.P.tolok = 1e-4`
- Use better initial guess: `L.P.initguess = {'vacuum', 'extrap', 'previous'}`
- Use tighter preconditioner: `L.P.ilu_droptol = 1e-5`
- Check for unphysical inputs (zero or negative `signeo`, unreasonable `Rp`)

### 6.3 "Conductivity of plasma has to be strictly positive"

Cause: `LX.signeo` contains zero or negative values.

**Fix**: Ensure `signeo` has a floor value:
```matlab
LX.signeo = max(LX.signeo, 1e3);  % floor at 1000 S/m
```

### 6.4 "Multiple magnetic axes but not doublet"

Cause: The nonlinear solver found a configuration with multiple O-points but `idoublet=false`.

**Fixes**:
- For a doublet plasma: set `P.idoublet = true`
- For an aberrant solution: check the initial condition and reduce `dt`
- Tighten convergence: `L.P.tolF = 1e-10`

### 6.5 Vertical Displacement Event (VDE)

Symptom: `LY.zA` rapidly drifts to the wall, followed by non-convergence.

The vertical growth rate for TCV is ~1000/s. Without vertical control, the plasma hits the wall in ~5 ms.

**Fixes**:
- Enable a vertical stabilization controller via `L.P.ctrlfct`
- Use a smaller `dt` to capture the instability dynamics
- Pre-compute the unstable mode: `[Ei, Vi] = fgeeig(L); fprintf('Growth rate: %.1f/s\n', real(Ei(1)))`
- For linearized simulation only: ensure the linear model includes the unstable eigenvalue correctly

### 6.6 "times must be equally spaced when passing multiple LX slices"

Cause: `LX.t` is not uniformly spaced.

**Fix**: Ensure `t = t0:dt:tend` with constant `dt`.

### 6.7 Linearization Fails for Some Time Slices

Symptom: `fgelinall` reports failed linearizations.

**Fix**: Use robust mode to skip failures:
```matlab
[Llin, LYlin, mask] = fgelinall(Llin, L, LY, 'robust', true);
fprintf('%d of %d linearizations succeeded\n', sum(mask), numel(mask));
```

Common cause: equilibrium with VDE or non-converged Newton at those time slices.

### 6.8 Memory Issues with Large Vessel Representation

Cause: `selu='v'` with TCV full vessel (nu~256) creates large matrices.

**Fix**: Use eigenmode truncation:
```matlab
[L, LX] = fge(tok, shot, t, 'selu', 'e', 'nu', 50);
% Uses only 50 dominant vessel eddy current eigenmodes
```

Check accuracy:
```matlab
[L_full, ~, LY_full] = fge(tok, shot, t, 'selu', 'v');
[L_eig,  ~, LY_eig ] = fge(tok, shot, t, 'selu', 'e', 'nu', 50);
% Compare LY_full.zA vs LY_eig.zA
```

### 6.9 "rzp is only a linear model"

Cause: Attempting to run nonlinear simulation with `code='rzp'`.

The RZIP model is purely linear; it can only be run with `P.lin=true`. FGE (code='fge') supports both linear and nonlinear modes.

### 6.10 Iterator Mode Errors

FGE supports an "iterator" mode for external control loops, where `fget` is called one step at a time. The iterator must pass back all state variables between calls:

```matlab
%% Init
[LYt, LYDt, xnl, xnldot, Prec, psstate, dstate, Tstate, nnoc] = ...
    fget(L, meqxk(LX, 1), 'dt', dt);

%% Each subsequent step
for ii = 2:nt
    [LYt, LYDt, xnl, xnldot, Prec, psstate, dstate, Tstate, nnoc] = ...
        fget(L, meqxk(LX, ii), ...
             'LYp', LYt, 'xnl', xnl, 'xnldot', xnldot, ...
             'Prec', Prec, 'dt', dt, 'psstate', psstate, ...
             'dstate', dstate, 'Tstate', Tstate, 'nnoc', nnoc);
    if nnoc > L.P.nnoc; break; end
end
```

Missing any of these state variables will cause errors or incorrect results.

---

## 7. Performance Tips

### 7.1 Choose the Right algoNL

| algoNL | Cost per step | Notes |
|---|---|---|
| `'all-nl'` | Medium | Default, robust, state is Iy |
| `'all-nl-Fx'` | Medium | State is Fx, similar cost |
| `'Newton-GS'` | Higher setup, faster iterations | Sparse Jacobian possible |
| `'Picard'` | Slow | No analytical Jacobian, avoid for FGE |

For TCV with a preconditioner: `'all-nl'` with `anajac=true` is typically fastest.

### 7.2 Vessel Modes

- Full vessel (`selu='v'`, `nu=256`): accurate, slow (256×256 Muu matrix)
- Eigenmode (`selu='e'`, `nu=20-50`): fast, captures 99% of dynamics for `dt > 0.1 ms`
- Single mode (`selu='e'`, `nu=1`): ultrafast, for rapid parameter sweeps

### 7.3 Preconditioner Update Frequency

```matlab
L.P.precupdate = 'once';   % compute once at start (default)
L.P.precupdate = 'time';   % update every time step (expensive but robust)
L.P.precupdate = 'newton'; % update every Newton iteration (very expensive)
```

For slowly-varying equilibria, `'once'` is usually sufficient.

### 7.4 Debug Levels

```matlab
L.P.debug = 0;  % Progress bar only (default)
L.P.debug = 1;  % Problem info + per-step summary
L.P.debug = 2;  % Detailed solver info per step
L.P.debug = 3;  % Full parameter display
L.P.debugplot = 1;  % Plot equilibrium per time step (very slow)
```

---

## 8. Typical Parameter Choices for TCV

```matlab
%% Standard TCV FGE run (open-loop, no CDE)
[L, LX] = fge('tcv', shot, t, ...
    'insrc',   'liu', ...      % Initial cond from LIUQE
    'selu',    'e',   ...      % Eigenmode vessel representation
    'nu',      50,    ...      % 50 vessel modes
    'algoNL',  'all-nl', ...   % Standard NL formulation
    'anajac',  true,  ...      % Analytical Jacobian
    'tolF',    1e-8,  ...      % Residual tolerance
    'mkryl',   200,   ...      % GMRES Krylov dimension
    'ssinit',  true,  ...      % Steady-state initialization
    'agcon',   {'Ip', 'bp', 'qA'});  % Profile constraints

%% TCV FGE run with 0D CDE
[L, LX] = fge('tcv', shot, t, ...
    'insrc',   'liu', ...
    'selu',    'e',   ...
    'nu',      50,    ...
    'cde',     'OhmTor_rigid_0D', ...
    'agcon',   {'bp', 'qA'});  % Ip replaced by CDE

%% For linearization / control design
[L, LX] = fge('tcv', shot, t, ...
    'insrc',   'liu', ...
    'selu',    'e',   ...
    'nu',      20,    ...      % Fewer modes for smaller state space
    'cde',     'OhmTor_rigid_0D', ...
    'agcon',   {'bp', 'qA'}, ...
    'lin',     true);          % Will run linearized model
```
