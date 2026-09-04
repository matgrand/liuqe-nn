# LIUQE Practical Usage Guide

Practical guide for using LIUQE in the MEQ (MatlabEQuilibrium) toolbox.
Companion to `liuqe_report.md` (technical reference).

---

## 1. Minimal Working Example

### 1.1 Basic TCV Reconstruction

```matlab
% Reconstruct a single time slice
shot = 65000;
t    = 0.5;       % [s]

[L, LX, LY] = liuqe(shot, t);
```

`liuqe` is simply a wrapper for `liu('TCV', shot, t)`. The three outputs:
- `L`: assembled parameter/geometry structure (can be reused for other slices)
- `LX`: loaded measurement data
- `LY`: reconstructed equilibrium

### 1.2 Multiple Time Slices

```matlab
shot = 65000;
t    = 0.3:0.01:1.0;    % 70 time slices

[L, LX, LY] = liuqe(shot, t);
% LY is now a structure where each field has an extra time dimension
% LY.Ip  → [1, 70] plasma current time trace
% LY.Fx  → [32, 27, 70] poloidal flux map at each time
```

### 1.3 Geometry-Only Call (No Reconstruction)

```matlab
% Get only L (parameters + geometry), without loading measurements
L = liuqe(shot);          % equivalent to liuqe(shot, [])
% or equivalently:
L = liu('TCV', shot, []);
```

Useful when you want to reuse the geometry for many shots/time slices.

### 1.4 Separate Geometry Loading and Reconstruction

```matlab
% Load geometry once (expensive: reads MDSplus static tree)
shot = 65000;
[L, LX] = liuqe(shot, t);

% Run reconstruction with different parameters without reloading geometry
LY1 = liut(L, LX, 'wreg', 1e-5);
LY2 = liut(L, LX, 'wreg', 1e-7);

% Or for a new shot with the same geometry (test only - geometry may differ):
[~, LX2] = liuqe(shot2, t2);
LY3 = liut(L, LX2);
```

### 1.5 Parameter Overrides

Any parameter from `liup.m` or `liuptcv.m` can be overridden:

```matlab
[L, LX, LY] = liuqe(shot, t, ...
    'wreg',    1e-5,  ...   % stronger regularization
    'itera',   100,   ...   % more iterations
    'psichco', 1e-5,  ...   % tighter convergence
    'bfp',     [2,3], ...   % more basis functions (5 total)
    'wFf',     ones(38,1), ... % restore flux loop 10-12
    'debug',   1);          % print iteration info
```

---

## 2. Setting Up Geometry (liug / liugtcv)

### 2.1 What Geometry Setup Does

The geometry setup is the expensive step (~seconds) that reads electromagnetic
coupling matrices from the TCV MDSplus static database. It is performed automatically
by `liuqe` / `liu` and results are cached in `L.G`.

For most users, the geometry is loaded transparently. The geometry depends on:
- Shot number (to select correct machine configuration, e.g., coil geometry versions)
- Parameters `P.sela`, `P.selu`, `P.selc`, `P.selx` (selection of coil, vessel,
  probe, and flux-surface configurations)
- `P.nz`, `P.nr` (grid size)

### 2.2 Reusing Geometry Across Shots

For efficiency when processing many shots from the same TCV period:

```matlab
% Load geometry from a reference shot
shotg = 50000;
L = liuqe(shotg, []);   % geometry only

% Process many shots efficiently
shots = [65000, 65001, 65002, 65003];
times = [0.5,   0.5,   0.5,   0.5  ];

for k = 1:numel(shots)
    % Load only measurements (fast)
    LX = liuxtcv(shots(k), times(k), L);
    LX = liux(L, LX);
    % Reconstruct
    LY = liut(L, LX);
    % ... process LY ...
end
```

**Caution**: the geometry (coil positions, probe locations, static matrices) can
change between TCV configuration periods. Using a reference shot from the same
campaign period is safer than always using shot 50000.

### 2.3 Key Geometry Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `P.nr`  | 27 | Number of R points in reconstruction grid |
| `P.nz`  | 32 | Number of Z points (must be power of 2) |
| `P.sela`| (from meqptcv) | Active coil selection string |
| `P.selu`| `'sx'` | Vessel representation: 'sx', 's', 'e' |
| `P.selx`| `'X'` | Flux surface selection for DML |
| `P.selc`| (from meqptcv) | Conductor selection |
| `P.idml`| 1 | Include diamagnetic flux measurement |

### 2.4 Inspecting the Geometry

```matlab
L = liuqe(shot, []);

% Grid
figure; imagesc(L.rx, L.zx, zeros(L.P.nz, L.P.nr));
hold on; plot(L.G.rl, L.G.zl, 'k-');  % limiter

% Coil positions
plot(L.G.ra, L.G.za, 'r*');           % active coils
plot(L.G.rm, L.G.zm, 'g+');           % magnetic probes
plot(L.G.rf, L.G.zf, 'bs');           % flux loops

fprintf('Grid: %d × %d\n', L.P.nz, L.P.nr);
fprintf('Coils: %d active, %d vessel segments\n', L.G.na, L.G.nu);
fprintf('Diagnostics: %d flux loops, %d probes\n', L.G.nf, L.G.nm);
```

For TCV standard configuration: 19 active coils, 38 vessel segments, 38 flux loops,
38 magnetic probes.

---

## 3. Loading Measurements (liux / liuxtcv)

### 3.1 Automatic Loading

Measurements are loaded automatically when you call `liuqe` or `liu` with a non-empty
time array. They are returned in `LX`:

```matlab
[L, LX] = liuqe(shot, t);
% LX.Ff    [38 × nt]  flux loop values [Wb]
% LX.Bm    [38 × nt]  probe values [T]
% LX.Ia    [19 × nt]  coil currents [A]
% LX.Iu    [38 × nt]  vessel currents [A]
% LX.Ip    [1  × nt]  plasma current [A]
% LX.Ft    [1  × nt]  diamagnetic flux [Wb]
% LX.Xt    [1  × nt]  diamagnetic flux (before Iy correction) [Wb]
% LX.rBt   [1  × nt]  R × Bt [m·T]
% LX.t     [1  × nt]  actual measurement times [s]
```

Note: `LX.t` may differ slightly from the requested `t` because measurements are
rounded to the nearest magnetics sample (100 µs grid for TCV).

### 3.2 Manual Measurement Loading

```matlab
% Load measurements separately
[L] = liuqe(shot, []);             % geometry only
LX = liuxtcv(shot, t, L);         % load from MDSplus
LX = liux(L, LX);                 % validate and augment

% Then reconstruct
LY = liut(L, LX);
```

### 3.3 Checking Measurement Data

```matlab
% Quick sanity check plots
figure;
subplot(3,1,1); plot(LX.t, LX.Ip/1e3); ylabel('Ip [kA]');
subplot(3,1,2); plot(LX.t, LX.Ff.');   ylabel('Flux loops [Wb]');
subplot(3,1,3); plot(LX.t, LX.Bm.');   ylabel('Probes [T]');
```

### 3.4 Real-Time Data (liuxrtc)

For shots where RT-LIUQE data was stored in the RTC tree MEMs:

```matlab
% For shots >= 65300 (Mark 2 SCD system):
L.P.xrtc = 1;                         % Enable RTC data source
[L, LX] = liuqe(shot, t);             % will call liuxrtc internally

% Alternatively:
LX = liuxrtc(shot, t, L);
```

### 3.5 Time Selection for Full Shot

The standard analysis time selection (matching `liuanasrv.m`):

```matlab
t = liuqe_times(shot);        % default: 'mixed' (1ms grid + Thomson times)
t = liuqe_times(shot, 'fast');     % 1ms grid
t = liuqe_times(shot, 'thomson'); % Thomson scattering times only
t = liuqe_times(shot, 'magnetics'); % all magnetics samples (10 kHz)
```

---

## 4. Understanding the Output LY Structure

### 4.1 Accessing Single vs. Multiple Time Slices

When multiple time slices are reconstructed, each field in `LY` has a trailing time
dimension:

```matlab
% Single slice:
LY = liut(L, meqxk(LX, 1));   % process only first slice
LY.Ip                          % scalar
LY.Fx                          % [nz × nr] matrix

% Multiple slices:
[L, LX, LY] = liuqe(shot, t);
LY.Ip                          % [1 × nt]
LY.Fx                          % [nz × nr × nt]
LY.Ia                          % [na × nt]
```

Use `meqxk(LY, kt)` to extract a single time slice from a multi-slice LY:
```matlab
LYt = meqxk(LY, 5);   % 5th time slice
```

### 4.2 Core Equilibrium Results

```matlab
% Poloidal flux map
psi = LY.Fx;        % [Wb/rad], (nz × nr) at each time

% Magnetic axis
R_axis = LY.rA;     % [m]
Z_axis = LY.zA;     % [m]
psi_axis = LY.FA;   % [Wb/rad]

% Boundary
psi_boundary = LY.FB;
is_diverted = LY.lB;  % 0=limited, 1=X-point bounded

% X-points (may be empty for limited plasmas)
R_xpoint = LY.rX;    % [m], multiple columns for multiple X-points
Z_xpoint = LY.zX;

% Plasma current distribution
Jtoroid = LY.Iy;    % [ny_z × ny_r × nt] toroidal current density
```

### 4.3 Profile Results

```matlab
% Basis function coefficients (raw output of reconstruction)
ag = LY.ag;         % [ng × nt], ng=3 for default bfp=[1,2]

% Safety factor profile
q = LY.aq;          % [noq × nt], noq=128 by default for TCV

% To get q profile on normalized flux grid:
psi_N = (0:noq-1)' / (noq-1);  % normalized flux ∈ [0,1]
plot(psi_N, LY.aq(:, kt));      % q profile at time kt

% Toroidal field profile T(ψ) = R*Bφ
T_prof = LY.TQ;   % [nthg × nt], nthg=129 default

% Global scalar parameters
beta  = LY.beta;   % toroidal beta
betap = LY.betap;  % poloidal beta
li    = LY.li;     % internal inductance li/2
Wp    = LY.Wp;     % stored thermal energy [J]
Vp    = LY.Vp;     % plasma volume [m³]
```

### 4.4 Fit Quality Assessment

```matlab
% Overall fit quality
chi = LY.chi;       % RMS of weighted residuals (should be ~1 for good fit)
conv = LY.isconverged;  % boolean: did reconstruction converge?

% Per-measurement residuals
% LY.werr: weighted residuals [nd × nt]
% nd = nf + nm + na + nu + 1 + 1 (flux loops + probes + coils + vessel + DML + Ip)
% Ordering: [flux loops, probes, coils, vessel, DML, Ip, regularization, constraints]

% Per-diagnostic comparison
% Reconstructed measurements:
Ff_fit = LY.Ff;     % [nf × nt] fitted flux loop values
Bm_fit = LY.Bm;     % [nm × nt] fitted probe values

% Compare with measurements:
figure;
subplot(2,1,1);
plot(LX.Ff - LY.Ff);    % flux loop residuals
title('Flux loop residuals [Wb]');
subplot(2,1,2);
plot(LX.Bm - LY.Bm);    % probe residuals
title('Probe residuals [T]');
```

### 4.5 Distinguishing Converged from Non-Converged Slices

By default, only converged slices are in `LY` (non-converged are removed). To include
all slices:

```matlab
[L, LX, LY] = liuqe(shot, t, 'LYall', true);
% LY now includes non-converged slices too
% Check convergence per slice:
mask = LY.isconverged;
```

---

## 5. Key Parameters to Tune

### 5.1 Basis Function Choice

The biggest impact parameter is `bfp` (and `bfct`). More basis functions = more
flexibility in the current profile, but risk of overfitting and oscillations.

```matlab
% Default TCV: 1 for p', 2 for TT' = 3 total
[~,~,LY] = liuqe(shot, t, 'bfp', [1,2]);

% More flexible: 2 for p', 3 for TT' = 5 total
% Use with caution: may not converge for weak plasmas
[~,~,LY] = liuqe(shot, t, 'bfp', [2,3]);

% Minimum (RT-like): 0 for p', 2 for TT' = 2 total (like bf3pmex)
[~,~,LY] = liuqe(shot, t, 'bfp', [0,2]);
```

### 5.2 Regularization

`wreg` controls smoothing of the current profile. Larger wreg = smoother profile but
potentially more biased.

```matlab
% Standard
LY1 = liut(L, LX, 'wreg', 2e-6);

% Stronger regularization (smoother, for noisy shots)
LY2 = liut(L, LX, 'wreg', 1e-5);

% No regularization (can produce oscillatory profiles)
LY3 = liut(L, LX, 'wreg', 0);
```

The adaptive regularization: `wregadapt=15` causes wreg to decrease for elongated
plasmas (κ > elomin=1.1). This makes sense as elongated plasmas are more uniquely
determined by the flux surface shape.

### 5.3 Measurement Weights

Turn off specific measurements (e.g., to test reconstruction quality):

```matlab
% Disable all probes
LY_noB = liut(L, LX, 'wBm', zeros(38,1));

% Disable DML (diamagnetic flux)
LY_nodml = liut(L, LX, 'idml', 0);

% Disable vessel currents
LY_noves = liut(L, LX, 'ivesm', 0);

% Increase weight of plasma current measurement
LY_highIp = liut(L, LX, 'wIp', 10);

% Disable specific flux loops (e.g., loops 5 and 6 if faulty)
wFf = ones(38,1); wFf([5,6]) = 0;
LY_masked = liut(L, LX, 'wFf', wFf);
```

### 5.4 Convergence Control

```matlab
% More iterations (for difficult shots)
LY = liut(L, LX, 'itera', 100, 'psichco', 1e-6);

% Accept less-converged result (faster)
LY = liut(L, LX, 'itera', 20, 'psichco', 1e-4);

% Relax chi acceptance threshold
LY = liut(L, LX, 'chimax', 5);
```

### 5.5 Vertical Stabilization

```matlab
% Default: dz is a free parameter (passive mode)
% L.P.zfct = @(P,rst,z0,zA,dz) dz;

% Enable PI feedback on fitted shift
[L,LX,LY] = liuqe(shot, t, 'zfct', @liudz, 'stab', 1, 'stabint', 1);

% Enable PID feedback on axis position
[L,LX,LY] = liuqe(shot, t, 'zfct', @liuz0, 'z0', 0.0, ...
    'stab', 1, 'stabint', 0.1, 'stabder', 0.01);

% Disable vertical stabilization entirely
[L,LX,LY] = liuqe(shot, t, 'stabz', false);
```

---

## 6. Common Failure Modes and Diagnostics

### 6.1 Reconstruction Fails to Converge

**Symptom**: `LY.isconverged` is false; `LY` is empty or contains few slices.

**Diagnosis**:
```matlab
% Run with debug output
LY = liut(L, LX, 'debug', 1);

% Check chi per time slice
plot(LY.chi);
yline(L.P.chimax);
title('Chi (dashed = chimax threshold)');

% Check number of iterations
plot(LY.niter);
title('Iterations to convergence');
```

**Common causes and fixes**:

| Cause | Symptom | Fix |
|-------|---------|-----|
| Low plasma current | Warning 'LowIp' | Check `Ipmin` setting; shot may be pre-breakdown |
| Bad measurements | High chi on specific channels | Check `LY.werr`, disable bad channels |
| FE fit failure | `chih` is large | Increase `iterh`, adjust `nelem` |
| Wrong vertical position | Large `dzg` | Enable vertical stabilization (`zfct`) |
| Profile oscillation | `wreg` too small | Increase `wreg` |
| Plasma hits limiter | lB switches between 0 and 1 | Check shot — may be a fast event |

### 6.2 Systematic Bias in Flux Loops

**Symptom**: Reconstructed flux loops (`LY.Ff`) consistently offset from measured
(`LX.Ff`) for specific channels.

**Diagnosis**:
```matlab
% Per-channel flux loop residuals at one time
kt = 5;
bar(LX.Ff(:,kt) - LY.Ff(:,kt));
title('Flux loop residuals [Wb]');
```

**Fixes**:
- If one channel is always a large outlier: disable it (`wFf(k) = 0`)
- If all channels are offset: check DML correction (try `idml=0`)
- Check measurement calibration offsets in the static tree

### 6.3 Wrong Current Profile (q(0) ≠ 1)

**Symptom**: `LY.qA` (safety factor at axis) is far from expected values.

This is expected for magnetics-only reconstruction — q(0) is poorly constrained
by external magnetic measurements alone. To constrain it:
- Add Faraday rotation measurements (requires `ifara=1`, not supported in standard MEQ)
- Use SXR or Thomson scattering constraints (external post-processing)
- Use a prior (regularization toward expected shape)

The default polynomial basis functions (`bfp=[1,2]`) produce smooth, monotonically
increasing q profiles. More complex q profiles require more basis functions.

### 6.4 Plasma Current Doesn't Match

**Symptom**: `|LY.Ip - LX.Ip| / LX.Ip > 10%`.

This is checked automatically in `liutcv_tests.m` (tolerance 10%). Causes:
- Rogowski coil miscalibration (affects `LX.Ip`)
- Too-strong vessel current regularization or wrong `Iuerr`
- Poor probe weights (`wBm`) affecting the plasma current estimator

```matlab
% Check Ip balance
figure;
plot(LX.t, LX.Ip/1e3, 'b-', LX.t, LY.Ip/1e3, 'r--');
legend('Measured Ip', 'Fitted Ip');
ylabel('kA');
```

### 6.5 Runtime Performance

```matlab
% Check computation time per time slice
plot(LY.cycle / 1e3); ylabel('ms'); title('Cycle time');

% If too slow, try:
% 1. Reduce basis functions
LY = liut(L, LX, 'bfp', [1,1]);   % 2 basis functions

% 2. Fewer iterations
LY = liut(L, LX, 'itera', 20);

% 3. Warm-start from previous (for time series)
LY = liut(L, LX, 'initguess', 'previous');

% 4. Parallelise time slices (see liuanasrv.m example)
parfor k = 1:np
    LY_{k} = liut(L, LX_{k});
end
LY = meqcat(LY_{:});
```

### 6.6 Shot Exclusion List

Some TCV shots have known issues and should not be used for testing:
- Shots 40042, 57942: faulty magnetic probe
- Shot 49823: stray shot  
- Shot 54001: DML issues at early times
- Shots 50000–50606: faulty magnetics sensors

See `liutcv_tests.get_excludelist()` for the current list.

---

## 7. Advanced Usage

### 7.1 Using LIUQE as a Real-Time Iterator

```matlab
% Configure for real-time mode
L.P.itert = 1;         % 1 iteration per step
L.P.initguess = 'previous';

% First step (with full FE fit)
LY1 = liut(L, meqxk(LX, 1));

% Subsequent steps (warm-start from previous)
x = L.LX2x(LY1);
u = [LY1.Ia; LY1.Iu; LX2Co(L, LY1)];
z = LY1.dz;

LY2 = liut(L, meqxk(LX, 2), 'x', x, 'u', u, 'z', z, 'rst', false);
```

### 7.2 Custom Base Functions

Using `bf3imex` / `bfli` with user-defined profiles:

```matlab
% Define 3 custom basis functions on a 41-point grid
n = 41;
FN = linspace(0,1,n)';

% Example: flat, peaked, hollow
gNg(:,1) = ones(n,1) .* (1-FN);      % linearly decreasing
gNg(:,2) = exp(-3*FN);                % peaked at axis
gNg(:,3) = FN.*(1-FN);               % hollow (zero at both ends)

% Compute primitives (use bfprmex or trapezoidal integration)
% Here using simple trapz:
IgNg = -cumtrapz(FN, gNg);            % note: integral from FN to 1

% Set up basis function parameters
bfp = struct('gNg', gNg, 'IgNg', IgNg, ...
             'fPg', [true false false], ...  % 1st function is p'
             'fTg', [false true true]);      % 2nd and 3rd are TT'

[~,~,LY] = liuqe(shot, t, 'bfct', @bfli, 'bfp', bfp);
```

### 7.3 Doublet Configuration

For TCV double-null divertor plasmas with two plasma domains:

```matlab
[~,~,LY] = liuqe(shot, t, 'idoublet', 1);
% LY.rA   [2 × nt]  R of both magnetic axes
% LY.zA   [2 × nt]  Z of both axes
% LY.IpD  [2 × nt]  current in each domain
```

### 7.4 Writing Results to TCV Tree

The standard production pipeline writes results using:
```matlab
% After running liuanasrv (which calls these internally):
meqwtcv(shot, 'RESULTS', '\TOP.EQUIL_1.RESULTS', L, LY, dep)
liuwtcv(shot, 'RESULTS', '\TOP.EQUIL_1.RESULTS', L, LY, dep)
```

`liuwtcv` writes LIUQE-specific data (currently just `DELTA_Z` = the vertical shift).
`meqwtcv` writes the full equilibrium output (Fx, aq, beta, li, etc.).

### 7.5 Comparing LIUQE Runs

```matlab
% Run with two different settings
[L, LX, LY1] = liuqe(shot, t, 'wreg', 2e-6);
LY2 = liut(L, LX, 'wreg', 1e-5);

% Compare q profiles at one time
kt = 10;
psiN = linspace(0, 1, L.P.noq);
figure;
plot(psiN, LY1.aq(:,kt), psiN, LY2.aq(:,kt));
legend('wreg=2e-6', 'wreg=1e-5');
xlabel('\psi_N'); ylabel('q');

% Compare plasma shapes
figure;
contour(L.rx, L.zx, LY1.Fx(:,:,kt), [LY1.FB(kt), LY1.FA(kt)], 'b');
hold on;
contour(L.rx, L.zx, LY2.Fx(:,:,kt), [LY2.FB(kt), LY2.FA(kt)], 'r');
plot(L.G.rl, L.G.zl, 'k-');
legend('wreg=2e-6', 'wreg=1e-5', 'limiter');
```

---

## 8. Function Quick Reference

| Function | Purpose |
|----------|---------|
| `liuqe(shot,t,...)` | TCV top-level: runs full reconstruction |
| `liu('TCV',shot,t,...)` | Generic top-level wrapper |
| `liuptcv(shot,...)` | Load TCV parameters from MDSplus |
| `liugtcv(shot,P)` | Load TCV geometry from MDSplus |
| `liuxtcv(shot,t,L)` | Load TCV measurements from MDSplus |
| `liup(P)` | Validate and fill defaults in P |
| `liug(G,P)` | Validate and fill defaults in G |
| `liux(L,LX)` | Validate and augment measurement structure |
| `liuc(P,G)` | Assemble full L structure (ancillary data) |
| `liut(L,LX,...)` | Run reconstruction time loop |
| `liuanasrv(shot,...)`| Production server: full shot analysis + write |
| `liusimc(L,LX,...)`| Prepare Simulink real-time model |
| `liuqe_times(shot,...)`| Get standard time array for a shot |
| `liuwtcv(shot,tree,branch,...)` | Write LIUQE results to MDS tree |
| `meqxk(LX,k)` | Extract k-th time slice from LX or LY |
| `meqcat(LY1,LY2,...)` | Concatenate multiple LY structures |

---

## 9. Typical Workflow for a New Shot Analysis

```matlab
%% 1. Quick look
shot = 65000;
[L, LX, LY] = liuqe(shot, 0.5);
fprintf('Ip = %.0f kA, chi = %.2f\n', LY.Ip/1e3, LY.chi);

%% 2. Check measurement quality
figure;
subplot(2,1,1); bar(LX.Ff - LY.Ff); title('Flux loop residuals [Wb]');
subplot(2,1,2); bar(LX.Bm - LY.Bm); title('Probe residuals [T]');

%% 3. Full shot reconstruction
t = liuqe_times(shot, 'fast');    % 1ms grid + Thomson times
[L, LX, LY] = liuqe(shot, t);

fprintf('%d / %d slices converged\n', sum(LY.isconverged), numel(t));

%% 4. Plot key traces
figure;
subplot(4,1,1); plot(LY.t, LY.Ip/1e3);    ylabel('Ip [kA]');
subplot(4,1,2); plot(LY.t, LY.kappa);      ylabel('\kappa');
subplot(4,1,3); plot(LY.t, LY.qA);         ylabel('q(0)');
subplot(4,1,4); plot(LY.t, LY.chi);        ylabel('\chi');
```
