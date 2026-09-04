# Hybrid Controller Practical Guide
## TCV Tokamak — Swiss Plasma Center, EPFL

**Repository:** `/home/grandin/repos/hybrid`  
**Companion document:** `hybrid_report.md` (architecture reference)

---

## 1. How to Initialize the Hybrid Model for a TCV Shot

### 1.1 Minimal Setup

The entry point for all offline controller work is `init_hybrid_model`. It requires a **Shape Design Parameters** (SDP) structure specifying the shot and time:

```matlab
% Add repo to path first
addpath(genpath('/home/grandin/repos/hybrid'))

shot = 78071;   % reference shot number
t    = 0.6;     % time instant (seconds) within the flat-top

% Step 1: Get default shape design parameters
FP  = get_shape_fixed_params();               % fixed timing, dimensions
SDP = get_shape_design_params('shot', shot, 't', t, 'ncC', 7);

% Step 2: (Optional) get shape control points from FBT
SC = get_shape_control_points(FP, SDP);

% Step 3: Build the linearized model + hybrid controller
[L, cpar, ssHyb, ssA, ssG, ssPlant, ssContr] = init_hybrid_model(SDP, SC);
```

After this:
- `L` is the FGEL equilibrium structure (plasma geometry, coil data, linearized model)
- `cpar` is the hybrid controller parameter structure with A, G, M matrices
- `ssHyb` is the closed-loop state-space model (plasma + controller)
- `ssA`, `ssG` are the A and G matrix state-space blocks (for named interconnection)
- `ssPlant`, `ssContr` are the open-loop plant and controller

### 1.2 The SDP Structure

`get_shape_design_params` accepts many optional key-value pairs. The most important:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `shot` | required | TCV shot number |
| `t` | required | Equilibrium time instant |
| `ncC` | 7 | Number of core boundary control points |
| `code` | `'fbt'` | Equilibrium code to use (`'fbt'`, `'mag'`) |
| `svdtol` | 4 | SVD truncation tolerance for shape control |
| `wcross` | 20 | Shape controller crossover frequency (Hz) |
| `PM` | 65 | Phase margin (degrees) |
| `wC` | 1 | Weight on core boundary control points |
| `wX1` | 1 | Weight on first X-point control points |
| `wL` | 0 | Weight on limiter points |
| `tuning` | `'pidtune'` | Controller tuning method: `'pidtune'`, `'mpc'` |
| `SCsrc` | `'LY'` | Shape control point source: `'LY'`, `'rz'`, `'LX'` |
| `fixip` | `false` | Fix plasma current in equilibrium (vacuum model) |
| `bp` | 0 | Poloidal beta to use in linearization |

For shots before shot #900000, the data is read from the `tcv_shot` MDS tree; for later shots from `pcs`.

### 1.3 What `init_hybrid_model` Does Internally

1. Calls `init_fgel_model(shot, t, fixip, code)` → gets the linearized equilibrium `L`
2. Calls `hybctrl('setup', [], [], L, [], 't', t, 'ton', t)` → initializes `cpar` with A, G, M matrices read from MDS
3. Optionally scales gains via `scale_hybrid_gains(cpar, t, rgain, zgain)`
4. Calls `build_hybrid_linmod(L, cpar, debug, SC, ...)` → builds all state-space blocks
5. Closes the feedback loop via `close_identity_feedback` → returns `ssHyb`

### 1.4 Vacuum (No-Plasma) Initialization

For testing without a plasma, use `fixip = true` or use a vacuum shot:

```matlab
SDP = get_shape_design_params('shot', shot, 't', t, 'fixip', true, 'code', 'mag');
[L, cpar, ...] = init_hybrid_model(SDP, []);
```

This uses only the magnetic (vacuum) model without plasma response.

---

## 2. How to Run a Closed-Loop Simulation

### 2.1 Using the Shape Controller (Full System)

The primary simulation script is `manual_tests/test_closed_shape.m`. The typical workflow:

```matlab
% Step 1: Run main_shape_ctrl.m to generate controller parameters
% (or load from file if already computed)
run('/home/grandin/repos/hybrid/scripts/main_shape_ctrl.m')
% This produces: SC, STP, ATP, MTP, LO structures

% OR load previously saved parameters:
load('/home/grandin/repos/hybrid/local_data/shapecontr', ...
     'shot','FP','SDP','ADP','STP','MTP','ATP','SC','LO','rgain','zgain');

% Step 2: Set up simulation time vector
dt   = 1e-4;
tsim = STP.ton(1)-0.05 : dt : STP.ton(1)+0.3;
tsim = round(tsim, 5);  % avoid floating-point issues

% Step 3: Load FGE model and reference data
[L, LX] = init_fgel_model(shot, tsim(1), SDP.fixip, SDP.code{1});
LX = fgex([], tsim, L, LX);  % expand to full time vector

% Step 4: Get LIUQE parameters for resistance values (needed for NN path)
Lliu = liu('tcv', shot, SDP.t);

% Step 5: Set up controller
SCref = SC;  % or modify to change references
L.P.ctrlfct = @shapectrl;
[~,~,cpar,~] = L.P.ctrlfct('setup', [], [], L, [], ...
    't', tsim, 'dt', dt, 'Gdead', 0, ...
    'SC', SCref, 'FP', FP, 'STP', STP, 'ATP', ATP, 'MTP', MTP, 'LO', LO, ...
    'Rs', Lliu.G.Rs, 'useliunn', false);  % set useliunn=true to use NN

% Step 6: Run FGE simulation
L.P.ctrlpar     = cpar;
L.P.debug       = 0;
L.P.stoponalarm = 0;
L.P.lin         = true;   % true = linear, false = nonlinear

LY = fget(L, LX, 'debugplot', 0);

% Step 7: Plot results
shape_ctrl_plot(L, LY, SCref, cpar, 'shape_0');
```

### 2.2 Using Just the Hybrid (RZIP) Controller

For testing the inner loop only, without shape control:

```matlab
[L, LX] = init_fgel_model(shot, tsim(1), 0, 'fbt');
LX = fgex([], tsim, L, LX);

L.P.ctrlfct = @hybctrl;
[~,~,cpar,~] = L.P.ctrlfct('setup', [], [], L, [], 't', tsim, 'Gdead', 0);

L.P.ctrlpar = cpar;
LY = fget(L, LX);
```

### 2.3 Using the FMag Controller

```matlab
[L, LX] = fge('tcv', shot, t_vec, ...
    'cde', 'OhmTor_0D', 'sela', {'E','F','OH','G'}, ...
    'agcon', {'bp','qA'}, 'selu', 'ex', 'nu', 100, 'izgrid', 1, 'nnoc', 2);

FP  = get_fmag_fixed_params();
LX1 = meqxk(LX, 1);  % design snapshot

% Build controller (R/Z/Ip + X-point example)
ctrl1 = build_RZIp_Xpoint_controller(FP, L, LX1, t_vec, 'useXPointControl', true);
[FTP, waveforms] = FMag_generate_mats(FP, {ctrl1}, 0, 1, t_vec, dt);

% Setup simulation
L.P.ctrlfct = @FMag_control;
[~,~,cpar,~] = L.P.ctrlfct('setup', [], [], L, [], 'FTP', FTP, 'waveforms', waveforms);
L.P.ctrlpar = cpar;
LY = fget(L, LX);
```

See `fmag_liuqenn_comparison.m` for a complete example including NN comparison.

### 2.4 Reading Simulation Results

After `fget`, the `LY` structure contains:
- `LY.t`: time vector
- `LY.Ia`: coil currents at each time step
- `LY.Va`: coil voltages
- `LY.Fx`: 2D flux map at each time step (if using nonlinear simulation)
- `LY.ctdebug`: controller debug outputs, including:
  - `LY.ctdebug.shapeErrs`: shape errors at each step
  - `LY.ctdebug.shapeFbck`: shape feedback values
  - `LY.ctdebug.shapeRefs`: shape references
  - `LY.ctdebug.cout`: coil voltage output
  - `LY.ctdebug.Arefs`, `LY.ctdebug.Aout`: RZIP references and measurements

---

## 3. How to Extract the Linear Model (`get_hybrid_linear_model.m`)

The script `scripts/get_hybrid_linear_model.m` demonstrates the full workflow for getting a named state-space model. Here is the condensed version:

```matlab
shot = 78071;
t    = 0.6;
ncC  = 7;

% Get design parameters
FP  = get_shape_fixed_params();
SDP = get_shape_design_params('shot', shot, 't', t, 'ncC', ncC);

% Get shape control points (needed for shape outputs in linear model)
SC = get_shape_control_points(FP, SDP, 1);  % debug=1 shows plot

% Build the linear model
[L, cpar, ssHyb, ssA, ssG, ssPlant, ssContr] = ...
    init_hybrid_model(SDP, SC, 1);  % debug=1 shows step responses
```

### 3.1 The ssHyb Model

`ssHyb` is the closed-loop state-space model from references/disturbances to controlled outputs. Its input/output names can be inspected:

```matlab
ssHyb.InputName   % reference channels (CP013IPA, CP012IP*ZA, etc.)
ssHyb.OutputName  % measured channels (same)
ssHyb.StateName   % state names (coil/structure currents, plasma position)
```

Key input/output names for position channels:
- `'CP013IPA'` or `'IPA'`: plasma current
- `'CP012IP*ZA'` or `'IP*ZA'`: vertical position (Z·Ip)  
- `'CP018_DELTA_R'` or `'DELTA_R'`: radial position

### 3.2 Using the Model for Analysis

```matlab
% Check closed-loop eigenvalues
max_eig = max(real(eig(ssHyb)));
fprintf('Max CL eigenvalue: %.3f\n', max_eig);

% Step responses (position channels)
figure
subplot(311); step(ssHyb('CP013IPA','CP013IPA')); title('Ip step')
subplot(312); step(ssHyb('CP012IP*ZA','CP012IP*ZA')); title('Z·Ip step')
subplot(313); step(ssHyb('CP018_DELTA_R','CP018_DELTA_R')); title('R step')

% Bode plots
bodeopts = bodeoptions; bodeopts.FreqUnits = 'Hz';
bode(ssHyb('CP013IPA','CP013IPA'), bodeopts)

% Shape outputs (if SC was passed)
shape_outputs = ssHyb.OutputGroup.sh;  % indices of shape channels
```

### 3.3 Including Optional Outputs

```matlab
% Include magnetics (Bm, Ff), G coil current, passive currents
[L, cpar, ssHyb, ssA, ssG, ssPlant, ssContr] = ...
    init_hybrid_model(SDP, SC, 0, [], [], ...
    true,   % includeMagnetics
    true,   % includeIG
    true);  % includePassives
```

### 3.4 Modifying Gains Before Building the Model

```matlab
% Scale radial and vertical gains (useful for gain sensitivity analysis)
rgain = 0.5;  % reduce radial gain to 50%
zgain = 1.5;  % increase vertical gain to 150%
[L, cpar, ssHyb, ...] = init_hybrid_model(SDP, SC, 0, rgain, zgain);
```

---

## 4. How to Modify Shape Control Points

### 4.1 Default Control Point Selection

The default behavior selects points using the FBT boundary reconstruction:

```matlab
SDP = get_shape_design_params('shot', shot, 't', t, ...
    'SCsrc', 'LY',  % use FBT boundary for point selection
    'ncC', 7,       % 7 points on core boundary
    'wX1', 1,       % include X-point control
    'wL',  0);      % no limiter points

SC = get_shape_control_points(FP, SDP, 1);  % debug=1 plots result
```

### 4.2 Control Point Sources

Three `SCsrc` options are available:

**`'LY'` (recommended):** Points are placed on the reconstructed boundary from FBT at intervals:
- Core boundary points: ncC points between midplane and X-point
- X-point vicinity points: additional points near the X-point
- Distribution controlled by `wC`, `wX1`, `wX1f`, `wX2` weights

**`'rz'`:** Pure (R,Z) position control — no isoflux points, just the plasma centroid. Use this for plasmas without a well-defined separatrix:
```matlab
SDP = get_shape_design_params('shot', shot, 't', t, 'SCsrc', 'rz');
SC = get_shape_control_points(FP, SDP);
% SC.nc = 0, SC.ref = [rIp; zIp; zeros...]
```

**`'LX'`:** Points based on FBT X-point tracking (least common, used for specific X-point control experiments).

### 4.3 Multi-Snapshot (Time-Varying) Shape Control

For shots where the plasma shape changes significantly during the control window:

```matlab
% Multiple time instants
tshape = [0.6, 0.8, 1.0];  % design times

SDP = get_shape_design_params('shot', shot, 't', tshape, 'ncC', 7, ...
    'activate', [true, true, true]);

% This generates one SC snapshot per time instant
SC = get_shape_control_points(FP, SDP);
% SC.t = tshape
% SC.rc (nc × 3), SC.zc (nc × 3)  ← different points at each time
```

The shape controller will smoothly interpolate between snapshots during the shot.

### 4.4 Fixing Control Points Across Snapshots

To keep the same control point locations across all snapshots (e.g., to test reference blending only):

```matlab
SC = get_shape_control_points(FP, SDP, 0, true);  % fixpoints=true
```

### 4.5 Manually Specifying Points

[GUESS: There is no direct API for manually specifying control points outside of `fbt`/`shapobs`]. To use custom points, you would need to manually construct the SC structure with the required fields:
- `SC.rc`, `SC.zc`: control point coordinates
- `SC.ref`: reference flux values at those points
- `SC.Mysyc`: sensitivity matrix (from `shapobs`)
- `SC.dims`: cell array of channel names
- `SC.nc`, `SC.ns`, `SC.t`

---

## 5. The liuqeml Workflow

### 5.1 Enabling the NN in a Simulation

To enable `liuqeml` inference in a closed-loop simulation with `shapectrl`:

```matlab
% Ensure the liuqe-ml repo is on the path
addpath('/home/grandin/repos/liuqe-ml/')
addpath('/home/grandin/repos/liuqe-ml/pretrained/')

% Set useliunn=true in the shapectrl setup
[~,~,cpar,~] = shapectrl('setup', [], [], L, [], ...
    't', tsim, 'dt', dt, ...
    'SC', SCref, 'FP', FP, 'STP', STP, 'ATP', ATP, 'MTP', MTP, ...
    'Rs', Lliu.G.Rs, ...
    'useliunn', true);   % <-- enable NN
```

The NN will automatically be called instead of `qintmex` in the shape controller step.

### 5.2 The Physics Input Vector

The NN takes a row vector `phys` of size `[1 × nphys]`:

```matlab
% In shapectrl.m (assembled at each step):
phys = [LY.Bm; LY.Ff; Ia; LY.Uf; LY.rBt];
%       ^^^^^^  ^^^^^  ^^  ^^^^^^  ^^^^^^
%       poloidal  flux  EF+G+OH  passive  r*Bt
%       probes   loops  currents  voltages
```

The ordering of `Ia` must match the LIUQE SCD algorithm ordering: `[E coils; F coils; G coil; OH coils]`.

### 5.3 Enabling the NN in FMag (X-point Observer)

```matlab
% In FMag simulation setup:
cpar.useliunn = true;
cpar.rX = rX_design;   % r-coordinates of X-points to query
cpar.zX = zX_design;   % z-coordinates of X-points to query
```

The FMag controller will then call `liuqenn_inference_lih` at each step to get the flux and its derivatives at the X-point positions. See `fmag_liuqenn_comparison.m` for a complete comparison example.

### 5.4 Comparing NN vs. Standard Observer

The `fmag_liuqenn_comparison.m` script runs both paths and saves results:

```matlab
run('/home/grandin/repos/hybrid/fmag_liuqenn_comparison.m')
% Saves to: local_data/liuqenn_comparison/fmag_lih_run.mat
%           local_data/liuqenn_comparison/fmag_lih_run.mat
```

Similarly, for shape control comparison:
```matlab
useliunn = 0;  % standard path
LY_std = fget(L, LX, 'debugplot', 0);

useliunn = 1;  % NN path
[~,~,cpar_nn,~] = shapectrl('setup', [], [], L, [], ..., 'useliunn', true);
L.P.ctrlpar = cpar_nn;
LY_nn = fget(L, LX, 'debugplot', 0);

% Compare shape errors
figure
plot(LY_std.ctdebug.t, LY_std.ctdebug.shapeErrs(1,:), 'b')
hold on
plot(LY_nn.ctdebug.t, LY_nn.ctdebug.shapeErrs(1,:), 'r')
legend('standard (qintmex)', 'NN (liuqeml)')
```

### 5.5 Running the SCDalgo_liuqeml Node

The `liuqeml` algorithm is also available as a standalone SCD node for real-time deployment. It runs at 100 µs and takes the full TCV measurement bus as input. To configure:

```matlab
obj = SCDalgoobj_liuqeml();
% Timing: starts at t=-0.5s, dt=1e-4s
% Input bus: SCDalgo_liuqeml_inBus (Bm, Ff, Ia, Ip, Ft, Iu, rBt, rn, zn)
% Output bus: SCDalgo_liuqeml_outBus (Fn, Brn, Bzn)
```

The input/output bus definitions are in `SCDalgo_liuqeml_loadbus.m`.

### 5.6 NN Weights and Model Versions

The NN weights are loaded from `data/3075260/weights.mat` (inside the `liuqe-ml` repo) at MEX compile time. The weight file contains all network parameters as `single` precision arrays. The architecture is fixed at compile time via `coder.const`. To use a different model version, change the `weights_path` in `net_forward2.m` and recompile the MEX.

### 5.7 MEX Compilation

The `net_forward2.m` function is designed for MATLAB Coder MEX compilation for real-time use:

```matlab
% [GUESS: approximate MEX compilation workflow]
cfg = coder.config('mex');
cfg.GenerateReport = true;
codegen net_forward2 -config cfg -args {coder.typeof(single(0),[1,nphys]), coder.typeof(single(0),[1,npts]), coder.typeof(single(0),[1,npts])}
```

The resulting MEX file (`net_forward2_mex`) is called in `shapectrl.m` as `liuqenn_inference` (which wraps the MEX). The pure MATLAB version (`net_forward2`) can be used for development/debugging.

---

## Troubleshooting

### Common Issues

**"Missing time vector" error in hybctrl setup:**  
Always pass `'t', tsim` to `hybctrl('setup', ...)`. The time vector is required to preload reference waveforms.

**NaN in control action:**  
Usually indicates an instability. Check that the gain matrices are correct for the shot number and time. Try reducing gains with `rgain` and `zgain`.

**Shape errors not converging:**  
- Check that `STP.ton` and `STP.toff` encompass the simulation time
- Check that `SCref.t` aligns with the shape control window
- Try increasing `wcross` (crossover frequency) — default 20 Hz may be too slow

**NN path gives worse results than standard:**  
- Check that `Rs` (passive structure resistances) is correctly passed from `Lliu.G.Rs`
- Check that the control point coordinates `ri`, `zi` are in the correct units (meters) and are all positive `ri > 0`

**MDS connection errors:**  
The controller reads A, G, M matrices from the MDS database. Ensure that `mdsopen` can connect to `pcs` tree for the shot number. For shots ≥ 900000, use `tree = 'pcs'`; for older shots, both `'pcs'` and `'tcv_shot'` may work.

**G matrix not found at the requested time:**  
`build_gmat(cpar.Gmat, t)` selects the G matrix active at time `t`. If `t` is before all scheduled G matrix times, it uses the first one. Check `cpar.Gmat.addr.time` to see the available G matrix time windows.
