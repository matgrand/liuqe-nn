# Hybrid Controller Technical Reference
## TCV Tokamak — Swiss Plasma Center, EPFL

**Repository:** `/home/grandin/repos/hybrid`  
**Purpose:** Plasma shape and position control system for the TCV tokamak  
**Primary language:** MATLAB  
**Real-time target:** Simulink/SCD framework running on TCV PCS (100 µs sampling, 10 kHz)

---

## 1. Overview of the Control Architecture

### 1.1 Hierarchy

The control system is structured in two nested layers:

```
Outer layer (shape/position objectives)
    shapectrl.m  ─── calls ──►  hybctrl.m  (inner layer)
         │                           │
         │ computes Δrefs, Δiarefs,  │  A-matrix observables
         │ Δffs from shape errors    │  G-matrix gains
         │                           │  M-matrix coil allocation
         ▼                           ▼
   Shape controller           Hybrid (RZIP) controller
   (2 ms updates)             (100 µs updates)
         │                           │
         └─────────────►  Power supply voltages  Va → coils
```

**Outer layer — shape/position control** (`shapectrl.m`)  
Runs at 2 ms intervals (every 20 hybrid steps). Computes corrections to the hybrid references (`deltaH`), coil current references (`deltaI`), and voltage feedforwards (`deltaV`) based on the error between measured isoflux shape control points and their references.

**Inner layer — hybrid (RZIP) controller** (`hybctrl.m`)  
Runs at 100 µs. Directly implements the PID loop that converts observable errors (as computed by the A matrix) into coil voltage commands (via the G and M matrices). The term "hybrid" historically refers to the combination of voltage-driven and current-driven coil control modes within the same controller, as well as the combination of different feedback types (R, Z, Ip and OH current).

### 1.2 What "Hybrid" Means

The "hybrid" designation has two distinct meanings in this codebase:

1. **Mixed actuation modes:** Some coils are driven in voltage mode (E1–E8, F1–F8, OH1–OH2) while some are driven in current mode (the G coil, `lhyb` flag). The power supply for current-mode coils runs its own internal PI current-tracking loop (`hybrid_ps_step`), and the hybrid controller bypasses these channels (ensures zero voltage contribution from its feedback path, relying on the PS internal loop instead).

2. **Combined-variable feedback:** The controller simultaneously controls:
   - Plasma current `Ip` (scalar)
   - Vertical position of the current centroid `Z·Ip` (vertical instability)
   - Radial position `ΔR` (isoflux radial displacement)
   - OH coil current difference `ΔI_OH` (inductive drive)
   - Shape control points (isoflux displacements along boundary contour normals)

### 1.3 Role of the Plasma Linear Model

The controller design and (optionally) execution rely on a linearized plasma model derived from a Grad–Shafranov equilibrium reconstruction. This model is described by the RZIP formalism (see Section 2). The key uses are:

- **Offline:** Computing the A, G, M gain matrices for a target shot from an FBT/FGE equilibrium reconstruction (functions `get_amat`, `get_gmat`, `get_mmat`).
- **Offline (design):** Building the full closed-loop linear model for PID tuning or MPC design (`build_hybrid_linmod.m`, `init_hybrid_model.m`).
- **Online (simulation):** The `fgess` / `fge` functions provide the LY structure (current equilibrium state) at each time step during closed-loop simulation.

---

## 2. The RZIP Linear Model

### 2.1 The State-Space Model

The plasma + passive structure + active coil system is described in the RZIP formalism. The governing equation is the circuit equation:

$$L \frac{d\mathbf{I}}{dt} + R \mathbf{I} = \mathbf{V} + G I_y$$

where:
- $\mathbf{I}$ is the vector of all circuit currents (active coil currents + passive structure eddy currents)
- $L$ is the inductance matrix (coil self- and mutual inductances, including plasma)
- $R$ is the resistance matrix (diagonal, coil and structure resistances)
- $\mathbf{V}$ is the applied voltage vector
- $G I_y$ is the back-EMF from plasma current $I_y$ (plasma acts as an additional inductive element)

The "RZIP" acronym reflects the key state variables: R (radial position), Z (vertical position), I (current), P (plasma). In practice the MATLAB code uses `L` as the structure holding all model matrices and `L.G` for the geometric/physical parameters.

### 2.2 The A, M, G Matrices — What They Represent

These three matrices form the core of the real-time hybrid controller signal chain:

**A matrix** (`get_amat`, stored at `\pcs::phys_mat_a` in MDS)  
Size: 24 × 120 (outputs × inputs)

Maps raw diagnostic measurements to 24 "control-oriented observables":
- Inputs (120 channels): flux loops (38), integrated poloidal field probes (38), PF coil currents (18), non-integrated probes (20), special signals (6: δI_OH, FIR, FPS current, interferometer, inner flux loop, TF current)
- Outputs (24 channels):
  - E1–E8 coil currents (or Vloop / OH corrections when in certain modes)
  - F1–F8 coil currents
  - δI_OH (OH current difference)
  - Radial displacement ΔR
  - Plasma current Ip
  - Vertical observer Z·Ip (or IP*ZA)
  - Density
  - Elongation
  - Fast vertical observer signals (differentiated / pre-differentiated)

This is a static gain matrix. It is a linear combination of sensor readings specifically designed to reconstruct the physically relevant control variables from raw measurements.

**G matrix** (`get_gmat`, stored at `\pcs::phys_mat_g1/g2/g3` in MDS)  
Size: 22 × 24, time-varying (up to 10 time windows)

The gain matrix from control-variable errors to coil current-rate requests. It has three components for a PID structure:
- `G2P` (22 × 24): Proportional gain
- `G1I` (22 × 24): Integral gain  
- `G3D` (22 × 24): Derivative gain

The outputs (22) correspond to the M-matrix input channels (E1–E8, F1–F8, OH1–OH2, Ip, ZIp, gas, FPS). The G matrix encodes the control allocation: which coils respond to which errors and with what gain. It can change at discrete time instants during a shot (e.g., different gains during ramp-up vs. flat-top).

**M matrix** (`get_mmat`, stored at `\pcs::phys_mat_m` in MDS)  
Size: 20 × 40 (outputs × inputs)

Maps G-matrix outputs + coil currents (for resistive compensation) to actual coil voltages:
- Inputs: G-matrix outputs (22 "GO" channels) + coil currents for resistive compensation (E1–E8, F1–F8, OH1–OH2)
- Outputs: E1–E8 voltages, F1–F8 voltages, OH1–OH2 voltages, gas valve, FPS voltage

The M matrix also encodes the resistive compensation (`R·Ia` feedforward) and can include an integral term (`R/s`) for improved steady-state decoupling. When coils are in "hybrid mode" (current-driven), the corresponding M matrix columns are zeroed out for those channels.

### 2.3 Signal Chain in the Real-Time Controller

At each 100 µs step (`hybctrl_step.m`):

```
LY (sensor readings)
    │
    ▼ LY2Ain()
Ain (120 × 1)   ← raw diagnostic vector
    │
    ▼ Amat.data * Ain
Aout (24 × 1)   ← control observables (ΔR, Z·Ip, Ip, ...)
    │
    ├── compare with Arefs (24 × 1) from wavegens
    ▼
Aerrs (24 × 1)  ← observable errors
    │
    ▼ ssPIDd (discrete-time PID state-space)
uhyb (24 × 1)   ← PID output
    │
    ▼ G = [G2P, G1I, G3D] (selected at current time)
Min[GO] (22 × 1) ← G-matrix output, enters M matrix "GO" input channels
    │
    ├── Min[Ia] ← actual coil currents (for R compensation)
    ▼
Min (40 × 1)    ← full M matrix input
    │
    ▼ Mmat.data * Min
cout (20 × 1)   ← coil voltage commands
    │
    + ffs (20 × 1)  ← voltage feedforwards from wavegens
    ▼
to power supplies
```

### 2.4 How Matrices Are Obtained from MEQ

The A, G, M matrices are stored in the TCV MDS+ database and retrieved per-shot by `get_amat`, `get_gmat`, `get_mmat`. They are programmed into the PCS (Plant Control System) before a shot by the physics team using the MGAMS (Matrix-based Gain And Mode Selection) interface.

For the **linearized plant model** used in offline control design, the key function is `fgess(L, t, outputs)` which takes the `L` structure from FGE/FGEL (the Grad–Shafranov equilibrium solver) and returns a continuous-time state-space model `(A, B, C, D)` of the plasma + passive structures, with:
- States: passive structure eddy currents + linearized plasma position/shape
- Inputs: active coil voltages `Va`
- Outputs: selectable from `{Ia, Bm, Ff, Fx, Ip, ...}`

The `fge` function returns the full nonlinear equilibrium at each time step (for simulation), while `fgel` / `init_fgel_model` returns the linearized model structure `L` around a given equilibrium.

### 2.5 The RZIP Naming Convention

In the code, "RZIP channels" specifically refers to the four position/current control channels in the G matrix:
- `DELTA_R`: radial displacement (R)
- `IP*ZA`: vertical centroid position times plasma current (Z)
- `CP023A`: fast vertical observer (also called "Zfast", used for derivative channel only)
- `IPA`: plasma current (I, for the "P" = plasma)

The function `inspect_RZIP_channels` and related utilities specifically analyze how the G matrix assigns gains to these four channels.

---

## 3. Shape Control

### 3.1 What Shape Control Points Are

Shape control points are a set of `nc` locations `(r_i, z_i)` on or near the plasma boundary (typically on the separatrix, between the midplane and the X-point) at which the isoflux condition is evaluated as a shape feedback signal.

At each point, the shape observer computes:
- `Fn(i)`: poloidal flux at `(r_i, z_i)`
- `Brn(i)`: radial magnetic field at `(r_i, z_i)`
- `Bzn(i)`: vertical magnetic field at `(r_i, z_i)`

The combination `SC2FB * [Fn; Brn; Bzn]` maps these raw values to "shape feedback" signals, which represent displacements of the boundary along the normal direction at each control point.

### 3.2 How Control Points Are Defined

Control points are computed by `get_shape_control_points(FP, SDP)`:

1. **FBT equilibrium reconstruction** (`fbt('tcv', shot, t, ...)`) is called to get the boundary geometry for a reference shot and time.
2. Control points are placed on the boundary using one of:
   - `SCsrc = 'LY'`: interpolation of the FBT output boundary (most common)
   - `SCsrc = 'rz'`: pure (R,Z) position control mode (no isoflux points)
   - `SCsrc = 'LX'`: based on FBT X-point structure
3. The function `shapobs(SC, Lfbt, LYfbt)` runs the shape observer to compute the sensitivity matrices (`Mysyc`) that linearize the shape feedback around the reference equilibrium.

The resulting `SC` structure contains:
- `SC.rc`, `SC.zc`: arrays of `(nc × nt)` control point coordinates (may vary with time for multi-snapshot)
- `SC.ref`: `(ns × nt)` reference shape feedback values
- `SC.Mysyc`: sensitivity matrix mapping `[Fn; Brn; Bzn]` → shape displacement
- `SC.dims`: cell array of channel names
- `SC.nc`: number of spatial control points
- `SC.ns`: number of shape feedback channels (up to `FP.nShapeFbck = 25`)

### 3.3 The SC (Shape Control) Vector and Error Computation

At each shape controller step (every 2 ms, `shapectrl_step.m`):

```
shapeRefs (ns × 1)   ← reference shape feedback from SC.ref
shapeFbck (ns × 1)   ← measured: STPt.SC2FB * [Fn; Brn; Bzn]
                       (interpolated at control point positions from LY.Fx)

shapeErrs = shapeRefs - shapeFbck
projErrs  = STP.Tpre * shapeErrs    ← project onto nShapeCtrl = 16 controlled directions
filtErrs  ← low-pass filtered projErrs
ctrls ← PID_step(filtErrs, ...)      ← diagonal PID in projected space
shapeCtrl = STP.Tpost * ctrls        ← back-project to nH = 24 hybrid references

deltaH = shapeCtrl                   ← added to hybrid refs  (changes refs for RZIP channels)
deltaI = STP.S2I * shapeCtrl         ← added to coil current references
deltaV = STP.S2V * shapeCtrl         ← added to voltage feedforwards
```

The `Tpre`/`Tpost` matrices are computed via SVD of the shape-to-coil sensitivity matrix, keeping only the most controllable singular directions. `svdtol` is the threshold for selecting which singular values to include.

### 3.4 How the Shape Controller Converts SC Errors to Coil Requests

The shape PID gains (`STP.Kp`, `STP.Ki`, `STP.Kd`) operate in the projected (SVD) space and are computed by `build_shape_tunable_params` using the linearized closed-loop model `ssHyb`. The design typically uses `pidtune` (MATLAB's automatic PID tuner) or a DC-gain-based approach. The resulting corrections are distributed back onto the 24 hybrid reference channels (`shapeCtrl → deltaH`), which modify the Arefs that the inner hybrid loop tracks.

For **MPC mode** (`SDP.tuning = 'mpc'`), the shape controller is replaced by `shapempc_step` which uses a finite-horizon quadratic optimization with a Kalman filter state estimator. The prediction model (`MTP.Akf`, `MTP.Ckf`) is derived from the linearized closed-loop model.

The **QP allocator** (`alloc_step_wrapper`) runs in parallel and redistributes the shape control action across coils to satisfy current limit constraints.

---

## 4. FMag Controller

### 4.1 What FMag Is

FMag (Fast Magnetics controller) is a second-generation, more modular controller that is designed to eventually replace or complement the classical hybrid controller. Key features:

- Operates at the same 100 µs rate as the hybrid controller
- Works with a **unified measurement vector** `X = [Ff; Bm; Ia; Itor; Uf; Bmdot]` (flux loops + field probes + coil currents + toroidal current + loop voltage + field derivatives)
- Uses **linear observers** (`ObsMat * X`) to compute all controlled quantities, with support for **nonlinear quadratic X-point observers**
- Supports up to `nctrl = 10` different controller configurations with smooth switching between them
- Implements bumpless transfer between controllers
- Outputs desired coil current rates `Iadot` and computes voltages as `Va = M·Iadot + R·Ia + feedforwards`
- Has a built-in **anti-windup** scheme based on voltage saturation detection

The entry point is `FMag_control.m` (meqctrl-compatible wrapper) which calls `FMag_step.m` at each time step.

### 4.2 Observer Design

The observer matrix `FTP.ObsMat` is a static matrix that maps the measurement vector `X` to linear observables. It is computed offline by `FMag_generate_mats` / `build_RZIp_Xpoint_controller` and packed into the FTP (Tunable Parameter) structure.

**Linear observers** cover:
- Radial and vertical positions (`DELTA_R`, `IP*ZA`)
- Plasma current `Ip`
- Boundary flux `FB` (max flux among boundary probes)

**Nonlinear X-point observers** are computed per-step in `FMag_step.get_observables`:
1. Linear part: `ObsMat * X` gives flux and derivatives at the X-point: `Fn`, `dFn/dr`, `dFn/dz`, `d²Fn/dr²`, `d²Fn/drdz`, `d²Fn/dz²`
2. Quadratic solve (`2×2 Hessian inversion`) finds the X-point displacement `(δr, δz)` from the Newton step
3. The X-point flux difference from the boundary flux is the "X-point control variable"

**Nonlinear shape observers** use a quadratic expansion of the flux at each control point and solve for the displacement along the control direction.

### 4.3 How FMag Uses the Linear Plasma Model

For **linearization and loop analysis** (`FMag_linearization.m`), the Jacobian of the observer chain with respect to the measurement vector is computed analytically (`FMag_observer_linearization.m`). This gives a linear state-space model of the FMag controller that can be connected with the linearized plant model (`fgess` output) for stability analysis.

The linear model of the FMag controller has the form:
```
A_ctrl = FTP.A_SS(:,:,i_ctrl)
B_ctrl = FTP.B_SS(:,:,i_ctrl) * (dobs/dX)
C_ctrl = M * Tpost * FTP.C_SS(:,:,i_ctrl)
D_ctrl = M * Tpost * FTP.D_SS(:,:,i_ctrl) * (dobs/dX) + diag(R) [feedthrough for R·Ia]
```
where `A_SS, B_SS, C_SS, D_SS` are the PID state-space matrices in the projected observable space, `Tpre` mixes errors, and `Tpost` distributes back to coil directions.

### 4.4 Connection to liuqeml

The FMag controller has an optional NN overwrite path (see Section 5). When `cpar.useliunn = true` and `cpar.rX/zX` are set, the FMag step calls `liuqenn_inference_lih` to replace the LIH-based X-point flux and derivatives with NN predictions:

```matlab
[Fn, drFn, dzFn, dr2Fn, dz2Fn, drzFn] = liuqenn_inference_lih(phys, rX, zX);
signals.nn = struct('Fn',Fn,'drFn',drFn,...);
[outputs,cstate] = FMag_step(FP, FTP, t, signals, cstate);
```

Inside `FMag_step`, when `signals.nn` is non-empty, the six NN outputs overwrite the corresponding `observables` entries for the X-point channels, and the quadratic solve is re-run with the NN inputs.

---

## 5. liuqeml Integration

### 5.1 What liuqeml Is

`liuqeml` is a machine-learning surrogate for the LIUQE equilibrium reconstruction code. It is implemented in a separate repository at `/home/grandin/repos/liuqe-ml/`. The MATLAB interface in the `hybrid` repo calls the pre-trained neural network to compute poloidal flux and field components at arbitrary query points, given the current magnetic measurements.

There are two distinct usage paths in the `hybrid` repo:

1. **Shape feedback (shapectrl):** Used in `shapectrl.m` step when `cpar.useliunn = true`. Replaces the full `qintmex` (flux grid interpolation) call with a NN inference to compute `[Fn, Brn, Bzn]` at the shape control points.

2. **X-point observers (FMag):** Used in `FMag_control.m` step when `cpar.useliunn = true`. Provides flux and its first and second derivatives at X-point positions for the nonlinear X-point observer.

### 5.2 What Inputs/Outputs the ML Model Takes

The NN takes:
- **Physics inputs** (`phys`): `[Bm; Ff; Ia; Uf; rBt]` — a concatenated vector of all available magnetic diagnostics and coil currents (physical measurements, not reconstructed quantities)
  - `Bm`: poloidal field probe signals (Teslas)
  - `Ff`: flux loop signals (Webers)
  - `Ia`: active coil currents (E, F, G, OH)
  - `Uf`: passive structure eddy currents estimated from voltage (Volts)
  - `rBt`: toroidal field proxy (R × Bt product)
- **Query points** (`r`, `z`): spatial coordinates where outputs are requested (meters)

The NN outputs:
- `Fn`: poloidal flux at each query point (Weber/rad)
- `Brn`: radial magnetic field at each query point (Tesla)
- `Bzn`: vertical magnetic field at each query point (Tesla)

### 5.3 The `net_forward2.m` Function

`net_forward2.m` implements the NN forward pass with a dual-encoder architecture:

1. **Points encoder** (`pts_enc`): encodes query positions `(r, z)` through 2 dense layers with Swish activation → `pts_embs` (embedding per query point)
2. **Physics encoder** (`phys_enc`): normalizes and encodes the full physics vector through 2 dense layers → `phys_emb` (single embedding for the whole shot state)
3. **Head network** (`head`): element-wise product `phys_emb × pts_embs` (FiLM-style conditioning), then 2 dense layers → `[Fn, Br, Bz]` per query point

The weights are loaded from `data/3075260/weights.mat` at MEX compile time via `coder.const`. The function is designed for MATLAB Coder code generation (MEX compilation). Typical inference times are ~55–150 µs for 1–25 query points.

### 5.4 How This Differs from Running Full LIUQE

Full LIUQE runs an iterative fixed-point solve of the Grad–Shafranov equation with current profile models (~5–50 ms per reconstruction). The `liuqeml` NN bypasses this entirely:
- No iterative solve
- No current profile model
- Single forward pass through ~4 dense layers per encoder
- Latency: ~55–150 µs (single call), well within the 100 µs controller cycle [GUESS: the MEX version would need to be called; Python inference would be too slow]
- Does not return a full flux map, only values at requested points
- Does not return equilibrium parameters (β, q-profile, etc.) — only local flux and field

The `liuqeml` SCD algorithm object (`SCDalgoobj_liuqeml.m`) is set up as a separate algorithm node running at 100 µs, with its own input bus (Bm, Ff, Ia, Ip, Ft, Iu, rBt, rn, zn) and output bus (Fn, Brn, Bzn). This allows it to be deployed as a standalone RTC node on TCV.

---

## 6. Gain Computation and Tuning

### 6.1 How Gains Are Computed

**Hybrid (RZIP) gains (G matrix):** The G matrix is not computed from scratch in the `hybrid` repo itself — it is retrieved from the TCV MDS database where it was previously programmed by physicists via the MGAMS interface. The functions `get_hybrid_gains` / `inspect_RZIP_channels` only read and interpret the stored gains.

For new gain design, the workflow would be:
1. Build the linearized closed-loop model via `init_hybrid_model` + `build_hybrid_linmod`
2. Analyze stability and performance using standard MATLAB control toolbox
3. Load new gains into MGAMS / the `G_time`/`G_order` tunable parameter arrays

**Shape controller gains (STP.Kp, Ki, Kd):** Computed by `build_shape_tunable_params.m` using one of:
- `pidtune` (MATLAB's automatic PID tuner) — most common
- Manual specification
- Optimization-based tuning (`SDP.tuning = 'optim'`)

The tuning is done in the SVD-projected space (after `Tpre`) and then the closed-loop poles of `ssHyb` (linearized hybrid + plasma) are used to set the crossover frequency (`wcross = 20 Hz` default) and phase margin (`PM = 65°` default).

**FMag gains:** Built by `build_RZIp_Xpoint_controller` using the linearized plant model and pole-placement or specified crossover frequencies.

### 6.2 QP Allocation — Handling Coil Current Limits

The QP allocator (`build_alloc_tunable_params`, `alloc_step_wrapper`) solves a quadratic program at each 2 ms shape controller step to redistribute the shape control action `shapeCtrl` across coils subject to current limit constraints:

```
minimize   ||u - shapeCtrl||²  (with weighting)
subject to:  satneg ≤ I_PF + u ≤ satpos   (current limits)
             zeroneg ≤ u ≤ zeropos         (zero-crossing avoidance)
```

where `u` is the correction to the EF coil currents. The allocator runs alongside the shape controller and adds its output (`allocH`, `allocV`, `allocI`) to the shape controller output.

The annihilator modes (`ATPt.annmode`) support:
- `'static'` (mode 1): precomputed static annihilator matrix
- `'quadprog'` (mode 2): runtime QP solve

### 6.3 The Hybrid Gain Structure

The complete gain parameter structure in tunable params (`get_hybrid_tunable_params`) includes:
- `tp.A_matrix` (24 × 120): the A matrix
- `tp.M_matrix` (20 × 40): the M matrix
- `tp.G1_matrix`, `tp.G2_matrix`, `tp.G3_matrix` (22 × 24 × 10): time-varying G matrix (I, P, D channels, up to 10 time windows)
- `tp.G_time` (50 × 1), `tp.G_order` (50 × 1): scheduling of G matrix transitions
- `tp.egain`, `tp.fgain`, `tp.ohgain`: global gain scalings for E, F, OH coil families
- `tp.shapectrl`: all shape controller parameters (SC points, STP, ATP, MTP) packed into a flat structure
- `tp.dfpskp`, `tp.dfpski`: FPS (fast power supply) proportional and integral gains

---

## 7. Simulink Integration

### 7.1 The HYBRIDlib.slx Simulink Model

[GUESS: The file `HYBRIDlib.slx` is referenced in the codebase but not directly readable. Based on the SCD framework structure:] The Simulink model implements the hybrid controller as an S-function or referenced model block. The bus definitions (`SCDalgo_hybrid_loadbus.m`) define the exact signal widths and types for the Simulink interface. The controller algorithm (`SCDalgo_hybrid`) is configured as an `SCD_algo` object in `SCDalgoobj_hybrid.m`.

### 7.2 SCD Framework Overview

The SCD (SPC Control and Data acquisition) framework is EPFL-specific middleware for real-time Simulink deployment. Key concepts:

- **Algorithm objects** (`SCD_algo`): wrap a Simulink subsystem with timing, parameter loading, and bus definitions
- **Fixed parameters** (`fp`): compiled into the Simulink model at build time, loaded by `SCDalgo_hybrid_loadfp`
- **Tunable parameters** (`tp`): can be updated between shots or even during a shot via MDS wavegens, loaded by `SCDalgo_hybrid_loadtp`
- **Buses**: Simulink bus objects defining input/output signal structures, defined by `SCDalgo_hybrid_loadbus`
- **Wavegens**: time-varying reference signals loaded from MDS database and played back during the shot

### 7.3 Timing

- **Standard hybrid** (`SCDalgoobj_hybrid`): runs at `dt = 1e-4 s` (10 kHz), starts at `t = -4.5 s`
- **FPS hybrid** (`SCDalgoobj_hybridFPS`): runs at `dt = 25e-6 s` (40 kHz) for faster power supply
- **Shape controller**: runs at `FP.shapeTs = 2e-3 s` (500 Hz) internally, triggered every 20 hybrid steps
- **liuqeml** (`SCDalgoobj_liuqeml`): runs at `dt = 1e-4 s` (10 kHz)

---

## 8. Data Flow Summary

### 8.1 Full Pipeline: TCV Diagnostics → Controller → Coil Voltages

```
TCV Hardware
    │
    ├── Flux loops (Ff, 38 channels)
    ├── Poloidal field probes (Bm, 38 channels)
    ├── Coil currents (Ia, 18 channels)
    ├── Special probes (interferometer, TF current, etc.)
    └── Loop voltages (Uf, from passive structures)
    │
    ▼ ADC (analog-to-digital conversion)
    │
    ▼ AIN BUS (SCDalgo_hybrid_Ain)
    │   (ctlint_vloop, ctlint_rvloop, ctlint_bpol_avg, etc.)
    │   preprocessed (integrated, averaged) sensor signals
    │
    ├──── [liuqeml path] ──────────────────────────────────────────┐
    │     liuqenn_inference(phys, rn, zn) → [Fn, Brn, Bzn]         │
    │     replaces qintmex in shapectrl and/or FMag                 │
    │                                                               │
    ▼ [standard path]                                               │
    A matrix (24×120)                                               │
    │                                                               │
    ▼                                                               │
    Aout (24 control observables)  ◄────────────────────────────────┘
    │
    + Arefs (from MDS wavegens: hyb_A_matrix_refs 1:24)
    │
    ▼ Error: Aerrs = Arefs - Aout
    │
    ▼ ssPIDd (discrete-time PID in 24-dim space)
    │
    ▼ G matrix [G2P, G1I, G3D] × PID output
    │
    + coil currents (for resistive compensation via M matrix)
    │
    ▼ M matrix (20×40) → 20 coil voltage requests
    │
    + voltage feedforwards (from MDS wavegens: hyb_FF_e/f/oh_u)
    │
    ▼ Power supplies
    │
    ▼ Coil voltages Va → plasma
```

### 8.2 Offline vs. Online Use of FBT/FGE

**Offline (before a shot):**
1. `fbt('tcv', shot, t)` reconstructs the reference equilibrium from the reference shot
2. `get_shape_control_points` computes the isoflux control points and sensitivity matrices
3. `init_hybrid_model` builds the linearized closed-loop state-space model
4. `build_shape_tunable_params` / `build_hybrid_linmod` compute the gain matrices
5. Results are saved and loaded into the TCV PCS via MDS wavegens and tunable parameters

**Online (during a shot):**
- The A, G, M matrices and references are static matrices/wavegens loaded at shot start
- The controller computes `Aout = A × Ain`, errors, PID, G×PID, M×Min at each 100 µs step
- The shape controller computes flux values at control points (via `qintmex` or `liuqeml`) every 2 ms

---

## 9. Issues for Surrogate Model Integration

### 9.1 What the Surrogate Model Needs to Replace

There are two distinct roles a surrogate model could play:

**Role A: Online shape feedback (replacing LIUQE/qintmex)**  
This is what `liuqeml` currently does in `shapectrl.m` and `FMag_control.m`. The surrogate replaces the call to `qintmex(L.rx, L.zx, LY.Fx, ri, zi, qpM)` which interpolates the FGE flux grid at the control point locations. The surrogate must produce `[Fn, Brn, Bzn]` at the shape control points given current measurements.

**Role B: Online equilibrium reconstruction (replacing full LIUQE)**  
A deeper integration that would replace LIUQE entirely for providing the equilibrium state `LY` (flux map `LY.Fx`, plasma position `LY.rIp`, `LY.zIp`, etc.). This is a significantly harder problem as the full flux map is a high-dimensional output.

**Role C: Offline equilibrium for control design (replacing FBT/FGE)**  
The surrogate could also be used offline in `init_hybrid_model` to generate the linearized model `L` and the shape control points `SC` without running FGE/LIUQE.

### 9.2 What Inputs the Controller Expects from the Equilibrium Code

The `LY` structure passed to the controller contains the following fields used by the controller:

| Field | Description | Used by |
|-------|-------------|---------|
| `LY.t` | time (scalar) | hybctrl, shapectrl |
| `LY.Ia` | active coil currents (na × 1) | hybctrl_step (resistive comp.), shapectrl |
| `LY.Va` | active coil voltages (na × 1) | FMag_control |
| `LY.Ff` | flux loop signals (nFf × 1) | shapectrl (useliunn path), FMag |
| `LY.Bm` | poloidal field probe signals (nBm × 1) | shapectrl (useliunn path), FMag |
| `LY.Uf` | passive structure voltages (nIu × 1) | shapectrl (useliunn path), FMag |
| `LY.rBt` | R×Bt product | shapectrl (useliunn path), FMag |
| `LY.Fx` | 2D flux grid (nzx × nrx) | shapectrl standard path (qintmex) |
| `LY.rIp`, `LY.zIp` | current centroid position | shapectrl (position control mode) |
| `LY.Ip` | plasma current | shapectrl (for normalization) |

The `L` structure (from `fgel`/`init_fgel_model`) is used for:
- `L.rx`, `L.zx`, `L.drx`, `L.dzx`: grid definition for flux interpolation
- `L.G.dima`: coil names and ordering
- `L.G.Ra`: coil resistances (for resistive compensation)
- `L.G.Iamax`, `L.G.Iamin`: coil current limits
- `L.G.na`: number of active coils
- `L.P.shot`: shot number

### 9.3 Speed Requirements

The real-time constraint is strict:
- **Hybrid controller step:** must complete in < 100 µs (10 kHz cycle)
- **Shape controller step:** must complete in < 2 ms (500 Hz cycle)
- **liuqeml inference:** measured at 55–150 µs for 1–25 query points (MEX compiled)
- **Full LIUQE:** ~5–50 ms — far too slow for real-time; currently runs offline or at much lower rate

The `net_forward2.m` function with the MEX compilation target (`%#codegen`) is designed to meet the 100 µs constraint for the online inference path. For the FMag X-point observer path, only 1–2 query points are typically needed, well within budget.

For **Role A (replacing qintmex)**: the surrogate must produce `[Fn, Brn, Bzn]` at up to `nShapeFbck = 25` points within 2 ms. The current `liuqeml` NN achieves this.

For **Role B (replacing LIUQE)**: the surrogate would need to produce the full flux map `LY.Fx` or at least all the quantities listed above, within the 100 µs cycle. This is significantly more demanding and would require a very compact model or a dedicated architecture.

### 9.4 Key Interface Functions for Surrogate Integration

The cleanest integration points are:

1. **In `shapectrl.m` (line ~180):** The `if cpar.useliunn` branch that calls `liuqenn_inference`. A surrogate model can be plugged in here by replacing `liuqenn_inference` with an equivalent function that takes the same `phys` vector and control point coordinates `ri`, `zi`.

2. **In `FMag_control.m` / `FMag_step.m`:** The `signals.nn` path for X-point observers. A surrogate providing flux and derivatives at X-point positions can be integrated here.

3. **In `LY2Ain.m`:** [GUESS] This function maps the `LY` structure to the A-matrix input vector. If the surrogate replaces LIUQE, it would need to produce outputs compatible with this function.

4. **In `init_fgel_model.m` / `fge`:** For offline design, a surrogate could provide the equilibrium `L` structure to `init_hybrid_model`.

The most practical near-term integration is Role A (already implemented via `liuqeml`) with possible extension to include first and second derivatives of the flux for improved shape observer accuracy.
