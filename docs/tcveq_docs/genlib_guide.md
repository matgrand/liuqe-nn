# genlib Usage Guide

**Purpose**: When to use each function, example call patterns, common workflows, and a dependency diagram.

---

## When to Use Each Function

### Green's functions (`greenem`, `green3d`, `ellipfast`)

Use `greenem` whenever you need:
- The mutual inductance matrix between sets of toroidal current filaments (coils, vessel, plasma)
- Magnetic field (Br, Bz) produced by filaments
- Derivatives of inductance or field (needed for Jacobians in equilibrium solvers)

Use `green3d` when:
- One or both conductor sets are non-axisymmetric (3D saddle coils, correction coils)
- Loops are given as piecewise-linear polygons in (r, z, phi) coordinates

Use `ellipfast` when:
- You need K(m) or E(m) directly with < 2e-8 accuracy
- You want the fastest possible scalar or vectorized evaluation

### Vessel modes (`vveig`)

Use `vveig` when:
- Reducing a high-dimensional vessel filament model to a smaller set of RL modes
- Building the reduced-order passive conductor model for real-time control
- Computing effective time constants of vessel current diffusion

### B-splines (`bspbase`, `bspsum`, `bsp2deval`, `csdec`)

Use when:
- Representing the plasma flux function ψ(r,z) as a 2D B-spline
- Computing ψ and its spatial derivatives at arbitrary (r,z) points
- Fitting equilibrium fields on a grid with smooth interpolation

### Cubic splines (`csmod`, `csder`)

Use when:
- Simple 1D interpolation with known edge conditions
- Working with piecewise polynomial data in MATLAB pp-form
- Derivatives are specified at knots (Hermite interpolation → `csder`)

### Signal processing (`sti`, `resamp`, `rpdcov`, `toepconv`)

Use when:
- Fitting a transfer function model to frequency response data (`sti`)
- Resampling diagnostic signals to a common time base (`resamp`)
- Converting rational fraction coefficient uncertainties to residue-pole uncertainties (`rpdcov`)

### Utilities

| Utility | When to use |
|---|---|
| `ifloor`, `iceil`, `iround` | Binary search in sorted arrays (faster than `find`) |
| `xmtimes` | Matrix multiply when inputs may be N-D arrays |
| `isort`, `imax`, `imin` | When you only need the index, not the value |
| `ifc(c,a,b)` | Inline ternary in expressions |
| `ppzero` | Finding roots of a fitted spline |
| `appendstruct` | Merging parameter structures |
| `simple_polygon` | Validating plasma boundary contours |

---

## Example Call Patterns

### Building a mutual inductance matrix (Mvv)

```matlab
% Vessel filaments: r_v [Nv x 1], z_v [Nv x 1]
Mvv = greenem('mut', r_v, z_v);        % [Nv x Nv], self on diagonal = 0

% For a rectangular conductor of dimensions a x b:
% greenem('mut', r_v, z_v, [a_v(:) b_v(:)])  % includes self-inductance on diagonal

% Coil-to-vessel mutual inductances:
Mcv = greenem('mut', r_c, z_c, r_v, z_v);  % [Nc x Nv]

% Magnetic field at sensor locations from vessel filaments:
Br  = greenem('br', r_s, z_s, r_v, z_v);   % [Ns x Nv]
Bz  = greenem('bz', r_s, z_s, r_v, z_v);   % [Ns x Nv]

% With post-multiplier (e.g., turn matrix):
Mce = greenem('mut', r_c, z_c, r_v, z_v, T_ve);  % [Nc x Ne]
```

### Multiple outputs in one call

```matlab
% Compute mutual inductance AND its derivative w.r.t. r2 simultaneously:
[M, dMdr2] = greenem({'mut','dr2m'}, r1, z1, r2, z2);
% M: [N1 x N2], dMdr2: [N1 x N2]

% All field components at once:
[Br, Bz] = greenem({'br','bz'}, r_sensors, z_sensors, r_coils, z_coils);
```

### Self-inductance only

```matlab
% Circular cross-section, radius a:
L_self = greenem('self', r_coil, [], a_coil);         % [Nc x 1]

% Rectangular cross-section, dimensions a x b:
L_self = greenem('self', r_coil, [], a_coil, b_coil); % [Nc x 1]
% or equivalently:
L_self = greenem('self', r_coil, [], [a_coil, b_coil]);
```

### Vessel eigenmode computation

```matlab
% Compute eigenmodes with resistance normalization (default):
[Tve, Le, Re, sym] = vveig(Mvv, Rv);  % Rv = [Nv x 1] resistance vector

% Transform mutual matrix to eigenmode basis:
Mce = Mcv * Tve;           % [Nc x Nv] * [Nv x Nu] = [Nc x Nu]
Mee = Tve' * Mvv * Tve;    % = diag(Le)

% Convert eigenmode currents back to filament currents:
I_vessel = Tve * I_eigen;  % [Nv x 1]

% Truncate to Nu dominant modes:
Tve_reduced = Tve(:, 1:Nu);   % take first Nu columns (longest time constants)
```

### B-spline flux function representation

```matlab
% Setup: knot sequences for r and z
Nr = 20; Nz = 30;
r_min = 0.6; r_max = 2.2;
z_min = -1.0; z_max = 1.0;
t_r = linspace(r_min, r_max, Nr);
t_z = linspace(z_min, z_max, Nz);

% Build the decomposition matrices (fit operators):
M_r = csdec(t_r, r_grid_points);    % [nb_r x Nr]
M_z = csdec(t_z, z_grid_points);    % [nb_z x Nz]

% Given flux values on a grid psi_grid [Nr x Nz]:
C = M_r * psi_grid * M_z';   % B-spline coefficients [nb_r x nb_z]

% Evaluate flux and derivatives at arbitrary points:
[psi, dpsi_dr, dpsi_dz, d2psi_drr, d2psi_drz, d2psi_dzz] = ...
    bsp2deval(taug_r, taug_z, C, r_pts, z_pts, 2, npar);
% Note: taug_r and taug_z are the augmented knot sequences returned by csdec
```

### Using bspbase and bspsum directly

```matlab
% Build basis matrix:
k = 4;  % cubic splines
t = linspace(0, 1, 20);  % interior knots (simple knot sequence)
% augment:
dt = t(2)-t(1);
taug = [(-3:-1)*dt, t, 1+(1:3)*dt];
x = linspace(0,1,100);
B = bspbase(taug, k, x);   % [nb x 100] where nb = length(taug)-k

% Given coefficients C [nb x 1], evaluate spline:
y = bspsum(taug, C, x);    % equivalent to C' * B but faster

% Evaluate derivative:
dy = bspsum(taug, C, x, 1);

% 2D with multiple outputs:
[f, gr] = bspsum(taug, C_matrix, x, 0);  % C_matrix: [nb x M]
```

### 1D cubic spline interpolation

```matlab
% Standard not-a-knot spline:
y_interp = csmod(x_knots, y_knots, 'a', x_eval);

% Natural spline (zero 2nd derivative at endpoints):
y_interp = csmod(x_knots, y_knots, 'n', x_eval);

% Get the pp-form for later evaluation:
pp = csmod(x_knots, y_knots);  % omit last argument
y_eval = ppval(pp, x_eval);

% Hermite spline (both values and derivatives known):
y_interp = csder(x_knots, y_knots, yp_knots, x_eval);
```

### Transfer function fitting (`sti`)

```matlab
% Experimental frequency response H at frequencies w [rad/s]:
% H: [Nfreq x Noutput] complex matrix
[A, B, RAB, Hfit, iter] = sti(w, H, 2, 2);
% A: denominator coefficients [1, a1, a2] (degree 2)
% B: numerator coefficients [b1, b2] (degree 1) per output

% Transfer function: H(s) = B(s)/A(s) = (b1*s + b2)/(s^2 + a1*s + a2)

% Convert to residue-pole form with uncertainty:
[r, p, d, Rrpd] = rpdcov(A, B, RAB);
```

### Signal resampling

```matlab
% Resample two signals to fastest common time base, linear interpolation:
[X1_r, X2_r, T_common] = resamp(X1, T1, X2, T2, 'fel');
% 'f' = fastest period, 'e' = extended time range, 'l' = linear interp

% Resample to a specific time vector:
[X1_r, T_common] = resamp(X1, T1, T_target, 'iel');
% 'i' = imposed time base, 'e' = extended, 'l' = linear
```

### Array index utilities

```matlab
% Binary search in sorted array:
t = 0:0.01:10;
x = [2.345, 5.678, 9.999];
k_floor = ifloor(t, x);   % index of largest t(k) <= x(i)
k_ceil  = iceil(t, x);    % index of smallest t(k) >= x(i)
k_near  = iround(t, x);   % index of nearest t(k) to x(i)

% Unsorted x: pass 3rd argument
k_floor = ifloor(t, x_unsorted, 1);

% Index of extremes:
k_max = imax(A);    % equivalent to [~, k_max] = max(A)
k_min = imin(A);
k_sorted = isort(A);  % equivalent to [~, k_sorted] = sort(A)
```

### 3D mutual inductance

```matlab
% One axisymmetric loop (phi1=[]) vs a 3D polygon:
% r2, z2, phi2 are [Npts x Nloops] matrices (NaN-padded for variable length)
M = green3d(r_circle, z_circle, [], r2_polygon, z2_polygon, phi2_polygon);
% Returns [numel(r_circle) x Nloops]

% Two 3D loops:
M = green3d(r1, z1, phi1, r2, z2, phi2);
% Uses numerical Neumann formula, ~O(N1*N2) cost
```

---

## Common Workflows in MEQ

### Workflow 1: Build complete inductance matrix set

Given coil positions (r_c, z_c), vessel filament positions (r_v, z_v), and plasma boundary filaments (r_p, z_p):

```matlab
% 1. Vessel-vessel mutual inductance matrix
Mvv = greenem('mut', r_v, z_v);                    % [Nv x Nv]

% 2. Coil-vessel mutual inductance
Mcv = greenem('mut', r_c, z_c, r_v, z_v);         % [Nc x Nv]

% 3. Coil-coil mutual inductance (with self on diagonal)
% For circular cross-section conductors:
Mcc = greenem('mut', r_c, z_c, [a_c, b_c]);       % [Nc x Nc], self on diagonal

% 4. Plasma-coil mutual inductance
Mpc = greenem('mut', r_p, z_p, r_c, z_c);         % [Np x Nc]

% 5. Vessel eigenmode basis
[Tve, Le, Re] = vveig(Mvv, Rv);

% 6. Transform to eigenmode basis:
Mce = Mcv * Tve;    % coil-to-eigenmode
Mpe = Mpc * Tve;    % plasma-to-eigenmode (via vessel)

% 7. Sensor response matrices (flux loops and Rogowski coils):
Bz_matrix = greenem('bz', r_s, z_s, r_c, z_c);   % [Ns x Nc]
Br_matrix = greenem('br', r_s, z_s, r_c, z_c);
```

### Workflow 2: Represent and evaluate plasma flux

```matlab
% Given: 2D grid of flux values psi_grid [Nr x Nz]
% on r_pts [Nr x 1], z_pts [Nz x 1]

% Build fitting operators
[M_r, taug_r] = csdec(linspace(r_min,r_max,Nr), r_pts);
[M_z, taug_z] = csdec(linspace(z_min,z_max,Nz), z_pts);

% Compute B-spline coefficients
C = M_r * psi_grid * M_z';   % [nb_r x nb_z]

% Evaluate at arbitrary points (e.g., for MHD equilibrium constraints):
r_eval = ...; z_eval = ...;
[psi, psi_r, psi_z, psi_rr, psi_rz, psi_zz] = ...
    bsp2deval(taug_r, taug_z, C, r_eval, z_eval, 2, npar);
```

### Workflow 3: Vessel time constant reduction

```matlab
% Start with fine vessel filament model
Mvv_fine  = greenem('mut', r_vessel_all, z_vessel_all); % [Nv x Nv]
Rv_fine   = resistance_vector;  % [Nv x 1]

% Compute eigenmodes (all Nv)
[Tve, Le, Re] = vveig(Mvv_fine, Rv_fine, 'R');

% Time constants in seconds:
tau = Le / Re;   % eigenmode time constants L/R [s]
fprintf('Longest time constant: %.3f s\n', max(tau));
fprintf('Shortest kept time constant: %.3f s\n', tau(Nu));

% Truncate to Nu modes (those with tau > tau_min):
Nu = sum(tau > tau_min);
Tve_u = Tve(:, 1:Nu);

% Reduced-order mutual inductance and resistance matrices:
Muu = Tve_u' * Mvv_fine * Tve_u;  % = diag(Le(1:Nu))
Ruu = Re * eye(Nu);                 % = scalar*I for type='R'

% Test: coil-to-mode mutual
Mcu = Mcv_fine * Tve_u;  % [Nc x Nu]
```

---

## Dependency Diagram

```
greenem.m
  └── greenemmex (MEX)
        └── greenem.c           [core computation]
              └── ellipK/E     [inline polynomial, Abramowitz 17.3.34]

green3d.m
  ├── greenem.m                 [axisymmetric case]
  └── greenxyz (inline)         [3D Neumann formula]

ellipfast.m
  └── ellipfastmex (MEX)
        └── ellipfast.c         [same polynomial as greenem.c]

vveig.m
  ├── isort.m                   [sort index utility]
  └── [no MEX, pure MATLAB]

bspbase.m
  └── bspbasemex (MEX)
        └── bspbase.c           [Cox-de Boor recurrence]

bspsum.m (bspsum is a MEX file itself — no .m wrapper beyond documentation)
  └── bspsummex (MEX, OpenMP)

bsp2deval.m
  └── bspsum.m (called twice per derivative order)

csdec.m
  ├── bspbase.m                 [build basis]
  └── bspsum.m (taug, M, xx)   [optional evaluation at xx]

cs2deval (MEX)
  └── cs2deval.c
        ├── csint               [cubic B-spline local basis, inline]
        └── adotc               [4x4 coefficient multiply, inline]

csmod.m                         [pure MATLAB, no MEX]
csder.m                         [pure MATLAB, no MEX]

sti.m                           [pure MATLAB, no MEX]

rpdcov.m
  └── rpdjacobian.m
        └── toepconv.m

resamp.m
  └── [interp1, pure MATLAB]

simple_polygon.m
  └── [pure MATLAB, sortrows]

ifloor.m, iceil.m, iround.m
  └── iflooriceiliroundmex (MEX, single shared binary)

isnans.m
  └── nansmex (MEX)

bspsum is sometimes listed as a MEX file directly (no .m wrapper) — Contents.m
entry "bspsum" is purely a comment block at the top of bspsum.m which acts as a
MEX stub (the actual computation is in bspsummex).
```

---

## Notes on Data Conventions

### Coordinate system

All `greenem` functions use cylindrical coordinates (r, z):
- r: major radius (meters), always positive
- z: height above midplane (meters), can be negative

The toroidal angle phi is integrated out (axisymmetric assumption). The result is per unit of the toroidal angle integral already included in the Neumann formula.

### Matrix indexing convention

`greenem('mut', r1, z1, r2, z2)` returns a [N1 x N2] matrix where:
- Row index i corresponds to filament i in (r1, z1)
- Column index j corresponds to filament j in (r2, z2)
- Element [i,j] = M(r1(i),z1(i), r2(j),z2(j))

When r2/z2 are omitted, the result is [N1 x N1] with zeros on the diagonal (mutual only, no self).

### Units

- Mutual inductance: Henry [H]
- Self-inductance formulas: Henry [H], requires conductor dimensions in meters
- Magnetic field: Tesla [T] = Weber/m^2 = H/m * A/m

`mu0 = 4e-7*pi = 1.256637061435917e-6` H/m (hardcoded in greenem.c)

### B-spline order convention

genlib uses **order k = degree + 1**:
- k=1: piecewise constant
- k=2: piecewise linear
- k=4: cubic (most common in MEQ)

The number of basis functions is `nb = length(t) - k`.

### Knot sequence for csdec

`csdec` automatically augments the user-provided knot sequence `t` by extending it 3 intervals on each side with uniform spacing. This creates `length(t) + 6` augmented knots, giving `length(t) + 2` cubic basis functions. The augmentation uses the spacing `dt = (t(end)-t(1))/(length(t)-1)`, so the knot sequence need not be uniform but the extension assumes uniform spacing at the ends.

### Resistance vector in vveig

`Rv` should be a column vector [Nv x 1]. If you have a resistance matrix `Rvv` (e.g., from MDSplus), use `diag(Rvv)` to extract the diagonal.

### Transfer function convention in sti

`sti` uses MATLAB `polyval` convention (descending coefficient order). The denominator polynomial A has leading coefficient 1 (monic). The transfer function is:
```
H(s) = polyval(B, s) / polyval(A, s)
```
with s = j*w in rad/s. After fitting, `sti` denormalizes by reversing the `w/wnorm` scaling.
