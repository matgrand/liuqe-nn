# genlib Technical Reference

**Library**: genlib — General Purpose Library for SPC/TCV
**Version**: 3.14.0
**Copyright**: 2022–2024 Swiss Plasma Center EPFL (Apache 2.0 License)
**Maintainers**: Antoine Merle, Federico Felici
**Role**: Mathematical foundation for the MEQ plasma equilibrium toolbox used at the TCV tokamak (EPFL Swiss Plasma Center)

---

## Section 1 — Overview

genlib is a MATLAB/C MEX library providing:

1. **Green's functions in toroidal geometry** — mutual inductance between current loops, magnetic field components, and all their partial derivatives. This is the core physical computation for building the inductance matrices (Mcc, Mcs, Mvv, Mpc, etc.) used throughout MEQ.
2. **Fast elliptic integrals** — K(m) and E(m) via degree-4 polynomial approximations, accurate to < 2e-8.
3. **Vessel eigenmode decomposition** — compression of vessel filament models into dominant RL-circuit modes.
4. **B-spline and cubic spline machinery** — basis evaluation, coefficient fitting, 2D tensor-product evaluation; used for representing flux functions ψ(R,Z).
5. **Signal processing utilities** — rational transfer function fitting (STI), resampling, Toeplitz convolution, residue-pole-covariance.
6. **Array and structure utilities** — typed index search, polygon geometry, etc.

The library consists of MATLAB wrappers calling compiled C MEX files (`*mex`). OpenMP parallelism is used in `bspsummex`. The mex files are built with `make CC=gcc MATPATH=...`.

---

## Section 2 — Green's Functions

### 2.1 Physical background

In toroidal (axisymmetric) geometry the Neumann formula for mutual inductance between two circular current loops reduces to:

$$M(r_1,z_1,r_2,z_2) = \mu_0 \sqrt{r_1 r_2} f(k)$$

where the modulus $k$ is defined by:

$$k^2 = \frac{4r_1r_2}{(r_1+r_2)^2 + (z_1-z_2)^2}$$

and $f(k)$ involves the complete elliptic integrals $K(k^2)$ and $E(k^2)$:

$$f(k) = \frac{(2 - k^2)K(k^2) - 2E(k^2)}{k}$$

In the C implementation the variable `m` is $k^2$ (the parameter, not the modulus), and the variable $x = d \cdot is = 1 - m$ is the complementary parameter (i.e., $x = 1 - k^2$). Using $x$ instead of $1-m$ is numerically more precise when the two loops nearly coincide ($k^2$ close to 1, $x$ close to 0).

### 2.2 `greenem.m` — Primary interface

**Signature**:
```matlab
[Y1,...] = greenem(mode, R1, Z1 [, R2, Z2] [, T])
```

**Arguments**:
- `mode`: string or cell array of strings selecting what to compute
- `R1, Z1`: coordinates of the first set of filaments (column vectors, N1 elements)
- `R2, Z2`: coordinates of the second set of filaments (column vectors, N2 elements); if omitted, R2=R1, Z2=Z1 (self-matrix)
- `T`: optional post-multiplier matrix; output becomes Y*T

**Output**: each output is [N1 x N2] (or [N1 x N2]*T if T given)

**Mode strings and their integer codes (used internally)**:

| Mode string(s) | Code | Meaning |
|---|---|---|
| `'mut'`, `'m'` | 0 | Mutual inductance M(r1,z1,r2,z2) |
| `'dr1m'`, `'dmdr1'` | 1 | ∂M/∂r1 |
| `'dz1m'`, `'dmdz1'` | 2 | ∂M/∂z1 |
| `'dr2m'`, `'dmdr2'` | 3 | ∂M/∂r2 |
| `'dz2m'`, `'dmdz2'` | 4 | ∂M/∂z2 |
| `'br'` | 5 | Br field at (r1,z1) from unit current at (r2,z2) |
| `'bz'` | 6 | Bz field at (r1,z1) from unit current at (r2,z2) |
| `'dr1r1m'` | 11 | ∂²M/∂r1² |
| `'dr1z1m'` | 12 | ∂²M/∂r1∂z1 |
| `'dr1r2m'` | 13 | ∂²M/∂r1∂r2 |
| `'dr1z2m'` | 14 | ∂²M/∂r1∂z2 |
| `'dz1z1m'` | 22 | ∂²M/∂z1² |
| `'dz1r2m'` | 23 | ∂²M/∂z1∂r2 |
| `'dz1z2m'` | 24 | ∂²M/∂z1∂z2 |
| `'dr2r2m'` | 33 | ∂²M/∂r2² |
| `'dr2z2m'` | 34 | ∂²M/∂r2∂z2 |
| `'dz2z2m'` | 22 | ∂²M/∂z2² (same as dz1z1m by symmetry) |
| `'dr1br'` | 51 | ∂Br/∂r1 |
| `'dz1br'` | 52 | ∂Br/∂z1 |
| `'dr2br'` | 53 | ∂Br/∂r2 |
| `'dz2br'` | 54 | ∂Br/∂z2 |
| `'dr1bz'` | 61 | ∂Bz/∂r1 |
| `'dz1bz'` | 62 | ∂Bz/∂z1 |
| `'dr2bz'` | 63 | ∂Bz/∂r2 |
| `'dz2bz'` | 64 | ∂Bz/∂z2 |
| `'self'`, `'s'` | -1 | Self-inductance (requires conductor size) |
| `'zero'`, `'0'` | 99 | Return zeros |

**Self-inductance formulas** (returned on diagonal or as standalone):

For a circular conductor of radius `a`:
```matlab
s = 4e-7*pi * r1 .* (log(8*r1./a) - 1.75)   % [H]
```

For a rectangular conductor `a x b`:
```matlab
s = 4e-7*pi * r1 .* (log(8*r1./(a+b)) - 0.50)  % [H]
```

**Parallelism**: the MEX function uses `npar = max(1, min(maxNumCompThreads, ceil(N1*N2/1e4)))` threads.

### 2.3 `greenem.c` — Core computation

The C function signature is:
```c
void greenem(int *mode, int nm, genlib_bool der1, genlib_bool der2,
             double r1, double z1, double r2, double z2, double **y, int o)
```

**Variable definitions** (immediately after entry):

$$r = r_1 + r_2 \quad \text{(sum of radii)}$$
$$h = z_1 - z_2 \quad \text{(axial separation)}$$
$$w = r_1 - r_2 \quad \text{(radial difference)}$$
$$s = r^2 + h^2 = (r_1+r_2)^2 + (z_1-z_2)^2$$
$$d = w^2 + h^2 = (r_1-r_2)^2 + (z_1-z_2)^2$$
$$is = 1.0/s$$

Then $r_1, r_2, r$ are doubled:

$$m = r_1 r_2 is = \frac{4r_1 r_2}{s} = k^2$$
$$x = d \cdot is = 1 - k^2 \quad \text{(complementary parameter)}$$
$$lx = \log(x)$$

**Elliptic integral computation** (inline polynomial in greenem.c, parameterized by $x = 1 - k^2$):

$K(m)$ where $m = 1-x$:
$$K = (((0.01451196212 x + 0.03742563713) x + 0.03590092383) x + 0.09666344259) x + 1.38629436112$$
$$- ((((0.00441787012 x + 0.03328355346) x + 0.06880248576) x + 0.12498593597) x + 0.50000000000) \log(x)$$

$E(m)$ where $m = 1-x$:
$$E = (((0.01736506451 x + 0.04757383546) x + 0.06260601220) x + 0.44325141463) x + 1.00000000000$$
$$- ((((0.00526449639 x + 0.04069697526) x + 0.09200180037) x + 0.24998368310) x) \log(x)$$

**Core Green's function building blocks**:

$$v_0 = (1.0 - 0.5m)K - E \quad \text{(the key combination)}$$
$$u_0 = \sqrt{s} = \sqrt{(r_1+r_2)^2 + (z_1-z_2)^2}$$

The mutual inductance is:

$$M = \mu_0 u_0 v_0 = \mu_0 \sqrt{(r_1+r_2)^2 + (z_1-z_2)^2} \left[\left(1 - \frac{2r_1r_2}{s}\right)K - E\right]$$

where $\mu_0 = 4\pi \times 10^{-7} = 1.256637061435917 \times 10^{-6}$ H/m, $\mu_1 = 2 \times 10^{-7}$ H/m.

Note: in the C code, r1, r2 have been doubled, so `s` after the doubling is `(2r1+2r2)^2 + h^2 = 4*(r1+r2)^2 + h^2`. The formula reconstitutes via u0=sqrt(s) and v0=(1-0.5*m)*K-E with m=4*r1orig*r2orig/sorig. The actual Neumann formula result is mu0 * sqrt(sorig) * vorig where sorig uses the original un-doubled radii. The C code doubles them so that `m = r1_doubled * r2_doubled * is = 4*r1orig*r2orig / (4*(r1orig+r2orig)^2+h^2)` — this matches k^2 exactly.

**Macro definitions for output**:
```c
#define mu0 1.256637061435917e-6   /* 4e-7*pi */
#define mu1 2e-7
#define G0(u0,v0)             (u0*v0)
#define G1(u0,u1,v0,v1,dx,kx) (u1*v0*dx + u0*v1*kx)
#define G2(u0,u1,u2,v0,v1,v2,dx,dy,kx,ky,kxy) \
    (2.0*u1*v0 + u0*v1*kxy + u2*v0*dx*dy + u1*v1*(dx*ky+kx*dy) + u0*v2*kx*ky)
#define GZ(...)  /* same as G2 but without the 2*u1*v0 leading term */
```

**First-derivative quantities** (computed when der1 = true):
```c
s2  = 2*h;
mr  = m*r;   // m*(2*(r1+r2)), r here is already doubled sum
m1  = (2*r2 - mr) * is;   // dm/dr1
m3  = (2*r1 - mr) * is;   // dm/dr2
m2  = -m * s2 * is;       // dm/dz1 (= dm/dh * dh/dz1 = -dm/dz2)
im  = 1.0/x;              // 1/(1-k^2)
v1  = 0.25*(E*im - K);    // dv0/dm
u1  = 0.5/u0;             // du0/ds * (ds/d(r1 or r2 or z1))... scaled
ir1 = 2.0/r1;             // 1/r1_original (r1 was doubled)
```

**Output for each mode** (key cases):
```c
case 0: M   = mu0 * G0(u0, v0)                           // mutual inductance
case 1: dMdr1 = mu0 * G1(u0,u1,v0,v1, s1=r, m1)         // ∂M/∂r1
case 2: dMdz1 = mu0 * G1(u0,u1,v0,v1, s2=2h, m2)        // ∂M/∂z1
case 3: dMdr2 = mu0 * G1(u0,u1,v0,v1, s3=r, m3)         // ∂M/∂r2
case 4: dMdz2 = -mu0 * G1(u0,u1,v0,v1, s2, m2)          // ∂M/∂z2 = -∂M/∂z1
case 5: Br    = -mu1 * G1(...) * ir1                      // Br = -(1/2pi*r1)*∂M/∂z1
case 6: Bz    =  mu1 * G1(...) * ir1                      // Bz = (1/2pi*r1)*∂M/∂r1
```

### 2.4 `dgreenem.c` — Legacy derivative implementation

This is an older C implementation (circa 1995, updated 2008) with a different interface and parameter convention. It takes individual `mode` integers:

```
mode: mut=0, br=1, bz=2, dbzdr1=3, dbzdz1=4, dbrdr1=5, dbrdz1=6,
      dbzdr2=7, dbzdz2=8, dbrdr2=9, dbrdz2=10, dmdr1=11, dmdz1=12,
      dmdr2=13, dmdz2=14
```

**Key difference in convention**: dgreenem uses k^2 directly as the modulus argument to `ellipfast(k2, kind)`, and computes:
```c
k2 = 4*r1*r2 / u2;   // u2 = (r1+r2)^2 + h^2
```

Mode 0 (mutual inductance):
```c
return ((2.-k2)*K-2.*E) * 6.28318530718e-07 * u;
// = mu0 * u * ((1 - k2/2)*K - E)  [note: mu0/2 = 6.28e-7 ??? ]
// Actually 6.28e-7 * 2 * ... wait:
// 6.28318530718e-07 = pi * 2e-7 = mu0/2
// so result = mu0/2 * 2 * u * ((1-k2/2)*K - E)
// but greenem.c has mu0 * u0 * ((1-m/2)*K - E) with u0 = sqrt(s)
// The factor of 2 difference is absorbed in the choice of u vs sqrt(s)
```

Note: dgreenem.c appears to be legacy code that is no longer called by the main greenem.m path (which calls `greenemmex`). The newer greenem.c is the active implementation.

### 2.5 `ellipfast.m` and `ellipfast.c`

**Source**: Abramowitz & Stegun, formula 17.3.34 (page 591)

**Interface**:
```matlab
f = ellipfast(m, kind)
% m = k^2 (the parameter, not the modulus k)
% kind: 1 for K, 2 for E  (or 'K' / 'E')
% Accuracy: error < 2e-8 for 0 <= m <= 1
```

**C implementation** — argument transformation: `y = 1 - m` (complementary parameter):

**K(m)** (kind=1), valid for 0 <= y <= 1:
```c
K(m) = (((0.01451196212*y + 0.03742563713)*y + 0.03590092383)*y + 0.09666344259)*y + 1.38629436112
      -((((0.00441787012*y + 0.03328355346)*y + 0.06880248576)*y + 0.12498593597)*y + 0.50000000000)*log(y)
```

**E(m)** (kind=2), with special case E(1) = 1 (y=0):
```c
E(m) = (((0.01736506451*y + 0.04757383546)*y + 0.06260601220)*y + 0.44325141463)*y + 1.00000000000
      -((((0.00526449639*y + 0.04069697526)*y + 0.09200180037)*y + 0.24998368310)*y)*log(y)
```

Note: `y = 1-m = 1-k^2`. At `y=0` (i.e., m=1, k=1), K diverges as `-0.5*log(y)`, and E=1. At `y=1` (m=0, k=0), K = π/2 ≈ 1.5707963, E = π/2 ≈ 1.5707963. The log(y) term carries the logarithmic singularity of K near k=1.

The MATLAB wrapper does `if ischar(kind), kind = 13.5 - double(kind)/6; end` which maps 'K'→1, 'E'→2.

In `greenem.c` the same polynomial coefficients appear but parameterized as `x = d*is = 1-k^2 = y`, so the two implementations are identical.

### 2.6 `green3d.m` — 3D mutual inductance

**Signature**:
```matlab
M = green3d(r1, z1, phi1, r2, z2, phi2)
```

**Cases**:

1. Both phi1 and phi2 empty: falls back to `greenem('mut',r1,z1,r2,z2)` (axisymmetric).

2. phi1 empty, phi2 provided: one loop is a circle, the other is a 3D polygon. Uses `green2d5`, which computes:
   ```
   M = greenem('mut', r_circle, z_circle, r2_polygon, z2_polygon) * dphi2 / (2*pi)
   ```
   (numerical integration of the axisymmetric Green's function over the 3D path)

3. Both phi1 and phi2 provided: full 3D Neumann formula via piecewise-linear segments:
   $$M = \frac{\mu_0}{4\pi} \sum_i \sum_j \frac{\mathbf{dl}_1^{(i)} \cdot \mathbf{dl}_2^{(j)}}{|\mathbf{x}_2^{(j)} - \mathbf{x}_1^{(i)}|}$$
   
   Implemented in `greenxyz`:
   $$M = \sum_i \sum_j \frac{dx_1 dx_2 + dy_1 dy_2 + dz_1 dz_2}{\sqrt{(x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2}} \times 10^{-7}$$
   
   where $10^{-7} = \mu_0/(4\pi)$.

**Input format**: rows of each matrix are the polygon vertices; NaN rows are used as padding for variable-length polygons.

---

## Section 3 — Vessel Eigenmodes (`vveig`)

### 3.1 Purpose

`vveig` compresses a vessel model consisting of Nv individual filaments into Nu dominant electromagnetic modes. This reduces the dimensionality of the passive conductor model for real-time control applications.

### 3.2 Algorithm

**Signature**:
```matlab
[Tve, Le, Re, sym] = vveig(Mvv, Rv, type)
```

**Inputs**:
- `Mvv`: [Nv x Nv] mutual inductance matrix of vessel filaments (symmetric, positive definite)
- `Rv`: [Nv x 1] resistance vector of each filament (in Ohm)
- `type`: `'R'` (resistance normalization, default) or `'M'` (inductance normalization)

**Outputs**:
- `Tve`: [Nv x Nv] transformation matrix; `IV = Tve * IE` converts eigenmode currents IE to filament currents IV
- `Le`: eigenmode inductances (array for type='R'; scalar=1 for type='M')
- `Re`: eigenmode resistances (scalar for type='R'; array for type='M')
- `sym`: [1 x Nv] logical vector, 1 for symmetric (even-parity) eigenmodes

**Mathematical procedure**:

1. Form the resistance-normalized matrix:
   $$MM = \frac{M_{vv}}{\sqrt{R_v} \sqrt{R_v}^T} = D_{R_v}^{-1/2} M_{vv} D_{R_v}^{-1/2}$$
   
   where $D_{R_v}^{-1/2} = \text{diag}(R_v)^{-1/2}$.

2. Symmetrize and diagonalize:
   $$[T_{ve,raw}, \Gamma_e] = \text{eig}\left(\frac{MM + MM^T}{2}\right)$$
   $$\gamma_e = \text{diag}(\Gamma_e)$$
   
   Eigenvalues $\gamma_e$ are the time constants $L_e/R_e$ of the RL circuits in the normalized basis.

3. Sort eigenvalues in descending order (longest time constants first):
   ```matlab
   k = flipud(isort(gammae))
   ```

4. Two normalization choices:

   **type='R'** (resistance normalization):
   ```
   Re = 1 / sum(1/Rv)         % parallel resistance = scalar
   Le = Re * gamma_e          % eigenmode inductances = array
   ```

   **type='M'** (inductance normalization):
   ```
   Re = 1 / gamma_e           % eigenmode resistances = array
   Le = 1                     % eigenmode inductance = scalar
   ```

5. Build the transformation matrix:
   ```matlab
   Tve = (Qve ./ sqrt(Rv)) * sqrt(diag(Re))
   ```
   where `Qve` is the eigenvector matrix sorted by descending gamma.

6. Parity detection: eigenmodes are flagged as symmetric (even) if the sum of squared values is larger in the upper half than the "antisymmetric" combination. Sign conventions are applied so that symmetric modes have Tve(1,:) > 0 and antisymmetric modes have Tve(2,:) > Tve(1,:).

### 3.3 Key identities verified by tests

For type='R':
```
diag(Le) == Tve' * Mvv * Tve      (eigenmode mutual inductance matrix)
Re * I   == Tve' * diag(Rv) * Tve (eigenmode resistance matrix)
```

For type='M':
```
Le * I        == Tve' * Mvv * Tve
diag(Re)      == Tve' * diag(Rv) * Tve
```

The eigenvalues Le/Re (type='R') equal the eigenvalues of diag(Rv)\Mvv sorted in descending order.

### 3.4 Truncation (not in vveig itself)

`vveig` always returns all Nv eigenmodes. Truncation to Nu < Nv modes is done externally by taking only the first Nu columns of Tve (corresponding to the Nu largest time constants).

---

## Section 4 — B-splines

### 4.1 `bspbase` — Basis matrix

**Signature**:
```matlab
B = bspbase(t, k, x [, d])
```

**Returns**: [nb x nx] matrix where `nb = numel(t) - k` is the number of B-spline basis functions, `nx = numel(x)`.
- `B(j, i)` is the j-th B-spline basis function of order k evaluated at x(i).
- d-th derivative is returned if d > 0 (default d=0).

**Requirements**: x must be sorted (ascending); the function sorts internally if needed. x must lie in `[t(k), t(end-k+1)]`.

**Convention**: the knot sequence t is of length `nb + k`. Boundary knots are NOT automatically repeated — the caller must provide a properly augmented knot sequence if needed.

**C implementation** (`bspbase.c`) — de Boor recurrence:

The algorithm finds the knot span j such that `T[j] <= X[i] <= T[j+1]`, then applies the standard Cox–de Boor recurrence:

```c
// Initialize: b[0] = 1 (only the j-th basis is nonzero at step 0)
b[0] = 1;
// Recurrence for order (builds up from order 1 to k-d)
for r = -1 down to d-k+1:
    for s = r to -1:
        o = (X[i] - T[j+s+1]) / (T[j+s+1-r] - T[j+s+1])  // or 0 if denominator=0
        b[s]   += (1 - o) * b[s+1]
        b[s+1] *= o
// Differentiation pass (for d > 0)
for r = d-k up to -k+1:
    for s = r to -1:
        o = kk / (T[j+s+1-r] - T[j+s+1])   // kk counts down from k-1
        b[s]   -= o * b[s+1]
        b[s+1] *= o
```

The result is stored column-major: `B[j + nbase*i]` is the j-th basis at the i-th point.

### 4.2 `bspsum` — Evaluate spline at points

**Signature**:
```matlab
y = bspsum(t, C, x [, d] [, flag] [, npar])
```

- `t`: knot sequence (length = nb + k, where nb = size(C,1))
- `C`: coefficient array, first dimension = nb (number of basis functions)
- `x`: evaluation points
- `d`: derivative order (default 0)
- `flag`: if 1, C is [nb x numel(x) x DIM] and result is [numel(x) x DIM]
- `npar`: number of OpenMP threads

**Output size**: `[size(x), size(C)(2:end)]` or `[numel(x), size(C)(2:end)]` for vector x.

**Spline order**: `k = numel(t) - size(C,1)`

The function evaluates `f(x_i) = sum_j C_j * B_j(x_i)` where B_j are B-splines of order k.

The MEX file (`bspsummex`) uses OpenMP for parallelism.

### 4.3 `bsp2deval` — 2D tensor-product B-spline

**Signature**:
```matlab
[f, gr, gz, hrr, hrz, hzz, d3rrr, ...] = bsp2deval(rk, zk, C, r, z, d, npar)
```

Evaluates a 2D function `f(r,z) = sum_{i,j} C_{ij} * B_i^r(r) * B_j^z(z)` and all its derivatives up to order d (max 6).

**Algorithm**: separable evaluation by applying bspsum twice — first along r for each z-slice of coefficients, then along z:
```matlab
cr0 = bspsum(rk, C, r, 0, 0, npar).'    % [nr x ncoefz] intermediate
f   = reshape(bspsum(zk, cr0, z, 0, 1, npar), sz)
```

Mixed derivatives like hrz = ∂²f/∂r∂z are obtained by differentiating the r-pass to get cr1, then z-passing cr1 with d=1.

For regular knot sequences, `cs2deval` (described below) is preferred for speed.

### 4.4 `csdec` — Cubic B-spline decomposition (fitting)

**Signature**:
```matlab
[M, taug, b] = csdec(t, x [, xx] [, ec] [, d])
```

Builds the matrix M such that `M * Y` gives the B-spline coefficients for a cubic spline that best fits (or exactly interpolates when t=x) the data values Y at points x for knot sequence t.

**Augmented knot sequence**: internally, csdec augments t with 3 knots on each side:
```matlab
dt = (t(end) - t(1)) / (nt - 1)
taug = [(-3:-1)'*dt + t(1); t(:); t(end) + (1:3)'*dt]
```
This creates a B-spline basis of order 4 (cubic) with `nt + 2` basis functions.

**Edge conditions** (ec argument):

| Code | Name | Constraint |
|---|---|---|
| `'a'` | not-a-knot (default) | f''' continuous at t(2) and t(end-1) |
| `'n'` | natural | f''(t(1)) = f''(t(end)) = 0 |
| `'p'` | periodic | f'(t(1)) = f'(t(end)), f''(t(1)) = f''(t(end)) |
| `'s'` | symmetric | f'(t(1)) = f''(t(end)) = 0 |
| `'t'` | [ill-conditioned, not recommended] | f'(t(1)) = f'''(t(1)) = 0 |
| `'u'` | over-constrained | f'(t(1)) = f'''(t(1)) = f''(t(end)) = 0 |
| `'v'` | over-constrained | f'(t(1)) = f'''(t(1)) = f'''(t(end)) = 0 |
| `'w'` | mixed | f'(t(1)) = 0, not-a-knot at t(end) |

**Constraint regression method**: the edge conditions T·c = 0 are enforced via QR decomposition. Let `b = bspbase(taug, 4, x)'` be the basis matrix. Then:
```matlab
[q, ~] = qr(T)
q = q[:, size(T,2)+1:end]   % null space of T
b = b * q                    % constrained basis
M = q / b                   % = q * pinv(b)  (least-squares fit matrix)
```

**2D splines**: use `M_r = csdec(t_r, x_r)` and `M_z = csdec(t_z, x_z)`, then the coefficient matrix is `C = M_r * Z_data * M_z'`.

**Optional evaluation**: if `xx` is provided (numeric), `csdec` also evaluates the spline at `xx` by calling `bspsum(taug, M, xx, d)`.

### 4.5 `cs2deval.c` — Fast regular 2D cubic spline evaluation

**Signature** (C):
```c
void cs2deval(double *V, double *Vx, ..., double *Tx, int nTx,
              double *Ty, int nTy, double *C, int nCx, int nCy,
              double *xp, double *yp, int np, int dmax, double idx, double idy)
```

For **regular** knot sequences (uniform spacing), this is faster than `bsp2deval` because it can compute the knot interval directly by division rather than search.

The helper `csint(c0,c1,c2,c3, rho, d, fac)` computes the 4 cubic B-spline coefficients for fractional position `rho` within a cell [0,1]:

- d=0 (evaluation): standard cubic B-spline basis in Bernstein form
- d=1 (1st derivative): derivative coefficients, scaled by `fac = idx`
- d=2 (2nd derivative): second derivative coefficients, scaled by `fac = idx^2`
- d=3 (3rd derivative): constant coefficients [1,-3,3,-1]*fac

The helper `adotc` computes the matrix-vector product: given row vector `[a0,a1,a2,a3]` and a 4x4 submatrix of C at (ix,iy), returns 4 intermediate column coefficients.

Evaluation at (xp[i], yp[i]):
```
ix = int((xp[i] - Tx[0]) * idx)    // cell index in x
iy = int((yp[i] - Ty[0]) * idy)    // cell index in y
dx = (xp[i] - Tx[ix]) * idx        // fractional position in cell [0,1]
dy = (yp[i] - Ty[iy]) * idy
V[i] = sum_{p,q} a_p * b_q * C[ix+p-3, iy+q-3]
```

---

## Section 5 — Cubic Splines (legacy)

### 5.1 `csmod` — Cubic spline interpolant (pp-form)

**Signature**:
```matlab
output = csmod(x, y [, ec] [, xx])
```

Computes the cubic spline interpolant to data (x, y) in MATLAB pp-form using the classical tridiagonal system for slopes. If `xx` is omitted, returns the pp-form; otherwise evaluates at `xx`.

**Algorithm**: sets up the tridiagonal linear system for the first derivatives (slopes) at each knot, with edge conditions applied at row 1 and row n:

```
A * b = B    (b = slopes at knots)
```

where A is tridiagonal: `A(i,i-1) = dx(i), A(i,i) = 2*(dx(i)+dx(i-1)), A(i,i+1) = dx(i-1)`.

After solving for slopes b, converts to pp-form:
```
d = (b(1:n-1) + b(2:n) - 2*dy) / dx
c = (dy - b(1:n-1)) / dx - d
```
giving the local cubic `p(s) = d(i)*s^3 + c(i)*s^2 + b(i)*s + y(i)` for s in [0, dx(i)].

**Edge conditions**:

| Code | Constraint |
|---|---|
| `'a'` | not-a-knot: 3rd derivative continuous at x(2) and x(n-1) |
| `'n'` | natural: 2nd derivative zero at endpoints |
| `'s'` | symmetric: slope=0 at x(1), not-a-knot at x(n) |
| `'p'` | periodic: f'(x(1)) = f'(x(n)), f''(x(1)) = f''(x(n)) |

Special cases: n=2 → linear, n=3 → parabola (single piece from x(1) to x(3)).

### 5.2 `csder` — Cubic spline from derivative data

**Signature**:
```matlab
output = csder(x, y, yp [, xx])
```

Builds a cubic spline interpolant with both function values y and derivatives yp specified at knots x. Returns pp-form or evaluates at xx.

**Formula**: given y, yp (slopes), computes:
```
dy = diff(y) / dx^2            (scaled differences)
d  = (yp(2:n) + yp(1:n-1)) / dx - 2*dy   (cubic coefficient)
c  = dy - yp(1:n-1) / dx - d              (quadratic coefficient)
d  = d / dx                               (normalize)
```

This is a Hermite cubic interpolation: each piece is uniquely determined by y and yp at both endpoints.

### 5.3 `pp2bsp` — pp-form to genlib B-spline format

Listed in Contents.m but not read in detail. [GUESS] Converts MATLAB's native pp-form (from spline/csapi) to the (t,C) knot+coefficient format used by bspsum.

---

## Section 6 — Other Utilities

### 6.1 `sti` — Frequency domain transfer function fitting

**Signature**:
```matlab
[A, B, RAB, Hfit, iter] = sti(w, H, nA, nB [, maxiter, tol, W])
```

Fits a rational transfer function `B(s)/A(s)` (numerator degree nB-1, denominator degree nA) to experimental frequency-response data H(w) using an iterative weighted least-squares algorithm (Sanathanan-Koerner / Steiglitz-McBride style).

**Algorithm**:
1. Normalize `w` by `wnorm = max(w)`.
2. Build regressor matrix phi from powers of `jw` (up to nA-1 for denominator, nB for numerator).
3. Iterate (up to maxiter=20 times):
   - Weight each row by `1/A(jw)`.
   - Solve `real(phi'*phi) \ real(phi'*phy)` for coefficients (uses real part to enforce real-valued result).
   - Update A.
4. Return B coefficients, covariance matrix `RAB`, and fitted response `Hfit`.

**Polynomial convention**: uses descending coefficient order (MATLAB `polyval` convention). E.g., nA=1, nB=1 gives `H(s) = B[1] / (A[2] + s*A[1])`.

**Output covariance**: `RAB = (residuals'*residuals)/(nw*nH) * inv(real(phi'*phi))`

### 6.2 `xmtimes` — Extended matrix multiply

**Signature**:
```matlab
z = xmtimes(x, y)
```

Generalization of matrix multiplication to higher-dimensional arrays. Contracts the last dimension of x with the first dimension of y (inner product):
```matlab
z = reshape(reshape(x, [prod(nx(1:end-1)), nx(end)]) * ...
            reshape(y, [ny(1), prod(ny(2:end))]), [nx(1:end-1), ny(2:end)])
```

Equivalent to `numpy.tensordot(x, y, axes=1)` or `np.einsum('...k,k...->...',x,y)`. "Extended" means the inputs can be N-D, not just 2-D matrices.

### 6.3 `rpdcov` — Residue-pole covariance

**Signature**:
```matlab
[r, p, d, Rrpd, rsd, psd, dsd] = rpdcov(A, B [, RAB])
```

Decomposes B(s)/A(s) into residues `r`, poles `p`, and polynomial division rest `d` (via MATLAB `residue`). If the coefficient covariance matrix `RAB` is provided, propagates it through the Jacobian to give the residue-pole-rest covariance `Rrpd`, using `rpdjacobian`.

### 6.4 `rpdjacobian` — Residue-pole-rest Jacobian

**Signature**:
```matlab
[J, Jinv] = rpdjacobian(A, B, r, p, d [, normalised])
```

Computes the Jacobian of the mapping `[A, B] → [r, p, d]` (polynomial coefficients to partial fractions). Uses Toeplitz convolution matrices built by `toepconv`. The Jacobian and its inverse are obtained by solving two linear systems.

### 6.5 `toepconv` — Toeplitz convolution matrix

**Signature**:
```matlab
T = toepconv(A, N)
```

Creates a `(N+1) x (length(A)+N)` Toeplitz matrix such that `c * T` gives the polynomial product of c (degree N) with A (degree length(A)-1). Equivalent to convolution implemented as matrix multiplication. Useful for building the Jacobian of polynomial operations.

### 6.6 `resamp` — Time signal resampling

**Signature**:
```matlab
[R1, R2, ..., TR] = resamp(X1,T1, X2,T2, ..., [T], [MODE])
```

Resamples multiple signal-time pairs onto a common time base. Mode is a 3-character string:
- Character 1: `'f'`=fastest (smallest period) or `'s'`=slowest common period
- Character 2: `'r'`=restricted (intersection), `'e'`=extended (union), `'p'`=prolonged with first/last value
- Character 3: `'n'`=nearest, `'l'`=linear, `'c'`=cubic, `'s'`=spline interpolation

Uses MATLAB `interp1` internally.

### 6.7 `simple_polygon` — Polygon self-intersection test

**Signature**:
```matlab
isSimple = simple_polygon(P)
% P: [n x 2] vertex coordinates, no repeated first/last vertex
```

Implements the **Shamos-Hoey algorithm** (O(n log n)) to detect whether a polygon has any self-intersections. Processes edges in sweep-line order; returns `true` if no intersections found.

Intersection test (`intersects_`): checks if two edges properly cross using the cross-product sign test.

### 6.8 Index search utilities

All three (`ifloor`, `iceil`, `iround`) call the same MEX file `iflooriceiliroundmex(t, x, int32(mode))` with mode=1,2,0 respectively. They require `t` to be monotonically increasing and `x` to be monotonically increasing (or pass a 3rd argument to sort x automatically).

| Function | Returns | Python equivalent |
|---|---|---|
| `ifloor(t, x)` | index of largest t(k) <= x | `np.searchsorted(t, x, side='right') - 1` |
| `iceil(t, x)` | index of smallest t(k) >= x | `np.searchsorted(t, x, side='left')` |
| `iround(t, x)` | index of nearest t(k) to x | `np.argmin(abs(t - x))` per element |

### 6.9 Remaining utility one-liners

| Function | Description | Python equivalent |
|---|---|---|
| `isodd(x)` | `logical(rem(x,2))` — true where x is odd | `x % 2 != 0` |
| `isort(x)` | returns sort indices only (second output of `sort`) | `np.argsort(x)` (note: 0-indexed) |
| `imax(x)` | index of maximum | `np.argmax(x)` (0-indexed) |
| `imin(x)` | index of minimum | `np.argmin(x)` (0-indexed) |
| `isnans(x)` | true for **signaling** NaN (not quiet NaN) | No direct numpy equivalent; MEX-based |
| `ifc(c,a,b)` | ternary: returns a if c else b; supports lazy lambdas | `a if c else b` |
| `ppzero(pp, z, d, tol)` | real zeros of pp-form polynomial minus z | scipy.interpolate equivalent not trivial |
| `appendstruct(S1,S2)` | merge two structs, S2 fields overwrite S1 | `{**d1, **d2}` for dicts |
| `struct2str(S, n)` | generate MATLAB code string defining a struct | `repr(dict)` roughly |
| `xmtimes(x,y)` | N-D inner product on last/first dims | `np.tensordot(x,y,axes=1)` |

---

## Section 7 — Python / NumPy Rewrite Notes

### 7.1 `greenem` → Python

The core mutual inductance computation. The exact polynomial approximations should be copied verbatim.

```python
import numpy as np

MU0 = 1.256637061435917e-6  # 4e-7*pi
MU1 = 2e-7                  # mu0/(2*pi)

def ellipK(m):
    """K(m) via Abramowitz 17.3.34. Input m = k^2, accuracy < 2e-8."""
    y = 1.0 - m
    return ((((0.01451196212*y + 0.03742563713)*y + 0.03590092383)*y
              + 0.09666344259)*y + 1.38629436112
            - ((((0.00441787012*y + 0.03328355346)*y + 0.06880248576)*y
                 + 0.12498593597)*y + 0.50000000000) * np.log(y))

def ellipE(m):
    """E(m) via Abramowitz 17.3.34. Input m = k^2, accuracy < 2e-8."""
    y = 1.0 - m
    if np.any(y == 0): return np.where(y == 0, 1.0, _ellipE_inner(y))
    return _ellipE_inner(y)

def _ellipE_inner(y):
    return ((((0.01736506451*y + 0.04757383546)*y + 0.06260601220)*y
              + 0.44325141463)*y + 1.0
            - ((((0.00526449639*y + 0.04069697526)*y + 0.09200180037)*y
                 + 0.24998368310)*y) * np.log(y))

def greenem_mut(r1, z1, r2, z2):
    """
    Mutual inductance [H] between circular filaments.
    r1, z1: [N1] arrays; r2, z2: [N2] arrays.
    Returns [N1 x N2] matrix.
    """
    r1 = np.asarray(r1).ravel()
    z1 = np.asarray(z1).ravel()
    r2 = np.asarray(r2).ravel()
    z2 = np.asarray(z2).ravel()

    # Broadcast to [N1, N2]
    R1 = r1[:, None]; Z1 = z1[:, None]
    R2 = r2[None, :]; Z2 = z2[None, :]

    rsum = R1 + R2
    h    = Z1 - Z2
    s    = rsum**2 + h**2        # (r1+r2)^2 + (z1-z2)^2
    d    = (R1-R2)**2 + h**2     # (r1-r2)^2 + (z1-z2)^2
    m    = 4*R1*R2 / s           # k^2
    x    = d / s                 # 1 - k^2 (complementary parameter)
    K    = ellipK(m)
    E    = ellipE(m)
    v0   = (1.0 - 0.5*m)*K - E
    u0   = np.sqrt(s)
    return MU0 * u0 * v0
```

**Pitfall**: MATLAB `greenem` returns a [N1 x N2] matrix where the *rows* index source 1 and *columns* index source 2. NumPy default broadcasting follows the same convention if you use `r1[:,None]` and `r2[None,:]`.

**Self-inductance**: use `MU0 * r * (np.log(8*r/a) - 1.75)` for circular conductors.

**Alternative**: `scipy.special.ellipk(m)` and `scipy.special.ellipe(m)` give K and E to full double precision, but are slower than the polynomial approximation for large batches.

### 7.2 `vveig` → Python

```python
import numpy as np

def vveig(Mvv, Rv, type_='R'):
    """
    Vessel eigenmode decomposition.
    Mvv: [Nv, Nv] mutual inductance matrix
    Rv:  [Nv] resistance vector
    type_: 'R' or 'M'
    Returns: Tve [Nv,Nv], Le, Re, sym [Nv]
    """
    nv = len(Rv)
    sqrtRv = np.sqrt(Rv)
    MM = Mvv / np.outer(sqrtRv, sqrtRv)
    MM = (MM + MM.T) / 2  # ensure symmetry

    gammae, Qve = np.linalg.eigh(MM)  # eigh for symmetric matrix
    # eigh returns ascending order; we want descending (longest time constants first)
    idx = np.argsort(gammae)[::-1]
    gammae = gammae[idx]
    Qve = Qve[:, idx]

    if type_ == 'R':
        Re = 1.0 / np.sum(1.0 / Rv)
        Le = Re * gammae
    elif type_ == 'M':
        Re = 1.0 / gammae
        Le = 1.0

    Tve = (Qve / sqrtRv[:, None]) * np.sqrt(Re)

    # Parity detection
    half = (nv + 1) // 2
    sym = (np.sum((Tve[:half, :] + Tve[nv-1:nv//2-1:-1, :])**2, axis=0) > 1e-6)
    fact = np.sign(Tve[0, :])
    fact[~sym] = np.sign(Tve[1, ~sym] - Tve[0, ~sym])
    Tve = Tve * fact[None, :]

    return Tve, Le, Re, sym
```

**Pitfall**: MATLAB `eig` does NOT sort eigenvalues; `np.linalg.eigh` sorts ascending. The `isort`/`flipud` step in MATLAB must be replicated explicitly. Use `eigh` (not `eig`) for symmetric matrices.

### 7.3 `bspbase` → Python

```python
from scipy.interpolate import BSpline
import numpy as np

def bspbase(t, k, x, d=0):
    """
    B-spline basis matrix. Returns [nb, nx] matrix.
    nb = len(t) - k, nb > 0 required.
    """
    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float).ravel()
    nb = len(t) - k
    B = np.zeros((nb, len(x)))
    for j in range(nb):
        c = np.zeros(nb); c[j] = 1.0
        spl = BSpline(t, c, k-1, extrapolate=False)
        if d == 0:
            B[j, :] = spl(x)
        else:
            B[j, :] = spl.derivative(d)(x)
    return B
```

**Pitfall**: scipy BSpline uses `order = degree = k-1`, while genlib uses `k = degree + 1 = order`. So `BSpline(t, c, k-1, ...)` is correct.

**Pitfall**: The knot sequence convention in genlib uses a linearly extended augmented sequence in `csdec`. When replicating csdec, you must build taug with the same 3-knot extension on each side.

### 7.4 `csdec` → Python

```python
import numpy as np
from scipy.linalg import qr, lstsq

def csdec(t, x, ec='a'):
    """
    Cubic B-spline decomposition.
    Returns M [nb x nx] such that M @ y gives B-spline coefficients.
    """
    t = np.asarray(t, dtype=float)
    nt = len(t)
    dt = (t[-1] - t[0]) / (nt - 1)
    # Augment knot sequence
    taug = np.concatenate([np.arange(-3, 0)*dt + t[0], t, t[-1] + np.arange(1, 4)*dt])
    order = 4
    b = bspbase(taug, order, x)  # [nb, nx]

    # Build constraint matrix T (edge conditions)
    # ... (see csdec.m for each ec case) ...
    # [GUESS] Implementation follows the same logic as the MATLAB code.

    Q, _ = qr(T_matrix)                # QR of constraint matrix
    q = Q[:, T_matrix.shape[1]:]       # null space
    b = b.T @ q                        # constrained [nx, reduced_nb]
    M = np.linalg.lstsq(b, q.T, rcond=None)[0]  # [reduced_nb, nx]
    return M, taug, b
```

**Pitfall**: MATLAB `qr(T)` without column permutation gives a full QR; the null space columns are `Q[:, rank:]`. In Python use `scipy.linalg.qr(T, mode='economic')` carefully.

### 7.5 `csmod` → Python

```python
from scipy.interpolate import CubicSpline

def csmod(x, y, ec='a', xx=None):
    """Wrapper: use scipy CubicSpline."""
    bc = {'a': 'not-a-knot', 'n': 'natural', 'p': 'periodic'}[ec]
    cs = CubicSpline(x, y, bc_type=bc)
    return cs if xx is None else cs(xx)
```

**Pitfall**: scipy `CubicSpline` with `bc_type='not-a-knot'` matches csmod's `'a'` mode. The `'s'` (symmetric) mode has no direct scipy equivalent.

### 7.6 `ellipfast` vs scipy

`scipy.special.ellipk(m)` and `scipy.special.ellipe(m)` accept `m = k^2` (parameter convention matching genlib) and are full double precision. They are interchangeable for correctness; the polynomial approximation in genlib is faster for vectorized batch computation where < 2e-8 accuracy is sufficient.

### 7.7 Common pitfalls

1. **MATLAB is 1-indexed; Python is 0-indexed**. All index functions (ifloor, iceil, iround, imax, imin, isort) return 1-based indices in MATLAB. Subtract 1 in Python.

2. **Column-major vs row-major**: MATLAB stores arrays column-major; numpy stores row-major by default. Beware when reshaping C coefficient arrays from bspbase/bspsum.

3. **`m` vs `k` vs `y`**: genlib uses `m = k^2` (the parameter). Abramowitz uses both conventions. scipy uses `m = k^2`. Be explicit.

4. **Self-inductance diagonal**: `greenem('mut', R, Z)` fills the diagonal with self-inductances automatically when conductor size A is provided. In Python you must compute the diagonal separately and insert it.

5. **`vveig` eigenvalue ordering**: MATLAB `eig` returns unsorted eigenvalues; `isort` then `flipud` sorts descending. `np.linalg.eigh` returns ascending — you must reverse.

6. **`bspsum` flag=1 mode**: this is a non-standard calling convention for applying coefficients element-wise; must be replicated carefully.

7. **Thread count in greenemmex**: the MEX uses `maxNumCompThreads` for auto-parallelism. In Python use `concurrent.futures` or numpy's implicit BLAS parallelism for the outer loop.

8. **`csmod` returns pp-form by default**: the pp-form is a MATLAB-specific piecewise polynomial format. The Python equivalent is `scipy.interpolate.CubicSpline` or `scipy.interpolate.PPoly`.
