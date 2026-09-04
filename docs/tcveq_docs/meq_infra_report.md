# MEQ Infrastructure Layer — Technical Reference

**Codebase**: MEQ (Matlab EQuilibrium) Toolbox v6.10, Swiss Plasma Center EPFL  
**Purpose**: Axisymmetric tokamak plasma equilibrium calculations  
**Solvers built on this layer**: LIUQE (reconstruction), FBT (free-boundary equilibrium), FGE (evolutive forward solver)

---

## 1. The MEQ Data Model

### 1.1 Top-level structures

MEQ organises all computation-time data into three structures. In most top-level solver calls these are named `L`, `LY`, and `LX`, but the actual field names are identical regardless of variable name.

**`L` — the "setup" or "ancillary" structure**  
Built by `meqc(P,G)`. Contains everything that does not change between time slices: grid metadata, Poisson solver coefficients, interpolation matrices, normalisation scalings, and base-function state. Its two principal sub-structures are:

- `L.G` — geometry and electromagnetic parameters (Green's function matrices). Built by `meqg(G,P)`.
- `L.P` — configuration parameters (algorithm switches, numerical tolerances, grid sizes). Built by solver-specific parameter functions (`liup`, `fbtp`, `fgep`).

Everything else in `L` is a pre-computed derived quantity for speed: grid spacings, inverse spacings, FE matrices, Poisson solver coefficients, etc.

**`LY` — the "state" or "output" structure**  
Contains the time-varying equilibrium result for one or more time slices. Fields are packed so that the last dimension is always the time index. Key fields:

| Field | Description | Dimensions |
|---|---|---|
| `LY.t` | time array | `[1, nt]` |
| `LY.shot` | shot number | scalar or `[1, nt]` |
| `LY.Fx` | poloidal flux on x-grid (Wb) | `[nzx, nrx, nt]` |
| `LY.Iy` | plasma toroidal current per cell on y-grid (A) | `[nzy, nry, nt]` |
| `LY.Opy` | plasma domain label map on y-grid (int8, 0=outside, 1..nD) | `[nzy, nry, nt]` |
| `LY.Ia` | active coil currents (A) | `[na, nt]` |
| `LY.Iu` | generalised vessel mode currents (A) | `[nu, nt]` |
| `LY.Iv` | vessel filament currents (A), `Iv = L.G.Tvu * Iu` | `[nv, nt]` |
| `LY.ag` | base-function coefficients for p' and TT' | `[ng, nt]` |
| `LY.rA`, `LY.zA`, `LY.FA` | magnetic axis position and flux | `[nD, nt]` |
| `LY.rB`, `LY.zB`, `LY.FB` | LCFS point and flux | `[nD, nt]` |
| `LY.lX` | true if domain is X-point diverted | `[nD, nt]` |
| `LY.lB` | true if LCFS was found | `[1, 1]` or `[1, nt]` |
| `LY.F0`, `LY.F1` | flux at axis (F0) and LCFS (F1) per domain | `[nD, nt]` |
| `LY.Ip` | total plasma current (A) | `[1, nt]` |
| `LY.rBt` | vacuum R×Bₜ (T·m) | `[1, nt]` |
| `LY.Wk` | kinetic (thermal) energy (J) | `[1, nt]` |
| `LY.Wp` | poloidal field energy (J) | `[1, nt]` |
| `LY.Wt` | plasma toroidal field energy (J) | `[1, nt]` |
| `LY.Vp` | plasma volume (m³) | `[1, nt]` |
| `LY.bp`, `LY.bt`, `LY.li`, `LY.mu` | betas, internal inductance, diamagnetism | `[1, nt]` |
| `LY.qA` | safety factor at axis | `[nA, nt]` |
| `LY.Bm` | simulated magnetic probe signals | `[nm, nt]` |
| `LY.Ff` | simulated flux loop signals (Wb) | `[nf, nt]` |
| `LY.PpQ`, `LY.TTpQ` | p'(ρ), TT'(ρ) profiles | `[nQ, nD, nt]` or `[nQ, nt]` when `nD=1` |
| `LY.PQ`, `LY.TQ`, `LY.iTQ` | p(ρ), T(ρ)=RBₜ, 1/T(ρ) profiles | same |
| `LY.aq` | flux-surface distance array (for contour post-processing) | `[noq, npq, nD, nt]` |
| `LY.aW` | wall-gap distances | `[nW, nFW, nt]` |
| `LY.nA`, `LY.nB`, `LY.nX` | number of axes, bounded domains, X-points | int32 |
| `LY.Brx`, `LY.Bzx` | Br, Bz on x-grid (optional) | `[nzx, nrx, nt]` |
| `LY.Fz` | flux on z-grid (optional extended post-processing) | `[nzz, nrz, nt]` |

**`LX` — the "input/measurement" structure**  
Identical field layout to `LY`, but contains only measured/external inputs. For LIUQE: `LX.Bm` (measured probe signals), `LX.Ff` (measured flux loops), `LX.Ia` (coil currents). For FGE: also `LX.Va` (coil voltages). The function `meqxk(LX,k)` extracts a single or subset of time indices from an LX or LY structure.

### 1.2 Variable naming conventions (from `meqhelp.m`)

The naming is systematic — every name is: optional prefix + root letter(s) + suffix/system letter.

**Root letters (capital = matrix/vector quantity)**:
- `B` magnetic field; `F` flux (ψ); `I` measured current; `J` fitted current; `M` mutual inductance; `R` resistance; `T` transformation matrix; `U` voltage; `Q` flux-surface integrals

**Suffix/system letter (identifies the grid or current type)**:
- `x` — x-grid (computational grid)
- `y` — y-grid (plasma current grid, interior of x-grid)
- `z` — z-grid (extended post-processing grid)
- `a` — active coils
- `u` — generalised vessel modes
- `v` — vessel filaments
- `m` — magnetic probes
- `f` — flux loops
- `b` — boundary of x-grid
- `l` — limiter
- `A` — magnetic axis
- `B` — LCFS
- `X` — X-point
- `D` — per-domain quantity

**Prefix letters**:
- `n` — number (count); `d` — mesh step; `i` — inverse; `dr`, `dz` — partial derivatives in R and Z; `l` — logical index (boolean); `k` — 1-based index; `j` — 0-based index; `T` — transformation matrix

**Dimension label arrays** (cell arrays of strings describing each element):
- `G.dima` — active coil labels
- `G.dimu` — generalised vessel mode labels
- `G.dimv` — vessel filament labels
- `G.dimf` — flux loop labels
- `G.dimm` — magnetic probe labels

---

## 2. Coordinate Systems and Grids

### 2.1 Tokamak cylindrical coordinate system

MEQ uses the standard tokamak cylindrical coordinates `(R, Z, φ)`:
- `R` — major radius (horizontal, metres). Symbol `r` in code.
- `Z` — vertical coordinate (metres). Symbol `z` in code.
- `φ` — toroidal angle (not explicitly gridded; axisymmetry assumed).

The poloidal flux function $\psi$ is stored as `Fx(Z,R)` — **MATLAB column-major order means the first index is Z and the second is R**. This is the conventional choice: `Fx(iz, ir)`. Every 2D grid array follows this `(nz, nr)` shape convention.

### 2.2 The x-grid (computational grid)

The x-grid is the primary rectangular grid covering the interior of the limiter (first-wall contour). It is uniform and equispaced.

| Field | Location | Description |
|---|---|---|
| `L.rx` | `L` | R coordinates, column vector, `[nrx, 1]` |
| `L.zx` | `L` | Z coordinates, column vector, `[nzx, 1]` |
| `L.nrx` | `L` | number of R grid points |
| `L.nzx` | `L` | number of Z grid points, must satisfy `(nzx-1) = 2^k` for FFT in Poisson solver |
| `L.drx` | `L` | R grid spacing |
| `L.dzx` | `L` | Z grid spacing |
| `L.idrx`, `L.idzx` | `L` | inverses (precomputed for speed) |
| `L.dsx` | `L` | cell area `drx * dzx` |
| `L.idsx` | `L` | inverse cell area |
| `L.nx` | `L` | total number of x-grid points `nrx * nzx` |
| `L.rrx`, `L.zzx` | `L` | 2D meshgrid arrays `[nzx, nrx]`, from `meshgrid(rx, zx)` |
| `L.lxy` | `L` | logical `[nzx, nrx]`, true for interior points (= y-grid embedded in x-grid) |

In `meqg.m`, `rrx` and `zzx` are built temporarily via `meshgrid(G.rx,G.zx)` for Green's function computation, then removed from `G` at the end (`G = rmfield(G,{'rrx','zzx'})`). They are stored permanently in `L` by `meqmesh`.

The grid is built from P parameters (`nz`, `ri`, `ro`, `zl`, `zu`) and if provided, the limiter contour `(G.rl, G.zl)`. The half-grid offset convention is used: the x-grid is constructed so that limiter points fall between grid cells.

### 2.3 The y-grid (plasma current grid)

The y-grid is strictly the interior of the x-grid — all points except the boundary row/column:

```
ry = rx(2:end-1)    [column vector, nry = nrx-2]
zy = zx(2:end-1)    [column vector, nzy = nzx-2]
```

The plasma current density `Iy` is defined only on the y-grid. This is because the Poisson solver (`gszrmex`) takes zero flux boundary conditions (handled separately), so the source term only exists on interior points.

| Field | Description |
|---|---|
| `L.ry`, `L.zy` | y-grid 1D coordinates |
| `L.nry`, `L.nzy` | dimensions |
| `L.ny` | total y-grid points `nry * nzy` |
| `L.iry` | `1./ry` (precomputed for flux-surface integrals that need 1/R) |
| `L.rry`, `L.zzy` | meshgrid for y-grid `[nzy, nry]`, equal to `L.rrx(2:end-1,2:end-1)` |

### 2.4 The z-grid (extended post-processing grid)

An optional extended grid covering coils and vessel conductors outside the x-grid, enabled by `P.izgrid = true`. Used for computing flux and fields at locations outside the limiter.

| Field | Description |
|---|---|
| `L.G.rz`, `L.G.zz` | z-grid 1D axis vectors |
| `L.rz`, `L.zz` | same, copied into `L` by `meqmesh` |
| `L.nrz`, `L.nzz` | dimensions |
| `L.drz`, `L.dzz` | spacings |
| `L.rrz`, `L.zzz` | meshgrid `[nzz, nrz]` |
| `L.G.lzx` | boolean `[nzz, nrz]`, true where z-grid overlaps with x-grid (when `P.izxoverlap=true`) |

When `P.izxoverlap=true`, the z-grid is constructed by extending the x-grid with additional rows/columns outward, so the x-grid becomes a contiguous sub-region and flux values can be overwritten directly.

### 2.5 Grid boundary representation

The x-grid boundary (perimeter of the rectangular grid, excluding corners for the normal derivative) is stored as:

```
[F(1:nzx,1)   F(1,2)  F(nzx,2)  ...  F(1:nzx,nrx)]
```

i.e., left column first (bottom to top), then top row (excluding corners), bottom row (excluding corners), then right column. This linear ordering is used by `gszrmex`, `nfdbmex`, and the boundary condition arrays.

- `L.rb`, `L.zb` — R,Z of all boundary x-grid points (size `nb = 2*(nzx+nrx)-4`)
- `L.lxy` — logical `[nzx, nrx]`, true for interior (y-grid) points
- `L.Tbc`, `L.kxc` — precomputed matrices/indices for computing the normal flux derivative approximation `(4F(x+dx) - F(x+2dx) - 3F(x=0))/(2dx)` on boundary (Lackner's trick; see Section 5)

### 2.6 Limiter contour

| Field | Description |
|---|---|
| `G.rl`, `G.zl` | limiter polygon vertices, column vectors, closed polygon, `[nl, 1]` |
| `G.rlh`, `G.zlh` | midpoints between consecutive limiter vertices |
| `G.nl` | number of limiter points |
| `L.Oly` | boolean `[nzy, nry]`, true where y-grid points are inside the limiter polygon |
| `L.kxl`, `L.clx` | bilinear interpolation indices/coefficients for flux at limiter points |

The limiter contour must define a simple (non-self-intersecting) polygon with consistent orientation. `meqg` asserts this with `simple_polygon([G.rl,G.zl])`.

---

## 3. Electromagnetic Structure Matrices

All matrices are stored in `L.G` and are built by `meqg.m` using Green's function calls from the external `genlib` library (`greenem.m`).

### 3.1 Current source types

MEQ distinguishes three types of conductors:

**Active coils (suffix `a`)**: Poloidal field coils driven by external power supplies. `na` circuits.
- Physical windings are stored as filament arrays `(G.rw, G.zw, G.ww, G.hw)` with connection matrix `G.Twa` mapping `nw` winding filaments to `na` active circuits: `Iw = Twa * Ia`.
- `G.dima` — cell array of circuit labels (e.g. `'E_001'`, `'F_003'`, `'OH_001'`).

**Vessel filaments (suffix `v`)**: The conducting vessel wall discretised into rectangular or circular filament elements. `nv` filaments.
- `G.rv`, `G.zv` — filament centres; `G.wv`, `G.hv` — dimensions (or `G.dv` for circular).
- `G.Rv` — resistances of individual filaments.

**Generalised vessel modes (suffix `u`)**: A reduced-order description of vessel eddy currents. `nu` modes.
- `G.Tvu` is the `[nv, nu]` matrix: `Iv = Tvu * Iu`.
- Three options selected by `P.selu`:
  - `'v'` or `'vx'`: `Tvu = eye(nv)`, i.e. filament currents directly (nu = nv)
  - `'s'` or `'sx'`: segment description (user-provided `Tvs` or `Tivs`)
  - `'e'` or `'ex'`: eigenmode decomposition via `vveig(Mvv, Rv)` from genlib; `Tvu = Tve(:,1:nu)`
  - `'n'`: no vessel (`nu=0`)
- `G.Muu = Tvu' * Mvv * Tvu` — the `[nu, nu]` vessel mutual inductance in modal coordinates
- `G.Ru` — resistances for generalised vessel (diagonal of `Tvu' * diag(Rv) * Tvu`)

### 3.2 Mutual inductance matrices

All "M" matrices are computed as `M = greenem('mut', r_target, z_target, r_source, z_source)` where `greenem` computes the toroidal mutual inductance between circular current loops. Units: Henry.

**Plasma grid mutuals (most expensive; only computed when `P.gsxe=3`)**:
- `G.Mxx` — `[nx, nx]` mutual inductance between all x-grid filaments. Built as: `greenem('mut', rrx(:), zzx(:), [drrx(:), dzzx(:)])` treating each grid cell as a finite rectangular filament. This is a dense `(nrx*nzx)² ` matrix — the most expensive object in the problem.
- `G.Mby` — `[nb, ny]` boundary-to-interior subblock of `Mxx`: mutuals from boundary x-grid points to interior y-grid points. Used in Lackner's boundary condition trick.
- `G.Mbc` — `[nb, nc]` mutuals from boundary x-grid points to near-boundary interior points (excluding corners); used to build `L.Tbc`.

**Cross-mutual inductances**:
- `G.Mxa` — `[nx, na]` x-grid to active coil. Built as `Mxw * Twa` where `Mxw = greenem('mut', rrx, zzx, rw, zw)`.
- `G.Mxu` — `[nx, nu]` x-grid to vessel modes. Built as `Mxv * Tvu`.
- `G.Mza` — `[nz, na]` z-grid to active coil (optional, only when `P.izgrid`).
- `G.Mzu` — `[nz, nu]` z-grid to vessel modes (optional).
- `G.Mzy` — `[nz, ny]` z-grid to y-grid (for computing flux on z-grid from plasma current).

**Active coil self/mutual**:
- `G.Maa` — `[na, na]` computed as `Twa' * Mww * Twa`; only loaded on request.
- `G.Mau` — `[na, nu]` active-to-vessel: `Mav * Tvu`.

**Vessel self/mutual**:
- `G.Mvv` — `[nv, nv]` vessel filament mutual inductance (dense). The self-inductance diagonal uses `greenem('self', rv, zv, dv)` for circular filaments.

### 3.3 Magnetic probe response matrices

The magnetic probe at position `(rm, zm)` with orientation angle `am` (measured CCW from +R direction) measures:

```
Bm = cos(am) * Br + sin(am) * Bz
```

The response matrices are:
- `G.Bma` — `[nm, na]` probe response to active coil currents. Built as `Tmm * (cos(am) .* Bmw_r + sin(am) .* Bmw_z) * Twa`.
- `G.Bmx` — `[nm, nx]` probe response to plasma grid currents. Built analogously with `greenem('br',...)` and `greenem('bz',...)`.
- `G.Bmu` — `[nm, nu]` probe response to vessel modes. Built as `Bmv * Tvu`.
- `G.Tmm` — `[nm, nm_ideal]` mixing matrix composing effective probe measurements as linear combinations of ideal ones.

The simulated probe signal is: `Bm = G.Bma*Ia + G.Bmx(:,lxy)*Iy(:) + G.Bmu*Iu`

### 3.4 Flux loop response matrices

Flux loop at `(rf, zf)` measures total poloidal flux passing through the loop:
- `G.Mfa` — `[nf, na]` flux loop to active coils. Built as `Tff * greenem('mut', rf, zf, rw, zw) * Twa`.
- `G.Mfx` — `[nf, nx]` flux loop to plasma grid. Built as `Tff * greenem('mut', rf, zf, rrx, zzx)`.
- `G.Mfu` — `[nf, nu]` flux loop to vessel modes. Built as `Mfv * Tvu`.
- `G.Tff` — `[nf, nf_ideal]` flux loop mixing matrix.

Simulated signal: `Ff = G.Mfa*Ia + G.Mfx(:,lxy)*Iy(:) + G.Mfu*Iu`

### 3.5 Magnetic field component matrices

- `G.Brxa` — `[nx, na]` radial component Br on x-grid due to active coils. `greenem('br', rrx, zzx, rw, zw) * Twa`.
- `G.Bzxa` — `[nx, na]` vertical component Bz. `greenem('bz', rrx, zzx, rw, zw) * Twa`.
- `G.Brxu` — `[nx, nu]` Br due to vessel modes. `Brxv * Tvu`.
- `G.Bzxu` — `[nx, nu]` Bz due to vessel modes.

### 3.6 Derivative matrices

For stability calculations and some solver variants, z-derivative Green's functions are available:
- `G.dzMxa` — `[nx, na]` ∂Mxa/∂z
- `G.dzMxu` — `[nx, nu]` ∂Mxu/∂z
- `G.dzBmx` — `[nm, nx]` ∂Bmx/∂z (used in `meqreprocess` to correct for vertical displacement `dz`)
- `G.dzMfx` — `[nf, nx]` ∂Mfx/∂z
- `G.dzdzMxa`, `G.dzdzMxu` — second z-derivatives (for stability margin)

These are all computed via `greenem('dz1m',...)` etc.

### 3.7 Wall gaps

- `G.rW`, `G.zW`, `G.aW`, `G.oW` — origin (R,Z), length, orientation of wall gap segments (`nW` gaps)
- `L.erW`, `L.ezW` — fractional grid indices of gap origins
- `L.crW`, `L.czW` — direction cosines of gaps

### 3.8 Dipole matrix

- `G.Dda` — `[nd, na]` matrix representing differences between adjacent E/F coil groups; `nd = na_PF - 1`.

---

## 4. Plasma Domain Detection

This subsection explains how the plasma boundary is found from the flux map `Fx(Z,R)`.

### 4.1 Finding critical points (O-point, X-points): `asxymex`

The C MEX function `asxymex` scans the 2D flux map `Fx` to find all **extrema** (O-points, i.e. magnetic axes) and **saddle points** (X-points). The MATLAB fallback `.m` file just throws an error — the real code is in C for speed.

Algorithm (inferred from `asxyJac.m` and the help text):
1. At every interior grid cell centre `(iz, ir)`, fit a local 6-point quadratic:
   ```
   f(x,y) ≈ f0 + a*x + b*y + (c/2)*x² + (d/2)*y² + e*x*y
   ```
   where `x = (r - r_cell)/dr`, `y = (z - z_cell)/dz`, and the 6 stencil points are the cell centre, left, down, lower-left, right, up neighbours.
2. The critical point within the cell is at `xe = (e*b - a*d)/h`, `ye = (e*a - b*c)/h` where `h = c*d - e²` (the Hessian determinant). If `|xe|, |ye| <= 1`, the critical point is inside the cell.
3. An extremum is detected if `h > 0` (Hessian positive or negative definite). A saddle is detected if `h < 0`.
4. The sub-grid location is refined further by the Gauss-Newton iteration in `asxycs.m` (when `P.icsint=true`), which uses bicubic spline interpolation of `Fx` to converge the critical point to machine precision.

Signature:
```matlab
[zA,rA,FA,dz2FA,dr2FA,drzFA,ixA, zX,rX,FX,dz2FX,dr2FX,drzFX,ixX, stat] = ...
  asxymex(Fx, zx, rx, DC, dzx, drx, idzx, idrx, L_mask, DIMW)
```
- `zA,rA,FA` — Z, R, flux at O-points (axes)
- `zX,rX,FX` — Z, R, flux at X-points
- `dr2FA,dz2FA,drzFA` — second derivatives of flux at axes (Hessian components)
- `ixA`, `ixX` — flattened x-grid indices of cells containing each critical point (used in `asxyJac`)
- `DC` — distance threshold to merge nearby extremum/saddle pairs
- `L_mask` (`L.Oasx`) — boolean search domain; the critical point search is restricted to cells where `L_mask` is true plus their immediate neighbours
- `DIMW` — maximum number of critical points to search (workspace size limit)

The search domain `L.Oasx` is built in `meqc.m` as the y-grid limiter interior (`L.Oly`) expanded by one cell on all four sides.

When `P.icsint=true`, `asxycs` refines using `cs2devalmex` (bicubic spline evaluation) via Gauss-Newton: move each critical point by `-(H^{-1} * grad F)`, clamped to one cell maximum, until convergence.

### 4.2 O-point (magnetic axis) selection

After `asxymex`, the axes are filtered in `meqpdom`:
- Axes with wrong Hessian sign are removed: for `Ip > 0`, valid axes must have `dr2FA + dz2FA < 0` (i.e. flux is a maximum at the axis).
- Axes are sorted by descending Z (convention: top axis first for doublets).

### 4.3 X-point polygon and domain membership (`bavxmex`)

The X-point defines a separatrix. Each X-point has an associated "outward normal" direction `(vrX, vzX)` computed from the Hessian eigenvectors:
```matlab
H1 = sIp*(dr2FX - dz2FX) * 0.5;
H0 = sqrt(H1^2 + drzFX^2);
vrX = sqrt(H0 + H1);
vzX = sqrt(H0 - H1) * sign(drzFX) * sIp;
```
`bavxmex` constructs the "X-point polygon" — a half-space defined by the direction pointing from the X-point toward the axis — and returns a logical array of which grid points (or limiter/other X-points) lie inside that polygon. This is used to restrict which X-points and limiter points can contribute to a given domain's boundary.

### 4.4 Finding the LCFS (`fbndmex`)

`fbndmex` selects the LCFS flux by comparing:
- Flux values at filtered limiter points `Fl` (interpolated from the x-grid to `(rl, zl)`)
- Flux values at candidate X-points `FX`

For `Ip > 0` (standard direction), `FB = max(valid Fl and FX candidates)`. If an X-point provides the maximum, `lX = true` (diverted plasma). Otherwise `lB` comes from a limiter point.

The marker value `FN = -Inf * sign(Ip)` flags invalid/excluded entries in `Fl` and `FX` arrays.

### 4.5 Plasma domain map (`Opy`)

Once `FB` is known for each domain, `Opy` is computed on the y-grid:
- `pdommex` sets `Opy(iz,ir) = domain_index` where `sign(Fx - FA) * sign(Fx - FB) <= 0` and the point is inside the X-point polygon.
- After masking with `L.Oly` (inside limiter), `Opy` is an `int8` array on the y-grid:
  - `0` = outside all plasma domains
  - `1..nD` = domain index

For simple single-null diverted or limited plasmas, `nD=1` and `Opy` is simply 1 inside the plasma and 0 outside. For doublets, `nD=3` (two inner domains + one outer mantle domain).

### 4.6 Limiter vs. diverted distinction

| `lX` | Condition | Plasma type |
|---|---|---|
| `false` | LCFS defined by a limiter point | **Limited plasma** |
| `true` | LCFS defined by an X-point | **Diverted plasma** |

The flux ratio `FR` characterises proximity to the divertor transition:
- `FR = 1`: X-point is exactly on LCFS (diverted)
- `FR > 1`: limited plasma, X-point exists but doesn't define LCFS; `FR` measures how much further the X-point flux is from the axis relative to the limiter flux range
- `FR < 1` (mantle domain in doublet): X-point is inside the mantle separatrix

---

## 5. Grad-Shafranov Poisson Solver

### 5.1 Equation being solved

The Grad-Shafranov equation is:

```
Δ* ψ = -μ₀ R j_φ
```

where `Δ* = R ∂/∂R (1/R ∂/∂R) + ∂²/∂Z²` is the Grad-Shafranov operator. In the finite difference form on the rectangular (R,Z) grid, this becomes:

```
(dz/dr) * [1/(r_i + dr/2) * (ψ_{i+1,j} - ψ_{i,j}) - 1/(r_i - dr/2) * (ψ_{i,j} - ψ_{i-1,j})]
  + (dr/dz) * (ψ_{i,j+1} - 2*ψ_{i,j} + ψ_{i,j-1}) = -cx_i * I_{i,j}
```

where `cx = π μ₀ (dz/dr) r` (see `gszrc`), and `I_{i,j} = Iy(iz, ir)` is the plasma current in the `(iz, ir)` cell.

This is solved by `gszrmex(Fb, Iy, cx, cq, cr, cs, ci, co, dz)`.

### 5.2 Finite difference discretisation and the FFT-ADI solver (`gszrc`, `gszrmex`)

The solver uses an **Alternating Direction Implicit (ADI) / Fourier-decomposition** approach, following [Hofmann, Comput. Phys. Commun. 48(2) 1988]:

**Step 1 — Fourier decomposition in Z**: Since the Z-boundary conditions are `ψ=Fb` (known external flux), the system decouples into `nry` independent tridiagonal systems after a discrete sine transform (DST) in Z. The constraint `(nzx-1) = 2^k` is required for the FFT. The `two` array in `gszrc` precomputes the DST eigenvalues.

**Step 2 — Tridiagonal solve in R**: For each Z Fourier mode, a tridiagonal system in R is solved. The precomputed arrays are:
- `cx` — `[nry, 1]` RHS scaling: `π μ₀ (dz/dr) r`  
- `cq` — `[nry, nzy]` precomputed LDL pivot values (inverse diagonal)  
- `cr` — `[nry, nzy]` off-diagonal coupling coefficients  
- `cs` — `[nry, nzy]` = `cw * cq` (lower off-diagonal)  
- `ci` — scalar, inner boundary coupling (left wall effect)  
- `co` — scalar, outer boundary coupling (right wall effect)

The `gszrmex` MEX function performs the solve in C. The `dz` argument shifts the solution by `dz * dzx` (used for vertical displacement in LIUQE stabilisation function).

### 5.3 Boundary conditions

The Poisson solver requires the flux `Fb` on the rectangular x-grid boundary. This is assembled from two contributions:

**External contribution** (coils + vessel):
```
Fb_ext = G.Mbe * [Ia; Iu]
```
where `G.Mbe = [G.Mxa(~L.lxy,:), G.Mxu(~L.lxy,:)]` — the subrows of the mutual inductance matrices corresponding to boundary grid points.

**Plasma contribution** (`meqfbp`):
Two options controlled by `P.ilackner`:
- `ilackner = 1` (Lackner's trick, standard): 
  ```
  Fb_plasma = L.Tbc * nfdbmex(gszrmex(zeros(2*(nzx+nrx-2),1), Iy, cx, cq, cr, cs, ci, co, 0))
  ```
  This solves the Poisson equation with `Fb=0` (internal problem), then uses the normal derivative on the boundary to approximate the external flux contribution from the plasma. `nfdbmex` computes `4*F(x+dx) - F(x+2dx)` on the boundary (excluding corners). `L.Tbc` converts this to the actual boundary flux using the self-inductance of the boundary segments.
  
  The self-inductance of each boundary segment is computed in `meqmesh` using the formula:
  ```
  Lcc = mu0 * rc * ((1 + ...) * log(8*rc/(drc+dzc)) + ...)
  ```
  
- `ilackner = 0` (explicit Green's function):
  ```
  Fb_plasma = G.Mby * Iy(L.lxy)
  ```
  where `G.Mby = G.Mxx(boundary, interior)` — the off-diagonal block of `Mxx`. Requires `Mxx` to be loaded (very expensive).

- `ilackner = 2`: Pre-computes `L.Mby` at setup time by solving the internal Poisson problem for each unit y-grid basis vector.

### 5.4 Plasma current source term

The source `Iy` is the toroidal current **per cell** (Amperes), related to current density by `j_φ = Iy / (dr * dz * R)`. For the Poisson equation, the RHS factor `rhsf` is precomputed:

```matlab
L.rhsf = (2.0 * dzx * idrx * pi * mu0) .* rry(:)
```

The relationship `gszrmex` uses internally is `Δ*ψ_cell = -cx_i * Iy_cell` where `cx = π μ₀ (dz/dr) r`.

The plasma current is parameterised via basis functions: `j_φ = R p'(ψ) + (1/R) T T'(ψ)` where `p'` and `TT'` are expanded in `ng` basis functions with coefficients `ag`. The base function code (`bf3pmex`, `bfabmex`, etc.) evaluates these profiles. The current density integrated over the y-grid cells gives `Iy`.

### 5.5 External currents inside the computational grid (`P.gsxe`)

When active coils or vessel filaments lie inside the x-grid boundary (unusual but possible), three options exist:

- `P.gsxe = 1`: Assign external filament current to the nearest y-grid cell using a sparse map matrix `L.Tye`.
- `P.gsxe = 2`: Use the mutual inductance to compute the effective y-grid source: `Iyie = -(dlst * Mxie) ./ rhsf` where `dlst` is the discrete Laplace-Shafranov operator matrix. [GUESS: This converts an external filament current into the equivalent Iy that would produce the same flux via the Poisson operator.]
- `P.gsxe = 3`: Use full Green's function approach: store `L.Mxe` and `L.Mxy` for direct flux computation.

### 5.6 Block solver (`bslvmex`)

Used to solve the coupled system `[ag; Ie]` (basis function coefficients and external currents):
```
[a; b] = [a1; b1] - [A, B] * [a; b]
```
solved iteratively by alternating:
```
a = a1 - A*b
b = b1 - B*a
```
starting from `b = b0` and running `NI` iterations. This is a fixed-point / Gauss-Seidel-style iteration on the two-block system. It is used in LIUQE to decouple the plasma-profile fit from the coil-current fit.

### 5.7 The DLST operator

`L.dlst = dlstc(nzx, nrx, drx, dzx, rry)` — the **Discrete Laplace-Shafranov Transform** operator as a matrix, mapping a y-grid flux distribution to its Δ* value. Used in `P.gsxe=2,3` to invert the Poisson problem analytically.

### 5.8 Finite element representations (`fetr`, `febi`)

For parameterising the plasma current profile using finite elements (used in FBT):

**Triangular FE (`fetr`)**: Constructs a set of piecewise-linear tent functions on a triangular mesh inscribed in the limiter contour. Returns `Tyh` — the FE basis function values evaluated at all y-grid points. The node positions `(rh, zh)` have one inner node at mid-width at top/bottom levels and two inner nodes at intermediate levels (schematic in the help text). Used when `P.bfct = @bffbt` or similar.

**Bilinear FE (`febi`)**: Bilinear (Q1) elements on the rectangular x-grid. [GUESS based on naming and usage: returns bilinear basis function values at grid points.]

---

## 6. Flux Surface Integrals

### 6.1 Flux surface parameterisation

After the plasma domain is identified, flux surfaces are parameterised in `(ρ, θ)` coordinates:
- `ρ` (rho) — normalised minor radius-like coordinate. `L.pq` = `[ρ₁, ..., ρ_npq, ρ_LCFS]` — the rho grid, `npq+1` points.
  - At LCFS: `ρ = pq(end)` ≈ `sqrt(2) - sqrt(1 - pqu²)` (upper end of `pq` grid).
  - On axis: `ρ = 0`.
  - The normalised flux `ψ_N = ρ²` (where `ψ_N = 0` on axis, 1 on LCFS), stored as `L.pQ = [0, L.pq]`.
- `θ` — poloidal angle. `L.oq = [0, 2π/noq, ..., 2π*(noq-1)/noq]`, `noq` equispaced points.
  - `crq = -cos(oq)`, `czq = sin(oq)` — direction cosines in R and Z from axis.

The flux surface at `(ρ, θ)` is traced by `rtcimex`, which iteratively solves for the distance `a(ρ, θ)` from the magnetic axis along direction `(crq, czq)` such that `Fx(rA + a*crq, zA + a*czq) = FQ(ρ)`. The flux surface position is then `rq = rA + aq*crq`, `zq = zA + aq*czq`.

### 6.2 Flux surface geometrical integrals (`fsgimex`)

`fsgimex` computes all major geometrical flux-surface integrals from the traced surfaces. Inputs: `a1q, a2q` — the Jacobian-related quantities `A1 = da²/dρ/ρ` and `A2`, computed from the `aq` array via the spline derivative matrices `L.M1q` (ρ-derivative) and `L.M2q` (θ-derivative).

Key outputs (all as functions of `(ρ, domain)` on the `pQ` grid):

| Output | Symbol | Formula |
|---|---|---|
| `Q0Q` | `<1/R>` | Flux surface average of 1/R |
| `Q1Q` | `-dψ/dV` | Inverse of volume gradient |
| `Q2Q` | `<1/R²>` | Flux surface average of 1/R² |
| `Q3Q` | `<|∇ψ|²/R²>` | Related to ∇ψ squared average |
| `Q4Q` | `<|∇ψ|²>` | |
| `iqQ` | `1/q` | Inverse safety factor |
| `ItQ` | `Ip(ρ)` | Enclosed poloidal current `(1/μ₀) ∮ Bp dlp` |
| `LpQ` | `∮ dlp` | Poloidal circumference |
| `rbQ` | `∮R dlp / ∮dlp` | Poloidal barycentre in R |
| `Q5Q` | `<|∇ψ|/2π>` | |
| `SlQ` | `2π ∮ R dlp` | Toroidal lateral area |

At the axis, limits are applied: `Q0Q = 1/rA`, `iqQ = 0` (regular), etc. At the X-point LCFS, special limiting values are used for `Q0Q, Q2Q, iqQ, ItQ` as detailed in the `fsgimex` help text.

The safety factor `q = 1/iqQ`.

### 6.3 Volume integrals (`vizrmex`)

`vizrmex(Fx, Opy, ry, iry, rBt, drx, dzx)` computes three volume integrals over the plasma:

- `Wp` — poloidal field energy: `∫_V Bₚ²/(2μ₀) dV`. Computed using finite differences of `Fx` to get `Bₚ = |∇ψ| / R`.
- `Ft0` — vacuum toroidal flux: `∫_S R Bφ₀ dR dZ` (using `rBt/R`).
- `Vp` — plasma volume: `∫_V dV = 2π ∫∫_plasma R dR dZ`.

The Jacobian matrix `vizrJac` computes `dWp/dFx` for Newton-based solvers.

### 6.4 Cumulative integrals on ρ grid (`meqintQ`)

`meqintQ` computes all integrals as cumulative functions of ρ (from axis to LCFS):
- `IpQ(iQ, iD)` — enclosed plasma current within ρ = `pQ(iQ)`, domain `iD`
- `WkQ(iQ, iD)` — kinetic energy within ρ
- `WpQ(iQ, iD)` — poloidal field energy within ρ
- `VpQ(iQ, iD)` — volume within ρ
- `FtQ(iQ, iD)` — toroidal plasma flux within ρ
- `Ft0Q(iQ, iD)` — vacuum toroidal flux within ρ

The total quantities are `IpQ(end,:)`, `WkQ(end,:)`, etc.

### 6.5 Global integral quantities (`meqint`)

From the profile integrals:
- `Ip` — total plasma current (A)
- `Wk` — kinetic energy (J): `1.5 * x * fPg` where `x = ITpDg .* ag`
- `Wt` — plasma toroidal field energy: `x * fTg`
- `Wp` — poloidal field energy (from `vizrmex`)
- `Vp` — plasma volume (from `vizrmex`)
- `WN = 1e-7 * π * r₀ * Ip²` — normalisation for beta
- `bp = (2/3) * Wk / WN` — poloidal beta
- `bt = (2/3) * Wk / Wt0` — toroidal beta
- `li = Wp / WN` — normalised internal inductance
- `mu = Wt / WN` — diamagnetism

The `bpli2 = bp + li/2` is the Shafranov parameter.

---

## 7. MEX Architecture

### 7.1 Why C MEX routines

The innermost computational loops are implemented in C MEX files for two reasons:
1. **Performance**: Operations like the Poisson solve (`gszrmex`), critical point scanning (`asxymex`), and flux surface contouring (`rtcimex`) must run in milliseconds for real-time operation (LIUQE runs on TCV's Real-Time Controller).
2. **Non-vectorisable control flow**: `asxymex` searches for extrema/saddles with early-exit logic. `rtcimex` traces flux surfaces with Newton iteration per point. These do not vectorise naturally in MATLAB.

### 7.2 The MATLAB fallback pattern

Every MEX file follows the same pattern:
```matlab
function varargout = FOOmex(varargin)
  warning('MEQ:MEX "%s" not compiled', mfilename)
  [varargout{1:nargout}] = FOOmexm(varargin{:});
end
```
The `FOOmexm` suffix denotes a pure-MATLAB fallback implementation. This means every MEX has a reference MATLAB implementation. The fallbacks are slower (often 10-100×) but are maintained for portability and testing.

### 7.3 MEX function catalogue

| MEX | Purpose |
|---|---|
| `gszrmex` | Poisson solver via FFT-ADI (Hofmann algorithm) |
| `asxymex` | Find all extrema and saddle points of 2D field (6-point quadratic fit) |
| `fbndmex` | Select LCFS from limiter and X-point flux candidates |
| `pdommex` | Flag plasma domain given flux map, LCFS flux, sign of Ip, and X-point polygon |
| `rtcimex` | Newton iteration to trace flux surfaces (isoflux contouring) |
| `fsgimex` | Geometrical integrals over traced flux surfaces |
| `nfdbmex` | Normal flux derivative approximation on rectangular boundary |
| `bavxmex` | Classify grid points inside/outside an X-point polygon |
| `bboxmex` | Compute bounding box of a 2D logical domain |
| `bslvmex` | Alternating block solver for `[ag; Ie]` system |
| `respmex` | Efficient `[A*B, C*D]` horizontal concatenation |
| `iatamex` | LDL decomposition inversion of `A1'*A1 + A2'*A2` (3 or 4 row A) |
| `uatamex` | Upper triangle of `A'*A` |
| `bintmex` | 4-point bilinear interpolation at arbitrary points |
| `cizrmex` | Cross-section integrals of arbitrary y-grid fields as functions of ρ |
| `vizrmex` | Volume integrals (Wp, Ft0, Vp) |

### 7.4 The block solver algorithm (`bslvmex`)

The MATLAB fallback code reveals the exact algorithm:
```matlab
b = b0;           % initial guess for b (vessel currents or Iu)
for k = 1:ni
  a = a1 - A*b;  % solve for ag given current Ie/Iu
  b = b1 - B*a;  % solve for Ie/Iu given current ag
end
```
where `a = ag` (basis function coefficients), `b = Ie` (external currents). This is a Gauss-Seidel block iteration. In the RT version `ni=1`, providing a single-sweep approximation.

### 7.5 LDL factorisation (`iatamex`)

`iatamex(A1, A2)` inverts `A1'*A1 + A2'*A2` where `A1` and `A2` have 3 or 4 rows. The restriction to small matrices means this is a special-purpose fast path: with only 3-4 rows, the normal equations matrix has dimension `ng × ng` where `ng` is the number of basis functions (typically 4-10). The LDL decomposition is unrolled at compile time in the C code.

---

## 8. TCV-Specific Geometry (`meqgtcv`, `mds2meq`)

### 8.1 TCV coil configuration

TCV (Tokamak à Configuration Variable) has:
- **E coils**: 8 individual shaping coils (`E_001` through `E_008`)
- **F coils**: 8 outer shaping coils (`F_001` through `F_008`)
- **G coil**: 1 gas-puff / glow coil (`G_001`)
- **OH coils**: 2 ohmic heating coils (`OH_001`, `OH_002`)
- **TOR coil**: Optional toroidal field coil supplement

The `sela` argument to `meqgtcv` selects which coil groups to include (default: all `{'E','F','G','OH'}`). Individual coils can also be specified directly.

`G.Twa` maps individual conductor windings to these circuit definitions. Its shape is `[nw, na]`.

The dipole matrix `G.Dda` is constructed to represent the difference between adjacent E and F coil currents (used for vertical stability control).

### 8.2 Vessel description

The TCV vessel is discretised into `nv` rectangular filament elements, loaded from MDSplus `STATIC("R_V")`, `STATIC("Z_V")`, etc.

Available vessel descriptions (`selu`):
- `'vx'` (default): use all `nv` filaments directly
- `'sx'`: segment description — filaments grouped into `ns` segments via `Tivs` (inverse: `Is = Tivs * Iv`)
- `'ex'`: eigenmode decomposition — `nu` dominant eddy current modes from `vveig(Mvv, Rv)`

The Green's function matrices `Mfu` and `Bmu` are loaded from MDSplus `STATIC` tree to avoid recomputation.

### 8.3 Magnetics geometry

From MDSplus:
- **Magnetic probes**: `G.rm`, `G.zm`, `G.am` from `STATIC("R_M")`, `STATIC("Z_M")`, `STATIC("ANG_M")`. Typically 38 probes (`selm = 1:38`).
- **Flux loops**: `G.rf`, `G.zf` from `STATIC("R_F")`, `STATIC("Z_F")`. Typically 38 loops (`self = 1:38`).

Pre-computed Green's function matrices `G.Mfa`, `G.Bma`, `G.Mfx`, `G.Bmx` are loaded directly from MDSplus rather than recomputed:
```matlab
G.Mfa = mdsdata('STATIC("MUT_"//$1//"_A")[$2,$3]', self, G.dimf, G.dima)
G.Bma = mdsdata('STATIC("BT_M_A")[$1,$2]', G.dimm, G.dima)
```
Note: `G.Bma(:, TOR_coil)` is zeroed out because TCV magnetic probe signal processing removes the TOR contribution.

### 8.4 Limiter contour

The limiter is stored as a piecewise polynomial (pp) in MDSplus:
```matlab
G.ppl = mkpp(mdsdata('STATIC("PP_T:BREAKS")'), mdsdata('STATIC("PP_T:COEFS")'), 2)
```
`128` equispaced limiter points are sampled from this pp representation.

### 8.5 Coil protection limits

From the `BASE` tree:
- `G.Iamax`, `G.Iamin` — per-coil current limits
- `G.Vamax`, `G.Vamin` — per-coil voltage limits
- `G.Ceq`, `G.Seq`, `G.Ilow`, `G.Ihig` — coil protection linear inequality constraints (`Ceq * Ia` in `[Ilow, Ihig]`)

### 8.6 MDSplus data flow (`mds2meq`)

The `mds2meq(shot, source, server)` function reads stored equilibrium results:

**For `source = 'LIUQE.M'`** (stored LIUQE results):
1. Opens `TCV_SHOT` tree
2. Calls `liuqe(shot, [])` to build `L` structure
3. Reads time array from `tcv_eq('TIME_PSI', source)`
4. Reads flux `Fx` (from `PSI` node), converted: `LY.Fx = permute(FxB + FB, [2,1,3])` — note the transpose to go from `(nrx, nzx)` MDSplus storage to MEQ `(nzx, nrx)` convention
5. Reads plasma current: `LY.Iy = permute(Jx(2:end-1,2:end-1,:), [2,1,3]) * L.dsx` — converts `j_φ` (A/m²) to per-cell current by multiplying by cell area
6. Reads all equilibrium scalars, profiles, coil currents

**For `source = 'MAG'`** (raw magnetics only):
1. Loads coil voltages from power supply measurements
2. Builds L structure without equilibrium calculation

### 8.7 x-grid naming convention in MDSplus

TCV stores pre-computed Green's function matrices in MDSplus for three grid variants selected by `selx`:
- `''` — standard grid
- `'X'` — X-point mesh (finer near X-point region)
- `'XF'` — extended X-point mesh

The MDSplus node names follow patterns like `STATIC("R_X:AXIS")`, `STATIC("MUT_X_A")[*,$2]`, etc.

---

## 9. Issues for Python Rewrite

### 9.1 Matrix building (Mxx, Mxa, etc.)

**Current approach**: All matrices are dense `float64`. The biggest is `Mxx` which is `(nzx*nrx) × (nzx*nrx)`. For TCV with `nzx=65, nrx=33`, that is `2145 × 2145` = 4.6M doubles ≈ 37 MB. Manageable as NumPy dense arrays.

**Python equivalent**: Use `numpy.ndarray` for all Green's function matrices. The `greenem` function from genlib must be reimplemented; it computes the Neumann function (mutual inductance between two coaxial loops) using complete elliptic integrals.

**Sparse consideration**: `Mxx` could [GUESS] potentially be stored sparse since mutual inductance between distant cells is small, but the code currently uses it dense. `G.Mby` and sub-blocks can also remain dense given the sizes involved.

### 9.2 MEX routines in Python

| MEX | Python equivalent |
|---|---|
| `gszrmex` | Custom FFT-based solver or `scipy.sparse.linalg.spsolve` on the Δ* operator assembled as sparse matrix |
| `asxymex` | Custom scan with `scipy.ndimage` for local extrema + Newton refinement |
| `fbndmex` | Simple `numpy` argmax/argmin operation |
| `pdommex` | Boolean numpy array operations |
| `rtcimex` | Newton iteration in `numpy`; can be vectorised over θ angles |
| `fsgimex` | Numerical integration (trapezoid, Simpson) along θ dimension |
| `nfdbmex` | Simple finite difference stencil extraction in numpy |
| `bavxmex` | Custom half-plane test; or `shapely.geometry.Polygon.contains` |
| `bboxmex` | `numpy.where` followed by `min/max` |
| `bslvmex` | Pure Python loop; already has MATLAB fallback |
| `iatamex` | `numpy.linalg.solve` or LDL via `scipy.linalg.ldl` |
| `uatamex` | `numpy.triu(A.T @ A)` |
| `bintmex` | `scipy.interpolate.RegularGridInterpolator` or manual 4-point bilinear |
| `vizrmex` | Finite difference gradient + masked summation in numpy |
| `cizrmex` | Masked cumulative summation in numpy |

The Poisson solver deserves special attention: `gszrmex` implements the Hofmann FFT-ADI algorithm which is algorithmically superior to a direct sparse solve. For a Python replacement, assemble the Δ* operator as a `scipy.sparse.csr_matrix` and use `scipy.sparse.linalg.spsolve` or, to match the original algorithm more closely, implement DST-based decoupling using `scipy.fft.dst`.

### 9.3 Array indexing conventions

**Critical difference**: MATLAB is 1-based column-major; Python/NumPy is 0-based row-major.

- All stored indices (`L.kxl`, `L.kxW`, `ixA`, `ixX`) are 1-based MATLAB linear indices. Subtract 1 when converting.
- The flux map `Fx` is `[nzx, nrx]` in MATLAB (Z is the "row" dimension, R is the "column"). In NumPy, this is the same `(nzx, nrx)` shape but the memory layout is different (C-order vs Fortran-order). The `Br` and `Bz` formulas use `Fx[iz+1,ir] - Fx[iz-1,ir]` for finite differences — these translate directly.
- `meshgrid(rx, zx)` in MATLAB produces `rrx[iz, ir] = rx[ir]` and `zzx[iz, ir] = zx[iz]` — this matches NumPy's `np.meshgrid(rx, zx)` (note the 'ij' vs 'xy' indexing: MATLAB `meshgrid(rx,zx)` is equivalent to `np.meshgrid(rx, zx)` with default `indexing='xy'`).

### 9.4 The struct-of-arrays pattern

MEQ packs time series as struct-of-arrays where the last dimension is always time. For example `LY.Ia` is `[na, nt]`.

**Python recommendation**: Use a Python dataclass with NumPy arrays. The `meqxk(LX, k)` extractor becomes:
```python
@dataclasses.dataclass
class MEQState:
    t: np.ndarray        # shape (nt,)
    Fx: np.ndarray       # shape (nzx, nrx, nt)
    Iy: np.ndarray       # shape (nzy, nry, nt)
    Ia: np.ndarray       # shape (na, nt)
    ...
    
def meqxk(LY: MEQState, k: Union[int, np.ndarray]) -> MEQState:
    # Slice last dimension of each array field
```

Alternatively, a dict of arrays works with less boilerplate. The `meqlpack` / `meqxk` pattern is essentially a time-series container; Pandas DataFrame or xarray Dataset could handle the 1D fields, but the 2D fields (`Fx`, `Iy`, `Opy`) would need special treatment.

### 9.5 The L structure

The `L` structure mixes two conceptually different things: precomputed scalars/matrices (suitable for a Python dataclass or frozen namespace) and mutable state (not appropriate). The cleanest split:

- `L.G` → a frozen `Geometry` dataclass (build once from `meqg` equivalent)
- `L.P` → a frozen `Parameters` dataclass (configuration)
- Remaining `L` fields → a `Setup` dataclass (built by `meqc` equivalent)

### 9.6 The `nzx` power-of-2 constraint

`meqmesh` asserts `(nzx - 1)` must be a power of 2. This is required by the FFT-based Poisson solver. Any Python implementation of `gszrmex` that uses DST must respect this.

### 9.7 The `Fx` storage sign convention

MEQ stores the poloidal flux `Fx` as an absolute value in Webers. The normalised flux `ψ_N = (Fx - FA) / (FB - FA)` ranges from 0 (at axis) to 1 (at LCFS). Note that for `Ip > 0`, `FA > FB` typically (flux is a maximum at the axis), so `FB - FA < 0`. The sign of `Ip` determines the direction of the flux gradient: `Bz = (1/(2πR)) * ∂Fx/∂R` should be positive at the outboard midplane for standard plasma.

### 9.8 Multi-domain (doublet) capability

For doublet plasmas (`P.idoublet=true`), `L.nD=3` domains are tracked. Each of `F0`, `F1`, `FA`, `FB`, `Opy`, `IpD`, etc. has a leading domain dimension. A Python implementation should support `nD >= 1` from the outset.

### 9.9 The `Txy` / `Txb` embedding matrices

`L.Txy` is a `(nx, ny)` sparse identity that embeds y-grid vectors into x-grid (zero on boundary). `L.Txb` is the complementary `(nx, nb)` sparse matrix for boundary points. These are used in Newton-type solvers to avoid indexing with `L.lxy` repeatedly. In Python, `scipy.sparse.csr_matrix` or simple boolean indexing with `lxy` would suffice for most uses.

---

## Appendix: Key Constants and File Locations

- `mu0 = 4π × 10⁻⁷` H/m — used as `mu0` throughout (MATLAB global or constant)
- `L.P.r0` — reference major radius for normalisation (TCV: ≈ 0.88 m)
- `L.P.b0` — reference toroidal field (TCV: ≈ 1.4 T)
- `L.Ip0` — nominal plasma current for scaling: `1e6 * al² / Rl * P.b0` (from limiter geometry, q95≈4)
- `L.Fx0` — nominal flux swing for scaling (from self-inductance at limiter major radius)

**File locations** (all in `/home/grandin/repos/meq/`):

| File | Role |
|---|---|
| `meqg.m` | Builds G: all Green's function matrices |
| `meqmesh.m` | Builds mesh metadata (rx,zx,ry,zy and derived) |
| `meqc.m` | Builds full L structure (calls meqmesh, gszrc) |
| `meqpdom.m` | Plasma domain detection (calls asxymex, fbndmex, bavxmex) |
| `meqpost.m` | Post-processing: integrals, profiles, simulated measurements |
| `meqint.m` | Volume integral scalars from profile integrals |
| `meqfbp.m` | Plasma boundary condition for Poisson solver |
| `gszrc.m` | Precompute Poisson solver coefficients |
| `gszrmex.m` | MEX wrapper for FFT-ADI Poisson solver |
| `fetr.m` | Triangular FE basis functions |
| `meqgtcv.m` | TCV geometry from MDSplus |
| `mds2meq.m` | Load stored TCV equilibrium into MEQ structures |
| `meqxk.m` | Extract time slices from LX/LY |
| `meqxconvert.m` | Convert LX between different L configurations |
| `meqreprocess.m` | Recompute post-processing from stored Fx, Iy |
| `tests/meqpdom_test.m` | Tests for all plasma domain configurations |
| `tests/meq_test.m` | Base test class with circular equilibrium fixture |

**Reference**: [J-M. Moret et al. Fus.Eng.Des. 91, 2015] — the primary algorithmic reference for MEQ.
