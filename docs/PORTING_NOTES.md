# TITAN 1.2 → Julia: porting notes

TITAN 1.2 (Gehmeyr & Mihalas 1994, Physica D 77, 320; release of Oct 1995, UIUC
permissive licence; the notice must be retained in derived works).

Evidence tags used below:
- **[run]**: verified by building and running the code here.
- **[read]**: verified by reading the source.
- **[inferred]**: my interpretation, not directly checked. Treat these with care.

Layout of `/home/claude/TitanPort/`:

| path | contents |
|---|---|
| `fortran/` | copy of `version_1.2/inner/*` plus top-level README/UPDATE/BENCHMARKS, patched minimally (see §2) and built in place (`titan.x`). Originals of the patched files are kept as `*.orig`. |
| `reference_runs/<P>/` | raw outputs (`*.out`, `*.dmp`, `*.nrg`, `*.sty`, `*.tty`), `time.txt`, `<P>_final.csv` (full double precision), `dumps/<P>_dNNN.csv` (every archived model) |
| `reference_runs/summary.csv` | the benchmark table from §3 |
| `scripts/` | `run_bench.sh` (one problem, 600 s cap), `run_all.sh` (reproduces everything), `read_dump.py` (reads the direct-access `.dmp` files) |
| `logs/` | gfortran build logs |

---

## 1. Architecture map

### 1.1 Control flow [read]
```
titan → opnfil (units; reads EOS tables if requested) → start → start{01..10}
  loop jstep: step → archive → writout → writsty ; stop when timen ≥ tmax
step:  peq (save old state; extrapolate r, m to t+dt with old drdt, dmdt;
            check monotonicity, retry up to `next` times; EOS/opac;
            Eddington factors; time-centred quantities)
       Newton loop (iter ≤ niter=15):
            matgen → matslv → maxdel → update(eos, opac, [eddfac if ltran=0])
            converged when dmax ≤ conv (1e-5)
       on failure: dt *= 0.05, reset, retry (ntry=6), then stop 'step01'
       on success: one extra "boost" Newton iteration (lboos=1), then timstp, totnrg
```
`peq`, `update`, `maxdel`, `timstp` and `matgen` each dispatch on `(lhydr, lrad, ltur)`
to `*h` (pure hydro), `*r` (radiation in a static medium), or `*rh` (radiation
hydrodynamics). Every `ltur ≥ 1` branch is `stop`.

### 1.2 Unknowns per zone and staggering [read]
Fortran index `k`. Zone `k` spans `r(k) … r(k+1)`. The real zones are `ngrs … ngre`.
There are 2 phantom zones on each side (`ngrs-2, ngrs-1, ngre+1, ngre+2`), and
`mgr = 330` is the maximum.

| eq. index | unknown | location | Newton correction |
|---|---|---|---|
| `ir` grid | r | interface | δ ln r |
| `im` mass | x_me = exterior mass / 4π (v1.2 change from interior mass) | interface | δ ln x_me |
| `id` continuity | ρ | centre | δ ln ρ |
| `iu` gas momentum | u | interface | δu / u_nom, u_nom = \|u\| + 1e-30 |
| `it` gas energy | T (equation is for e_gr = e_gas + E_r/ρ in the rh case) | centre | δ ln T |
| `ie` radiation energy | E_r | centre | δ ln E_r |
| `if` radiation momentum | F_r (comoving) | interface | δF / F_nom, F_nom = \|F\| + 1e-30 |
| (8th slot, `ic`/turbulence) | unused | | |

- rh: `neqn = 7`.
- Hydro: `neqn = 5` (r, m, ρ, u, T).
- Static radiation: `neqn = 4` (r, T, E, F).
- Masses and volumes carry no 4π: `dvol = (r(k+1)³ − r(k)³)/3`, interior mass `xm = xmtot − xme`.

### 1.3 Equations (conservative finite volume on a moving mesh) [read]
- **Mass**: the cumulative-mass relation `xme(k-1) − xme(k) = ρ ΔV`, plus continuity in
  flux form with relative velocity `urel = u − dr/dt`. Interface values are advected
  with a van Leer–type monotone interpolation (`advectc` for centred quantities,
  `advecti` for interface quantities; `cadv = 2` gives van Leer and `cadv = 0` gives donor cell).
  Optional artificial mass and energy diffusion (`diffuse`; `sigd`, `sige`).
- **Gas momentum** (`gasmom{h,rh}`): pressure gradient, gravity
  `2μπG·xm/r^μ` (equal to GM/r² in spherical geometry), or constant `g` in planar geometry.
  The rh variant adds the radiation force `−χ_F F/c` and the viscous-stress divergence.
- **Gas energy**:
  - `gasnrgh` is the internal-energy equation.
  - `gasnrgrh` is the equation for gas internal energy plus radiation energy. It has the
    flux divergence of F, P dV work with `(p + f E)`, a `(1−3f)` geometric term, and
    viscous heating `q ∂u/∂r`.
  - `gasnrgr` is the static-medium version.
- **Radiation energy** (`radnrg*`): comoving O(v/c) form, with coupling
  `cκ_E E − 4πκ_P B` and `B = σT⁴/π`.
- **Radiation momentum** (`radmom*`): `D(F/c²)/Dt`, plus `∇(fE)`, plus `(3f−1)E/r`, plus
  `ρχ_F F/c`. The interface opacity is the harmonic mean of zone values.
  - `ltran ≥ 2` replaces this equation with a diffusion law `F = −c∇E/(3χρ)`, and sets f = 1/3, g = 1/2.
  - Within that diffusion branch, a flux limiter (`lam = 1`, `xlam = 2/3·lam`) is available.
  - `ltran ≥ 3` additionally replaces the radiation energy equation with `E = aT⁴`
    (equilibrium diffusion).
  - The benchmarks use none of these options; they all use `ltran = 1`.
- **Grey (frequency-integrated)**: there are no frequency groups. [read]
- **Time discretisation**: θ-method with `thet` (0.55 for the hydro tests, 1.0 for the
  radiating shocks). Time-centred values are `x = θ x_new + (1−θ) x_old`.
- **Gravity is off in every benchmark** (`cgrav = 0`). [read]

### 1.4 Jacobian [read]
- **Fully analytic**. Every equation routine writes the residual into `rhs(i,k)` and
  hand-coded partial derivatives into 5 diagonal blocks `em2, em1, e00, ep1, ep2`
  (`meqn×meqn` each, per zone). The derivatives are with respect to the log or
  normalised variables listed above.
- The equation numbers in the comments (c13, gm48, ag42, ef51, …) refer to the
  *TITAN Code Reference Manual*, which is **not** in this distribution.
- Some derivatives are knowingly dropped. Examples:
  - `opac` with `lopac = 5`: "eoneb derivative ignored".
  - Eddington factors are frozen during the Newton iteration when `ltran = 1`.
  - Upwind switches are frozen. [read]

### 1.5 Linear system and solver [read]
The system is block-pentadiagonal: each zone couples to k±2. It is solved by `matslv`:
- **`lband = 0`** (used by every benchmark): custom block-pentadiagonal Gaussian
  elimination (`sgemm`, `sgemv`). Each diagonal block is factored by `sgefad` and
  solved by `sgesld`. These are **LU without pivoting** ("factorizers and solvers have
  no pivoting").
- **`lband = 1`**: maps the system to a LINPACK band matrix (`sgbfa`/`sgbsl`,
  partial pivoting, `ml = mu = 3·neqn − 1`).
- **`lcray = 1`**: CAL-coded Cray solvers (`sgefaj`, `sgeslmj`, `bglsdc`, `bglssl`).
  Their kernels `bglsrd` and `trslbl` are **empty stubs** in `slvpack.f`/`craypack.f`,
  so `lcray = 1` cannot work off a Cray. [read]
- The bundled `blas.f` and `linpack.f` are REAL*8 code under single-precision names
  (`sgemm`, `isamax`, …). This is a Cray naming convention. Do not link them against a
  system BLAS.

### 1.6 Adaptive grid (`grid.f`, `gridinit.f`) [read]
This is a Dorfi–Drury (1987) type grid equation.
- **Point concentration**:
  - `n_k = xscale/Δr_k` (`ladx = 1`), or
  - `n_k = (r_{k+1}+r_k)/(2Δr_k)` (`ladx = 2`).
- **Spatial smoothing**: `ñ_k = n_k − α(α+1)(n_{k−1} − 2n_k + n_{k+1})`, with `alph = α`
  (1.5 or 3 in the tests).
- **Temporal smoothing**: `n̂ = (1+γ)ñ − γñ_old`, with `γ = (τ/Δt)^ibet`. `tau` ranges from
  1e-20 to 1e3 in the tests. `ibet = 1`.
- **Resolution function**: `R_k = sqrt(1 + n_k² Σ_l w_l cs_l²)`. `cs_l` is the difference
  of variable `l` across zone k:
  - linear: `Δy/yscl`
  - logarithmic: `Δy/(y_{k+1}+y_k)`
  - harmonic: `y_{k+1}/y_k − y_k/y_{k+1}`
- **Selectable variables** `lady(l)`: 1 interior mass, 2 ρ, 3 T, 4 E_r, 5 p_g, 6 e_g,
  7 χ, 8 viscous pressure, 9 F_r, 10 `bri` (luminosity or velocity), 11 exterior mass.
  Each has a weight `wt(l)`.
- **Grid equation**: `n̂_k/R_k = n̂_{k−1}/R_{k−1}` for `k = ngrs+2 … ngre−1`. The equations
  next to the boundaries just set `n_k = n_{k−1}`.
- **Alternatives**: `lgrid = 2` (Eulerian, r fixed) and `lgrid = 3` (Lagrangean, dr/dt = u).
- **`gridinit.f`**: problem-specific grid relaxers (`grid3nit` for SOD, `grid4nit`,
  `grid5nit`, `grid6nit`, `grid9nit`). They iterate the grid equation against analytic
  initial profiles, with damping `min(1, dtol/dmax)` and `min(1, ctol/cmax)`.

### 1.7 Eddington factors (`eddfac`, `eddf`, `eddg`) [read]
This is a static (no Doppler or aberration) grey formal solution of the transfer equation.
- **Source function**: `S = [(χρ − σ_T n_e)B + σ_T n_e cE/4π]/(χρ)`. This includes
  coherent Thomson scattering using the *current* E_r.
- **Within each zone**: S is piecewise linear in optical depth, with harmonic-mean slope
  limiting (monotonised).
- **Spherical geometry**: tangent-ray (impact-parameter) method on an augmented mesh of
  zone edges and midpoints. It uses `ncor` core rays (20 in the tests) plus one ray tangent
  to each augmented radius. Rays are integrated analytically, inward and then outward.
  - The moments J and K come from trapezoid-type quadrature in μ (`eddf`).
  - `f = K/J` is clipped to [0, 1].
  - The boundary flux factors `geddl` and `geddr` (`eddg`) are clipped to [0.5, 1].
- **Planar geometry**: a fixed 23-point μ grid (`titan.dat`).
- **`ltran = 1`** (all radiation benchmarks): f is computed once per time step in `peq*`
  from the old state.
- **`ltran = 0`**: f is recomputed in every Newton update.

### 1.8 EOS and opacities [read]
- **`leos = 2`**: ideal gas (γ, μ). Every benchmark uses this.
- **`leos = 1`**: bicubic-Hermite interpolation (`intrp`/`intrc`) in the MHD
  (Mihalas–Hummer–Däppen) tables `*.pg`, `*.eg`, `*.pe`.
  - Table contents: ln p_gas, ln e_gas, ln p_e on ln ρ ∈ [−27.63, −4.61] (51 points;
    ρ = 1e-12 to 1e-2) × ln T ∈ [8.06, 16.12] (176 points; log T = 3.5 to 7.0).
  - Six mixtures: smr2, smr1, solar, pop2a, pop2b, pop2c. All have He/H = 0.1 by number,
    and all contain metals. **There is no Z = 0 table.**
  - Gas only. Radiation is a separate variable.
  - Not exercised by any benchmark. Suspected problems are listed in §1.11.
- **`leos = 3`**: Stellingwerf (1982/84) analytic fit (`eos3`).
- **Opacities**:

| `lopac` | description | notes |
|---|---|---|
| 1 | log-table Rosseland mean from `<name>.cf` | **the `.cf` tables are not distributed**; the Planck-table branch is disabled with `if(.false.)`, so κ_E = κ_P = χ_R − σ_T n_e/ρ |
| 2 | Stellingwerf (1975) fit | |
| 3 | constant χρ | SUBC, SUPC, HEAT, COOL |
| 4 | Thomson | RLST |
| 5 | free-free + Thomson + Compton + crude bound-free | |
| 6 | Kramers + Thomson | |

### 1.9 Artificial viscosity, time step, boundary conditions [read]
- **Artificial viscosity** (`viscous`): tensor (Tscharnuter–Winkler) form.
  - Viscous pressure `q = −cqvis·ρ·μ_q·(∂u/∂r − (μ/4)(u_k/r_k + u_{k+1}/r_{k+1}))`.
  - `μ_q = l(cq1·c_s − cq2·l·min(∇·u, 0))`, with `l = ql0 + ql1·r̄`.
  - It enters momentum as `r^{−μ/2}∂(r^{3μ/2}q)`, and energy as `q ∂u/∂r`.
- **Time step** (`timstp`):
  - `smax` = max relative change over the step of r, m, ρ, T, E. Velocity and flux are
    excluded.
  - `x = stol/smax`; `φ = √x` if x > 1, otherwise x².
  - `tfac = clip(√(tfac·clip(φ, 0.3, 3)), 0.1, 10)`; `dt ← tfac·dt`.
  - There is no CFL limit.
  - The initial dt is tiny (1e-25 s in the shock tests).
- **Hydro boundary conditions**:
  - Eulerian inner (`leibc`) and outer (`leobc`):
    - 1: zero flux
    - 2: prescribed inflow or outflow state (ρ, T, u, mass flux `phi*0`)
    - 3: transmitting (outer only)
  - Lagrangean (`llibc`, `llobc`): 1 is a piston or specified velocity; `llobc = 2` is a
    specified external pressure `pextr` with `omega`.
- **Radiation boundary conditions** `lribc`/`lrobc`:
  - 1: transmitting, with incident intensity `xipl`/`ximr`
  - 2: reflecting
  - 3: imposed net flux (this is how a luminosity is injected, e.g. HEAT)
- **Gravity**: the interior mass `xm = xmtot − xme` can include a central mass. The UPDATE
  notes describe this explicitly. No benchmark uses gravity.

### 1.10 Convection [read]
**There is none.**
- `titan.com` declares turbulence arrays (`et`, `ft`, `cmlt`, `ctvis`, `cdrag`, `cdiss`, `ltur`, `ntur`).
- `titan.flow` lists `turnrg`, `turscl` and `turdrv`, but **those routines do not exist**.
- Every `ltur ≥ 1` branch executes `stop`.
- [inferred] The stubs were probably meant for the Gehmeyr & Winkler (1992, A&A 253)
  time-dependent turbulent-convection model.

### 1.11 Machine-specific hacks, known and suspected bugs
**Machine-specific code** [read]:
- `lcray` paths; Cray `cvmg?` merge intrinsics emulated in `craypack.f`.
- Unit numbers: `itty = 4` opens `/dev/tty`, `iout = 0`; comments say "6 on SGIs".
- Direct-access dump `recl` is in **bytes**. This is right for gfortran and Sun; SGI f77
  would need `-bytereclen`.
- `dimension x(1)` dummy arrays in `eddf`/`eddg`/LINPACK. These trip `-fcheck=bounds`.

**Disabled safety checks** [read]:
- `step.f`: the "changes over timestep too large" rejection is effectively disabled
  (`xchek = 1.0D6`; the original `cvmgt(2.0, 1e3, …)` is commented out).
- `peq*`: `frnom = |F| + 1e-30` (and `unom` likewise) overrides the physically scaled
  nominal flux ("that's the choice for all tests").

**Known bugs from UPDATE**: ISAMAX stride in `maxdel`/`gridinit` and the TOTNRG energy
terms, both fixed in 1.2. The README admits that the identical source gives different
results on different machines (unresolved).

**Found here**:
- `pardoc.f` includes `titan.par` before `titan.imp`, so its PARAMETERs are typed REAL*4.
  This is a hard compile error in gfortran. Fixed (§2). [run]
- `eddfac` planar branch with `lribc = 2`: loop `do 9 kp = ngre, ngrs, -1` uses the stale
  index `k` (`dtau(k)`, `asf(k)`, `bsf(k)`) instead of `kp`. [read, suspected bug; affects
  only planar geometry with a reflecting inner boundary, which no benchmark uses]
- `eos` with `leos = 1`:
  - The tables hold the energy **per unit volume**. README_eos says so, and the table
    columns step by exactly 1 in ln e per 1 in ln ρ.
  - TITAN uses `egn` as **specific** energy (erg/g), e.g. `egrn = egn + ern/dn`, and there
    is no division by ρ.
  - `gamman` is set to the input `gam` instead of Γ₁, and `cpn`, `delsn`, `betan` are not set.
  - Conclusion: [read + inferred] the table EOS path looks unusable as shipped.
- `lopac = 1` needs `.cf` files that are not provided (`opnfil` reads only `<name>.cf`). [read]
- The hydro runs end with gfortran's "IEEE_INVALID_FLAG signalling" note. A build with
  `-ffpe-trap=invalid,zero -O0` ran SOD to completion without trapping, and no
  NaN or Inf appears in any output.
  - [run] It is benign as far as these runs show.
  - [inferred] The flag is probably raised inside a library or runtime routine.

---

## 2. Build notes [run]
- **Compiler**: gfortran 13.3.0 (Ubuntu 24.04), installed with `apt-get install -y gfortran`.
- **Precision**: all code goes through `titan.imp` (`implicit real*8 (a-h,o-z)`), and
  literals use `D` exponents. **`-fdefault-real-8` is not needed and should not be used.**
- **Changes needed** (originals saved as `*.orig`):
  1. `Makefile`: `CC = f77` → `gfortran`; `FFLAGS = -O -cg89` → `-O2 -std=legacy $(XFLAGS)`.
     Build with `make codes XFLAGS="-fno-automatic -finit-local-zero"`.
     The default target `codes` links all objects. `make clean` errors harmlessly if no `.a` files exist.
  2. `pardoc.f`: swapped the order of `include 'titan.par'` and `include 'titan.imp'`.
     This was a compile error.
  3. `opnfil.f`: `/dev/tty` → `/dev/stdout`, because there is no controlling tty in batch
     and the original stops with `opnfil1`.
     - Caveat: tty-unit output still does not show up on stdout.
     - [inferred] This is an interaction between unit 4 and the reopened unit 6, which both
       point at fd 1.
     - Workaround: everything is also in `<P>.out` (`ltalk ≥ 1`), and inputs with an `rtty`
       line write it to `<P>.tty`.
- **Unaffected file**: `code.f` (a monolithic `include` build) is not part of the Makefile
  target and was left untouched.
- **Compiler warnings**: 74 warnings remain. They are `(1)`-dimension bound warnings and
  REAL*4/8 mismatches in the dead Cray stubs and LINPACK `srotg`.
- **Robustness checks** [run]:
  - `-fno-automatic -finit-local-zero` versus automatic storage: identical `SOD.out`.
  - `-O0` versus `-O2`: SUBC and SUPC gave **bit-identical** final dumps.
  - HEAT gave the same jstep and time to every printed digit across `-O0`, `-O2` and
    automatic storage.
  - Not tested: `-march=native`, which allows FMA contraction and could change the bits.

---

## 3. Benchmark results [run]
The problem is selected by the third input line (`lprob` = 3…10), which calls `start{lprob}`.
`titan.inNN` is copied to `titan.in`.
- **COOL** restarts from HEAT's dump 4 (`linit = 4`). `HEAT.dmp` is copied to `COOL.dmp`.
- **RLST** runs in two legs:
  - leg 1: `titan.in10`, 2000 steps;
  - leg 2: `titan.in`, which restarts from dump 19 / jstep 1801 and runs to tmax.
- RLST's `start10` first relaxes a static-radiation model (`lgrid = 3`, `lhydr = 0`), then
  resets t = 0 and deposits `ecnt = 1e47` erg.

Wall times are from a 2-core sandbox with two runs going at once. The 1995 column is user
time on an SGI Indigo2 (RLST: Power Challenge).

| problem | final t (s) | jstep | 1995 t / jstep | match? | ends by | wall (s) | 1995 user (s) |
|---|---|---|---|---|---|---|---|
| SOD | 1.1168 | 454 | 1.12 / 454 | **yes** | tmax | 2.5 | 377 |
| BLST | 5.0026e-2 | 572 | 5.00e-2 / 572 | **yes** | tmax | 10.3 | 587 |
| SDOV | 2.0003e5 | 1097 | 2.00e5 / 1097 | **yes** | tmax | 9.0 | 801 |
| HEAT | **2.9609e15** | 151 | **3.23e15** / 151 | steps yes (trivial: nstep limit); **time 8% low** | nstep | 2.1 | 196 |
| COOL | 2.1113e17 | 551 | 2.11e17 / 551 | yes (nstep-limited; time matches to 3 digits) | nstep | 4.4 | 323 |
| SUBC | 8.7003e4 | 668 | 8.70e4 / 668 | **yes** | tmax | 11.7 | 836 |
| SUPC | 1.5093e4 | 714 | 1.51e4 / 714 | **yes** | tmax | 11.6 | 810 |
| RLST | 6.8095e3 | 3350 | 6.81e**4** / 3350 | **yes** (see note) | tmax | 30.4 + 12.6 | 3183 |

Notes:
- **Step counts**: every tmax-terminated run reproduces the 1995 jstep **exactly**.
  The speed-up is roughly 60–150×.
- **RLST time**: the input has tmax = 6.8e3 and the step count matches exactly.
  [inferred] The table's 6.81E+04 is a typo for 6.81E+03.
- **HEAT time**: the step count is fixed by `nstep`, so it does not test anything. Our time
  after 151 steps is 2.96e15 against 3.23e15 in 1995.
  - The time step oscillates in this problem; dt goes 5.1e14 → 1.3e14 near the end.
    [inferred] Platform round-off probably shifted the dt history. The README warns that
    results differ across machines.
  - COOL restarts with t = 0, so its time does not inherit this difference.

**Regression targets**:
- `reference_runs/<P>/<P>_final.csv`, and `dumps/*.csv` for every archived model.
- Full double precision, taken from the binary dump rather than the 5-digit printout.
- Columns: `k, r, xme_over_4pi, u, Fr` (interfaces) and `rho, T, eg, pg, Er, fedd, chiF,
  kapE, kapP, cs, r_center, Trad` (zones; the last row is a NaN phantom).
- Header: `jstep, timen, dtime, geddl, geddr, xmtot, tetot`.

Caveats on these targets:
- `fedd` in a dump is the value used during the last step. With `ltran = 1` it was computed
  from the *start-of-step* state.
- RLST_leg2 dumps 1–19 come from leg 1. COOL dumps 1–4 are HEAT's.
- The time series of energy terms is in `*.nrg`.

---

## 4. Julia port plan (module by module)

**Size**: about 13k non-comment Fortran lines, excluding BLAS/LINPACK (~1k lines replaced by
the standard library). Of those, ~4.6k lines are problem setups (`start*.f`, `gridinit.f`).
The physics and numerics core is ~5k lines.

**Guiding choice: residual-first.**
- Port each `rhs(...) =` residual expression exactly.
- Get the Jacobian blocks from **ForwardDiff on local stencils**, i.e. per zone the
  dependence on zones k−2…k+2. This gives exactly the block-pentadiagonal pattern.
  Use the same unknowns: ln r, ln m_e, ln ρ, ln T, ln E, u/u_nom, F/F_nom. Matching these
  keeps the Newton iteration counts comparable.
- Keep the Fortran analytic Jacobian **only as a test oracle**. Add a tiny Fortran hook that
  dumps `em2…ep2, rhs` after `matgen` for one step, and compare to the AD Jacobian.
  Expect deliberate differences where TITAN drops terms (§1.4).

| Fortran | Julia module | keep / replace |
|---|---|---|
| `titan.com`, `.par`, `.eq*`, `.imp`, `.dat` | `State` and `Params` structs (new/old/centred arrays as fields; constants in `Const`) | replace; remove EQUIVALENCE aliasing |
| `titan.f`, `step.f`, `peq*.f`, `update*.f`, `reset.f`, `maxdel.f` | `Integrator`: predictor, Newton loop, boost iteration, retry logic | keep logic exactly (convergence measure, dt × 0.05 on failure); merge the h/r/rh triplicates using multiple dispatch on a `Physics` type |
| `matgen.f` | `assemble!`: loops zones and calls the equation kernels | replace with AD-assembled blocks |
| `mass`, `contin`, `gasmom*`, `gasnrg*`, `radnrg*`, `radmom*`, `viscous`, `diffuse`, `advectc`/`advecti` | `Equations`: one residual function per equation, on a stencil view | keep the discretisation verbatim; drop the hand derivatives (~60% of these files) |
| `grid.f` | `AdaptiveGrid`: resolution functions plus the Dorfi–Drury equation | keep; table-driven resolution variables instead of 11 copy-pasted blocks |
| `gridinit.f` | `GridInit`: one generic relaxer (damped Newton on the grid equation with fixed profiles) | replace the 5 copies with one |
| `matslv.f`, `slvpack.f`, `linpack.f`, `blas.f`, `craypack.f` | `LinearSolve` | replace. First choice: BandedMatrices.jl, LAPACK `gbtrf`/`gbtrs` with pivoting (N ≲ 2500, bandwidth ≲ 23, negligible cost). Optional later: a block-pentadiagonal solver with StaticArrays and **pivoted** block LU. The Fortran default has no pivoting. |
| `eddfac`, `eddf`, `eddg` | `Transfer`: tangent-ray formal solver returning f(r), g_l, g_r | keep the algorithm; add a planar fix for the stale-`k` bug; unit-test on analytic cases (isothermal sphere, f → 1/3 at depth) |
| `eos.f`, `readeos.f`, `intrp.f`, `intrc.f` | `EOS` interface: `eos(model, ρ, T) → (p, e, ∂ln p/∂ln ρ, ∂ln p/∂ln T, ∂ln e/∂…, Γ₁, n_e)` | ideal gas: keep. MHD tables: port the monotone bicubic Hermite only if wanted, and fix the per-volume/per-mass bug. Stellingwerf: optional. |
| `opac.f`, `readopac.f` | `Opacity` interface: `(χ_F, κ_E, κ_P, derivatives)` | keep lopac 3–6; replace 1 and 2 with real tables (§5) |
| `timstp.f` | `TimeStep` | keep exactly (needed for jstep regression) |
| `archive`, `writout`, `writsty`, `totnrg`, `opnfil`, `clzfil`, `pardoc`, `hdfinit`, `hstore` | `IO`: JLD2/HDF5 snapshots, CSV, an energy log, and a reader for the Fortran `.dmp` (port `read_dump.py`) | replace; keep `totnrg` as an energy diagnostic |
| `start03`–`start10` | `problems/*.jl`, one file per benchmark | port each as a test case |
| `start01`, `start02` | — | blank templates; drop |

**Test strategy**: see §6, the three-tier verification policy. Bitwise parity is a debugging aid only.

**Effort estimate** (honest; one experienced developer, AI-assisted; ±50%):
| stage | time |
|---|---|
| hydro core plus grid, reproducing SOD, BLST, SDOV | 2–3 weeks |
| radiation (transfer, static and rh variants), reproducing HEAT, COOL, SUBC, SUPC, RLST | 2–3 weeks |
| quasi-star additions excluding convection (§5) | 3–5 weeks |
| time-dependent convection and its validation | 4–8 weeks (the largest uncertainty) |
| **total to a validated quasi-star-capable code** | **~3–4 months** |
| reproduction of the existing benchmarks alone | **~1–1.5 months** |

**Main risks**:
1. **Untested code paths**:
   - No benchmark exercises gravity, a central mass, an Eulerian inflow boundary
     (`leobc = 2`), a non-ideal EOS or tabulated opacities. These are exactly the paths the
     quasi-star needs.
   - "Reproducing TITAN" will validate the solver, not that physics.
   - New tests are needed: a hydrostatic polytrope or envelope, a Bondi/free-fall accretion
     column, Saha ionisation fronts, a Marshak wave.
2. **Newton robustness depends on TITAN's conventions**:
   - log variables;
   - |x|+1e-30 normalisation of u and F;
   - convergence judged without u and F;
   - the boost iteration;
   - grid parameters (α, τ, weights) tuned per problem.
   Deviating from any of these may change convergence behaviour.
3. **Bitwise trajectory divergence** (pivoted versus unpivoted LU, AD round-off) leads to
   different dt histories. The HEAT case shows this sensitivity even between machines in 1995.
4. **Missing manual**: the Reference Manual is not in the tarball, so the equation numbers
   cannot be checked against the derivations.
5. **Physics limits**: grey, O(v/c), static formal solution for f. This is adequate for
   optically thick envelopes. It is questionable where the accretion shock radiates into an
   optically thin, fast inflow (v/c ~ 0.1 near a massive seed).

---

## 5. Physics to ADD for a quasi-star accretion-shock problem
| need | status in TITAN 1.2 | what to add |
|---|---|---|
| Spherical geometry | **present** (`lgeom = 2`) [read, run] | none |
| Central point mass | **partly present** [read]: gravity uses `xm = xmtot − xme`, so a core mass enters by offsetting `xmn(ngrs)`. Never exercised (cgrav = 0 in all tests). | make the core mass explicit (and possibly time-dependent via the accretion rate through the inner boundary); add a hydrostatic-equilibrium test |
| Central luminosity source | **partly present**: `lribc = 3` holds the inner flux at its initial value (HEAT uses `xlum`) | L(t) or an L tied to accretion onto the core; this changes the inner-boundary F equation and its Jacobian entries |
| H/He Saha EOS with radiation | **absent**: ideal gas only in practice; MHD tables have metals and the suspected bug | Saha for H, He, He⁺ (optionally H₂ and H⁻ at low T) with ionisation energy in e and consistent Γ₁. Radiation stays the separate E_r variable, but provide an LTE (T_r = T) consistency check. Derivatives with ForwardDiff through the Saha solve (implicit-function rule, or a Newton loop with AD). |
| Z = 0 Rosseland and Planck means | **absent** | Tables in (ρ, T). Candidate sources, to verify: OPAL GN93 X = 0.75, Z = 0 for log T ≳ 3.75; low-T primordial means (e.g. Mayer & Duschl 2005; Lenzuni, Chernoff & Salpeter 1991). Separate Planck (κ_P, κ_E) from Rosseland (χ_F), which TITAN supports structurally but its table branch ignores. Add a monotone bicubic interpolator. |
| Free-fall inflow outer boundary | **mechanism present, untested**: `leobc = 2` (prescribed dr, ur, tr, mass flux `phir0`, momentum and energy fluxes `phir1`, `phir2`) plus `lrobc = 1` with incident `ximr` | wire u = −√(2GM/r), ρ from Ṁ = 4πr²ρ\|u\|, kinetic plus enthalpy energy flux; check sign conventions of `phir*` against `dmdt`; test on a pressureless or Bondi column |
| Time-dependent convection | **absent**: stubs only (`ltur`, `et`, `turnrg`/`turscl`/`turdrv` missing) | Add a turbulent kinetic-energy (or convective-flux) equation as an extra block variable; `meqn = 8` already reserves a slot. Candidates: Kuhfuß (1986) one-equation model; Gehmeyr & Winkler (1992), likely what the TITAN stubs were for [inferred]; Wuchterl & Feuchtinger (1998) / Feuchtinger (1999), implemented in an implicit adaptive-grid RHD code of the same Dorfi–Drury lineage. Needs a turbulent flux in the energy equation, turbulent pressure and viscosity in momentum, a source/dissipation balance, and new resolution weights in the grid equation. |
| Large dynamic range in mass | **present**: the v1.2 exterior-mass variable exists for exactly this case (accretion onto a compact object, steep envelopes) [read] | none |

---

---

## 6. Verification policy (three tiers)
Tom's decision (phase 2): the Julia port does **not** need to be bitwise identical to TITAN.
It may use Julia's own math, `muladd`/FMA, AD Jacobians, restructured loops and better
linear solvers.
- The bitwise machinery of Phase 1 is kept only as a **translation-debugging aid**: glibc
  parity mode, byte-compare of `.dmp`/`.out`, the trace executable.
- Acceptance is defined by the three tiers below.
- Every tolerance is a *proposal* justified from TITAN's own numerics. None is a physical
  requirement.

### (a) Translation check, per routine
- **Test**: same inputs (a COMMON snapshot from a benchmark state) into the Fortran routine
  and the Julia routine; compare every output array.
- **Pass**: max relative difference ≤ 1e-12 on outputs.
- **Exception for differences of large terms**: here the difference is measured relative to the
  largest term. Example: the continuity residual is O(1e17), formed from terms of O(1e22).
- **Debugging**: run parity mode (glibc `exp`/`log`/`pow`, gfortran multiplication order) only
  when a check fails and the cause has to be localised.
- **Tools**: `phase2_run.jl … tierb` (`tier_a`), `contin_ad.jl`, `harness_eos.jl`.

### (b) Per-step check
- **Test**: from an *identical* full COMMON snapshot, take one step with the reference
  routines and one step with the candidate routines (`TD.snapshot()` / `TD.restore!()`).
  Compare the new-time-level state, the Newton iteration count and the next dt.
- **Pass**:
  - max relative difference in r, m, ρ, T, E_r ≤ `conv` (= 1e-5 in every benchmark);
  - |Δu|/max|u| and |ΔF|/max|F| ≤ `conv`.
- **Why this tolerance**:
  - Newton is declared converged when the largest relative correction to (ln r, ln m, ln ρ, ln T,
    ln E) is ≤ `conv`. One extra boost iteration follows.
  - Two correct implementations of the same discrete equations therefore agree to within the
    accuracy that TITAN itself guarantees: `conv` in those variables.
  - u and F are excluded from TITAN's convergence measure. They are normalised by their own
    maxima here.
- **Iteration counts**: report them. A difference in count is allowed but flagged.
- **Expectation**: round-off-level changes, such as reassociated EOS arithmetic, should give
  ~1e-12; anything near 1e-6 indicates a real discrepancy. This is inferred from the
  quadratic convergence of Newton near the solution.

### (c) Benchmark acceptance on physical observables (not bits)
Observables are computed by `julia/Metrics.jl` from per-step snapshots. They are linearly
interpolated in time to fixed target times, because two runs never land on the same times.
Proposed thresholds:

| observable | threshold |
|---|---|
| positions (shock, contact, rarefaction head/tail, Zel'dovich-spike position, precursor extent) | \|Δ\| ≤ 0.5% of a length scale L: L = 1 (box) for SOD; L = shocked-layer width (shock − piston) for SUBC/SUPC |
| temperatures, T_rad and plateau values (T_spike, T_post, T_rad at the shock, T_pre, p\*, u\*, ρ\*) | ≤ 1% relative |
| energy-error diagnostic E_err = (tetot − tetot₀)/max(\|tee\|, \|tew\|) | \|Δ\| ≤ 1e-3 absolute |
| step count to tmax | within 3% |

Why these thresholds:
- The time-step history is chaotic at round-off level. HEAT's dt history already differed
  between 1995 machines.
- Step counts and late-time profiles can therefore drift without any error, while shock
  positions and temperatures are robust.
- The SOD reference itself deviates from the exact Riemann solution by about 1e-3 in
  shock/contact position and 0.03% in p\*. That bounds what "the same answer" can mean.

**Observable definitions**:
- **Shock (radiating shocks)**: steepest outward density drop per unit length, refined to the
  half-jump crossing.
- **SOD shock**: steepest pressure drop per unit length.
- **SOD contact**: maximum of (Δln ρ − Δln p)/Δx left of the shock, i.e. a density jump
  without a pressure jump.
- **Rarefaction head**: first zone with ρ < 0.99 ρ_L. **Tail**: first zone with ρ ≤ 1.01 ρ\*_L.
- **p\*, u\***: values midway between contact and shock.
- **T_spike**: max T.
- **T_post, T_rad,post**: values at r_sh − 0.1·W, where W is the shocked-layer width.
- **T_pre, T_rad,pre**: values at r_sh + 0.05·W.
- **T_rad at the shock**: (E/a)^{1/4} at r_sh.
- **Precursor extent**: outermost zone with T > 170 K (twice the hottest initial gas
  temperature), minus r_sh.

Limitations of these definitions:
- E_err takes the state after step 1 as its zero point.
- The precursor metric saturates once the whole domain exceeds 170 K. In SUBC that is the case
  at the final time, so the target times stop at 6e4 s.
- The SUPC target times (4e3, 7.5e3, 1.3e4 s) are my choice. The comment in `start09.f` lists
  the SUBC times and looks copy-pasted.

---

## Phase 1: shared library + harness

New directories:
- `lib/`: `Makefile`, `src/titan_lib.f`, `src/titan_trace.f`, the build outputs,
  `common_map.txt` and the verification logs.
- `julia/`: `TitanFortran.jl`, `verify_phase1.jl`, `harness_eos.jl`.

Nothing in `fortran/` was changed for this phase.

### P1.1 Build [run]
```
cd /home/claude/TitanPort/lib && make -j1          # libtitan.so, titan.x, titan_trace.x
# flags: -O2 -std=legacy -fPIC -ffp-contract=off -finit-local-zero   (no -march, no -fno-automatic)
# make STATIC=1 adds -fno-automatic (not needed)
```
- **Sources**: all objects are compiled once from `../fortran/*.f` into `lib/obj/`.
  - `libtitan.so` = core objects + `src/titan_lib.f`.
  - `titan.x` = core objects + the **unmodified** `fortran/titan.f`. The original
    `fortran/Makefile` also still works.
  - `titan_trace.x` = core objects + `titan_trace.f`. It is a copy of `titan.f` that also
    streams COMMON `/new/` after each of the first 50 steps.
- **Entry points** (`titan_lib.f`; gfortran symbols `titan_init_`, `titan_step_`, `titan_finalize_`).
  The bodies are the `titan.f` call sequence, split up:
  - `titan_init(istat)`: `opnfil`, then `start`, then `jstep = jsteps`.
  - `titan_step(istat)`: `step`, `archive`, `writout`, `writsty`; then the tmax test and
    `jstep += 1`. Returns istat 0 (continue), 1 (t ≥ tmax) or 2 (jstep > jstepe).
  - `titan_finalize()`: closes the remaining units.
  - The `titan.dat` DATA statements moved into a `BLOCK DATA titanbd`.
  - The problem is selected by the input deck, exactly as before: the wrapper copies
    `titan.inNN` to `<rundir>/titan.in` and `cd`s into `rundir`. File I/O is unchanged
    (`.out`, `.dmp`, `.nrg`, `.sty`).
- **SAVE semantics** [run]. There are no SAVE statements; DATA appears only for read-only
  constants (`eos3`, `calcei`). Built *without* `-fno-automatic`, all 9 benchmark legs
  (`scripts/check_exe.sh`) gave **byte-identical `.dmp` files**. The only difference was in
  `.out`: `writout` prints a local array element, `odepth(ngrs-1)`, that is never assigned.
  `-finit-local-zero` makes that printout 0 again, and then `.out` is byte-identical too.
  So `-fno-automatic` is not needed. `-finit-local-zero` is kept only for printout parity.
- **FMA**: `-ffp-contract=off`, with no `-march`.
  - On baseline x86-64 there is no FMA instruction anyway, so this is a guard.
  - The new flags (PIC, contract-off, automatic storage) still reproduce every
    reference `.dmp` and `.out` byte for byte. [run]

### P1.2 COMMON-block map [run]
`TitanFortran.parse_commons` computes the layout by parsing `titan.com` and `titan.par`.
- Types follow `titan.imp`: i–n is `Int32`, everything else `Float64`, and `header` is
  character*80.
- Layout is packed, in declaration order; Fortran column-major order maps directly onto
  Julia arrays.
- All 30 blocks match the ELF symbol sizes from `nm -S`:
  - 25 match exactly;
  - `agrid`, `geom`, `integrat`, `turb2` and `matrix` end with integers and carry only
    trailing padding to a multiple of 8.
- The full offset table is in `lib/common_map.txt`.
- Access is zero-copy: `TitanFortran.var(:rn)` returns `unsafe_wrap` over `dlsym(lib, :new_)`
  plus the offset.

Main blocks:

| block | contents |
|---|---|
| `/new/` | 41 arrays(330): rn, rmun, rmup1n, rmum1n, rcn, xmn, xmen, xmcn, dvoln, dn, dsn, un, usn, tn, asn, egn, egsn, pgn, xnen, ern, ersn, egrn, egrsn, etn, etsn, egtn, egtsn, frn, frsn, egrtn, egrtsn, avchin, betan, delsn, cpn, chifn, xken, xkpn, plfn, bri, gamman |
| `/old/` | the same quantities with suffix `o` |
| `/mid/` | time-centred values, plus fedd, geddl, geddr |
| `/matrix/` | em2, em1, e00, ep1, ep2 (8,8,330); rhs (8,330); band matrix bm; ngrs, ngre, neqn |
| `/integrat/` | sx, smax, stol, timeo, timen, dtime, tfac; jstep, jsteps, jstepe, jext, next, jback, nback |
| `/iterate/` | dx, cx, dmax, dtol, cmax, ctol, conv; kx, iter, niter, jtry, ntry |
| `/logic/` | every switch (`lprob`, `lgeom`, `leos`, `lopac`, …) |
| `/const/` | physical constants |
| `/eostab/`, `/opactab/` | tables and derivatives |

`get_state()` returns:
- at interfaces ngrs..ngre+1: r, xme, u, F;
- in zones ngrs..ngre: ρ, T, e_g, p_g, E_r, f, χ_F, κ_E, κ_P, c_s;
- plus jstep, iter, time, dtime.

### P1.3 Bitwise verification through Julia [run]
Run with `julia verify_phase1.jl SOD 3`, `… SUBC 8` and `… SUPC 9 50`, each in its own
process. The log is `lib/verify_phase1.log`.

| problem | final jstep (ref) | live COMMON vs reference dump (63 quantities, full 330-length arrays, compared as UInt64) | files written by the Julia-driven run vs reference |
|---|---|---|---|
| SOD | 454 (454) | all identical | `.dmp`, `.out`, `.nrg` byte-identical |
| SUBC | 668 (668) | identical except `fedd[4]` | `.dmp`, `.out`, `.nrg`, `.sty` byte-identical |
| SUPC | 714 (714) | identical except `fedd[4]`, `fedd[105]` | `.dmp`, `.out`, `.nrg`, `.sty` byte-identical |

- **The `fedd` exceptions are expected, not drift.** They are the phantom entries
  `ngrs−1` and `ngre+1`. `writout` assigns `fedd(ngre+1) = geddr` and
  `fedd(ngrs−1) = geddl` *after* `archive` has written the dump. [read + run]
- **Per-step SUPC check**: for the first 50 steps, the Julia-driven run and the standalone
  `titan_trace.x` agree bitwise on jstep, iter, timen, dtime and the whole 108,240-byte
  `/new/` block: **50/50 identical**.

### P1.4 Per-routine harness [run]
`harness_eos.jl` does the following:
1. Runs SUBC to the end.
2. Zeros a routine's outputs and calls the single Fortran subroutine with `ccall` on that
   state.
3. Runs a hand translation in Julia with the Fortran operation order copied.
4. Compares bitwise over k = ngrs−1..ngre+2 (103 entries).
5. The switch `parity = true` routes `sqrt` and `pow` through glibc (`ccall((:pow, "libm.so.6"), …)`).

The log is `lib/harness_eos.log`.

| routine | input | parity = false (Julia math) | parity = true (glibc) |
|---|---|---|---|
| `eos`, leos = 2 (20 outputs) | SUBC final state | identical | identical |
| `eos`, leos = 2 | T perturbed by ±37% | identical | identical |
| `opac`, lopac = 6 (12 outputs; extra libm demo, lopac forced to 6 on the SUBC state) | snapshot | **χ, κ_E, κ_P, dlnχ/dlnT differ in 2 of 103 zones** | identical |
| `opac`, lopac = 6 | T perturbed | 3 outputs differ in 2 zones | identical |

The ideal-gas `eos` is bitwise reproducible without parity mode because its only
non-arithmetic call is `sqrt`. gfortran inlines that as `sqrtsd`, and it is correctly rounded
in both Julia and glibc.

**Which calls each object file makes** (from `objdump` relocations) [run]:

| call | where |
|---|---|
| `exp` | eddfac (8), eos (8, table path), gridinit (24), opac (9) |
| `log` | eos (4), opac (4) |
| glibc `pow` (real exponents) | gasmom{h,rh} (`r**(1.5·xmu)`), opac (10), peq* and update* (`**rxm1`), several start routines |
| libgcc `__powidf2` (**integer-variable** exponents, e.g. `r**mup1`, `r**mu`) | grid, peq*, update*, totnrg, start* |
| inline multiplications (**constant** integer powers) | `x**3 = x·(x·x)`, `x**4 = (x·x)·(x·x)`, `x**5 = x³·x²` |
| inline `sqrtsd` | sqrt |

**Micro-checks over 10⁶ arguments** (number of results differing in any bit) [run]:

| comparison | differing results |
|---|---|
| exp: Julia vs glibc | 9,730 |
| log: Julia vs glibc | 83 |
| `pow(x, −3.5)`: Julia vs glibc | 9,770 |
| `pow(x, 0.25)`: Julia vs glibc | 9,757 |
| sqrt | 0 |
| Julia `x^n` (runtime n = 2 or 3) vs `__powidf2` | 0 |
| Julia literal `x^3` vs `x*(x*x)` | 0 |
| **Julia literal `x^4` vs gfortran's `(x*x)*(x*x)`** | **499,739 (50%)** |

Rules that follow for a bitwise port. Since phase 2 these apply only in debugging/parity mode (§6):
1. In parity mode, call glibc `exp`, `log` and `pow` through `ccall`. Julia's `exp`, `log`
   and `^(Float64, Float64)` are its own implementations and differ in about 1% of cases.
2. Write constant integer powers as gfortran's multiplication chains, e.g.
   `t2 = t*t; t2*t2` for `t**4` in `plfn = σT⁴/π`. Never use `x^4`; Julia evaluates it with
   extra precision.
3. Integer-variable powers can use `x^n` for the n values that occur (1, 2, 3), or call
   `__powidf2` directly.
4. Keep the Fortran left-to-right association, e.g. `ck/(cm0*gmu)*dn*tn`.
5. Keep `-ffp-contract=off` on the Fortran side. On the Julia side, never use `@fastmath`
   or `muladd` in parity code.
6. [inferred] GPU kernels (KernelAbstractions, phase 3) will not reproduce glibc bits for
   exp/log/pow, and GPU compilers contract to FMA by default. Plan GPU verification as
   tolerance-based against the bitwise-verified CPU Julia version.

### P1.5 Re-entrancy and limits [read + run]
- **One model per process on the Fortran side.** All state lives in global COMMON blocks.
  Fortran I/O units are global, with fixed numbers 4–10 taken from the input deck. COMMON
  is not re-zeroed between models, so a second `init` in the same process is
  unsupported. The verification scripts therefore start one Julia process per problem.
- **A Fortran `STOP` terminates Julia as well.** Examples: non-convergence `step01`,
  `opnfil*`, `start*`. Mitigation: run each Fortran model in a subprocess or a Distributed
  worker.
- **Working directory**: `init` changes the process working directory to `rundir`, and
  `finalize` restores it. Not thread-safe.
- **The hooks overwrite COMMON in place.** After using `setscalar!` or `var(...) .=` (as
  the harness does), the live model is no longer the reference trajectory.
- **Not done in phase 1**: no routine has been replaced by Julia inside a running model. The
  harness compares a routine in isolation on a snapshot. Phase 2 needs a Julia-callable
  hook to substitute routines inside `step`; the obvious approach is that Julia writes its
  outputs into COMMON and the Fortran call is skipped.

---

## Phase 2: Julia step driver, swapped EOS/opacity, AD Jacobian

### Files
| file | contents |
|---|---|
| `julia/TitanDriver.jl` | `step!()`: `step.f` re-implemented in Julia (predictor `peq`; Newton loop `matgen → matslv → maxdel → update`; convergence test `dmax ≤ conv`; divergence retry with dt × 0.05 up to `ntry`; boost iteration; `timstp`; the "changes too large" branch; `totnrg`) |
| | `titan_step!()`: the `titan.f` loop body (`archive`, `writout`, `writsty`, tmax/jstepe tests) |
| | swap points: `ROUTINES[:matgen] = f`, etc. for routines called by step; `install_hooks!(eos = f, opac = g)` for routines called *inside* Fortran |
| | `snapshot()` / `restore!()` copy all 30 COMMON blocks |
| `lib/src/hooks.c` | the deep-level hook mechanism. In `libtitan.so`, `eos_`/`opac_` are C shims that call a Julia `@cfunction` if one is set, otherwise the original Fortran, which `objcopy` renamed to `eos_f_`/`opac_f_`. `titan.x` is unaffected. |
| `lib/src/port_msg.f` | `step.f`'s WRITE/STOP statements, so `.out` stays identical |
| `julia/JuliaPhysics.jl` | Julia-native `eos!` (leos = 2) and `opac!` (lopac = 3, 4, 6), written for speed and clarity rather than parity: precomputed c_v and R, `c_s = sqrt(γRT)`, ρ cancelled analytically in n_e, `T^-3.5` as `1/(T³√T)`. Other options fall back to Fortran. |
| `julia/Metrics.jl`, `julia/tierc.jl` | tier-(c) observables and acceptance |
| `julia/phase2_run.jl` | runs (modes `fortran`, `julia`, `tierb`) |
| `julia/contin_ad.jl` | continuity residual and its ForwardDiff Jacobian |
| `julia/env/` | project environment containing ForwardDiff 1.4.6 |
| `reference_runs/metrics_reference.csv` | reference observables at the target times |
| `lib/p2_runs/` | logs and per-step snapshots |

### P2.1 Scaffold check: Julia driver with all-Fortran routines [run]
SOD, SUBC and SUPC produce `.dmp`, `.out`, `.nrg` and `.sty` files **byte-identical** to the
standalone references. Step counts are 454, 668 and 714.
- **Exercised**: the divergence/retry path, which appears 14–32 times per run in `.out`.
- **Not exercised**: the "changes too large" (`jback`) branch, which TITAN disables
  (`xchek = 1e6`), and the STOP paths.

### P2.2 Julia EOS/opacity swapped into the running model [run]
The Julia `eos!` and `opac!` are active from t = 0, including the initial-model setup.

**Tier (a)**:
- `eos`: max relative difference 2.9e-16 across 8 outputs × 103 entries, at 35 checkpoints.
- `opac` (lopac = 3): 0, across 6 outputs at 25 checkpoints (the SUBC and SUPC ones).

**Tier (b)** at 35 checkpoints (SOD 10, SUBC 12, SUPC 13), one step from an identical snapshot:
- all pass;
- max relative difference in (r, m, ρ, T, E) = 8.5e-12;
- |Δu|/max|u| ≤ 7.6e-13; |ΔF|/max|F| ≤ 3.7e-12; Δdt/dt ≤ 9.2e-11;
- **identical Newton iteration counts at every checkpoint**.

**Tier (c)**: all three problems **ACCEPT**. Each cell shows the maximum deviation of the
Julia-swapped run from the Fortran reference over all target times.

| problem | steps ref / Julia | positions \|Δ\|/L | temperatures and plateau values (relative) | E_err \|Δ\| |
|---|---|---|---|---|
| SOD | 453 / 452 (−0.2%) | 8.0e-9 | 2.5e-10 (≤ 2.4e-8 for p\*, u\*) | 1.7e-12 |
| SUBC | 667 / 667 | 6.5e-11 | 9.0e-11 | 8.2e-13 |
| SUPC | 713 / 713 | 8.1e-12 | 8.1e-13 | 1.8e-12 |

SOD needs one fewer step. Through t = 0.2 the trajectories agree to ~1e-8. They separate later,
after the shock has reflected, and the final step overshoots tmax by a different amount
(t_end = 1.1245 versus 1.1168).

**Reference observables** (Fortran TITAN; also in `reference_runs/metrics_reference.csv`).

SOD (planar, γ = 1.4). Exact Riemann values: shock 1.6752 / 1.8504; contact 1.5928 / 1.6855;
head 1.3817 / 1.2634; tail 1.4930 / 1.4860; p\* 0.30313; u\* 0.92745; ρ\*_L 0.42632; ρ\*_R 0.26557.

| problem | t | shock | contact | head | tail | p\* | u\* | ρ\*_L | ρ\*_R | E_err |
|---|---|---|---|---|---|---|---|---|---|---|
| SOD | 0.1 | 1.6762 | 1.5932 | 1.3737 | 1.5043 | 0.3031 | 0.9273 | 0.4261 | 0.2656 | −6.18e-05 |
| SOD | 0.2 | 1.8514 | 1.6860 | 1.2499 | 1.5102 | 0.3030 | 0.9273 | 0.4260 | 0.2656 | −1.17e-04 |

Reading the SOD comparison with the exact solution:
- TITAN matches it to ~1e-3 in shock and contact position, and to 0.03% in p\*, u\*, ρ\*.
- The rarefaction head and tail differ by ~0.01–0.02. This is mostly an artefact of my 1%
  thresholds on a smoothly smeared fan.

Radiating shocks:

| problem | t (s) | r_shock (cm) | W (cm) | T_spike | T_post | T_rad,post | T_rad,shock | T_pre | T_rad,pre | precursor (cm) | E_err |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SUBC | 5.4e3 | 8.042959e11 | 1.046e9 | 860.3 | 848.6 | 655.2 | 627.1 | 130.0 | 606.9 | 2.36e7 | −0.069 |
| SUBC | 1.7e4 | 8.138441e11 | 3.620e9 | 885.9 | 861.7 | 744.3 | 696.6 | 200.0 | 654.4 | 4.49e8 | −0.059 |
| SUBC | 2.8e4 | 8.228424e11 | 6.018e9 | 868.5 | 841.4 | 755.9 | 697.1 | 195.9 | 639.9 | 4.73e8 | −0.040 |
| SUBC | 3.8e4 | 8.309962e11 | 8.156e9 | 847.1 | 819.3 | 752.7 | 687.0 | 177.5 | 618.9 | 3.52e8 | −0.044 |
| SUBC | 6.0e4 | 8.487553e11 | 1.269e10 | 795.4 | 771.7 | 729.8 | 654.3 | 129.5 | 568.0 | 3.23e7 | −0.041 |
| SUPC | 4.0e3 | 8.091270e11 | 1.120e9 | 6206.9 | 3985.5 | 3826.3 | 3924.6 | 3843.7 | 3846.1 | 1.84e10 | −0.023 |
| SUPC | 7.5e3 | 8.173132e11 | 2.299e9 | 6363.0 | 4197.8 | 4194.9 | 4296.3 | 4205.6 | 4213.0 | 2.79e10 | −0.024 |
| SUPC | 1.3e4 | 8.303217e11 | 4.299e9 | 6484.2 | 4540.3 | 4542.8 | 4634.4 | 4549.9 | 4552.7 | 3.77e10 | −0.021 |

Temperatures are in K.

Qualitative physics checks (inferred; not compared numerically with Ensman 1994 or other literature):
- **SUPC** shows the supercritical signature T_pre ≈ T_post, and a Zel'dovich spike at
  ~6200–6500 K.
- **SUBC** shows T_pre ≪ T_post.

**TITAN's own energy diagnostic** is not small:
- SUBC: E_err of 4–7% of the piston work.
- SUPC: about 2%.
- SOD: about 1e-4.

Whether this is real non-conservation of the internal-energy formulation or an imperfect
`totnrg` diagnostic is **not determined**. UPDATE notes earlier `totnrg` bugs. Tier (c)
therefore compares E_err *against the reference*, not against zero.

### P2.3 Julia continuity residual with a ForwardDiff Jacobian [run]
`contin_ad.jl` ports, in Julia:
- the continuity residual (`contin.f`) for zones ngrs..ngre, including both boundary rows;
- the van Leer interface densities (`advectc.f`);
- the mass-diffusion flux (`diffuse.f`, zet = 0).

Unknowns: ln r, ln ρ and u/u_nom on k = ngrs−2..ngre+2, i.e. TITAN's own Newton variables.

**Test state**: run the model to step N; do one Newton iteration of step N+1 so the state is
consistent but unconverged; then call `matgen_` for the Fortran residual and blocks
`em2…ep2`. The Julia residual and `ForwardDiff.jacobian` are evaluated on the same state.

| case | residual \|ΔR\| / largest term | Jacobian entries compared | max \|ΔJ\| / row max | AD ≠ 0 but Fortran = 0 | Fortran ≠ 0 but AD = 0 |
|---|---|---|---|---|---|
| SUBC, step 300 | 3.9e-13 | 791 | 1.8e-11 | 0 | 0 |
| SUPC, step 400 | 2.9e-13 | 791 | 7.3e-11 | 0 | 0 |
| SOD, step 100 (planar, sigd = 1e-5) | 1e-18 | 991 | 2.0e-12 | 0 | 0 |

- **Result**: TITAN's hand-coded continuity Jacobian is complete; it drops no terms.
- **Residual normalisation**: relative to the residual itself the agreement is 1e-8 to 2e-7.
  The residual is a near-cancellation of O(1e22) terms, so the largest term is the right scale.
- **Frozen switches**: the upwind switch `s(k)` and the limiter's sign test are piecewise
  constant, and AD and the Fortran treat them identically.
- **Not covered**: the phantom-zone boundary rows of the id equation.
- **Dropped terms I know of from reading** [read; not yet measured]:
  - lopac = 5 ignores the derivative of E₁;
  - with ltran = 1, Eddington factors are frozen during the Newton iteration (a cross-routine
    dependence);
  - with ltran = 0 their derivatives are absent even though f is updated each iteration.
  - These will show up when `radmom*` and `radnrg*` are ported the same way.

### P2.4 What remains
- **Port the remaining residual routines** the same way: mass, gasmom*, gasnrg*, radmom*,
  radnrg*, viscous, grid. Then assemble the full block-pentadiagonal Jacobian with
  ForwardDiff (or a sparsity-colored AD over 5-zone stencils), and replace `matgen`.
- **Replace `matslv`** with a pivoted banded LU (BandedMatrices.jl / LAPACK `gbtrf`).
  The Fortran block solver does not pivot.
- **Port in Julia**: peq*/update* (they call eos/opac and set the time-centred quantities),
  `timstp`, `eddfac` (tangent rays), and `archive`/`writout`.
- **Remove all COMMON dependence** once every routine is in Julia. Then multiple models per
  process and KernelAbstractions become possible.
- **Tolerance thresholds**: the tier (b)/(c) thresholds are proposals. They should be revisited
  when a change is expected to move results by more than round-off, e.g. a different
  linear solver or AD terms TITAN drops.
- **Not verified**:
  - the hook mechanism with routines other than eos/opac;
  - problems other than SOD, SUBC and SUPC under the Julia driver (BLST, SDOV, HEAT, COOL, RLST);
  - any GPU path.

---

## Phase 3: energy diagnostic, Jacobian completeness, Julia solver/assembly

New files:
- **Julia scripts** (`julia/`):
  - `energy_ledger.jl`, `energy_table.jl`, `energy_decompose.jl`: independent energy budget;
  - `fd_jacobian.jl`: coloured FD of the Fortran residual vs the analytic blocks;
  - `eddington_newton2.jl`: Newton effect of the frozen Eddington factors;
  - `lopac5_check.jl`;
  - `JuliaSolver.jl`: banded solver plus coloured-AD continuity and mass rows;
  - `phase3_run.jl`;
  - `driver_all.jl`.
- **Library additions** (`lib/src/`): `start03/08/09.f` copies with an optional `TITAN_NCELL`
  override (no effect when unset; defaults verified byte-identical); `port_env_ncell` in
  `port_msg.f`.
- **Logs**: `lib/p3_runs/`.

### P3.0 Load-order bug in the harness (found and fixed) [run]
- **Symptom**: a Julia-driven run diverged from the reference in the last bits whenever any
  library symbol was looked up before `init`.
- **Cause**, shown with `LD_DEBUG=bindings`: libtitan's undefined `pow`/`exp`/`log` were bound
  either to glibc libm or to **Julia's `libopenlibm`**, depending on whether OpenLibm was
  already in the global symbol scope at `dlopen` time.
- **Fix**: `TitanFortran.load()` now uses `RTLD_DEEPBIND` and runs in `__init__`, so libtitan
  always binds its own dependency (glibc libm).
- **Re-verification after the fix**: SOD/SUBC/SUPC through the Julia driver are again
  byte-identical (`.dmp`/`.out`/`.nrg`).
- **Phase 1 results** (`verify_phase1.jl`, `harness_eos.jl`) were obtained in the safe load
  order (glibc); they stand.
- **Phase 2 `julia` mode runs** installed hooks before `init`, so their *Fortran* routines used
  OpenLibm's pow/exp/log. Their tier-(c) results (≤ 1e-8) were obtained in that configuration.
  This makes them only a looser test; the conclusions stand.

### P3.1 Energy-diagnostic investigation [run unless marked]
**Conclusion 1: TITAN's own `tetot` diagnostic is wrong; the reported 4%/2.6% errors are almost
entirely bookkeeping.**

Evidence:
- **Kinetic + internal + radiation energy agrees exactly.** My independent ledger (from the
  per-step COMMON state, using the dt actually taken) reproduces TITAN's `tee` to all digits.
- **The discrepancy is entirely in `tew`, the boundary-work term.** From step 1 to the end,
  TITAN's `tew` exceeds the independently integrated piston work by 3.7% (SUBC) and 2.6% (SUPC).
- **Mechanism, verified per step.** `step.f` calls `timstp` *before* `totnrg`, so `totnrg`
  multiplies every step's fluxes by the dt *planned for the next step*.
  - In normal steps the per-step ratio TITAN/true equals dt_next/dt_used to 4 digits.
  - That shift alone is harmless: it changes W by about 1e-4.
  - The damage is done at **divergence retries**. The planned dt is booked, but the step
    actually taken has 0.05 × that dt. In SUBC (14 retries) this accounts for essentially
    the whole excess.
- **Second, independent bug: the clock.** In the divergence branch, `step.f` sets
  `timen = timeo + dtime` and then calls `reset`, whose first statement is `timen = timeo`.
  The retried step is integrated over dtime, but the clock does not advance.
  - SUBC loses 170 s (0.2% of t_end) in 14 retries; SUPC loses 11 s in 6; SOD loses 0.0054
    (0.5%) in 113 retried steps.
  - Consequences: all TITAN times after a retry are early by the lost amount. This includes
    the 1995 BENCHMARKS times and the observable times in §6, which use TITAN's clock.
  - The Julia driver has an opt-in fix, `TD.FIX_CLOCK[] = true` (env `TITAN_FIX_CLOCK=1` in
    `phase3_run.jl`), off by default so that TITAN is reproduced bitwise.
  - With the fix, SUBC reaches tmax in 587 instead of 667 steps (t_end 87009 s). Observables
    at fixed times with the corrected clock were not yet compared.
- **Radiation term**: TITAN's `tel` agrees with my radiative boundary losses to within the same
  next-dt booking (SUBC: 2.72e-3 vs 2.87e-3 of W).

**Conclusion 2: the scheme's genuine non-conservation is about 10× smaller and is a
spatial-discretisation error that converges with resolution.**

Independent residual R = ΔE − W + L at the final time:

| problem | residual | normalisation |
|---|---|---|
| SUBC | −2.8e-3 | piston work W |
| SUPC | −3.1e-3 | piston work W |
| SOD (closed box) | −1.5e-3 | E₀ |

Same quantity at fixed physical times (`lib/p3_runs/energy_table.txt`):

| run | SUBC t=5.4e3 | SUBC t=3.8e4 | SUBC t=8e4 | SUPC t=4e3 | SUPC t=1.3e4 | SOD t=1.1 |
|---|---|---|---|---|---|---|
| base (100 zones) | −2.2e-2 | −4.1e-3 | −2.7e-3 | −9.8e-3 | −3.5e-3 | −1.5e-3 |
| conv 1e-5 → 1e-7 | unchanged | unchanged | unchanged | unchanged | unchanged | unchanged |
| stol 0.1 → 0.05 | −1.9e-2 | −3.8e-3 | −2.6e-3 | −8.7e-3 | −3.2e-3 | −1.4e-3 |
| stol → 0.025 | −1.8e-2 | −3.6e-3 | −2.5e-3 | −7.9e-3 | −2.9e-3 | −1.3e-3 |
| 200 zones | −1.2e-2 | −2.0e-3 | −1.1e-3 | −5.4e-3 | −1.8e-3 | −1.8e-4 |
| 300 zones | −8.6e-3 | −1.4e-3 | −7.5e-4 | −3.9e-3 | −1.3e-3 | n/a |
| θ = 0.75 (SUBC) / θ = 1 (SOD) | −2.0e-2 | −3.8e-3 | −2.6e-3 | n/a | n/a | −3.5e-3 |
| Lagrangian grid (lgrid = 3) | −0.11 | −0.16 | −0.17 | −5.2e-2 (t = 4e3) | n/a | −1.2e-2 |

What the table shows:
- **Not the Newton tolerance**: conv has no effect.
- **Weak dependence on dt**: halving `stol` twice reduces the residual by only ~20%.
- **Close to first-order convergence in the zone count**: exponent ≈ 0.8–0.9 for SUBC/SUPC;
  SOD's ×8 drop from 100 to 200 zones is faster.
- **The adaptive grid helps**: the Lagrangian grid, coarse at the shock, is far worse. So the
  adaptive remap is not the culprit.
- **Localised in time**:
  - SUBC/SUPC: 40% / 68% of the final residual accrues in steps 8–13 (t ≈ 1e-22–2e-21 s), the
    impulsive piston start and shock formation.
  - SOD: it accrues during the wall reflection (t ≈ 0.285, dt ≈ 1e-11).
  - Spatial localisation was not computed. [not done]

Decomposition of the residual (`energy_decompose.jl`; `lib/p3_runs/energy_decompose.txt`):

| component | finding |
|---|---|
| **artificial viscosity** | Heating −Σ q (∂u/∂r) dV and viscous work Σ u r⁻¹ Δ(r³ q), in TITAN's discrete forms, balance to ≤ 4e-9 (1e-11 for SUBC/SUPC) of the scale. **The viscous heating is not mis-accounted.** |
| **time integration** | The θ-scheme kinetic-energy term −(θ−½) Σ m̄ (Δu)² gives 20% (SUBC), 25% (SUPC), 22% (SOD, θ = 0.55) of the residual, and 75% for SOD at θ = 1. |
| **remaining 60–80%** | Spatial discretisation of the non-conservative internal-energy formulation (advection/shock capturing). Not decomposed further. **[inferred]** from the zone-count convergence and its presence in pure-hydro SOD. |
| **O(v/c) radiation terms** | Not separately isolated. **[inferred, magnitude bound]** v/c ≤ 7e-5 in SUBC/SUPC, and radiative losses are 0.3% / 3% of W, so any O(v/c) mis-accounting is ≤ 1e-6 of W. The same residual also appears in radiation-free SOD. |

**What the port needs**:
1. **Fix the diagnostic**: book each step's fluxes with the dt *actually used*, i.e. before
   `timstp` and after a retry; keep the independent ledger as a regression observable.
2. **Fix the clock** (`FIX_CLOCK`).
3. **For the quasi-star problem** (energy-dominated, strong accretion shock):
   - TITAN's internal-energy formulation loses ~1% of the shock-processed energy at 100–300
     zones near shock formation.
   - Either require shock resolution (zone convergence of the energy residual as an
     acceptance test), or move to a total-energy-conserving formulation.
   - **[recommendation; not implemented]** Two candidate formulations: evolve total energy
     with internal energy recovered; or add the discrete kinetic-energy defect back as heat.

### P3.2 Completeness of TITAN's analytic Jacobian [run]
**Method**:
- Central finite differences of the *actual Fortran residual*: perturb unknowns in COMMON,
  call `update_` with a zero correction to recompute EOS, opacity and time-centred quantities,
  then call `matgen_`.
- 5-colour stencil.
- Step sizes chosen per variable (r: 1e-10, m: 1e-11, others: 1e-6). Larger steps in r move an
  interface by a sizeable fraction of a zone and are useless.
- Checked at two step sizes; compared over rows k = ngrs+2..ngre−2.

Results (`lib/p3_runs/fd_jacobian.log`; SUBC step 300, SUPC step 400, SOD step 100):

| equation | max \|J_FD − J_an\| / row max |
|---|---|
| grid, mass, gasmom, gasnrg, radnrg, radmom (at frozen f) | ≤ 1e-6, all variables. This is FD round-off: errors fall as h grows. |
| contin vs r | 8e-4 / 5e-3. Also round-off: it falls with larger h. Coloured AD confirms contin to ~1e-11 at these states. |
| contin vs ρ (all problems); SOD gasmom vs u, gasnrg vs T | 1e-6 to 0.1, and h-dependent. **[inferred]** These are kinks of the van Leer limiter (FD across a switch). For contin, AD at the same SOD state shows exact agreement; gasmom and gasnrg were not ported, so for them this is not proven. |

So **at frozen Eddington factors, TITAN's hand-coded Jacobian is complete** for all seven
equations in SUBC/SUPC and five in SOD, apart from the limiter-kink caveat above.

**Dropped term 1: the Eddington-factor response ∂f/∂x** (uncoloured FD, because the formal
solution is nonlocal; `eddington_newton.jl`, SUBC step 300):

| equation | ∂f/∂x terms relative to row max |
|---|---|
| radmom | up to 1.9e2 inside the ±2-zone band, 17 outside it (11% of the term's ℓ1 mass lies outside the band) |
| radnrg | 1e-5 |
| gasnrg | 4e-8 |
| others | 0 |

The frozen-f Jacobian itself has nothing outside the band (4e-10). The effect on Newton,
measured in TITAN's full block system (dense pivoted LU with row equilibration, corrections
applied by Fortran `maxdel`/`update`; `eddington_newton2.jl`):

| state | C: ltran = 1 (as the benchmarks) | A: ltran = 0, TITAN J | B: ltran = 0, J + ∂f/∂x |
|---|---|---|---|
| SUBC step 100 | 2 | 3 | 2 |
| SUBC step 300 | 4 | 5 (convergence turns linear) | 4 |
| SUPC step 400 | 3 | 5 | 4 |

- With `ltran = 1` the dropped term does not matter, because f is part of the frozen-f
  discrete system and the Jacobian is consistent with it.
- With `ltran = 0`, dropping ∂f/∂x costs 1–2 extra iterations; including it restores
  quadratic convergence.
- The term is dense, so including it needs either a dense or low-rank update, or keeping its
  ±2-zone banded part (89% of its ℓ1 mass). The banded-part variant was not tested.

**Dropped term 2: lopac = 5 opacity derivatives** (`lopac5_check.jl`; analytic d ln κ/d ln T vs
FD of the same Fortran routine; ρ = 1e-9 g cm⁻³; T = 1e4–1e7 K):

| T (K) | d ln κ_P/d ln T, analytic | FD |
|---|---|---|
| 1e4 | +232 | −3.6 |
| 3e4 | +27 | −3.7 |
| 1e5 | −0.1 | −3.8 |
| ≥ 3e5 | ≈ correct (−3.3 to −3.7) | |

- The analytic bound-free derivative omits exactly the E₁′ term (the "eoneb derivative
  ignored" comment). At large b = 1.6e5/T, the remaining term no longer cancels and the error
  becomes enormous. Shown by hand differentiation [read] plus the numbers above [run].
- The flux-mean derivative has the same error. κ_E's derivative is off by up to 3.0. From
  reading the code: the Compton derivative `dxkccdlt = 2.676e-10` lacks the factor T, and κ_E
  ignores the bound-free derivative.
- No benchmark uses lopac = 5, so the iteration-count effect was not measured. The port should
  differentiate opacities with AD.

### P3.3 Port status and solver [run]
**Routine-by-routine status**:

| Fortran routine | Julia status | verification |
|---|---|---|
| step.f control flow, titan.f loop | ported (`TitanDriver.step!`) | byte-identical on **all 9 benchmark legs** (P3.4) |
| eos (leos = 2) | ported, Julia-native math | tier a 2.9e-16; tier b/c pass |
| opac (lopac 3, 4, 6) | ported | tier a 0 (lopac 3); 3 and 6 used in runs |
| matslv | **replaced**: pivoted banded LU (LAPACK gbtrf/gbtrs) with row equilibration | SOD/SUBC/SUPC: same step counts, tier c passes |
| contin + advectc + diffuse (mass branch) | ported; residual + stencil-coloured ForwardDiff Jacobian (15 partials, one sweep), assembled into the running system | residual 4e-13 of the largest term; J vs analytic ≤ 7e-11; tier b passes |
| mass | ported, coloured AD | J vs analytic ≤ 4e-16; residual identical (SUBC/SUPC/SOD) |
| gasmom*, gasnrg*, radnrg*, radmom*, grid, viscous, advecti | **not ported** | analytic Jacobians verified complete by FD (P3.2), except the Eddington-factor response in radmom |
| peq*, update*, timstp, maxdel, reset, eddfac, totnrg, archive/writout | Fortran | n/a |

**Stability fixes for the Julia solver**:
- The first version (no equilibration) hit an *exactly* singular pivot.
- The matrix rows span about 1e-30 to 1e30. The 1e-30 comes from `unom = |u| + 1e-30`
  at u = 0.
- Row equilibration fixes it; the dense LU needed the same treatment.

**Timing**: Julia driver, `nice`d, 2 cores shared with other jobs. Treat these as ±20%.

| routine | Fortran | Julia |
|---|---|---|
| matslv | 0.15–0.30 ms/call | 0.42–0.63 ms/call (band LU, n = 525–728) |
| matgen | 0.65 ms/call | 3.7 ms/call with Julia contin + mass rows (2.6 ms with contin only) |

- The extra matgen time is mostly the Julia continuity residual (2 ms/call). It allocates
  mgr-length Dual arrays per call; it is not yet optimised.
- **Per step**: Fortran 6.0–6.7 ms; `jsolve` 7.7–9.3 ms; `jall` (solver + eos/opac + contin
  and mass rows) 22–27 ms.
- **The Julia version is currently about 4× slower per step.** Expected, since the ported
  residuals are unoptimised and allocate per call; no optimisation has been attempted yet.
- **Profile of the all-Fortran run** (SUBC): matgen 49%, matslv 22%, peq 10% (Eddington
  factors), driver and I/O 16%.

**Tier results** (`jall` = Julia solver + Julia eos/opac + Julia contin and mass rows):
- **Tier (b)**:
  - with mass rows included: SUBC 9/9, SUPC 10/10;
  - SOD 8/8, run with contin rows only (before the mass rows were added);
  - max relative state difference 1.2e-9; identical iteration counts at every checkpoint.
- **Tier (c)**: SUBC and SUPC ACCEPT with identical step counts.
  - SUBC: positions ≤ 1.9e-5 of the layer width, temperatures ≤ 2.6e-5.
  - SUPC: both ≤ 5e-8.
- **SOD**: all observables ACCEPT (positions ≤ 9e-4, plateau values ≤ 1.2e-4), but the **step
  count fails the 3% criterion**.

SOD step count for variants that pass tier (b), compared with round-off perturbations of the
Fortran model:

| run | steps to tmax |
|---|---|
| reference | 453 |
| jsolve | 453 |
| jeos | 423 |
| jall | 503 |
| jcontin | 538 |
| Fortran + one-time 1e-14 T perturbation (6 seeds) | 451–468 |
| Fortran + 1e-14 T noise every step (6 seeds) | 448–453 |

- The count is governed by the number of divergence retries (27–47) during the wall
  reflection.
- **[open]** Why the Julia variants spread more widely than random noise was not resolved.
  Per-step agreement is ~1e-12, so no defect is indicated.
- **Proposal**: drop the step-count criterion for SOD, or replace TITAN's crude retry (dt × 0.05
  after 15 failed iterations) with a damped / line-search Newton in the port, and re-evaluate.
- **Metric fix**: the rarefaction head/tail metrics were changed to interpolated crossings.
  They had jumped by a whole zone.

### P3.4 Complete reference set under the Julia driver [run]
BLST, SDOV, HEAT, COOL (restarting from the HEAT dump), RLST leg 1 and RLST leg 2 all run
through the Julia driver with all-Fortran routines. Every TITAN output file (`.dmp`, `.out`,
`.nrg`, and `.sty`/`.tty` where present) is **byte-identical** to the standalone references.
With SOD, SUBC and SUPC, all 9 benchmark legs are now covered (`lib/p3_runs/driver_all.out`).

### P3.5 What remains before quasi-star physics
1. **Port the rest**: residuals gasmom*, gasnrg*, radnrg*, radmom*, grid, viscous, advecti
   with coloured AD; then peq/update, timstp, eddfac. Optimise the Julia residuals so they
   allocate nothing per call.
2. **Fix TITAN's bugs** in the port: clock and energy booking (P3.1), lopac = 5 derivatives
   (P3.2), planar `eddfac` stale index (§1.11). Decide on `ltran` (frozen f) vs live f with
   ∂f/∂x.
3. **Gravity with a central point mass**:
   - The Fortran path exists, but no benchmark exercises it (cgrav = 0 everywhere).
   - New test: a **hydrostatic envelope** (polytrope or radiative envelope) held static over
     many dynamical times. Acceptance: |u|/c_s and the energy residual stay small.
4. **Central luminosity source**:
   - `lribc = 3` holds the inner flux fixed; a time-dependent or accretion-powered L is needed.
   - Test: a radiative-equilibrium envelope reaching the analytic L(r) = const.
5. **Z = 0 EOS and opacity from `/home/claude/QuasiStarInfall/src`** [read]:
   - The EOS (`eos.jl`) is ideal H/He gas, X = 0.76, Y = 0.24, with Saha ionisation. Its
     charge balance is solved by bisection, and its thermodynamic derivatives are FD.
   - For TITAN it must also provide gas specific *internal energy*, including ionisation
     energy (not present in `eos.jl`), plus ∂ln p/∂ln ρ, ∂ln p/∂ln T, ∂ln e/∂…, Γ₁ and n_e.
     Radiation stays TITAN's separate E_r.
   - The derivatives should come from the implicit-function theorem on the charge balance.
     ForwardDiff through a bisection returns zero derivatives.
   - The opacity (`opacity.jl`) is a Rosseland mean only, from continuum sources without
     lines (stated there as a lower bound around 1e4–1e5 K), interpolated bilinearly in log.
     TITAN additionally needs Planck/absorption means (κ_P, κ_E) and a χ_F that includes
     scattering.
   - Bilinear interpolation gives derivatives that jump at cell edges. Use a monotone bicubic
     (port TITAN's `intrc`/`intrp`) or smooth splines for Newton.
6. **Free-fall inflow outer boundary**:
   - `leobc = 2` (prescribed ρ, T, u and mass flux `phir0`, plus momentum/energy fluxes
     `phir1`/`phir2`) exists but is untested.
   - Wire it to v_ff = √(2GM/r), with ρ from Ṁ, matching `boundary.jl`'s
     `rho_infall`/`vff`, and check the sign conventions of `phir*` against `dmdt`.
   - Test: a **free-fall column**, pressureless or Bondi-like, steady state against the
     analytic ρ ∝ r^−3/2, with mass and energy flux conservation through the domain.
7. **Marshak wave test** for the radiation coupling and the diffusion limit: analytic
   front position ∝ t^½ for constant opacity.
8. **Re-validation**: re-baseline the energy ledger (P3.1) as an acceptance observable for
   every new test.

## Phase 4: bookkeeping fixes, total-energy prototype, full Julia residuals, solver, damped Newton
(Work done 15–16 Oct 2026. The working agent was interrupted before writing this section, so it was written from the logs in
`lib/p4_runs/` and re-runs from the git repository. **[run]** unless marked.)

### P4.1 Bookkeeping fixes: ON by default
- `TitanDriver.FIX_CLOCK` and `FIX_ENERGY` are on by default. `TITAN_BITWISE=1` turns both off and reproduces TITAN 1.2 bitwise.
- Neither fix changes the solution of any step. They change the time labels, so the step at which t_max is reached, and the
  energy diagnostic.
- With the fixes, SUBC reaches t_max in 587 steps instead of 667, and SOD in 452 instead of 453.
- New reference observables for all 9 legs: `reference_runs/metrics_reference_fixed.csv` (from `legcompare.jl`). The originals are
  in `metrics_reference_orig_all.csv`.
- `tierc.jl` takes `TIERC_REF=std` to compare against the fixed all-Fortran run (Julia driver). The default remains the bitwise
  `fortran` run, and each reference writes to its own CSV.

### P4.2 Total-energy formulation (`TitanResidFast.ETOTAL[] = true`)
The gas-energy row is replaced by a conservative total-energy row: gas + kinetic + radiation, with fluxes in finite-volume form
on the adaptive grid. Internal energy is recovered by subtraction.

| run | steps | Newton it/step | cumulative energy imbalance / scale |
|---|---|---|---|
| SOD std / etot | 471 / 497 | 5.77 / 5.82 | −1.6e-3 / −6e-14 |
| SUBC std / etot | 667 / 777 | 4.87 / 4.87 | −2.8e-3 / +2.5e-13 |
| SUPC std / etot | 713 / **1699** | 4.29 / 3.66 | −3.2e-3 / +5.1e-13 |

- **Conservation:** to round-off.
- **SOD wall heating removed:** at the reflecting wall (t = 0.32), the wall-zone ρ/ρ_exact goes from 1.033 to 0.999 and
  T/T_exact from 0.980 to 1.0008.
- **Accuracy at fixed resolution** (`etot_resolution.txt`): at 100 zones, SUBC observables are closer to the 200-zone answer.
  For example, at t = 3.8e4 the shock position error is 1.6e-2 (std) vs 3.2e-3 (etot) of the width, and T_post is 2.6e-3 vs 8e-5.
  At 200 zones the two schemes agree to ≤ 6e-3.
- **Tier (c) against the internal-energy reference fails** (positions 1–5 % of L). This is expected for a *different discretisation*.
  The criterion for etot is convergence with resolution, not agreement with the old scheme.
- **Cost:** SUPC needs 2.4× the steps (2× the Newton iterations) through shock formation. The time step is still limited by the
  relative change in T, and in cold, kinetic-energy-dominated zones T comes from the small difference E − KE.
- **Dual-energy switch** (`ETASW`, internal-energy row where e_int/(e_int + KE) < η): with η = 0.3 the step counts are normal
  (SUBC 697, SUPC 727), but the switched zones give back the conservation (SUBC −9.8e-4; SUPC at η = 0.1: −1.9e-3).
- **[recommendation]** Keep the pure total-energy row and change the time-step control instead (limit the relative change in
  E_tot, or in e_int only where e_int/KE is not small). For the quasi-star accretion shock, where energy is the whole point, this
  is the preferred direction. Not yet implemented.

### P4.3 All residuals in Julia (`TitanResidFast`) and static block-pentadiagonal solver (`BlockSolver.Penta`)
- **Residuals:** all 7 equations (grid, mass, continuity, gas momentum, gas energy, radiation energy, radiation momentum), with
  one stencil-coloured ForwardDiff sweep.
  - `ncol = 6` is the default. With 5 colours, the single ±3-block dependence (radmom: ρ(k−3) through the F/ρ interface
    interpolation) aliases into the +2 entry, giving a Jacobian error of 1.3e-3 of the row max on SOD.
  - Jacobian vs the Fortran analytic one: ≤ 2.3e-11 of the row max. RHS ≤ 3e-18.
  - Side outputs (slopes etc. left in COMMON): identical to Fortran.
- **Solver:** a static 7×7 block-pentadiagonal LU with row equilibration. Agrees with Fortran matslv to 7e-11 and with LAPACK
  gbtrf to 1e-10. Takes 0.20 ms/call vs 0.30 for Fortran.
- **"jfast" mode** (residuals + solver + EOS/opacity in Julia; everything else Fortran), against `std` (all-Fortran, fixes on):
  - SOD, SUBC and SUPC have identical step counts and final times.
  - Tier (c) observables agree to ≤ 1e-12.
- **Speed** (steady state, excluding JIT): Julia 4.1–5.1 ms/step vs Fortran 6.4–6.7.
  - First-call compilation is about 11 s for matgen and 3.5 s for matslv. A PackageCompiler sysimage, or precompiling
    a package, would remove that.
- **`TitanResidRef` (the slower reference implementation) is deprecated.** Its side outputs differ from Fortran by up to 1e-2
  (Jacobian and RHS agree); `TitanResidFast` is exact. Not investigated further.

### P4.4 Damped / line-search Newton (`DampedNewton`, experimental, OFF)
- **SOD** (fixed bookkeeping, full Julia stack): 428 steps vs 470 with TITAN's dt/20 retry. Divergence retries go from 32 to 25,
  and Newton iterations from 2709 to 2495. Every tier-(c) observable passes; only the step-count criterion flags (−5 %).
- **SUBC and SUPC:** abort with `SingularException` in the static block LU at a line-search trial state, where a 7×7 diagonal
  block is exactly singular. Not resolved.
- The end-time spread between "seeds" (1e-14 perturbations) is the final-step overshoot past t_max, not a physics difference.

## Phase 5: removing the remaining Fortran (in progress)
New files: `julia/TitanCore.jl` (translations), `julia/PortHarness.jl` (in-flight checks), `julia/phase5_check.jl`,
`scripts/phase5_jcore.sh`.

### P5.1 reset, update{h,rh}, maxdel, timstp, totnrg → Julia [run]
- **Method:** the harness wraps a routine so that on every call it runs the Fortran version and the Julia version from the same
  COMMON snapshot and compares **every** COMMON variable (Float64: max |Δ|/max|x|; integers: exact). The run continues on the
  Fortran state, so every call is an isolated test on a real state.

| routine | SUPC (rad-hydro), 713 steps | SOD (hydro), 470 steps |
|---|---|---|
| maxdel | 3057 calls, identical | 2709 calls, identical |
| timstp | 713 calls, identical | 470 calls, identical |
| totnrg | 713 calls, identical | 470 calls, identical |
| reset | 6 calls, identical | 32 calls, identical |
| update | 3057 calls; all variables ≤ 1e-13 except urel 6.6e-7, drdt 3.2e-7 | 2709 calls, identical |

- **The two update outliers are one ulp, amplified by cancellation.** The boundary-cell radius uses a real-exponent power
  (`rn^rxm1`), where Julia's `pow` differs from glibc's in the last bit. drdt = (rn − ro)/dt and urel then amplify that.
  Verified: with `PARITY=1` (glibc pow through ccall, `TitanCore.PARITY[]`) every variable is ≤ 4e-16.
- **"jcore" mode** (jfast + these five routines in Julia; only peq, I/O and initialisation remain Fortran), against `std`:
  - SUBC 587/587 steps, SUPC 713/713 steps: tier (c) ACCEPT.
  - SOD: observables ACCEPT; step count 470 vs 452, the same as jfast (known SOD retry sensitivity, §P3.3).
- **Speed:** jcore runs at 4.4–5.9 ms/step vs Fortran 6.4.
- **Lesson:** the first version was 30× slower in `update!`. The cause was `Int(TD.g(s))` on an untyped COMMON scalar, which is
  not inferred as `Int`, so every `rhs[ir, k]` was a dynamic dispatch. Concrete accessor types (`gi(s)::Int`, `vf(s)::Vector{Float64}`)
  fixed it. Keep this pattern for all COMMON access.

### P5.2 Remaining Fortran in a run
1. **`peq` → `peqh` / `peqrh`:** old ← new copies, extrapolation of radius and mass with the monotonicity retry, boundary values,
   EOS/opacity, and the Eddington-factor call.
2. **`eddfac`:** the formal solution; only called when `ltran = 0`, which no benchmark uses.
3. **I/O:** `archive`, `writout`, `writsty`, `clzfil`.
4. **Initialisation:** `start*`, `gridinit`, `readeos` / `readopac`, `opnfil`.

After (1) and (3), a benchmark run needs Fortran only for initialisation. Initialisation can then be replaced per problem by Julia
setup functions, which is also where the quasi-star problem setup belongs.

### P5.3 peq{h,rh} and eddfac → Julia: the whole time step runs in Julia [run]
- **Refactor:** `update!` is split into the Newton correction plus `_ghost_and_derived!(mode)`, which `peq!` shares. The Fortran
  peq and update tails differ in only three places, all reproduced:
  - peq always calls eddfac;
  - peqrh computes `plfn = csigr/cpi*T^4` where updaterh computes `csigr*T^4/cpi`;
  - peq also zeroes the `egt*` arrays (and peqh zeroes `xne*`).
- **`peq!`:** new→old copies; extrapolation of r and m with TITAN's halving retry (the retry counter `jext` is kept in COMMON,
  as in the Fortran); u and F boundary values; then the shared tail.
- **`eddfac!`:** planar (μ = 0) and spherical (μ = 2) formal solutions, with `eddf`/`eddg` moments.
  - Translated faithfully, including two TITAN 1.2 quirks. Planar lribc = 2: the return sweep uses the stale index k = ngre+1
    instead of kp; `FIX_EDDFAC_PLANAR[]` fixes it and is off by default. Spherical lribc = 3: uses the planar quadrature
    `ang(j)`; not changed.
  - SUBC and SUPC are spherical with lribc = 2, so neither quirk affects them.
- **In-flight checks**, every call over the whole run:

| routine | SUPC | SUBC | SOD |
|---|---|---|---|
| peq | 713 calls, ≤ 2e-11 (urel/drdt as in P5.1, the rest ≤ 2e-13) | — | 470 calls, identical |
| eddfac | 713 calls, ≤ 9e-13 (xmom2) | 587 calls, ≤ 4e-13 | n/a (hydro) |

- **"jstep" mode:** every routine of the time step in Julia (residuals, Jacobian, solver, EOS, opacity, peq, eddfac, update,
  maxdel, timstp, totnrg, reset). Only output (`archive`/`writout`/`writsty`) and initialisation are still Fortran.
  - Against `std`: SUBC 587/587 steps and SUPC 713/713 steps, tier (c) ACCEPT.
  - SOD: observables ACCEPT; 470 vs 452 steps (§P3.3).
  - Speed: 4.5 (SOD), 6.1 (SUPC) and 6.3 (SUBC) ms/step vs Fortran 6.4 (SUPC), excluding JIT.

## Phase 6: no Fortran at run time [run]
- **State backend** (`TitanFortran`): `BACKEND[] ∈ (:common, :julia)`.
  - On `:julia`, `var(name)` returns a Julia-owned array with the same name, type and shape as the COMMON variable. Every
    module reaches state through `var`, so the routines run unchanged on either backend.
  - `to_julia!()` / `to_common!()` copy between the two.
  - Modules that cache views register a clearer in `TF.ON_SWITCH`: `TitanCore.VC`, `JuliaPhysics.V`, `TitanDriver.SV`.
  - `snapshot` / `restore!` work on both backends.
- **State files:** `save_state` / `load_state!` use a simple Julia-version-independent format ("TITANST1": name, type, dims,
  raw data), gzipped. A benchmark's complete initial state is 35–63 kB. `TITAN_NOLIB=1` stops libtitan being opened at all.
- **Driver:**
  - `TitanDriver.run_julia!` is titan.f's loop without archive/writout/writsty.
  - `msg()` replaces port_msg on the Julia backend; the fatal codes raise Julia errors instead of STOP.
  - `TD.EOS` / `TD.OPAC` are the EOS and opacity entry points.
- **`phase6_run.jl`:**
  - `capture NAME NN` runs TITAN's Fortran set-up once and stores the state in `reference_runs/init/NAME.tstate.gz`.
  - `pure NAME` loads that state and runs every routine in Julia.

| problem | pure Julia | fixed Fortran reference | tier (c) | ms/step (steady) |
|---|---|---|---|---|
| SOD | 455 steps | 452 | ACCEPT | 2.4 |
| SUBC | 587 | 587 | ACCEPT | 4.7 |
| SUPC | 713 | 713 | ACCEPT | 4.3–5.0 |

- **Why SOD changed** (455 here vs 470 for jstep): in jstep the Julia EOS hooks were already installed during TITAN's set-up,
  whereas the captured state was set up with the Fortran EOS. SOD's step count is sensitive to such last-bit differences (§P3.3).
- **Remaining Fortran:** problem set-up (`start*.f`, `gridinit.f`, `readeos`/`readopac`) and the original text/dump outputs. The
  benchmarks keep using the captured states. New problems get Julia set-up functions that fill the same state; the quasi-star
  accretion problem is the first.

## Phase 7: quasi-star physics (in progress)
Design decisions (Tom, 16 Oct 2026):
- **Convection:** a local, time-independent MLT convective flux (α = 0.8, as in MESA-QUEST and the steady model) in the
  energy equation.
- **Domain:** an outer shell (inner boundary at ~0.3 R with imposed L and the steady-model state) through the photosphere into
  the infall.

### P7.1 Pluggable materials; Z = 0 Saha gas + tabulated opacities [run]
- **`Z0Physics.jl`:** the AD-ready Z = 0 EOS and opacity from QuasiStarInfall.jl v2 (validation record in
  `docs/Z0_EOS_OPACITY.md`). The opacity table is stored as gzipped raw Float64 (`julia/z0/kaptable_Z0.bin.gz`, 7 MB),
  independent of the Julia version; the round trip is bit-identical.
- **`TitanMaterials.jl`:** `Builtin` (TITAN's ideal gas + lopac 3/4) and `Z0Gas`.
  - `Z0Gas` provides gas-only pressure and internal energy (translational + ionization; radiation stays in TITAN's separate
    E_r), κ_flux = Rosseland (incl. scattering), κ_E = κ_P = Planck (absorption), and n_e for the eddfac scattering term.
  - a_s (used only in the artificial viscosity) is the monatomic estimate √(5p/3ρ). **[approximation]**
  - `use!(m)` installs the material in the residual (`TitanResidFast.MATERIAL`, generic in the Dual type) and in the state
    fills used by peq!/update!.
- **Nested AD:** the residual's hand-made tag `JTag` is declared ≺ every ForwardDiff.Tag, so the Newton polish inside the Saha
  solve differentiates correctly within the residual sweep. Unused ghost zones (ρ = T = 0) are skipped.
- **Checks:**
  - Z0 coloured AD Jacobian = dense AD Jacobian exactly (728 unknowns); dense AD vs central FD ≤ 1.5e-6 of the row max.
  - Assembly 0.9 ms/call (Builtin 0.5).
  - Benchmarks through the material interface are unchanged (SUPC pure: 713 steps, same final time).

### P7.2 MLT convective flux in the energy equation [run]
- **Where it enters:** with `TitanResidFast.CONVECTION[] = true` (α from `ALPHA_MLT`, default 0.8), the energy row gets
  r^μ F_conv at every interior interface. That is the gas+radiation energy row for radiation runs, the gas-energy row otherwise,
  or the total-energy row. F_conv = 0 at the two boundary interfaces.
- **Flux:** F_conv = ρ c_P T √(gδ) ℓ² ζ³ / (4√2 H_P^{3/2}), with ζ = √(U² + ∇ − ∇_ad) − U.
  - This is Kippenhahn & Weigert's MLT, solved for the flux at a given actual gradient, which is what a time-dependent code
    needs. It is smooth (∝ ζ³) at the onset ∇ = ∇_ad.
  - ∇ = Δln T / Δln P from the neighbouring zones at the n level, with P = p_gas + f E_r.
  - g is the same expression as in the momentum equation; no gravity means no convection.
- **Lagged coefficients:** ∇_ad, c_P, δ (gas + equilibrium radiation, `material_thermo`) and κ_R come from the current Newton
  iterate and are held fixed during the derivative sweep. The dependence on T, ρ, P and ∇ is differentiated exactly.
  **[design choice]** Newton convergence may become linear when the coefficients change strongly between iterates; not yet observed.
- **Check** (`CONV=1 z0_jac_check.jl`, a superadiabatic test state with all 99 interfaces convecting): coloured AD = dense AD
  exactly, and vs central FD ≤ 1.5e-6 of the row max. Assembly 1.07 ms/call.
- **Validated** against the steady model's MLT on the quasi-star profile in P7.3 (agreement to ≤ 0.002 L).

### P7.3 Quasi-star outer shell: set-up and hydrostatic test [run]
- **Profile export:** `QuasiStarInfall/scripts/export_profile.jl fBH Mdot out.csv` writes the steady envelope (r, m, P, T, ρ,
  ∇, ∇_ad, ∇_rad, …). Used here: M = 10⁶ M☉, f_BH = 0.1, L = L_Edd(κ_es) = 1.42e44 erg/s, Ṁ = 0, R = 8.435e15 cm.
- **`julia/QuasiStarSetup.jl` `setup!`** builds the TITAN state from it, starting from the SUBC state as a parameter template.
  - Domain: r_in = 0.4 R to R, plus either a grey hydrostatic atmosphere (Ṁ = 0) or a free-fall inflow region (Ṁ > 0).
  - Grid: 140 envelope zones uniform in ln P_tot; the outer region uniform in r (atmosphere) or geometric from R (inflow).
  - Masses: per steradian. The mass inside r_in, m(r_in) = 1.565e39 g (BH plus inner envelope), acts as a point mass;
    the shell holds 4.2e38 g.
  - Constants: TITAN's cgrav, cc, csigr and car are overwritten with those of the Z0 physics. With Γ ≈ 1 the net force on
    the gas is about 1 % of gravity, so TITAN's 2·10⁻⁴ differences would matter.
  - Boundaries: walls (leibc = 1, leobc = 1) with an Eulerian grid (lgrid = 2). The inner boundary has an imposed luminosity
    (lribc = 3, xlum = L), all of it radiative at the wall; the outer radiation boundary is transmitting (lrobc = 1, ximr = 0).
    Inflow runs use leobc = 2: a fixed external state ρ_ff, T_Lucy, u = −v_ff, with matching mass flux (phir0) and gas
    internal-energy flux (phir2).
  - Material: `Z0Gas` (optionally `grey = true`: κ_R for absorption too), MLT on with α = 0.8 (the steady model's value).
- **MLT optical-depth cut** (new in `TitanResidFast.convcoef`): F_conv = 0 where the Rosseland optical depth from the outer
  boundary (on the current iterate) is below `TAU_MLT[]` (default 1).
  - Without it, the optically thin atmosphere convected spuriously, because there P ≈ fE_r is nearly constant and ∇ blows up.
  - It also skips unphysical iterates (ρ or T ≤ 0); the step is then rejected by the usual divergence logic.
- **Driver:** `julia/qstar_run.jl setup|relax` (environment variables QS_NAME, QS_MDOT, QS_GREY, QS_LGRID, …).
  `julia/qstar_compare.jl` compares a snapshot with the steady profile.
- **Hydrostatic test, Ṁ = 0** (2000 steps, about 0.4 dynamical times; t_dyn = √(R³/GM) ≈ 6.7e7 s).
  Table: `reference_runs/qstar_hydrostatic_compare.txt`.
  - **Interior (0.4–0.98 R):** TITAN keeps the steady structure to |ΔT/T| ≤ 0.5 % and |Δρ/ρ| ≤ 1 %. Its MLT convective
    fraction matches the steady model's 1 − ∇/∇_rad to ≤ 0.002 of L at every point (0.191 vs 0.191 at 0.42 R; 0.0165 vs
    0.0163 at 0.85 R). Total L(r) is constant to ±0.7 %. **[run]** This validates both the Z0 material and the MLT flux.
  - **Surface (r > 0.993 R):** it does not hold, and the reason is physical, not numerical. **[run + analysis]**
    - The steady model's outer boundary is the Eddington grey T(τ) with τ = 2/3 at R. Its Γ_κ = 0.91 there comes from
      partly ionized H at 7300 K.
    - Integrating the atmosphere above R hydrostatically gives only τ ≈ 0.31 in 12 density e-folds. The column needed for
      τ = 2/3 is not there.
    - TITAN's radiative equilibrium then cools the top to ≈ 5700 K (non-grey) or ≈ 5000 K (grey). H recombines; at Z = 0
      the opacity drops by orders of magnitude (the "opacity cliff"); the radiative support vanishes. The atmosphere falls at
      30–60 km/s (close to free fall) and drains onto a new photosphere near 0.995 R.
    - The Eulerian top zone then empties below the table's range (ρ ~ 1e-21) and the run stops at t ≈ 2.8e7 s (non-grey) or
      3.5e7 s (grey).
    - Implication for the quasi-star question: a Z = 0, Γ ≈ 1 envelope has a sharp photosphere pinned at the H-recombination
      edge. An extended ~6000 K atmosphere is not hydrostatically supported. This agrees with the steady branch analysis.
- **Infall runs (Ṁ = 0.3 and 3e-5 M☉/yr), first attempts: not yet working.** **[run]**
  - Eulerian grid: the run stops at t ≈ 100 s, with time steps of 0.04–5 s. The accretion shock (v_ff = 1.8e8 cm/s, Mach
    ≈ 200) sits in one zone at T ~ 1e6 K, with a 1300 K undershoot ahead of it (below the 2000 K table edge). The same
    happens at Ṁ = 3e-5, where the inflow density 1e-20 is also below the table.
  - Adaptive grid (lgrid = 1, lady = log m, ρ, T): the first step does not converge even at dt = 1e-5 s. The density
    correction at the first inflow zone is −870× while the grid equation is far from equidistribution.
  - Next to try:
    (1) Start from an equidistributed grid: relax the grid equation on the static state first, as start08 effectively does
        from dt = 1e-25.
    (2) Ramp the inflow velocity and density at the outer boundary from zero, so the shock forms gradually.
    (3) Extend the opacity table below 2000 K and below log R = −12 (cold, rarefied inflow).
    (4) Put a smooth settling layer below the initial contact.

### P7.4 Infall: next steps (1)–(3) [run]
- **(3) Extended opacity table** `julia/z0/kaptable_Z0_ext.bin.gz` (`Z0_TABLE=ext`): log T 3.0–8.0, log R −20–1, built by
  `steady/scripts/build_kap_ext.jl` (754 s on 2 threads). On the v2 nodes it is identical (max |Δ log κ_R| = 0). No molecules:
  below ~2000 K κ is a lower bound (κ_R = 1.6e-9 at ρ = 1e-20, 1500 K). Outcome: removes the table edge but does not change the
  failure (run E1 vs R1: identical failure time 2.3e4 s).
- **(2) Ramped inflow** (`QS_RAMP_FAC`, `QS_RAMP_T` in `qstar_run.jl relax`): the external density and the mass and energy
  fluxes through r_out grow as exp(t/T_ramp) up to the target. Together with a **smoothed contact**
  (`QS_WCONTACT`: ln ρ, ln T log-linear and u linear across [R, R(1+w)]), the Eulerian run takes 200–400 s steps instead of 0.04 s
  and survives to 2.3e4 s (vs ~100 s). It still fails where the infall stops: a single-zone shock at ~7e4–1e5 K with a cold
  zone ahead.
- **(1) Adaptive grid.** Root cause of the first-step failures found and fixed: TITAN's grid equation sets
  ν(ngrs+1) = ν(ngrs) and ν(ngre) = ν(ngre−1) (no time smoothing there). An initial grid that violates these two rows is
  rearranged within one Newton step at any Δt, which moves interfaces by ~1 % of r and wrecks the mass coordinate
  (δ ln m_e of order unity at the contact). `QuasiStarSetup.setup!` now makes those two pairs equal for lgrid = 1. Newton then
  converges quadratically from step 1 (dmax 6.5e-4 → 4.5e-6 → 2.7e-9 at Δt = 1e-3 s). Also added: `equidistribute` (initial grid
  of equal arc length in (r/xscale, ln ρ, ln T)), `QS_LADY`, `QS_GRIDTAU`.
  Outcome: run A1 (adaptive, τ_grid = 1e3 s) advances with Δt = 0.04–0.1 s and fails after 7.9 s at the same place.
- **Status:** the Mach-200 radiating accretion shock is not yet handled. Candidates: larger grid τ so the grid follows the shock;
  artificial-viscosity parameters; the stiffness of gas–radiation coupling behind the shock (~1 s); total-energy row.
- **Steady-model tables regenerated** (16 Oct): the browser port of the steady model showed that `steady/output/branches.csv` of
  15 Oct came from an older code (R off by up to 4.6e-3). `branches.jl`, `fold.jl`, `validate2.jl`, `run_grid.jl` were rerun;
  no root appears or disappears. The book (`report/`) uses the regenerated numbers.

### P7.5 Infall: stable dynamical runs [run]
Four defects, found one after the other by stepping the stalled run E1 (`julia/stepdiag.jl`: per step the limiting variable
and zone; `julia/newtrace.jl`: per Newton iteration the largest correction, and a scan of the residual along the Newton
direction). None of them is in TITAN's own test problems; all are triggered by thin gas in a strong radiation field.

1. **Artificial viscosity too short.** The quasi-star set-up had `ql1 = 1e-5` (viscous length 1e-5 r ≈ 0.2 zones at R).
   A planar Noh test (`julia/infall_tests.jl noh1d`, Mach 200, 100 zones, `ProblemSetup.jl` builds the state): with a
   0.2-zone length the shock lags 10–12 zones, the post-shock density is 75 % off and velocities oscillate to 0.13 v0; with
   2 zones (`ql0 = 2Δx`, `cq2 = 1`) the shock is within 0.5 zones and ρ within 3.7 %, also at Mach 2000 and on the adaptive
   grid. Infall runs now use `QS_SET=ql1=1e-4,cq2=1` (≈2 zones).
2. **Gas energy by cancellation.** `Z0Physics.gas_e` and the material's e_g were e_int − aT⁴/ρ. Behind the shock aT⁴/ρ is up
   to 1e11 times e_g, so e_g kept ~5 digits; the residual is a staircase in T with 1e-5 steps and Newton 2-cycles at
   dmax = 1.2e-5/2.0e-5 just above `conv = 1e-5` (residual scan: piecewise constant with jumps of the full residual size).
   Every such step was cut by 20, so dt stayed at 1e-3 s. Fixed by computing e_g directly (`Z0._egas`). dt then grows to 200 s.
3. **Eddington-factor source function.** eddfac used thermal absorption = χ_F ρ − σ_T n_e (TITAN's grey opacities make
   that exact). With a tabulated Rosseland χ_F that difference is table noise, and in 1e6 K gas it multiplies B ≈ 1e9 J:
   f came out 0 or 1 in alternate zones and E decoupled zone by zone. For non-built-in materials the source is now
   (κ_P ρ B + σ_T n_e J)/(κ_P ρ + σ_T n_e) (`TitanDriver.EDD_MATERIAL_SOURCE`, set by `TitanMaterials.use!`).
4. **TITAN's gas + radiation energy row.** TITAN's energy row is for e_g + E/ρ; the gas temperature follows from its
   difference with the radiation row. In the infall E/ρ ≈ 3e5 e_g, and the two rows advect e_g + E/ρ and E/ρ with separate
   van Leer limiters, so the gas energy is the small difference of two discretisation errors: zones ahead of the shock
   cooled to 400 K inside a 6000 K radiation field and the run stopped (6 dt cuts). `TitanResidFast.GASROW[] = true`
   uses the gas-only row Δ(m e_g)/dt + advection(e_g) + p∇·u + Q_visc + MLT = −(cκ_E E − 4πκ_P B)ρ dV. Algebraically equal
   and still conservative (the exchange term is the radiation row's with opposite sign). Check (`gasrow_check.jl`): SUBC
   587/587 steps, profiles equal to 7e-8; SUPC 713 vs 693 steps, median |ΔT/T| 6e-4 (different end time, 0.6 %).
   Default stays false (TITAN).

Physics added on the way:
- **Compton exchange** (`TitanResidFast.COMPTON`): Λ_C = 4σ_T n_e c E k(T − T_r)/(m_e c²), in the radiation row and the gas
  row, fully implicit. Behind a 1e8 K shock in E ≈ 10 erg/cm³ it cools in t_C ≈ 1.5e6 s, far faster than free-free.
  `compton_box.jl`: relaxation of 1e8 K gas in an optically thick box against an RK4 reference including free-free;
  max |ΔT/T| 0.32 / 0.085 / 0.027 for stol 0.1 / 0.03 / 0.01 (first order in dt, as the scheme).
- **Opacity below the table** (`Z0Physics.EXTRAP_LOWR`, default on): log κ extrapolated linearly in log R instead of
  clamped. The corona (1e7–1e8 K, 1e-17 g/cm³) lies below log R = −20; clamping overstated κ_P, and free-free cooling, by
  30–100×. The extrapolation reproduces κ_P ∝ ρ of the in-table values exactly.

**Runs** (`qstar_run.jl relax`, all with `QS_SET=ql1=1e-4,cq2=1 QS_GASROW=1 QS_COMPTON=1 QS_RAMP_FAC=1000 QS_RAMP_T=2e6`,
Eulerian grid, `Z0_TABLE=ext`; analysed with `qstar_analyze.jl`):
- G0 / G1 (outer boundary 2R, ramp ×100, before the Compton and opacity changes): G0 with TITAN's row fails at 2.5e4 s
  as E1; G1 with the gas row runs to 1e7 s in 3336 steps.
- **R2** (outer boundary 1.2R so that the ramped inflow reaches the star within the run; Ṁ from 3e-3 to 3 M⊙/yr):
  13 942 steps to 6e7 s, typical dt 1e3–3e4 s, no failure. Two regimes:
  * low Ṁ (≲0.03 M⊙/yr at the shock, t ≲ 1.5e7 s): the infall shocks to 2e7–1e8 K and does not cool (free-free time ≫
    run; Compton ~1e6 s). A hot, tenuous layer (ρ ~ 1e-17) grows above the photosphere; the shock moves out to 1.010 R;
    the photosphere stays at 0.995–0.999 R with T_eff ≈ 7300 K, L ≈ L0.
  * Ṁ ≈ 3 M⊙/yr at the shock (t ≳ 2.8e7 s): the shock collapses onto the photosphere and becomes radiative (isothermal:
    pre-shock gas at T_r ≈ 6700 K, compressed 2600–6800×), with intermittent hot episodes (T_peak up to 2.8e7 K).
    The photosphere contracts steadily: r_ph/R = 0.9945 (2.8e7 s), 0.9898 (4.3e7 s), 0.9868 (5.8e7 s), T_eff 7260–7450 K,
    L_ph = 0.98–1.07 L0. 4.4 M⊙ accreted onto the grid by 6e7 s. Not yet stationary.
  * For comparison the steady model's large-R root at Ṁ = 3 (f_BH = 0.1) has R = 0.9876 R0, T_eff = 7359 K: the run is
    passing through it, not settling on it (yet). [interpretation: the contraction is the outer envelope adjusting to the
    ram pressure on a time longer than the run; continuation R2c to 2.5e8 s is running] — superseded: see P7.6, the
    contraction is not caused by the missing sink or the time step, and the adaptive grid does not show it.
- R2L (half resolution, 70 + 40 zones): same sequence; at Ṁ ≈ 3 r_ph/R = 0.9937 (3.6e7 s), 0.9891 (5.4e7 s), T_eff 7270–7400 K.
- R2A (adaptive grid): still limited to dt ≈ 0.045 s by ρ changes of 10 % per step in the zone at the envelope–infall
  contact as the grid moves through it; stopped. Open.

### P7.6 Mass sink, second-order time stepping, adaptive grid [run]
**Inner mass sink** (`TitanResidFast.INNER_SINK`, `qstar_run.jl` `QS_SINK_TAU`). The inner wall becomes porous to mass:
phil0 < 0 carries mass onto the central point mass, with the gas enthalpy phil0·(e_g + p_g/ρ) in phil2 and the
radiation enthalpy (1 + f)E/ρ·phil0 in the radiation row. A controller holds the mass below the shock at its initial
value: Ṁ_sink = max(0, M_s − M_s0)/τ. The inner wall itself stays at rest.
- **TITAN 1.2 bug found:** `mass.f` books the inner boundary mass flux as `+ rmu·phil0·dt`, the opposite sign to
  `contin.f` (inflow rmu·phil0 into zone ngrs) and to `peq` (dmdt(ngrs) = −r^mu·phil0). Never exercised (phil0 = 0 in
  every TITAN test problem). With mass leaving, m(r) lost what the core should gain (first sink run: "point mass gained
  −8.5 M⊙"), and every advective flux derived from dm_e/dt was off. Fixed (`FIX_PHIL0`, default on; SUBC unchanged,
  587/587 steps).
- Runs S6 (τ = 1e6 s) and S7 (τ = 1e7 s), otherwise as R2 to 1.2e8 s: the excess mass below the shock levels off at
  1.2 and 2.0 M⊙ (R2/R2c: ~14.5 M⊙ by 1.66e8 s); 8.8 / 7.9 M⊙ moved onto the point mass. **The photosphere contracts exactly as
  without the sink** (r_ph/R: 0.9868 at 5.7e7 s and 0.9841 at 1.16e8 s in S6; 0.9868 and 0.9841 in R2/R2c). The
  accumulated mass was therefore not the cause.

**Second-order time stepping** (`julia/Trbdf2.jl`, `QS_TRBDF2=1`): TR-BDF2 built from TITAN's θ-scheme. Stage 1: θ = ½
over γΔ (γ = 2 − √2). Stage 2: BDF2, written as backward Euler over (1−γ)/(2−γ)Δ from y* = A y_γ − B y_n, with y* formed
on the conserved quantities the residual reads (m, m e_g, m(e_g + E/ρ), E dV, mass-weighted momentum, the flux-row
inertia F·2M/(ρ₋ + ρ₊), exterior mass). Eulerian grid only. Needed on the way: the Compton term was fully implicit,
which made even θ = ½ first order; it is now θ-weighted like the other exchange terms.
- `order_test.jl` (fixed steps, error against N = 4096 or 16384; outputs stored in `reference_runs/order_test/`): Compton box: TR-BDF2 2.00, θ = ½ 2.00, backward Euler
  1.00 (errors at N = 128: 2.5e-5, 6.2e-5, 1.6e-2). Standing sound wave, 64 zones, amplitude 1e-3 c_s: TR-BDF2 3.0 → 2.3
  down to an absolute floor ≈ 1e-9 c_s (not identified; Newton tolerance and van Leer slopes ruled out); backward Euler
  0.76 → 0.97. Radiating wave (grey, τ = 1 per zone): pre-asymptotic up to N = 2048: TR-BDF2 1.4–1.7, θ = ½ 1.7–1.85, backward
  Euler 0.3–0.7; at N = 2048 TR-BDF2 is 80× more accurate than backward Euler.
- The implicit-midpoint stage fails Newton at an under-resolved shock where backward Euler converges. Fallback: the same
  Δ is retried with backward Euler, and the next 10 steps stay with it (`FALLBACK_BE`, `STICKY`). Those steps are first
  order.
- Run T2 (as S6, TR-BDF2, to 4e7 s): r_ph on the same zone edge as S6 at every 5e6 s (i.e. equal to the zone width,
  ~1e-3 R near the photosphere).
- **Courant cap** (`QS_CFL`): run K1 (as S6, Δt ≤ min Δr/|u|, about 2.5e3 s instead of up to 3e5 s): r_ph(t) within one zone
  (≤ 9.6e-4 R) of S6 up to 9e7 s; S6 has no snapshot between 7.6e7 and 1.16e8 s (K1 failed at 9.5e7 s; not investigated). On the Eulerian grid the result does not
  depend on the time step.

**Adaptive grid.** Diagnosis of the 0.05 s steps: the zones sweep across the envelope–infall contact, or later across
the single-zone radiative shock (a 600× density jump), so their contents change by 10 % per step whatever Δt. With the
grid frozen (τ_grid = 1e30), the adaptive code reproduces the Eulerian steps exactly, so the lgrid = 1 path is
consistent; the problem is how fast the grid moves. With τ_grid = 1e3 s the grid equation's imbalance at the jump
(x_ν² S ~ 1e6) moves interfaces at up to 1e10 cm/s. With τ_grid = 1e7 s and weights (exterior mass 1, ln ρ 1e-2, ln T 0)
the grid concentrates zones gradually and the run proceeds (AD1: median step ~830 s, 10 % of steps < 15 s, maximum
  2.4e5 s). The diagnosis numbers above (0.05 s steps, 10 %/step, 1e10 cm/s, frozen grid) are from interactive
  `stepdiag.jl` sessions, not stored runs.
- `to_adaptive.jl` switches a saved Eulerian state to lgrid = 1 (first and last zone pairs re-split to equal width,
  conserving mass and energies, keeping the profile shape: equal values in the two zones made the MLT gradient 0/0).
- **AD1** (R2 at 6e7 s → adaptive, to 1.2e8 s, 15 150 steps) and **AD2** (same, ln ρ weight 0.1): the grid shrinks the
  zones to 0.4–1.5e-4 R at the velocity jump and to ~1e-6 R (minimum 3.4e-7 R) at the bottom of the cooling layer, and
  resolves **a hot post-shock cooling layer**: 3.6e7 K right behind the shock (the
  strong-shock value 3/16 μm_H v²/k ≈ 3.8e7 K), cooling to 6900 K over ~1.3e-3 R. The Eulerian zones there are 8e-4 R;
  the Eulerian runs collapse this layer into one cool zone. [estimate: Compton cooling time behind the shock ≈ 2e6 s,
  times v/4 ≈ 4e7 cm/s gives ~1e-2 R; the layer is shorter because the gas compresses and free–free cooling grows as it
  cools]
- With the layer resolved **the photosphere turns round**: right after the switch it contracts faster than the Eulerian
  run (min 0.9824 at 7.28e7 s in AD1, 0.9826 at 6.74e7 s in AD2), then rises from ~7.5e7 s to 0.9894–0.9898
  (1.07–1.2e8 s), with L_ph = 0.91–0.92 L0 and T_eff = 7220–7240 K. Over 1.07–1.2e8 s AD1 and AD2 agree to 1e-4 R,
  0.004 L0, 8 K (from 6.5e7 s on: up to 9e-4 R, 0.046 L0, 85 K); their peak temperatures differ (3.4e7 vs 1.4–1.8e7 K, the jump itself is resolved
  differently). The Eulerian runs at the same times: r_ph 0.984, L_ph 0.97–1.02 L0 (R2c, S6).
- **R2A** (adaptive from the set-up: `QS_LGRID=1 QS_GRIDTAU=1e7 QS_LADY=2,2,0 QS_GRIDW=1,1e-2,0`, otherwise R2; 18 295
  steps to 1.22e8 s): r_ph/R has a minimum of 0.9873 at 6.6e7 s and is back at 0.9910 at 1.22e8 s (T_eff 7244 K,
  L_ph 0.97 L0). AD1 at the same time: 0.9894, 0.92 L0, T_eff 7243 K; Eulerian S6: 0.9831. The adaptive runs agree on
  the behaviour (minimum, then re-expansion; its depth depends on the start) but at 1.2e8 s still remember their start.
  The steady root at Ṁ = 3 has L = 1.02 L0 (R²T_eff⁴, accretion luminosity included), so the adaptive deficit is ~10 %.
- **R2AH** (adaptive, 200 + 114 zones): stalled at 2.53e7 s (Δt → 1e-3 s) with the MLT flux of one zone just below the
  shock (1.2e4 K, Δr = 7e-5 R) oscillating up to 2 L0. MLT coefficients are lagged (from the Newton iterate, not in the
  Jacobian), the likely cause [interpretation]. Stopped; convergence in zone count still open.
- Verification: an independent pass checked every number of this section and of chapter 7 against the snapshots
  (11 corrections applied, 17 Oct).
- **Conclusion so far:** the late, high-Ṁ evolution of R2 is set by an under-resolved shock. The adaptive runs are the
  better approximation; they are not yet a converged answer (one zone count, mgr = 330 caps the number of zones).
