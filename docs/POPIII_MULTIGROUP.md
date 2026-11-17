# Implicit multigroup radiation for the Pop III protostar (`MGROUPS`)

Branch `multigroup` (from `popiii-s6` 43269f6), 16 Nov 2026. Code: julia/MGRad.jl (included into TitanResidFast: unknowns, rows,
Newton hooks, Planck group fractions), julia/MGOpac.jl (included into PopIIIMaterial: group-opacity ratios, group photo-rates),
julia/G2Opacity.jl (`group_fast`, `group_index`), hooks in TitanResidFast.residual!, Trbdf2.jl, TitanFortran.jl (meqn), popiii_m2c.jl
(`mg_install!`, checkpoint, `mg_record`). Unit test julia/mg_test.jl. Test runs `D/t_mg1N`, D = /zpool/nvme/data/tabel_scratch/titan_popiii/runs,
from `git archive` trees in /zpool/nvme/data/tabel_scratch/titan_popiii/mg_trees, restarts of copies of z100mm's checkpoint.

Status of this note: §1–§9 are the design; §10 is what the first coupled prototype (stage 1) implements and what it measured; §11 the
plan of the remaining stages. Where a section says "stage 1", the prototype does exactly that; everything else is design.

Current state (17 Nov 2026, branch `mg2b`): §15–§17 Rayleigh/Raman scattering and stage 2b; §18 stage 2c (photo-rates and photo-heating
from the group fields inside the solve); §19 the switch audit: `MGROUPS` is the only switch, and the switches that §5–§14 name
(`MGVEF`, `MGSOLVER`, `MGMILNE`, `MGN2BF`, `MGMILNET`, `TIESMOOTH`) are history; §20 validation of the final configuration (bce234e): from a grey checkpoint after protostar formation the multigroup
run is as stable as grey and adds no energy error of its own; two gates are **not** met — the IonFS comparison (§20.3) and the run
from the initial conditions, which does not reach the protostar (§20.6).

## 1. Why, and what multigroup has to deliver

The run carries one radiation field (E_r, F_r) with the grey means κ_E (absorption), κ_P (emission) and χ (flux mean), and the
network's photo-processes are fed by patches that reconstruct a spectrum from that one field: STARLIGHT (J_ν = W B_ν(T_c), one colour
temperature and one dilution for every process), SLBLEND/SLEB (T_c from T_eff or T_r at τ + 2/3), SLWFLUX (W from the flux), LOCALJ (the
network's own photo-processes in that J_ν), SLIMP (W, T_c at the Newton iterate), SLSHIELD (one LW frequency, 12.87 eV, times a column
shielding fit), LYABB. Since the G2 material (docs/POPIII_G2_MATERIAL.md) the photosphere and envelope are network gas whose continuum
opacity is computed per frequency (G2OPAC=2, 484 nodes, exact AD), and the table is confined to the non-ideal interior: the grey means
are now a reduction of frequency-resolved information the code already has.

The photosphere is a transparency front about one photon mean free path thick: H⁻ and H bound–free opacity rise steeply with T, so
the IR (H⁻ free–free), the H⁻ bound–free band and the Balmer continuum leave from different depths. A grey field must put all of
them at one depth (the grey cliff of POPIII_SLABFLUX.md: one interface carries χ ratios of 1e3–1e5). Multigroup is the physical cure
of that cliff, and it is what the photo-rates need: H⁻ photodetachment (≥ 0.755 eV), H(n = 2) photoionisation (≥ 3.40 eV),
Lyman–Werner dissociation (11.2–13.6 eV), H, He, He⁺ photoionisation (≥ 13.6, 24.6, 54.4 eV) each see their own band, transported
through the absorbers that matter for that band, photon-conserving.

Requirements (the stack's, docs/POPIII_STACK.md): fully implicit (the radiation is a Newton unknown with exact derivatives, as the grey
field is now), TR-BDF2 with its error estimate, ~100 zones, stable over very long runs; bit-identical when off.

## 2. Group structure

Edges are the code's own thresholds (MGOpac.jl `mg_thresholds`, all from constants in the code, none from memory):

| edge | value [eV] | source in the code | why an edge here |
|---|---|---|---|
| H⁻ photodetachment threshold | 0.7551 | λ₀ = 1.6419 μm, ContinuumOpacity.σ_Hm_bf (John 1988, 1988A&A...193..189J; zero for λ ≥ λ₀) | below it only free–free (H⁻ ff, H⁺ ff), CIA, Rayleigh: the IR group; above it H⁻ bound–free dominates the photosphere |
| H n = 3 (Paschen) edge | 1.511 | χ_H/9, Z0Physics.χH | (8-group set) H(n=3) bound–free in thick layers |
| H n = 2 (Balmer) edge | 3.400 | χ_H/4 | H(n = 2) photoionisation of the three-level atom (β), Balmer continuum |
| Lyman–Werner band, lower end | 11.2 | the RADGROUPS trial's band (912–1107 Å, Draine & Bertoldi 1996, 1996ApJ...468..269D) | H₂ dissociation by line absorption (the code's rate is normalised at 12.87 eV, Abel et al. 1997, 1997NewA....2..181A) |
| H ionisation | 13.598 | Z0Physics.χH (the G2 continuum's threshold; ChemNet.χH = 13.598) | H(1s) bound–free |
| He I ionisation | 24.587 | Z0Physics.χHe1 = ChemNet.χHe1 (NIST ASD) | He I bound–free (Verner et al. 1996 in G2) |
| He II ionisation | 54.418 | Z0Physics.χHe2 = ChemNet.χHe2 | He⁺ bound–free |

Default sets (`MGROUPS=n`): 3 = IR | optical | ≥ 11.2 eV; **4 = IR | optical (H⁻ band, Balmer, Lyα) | LW | ionising** (stage 1);
6 = IR | H⁻ (0.755–3.40) | Balmer + Lyα (3.40–11.2) | LW | H I (13.6–24.6) | He (≥ 24.6); 7 = 6 with He I | He II split at 54.4;
8 = 7 with the H⁻ band split at the Paschen edge. `MGROUPS=e1,e2,…` gives any edge list in eV.

Lyα (10.2 eV = 3χ_H/4 from the code's E21) sits in the Balmer/optical group: it is pumped and scattered there and its energy comes
from the n = 2 cycle (LYABB, the three-level atom); a separate Lyα group needs line transfer (LyaTransfer, not grey) and is not a
continuum group. The 11.2 eV edge is a convention for the band, not a threshold of the continuum: the continuum between 11.2 and 13.6 eV
is small in network gas (no H⁻, which vanishes below 0.125 μm in σ_Hm_bf; H(n ≥ 2) bound–free, free–free, Thomson, capped Rayleigh), so
the LW group transports nearly freely until H₂ absorbs it (§6).

Group count: 6–8 groups resolve every process; 4 is the smallest set that separates the processes of the present stack (H⁻ and
three-level in one group since they use the same W in STARLIGHT now). The in-group spectral shape is the one assumption left: it is
B_ν(T_s) with T_s the local colour temperature (§4), so narrow groups at the thresholds are what makes the answer insensitive to it.

## 3. Equations per group

The grey rows of TitanResidFast (comoving frame, frequency-integrated moment equations, Mihalas & Mihalas 1984, 1984frh..book.....M;
TITAN's discretisation, Gehmeyr & Mihalas 1994, 1994PhyD...77..320G) are integrated over each group — the multigroup moment form of
Shestakov & Offner (2008, 2008JCoPh.227.2154S), Vaytet et al. (2011, 2011JQSRT.112.1323V), Zhang et al. (2013, 2013ApJS..204....7Z),
González et al. (2015, 2015A&A...578A..12G) — closed per group by a variable Eddington factor (Stone, Mihalas & Norman 1992,
1992ApJS...80..819S). Per steradian, zone k (volume ΔV_k), interface k, group g, θ-level quantities as the grey rows (MGRad `mg_rows!`):

* **energy** (zone): Δ(E_g ΔV)/Δt + r²F_g|_k^{k+1} + f_g E_g Δ(r²u) + (1 − 3f_g) E_g ½μ(u/r) ΔV
  + x_cap,g (c κ_E,g E_g − 4π κ_P,g B_g(T)) ρ ΔV + Q_g ΔV + advection of E_g/ρ with the step's mass fluxes = 0;
* **momentum** (interface): inertia (F_g/ρ̄ c²) + r²[q f_g E_g]_{k−1}^{k} (Auer sphericality form, VEF_SPH) + ½ r² (F_g/c²) Δu
  + ΔM χ̄_g F_g/c = 0, blended with the FLUXTIE closure F_g = h_g c E_g (donor-cell) where the group's depth from the outer boundary is small;
* **outer boundary**: F_g = c g_g E_g,s − π(2g_g + 1) I_cmb b_g(T_cmb) (E_g,s from the zone centre by the sphericality factor over the half zone);
* **inner boundary**: F_g = 0 at the regular centre; ghost rows as the grey ones.

B_g(T) = b_g(T) σT⁴/π, b_g the Planck fraction of the group: b_g = P(x_hi) − P(x_lo), P(x) = (15/π⁴)∫₀ˣ t³/(eᵗ − 1)dt, with the
Bernoulli series for x < 1 and the tail Σₙ e^{−nx}(x³/n + 3x²/n² + 6x/n³ + 6/n⁴) for x ≥ 1 (both differences of tails where the lower
edge has x ≥ 1, no cancellation), generic in the number type (MGRad `mg_bfrac`, the series of the RADGROUPS trial). Unit test: Σ b_g = 1
to 1e-16, against direct quadrature 2e-10, AD derivative against a 4th-order difference 1.6e-10. b_g is floored at 1e-200 (an energy of
1e-200 aT⁴: none) so that ln E_g stays finite where a group's Planck fraction underflows (the ionising groups at T ≲ 300 K).

**What the hydro sees.** The grey rows become E_r = Σ_g E_g (per zone) and F_r = Σ_g F_g (per interface and at the outer boundary):
algebraic rows that keep E_r, F_r as unknowns. Every term that reads E_r, F_r is then the group sum without change: the total-energy row
and the gas + radiation row (their flux r²F_r, the work f E_r ∇·u and (1 − 3f)E_r u/r), the radiation pressure f E_r of MLT and the
dual-energy switch, the formal solution (f, g, h), the Δt control. The gas–radiation exchange is not in those rows (it is internal to
gas + radiation), so the gas temperature follows from the total-energy row minus Σ_g of the group energy rows: energy is conserved by
construction, whatever the group opacities are. The one term that needs the groups is the **radiation force** on the gas,
Σ_g χ̄_g F_g/c in the momentum row instead of χ̄ F_r/c, and its work in the total-energy row needs nothing (that row has no force term:
the force is internal between gas and radiation momentum). With a group-dependent f_g (stage 2) the work terms become Σ_g f_g E_g.

Left out per group (as the RADGROUPS trial): the Doppler terms that move photons across group edges in compressing or expanding flow,
O(u/c) ≲ 3e-4 here (u ≲ 100 km/s, free fall onto the star); Vaytet et al. (2011) keep them. They matter only for narrow groups next to a strong edge in a fast flow.

## 4. Group opacities and emission from the network populations

*(Evaluation level: stages 1 and 2a evaluated the ratios below at the start of the step. Since stage 2b they are evaluated at the Newton iterate and the lagged path is deleted; §17. Lines per group: §17.)*

**Network zones** (the photosphere, the envelope, the precursor): the G2 continuum on the zone's own populations
(H, H⁻ from the network, n = 2 from the three-level atom, H₂, He, He⁺, He²⁺, e⁻; HM88 occupation probabilities; Gaunt factors;
dissolved levels; docs/POPIII_G2_MATERIAL.md §6.2) on the 484 nodes of G2OPAC=2, reduced per group with the weights of the grey means:

* κ_P,g = ∫_g α_P B_ν(T) / ∫_g B_ν(T) — emission (H⁻ in Saha with the zone's H and e⁻: the Kirchhoff–Milne emission of NLOPAC);
* κ_E,g = ∫_g α_E B_ν(T_w) / ∫_g B_ν(T_w) — absorption (the network's own H⁻), in-group spectrum of the absorbed light (T_w = T_c of the
  star light, NLOPACJ);
* χ_g = ∫_g ∂B_ν/∂T(T) / ∫_g α_ext⁻¹ ∂B_ν/∂T(T) — extinction (H bound–free and scattering included), the group's Rosseland-type flux mean.

With these weights Σ_g b_g(T) κ_P,g = κ_P, Σ_g b_g(T_w) κ_E,g = κ_E, Σ_g b'_g(T)/χ_g = 1/χ (b'_g the group fraction of ∂B/∂T) up to the
node quadrature: unit test 2 (four states from 1500 K molecular to 1.2e4 K ionised gas) 1.000–0.999 (κ_P, κ_E) and 1.0000–1.0071 (χ),
and the grey totals of `group_fast` equal `means_fast` bit for bit. The group means differ from the grey ones by large factors, which is
the point: in 5000 K photosphere gas κ_P,g/κ_P = 0.008 (IR) | 1.16 (optical) | 4e-9 (LW) | 1e-4 (ionising, without H(1s)), and the
H(1s) absorption of the ionising group is 2.5e6 cm² g⁻¹ against a grey κ_E of 0.11.

The in-group weights are those of the grey means at the gas T (emission, extinction) and T_w (absorption), so the identities hold at
every T; a group whose weight underflows (the ionising groups below ~300 K) is taken as grey (its Planck fraction is at the floor anyway).
Only the H(1s) absorption of the ionising groups uses the shape B_ν(max(T_w, 3000 K)) (MGTMIN), the star's light. (Stage 1 first floored
all in-group shapes at 3000 K: that broke Σ_g b_g κ_P,g = κ_P below 3000 K, and the 200–1100 K gas of the outer cloud over-emitted into
the optical group — 4.5 L☉ through the outer boundary after 400 steps with the grey f, 5e8 L☉ within 50 steps with per-group f (t_mg27).
Fixed in G2Opacity.group_fast; the stage-1 numbers of §10 are from before this fix, with the outer optical flux still ≤ 4.5 L☉.)

**Design: exact group opacities at the iterate.** κ_X,g(ρ, T, x) directly from the per-node kernels with Duals (G2OPAC=2 already does
this for the grey means at ~50 µs per call with a 2-partial Jacobian; a group reduction adds G accumulators to the same node loop).
**Stage 1 (prototype):** the run's grey means at the Newton iterate (exact AD, including the NLOIT LTE line calibration) times lagged
group ratios r_X,g = κ_X,g/κ_X from the start-of-step composition, κ_X,g = κ_X(iterate) r_X,g (PopIIIMaterial `mg_ratios!`, ~6 ms per step
at 100 zones). The ratios keep the run's calibration (lines in κ_R above 6300 K through NLOIT's LTE factor, OPAL/table where used) and
distribute it over the groups like the continuum; the T-dependence of the ratio is lagged by one step (Jacobian inexact in that piece only).

**Lines.** The table's κ_P above 6300 K is an LTE line Planck mean ×30–2000 over the continuum (G2 doc §7.3). With groups it would go
where the lines are: Lyα and Balmer lines in the optical/Balmer group, the Lyman series in the LW group. Stage 1 distributes it with the
continuum ratio (a limitation; most relevant in the 6300 K–5e4 K layer below the photosphere, which is thick in every group).

**Scattering.** χ_g includes Thomson and Rayleigh (Rayleigh capped at its 1300 Å value as in G2); κ_E, κ_P do not. The formal solution
(EDDSCAT) treats χ − κ_P as scattering; per group the same with χ_g, κ_P,g (stage 2).

**Table interior (non-ideal zones, w = 1 in G2BLEND; the HOT branch):** r_X,g = 1 — the grey field split by Planck fractions,
E_g → b_g(T) aT⁴ and F_g → −(c/3χρ)∂E_g/∂r = b'_g F_r in the diffusion limit: the groups sum exactly to TITAN's grey diffusion there
(the asymptotic diffusion limit of Larsen, Morel & Miller 1987, 1987JCoPh..69..283L, per group). The table has no frequency information
(OPAL κ_R, Z0 line κ_P); in the interior the radiation is thermalised and diffusive in every group, so nothing is lost.

## 5. Moment closure per group: VEF and the flux tie (SLABFLUX/SLABTIE)

* **Interface opacity:** χ̄_g = SLABFLUX's χ_n(χ_g,k−1, χ_g,k) (`sf_mean`, n = ∞ default; SUBZONE on top if set), else the harmonic mean
  — the slab argument of POPIII_SLABFLUX.md §2 holds per group (F_g constant across the dual cell, the exchange in the zone rows).
  The grey cliff can then split into G fronts at different interfaces (each group's χ_g ratio differs by up to 1e5 across the photosphere
  zone, e.g. χ_LW/χ = 0.14 below and 90–7e5 above it in t_mg14): the multigroup cure of the cliff; whether the group photospheres resolve
  it needs the per-group closure of stage 2a. **Stage 1: implemented** (`mg_theta!`).
* **Flux tie:** F_g = h_g c E_g with SLABTIE's donor h at the group's slab depth, `st_hdonor(h_FS, f, χ_g,k−1, χ_g,k, M₋, M₊)`.
  **Stage 1:** the grey formal solution's h and f with the group's χ_g in SLABTIE (the donor depth is per group).
* **Eddington factors f_g, g_g, h_g from a per-group formal solution** (TitanCore.eddfac! with the group's χ_g, κ_P,g, source
  (κ_P,g B_g + (χ_g − κ_P,g) c E_g/4π)/χ_g, boundary intensity I_cmb b_g): stage 2. Cost: G formal solutions per Newton iterate
  (VEFIMP). Where it matters: the thin layers where different groups have different photospheres (f_g → 1 in a group whose photosphere
  lies inside, while the grey f of the zone is still < 1/2).

**The two Newton-convergence fixes from the slab-flux work** (POPIII_SLABFLUX.md §4.5: 20–22 failures/100 steps with SLABFLUX+SLABTIE,
slow, not divergent, largest corrections in E_r and the front zone's species):

1. **Smooth the max(·, 0) kinks of the tie.** The donor-cell tie uses max(h_k−1, 0) E_k−1 + min(h_k, 0) E_k and SLABTIE's
   h = max(h_FS + h_slab(τ_i) − h_slab(a), 0). Replace max(h, 0) by the C¹ ramp s(h) = ½(h + √(h² + ε²)) − ε/2 with ε = 1e-3 (h is
   O(1/3) where it matters, O(1e-5) only deep in the photosphere where the tie weight is 0): s(0) = 0 (no inward light from the
   upwind side), ds/dh continuous through 0; the error against max is ≤ ε/2 · (ε/|h|) for |h| ≫ ε. Same in every group. A switch,
   evaluated per group, so it can be tested on the grey run first (`TIESMOOTH=ε`, stage 2).
2. **Exact derivatives of f, g, h.** Now f, g, h come from the formal solution of the iterate (VEFIMP) but are frozen inside the Jacobian
   (∂f/∂x is dense and was dropped: docs/PORTING_NOTES.md "Dropped term 1"), so the Newton contracts linearly where f depends strongly on
   the state (the photosphere: f from 1/3 to 1 across two zones). Two options:
   * **Jacobian-free Newton–Krylov** (Knoll & Keyes 2004, 2004JCoPh.193..357K): GMRES on J δ = −R with J v from a difference of the full
     residual including a fresh formal solution, J v ≈ [R(x + εv) − R(x)]/ε, preconditioned by the present banded AD Jacobian (f frozen).
     The preconditioned operator is I + (frozen-f error), whose spectrum clusters at 1 except for the few photosphere modes, so 2–5 Krylov
     iterations per Newton step are expected. Cost per Krylov iteration: one Float64 residual + one formal solution (at 100 zones a Float64
     residual is ~1/50–1/100 of the AD Jacobian assembly, the formal solution a few ms), i.e. ~10–20 % of a Newton iteration's present cost,
     and the convergence becomes quadratic again. The banded LU is reused (only its solve per Krylov step). Matches the code's structure:
     the matslv hook takes the preconditioned GMRES, everything else unchanged.
   * **Truncated local finite-difference derivatives of f, g, h** (∂f_j/∂y_k for |j − k| ≤ 2, by a coloured perturbation of y_k = χ_k, κ_P,k,
     E_k, T_k in groups of every 5th zone, 5 colours × 4 inputs = 20 formal solutions per Jacobian). This keeps the banded direct solve,
     but the truncation drops exactly the non-local part (f at the photosphere depends on the whole thin envelope through the rays), and the
     finite differences of a piecewise formal solution (limiter in `slopes!`, the μ-quadrature) are noisy.
   **Recommendation: JFNK with the banded Jacobian as preconditioner** (exact in the limit, no truncation, cheaper than 20 formal solutions
   per iteration). With groups it is the only affordable route to exact f_g, g_g, h_g (G formal solutions per J v product, against 20·G).
   Test plan: grey run first (z100mm restart, SLABFLUX + SLABTIE): Newton iterations per step and failures/100 against the 20–22 now.

## 6. Gas–radiation exchange per group, photo-rates and photo-heating from the group fields

**Exchange.** Per group, x_cap,g (c κ_E,g E_g − 4π κ_P,g B_g) ρ ΔV in the group's energy row (EXCAP per group: the factor
min(1, EXCAP/(c κ_E,g ρ Δt)) of the grey rows with the group's κ_E,g; same factor in absorption and emission, so the group conserves energy).
The gas gets −Σ_g of it through the total-energy row (§3).

**Photo-rates (stage 1).** In every lit network zone (starlight_setup `rates(k)`), each group has its own dilution

  W_g = c u_g / (4σT_c⁴ b_g(T_c)),  c u_g = c(E_g − E_cmb b_g(T_cmb))  (SLWFLUX: at least ½(F_g,k⁺ + F_g,k+1⁺)),

capped at 1 where the grey W is (SLBLEND), 1 in thermalised zones (τ ≥ 3), and the network sees J_ν = W_g B_ν(T_c) for ν in group g
(`mg_J`; T_c is now only the in-group shape). With W_g = W for every g this is today's field exactly. Rates: H⁻ photodetachment
(ChemNet.photorate of the piecewise J), the three-level β = W_Balmer β_2s(T_c) and P12 from the Lyα group's W, LW
k = 1.1e8 · 4π W_LW B_ν(12.87 eV, T_c) f_sh (the code's rate), and **new: H(1s) photoionisation**
k_1s = Σ_{g ≥ χ_H} W_g 4π ∫_g σ_1s(ν) B_ν(T_c)/(hν) dν (σ_1s the G2 cross section with its Gaunt factor) in lit, not deep, zones.
SLIMP: the iterate factor q = W(iterate)/W_start is taken from the H⁻-band group (`sl_qt_mg`: its E_g, F_g at the iterate), and scales
every star-light term of the coupled zone (one q per zone in this stage).

**Photo-heating** is not a separate term: the energy the group fields lose by absorption (c κ_E,g E_g, including the H(1s) part) goes to
the gas through the total-energy row, and the chemistry takes the ionisation/detachment energy from u (ChemNet's u includes the binding
energies). Photon conservation then requires that the absorption that heats is the absorption that ionises: for H(1s) the group's κ_1s
(absorption weighted with B_ν(T_c) over the same nodes, `group_fast` κ1s) and k_1s use the same cross section and in-group shape, so
c ρ κ_1s,g E_g = n_1s k_1s,g ⟨hν⟩_g in the steady state. Case B is kept: the network has no ground-state recombination, so the ionising
groups have no H(1s) recombination emission (on-the-spot), and H(1s) absorption is added only where H bound–free is not already in κ_E
(not in the deep zones of NLOHBFDEEP).

**What retires (design target, after stage 2):** STARLIGHT's single W and T_c for the photo-processes (W_g from the groups; T_c as the
in-group shape only, which narrow groups make irrelevant), SLWFLUX (the group flux is the group's dilution), LOCALJ (the network's own
photo-processes with the piecewise J_ν — stage 1 still uses LOCALJ's tabulated rates with the Balmer group's W, exact for 4 groups where
the H⁻ band, Balmer continuum and Lyα share one group, approximate for 6–8 groups), SLIMP (the group fields are Newton unknowns: the
rates become functions of E_g at the iterate, one factor per process instead of one per zone), SLSHIELD (below).

**LW self-shielding, photon-conserving (design).** The code's rate is k = k₀ f_sh(N_H₂, T) f_HI(N_HI) with the Wolcott-Green et al. (2011,
2011MNRAS.418..838W, eq. 10, α = 1.1) and Wolcott-Green & Haiman (2011, 2011MNRAS.412.2603W, eq. 13) fits on Sobolev-limited columns.
That is a column fit applied to a field that does not lose the photons. With a LW group the photons must be removed where they are
absorbed. The band-mean H₂ absorption opacity that reproduces the fit's attenuation along a column is
σ_eff(N) = −d ln f_sh/dN per H₂ (the local logarithmic slope of the shielding function: a beam through dN loses the fraction σ_eff dN), so

  κ_LW,H₂ = n_H₂ σ_eff(N_loc)/ρ   in the LW group's κ_E (absorption, no emission: the absorbed photons fluoresce mostly back into the
  band and the IR, a scattering-like part kept in the group, and a fraction p_d dissociates),

with N_loc the zone's own Sobolev column (the line saturation inside the zone; the attenuation over the flow then comes from the transport).
The dissociation rate per H₂ is then k = p_d c σ_eff E_LW / ⟨hν⟩_LW from the group field at the iterate, and f_sh is not applied a second
time. p_d follows from the thin limit: k(N → 0) must equal the code's unshielded k₀ = 1.1e8 · 4π J_ν(12.87 eV), which fixes p_d σ_eff(0) per
unit band energy density (p_d itself — 0.15 in the code's IonFS post-process, ionfs_post.jl `PD`, after Draine & Bertoldi 1996 — only sets
how much energy is fluorescence vs dissociation). The gas heat
per dissociation stays the code's 4.48 + 0.4 eV (lwheat!), now taken from the LW group (stage 1: from the optical group, with GASQ).
This is photon-conserving by construction and needs no column integration in the chemistry (the columns become local).

**H₂ line cooling** (H2Lines ALI/EP, LINEIMP at the iterate) stays as it is: the H₂ rovibrational lines are in the IR (≲ 1 eV) and are
transported by their own line transfer; their cooling enters the gas energy (cpL), not a group. Interplay: the line photons that escape
leave the domain; those absorbed by the continuum would belong to the IR group — the ALI does not hand them over (as now in the grey run).

## 7. Newton structure and cost

Unknowns per zone: r, x_me, ρ, u, T, E_r, F_r, x_H⁺, x_H₂ (9) + E_g, F_g (2G): neqn = 9 + 2G (13 for G = 2 … 25 for G = 8), meqn ≥ neqn
(TitanFortran sets it from MGROUPS). E_g are logarithmic variables (E_g = E_g,0 e^δ, positive; correction clamped to ±3 per iteration), F_g
linear with a per-interface nominal scale max(|F_g|, 1e-3 c Ē_g). The convergence measure includes |δ ln E_g| weighted by
E_g/(E_g + 1e-3 E_r) (a group that carries < 1e-3 of the field does not hold the Newton), the TR-BDF2 error norm likewise (with the grey
E floor). TR-BDF2's y* combination, save/restore and the filtered estimate carry E_g ΔV and the flux-row inertia of F_g like E_r, F_r
(Trbdf2.jl). The Jacobian is the analytic AD of the residual (5-colour stencil, ForwardDiff): the group rows are linear in E_g, F_g at
fixed (ρ, T, r, u), so their part is exact; the lagged ratios and the frozen f are the only inexact pieces.

Cost per Newton iteration (100 zones) has three parts: (a) the AD assembly — the Dual width is ncol × neqn (6 × 9 = 54 grey, 6 × 17 = 102
for G = 4, 6 × 25 = 150 for G = 8): the material sweep (EOS, NLOIT opacity, coupled chemistry lifted from cached local Jacobians) carries
all partials, so its cost grows ~linearly with neqn although the material does not depend on E_g (except through q); (b) the linear solve:
banded LU of n = neqn(N + 4) with bandwidth 3 neqn − 1, O(n·neqn²) = O(N neqn³), 1768 × 50 for G = 4 (pivoted LAPACK gbtrf, ms);
(c) the formal solution(s). Measured for stage 1 in §10. **Design for the cost:** colour the Jacobian in two passes — the hydro/material
part with Duals over the 9 base unknowns (as now) and the group rows (which read the material only through χ, κ_E, κ_P, B and T) with
Duals over (ρ, T, x, E_g, F_g of the stencil) — or write the group rows' Jacobian analytically (they are linear in E_g, F_g); either keeps
the expensive material sweep at 54 partials. With that the per-step cost over grey is ≈ the group rows (cheap) + the wider LU (O(G³)) +
the G formal solutions of stage 2.

## 8. Boundary conditions

Outer edge (143 pc, the CMB at z = 19, T_CMB = 54.5 K): incoming isotropic intensity I_cmb b_g(T_cmb) per group (TITAN's
F = c g E_s − π(2g + 1)I_in with I_in → I_in,g and the VEF_SPH half-zone carry of E_g − E_b,g, `mg_esurf`). At 54.5 K the CMB is entirely in the
IR group (b_opt = 1e-64, the LW and ionising groups at the floor). Inner boundary: the regular centre (F_g = 0). Restart from a grey
checkpoint: E_g = E_cmb b_g(T_cmb) + (E_r − E_cmb)⁺ b_g(T_s), T_s = T_r inside the photosphere (τ ≥ 2/3) and the star's T_eff outside,
F_g = F_r b_g(T_s), renormalised to the exact sums.

## 9. Grey quantities the hydro needs as sums over groups

| quantity | grey run | multigroup |
|---|---|---|
| E_r, F_r (total-energy and gas+radiation rows, ∇·F, work, MLT, Δt control) | unknowns | unknowns, = Σ E_g, Σ F_g (algebraic rows) |
| radiation force on the gas | χ̄ F_r/c | Σ_g χ̄_g F_g/c |
| radiation pressure P_r | f E_r | Σ_g f_g E_g (= f E_r with the grey f, stage 1) |
| gas–radiation exchange | c κ_E E_r − 4π κ_P B in the E_r row | Σ_g in the group rows; the gas sees −Σ through the total row |
| chemistry energy moves (RAYQ, GASQ) | from E_r | from their group (stage 1: the optical group MGQG) |
| Eddington factors f, g, h | formal solution of the grey field | grey (stage 1); per group (stage 2), P_r = Σ f_g E_g |
| L_em, T_eff, photosphere | from E_r, F_r | the same sums; plus L_g per group (the emergent spectrum) |

## 10. Stage 1: the first coupled prototype (`MGROUPS=4`)

**What is built** (commits c2c41e8 … a0cf74b on `multigroup`; every piece default off):

* 4 groups (IR < 0.7551 eV | optical 0.7551–11.2 | LW 11.2–13.598 | ionising ≥ 13.598), E_g and F_g as Newton unknowns 10–17 (neqn 17,
  meqn 17 set by TitanFortran when MGROUPS is set), the grey E_r, F_r rows as sums, the group energy and momentum rows of §3 with the grey
  f, g, h for every group, SLABFLUX/SLABTIE per group, the outer boundary per group with the CMB at T_cmb = 54.5 K, the radiation force
  Σ_g χ̄_g F_g/c, and the group exchange both in the group rows and in the GASROW gas row of the dual-energy zones (the stack's
  internal-energy row has the exchange explicitly: it becomes Σ_g there).
* group opacities = grey means at the iterate × lagged G2 group ratios (§4), table zones 1; H(1s) absorption added to the ionising group
  of case-B network zones; in-group shape floor 3000 K.
* photo-rates from group dilutions W_g (§6): H⁻ detachment (piecewise J_ν), three-level β and P12 (optical group), LW (LW group),
  H(1s) photoionisation (new, ionising group); SLIMP's iterate factor from the optical group; RAYQ/GASQ from the optical group.
* E_g updated in the logarithm (clamp ±3, floor 1e-150 E_r), convergence weight E_g/(E_g + 1e-3 E_r), TR-BDF2 y*/save/restore/error norm
  for the groups; the pivoted banded LAPACK solver (JuliaSolver) for the 17-wide blocks; group fields in the checkpoint (.aux, MG only);
  restart from a grey checkpoint partitions the grey field (§8). Records `mg_rec.jls` every SNAPEVERY steps (MG only).
* Diagnostics: `EDIAG=1` (total-energy rows of every converged stage → ediag.csv), and `TIESMOOTH=ε` (§5 fix 1; test 9).

**Tests** (100 zones; restarts of copies of z100mm's checkpoint `z19z100mmS5_auto0`, step 14400, t = 237841.0 yr, M★ = 0.849 M☉,
photosphere at zone 26, r_ph = 116 R☉, T_eff ≈ 4390 K; z100mm's own phaseB.env + run.env (SLABFLUX=1 SLABTIE=1 SUBZONE=0); 8 threads,
cores 24–31 / 32–39):

1. **Unit test** (julia/mg_test.jl): all passed — Σ_g b_g = 1 (≤ 1.1e-16), b_g against quadrature 2.2e-10, d b_g/dT AD vs 4th-order
   difference 1.6e-10; grey identities of the group means 0.999–1.000 (κ_P, κ_E) and 1.000–1.007 (χ); grey totals of `group_fast` =
   `means_fast` exactly.
2. **Bit-identical off.** 30 steps, `t_mg10` (popiii-s6 43269f6) vs `t_mg11` (c2c41e8): 591/591 state arrays and the .aux identical
   (audit/ckcmp.jl). Repeated for the final tree a0cf74b (with the GASROW, TIESMOOTH, EDIAG and photo-rate hooks): `t_mg20` vs `t_mg10`,
   591/591 + .aux identical.
3. **Grey vs multigroup over 400 steps** (`t_mg13` grey, `t_mg14` MGROUPS=4; steady window 14500–14800, the first 100 steps are the restart
   transient of the group partition):

| | grey `t_mg13` | MGROUPS=4 `t_mg14` |
|---|---|---|
| s/step (steady, 100-step rows) | 1.65–1.71 | 2.45–2.55 (×1.48) |
| Newton failures (TR + BDF) /100 steps | 18–19 | 21–22 |
| model time per 100 steps [yr] | 0.024–0.027 | 0.020–0.022 |
| Δt median [s] | 1.15e4 | 9.5e3 |
| L_em at 1.5 r_ph, median (rms ln) [L☉] | 6507 (0.7 %) | 6323 (1.3 %; one −9 % dip at 14510, after a 1.7e4 s step) |
| L_ph, T_eff (median) | 4467 L☉, 4387 K | 4260 L☉, 4336 K |
| photosphere | zone 26, 115.5–116.3 R☉ | the same |
| max \|Σ_g E_g/E_r − 1\| | — | ≤ 2e-12 (5e-8 once) |
| L through the outer boundary | −1.0e5 L☉ (CMB boundary artefact of the grey run, ~3e-7 of the CMB flux) | −1.0e5 L☉, all in the IR group |

   Sum of the groups against the grey run: the same photosphere zone and radius, L_em 2.8 % and T_eff 1.2 % lower. The emergent spectrum
   at 1.5 r_ph is IR 16.7 % | optical 83.3 % (a 4390 K blackbody has 18.1 % below 0.7551 eV); at r_ph IR 21 %.
4. **Energy conservation** (`EDIAG=1`, 200 steps each, `t_mg15` grey, `t_mg16` MG; the trajectories are those of t_mg13/t_mg14 bit for bit:
   the diagnostic is passive). The total-energy rows evaluated at every converged stage leave an imbalance of **3741 L☉ (grey; p10–p90
   3737–3745) and 3754 L☉ (MG; 3739–3768)** as a luminosity. It is a property of the grey run (constant, ~3.7 % of the outer-boundary
   flux artefact; the dual-energy / flux-tie zones), not of the groups: multigroup adds 13 L☉ (0.35 %) to it, with the gas–radiation
   exchange carried entirely by the group rows.
5. **Photo-rates against STARLIGHT/LOCALJ** (t_mg16, group W_g over the grey W of the same zone, median of 40 records): H⁻ detachment and
   the three-level rates (optical group) **0.97–1.03 × grey** from the photosphere to 1e5 R☉, 1.05–1.40 at 1e5–4e5 R☉; the IR group 0.87
   (not used by any rate); **LW 1.2 × at 116 R☉ rising to 1e3 × beyond 1e4 R☉** (the LW band is only Rayleigh-scattered in the molecular
   envelope, κ_E,LW/κ_E ~ 1e-7, while the grey W falls as the optical light is absorbed by H⁻: the effect multigroup is for, but see 6);
   H(1s) photoionisation (new) W_ion/W ≤ 4e-12 outside the photosphere: negligible at T_eff 4400 K.
6. **Emergent spectrum against the IonFS formal solution** (branch ionfs-v2, julia/ionfs_post.jl, on the same snapshots 14600, 14700,
   14795 of both runs): IonFS gives Qem_LW = 1.1–1.35e41 s⁻¹ and Qem_H = 0.9–1.3e38 s⁻¹ leaving r_ph, i.e. L_LW/L ≈ 1.5e-7 and
   L_ion/L ≈ 1.5e-10 (×12.4 and ×14.6 eV); its T_eff 4332–4339 K (MG snapshots) and 4385–4395 K (grey) match the runs'. **MGROUPS=4 at
   r_ph: L_LW/L = 1.1–1.5e-10 and L_ion/L ≈ 1e-19** — 1e3 and 1e9 below the formal solution. The major groups are right (the IR/optical
   split and the luminosity), the minor groups are not: the LW and ionising photospheres lie elsewhere than the grey one (LW: Rayleigh-thick
   where the grey field is thin; ionising: H(1s)-thick within ~1e3 cm above the photosphere), and these groups are closed with the grey
   f and the grey flux-tie h, which fix their free-streaming geometry at the grey photosphere. Their luminosities, and with them the
   LW and H(1s) rates, are not usable before stage 2a (per-group Eddington factors). Also tried (t_mg17, t_mg18): Newton in E_g instead of
   ln E_g with every group converged above 1e-12 E_r — the ionising group then flip-flops between ~1e-22 and ~1e-49 of E_r at the
   photosphere from step to step and Δt fell to 446–2228 s; reverted.
7. **Cost** (CPLPROF, 50-step blocks after the restart): grey matgen 64.5 s per 50 steps (1514 Newton iterations, 42.6 ms each), matslv
   (block pentadiagonal, static 9×9) 0.6 s; MG matgen 84.7 s (1588 iterations, 53.3 ms each: ×1.25 per iteration from the 102- instead of
   54-partial Duals, ×1.05 more iterations) and ~16–18 s of banded LAPACK LU (n = 1768, bandwidth 50; not timed by CPLPROF since the
   solver is replaced after the timers are installed: from the difference of the totals). The ×1.48 per step splits into ×1.31 Jacobian
   assembly and ~0.17 for the wider linear solve; both have the remedies of §7/§11. The run's static block-pentadiagonal solver
   (BlockSolver.Penta, pivoting only inside a block) is not usable for neqn = 17: `t_mg23` (MGSOLVER=penta) diverged in its first step
   (corrections 1e10–1e13, all coupled zones decoupled, Δt → 1e-5 s); the pivoted banded LU is required as it stands.
8. **Final tree** a0cf74b (log update with the 1e-150 floor): `t_mg21`, 200 steps — 22 failures/100, 2.39 s/step, L_em(1.5 r_ph) median 6289
   L☉ (rms 0.8 %), T_eff 4337 K: the same as t_mg14 within its scatter (not bit-identical to it: the floor acts in the cold outer zones).
9. **TIESMOOTH=1e-3 on the grey stack** (`t_mg22`, 400 steps): failures 18, 19, 18 per 100 steps against 18, 19, 19 (t_mg13), s/step and
   Δt the same: SLABTIE's max(h, 0) does not bind at this checkpoint (h stays far from 0), so fix 1 is not what limits the Newton here.

## 11. Plan of the remaining stages

1. **Stage 2a — per-group Eddington factors (needed for the minor groups, test 6).** TitanCore.eddfac! on group arrays (χ_g, κ_P,g,
   source with E_g and EDDSCAT, boundary I_cmb b_g), f_g, g_g, h_g at every iterate (VEFIMP), radiation pressure Σ f_g E_g in the total
   rows and MLT, the tie with h_g. Test: L_g at r_ph against ionfs_post (LW 1.5e-7, H 1.5e-10 of L at this checkpoint); grey vs MG L_em.
2. **Stage 2b — exact group opacities at the iterate** (group_fast generic in the number type inside NLOIT's material sweep; the lagged
   ratio goes) and **lines per group** (the LTE line factor of NLOIT assigned to the groups that hold the lines).
3. **Stage 2c — photo-rates fully from the groups**: one iterate factor per process (q_g of its group) instead of SLIMP's one per zone;
   LOCALJ with the piecewise J_ν (direct integrals or per-group tables) so that 6–8 groups are exact; **photon-conserving LW** (§6:
   σ_eff = −d ln f_sh/dN in the LW group's κ_E, the dissociation rate from E_LW at the iterate, the LW heat from the LW group); He I/He II
   photoionisation from the He groups (6–8-group sets). This retires STARLIGHT's W, SLWFLUX, SLIMP and SLSHIELD for these processes.
4. **Stage 3 — cost and Newton:** (a) a block solver with pivoting across neighbouring blocks (e.g. BlockSolver's zone-pair tridiagonal form
   with partial pivoting inside the 34×34 pair blocks) instead of banded LAPACK (−0.17 of the ×1.48; the unpivoted Penta fails, test 7); (b) two-pass colouring or an analytic Jacobian of the group rows so that the material sweep stays
   at 54 partials (target ≤ ×1.2 over grey for 4 groups); (c) the two Newton fixes of §5: TIESMOOTH (built; no effect at this checkpoint, test 9 — to be
   repeated where the photosphere is in the dark/bright switching of POPIII_SLABFLUX.md) and JFNK with the banded Jacobian as preconditioner for exact f, g, h derivatives (grey first).
5. **Stage 4 — validation and production tests:** 6- and 8-group runs; L_em steadiness over ≥ 1e4 steps; a from-IC 100-zone start through
   protostar formation (IC=abn02red on the z100mm mesh) with MGROUPS from the first step; the IonFS comparison at the break-out phase.

## 12. The grey stack's total-energy imbalance (~3.7e3 L☉ at z100mm step 14400)

Question (PI): the conservative total-energy rows, evaluated at every converged stage (`EDIAG`), leave 3741 L☉ in the grey run — 60 % of
the emergent L. What is it? Diagnostic `EDIAG=2` (TitanResidFast `zone_diag`: per zone the total-energy row residual, which rows the solve
used, the FLUXTIE weights, the radiation momentum row in its moment-equation form, the advection-limiter mismatch and ΔK). Runs `t_mg24`
(40 steps), `t_mg26` (30 steps) from the same checkpoint, grey stack, medians over the BDF stages (scripts in the scratchpad mg/ezone*.jl).

| where | imbalance [L☉] | of which explained by |
|---|---|---|
| zones 1–26 (the star inside the photosphere, total-energy rows) | −0.4 (Σ\|zone\| 11) | — (solved rows: 0) |
| zones 27–50 (r_ph = 116 R☉ … 3.6e4 R☉, the 10–49 km/s infall; internal-energy rows) | **+442** (Σ\|zone\| 1704, alternating ±150…580 zone to zone around the infall maximum) | advection-limiter mismatch −0.2; u·R_mom 0.1; the rest is the kinetic-energy/work discretisation of the dual-energy rows (ΔK/dt of these zones 298) |
| zones 51–87 (3.6e4 R☉ … 1.1 pc) | +19 | u·R_mom 23 |
| zones 88–100 (1.1 … 143 pc, CMB-dominated cloud; internal-energy rows, tie weight 1) | **+3280** | **u·R_mom = 3391** (zone by zone, e.g. 394 vs 375, 705 vs 555, 682 vs 875 L☉) |
| total | 3740 | — |

* It is **entirely the non-conservation of the dual-energy zones** (74 of 100 zones, every zone outside the photosphere, use the
  internal-energy/GASROW row: ETASW 0.6, ETASWTH, ETASWTIE); the zones on the total-energy row balance to −0.4 L☉. It is not accretion
  luminosity at the inner boundary or the moving mesh (zones 1–26: −0.4 L☉), not the chemistry or H₂ line sources (they are in both rows),
  not the table/network blend (interior zones balanced), not the separate van Leer limiters of e_g and E/ρ (−0.2 L☉).
* **88 % (3280 L☉) is in the outer cloud at 1–143 pc** and is, to 3 %, u × (residual of the radiation momentum equation) at the FLUXTIE
  interfaces: there the radiation row is the tie F = h c E, the gas feels χ̄F/c, but the energy rows do the compression work f E ∇·u of a
  field that is the CMB (E ≈ a T_CMB⁴), so total − kinetic − radiation rows leave −u·(∇·P + χρF/c) (the effect ETASW_TIE's docstring
  describes). It is a bookkeeping inconsistency of the CMB field at the 1e-2 level of its compression work, not an energy flow from or to
  the star or the gas (the gas row is solved; the radiation force of the tie flux there is ~1e-7 of gravity).
* **12 % (442 L☉ net, 7 % of L_em) is near the star**, in the supersonic infall just outside the photosphere, alternating in sign zone to
  zone: the internal-energy row plus the momentum row do not conserve e + K discretely where K/m ≫ e (kinetic-energy changes of the
  momentum discretisation are not returned as heat, or vice versa). This one can matter for the emergent luminosity.
* **Is it a bug?** No coding error found: total-energy rows everywhere (ETASW = ETASWTH = ETASWTIE = 0, `t_mg25`) conserve to −1.1 L☉ but
  the run fails after 23 steps (the photosphere jumps from zone 26 to 31; 66 Δt rejections against 8), as the dual-energy switches were
  introduced to prevent. The imbalance is the price of the dual-energy formulation and of the tie closure. Remedies (not built): (a) for the
  ledger, book u·R_mom of the tie zones explicitly (it accounts for 3391 of the 3740 L☉); (b) for the near-star 442 L☉, a dual-energy
  synchronisation that puts the kinetic-energy discretisation error into the gas row (e.g. the total-energy row's ΔE − ΔK − ΔE_r as a source
  of the internal-energy row where η is near the switch), to be tested for stability against t_mg25's failure.

## 13. Stage 2a: Eddington factors per group (then `MGVEF=1`; the only behaviour since §19)

**Built** (4cec045 … 691cbfe): after the grey formal solution, the same formal solution (TitanCore.eddfac!) once per group on the group's
arrays swapped into the state it reads (χ_g ≥ κ_E,g; absorption κ_E,g; emission κ_P,g B_g; the rest of χ_g scattering; E_g, F_g; the
background I_cmb b_g(T_cmb)) at every Newton update (VEFIMP): f_g, h_g, the tie weights of the group's own depth, g_g; the grey `fedd` then
holds Σ f_g E_g/Σ E_g (radiation pressure of the total rows). Diagnostic H⁺ per interface (TitanCore.HPON/HPLUS) for the IonFS comparison.

**Bugs found and fixed on the way** (each was in the stage-1 numbers of §10 as well, mostly invisible with the grey closure):
* in-group shape floor 3000 K broke Σ b_g κ_P,g = κ_P below 3000 K (§4; cold cloud over-emitted into the optical group);
* a group whose weight underflows (cold gas) fell back to the grey κ — the LW band then had the cold gas's H₂ CIA absorption and died within
  a few zones; now the in-group shape at max(T, 3000 K) for that group only (550e139, 7dcc6ae);
* restart partition: the grey excess over the CMB was given T_eff's shape out to 143 pc; now only in one star-lit region from the photosphere
  (2498992, de0b1ae). A residual: the outer cloud's own excess at 450–1100 K puts ~1e-3 of E into the optical group there, which streams out
  with f ≈ 1 (up to ~3e7 L☉ through the outer boundary; emission and absorption there are ~1e-5 L☉: a restart transient of the CMB-dominated
  cloud, decoupled from the gas);
* ionising group: H(1s) Kirchhoff–Milne recombination emission (K1 × n_1s*/n_1s, above 2500 K only — below, the floor of b_g broke the
  e^{χ/kT} e^{−hν/kT} cancellation, t_mg32) and the matching ground-state recombination a1 x_e x_H⁺ in the network where the group
  photoionisation is on (b1e0c5b, de0b1ae; `MGMILNE`).

**Runs** (z100mm checkpoint, 4 groups): `t_mg36` (300 steps, EDIAG), `t_mg38` (120 steps, H⁺ diagnostic), against grey t_mg13/t_mg15.

| | grey | stage 2a (t_mg36), steps 14500–14700 |
|---|---|---|
| s/step | 1.65–1.71 | 2.7–2.8 (×1.65: 5 formal solutions per update) |
| Newton failures / 100 steps | 18–19 | 10 |
| model time / 100 steps [yr] | 0.024–0.027 | 0.019–0.021 |
| L_em at 1.5 r_ph, median (rms) | 6507 L☉ (0.7 %) | 6493 L☉ (7.5 %; 6221–7361) |
| T_eff | 4387 K | 4498 K |
| total-energy imbalance (EDIAG) | 3741 L☉ | 3801 L☉ |

**Gate** (L_LW/L, L_ion/L at r_ph against ionfs_post.jl on the same snapshots 14450–14515 of t_mg38; IonFS's Qem × the band's mean photon
energy in B_ν(T_eff), over its L_run; MG: the outgoing L⁺_g = 4π r² 4π H⁺_g of the group's formal solution at the same interface — the light
leaving the layers inside, as IonFS's vacuum-above formal solution — and the net L_g):

| step | IonFS L_LW/L | MG L⁺_LW/L (net) | IonFS L_ion/L | MG L⁺_ion/L (net) |
|---|---|---|---|---|
| 14450 | 9.7e-9 | 2.3e-5 (3.8e-8) | 2.1e-11 | 3.3e-20 (3.5e-20) |
| 14475 | 1.5e-8 | 1.9e-5 (4.3e-8) | 8.2e-12 | 2.4e-12 (1.5e-20) |
| 14500 | 3.0e-9 | 1.9e-5 (3.0e-8) | 2.5e-11 | 5.5e-12 (2.5e-20) |
| 14515 | 3.6e-9 | 1.8e-5 (3.8e-8) | 2.1e-11 | 3.7e-12 (2.1e-20) |

* **Ionising group: outgoing flux 0.18–0.29 × IonFS in 3 of 4 snapshots** (×3.5–5.5 low; 1.6e-9 × in the fourth): from 1e-9 × before the
  Milne emission. The net flux at r_ph is ~1e-20 of L because the first zone above (neutral H, τ_1s ~ 1e9) absorbs and re-emits at nearly the
  same T — IonFS's own ray finds Qesc_H = 0 as well. Gate (factor 2): not met, close.
* **LW group: not met.** The outgoing LW at r_ph is 1e3–6e3 × IonFS's, the net 3–10 ×. The LW band above the photosphere is
  Rayleigh-thick (χ_LW ≈ 0.36 cm² g⁻¹, f_LW = 1/3 out to ~1e3 R☉), so outgoing ≈ inward and the net is small; the outgoing excess means the
  LW source function inside r_ph is far above IonFS's. Suspect: the H(n ≥ 2) bound–free absorption at 11–13.6 eV, which the G2 continuum
  takes from the three-level atom's n = 2 while IonFS puts n ≥ 2 in LTE relative to the continuum (n_e n_p) — a lower n = 2 makes the
  LW band thermalise deeper and hotter. To check next: the group's κ_E,LW/κ_P,LW profile against IonFS's coeffs3 on the same snapshot.
  Also note IonFS's LW itself varies ×5 between snapshots 25 steps apart (9.7e-9 → 3.0e-9), and was 1.5e-7 on the stage-1 snapshots
  (different photosphere structure).
* Bit-identity off for the final tree 691cbfe (EDIAG=2, zone_diag, HPON, starlight_terms' a1, all stage-2a code): `t_mg39` vs `t_mg10`,
  30 steps, 591/591 state arrays + .aux identical.

## 14. Stage 2a continued: the LW and ionising groups, the outer-cloud optical reservoir

**H n = 2 populations by regime** (snapshot t_mg44 step 14500; three-level n₂ = the network's own atom, PopIIIMaterial.n2_tl3, recorded per zone;
scratchpad mg/lwpop.jl):

| zones | T [K] | n₂ three-level | n₂ Boltzmann (LTE to 1s) | n₂ Saha (LTE to the continuum, IonFS) |
|---|---|---|---|---|
| 18–25 (τ ≫ 1, below the photosphere) | 37 400–7 400 | 1.7e14 … 1.6e9 | 1.6e14 … 1.4e9 | 1.1e14 … 1.6e9 (all three within 1–30 %) |
| 26–27 (photosphere) | 4 830–4 800 | 8.4e4, 3.0e4 | 9.2e4, 4.1e4 | 5.8e4, 1.1e4 |
| 28–34 (envelope 120–260 R☉) | 3 320–1 825 | 6.1e2 … 7.4e-3 | 0.43 … 1.2e-15 | 1.1e4 … 52 |

The three-level atom is LTE where it should be (the thick layers), between the two LTE limits in the photosphere, and in the irradiated,
over-ionised envelope 10–1e4 × below IonFS's LTE-to-continuum value (the star's Balmer-continuum photoionisation β ≫ collisional and
recombination rates empties n = 2; b₂ = 1 relative to the continuum ignores that) and 1e3–1e13 × above Boltzmann (recombination and Lyα
pumping fill it). The three-level value is the physically right one in every regime; IonFS's assumption is not, above the photosphere.

**But n₂ does not set the LW band.** The H n ≥ 2 bound–free at 11.2–13.6 eV is ≤ 1e-5 cm² g⁻¹ above the photosphere with *any* of the three
n₂ (G2 group means, the same snapshot), against the band's Rayleigh extinction 0.36 cm² g⁻¹ — so the n ≥ 2 absorption, now added to the LW and
ionising groups of case-B zones with the three-level n₂ (0998ca3, `MGN2BF`), changes nothing measurable (t_mg45 with it off has the same perf record as
t_mg36 to the printed digits). What did set the LW band was a **numerical** absorber: the ground state's dissolved Lyman pseudo-continuum (G2Opacity's DIS window,
q ≤ 1/1600, 13.590–13.598 eV) fell into the LW group, and in a Planck (arithmetic) mean that 6e-4-wide sliver of ground-state absorption raised
κ_E,LW ×1e5 in the LTE zones below the photosphere (zone 25: 9.5e3 against 0.07 cm² g⁻¹), thermalising the whole band at 7400 K one zone
deeper and 2500 K hotter than the front. The RADGROUPS trial met the same effect with the H(1s) continuum. Fix: the LW/ionising edge at the
onset of the dissolved window, χ_H(1 − 1/1600) (34318f8); the sliver is ground-state absorption and belongs with the ionising group.

**Gate after the fix** (t_mg47, IonFS on the same snapshots 14420–14450; MG outgoing L⁺ of the group formal solution at r_ph):

| step | IonFS L_LW/L | MG L⁺_LW/L | MG net L_LW/L | IonFS L_ion/L | MG L⁺_ion/L |
|---|---|---|---|---|---|
| 14420 | 1.9e-8 | 1.6e-6 (×82) | 6.2e-8 | 6.2e-12 | 2.3e-11 (×3.7) |
| 14430 | 2.8e-8 | 5.7e-7 (×20) | 1.2e-8 | 4.6e-11 | 2.5e-11 (×0.54) |
| 14440 | 2.0e-8 | 2.1e-6 (×104) | 5.5e-8 | 6.4e-12 | 3.9e-12 (×0.61) |
| 14450 | 2.0e-8 | 1.8e-6 (×91) | 5.6e-8 | 7.9e-12 | 6.6e-12 (×0.83) |

* **Ionising group: within a factor 2 of IonFS in 3 of 4 snapshots** (0.54–0.83; ×3.7 in the fourth), from 1e-9 × before the Milne emission and
  0.2–0.3 × before the edge fix. The quantity is exponentially sensitive (hν/kT ≈ 30 at the photosphere: a factor 2 is 2–3 % in the
  temperature where τ_1s = 2/3), so this is as close as one zone across the photosphere allows. The net ionising flux at r_ph stays ~1e-20 L:
  the first zone above absorbs it (IonFS's ray: Qesc_H = 0 as well).
* **LW group: outgoing 20–100 × IonFS** (from 1e3–6e3 ×), net 0.4–3 × IonFS. Remaining difference identified, not settled: the band is
  Rayleigh-scattering-dominated above the photosphere, and the two codes' Rayleigh cross sections differ there — G2 holds σ_Ray at its
  1300 Å value at shorter wavelengths, IonFS (ContinuumOpacity.alpha) uses the λ⁻⁴…λ⁻⁸ fit, 4× larger at 1050 Å and diverging toward Lyα;
  neither is valid between the Lyman lines (resonant scattering). More scattering above the front lowers the outgoing flux. IonFS used no
  sub-zone points on these snapshots (its front-zone opacity scaling inactive, checked), and its own L_LW varies ×1.5–5 between snapshots.
  The LW group is therefore fixed to within its scattering opacity; the dissociation rates need 2c's photon-conserving LW absorption anyway.

**Outer-cloud optical reservoir: a restart transient, not a decoupling.** In t_mg36 the optical group's energy in zones 94–100 (1.1–143 pc),
U = 7.5e49 erg, falls at dU/dt = −0.97 × its outflow through the outer boundary (4.0–4.8e5 L☉), while its emission − absorption there is
−1.4e-4 L☉ (scratchpad mg/optout.jl): it is the initial partition draining at c over r/c ≈ 1e3 yr, decoupled from the gas. The grey field's
excess over the CMB out there cannot be the cloud's own emission (1e-4 L☉): it is the grey CMB solution's residual. Re-partitioning it at
T_CMB, or as starlight capped at L/(4πr²c) plus CMB, removes the reservoir but fails the first step (t_mg40–t_mg43: the logarithmic update
cannot move a group by 200 decades in one step); kept the working partition (14a45a9). From the IC the reservoir does not arise.

**Stage-2a state** (t_mg47, 34318f8, steps 14420–14480): L_em(1.5 r_ph) median 6987 L☉ (rms 4.4 %), T_eff 4502 K, Δt median 9.8e3 s.

## 15. Round 3 (branch `mg2b`, from 4f1fe93): grey bit-identity re-verified

`t_mb01` (tree 4f1fe93, MGROUPS unset, 30 steps from the step-14400 checkpoint, cores 32–39) against `t_mg10` (popiii-s6 43269f6, the same
30 steps): the final checkpoint's `.aux` is byte-identical (md5 d540837f…), the decompressed `.tstate` is byte-identical (md5 246c3b43…; the
`.gz` files differ only in the gzip header time), and the run logs (121 lines, every step record) and segment outputs (161 lines) differ in
two wall-time lines only (`EP call 0.1 s`, the NLOIT timing summary); the NLOIT counters are equal (1 014 006 ratio evaluations, 88 779 slope
refreshes). The step-14400 checkpoint used here and below is the copy kept in `runs/t_mg10/ck` (md5 aux 705f40fe…, tstate.gz f9ff2599…):
z100mm's own `z19z100mmS5_auto0` rotates and was overwritten by the production run on 16 Nov 21:08.

## 16. H(1s) Rayleigh and Raman scattering between Lyα and the Lyman limit (46fe7b4)

**Source.** Kokubo 2024, MNRAS 529, 2131 (bibcode 2024MNRAS.529.2131K, refereed; found with `scix`, read in full): the Kramers–Heisenberg
cross sections of H(1s) with analytic radial integrals over the bound and continuum intermediate p states — Rayleigh σ_1s(ν) and the Raman
branches 1s → 2s, 3s, 3d, 4s, 4d. `julia/HScatter.jl` evaluates his sums (1000 bound states + the continuum integral, the tail of the bound
series analytic), every resonance regularised with the radiative width of np → 1s (his footnote). Branches to n_B ≥ 5 (open above 13.05 eV)
are not in the source and are left out. Checks in a limit (`hscatter_check`, scratchpad `mb/t1/hs1.jl`): the Thomas–Reiche–Kuhn sum of the
oscillator strengths used is 1.0000000018; σ_1s/(σ_T ε⁴) → 1.26563 = (9/8)² for ε = hν/χ_H → 0 (static polarisability 9/2 a₀³).

**What it shows for the LW band (11.2 eV … χ_H(1 − 1/1600)).**

* The low-energy expansion in the code (ContinuumOpacity `σ_ray_H`, Z0Physics `σRay`; G2 holds it at its 1300 Å value below 1300 Å) is not
  valid there: ratio to the exact σ_1s 0.855 at 1900 Å, 0.538 at 1520 Å, 0.124 at 1302 Å; between the Lyman lines σ_1s swings between
  resonances and zeros (zeros at 11.682, 12.604, 12.986, 13.182, 13.297, 13.370 eV).
* **There are no escape windows.** At the zeros of σ_1s the Raman cross section is 3–3.9e-25 cm²: where the atom does not scatter the photon
  it destroys it (1s → 2s + an optical photon at ν − ν_Lyα; the atom returns by two-photon decay). Band minima: σ_1s + σ_Ram ≥ 2.4e-25 cm²,
  effective (thermalisation) cross section k_eff = √(3σ_Ram(σ_Ram + σ_1s)) ≥ 2.3e-25 cm²; medians σ_Ram 2.8e-25, σ_1s 1.2e-24, k_eff 7.6e-25 cm².
  Band-averaged transmission ⟨e^{−k_eff N}⟩: 0.36 at N_HI = 1e24, 6.0e-3 at 1e25, 3e-13 at 1e26 cm⁻².
* On the checkpoint the H I column above the photosphere is 7e26–1e27 cm⁻², and every zone out to 3e4 R☉ has N_HI ≥ 1e25 cm⁻² on its own:
  **no LW photon of the photosphere reaches the envelope**; each zone is LW-thick in the best window.

**In G2 and the group rows** (the LW group only; the grey means are untouched — `group_fast` takes the sub-grid only when called with group
edges):

* flux mean χ_LW: harmonic mean of (smooth extinction of the G2 nodes, interpolated) + n_1s(σ_1s + σ_Ram) over a sub-grid of the band, weight
  ∂B_ν/∂T at the colour temperature of the light. The sub-grid is 221 cells of 32 table cells each with the cell's harmonic H(1s) cross
  section (its transmission) and Raman share; against the full 7044 cells the group means agree to 0.7 % where H(1s) dominates the extinction
  (smooth part ≤ 1e-25 cm² per atom) and χ_LW to 2 % everywhere (`mb/t1/compress.jl`; per atom χ ≈ 7.0–9.3e-25 cm², κ_Ram ≈ 2.1–2.9e-25 cm²
  for 3000–30 000 K);
* energy mean: Raman is absorption in the LW group's energy row, κ_Ram = χ_LW ⟨σ_Ram n_1s/χ_ν⟩ (in-group spectrum ∝ B_ν/χ_ν of the harmonic
  mean: photons sit where the gas is transparent, and the unresolved line cores, where they are trapped, do not count);
* the destroyed LW photon's energy goes to the group of its Stokes partner (optical), and the inverse process (H(2s) + optical photon → LW
  photon, anti-Stokes) is the LW group's source: per mass c κ_Ram (C E_opt − ω E_LW) into the LW group and the same out of the optical one
  (`M.xfer`; energy moves between two groups, the gas row sees nothing). By detailed balance of the Kramers–Heisenberg rate (g_1s = g_2s)
  the inverse-Raman source function is S_ν = (2hν³/c²)(n_2s/n_1s) n̄(ν − ν_Lyα), so C = (n_2s/n_1s) e^{E21/kT_w} b(pair)/b_opt(T_w) with the
  three-level atom's n_2s = n_2/4. LTE check on the run: zones 18–25 (T 7500–37 000 K) have C E_opt/(ω E_LW) = 0.99–1.05 with
  E_LW = (1.000 … 1.02) aT⁴ b_LW(T), and Σ_g of the transfer is zero by construction.

**What LW field the photosphere really emits** (t_mb02, 60 steps from the checkpoint; scratchpad `mb/t1/ram.jl`, `n2var.jl`):

* Out of r_ph: L⁺_LW/L = 1.7–2.4e-9, i.e. 0.86 × the Planck share b_LW(T_eff) = 9.8e-10 of the outgoing half flux (Q ≈ 2.0–2.6e39 s⁻¹). The
  net LW flux through r_ph is 3e-11 L: the zone above takes it (Raman conversion time 1/(c κ_Ram ρ) = 0.1–0.4 s in zones 26–28).
* Above that the LW field is not photospheric light at all: it is the local inverse-Raman source of the optical field on excited hydrogen,
  E_LW/E_r ≈ 1e-13 in zones 29–37 (126–740 R☉), with the transfer in balance to 1e-3 (C E_opt/ω E_LW = 1.000). In units of the
  undiluted photospheric Planck function, S_LW/B_LW(T_c) = 1.3e-4 (zone 29), 5.7e-6 (zone 36), 1e-5 (zones 40–41, 2000–3000 R☉),
  1e-6 (zone 44), 3.6e-11 (zone 49, 2.8e4 R☉) — against the geometric dilution W = 0.47, 1.7e-2, 1e-3–5e-4, 9e-5, 6e-6 that STARLIGHT's
  W B_LW(T_c) assumes: **lower by 3.6e3, 3e3, 50–100, 90 and 1.7e5.**
* **Uncertainty.** The field is proportional to n_2s. (i) The three-level atom pumps 1s → 2s with Λ_2s n̄(Lyα), one power of the dilution;
  two-photon absorption is quadratic in the occupation numbers (W² in the Wien limit). With W² the source drops by 1.7–3× in zones 29–36,
  not at all in zones 38–45 (there recombination, α_B n_e n_p, feeds n = 2: 77–100 % of the rate) and by 50–160× in zones 46–49. With no
  radiative pumping at all (collisions + recombination) the numbers are the same as with W² to 10 % outside zone 33 and 1.4–5× below the
  code's inside. (ii) Where recombination feeds it, n_2s is linear in n_e n_p, so it carries the envelope's ionisation fraction (4e-5) with
  it. (iii) The Sobolev escape factor of Lyα in the atom has no meaning at N_HI ≥ 1e25 (damping wings thick); the Lyα term is 1e-5 of the
  pumping, so this does not matter for n_2 here. (iv) 2s : 2p is taken statistical. **So: photospheric LW is known to ~20 % (two codes) and
  irrelevant beyond one zone; the envelope's LW field is 2–5 decades below the diluted Planck field, and uncertain by a factor ≤ 5 inside
  1e4 R☉ and by two orders of magnitude beyond**, where it is too small to matter either way (S/B ≤ 1e-9). The grey stack's `n2_tl3` is not
  changed in this round (it belongs to the chemistry); its one-power pumping term should become W n̄ when the three-level atom reads the
  group fields.

**IonFS with the same physics, and the gate redone.** `mb/ifs/ionfs_post_ram.jl` (a copy of tree 3499e3a's `ionfs_post.jl`): 440 extra
frequencies on the HScatter cells in the LW group, exact σ_1s in place of the expansion, Raman as absorption with the LTE source (the layers
inside r_ph, where the run's n_2 is within a factor 2–4 of LTE). Its LW photon rate out of r_ph falls from 1.1–1.7e40 s⁻¹ (control with
the old opacity on the same frequencies, `IFSRAM=0`) to 2.0–2.6e39 s⁻¹. Gate on t_mb02's snapshots (MG outgoing L⁺_g of the group formal
solution at r_ph over the net L):

| step | IonFS L_LW/L | MG L⁺_LW/L | ratio | IonFS L_ion/L | MG L⁺_ion/L | ratio |
|---|---|---|---|---|---|---|
| 14430 | 1.99e-9 | 1.73e-9 | 0.87 | 7.5e-12 | 3.4e-12 | 0.45 |
| 14440 | 2.25e-9 | 1.85e-9 | 0.82 | 8.9e-12 | 6.0e-12 | 0.67 |
| 14450 | 2.94e-9 | 2.36e-9 | 0.80 | 1.34e-11 | 8.0e-12 | 0.59 |
| 14460 | 2.15e-9 | 1.96e-9 | 0.91 | 8.3e-12 | 7.2e-12 | 0.87 |

**LW: 0.80–0.91 of IonFS on all four snapshots (was ×20–104).** The whole previous discrepancy was the scattering opacity: with the exact
cross sections and the Raman sink in both codes they agree to 20 %. Ionising: 0.45–0.87 (one snapshot 10 % outside a factor 2; §14's
sensitivity argument applies). IonFS's carrier beyond r_ph has no H I scattering or Raman (its `Qesc_LW` ≈ 0.94 Q_em is therefore not an
escape rate; with the band transmission above nothing escapes).

Newton failures and Δt of t_mb02 are those of t_mg47 over the same 60 steps (37 against 38 rejected attempts in the restart transient).

## 17. Stage 2b: group opacities, emissivities and lines at the Newton iterate (e88345e)

**What changed.** The lagged group ratios of stages 1/2a (`mg_ratios!`, `MGMILNET`; §4, §6 describe them and are superseded by this
section) are deleted, not switched off. Each zone carries a record of 7 numbers per group — the ratios r_E, r_P, r_R of the group's
κ_E, κ_P, χ_R to the grey means, the group's continuum κ_E, the ground-state photoionisation opacity K_1, its Milne emission A_1, the
Raman opacity — plus two per Raman pair and (n_2, T_w): `mg_nv() = 7G + 2·pairs + 2` (`MGRad.jl`). The record is evaluated from the
iterate's own composition in the same NLOIT call that gives the grey ratios (`material_opacg!`), and `mg_level` turns it and the grey
(χ, κ_E, κ_P) of the same level into (κ_E,g, χ_g, κ_P,g, A_1). The old time level's record is the converged record of the previous
step (`mg_fill!` → `MGVN` → `MGVO`).

* **Deviation from the plan ("group_fast generic in the Dual type").** The record's Jacobian is the same finite-difference log-slope set
  as the grey ratios' (`mg_slope`, from the same `nlo_ratio` calls), not an AD pass through `group_fast`: the grey NLOIT ratios are
  linearised that way, and one linearisation for both keeps Σ_g rows = grey row exact in the Jacobian as well as in the residual.
* **Lines per group.** The table opacity's line part is distributed with the Z0 line list's Planck-weighted shares per group
  (`mg_linetab!`, `mg_lines`) instead of the continuum ratio. The Lyman lines inside the LW band are left out of the shares (they are
  the Lorentz wings of §16, already in the continuum on the fast grid).
* **Defect of stage 2a found on the way.** 2a formed κ_E,g = (lagged continuum ratio) × (line-including grey κ_E) and floored
  χ_g ≥ κ_E,g. At step 14450: zone 23 κ_E,LW 2.27e3 cm² g⁻¹ (2a) against 49.5 (2b) with χ_LW ≈ 90; zone 26 0.381 against 0.119;
  ionising group zone 26 8.0e6 against 2.7e6, zones 22–25 30–60× lower in 2b. In the interior the floor was 3–400× the Rosseland
  mean, so the groups did not sum to the grey diffusion flux. Replaced by `mg_chif` (Rosseland mean where the zone is thick,
  absorption mean with weight e^{−χ_R m} where it is thin; docstring in `MGRad.jl`).
* **Milne emission** A_1 from the same K_1 at the iterate. Unit check of the implied α_1 against α_A − α_B: 1.065 (5000 K),
  1.018 (1e4 K), 1.006 (2e4 K), 1.011 (4e4 K).
* **Line cooling** of the coupled zones (H₂ rovibrational, HD, atomic; LINEIMP at the iterate) is emitted into the IR group (`Mg.lcl`,
  group `MGLG`) with the gas row's time weights, so Σ_g rows = grey row holds to round-off (as committed in 2b only where the gas row
  is the internal-energy row; made consistent for every row type in §20.4). *Merge note for energy-compat:* this is one
  added source term in the group energy row; the tie and boundary rows are untouched.
* **`G2Opacity.group_fast`**: function barrier and task-local scratch (3.5 ms and 1.8 MB per call → 0.1 ms and 5 kB); the LW band is
  carried on a compressed sub-grid (`LWSTRIDE = 32`, §16).

**Sum identities.** By construction E_r = Σ E_g (max |ΣE_g/E − 1| ~ 1e-12 in t_mb06) and the gas row's exchange is Σ_g. The record's
own Σ_g r_E b_g − 1 ≈ 0.9–1.1e-2 and Σ_g r_P b_g − 1 = 3.6e-2 are *not* a multigroup error: they are the Planck-quadrature error of the
484-node fast grid against the analytic band fractions b_g (bquad.jl: optical +6.4e-2 at 186.5 K; LW −1.0…−6.6e-2; ionising +0.71 at
1000 K falling to +5e-4 at 7.5e4 K). With the grid's own weights the identity holds to ≤ 1e-15 (unit test). The group means are
normalised to the grey mean of the same grid, so the error does not enter the rows.

**Grey path** (MGROUPS unset) bit-identical to 4f1fe93: t_mb05 against t_mb01, logs and final state md5
246c3b436a9e423b0c6fafc6710cfec4.

**Gate: 200 steps from the step-14400 checkpoint, 100 zones, MGROUPS=4 (t_mb06) against stage 2a (t_mg47).**

| | 2b (t_mb06) | 2a (t_mg47) |
|---|---|---|
| Newton failures / 200 steps (per 50: ) | 59 (31, 9, 10, 9) | 44 (31, 11, 5, 5) |
| error rejects | 0 | 8 |
| Newton failures after step 14450 | 28 | 13 (+8 error rejects) |
| median Δt of the failing attempts | 1.95e4 s | 1.33e4 s |
| simulated time per step, steps 14460–14600 | 8.7e3 s | 7.9e3 s |
| L(1.5 r_ph) median, rms | 6993 L☉, 0.47 % | 6955 L☉, 5.45 % |
| r_ph median | 115.64 R☉ | 116.34 R☉ |
| T_eff from (L, r_ph) | 4908 K | 4887 K |
| wall time per step, steps 14500–14600 | 2.96 s | 2.99 s |

*Verdict, not weakened:* the Δt gate is met in steady state (+11 % simulated time per step, +17 % per wall second); the luminosity noise
is 10× smaller; cost equals 2a (grey: 1.65–1.71 s/step). **The Newton-failure count is higher (59 against 44; all rejects 59 against
52).** The pattern is the step controller (limited by x_H2 in 171 of 200 accepted steps) growing Δt until Newton fails and cutting by
4; 2b's failure threshold lies 1.5× higher, and the controller probes it more often. Counting all rejected attempts per simulated time after step 14450: 2b 28 per
1.22e6 s, 2a 21 per 1.10e6 s — still 20 % more. Suspected obstruction: the controller has no memory of the last
failing Δt; this is the grey stack's controller and was not changed.

**IonFS gate on stage 2b (t_mb06, snapshots 14450/14500/14550/14600, `mb/ifs06`): NOT met.**

| step | IonFS L_LW/L | MG L⁺_LW/L | ratio | IonFS L_ion/L | MG L⁺_ion/L | ratio |
|---|---|---|---|---|---|---|
| 14450 | 3.68e-9 | 1.32e-9 | 0.36 | 1.74e-11 | 2.76e-12 | 0.16 |
| 14500 | 3.64e-9 | 1.33e-9 | 0.36 | 1.72e-11 | 3.01e-12 | 0.17 |
| 14550 | 3.70e-9 | 1.35e-9 | 0.36 | 1.77e-11 | 3.02e-12 | 0.17 |
| 14600 | 3.47e-9 | 1.45e-9 | 0.42 | 1.63e-11 | 3.47e-12 | 0.21 |

Suspected reasons, partly quantified, none proven:
1. *2a's agreement (0.80–0.91; §16) was partly fortuitous.* It rested on the κ_E,g defect above: an absorption opacity 3–50× too high
   in the LW and ionising groups raises the thermal source at the group photosphere.
2. *One albedo per group.* ε = κ_E/χ in the photospheric zone is 0.22 in 2b (0.67 in 2a): surface source factor √ε/(1 + √ε) ≈ 0.32
   against 0.45. Frequency-resolved transfer lets LW photons out through the windows between the Lyman wings where a single group mean
   cannot. The sensitivity to scattering is large in IonFS itself: without Rayleigh/Raman (`IFSRAM=0`, `mb/ifs06n`) its L_LW/L is
   1.4–1.6e-8, ×3.9–4.5 above the result with scattering (L_ion unchanged).
3. *The group photospheres lie in the unresolved temperature jump* between zone 26 (7525 K) and zone 27 (4861 K): B_LW changes by
   e^{−10.5} across one interface; d ln B/d ln T = 19 (LW), 23.6 (ionising). A factor 2.7 is 5.2 % in T (390 K of the 2660 K jump), a
   factor 5–6 is 7.6 %.

General remedy (not built): opacity-binned sub-groups (window/wing) in the LW band, which give the band two albedos; not a tuned
floor. What the shortfall means for the science at this epoch: nothing (T_eff ≈ 4500–4900 K; L_LW/L ~ 1e-9, and §16 shows the envelope
absorbs it within one zone). It matters when T_eff rises; the gate is repeated with the final configuration in §20.

## 18. Stage 2c: photo-rates and photo-heating from the group fields inside the solve (672c2fd)

Branch `mg2b`. With `MGROUPS` the chemistry's radiation field **is** the group fields at the Newton iterate. The grey star-light paths
(STARLIGHT's dilution W, SLWFLUX, the SLIMP variants, SLSHIELD, LOCALJ) no longer act in an MGROUPS run; they are still in the code for the
grey stack and are listed as removable in §18.8. Grey (MGROUPS unset) stays bit-identical to 4f1fe93 (§18.6).

### 18.1 What the solve sees

* **Star light of a zone** = c u_g = c(E_g − E_cmb b_g(T_cmb)) of every group (`mg_cu`; floor MGCUFL = 1e-60 erg cm⁻² s⁻¹ so that ratios
  exist where a group has underflowed; 1e-240 gave 0·∞ in the rate partials). No free-streaming floor, no cap at W = 1, no W = 1 rule for
  thermalised zones: the group energies carry all of that. The in-group spectral shape is B_ν(T_c) with the colour temperature of the
  SLBLEND/SLEB rule at the iterate (`mg_tc`).
* **Local variables of a coupled zone** (`cpl_zone_qg`, `cpl_g`, `cpl_tqg`): (ρ, T, x_H+, x_H2, q_1 … q_G, s, t_c) with q_g = c u_g / c u_g(start
  of step), s the LW shielding factor (§18.2) and t_c = T_c/T_c(start). All are functions of the iterate's unknowns, so the Jacobian has the
  derivatives of every photo-rate with respect to every group energy of the zone. Nothing is lagged.
* **Photo list** (H⁻ detachment, H₂⁺, He, He⁺, Li, the excited levels): `ChemNet.photorates_mgq`: R = R_static(CMB) + Σ_levels (g_l/Z)
  e^{−E_l/kT} Σ_g c u_g Q_{l,g}(T_c), Q_{l,g} = ∫_g σ_l B_ν(T_c)/(hν) dν / ∫_g B_ν(T_c) dν tabulated in T_c per group (`lj_gtables`; in the
  limit c u_g = 4σT⁴ b_g(T) the sum is the blackbody rate: thermalised zones k ≤ 24 give group/grey = 0.994–1.001, §18.5). The stimulated
  part of the Milne recombination takes the dilution of the group that holds the threshold.
* **Three-level H atom**: Balmer-continuum photoionisation β and Lyα pumping P12 from the dilution of the groups that hold those bands
  (`MGROLE`), at the iterate.
* **H(1s) ionisation and ground-state recombination**: the photon-conserving pair of §18.3 (two network reactions, rates at the iterate's T).
* **H₂ Lyman–Werner dissociation**: the photon-conserving closure of §18.2.
* **Photo-heating**: there is no separate heating term. The gas row is the total-energy row; each group row removes the absorbed energy from
  the group and returns the emitted energy, and the chemistry carries the binding energies. `lwheat!` (the grey stack's explicit LW heating)
  is off under MGROUPS; the LW group's absorbed energy is split into the share that stays in the gas (binding + 0.4 eV kinetic) and the
  fluorescence share that goes to the group below the LW band.
* **NLOIT opacity composition** is taken at the iterate's group fields (the same `cpl_tqg`), so opacity, emissivity and rates of a zone are
  functions of one state.
* `MGROUPS` now requires `SLIMP=2` (the chemistry's radiation at the iterate) and stops otherwise.

### 18.2 H₂ Lyman–Werner band: one closure for the rate and the group's energy

Sources (each read from the paper, `mb/lit/`): rate per molecule k = 1.1e8 F_ν(12.87 eV) s⁻¹ with F_ν = 4πJ_ν (Abel, Anninos, Zhang & Norman
1997, 1997NewA....2..181A, their k27); shielding f_sh(N, T) = 0.965/(1 + x/b_5)^1.1 + 0.035 e^{−8.5e-4 √(1+x)}/√(1+x), x = N/5e14 cm⁻²
(Wolcott-Green, Haiman & Bryan 2011, 2011MNRAS.418..838W, eq. 10: Draine & Bertoldi 1996, 1996ApJ...468..269D, eq. 37 with exponent 1.1);
dissociation probability per absorption ≈ 0.15 and ⟨hν⟩ = 12.4 eV (DB96).

The brief asked for σ_eff = −d ln f_sh/dN in the LW group's energy opacity. That form does not conserve photons: the energy a column
removes must equal the dissociations it causes. *Deviation, reported:* with J(N; T) = ∫_N^∞ f_sh dN′ (closed form, `mg_fshJ`) the band
energy that an infinite column can still use is ∝ J, so

    d ln(c u_LW)/dN = −f_sh/J ≡ −σ_eff,      c u(N) = c u(0) J(N)/J(0),      k = 1.1e8 (c u φ) f_sh J(0)/J ≡ k_unshielded(c u) · s

with φ = B_ν(12.87 eV, T_c)/∫_LW B_ν dν. The group value is the attenuated field and the rate is computed from it; every dissociation removes
e_d = Δν_LW/(1.1e8 J(0)) of band energy and nothing else removes any. κ_H2 = x_H2/μ_H · φΔν f_sh/J enters the LW group's energy opacity
(and χ); the share f_gas = min(1, 4.88 eV/e_d) of the absorbed energy stays in the gas, the rest (fluorescence) is emitted into the group
below the LW band. The column N is the smaller of the static column from the photosphere and the Sobolev column n_H2 b/|du/dr|; only the
zone's own x_H2 carries derivatives (the column is non-local: an inexact Newton column, exact residual).

Limiting cases (`mb/p2c/lwtest.jl`):

| check | result |
|---|---|
| −dJ/dN ÷ f_sh (finite difference) | 1.000000 |
| σ_eff(N → 0) | 2.19e-17 (100 K), 1.82e-17 (1000 K), 1.24e-17 cm² (8000 K); 1.0e-27 at N = 1e26 |
| s = f_sh J(0)/J at 100 K | 0.29 (N = 1e15), 5.2e-3 (1e17), 2.6e-5 (1e21), 4.6e-11 (1e26) |
| e_d, f_gas | 75.3 eV, 0.065 (30 K); 72.0, 0.068 (100 K); 59.6, 0.082 (1000 K); 50.3, 0.097 (3000 K); 40.8, 0.120 (8000 K) |
| DB96: ⟨hν⟩/p_diss | 12.4/0.15 = 83 eV: e_d of the closure is 0.5–0.9 of it (the fit's J(0) depends on b) |
| slab, T = 1000 K: energy removed per dissociation | 59.619 eV = e_d (identity) |
| slab attenuation, zone-centred/zone-edge vs exact J(N)/J(0) | N = 5e20: 0.443/0.405 vs 0.390; N = 1e24: 0.0786/0.0441 vs 0.0329 |

Caveat: the Sobolev column needs the velocity gradient of the iterate; in a static atmosphere it is the static column. WG11's fit is for
a static slab; no Doppler-shift correction of the shielding beyond the Sobolev cut is made.

### 18.3 H(1s): a photon-conserving ionisation / ground-state-recombination pair

The group rows exchange *energy* with H(1s): absorption c E_g n_1 S1/SB and Milne emission 4π n_e n_p λ_e³ s1_g(T) (σT⁴/π) X_g(T). The network
needs photon *numbers*. With the sums over the opacity grid's nodes of group g (σ_1s the grid's own cross-section, stimulated factor of the
gas temperature, in-group shape B_ν(T_1s), T_1s = max(T_c, MGTMIN)),

    S1 = Σ σ_1s (1 − e^{−hν/kT}) B_ν w,   SN = Σ σ_1s (1 − e^{−hν/kT}) B_ν w/(hν),   SB = Σ B_ν w,   ε̄_g = S1/SN,

both sides are the row's energy terms divided by the same mean photon energy ε̄_g:

    k_ion = Σ_{g ≥ H} c u_g SN/SB,          α_1 = Σ_{g ≥ H} 4π λ_e³ s1_g(T) (σT⁴/π) X_g(T) SN/S1

(`mg_kion`, `mg_arec`, reactions `mg_h1s_reactions` appended to the three-level network, switched per zone by `env.mgq.h1`: star-lit, not
thermalised, not TL3DEEP-deep). Then (ionisations − ground-state recombinations) × ε̄_g is the group's net energy exchange with H(1s),
identically. Where the group is opaque and its row is in on-the-spot balance the pair cancels exactly whatever T and T_c are, and each net
ionisation costs the transported field ε̄_g. It replaces the first 2c form (photon count with B(T_c) on the grid over the analytic band
fraction; recombination coefficient frozen at the start-of-step T), whose two sides were counted with different weights.

Limiting cases (`mb/p2c/h1test.jl`; 138 grid nodes in the ionising group; network H(1s) threshold 13.598 eV):

| check | result |
|---|---|
| α_1(T) at T_c = T against the Milne integral of the network's σ_1s (fine quadrature) | 0.980 (3000 K), 0.987 (5000), 0.992 (8000), 0.9936 (1e4), 0.9968 (2e4), 0.9988 (5e4), 0.9997 (1e5) |
| the same against the network's case A − case B | 1.018, 1.0012, 0.9977, 0.9973, 0.9976, 0.9985, 0.9991; α_1(1e4 K) = 1.580e-13 cm³ s⁻¹ |
| detailed balance in J = B_ν(T): k_ion λ_e³ e^{χ/kT}/α_1 | 1.00000000 for 3000–1e5 K |
| star light T_c on gas T: k_ion per c u against fine quadrature | 0.9869 (300 K, 5000 K), 0.9912 (5000, 7500), 0.9980 (8000, 3e4), 0.9998 (1e4, 1e5) |
| mean photon energy ε̄, mean cross-section | 13.98, 14.19, 15.91, 19.92 eV; 5.63e-18, 5.44e-18, 3.62e-18, 1.10e-18 cm² |
| closure: (absorbed − emitted)/ε̄ against the network's net rate, four states | relative difference ≤ 1.7e-16 |
| d ln α_1/d ln T (ForwardDiff) | −0.56 (5000 K), −0.67 (1e4 K) |

The 1–2 % at 3000–5000 K is the opacity grid's quadrature of the steep Wien tail (§17); it is the same on both sides of the pair. The sums
are scaled by e^{hν_L/kT_s} (ForwardDiff's quotient rule divided by an underflowing square otherwise: NaN at step 0, first attempt of t_mb14).

Approximation: one photon energy per group. The gross α_1 differs from the Milne value by ε̄_em(T)/ε̄_g(T_1s) (a few % for |T − T_c| ~ 1e3 K);
this cancels on the spot. He and He⁺ ionisation go through the photo list and have no such pair: recombinations to the He ground states
do not emit into the He groups (case B for He, as in the grey stack).

### 18.4 The first-step transient, and its cause: the opacity state was not restored at restart

The first 2c runs produced a burst in the first step after the restart. Record of what was seen and what it was:

| run | tree | first accepted Δt | zones 27–31 after the first steps |
|---|---|---|---|
| t_mb11 | 2c, frozen α_1 | 3427 s | T = 242–300 K, x_H2 0.28–0.47, W_ion 3.9e4; photosphere index 30 → 34, L 2.4e37 → 1.6e37 erg/s |
| t_mb12 | same, first Δt capped at 25 s | 25 s | zone 27: T 3278 K, x_H+ 2.2e-3, W_ion 3.1e4; F_ion/(cE_ion) = 1.3–12 in zones 27–33; 12 steps, 5 Newton failures |
| t_mb10 | same, 120 steps | | no recovery: 32 Newton failures; limiting variables x_H2 39, E_opt 33, x_H+ 23, E_LW 19 |
| t_mb13 | same, first Δt ≤ 0.05 s | | Newton diverges down to 7.6e-7 s |
| t_mb14 | + the pair of §18.3 | 857 s | zone 27: T 2295 K, x_H2 0.23, x_H+ 1.0e-3, W_ion 4.0e4: the burst persists; 30 steps, 8 Newton failures, L rms 60 % |
| t_mb17 | + restart fix | 857 s | zone 27: T 4724 K, x_H+ 8.2e-5, x_H2 2.8e-5, W_ion 6.9, W_LW 7.6: no burst; L rms 0.71 % over 30 steps |

Diagnosis (temporary print of the group formal solution's inputs in the first calls, t_mb14): in the first step the ionising group had
χ = 531, 121, 44.5, 0.86, 1.4e-3, 6.5e-4 cm² g⁻¹ in zones 26–31 — the *grey* opacity; the true n_1 σ_1s/ρ is 2.6e6 — K_1 = 0, tie weights
0.5–1 and h = 0.76–0.87 from zone 30 out: the ionising and LW groups were transparent and streaming, in an atmosphere where their mean
free path is 1e-9 of a zone. Cause: `nlo_it` returns the table opacity and a default group record (every group = the grey opacity)
for a zone without an NLOIT state `NLOST[k]`. That state is written by the end-of-step chemistry and, at a restart, by `nlost_restart!`,
which restored it only where `NLOPC[k] !== nothing` — and under NLOIT nothing ever assigns NLOPC. So after **every** restart the first
step ran all network zones on table opacities. With the group fields feeding the chemistry inside the solve, a transparent ionising group
(E_ion carried from the photosphere at c) ionises and dissociates the atmosphere in one step.

Fix (MGROUPS only): `nlost_restart!` restores the state of every network zone that has a three-level record with its chemistry part; the
opacity cache is reset when the groups are installed (the restart evaluates the opacities once before the groups exist; that evaluation
takes the grey means). With it the first group formal solution has χ_ion = 2.0–2.7e6 cm² g⁻¹ and tie weight 0 in zones 26–33 from call 0.

**The same defect is in the grey stack and in production** (the first step after every restart, i.e. every SEG = 1000 steps of z100mm,
uses table opacities in the network zones). It is not changed here because the grey path must stay bit-identical in this round; the fix
is the one condition in `PopIIIMaterial.nlost_restart!` (drop the `NLOPC[][k] === nothing` test). Reported to the PI.

What is *not* separated: the pair of §18.3 without the restart fix does not remove the burst (t_mb14); pair + restart fix does (t_mb17);
the restart fix with the first 2c rate form was not run. The earlier explanation of the burst by a loop gain ≠ 1 of the on-the-spot
cycle (the rate form's inconsistency is real: 3–10 % between the two sides) is therefore not shown to have been its cause.

### 18.5 Group-field rates against the paths they retire

One record of t_mb17 (step 14429, settled; `mb/p2c/ratecmp.jl`). Grey path: W = max(c(E − E_cmb), F)/4σT_c⁴ capped at 1 (STARLIGHT with
SLWFLUX) times B_ν(T_c), LW with `SLSHIELD=sobolev`. Group path: as above. W_g/W is the group's dilution over the grey one.

| zones | T [K] | W_g/W: IR, optical, LW, ionising | H⁻ detachment group/grey | H(1s) ionisation [s⁻¹] group / grey | H₂ dissociation [s⁻¹] group / grey |
|---|---|---|---|---|---|
| 15–24 (thermalised) | 7.6e4–1.1e4 | 1, 1, 1, 1 | 0.994–1.001 | 8.148e6 / 8.148e6 … 7.466e2 / 7.471e2 | — |
| 27 (first zone above the photosphere) | 4799 | 0.93, 1.0, 6.2, 5.8 | 1.012 | 1.6e-5 / 2.8e-6 | 0 / 0 |
| 30–39 | 2060–1630 | 0.81–0.85, 0.95–1.0, 4e-4, 0.7–1.2e-4 | 0.95–1.04 | 2.4e-13 … 3.2e-15 / 3.5e-9 … 2.6e-11 | 5.9e-22 … 1.2e-23 / 5e-31 … 3e-27 |
| 42–66 | 1461–673 | 0.66–1.2, 0.87–1.0, 1.4e-5–4.2, 7e-7–2e-4 | 0.87–1.00 | 6e-17 … 9e-21 / 3.7e-12 … 4.8e-17 | 8e-25 … 3e-20 / 1e-32 … 1e-69 |
| ≥ 69 (grey E_r < E_cmb) | 546–190 | optical 5e-2 → 1e-13 | 0.05 → 0 | 1e-20 / 3e-18 … 1e-15 | 1e-24 … 1e-45 / 3e-24 … 4e-18 |

Reading: (i) where the field is thermal the group rates equal the grey ones to the tables' accuracy. (ii) H⁻ detachment in the atmosphere
agrees to 5–13 %: the optical group is the grey field there. (iii) The grey path assumed the *undiluted-shape* Planck tail of T_c in the
ionising and LW bands (W of the bolometric field): it overestimates the H-ionising flux in the atmosphere by 4–7 orders of magnitude
against the transported field (W_ion/W ~ 1e-4–1e-6: the photons are absorbed in the first zone) and puts W = 1 where zone 27 really has
W_ion ≈ 6 (the gas is hotter than the radiation's colour). (iv) H₂ dissociation: both are negligible at this epoch (≤ 3e-20 s⁻¹); the group
field gives a small, transported LW field with the photon-conserving shielding, the grey path the Sobolev-shielded Planck tail. (v) In the
zones where the grey solution has E_r < E_cmb the grey path's SLWFLUX floor (W from the flux) kept a star-light rate (H⁻ 1e-6–0.1 s⁻¹)
that the group fields do not have: the optical group there is the draining restart reservoir of §14, not star light. This is the
outer-cloud reservoir defect that the energy-formulation work owns; from the IC it does not arise.

### 18.6 Behaviour

200 steps from the checkpoint with 672c2fd (`t_mb18`, MGROUPS=4) against stage 2b (`t_mb06`):

| | 2b `t_mb06` | 2c `t_mb18` |
|---|---|---|
| Newton failures per 50 steps (sum) | 31, 9, 10, 9 (59) | 12, 29, 32, 27 (100) |
| L(1.5 r_ph): median, rms | 6993 L☉, 0.47 % (20 records) | 6219 L☉, 0.32 % (4 records) |
| r_ph median | 115.64 R☉ | 115.68 R☉ |
| T_eff(L, r_ph) | 4908 K | 4766 K |
| Δt median of the records | 7.2e3 s | 1.07e4 s |
| wall time per step, second hundred | 2.96 s | 10.75 s |
| allocation per 100 steps | 430 GB | 1852 GB |

The first-step burst is gone (12 failures in the first 50 steps against 31), the photosphere and the luminosity are as quiet as in 2b.
The two numbers that are worse have causes outside the physics of this stage; both were found on this run and are repaired in the next
commit (§19.3, §20.1):

* *Cost.* The H(1s) pair of §18.3 summed over the opacity grid through an untyped reference to it (`G2Opacity.FAST` is a `Ref{Any}`):
  137 µs per rate, in every network evaluation of every star-lit zone. Behind a function barrier and over the group's own node range
  the same sum takes 2.9 µs (§19.3).
* *Newton failures.* The error estimate of the time step was zero in 8 of this run's 10 trace lines (`err` 1e-16 … 1e-13 against
  0.005–0.014 when it works): without an estimate the controller raises Δt five-fold after every accepted step until the Newton
  iteration fails, here every second step. Cause and repair: §20.1. Stage 2b had the same defect in 2 of its 10 trace lines.

Six groups (`t_mb20`, 672c2fd, `MGROUPS=6`, edges 0.7551, 11.2, 13.59, 24.59, 54.42 eV; stopped by hand after the record of step 14430
to free its cores): runs; 18 Newton failures in its ≈ 35 steps (four groups, `t_mb17`: 8 in 30); L(1.5 r_ph) 6137 L☉, r_ph 116.1 R☉,
T_eff 4741 K. At a 4540 K photosphere the helium groups hold nothing: E_g/E_r = 6.6e-23 (He I) and 4.4e-55 (He II) in zone 27, He II
at the floor (1e-150) in the outer cloud. The six-group configuration is exercised by this run, not validated: that needs a photosphere
above ≈ 2e4 K.

Grey bit-identity of the 2c tree: `t_mb19` (MGROUPS unset, 30 steps) tstate md5 246c3b436a9e423b0c6fafc6710cfec4, aux d540837fc3007163e0317e84e6cc0c24: identical to the reference (`t_mb15`, `t_mb05`).

### 18.7 Limitations of this stage (all stated, none tuned)

1. Group-field rates act in the zones that have a star-light record (`SLSTAR`, outside the τ = 2/3 zone of the step's setup); LW shielding
   columns start at the start-of-step photosphere index, and the colour temperature rule is the grey stack's SLBLEND/SLEB. Lyα pumping and
   the Balmer-continuum rate are exponentially sensitive to T_c (the in-group shape).
2. The column of the LW closure enters the Jacobian as a value; zones whose species are not unknowns use the start-of-step x_H2.
3. Photons of the LW and ionising groups absorbed by the n ≥ 2 bound–free continuum are removed as energy but are not counted as
   ionisations of those levels beyond the three-level atom's Balmer term. All LW fluorescence goes to the group below the band; the LW
   group has no H₂ emission term.
4. He: no ground-state recombination emission into the He groups.
5. Table (equilibrium) zones carry the grey opacity in every group (default record). They are optically thick in every group.
6. The emission mean switches its weight from B_ν(T) to B_ν(max(T, MGTMIN)) where the group's Planck sum underflows (T ≲ 255 K for the
   ionising group): a few-% discontinuity of a term that is ~1e-200 there. Below MGTMIN = 3000 K the H(1s) absorption shape is B_ν(3000 K);
   the row's LTE balance for the ionising group is then off by a few % at T < 3000 K, where b_ion ~ e^{−52}.
7. The restart partition of a grey checkpoint is far from the on-the-spot balance in the first step (first Δt 857 s, then growth).

### 18.8 Removable after stage 2c (not deleted in this round: outside the multigroup files)

With MGROUPS as the production path these grey star-light parameters have no function; evidence in §18.5:

| switch | what it was | status under MGROUPS |
|---|---|---|
| `STARLIGHT` dilution W, `SL3`, `SLBLEND`, `SLEB` | bolometric W × B_ν(T_c) for every band | W unused; T_c rule still used for the in-group shape (keep `SLBLEND/SLEB` until the shape comes from the groups) |
| `SLWFLUX` | W from max(c(E − E_cmb), F) | unused (group energies) |
| `SLIMP=0,1` | lagged / partially implicit W | unused; `SLIMP=2` is required, becomes the only path |
| `SLSHIELD` (`sobolev`, …), `lwheat!` | LW shielding and heating of the grey path | unused (§18.2) |
| `LOCALJ`, `LOCALJWMIN` | local J for the photo list | unused (`photorates_mgq`); the single `ljW` per zone was one dilution for all bands |
| `NLOPC` test in `nlost_restart!` | restart of the NLOIT state | defect (§18.4): drop the test |
| `W.cpQ`, `W.cpY` | workspace of an earlier coupling | stale (no reader) |

## 19. Switch audit: one code path (8076ce0)

Every switch that the multigroup stages introduced was read and given one of three verdicts: it is the only behaviour (switch and dead
branch deleted), it is a documented physical or numerical parameter, or it was a failed experiment and is gone. After this commit the
multigroup files (`MGOpac.jl`, `MGRad.jl`, the `mg_*` functions of `popiii_m2c.jl`) read exactly one environment variable, `MGROUPS`.

### 19.1 Verdicts

| switch | verdict | evidence | what the code does now |
|---|---|---|---|
| `MGROUPS` (unset or 0) | stays: the one switch | grey bit-identity (19.3) | no multigroup call is made |
| `MGROUPS=4`, `=6`, `=e1,e2,…` | stay | 4: every run of §10–§20; 6: `t_mb20` (§18.6) | 4 = IR, optical, LW, ionising; 6 = + He I, He II edges; an edge list in eV is checked by `mg_roles!`: it must contain 11.2 eV and the H edge, with a group below the LW band |
| `MGROUPS=3` | removed | no separate ionising group: the LW closure (§18.2) and the H(1s) pair (§18.3) each need a group of their own | error message |
| `MGROUPS=7`, `=8` | removed | never run; Balmer and Paschen edges can be given as an edge list | error message |
| `MGVEF=0` | removed (`MGVEF=1` is the only behaviour) | test 6 (§10) and §13: the minor groups need the factors and tie weights of their own optical depth | per-group formal solutions always; `MGVEFOK` only says whether the first one has run |
| `MGSOLVER=penta` | removed (band is the only solver) | test 7 (§10): `t_mg23` diverged in its first step; the block solver pivots only inside a block | `JuliaSolver.matslv!` always |
| `MGTMIN` | numerical parameter, now a constant (3000 K) | §4; limitation 6 of §18.7 | `const MGTMIN = 3000.0` |
| `MGMILNE`, `MGN2BF` | only behaviour since 672c2fd | §13, §14 | no switch |
| `MGMILNET`, `mg_ratios!` (lagged group ratios) | removed in stage 2b (e88345e) | §17 | group means at the Newton iterate |
| `TIESMOOTH=ε` (smooth ramp for SLABTIE's max(h, 0)) | removed: no effect | test 9 (§10): `t_mg22`, Newton failures 18, 19, 18 per 100 steps against 18, 19, 19 | `st_hdonor` has the one form (the eight deleted lines in `TitanResidFast.jl` are the switch and its ramp; the row is unchanged) |
| restart-partition variants (excess at T_CMB; continuous split capped at L/4πr²c) | were removed as code in e8d0b3b | §14: both failed the first step (`t_mg40`–`t_mg43`) | one partition; comment points here |
| "Newton in E_g" update with `MGCFL = 1e-12` | was reverted in stage 1 | §10 (`t_mg17`, `t_mg18`) | unknown is ln E_g; the comment block is replaced by a pointer |
| `MGCLAMP`, `MGCFL`, `MGEFL` as `Ref`s | constants | never set anywhere | `const` (19.2) |
| `MGNEWTON` | deleted | no reader, no writer | — |
| `EDIAG=1`, `=2` | stays: passive diagnostic | trajectories with and without it are bit-identical (test 4, §10) | total-energy rows of every converged stage → `ediag.csv` (and per zone → `ediagz.jls`) |
| `HPON`, `MGHPLUS` | stays: passive diagnostic | needed for the outgoing half-moment flux of the IonFS gate (§20) | stored by the group formal solutions, written to the record |
| record `mg_rec.jls` (`SNAPEVERY`) | stays | every analysis of §13–§20 reads it | one record per snapshot |

`refilter_mg!` (`Trbdf2.jl`, §20.1) is not a switch: with MGROUPS the second filter of the error estimate takes the exact image of
the estimate; the grey path keeps its finite difference, bit-identical.

Sections 5, 10, 11, 13 and 14 above describe `TIESMOOTH`, `MGSOLVER`, `MGVEF`, `MGMILNE` and `MGN2BF` as switches: that is the history; the
switches no longer exist.

### 19.2 Numerical parameters that remain (all constants, none read from the environment)

| constant | value | role | what happens beyond it |
|---|---|---|---|
| `MGCLAMP` | 3 | largest correction of ln E_g per Newton iteration | larger steps of a group that is 1e-20 of E_r overshoot (§10) |
| `MGEMIN` | 1e-150 | floor of E_g/E_r after an update | an underflowed group makes the matrix singular |
| `MGCFL` | 1e-3 | weight E_g/(E_g + MGCFL·E_r) of a group's correction in the convergence norm | a group below 1e-3 of the total does not decide convergence |
| `MGEFL` | 1e-3 | the same floor in the TR-BDF2 error norm | a group below 1e-3 of the total does not set Δt |
| `MGBFLOOR` | 1e-200 | floor of a group's Planck fraction | keeps ln E_g finite at the initial partition |
| `MGCUFL` | 1e-60 erg cm⁻² s⁻¹ | floor of c u_g in the ratios q_g = c u_g/c u_g(start) | 1e-240 gave 0·∞ (§18) |
| `MGTMIN` | 3000 K | floor of the in-group shape temperature of the absorption means | below it the ionising group's Planck sum underflows on the grid |
| `G2Opacity.LWSTRIDE` | 32 | every 32nd node of the LW band's Rayleigh-wing grid | §16 |

`MGCFL` and `MGEFL` are the two that touch the science: a transition in a group that carries less than 1e-3 of the local radiation energy
(the ionising and LW groups outside the photosphere, today) neither blocks convergence nor limits the step through the *radiation*
norms. It still does both through the species rows: the photo-rates feed x_H+ and x_H2, which are unknowns with their own tolerance, and
those are what limit Δt in the runs of §18 and §20.

### 19.3 Verification

* *The audit changes nothing with MGROUPS=4.* `t_mb21` (audited tree, 30 steps from the checkpoint) ends on the same state as `t_mb17`
  (672c2fd): aux md5 3d3839d04242a4e40c7cf59aaf5d68bf, uncompressed tstate md5 d82e3f6fb12bea2eb25d0cf770f9d518 in both (the `.gz`
  files differ in their headers only).
* *Grey bit-identity* (MGROUPS unset, 30 steps, `t_mb24`, tree of this commit): aux md5 d540837fc3007163e0317e84e6cc0c24, uncompressed tstate md5
  246c3b436a9e423b0c6fafc6710cfec4: identical to the reference (`t_mb19`, `t_mb15`, `t_mb05`).
* *Cost of the H(1s) pair.* `G2Opacity.FAST` is a `Ref{Any}`: the sums of `mg_h1s` over the grid's 484 nodes ran through dynamic
  dispatch, 137 µs per `mg_kion` call, in every network right-hand side and Jacobian of every star-lit zone; 672c2fd was 3.6 times
  slower per step than 2b for that reason (§18.6). The sums are now behind a function barrier (`_mg_h1s`) and run over the group's own
  contiguous node range (the grid ascends in ν; per group 210, 127, 9, 138 nodes): 2.9 µs per `mg_kion`, 5.7 µs per `mg_arec`
  (8.0 and 16.3 µs with 12 partials). Same terms in the same order: the table of §18.3 is reproduced to every printed digit and the
  30-step state above is bit-identical. Measured cost of the final code: §20.5.

### 19.4 Removable outside the multigroup files (listed, not deleted)

The table of §18.8, plus `SLIMPDBG` (`TitanResidFast.jl`, `PopIIIMaterial.jl`: debugging prints of the SLIMP=2 bring-up, off by default).
All of them are grey star-light parameters or leftovers that MGROUPS does not read; deleting them changes the grey stack and belongs with
the decision to make MGROUPS the production path. One more, found in §20.1: `CONV=1e-3` in the runs' `run.env` is read by
`popiii_m1.jl` and `s5_run.jl` only, not by `popiii_m2c.jl` — the production runs iterate to TITAN's default 1e-5 and the entry
suggests otherwise; remove it from the env files (or read it, which changes every run).

## 20. Validation of the final configuration (8076ce0, bce234e)

Tree `mg_trees/wip20` = commit 8076ce0 (§20.1–20.3, 20.5, 20.6); tree `mg_trees/wip25` = commit bce234e, which adds the line-cooling
repair of §20.4 (at this checkpoint it changes the solve by < 0.1 L☉ of line emission in zones 1–27 and the diagnostic by 640 L☉) and
with which the 500 steps were repeated (§20.2). All checkpoint runs: z100mm checkpoint of step 14400 (t = 237 841 yr, 100 zones, photosphere zone 26),
production physics, `MGROUPS=4` and nothing else, 8 threads.

### 20.1 The time-step error estimate vanished under MGROUPS: cause and repair

*Symptom.* Stage 2c as committed (672c2fd, `t_mb18`) had 100 Newton failures in 200 steps, and `err` of the accepted steps was
1e-16 … 1e-13 in 8 of the 10 trace lines. A diagnostic copy of the tree (`wip19d`, a print in `refilter!`; run `t_mb22`, 100 steps, every
step traced): the estimate is below 1e-9 in **53 of 100 steps**, from the first step on; 40 Newton failures.

*Cause.* TR-BDF2's error estimate is filtered twice (Hosea & Shampine 1996, Appl. Numer. Math. 20, 21, §6; `Trbdf2.jl`): Est =
(I − hdJ)⁻¹ est, then EST = (I − hdJ)⁻¹ Est. The second solve needs the image of Est in the conserved combinations, ΔQ = Q_y·Est. The
code took it by a finite difference, (Q(y + ε Est) − Q(y))/ε with ε = min(1, 1e-3/max|Est|), max over all unknowns. The group unknowns
are corrections to ln E_g. In the molecular cloud (zones 48–74 at this checkpoint) H₂ absorbs the LW group, and H the ionising group, down
to the floor 1e-150 E_r; there the first-filter estimate of δ ln E_g is 1e43 … 1e133 (a floor value has no meaningful logarithmic
error; the weighted norm gives these entries zero weight, so they are harmless *in the norm*). But they set ε = 1e-47 … 1e-137: y + ε Est
equals y in floating point, the difference is exactly zero, EST = 0, and the controller is told that the step has no error. It then
raises Δt five-fold per step until the Newton iteration fails, and quarters it. The same happened in 2 of the 10 trace lines of stage 2b
(`t_mb06`), i.e. the defect is as old as the logarithmic group update (stage 1); the Newton-failure counts of §10, §13 and §17 contain it.

*Repair* (`refilter_mg!`, used whenever MGROUPS is on): the image is taken exactly, with one evaluation of the residual's state in dual
numbers that carry Est as the single direction, ΔQ_k = ∂Q_k/∂y · Est. No step size, no cancellation; `combos` is generic in the number
type for this. This is the general form of what the finite difference approximated, not a guard on the group entries. The grey path
keeps its finite difference so that it stays bit-identical to production (`t_mb24`: 30 grey steps, final state identical to the
reference). Recommendation to the PI: let the grey path use the exact image as well (one code path; the grey trajectory then changes at
round-off level, which is why it was not done here).

*Result* (`t_mb23`, 200 steps, every step traced):

| | 672c2fd `t_mb18` | 8076ce0 `t_mb23` | grey (`t_mg13`) |
|---|---|---|---|
| steps with err < 1e-6 | 8 of 10 traced | **0 of 199** (first step after a restart has no estimate) | — |
| err of the accepted steps (p10, median, p90) | — | 0.0076, 0.0117, 0.018 | — |
| limiting variable | — | x_H2 192 steps (zone 28), x_H+ 7 | — |
| Newton failures per 50 steps | 12, 29, 32, 27 (100) | **10, 9, 9, 9 (37)** | 18, 19, 19 per 100 |
| Δt (p10, median, p90) | median 1.07e4 s | 5.7e3, 9.3e3, 2.0e4 s | — |
| L(1.5 r_ph) median, rms | 6219 L☉, — | 6208 L☉, 0.39 % | — |
| r_ph, T_eff(L, r_ph) | 115.68 R☉, 4766 K | 115.66 R☉, 4764 K | — |

The multigroup run now fails as often as the grey run (18.5 against 18–19 per 100 steps). The first-step burst of stages 2a/2b (31
failures in the first 50 steps of `t_mb06`) is gone as well: 10.


*What limits Δt now: the Newton iteration, in three separable ways.* The iteration stops at `dmax ≤ 1e-5` within 25 iterations
(`newton!`, `Trbdf2.jl`; the `CONV=1e-3` of the runs' `run.env` is not read by `popiii_m2c.jl` — only `popiii_m1.jl` and `s5_run.jl`
read it — so the runs iterate to TITAN's default 1e-5, as the failed iterations at 1.4e-5 below show; a dead setting, §19.4). The measure
is the base stack's max(|δr/r|, |δm/m|, |δρ/ρ|, |δT/T|, |δE_r/E_r|) plus the species' plus the groups' (`mg_dmax!`); u and F_r are not
in it. A copy tree with both parts printed separately (`wip20d`, runs `t_mb27`, `t_mb33`, `t_mb36`) gives:

1. *At the checkpoint: the base stack's photospheric stall, as in grey.* 100 steps of `t_mb27`: 17 iteration-limit failures; in all 17
   the base measure is the larger one and the group measure is 0.33–0.345 of it (e.g. 8.1e-5 / 2.8e-5 at Δt = 1.85e4 s; 1.4e-5 /
   4.9e-6). The stalled unknown is x_H+ of zone 26 (the photosphere) in 34 of 37 failures of the multigroup run and 74 of 75 of the
   grey run; the corrections oscillate at ≈ 3e-5 against the 1e-5 tolerance. The groups add nothing to it: 18.4 failures per 100
   steps against 18.75 grey (§20.2).
2. *At a restart from a grey checkpoint and at the first step from the initial conditions: an exactly singular matrix.* 2 failures in
   `t_mb27` and 1 in `t_mb33` are `LAPACKException(1766)` / `(1768)` from `gbtrf`: zero pivots in the F_3 and F_4 columns of the last
   ghost block. The matrix itself is well formed (dump in `t_mb35`: ghost rows −1, +1), but the groups are empty there: E_3 = 6.7e-158
   erg cm⁻³, nom(F_3) = 1.1e-50 against cE_3 = 2e-147, so within the same group's rows the ln E_g and F_g columns differ by ~100
   orders of magnitude (row scales 1.5e157 and 9.0e49), and an entry of column 1768 two blocks further in is 1.3e-36 of its row's
   maximum. Partial pivoting then meets pivots at round-off level or exactly zero. The step is repeated at Δt/4 and Δt/16 and
   then passes (1.37e4 → 3.4e3 → 8.6e2 s at the checkpoint), and the groups relax within the first steps.
   *Tried and not in the code* (negative result): at a grey restart the partition gives the empty groups F_g = F_r b_g(T), i.e.
   |F_g| ≫ cE_g by 97 orders; limiting every group's flux to c·½(E_g,k−1 + E_g,k) and giving the cut part to the group that
   holds the energy (copy tree `wip26`, run `t_mb40`) leaves the two failed attempts exactly as they were (divergence at the first
   iteration, the same 1.37e4 → 3.4e3 → 8.6e2 s; last corrections 1e13–1e17 against 1e9–1e12 without). The unrealisable fluxes are
   therefore not the cause; what remains is the emptiness itself (E_g 150 orders below E_r in the outer zones, where the grey
   field has no excess over the CMB). Cost: two rejected attempts once per restart from a grey checkpoint; open.
3. *From the initial conditions: the centre of the collapsing core, in grey and multigroup alike* (§20.6). In the early collapse
   (`t_mb33`, `t_mb36`, Δt = 1e10–1e11 s) 25 of 25 iteration-limit failures end with the group measure above the base measure, the
   base measure already below 1e-5 in 22 of them (2.3e-9 / 2.5e-5; 4.9e-11 / 1.6e-5); the worst group is always group 1 — the one
   that holds E_r — in zones 6–12, with δ ln E_1 = 1.2e-5 … 5.8e-5. The simplest reading (not verified iteration by
   iteration) is that this is not a separate failure mode: E_r is updated linearly and E_1 by exp(δ), so a correction δ of ln E_1
   leaves ½δ² in the linear sum row E_r − Σ_g E_g, which the next iteration removes from ln E_1 alone; a last correction of
   1e-5–1e-4 with E_r converged would then mean that the iteration before it still moved E_r by ~1e-2 — the iteration had not
   finished in 25 iterations and the logarithmic group variable is one iteration behind the base variable. Two repairs of the
   linear solve, tried on the assumption that it was solve error instead, changed nothing, and neither did 50 iterations instead
   of 25 (§20.6): the iteration stalls, it is not slow.

### 20.2 500 steps from the checkpoint

`t_mb26` (`MGROUPS=4 EDIAG=1 SNAPEVERY=10`, steps 14400–14898, 5.21e6 s of model time), against the grey run and stage 2b:

| | grey `t_mg13` (400 steps) | 2b `t_mb06` (200 steps) | final `t_mb26` (500 steps) |
|---|---|---|---|
| Newton failures per 100 steps | 18.75 | 29.5 (19 after the first 50 steps) | **18.4** (19, 18, 19, 18, 18) |
| steps rejected by the error estimate | — | — | 0 |
| variable with the largest error estimate | — | — | x_H2 in 470 steps (zone 28), x_H+ in 29 (err as in §20.1) |
| Δt median over the records | 1.15e4 s | 7.2e3 s | 1.15e4 s |
| Δt of the failed attempts (median; range) | 1.71e4 s; 0.58–2.97e4 | — | 2.58e4 s; 2.16–3.64e4 |
| L(1.5 r_ph) median; rms | 6507 L☉; 0.7 % | 6993 L☉; 0.47 % | **6233 L☉; 0.71 %** (6124–6336, 50 records) |
| r_ph median | 115.5–116.3 R☉ | 115.64 R☉ | 115.87 R☉ |
| T_eff(L(1.5 r_ph), r_ph) | — | 4908 K | 4764 K |
| wall time per step | 1.64–1.67 s | 2.96 s | 4.19 s (`t_mb23`; 4.35–4.64 s with EDIAG and 5× the records) |

*Repeated with the committed final code* (bce234e, `t_mb42`, same settings, 500 steps): 18.6 Newton failures per 100 steps (19, 18,
19, 19, 18, two of the 93 being the restart's first step, §20.1), L(1.5 r_ph) = 6238 L☉ with rms 0.70 % (6109–6333, 50 records),
r_ph = 115.91 R☉, T_eff = 4765 K, Δt median 1.18e4 s, 4.36–4.44 s per step, total-energy imbalance 3671 L☉ (§20.4).

No step was rejected, no run stopped, the luminosity and the photosphere are as quiet as in the grey run. The final configuration is
11 % less luminous than stage 2b at the same epoch and 4 % less than grey: with the photo-rates taken from the group fields the zones
just above the photosphere are cooler and molecular (step 14500, 2b → final: zone 27 T 4803 → 4559 K, x_H2 1.8e-4 → 2.1e-2; zone 28
3721 → 2767 K; zone 25 7579 → 8239 K). Which of the two is right cannot be decided on this grid (§20.3: zones 26/27/28 are at
4877/4559/2767 K, one zone per factor 1.7 in T).

### 20.3 IonFS gate: not met, and why it cannot discriminate at this epoch

Gate (task brief): on ≥ 4 snapshots the outgoing LW and ionising photon luminosities at the photosphere within a factor 2 of the
frequency-resolved formal solution (`ionfs_post_ram.jl`, static spherical Feautrier/ALI with scattering on the run's own structure; MG
side L⁺_g = 4πr²·4π·H⁺_g at the interface nearest r_ph, script `mg/gate.jl`). Five snapshots of `t_mb26`:

| step | IonFS L_LW/L | IonFS L_ion/L | MG L⁺_LW/L | MG/IonFS | MG L⁺_ion/L | MG/IonFS | IonFS L_LW/L without its sub-zone model | MG/IonFS |
|---|---|---|---|---|---|---|---|---|
| 14500 | 6.32e-13 | 1.28e-10 | 2.77e-10 | 439 | 4.61e-12 | 0.036 | 5.09e-9 | 0.054 |
| 14600 | 6.19e-13 | 2.30e-10 | 4.37e-10 | 705 | 1.09e-19 | 4.8e-10 | 8.50e-9 | 0.051 |
| 14700 | 6.23e-13 | 1.87e-10 | 3.74e-10 | 601 | 8.67e-20 | 4.6e-10 | 7.07e-9 | 0.053 |
| 14800 | 6.19e-13 | 1.26e-10 | 2.78e-10 | 449 | 5.80e-12 | 0.046 | 5.00e-9 | 0.056 |
| 14890 | 6.25e-13 | 1.31e-10 | 2.84e-10 | 454 | 4.67e-12 | 0.036 | 5.16e-9 | 0.055 |

**The gate is not met in either band on any snapshot** (stage 2b, §17: LW 0.36–0.42, ionising 0.16–0.21 — not met either). Nothing was
tuned to change that. The quantified reasons:

1. *The reference is undetermined to 4 dex in the LW band.* IonFS replaces the front zone by sub-points on a grey T(τ) =
   T_eff (¾(τ + ⅔))^¼ profile when that zone's half optical depth exceeds 1 (`if Lrun > 0 && τzc > 1`), and otherwise uses the zone
   temperature. In the 2b snapshots the branch was not taken (Q_LW = 3.7e39 s⁻¹); in the final-configuration snapshots it is, and the LW
   photons then leave from the τ → 0 skin at 0.84 T_eff ≈ 3640–3750 K: Q_LW = 5.3–5.5e35 s⁻¹. The same snapshots with the branch
   disabled (`mb/ifs/ionfs_post_ram_nosub.jl`, output `mb/ifs26n`): Q_LW = 4.4–7.3e39 s⁻¹, **8000–14000 times more**. The Wien-tail
   sensitivity d ln B_LW/d ln T is 19–33 across the temperature cliff of zones 26/27/28 (4877/4559/2767 K), which the grid does not
   resolve. The multigroup value lies between the two references (0.05× and 440–700×): L⁺_LW ≈ 2.3e38 photons/s, 0.4 of the blackbody
   half-flux at T(zone 27). A factor-2 gate against a reference with a 4-dex switch tests the switch, not the transport.
2. *The half-flux diagnostic of an opaque group at one interface is not a robust quantity.* The record's passive H⁺_g is not consistent
   with the moment solution it is derived from: in the LW group L⁺ (2.77e-10 L) is smaller than the net L_LW (3.36e-10 L) at r_ph, and
   cE_g/4 would give 7e-9 L (LW) and 2–4e-11 L (ionising). The ionising L⁺ at r_ph jumps between 4.6e-12 and 1e-19 of L from snapshot
   to snapshot while E_4/E_r in the adjoining zones is unchanged (5.4e-11 in zone 26, 2.7e-11 in zone 27): the interface nearest r_ph
   sits on either side of the zone where the group is absorbed. The gate quantity should be the net group luminosity through a fixed
   optical-depth surface, and the reference should come from a structure that resolves the photosphere; neither exists at 100 zones.
3. *Ionising band:* IonFS emits 9.3e37–1.6e38 photons/s from the photosphere model and all of them are absorbed in the first zone
   outside (Q_esc = 0, r_IF = r_ph in all snapshots); the multigroup net ionising luminosity at r_ph is 5–12e-20 L. Both say the same
   thing about what leaves; the ratio at the interface is dominated by point 2.

**What both calculations agree on, and what the science question needs:** nothing gets ahead of the photosphere at this epoch. The
multigroup net LW luminosity falls from 6.9e-3 L (zone 25) through 6.2e-6 (26), 3.4e-10 (27, r_ph), 5.8e-11 (28) to 7.4e-15 L (zone 29),
where x_H2 = 0.41–0.44; E_LW/E_r = 3.4e-8, 8.5e-9, 1.3e-12, 7.5e-14 in zones 26–29. The net ionising luminosity is 5.6e-20 L at r_ph
and 1.7e-25 L at zone 29. IonFS: Q_esc,H = 0, ionisation front at r_ph. The gate has to be repeated when the photosphere is hot enough
for the LW and ionising bands not to sit on the Wien tail of an unresolved cliff (T_ph ≳ 1e4 K), where the reference is well defined.
### 20.4 Conservation

**Sum identities** (all 50 records of `t_mb26`): max over zones of |Σ_g E_g/E_r − 1| = 1.1e-13 … 1.9e-12, of |Σ_g F_g − F_r|/|F_r|
= 2.1e-16 … 1.5e-15. The gas row's exchange is the sum of the group rows' exchange terms by construction (`Mg.xch`).

**Budget of each group row** (`mb/p2c/budget.jl`, step 14800 of `t_mb26`, zones 27–60, L☉; the script evaluates emission, absorption,
transfers, ΔL and dE/dt from the record and has no velocity-dependent terms, so the residual column is what advection and
compression work carry plus the record's rounding):

| group | emission | absorption | line cooling in | transfers | net source | ΔL across the shell | dE/dt | residual |
|---|---|---|---|---|---|---|---|---|
| 1 IR | 1030.8 | 692.5 | 633.1 | — | 971.4 | 971.2 | 2.1 | −1.9 |
| 2 optical | 7874.0 | 6199.7 | — | −3.8e-7 | 1674.3 | 1721.6 | −34.6 | −12.6 |
| 3 LW | 3.4e-12 | 2.24e-6 | — | +3.8e-7 | −1.86e-6 | −1.43e-6 | — | −4.3e-7 |
| 4 ionising | 415.14 | 415.12 | — | — | 0.021 | −2e-16 | — | 0.021 |
| Σ | 9319.9 | 7307.3 | 633.1 | 0 | 2645.7 | 2692.8 | | |

The LW group absorbs what flows in (1.44e-6 L☉ at zone 27) and hands 17 % of it to the group below as fluorescence; the ionising
group is in on-the-spot balance to 5e-5 (it emits and reabsorbs 415 L☉ in place and transports nothing).

**Total energy (EDIAG).** `EDIAG` evaluates the conservative total-energy row at every converged stage and reports what the rows the
solve actually used left over. Grey: 3741 L☉ (§12; bands: zones 1–26 −0.4, 27–50 +442, 51–87 +19, 88–100 +3280 — the last is u × the
radiation-momentum residual at the FLUXTIE interfaces and belongs to the energy-formulation work). Stage 1: 3754; 2a: 3801. The final
configuration at 8076ce0 gave **4315 L☉** (`t_mb26`, 490 steps, p10 4309, p90 4322), 574 more than grey, and the per-zone table
(`EDIAG=2`, `t_mb34` against grey `t_mg24`) put all of the difference in zones 36–57 (T = 900–1700 K, molecular infall):

| zone | 36 | 38 | 40 | 42 | 44 | 46 | 48 | 50 | 52 | 54 | 56 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| grey imbalance [L☉] | 38.7 | 13.4 | 6.1 | 3.9 | 2.6 | 1.9 | 1.5 | 1.1 | — | — | — |
| multigroup, 8076ce0 | 49.0 | 33.2 | 34.0 | 41.7 | 54.7 | 65.2 | 57.2 | 22.0 | 9.0 | 7.3 | 4.5 |
| line cooling emitted into the IR group, `Mg.lcl` ΔV | 10.4 | 19.7 | 28.1 | 38.0 | 52.3 | 63.3 | 55.6 | 21.0 | 8.5 | 6.7 | 4.2 |

The difference is the H₂ line cooling, zone by zone (Σ over zones 27–56: 621.5 L☉). *Cause:* stage 2b emits the coupled zones' line
cooling into the IR group, but only where the gas row is the internal-energy row; a row that holds gas + radiation (the conservative
total row, the `lrad` row without `GASROW`, and therefore EDIAG's all-zones total row) still subtracted the whole chemistry cooling
`cpL`, as in grey where the line photons leave the system. Two consequences: the diagnostic counted the line emission, which the solve
had moved from the gas into E_1, as lost; and — a real inconsistency — in the zones the solve closes with the total row (zones 1–26
here) the line photons were removed from the gas and not given to any group. At the checkpoint the second is harmless (line cooling
of zones 1–27: < 0.1 L☉); in a prestellar core closed with total rows it would not be.

*Repair:* the group row gets the line emission in every coupled zone, and a gas row that holds gas + radiation takes the same term
back (`mg_rows!` in `MGRad.jl`; one added term after the chemistry cooling in the gas row of `TitanResidFast.jl`, MGROUPS only —
*merge note for energy-compat:* a source term, no change of the row's discretisation). With it (`t_mb38`, `EDIAG=2`, 100 steps):

| band | grey `t_mg24` | multigroup 8076ce0 `t_mb34` | multigroup, repaired `t_mb38` |
|---|---|---|---|
| zones 1–26 | −0.4 | +0.4 | +0.4 |
| zones 27–50 | +442 | +1016 | **+438** |
| zones 51–87 | +19 | +85 | **+19** |
| zones 88–100 | +3280 | +3214 | +3213 |
| all | 3741 | 4315 | **3670** (per-step median 3676, p10 3671, p90 3681, 90 steps) |

Over 500 steps with the committed code (`t_mb42`): median 3671 L☉, p10 3665, p90 3677 (490 steps; 40 steps `t_mb29`: 3677).
The multigroup rows now leave the same total-energy imbalance as the grey rows, band by band: the group rows add no energy error of
their own, and what remains is the grey stack's (§12), owned by the energy-formulation branch.
### 20.5 Cost

Measured, 8 threads, 100 zones, from the checkpoint (protostar with a photosphere; 31 Newton iterations per step in all stacks):

| | grey `t_mg13` | stage 1 `t_mg14`, `t_mg16` | 2b `t_mb06` | 2c 672c2fd `t_mb18` | final 8076ce0 `t_mb23` |
|---|---|---|---|---|---|
| wall time per step | 1.64–1.67 s | 2.45–2.55 s | 2.96 s | 10.75 s | **4.19 s (×2.5 grey)** |
| allocation per 100 steps | 275–278 GB | — | 430 GB | 1852 GB | 717 GB |
| garbage collection per 100 steps | 24–25 s | — | — | — | 70 s |
| Jacobian assembly (`matgen`) per Newton iteration | 42.6 ms | 53.3 ms | — | — | 113 ms (`t_mb27`, CPLPROF) |
| `update` per Newton iteration | — | 7.5 ms | — | — | 18.7 ms |
| banded LU (n = 1768, bandwidth 50) per 50 steps | 0.6 s (block solver) | 16–18 s | — | — | not timed separately |

In the final code 85 % of a step is the Jacobian assembly and 13 % the update (which re-evaluates the group records at the new
iterate). The assembly went from 42.6 ms (grey) to 53.3 ms with the wider dual numbers of stage 1 (102 partials instead of 54) and to
113 ms with the group physics of 2b/2c evaluated at the iterate — the no-lag requirement: in 500 steps of `t_mb26` 1.9e7 evaluations of
the group opacity ratios (3562 thread-seconds, 7.1 per step) and 6206 thread-seconds of correlated-k expansions. The untyped-grid defect
that made 672c2fd cost 10.75 s is repaired (§19.3).

In the prestellar collapse from the initial conditions the multigroup run costs the same per step as the grey one (1.84–1.99 s against
1.96 s, §20.6): the group-field rates and the opacity ratios run only in zones with a record of the star's light.

What the plan's remedies (§7, §11) would buy, from these measurements — none was built:

* *Analytic group-row Jacobian or two-pass colouring* (the material sweep keeps 54 partials, the group rows get their own narrow
  pass): removes the width cost. Stage 1 measured that cost on the material alone: 53.3 − 42.6 = 10.7 ms per iteration, 0.33 s per
  step. Of the 60 ms per iteration that 2b/2c added, the part that is dual arithmetic at width 102 would shrink by up to 102/54; if all
  of it scaled, 0.9 s per step more. Bound: 0.3–1.2 s of the 4.19 s. It does not bring the cost to the plan's ×1.2 target by itself.
* *What would:* the opacity ratios and their correlated-k expansions at the iterate are the largest single item (≈ 20 thread-seconds
  per step, 2.5 s of wall time if perfectly parallel; the two timers of the previous paragraph may overlap, and part of the
  expansions is the grey opacity's). They are evaluated once per colour per iteration in dual numbers; with the
  group rows in a second pass they would be needed once per iteration in plain numbers plus their T, ρ, x derivatives (three more
  columns). That is the same restructuring as the first item and is where its real gain lies, but the number is an expectation, not
  a measurement.
* *Pivoted block solver in place of banded LAPACK:* the LU is 0.33 s per step (8 %); a block-pentadiagonal elimination with 17 × 17
  blocks scales as (17/9)³ of the grey solver's 0.012 s per step, ≈ 0.08 s.
* *Jacobian-free Newton–Krylov with the banded matrix as a lagged preconditioner:* replaces most assemblies by residual evaluations.
  Not measured; it changes the convergence behaviour of an iteration that already fails 18 times per 100 steps for reasons outside the
  group rows (§20.1), so it should not be tried before that is resolved.
* *Allocation:* 7.2 GB per step against 2.8 GB grey, 0.45 s per step more garbage collection; typed, preallocated group work arrays
  for the update would remove most of it.

For correctness nothing more was needed: the 2.5× is affordable for the 100-zone production runs (z100mm: 1e4 steps in 4.6 h grey,
11.6 h multigroup).
### 20.6 From the initial conditions: the multigroup run does not reach the protostar

Set-up as `docs/POPIII_STACK.md`: abn02red initial conditions, 100 zones, production physics (env files of `z100mm`), `MGROUPS=4` from
step 1, 8 threads. A grey control was run from the same tree and the same files, because `z100mm` is not a clean reference (it was
restarted through several code versions).

| | grey control `t_mb32` | multigroup `t_mb25` (8076ce0) |
|---|---|---|
| log n_c ≥ 8 | step 94, t = 2.2774e5 yr | step 147, 2.2827e5 yr |
| ≥ 9 | step 184, 2.3983e5 yr | step 346, 2.4028e5 yr |
| ≥ 9.5 | step 190 | step 451, 2.4188e5 yr |
| ≥ 10 | step 214, 2.4201e5 yr, Δt = 2.2e9 s | step 900, 2.4264e5 yr, Δt = 4.4e7 s |
| ≥ 22 (protostar) | **step 273**, 2.42760e5 yr, T_c = 4.6e4 K | **not reached**: stopped at step 1000, log n_c = 10.16, Δt = 1.9e7 s |
| Newton failures up to step 100 / 200 / 300 / 400 / 500 | 21 / 46 / 53 / 57 / 57 | 24 / 52 / 83 / 116 / 151 (333 at step 1000) |
| wall time per step | 1.2–2.2 s | 1.84–1.99 s |

**The multigroup stack does not pass this test.** Up to log n_c ≈ 8.5 the two runs fail equally often; then the grey run's failures
stop and it collapses to the protostar in 90 more steps, while the multigroup run keeps failing 33–37 times per 100 steps with Δt
held 50 times smaller at the same central density, and creeps. The run was stopped at step 1000 (0.1 dex in 100 steps); all runs of
this section are stopped.

What the failures are, as far as measured:

* In both runs the failing iterations end with the same pattern in the base unknowns (print of the last corrections per equation):
  δr/r = 1e-5 … 7e-4 and δm/m = 3 δr/r at zones 6–8, δE_r/E_r = 1e-4 … 5e-3 at zone 2, δF_r at zones 2–4, changing sign from one
  failure to the next. These are the innermost zones of the fixed logarithmic grid (Δr = 2–7e9 cm, optical depth per zone 1e-13,
  F/(cE) = 1e-11 … 1e-10, cΔt/Δr = 1e8 … 5e8 at Δt = 3e7 s): a free-streaming, isotropic radiation field crossed 1e8 times per
  step. The mode exists in the grey stack; the grey run gets through it.
* In the early collapse of the multigroup run (`t_mb33`, `t_mb36`) the failing iterations end with the base measure converged and
  ln E_1 of zones 6–12 still moving by 1e-5 … 1e-4: most simply the logarithmic group variable being one iteration behind
  (§20.1, item 3), i.e. the same slow iterations seen through the group measure.
* *Two repairs of the linear solve were tried and are not in the code* (negative results; the hypothesis was that the banded LU
  loses the central zones' corrections to round-off, cΔt/Δr being 1e8–1e11). (a) Column equilibration in addition to the row
  scaling (tree `wip22`, run `t_mb37`): failures 27 / 55 / 84 at steps 100 / 200 / 300, log n_c ≥ 8.5 at step 243 — as without
  (24 / 52 / 83; step 227). (b) Column equilibration plus two steps of iterative refinement (tree `wip24`, run `t_mb39`): the
  first solve is accurate to 1e-8 … 1e-10 of the largest correction and the refined one to 1e-14, and the failures are 22 / 51 /
  83 at steps 100 / 200 / 300, and step 300 is at t = 2.3926e5 yr with log n_c = 8.89, the same time as step 300 of `t_mb25`.
  **The linear solve is not the obstruction.**
* *The iteration limit is not the obstruction either* (`t_mb41`, bce234e, `NITER=50` instead of 25): failures 26 / 54 / 85 / 118
  at steps 100 / 200 / 300 / 400 (`t_mb25`: 24 / 52 / 83 / 116), log n_c ≥ 8.5 at step 189, 8.81 at step 300. After 50 iterations
  the last corrections are the same as after 25 (δr/r = −4e-5 and δm/m = −1.2e-4 at zone 7, δE_r/E_r = −3.4e-4 and δF_r at zone 2):
  the iteration does not converge slowly, it stalls in an oscillation of the centre zones.

**What is not known:** why the multigroup run keeps failing after log n_c ≈ 8.5 where the grey run stops failing. The failing
unknowns are the base stack's in both. Candidates that the measurements above do not separate: (i) the group pair (ln E_1,
F_1) duplicates (E_r, F_r) in zones where nothing absorbs or emits, and the pair of rows changes the Newton path of the
oscillating centre mode — the grey iteration leaves the oscillation, the multigroup one does not; (ii) the group tie and boundary
rows at the centre. (ii) is the energy-formulation branch's area and was not touched. The next measurement would be the per-
iteration corrections of one failing step in both stacks from the same state (a grey checkpoint at log n_c ≈ 9 restarted grey and
multigroup). Until this is resolved **a multigroup run has to start from a grey checkpoint after protostar formation** — the
configuration validated in §20.1–20.4 — and the science question at the earliest epoch (does LW or ionising light get ahead of
the photosphere during protostar formation) is answered only from step 14400 on.
