# Implicit multigroup radiation for the Pop III protostar (`MGROUPS`)

Branch `multigroup` (from `popiii-s6` 43269f6), 16 Nov 2026. Code: julia/MGRad.jl (included into TitanResidFast: unknowns, rows,
Newton hooks, Planck group fractions), julia/MGOpac.jl (included into PopIIIMaterial: group-opacity ratios, group photo-rates),
julia/G2Opacity.jl (`group_fast`, `group_index`), hooks in TitanResidFast.residual!, Trbdf2.jl, TitanFortran.jl (meqn), popiii_m2c.jl
(`mg_install!`, checkpoint, `mg_record`). Unit test julia/mg_test.jl. Test runs `D/t_mg1N`, D = /zpool/nvme/data/tabel_scratch/titan_popiii/runs,
from `git archive` trees in /zpool/nvme/data/tabel_scratch/titan_popiii/mg_trees, restarts of copies of z100mm's checkpoint.

Status of this note: §1–§9 are the design; §10 is what the first coupled prototype (stage 1) implements and what it measured; §11 the
plan of the remaining stages. Where a section says "stage 1", the prototype does exactly that; everything else is design.

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

## 13. Stage 2a: Eddington factors per group (`MGVEF=1`, default with MGROUPS)

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
