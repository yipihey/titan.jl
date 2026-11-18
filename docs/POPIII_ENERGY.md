# Energy conservation of the Pop III accretion stack: measurement, derivation, compatible rows

Branch `energy-compat` (from 4f1fe93, head of `multigroup`). Commits: 774c85b (zone_diag split), a37f48b (TIEFRAME), be38295 (ECOMPAT),
50f1609 (TIEFRAME for the multigroup rows), e0e6b2e (sub-zone kinetic-energy transport, TR-BDF2 kinetic reference); follow-up
(§14–§17): 39923e1 (ECGRAV), cd9922c (TIEX, ERESTOL), bb57e3e (GEQGP); round 3 (§18–§21): d676f9b (ERLAB, GFV), 6642493 (GVAR).
All runs: 100 zones (except the zone-number convergence tests of §21, at 50 and 200 zones; production stays at 100), restart of z100mm's checkpoint `z19z100mmS5_auto0` (step 17398, M* = 0.88 M☉ in perf.csv, photosphere zone 26, 111 R☉) unless
stated. Luminosities are 4π × (per-steradian value) / L☉. Sign convention of every table: a positive row residual means the row would
create energy; a negative term added to the gas row heats.

## 0. Answer

*Is a single energy-conserving formulation possible?* For internal + kinetic + radiation energy against the momentum row's forces:
yes, in every zone, and it does not need the E − K subtraction. **Including the gravitational binding energy: no, not to round-off
on this mesh** — the follow-up (§14–§17, summary in §14.7) measured and derived this; the bullets below are the first round's,
corrected where the follow-up changed them (marked "§14", "§15", "§16", "§17"). Round 3 (§18–§21: the lab-frame radiation row
alone, the gravity weight of the momentum row, ECGRAV=1 from the IC, zone-number convergence) is summarised in the last four
bullets of this section.

* The gas energy row can be the first law in **every** zone (star, shock, cold supersonic infall with K/e ~ 1e3–1e4) and still conserve
  internal + kinetic energy to round-off, if it is written as the exact discrete adjoint of the momentum row: the kinetic energy that the
  momentum row's discretisation does not keep is a sum of explicit products of velocity *differences* (never K_new − K_old), and it is
  returned to the gas as heat. On this mesh there are four such terms (§4): time discretisation, viscous-work form, upwinded momentum
  flux, and the kinetic energy of the TR-BDF2 reference state. Measured: per-zone identity residual 0.000 L☉ in zones 1–87
  (max over 1000 stages), |Σ| ≤ 1 L☉ in zones 88–100 in 90 % of the stages (round-off of the lab-frame radiation row there, §5).
* The compatible row books the gravity power the way the momentum row's force does it (Σ_j u_j g_j ½(m_j−1 + m_j)). The total-energy
  rows of the current stack book it as (mass flux r² u ρᵇ through the interface, ρᵇ the advected interface density) × g × Δr; the two
  differ per zone by G5 = +1.0e3 L☉ (first 0.2 yr of the restart) to +1.9e3 L☉ (time-weighted over 12.6 yr) in the star and
  +0.7…0.85e3 L☉ in zones 27–50: 25–30 % of L_em in all. It is the largest term of the dual-energy stack's ledger (§2) and the
  reason the star's structure depends on which zones are flagged (§7). *Corrected by §14:* the first round called the momentum
  row's form "the gravity power of this discretisation". It is not. An explicit binding energy W follows the flux form to 0.02–0.8 %
  and the momentum row's power is off by 2.6e3–3.2e3 L☉ (§14.2); the infall of every run arrives at the shock with up to 11 % more
  kinetic energy than free fall allows (§14.4). No state-only force conserves K + U + W + E_r to round-off on a mesh that both
  moves and lets mass through with an upwinded interface density (§14.3); putting the difference into the gas row conserves W
  and fails in the cold infall (§14.5), which is also why "total rows everywhere" failed (t_mg25, §6.3).
* The radiation energy row must be the conservation law of the *lab-frame* energy wherever a closure (FLUXTIE) replaces the moment
  equation of the flux: the comoving row's work term is the adjoint of the moment equation, not of the closure. This removes the
  3280 L☉ created in zones 88–100 (§3). The frame-consistent closure (H2) is correct and needed but does not remove it by itself.
* The kinetic energy the upwinded momentum flux dissipates (190 L☉ time-weighted, 3 % of L_em, in the AV-smeared accretion shock,
  Mach 3–11) is returned to the gas as heat. It has to be returned in the right place: put into the zone where it is computed it is a
  heat dipole across the photosphere and the run crawls (t_en08: 0.20 yr in 500 steps); transported with the sub-zone kinetic energy
  (§4.4) it is sign-definite and of second order in smooth flow, heats nothing in the cold infall (|ΔT/T| < 1e-3 from zone 31 outwards,
  §5 gate d) and costs nothing: 500 steps cover 12.65 yr with 74 Newton failures, against 11.9 yr / 89 when the same energy is merely
  booked (t_en20) and 8.04 yr / 88 in the current stack (§6). Thermalise, not book.
* What is **not** shown: 1 L☉ at *every* stage. The identity closes per zone to 0.000 L☉; what is left (|Σ| > 1 L☉ in a third of the
  stages, 4e-6 of L_em time-integrated) is the residual the Newton iteration leaves in the solved rows, a stopping-criterion matter
  (§5a). A residual-based stopping test (ERESTOL, §17) removes five sixths of those stages at 8 % of the solver's wall time; 48 of
  1000 stay above 1 L☉ (largest 11 L☉). And the new rows give a *different star* than the current stack (cooler, larger; from the
  IC a central temperature lower by 27–29 % at equal mass, which the control on the same tree attributes to the rows, §16): by a
  controlled run two thirds of the difference is the gravity term G5 of the star's zones (§14.5), and at 100 zones neither
  star is the converged one (§14.7). *Round 3 (§21):* with the zone number the two stars approach each other, towards the
  stack's; last bullet of this section.
* **For the decision on the default (§14.7):** counted against an explicit binding energy, ECOMPAT gives the gas of the star
  (zones 1–26) 2.0e3 L☉, about 30 % of L_em, that the binding energy does not lose; the current stack's total row is open there
  by 0.15e3 L☉. In the infall both are open by 0.6e3–0.75e3 L☉ (the momentum row's error, in every variant); in the outer cloud
  the stack creates 3.3e3 L☉ and ECOMPAT 3 L☉. Whole-mesh ledger against the binding energy: stack 3.8e3–4.1e3, ECOMPAT
  2.6e3–2.9e3, ECOMPAT with the gravity term put back in the star's zones (ECGRAV=1) 0.7e3–1.0e3 L☉ — but that last one needs
  the dual-energy flag again, and without the flag it fails (t_en32). So ECOMPAT as it stands is not the more conserving
  formulation for the star, and the measurements speak against making it the default unchanged. A mesh weight on density jump ×
  gravity power does not reduce the term at 100 zones (§14.6). *Round 3:* the stack's 3.3e3 L☉ of the outer cloud go with
  the lab-frame radiation row alone (`ERLAB`, §18); first bullet of round 3 below.
* A radiation energy density below the bath's (first round: open point 10) is the flux-tie closure's defect, not ECOMPAT's: it is
  in the control and, deeper (E down to 0.02–0.3 % of a T_CMB⁴), in the production run (§15). Not fixed; §15.5 says what a fix needs.
* The isolated-star line (S5, 100 M☉, 100 zones) runs with the same rows at the same Δt, Newton cost and ledger closure as with its
  own combined row, over the first 170 yr of the contraction and only with the exchange cap EXCAP=1e8 that the accretion stack
  already uses: a gas row that carries the gas–radiation exchange explicitly is singular to round-off in an opaque interior
  (§7.3). Not run to the ZAMS.

**Round 3 (§18–§21), the four measurements asked for; where they change a bullet above, the bullet says so.**

* **The lab-frame radiation row alone (`ERLAB`, §18, §20.4).** With the stack's gas rows unchanged it gives the stack's star
  with the outer cloud's defect removed. From the checkpoint: r_ph within 0.4 %, T_eff within 8 K, T_c within 0.05 % of the
  stack over 7.4 yr, the residual of zones 51–100 from 3296 to 1.9 L☉, the same steps and failures; binding-energy ledger
  527 / 479 / 741 L☉ (W i / ii / iv) against 3855 / 3805 / 4062 (stack), 2614 / 2602 / 2886 (ECOMPAT) and 700 / 665 / 957
  (ECOMPAT + ECGRAV=1). From the IC: the control's star at equal mass (R\* within 1 %, T_c within 0.3 %), cloud residual 4.8
  against 3407 L☉, ledger −375 / +313 L☉ against +3217 / +3974 L☉ (M\* 0.10–0.20). At 200 zones the run with the comoving row
  does not pass the collapse and the run with the lab-frame row does (§21.1). **This is the minimal change to the production
  stack and the configuration the measurements favour.**
* **The gravity weight of the momentum row (`GFV`, `GVAR`, §19).** The finite-volume weight does in the infall what §14.4
  estimated (G5 there 885–917 → −18…−36 L☉, d(u²/2)/d(GM/r) 1.09–1.12 → 1.00 ± 0.03) and from the checkpoint closes the ledger
  of the stack + `ERLAB` to −137 / −202 / +105 L☉. It has no binding energy (on the Lagrangian S5 star the ledger opens from
  3.5e-7 to 2.6e-3 of the energy radiated), leaves the density-jump part of G5 (21–30 % of the release rate at 100 zones), and
  at 100 zones from the IC it makes a star whose centre is 0.33–0.35 dex less dense and 19–21 % cooler, further from the
  200-zone star than the node weight's. At 200 zones that error is gone (centre within 0.02–0.03 dex and 2 % of the
  node-weight star, G5 0.075 of the release rate against 0.070, row residual −23 against 198 L☉, §21): it is a weight for a
  mesh that resolves the centre, not for the 100-zone mesh. A form with both limits, continuous in the mesh velocity and in u and without a
  flag, does not exist (§19.2). The variational weight keeps a binding energy, overshoots in the infall and changes the centre
  most (from the IC at 100 zones: centre 0.43–0.46 dex less dense and 36–39 % cooler than the node weight's; ledger against
  its own energy W (iv) −187 L☉, against W (i) −952 L☉; no 200-zone run). Neither is proposed as a default at 100 zones.
* **ECGRAV=1 from the IC (§20).** It restores the control's central temperature (+1.4…+3.0 % at equal mass; ECOMPAT alone
  −27…−28 %) and most of its central density (−0.05…−0.064 dex; ECOMPAT alone −0.26…−0.28 dex), and the radius to +3…+4 %
  (ECOMPAT alone +11…+18 %); at 50 zones likewise (log n_c 22.67 against the stack's 22.71 and ECOMPAT's 22.18). At 200 zones (to M\* 0.15) it is
  the stack's star to 0.008 dex in n_c, 1.6 % in T_c and 0.5 % in R\* (ECOMPAT alone −0.11 dex, −12 %, +5 %; §21.6). So
  the difference between the two stars is the G5 term of the star's zones, measured. Over M\* 0.10–0.335 the ledger of
  ECOMPAT + ECGRAV=1 is −424 / +119 L☉ (W i / iv), ECOMPAT's +1610 / +2023 L☉. Radius at M\* 0.33: control 44.1, ECOMPAT 48.9,
  ECOMPAT + ECGRAV=1 46.0, z100mm 55.6 R☉; z100mm's is larger because its early settings (before SLABFLUX / SLABTIE) kept its
  photosphere dark (1.07e44 erg radiated since core formation against 6.45e44 erg), not because of the energy rows (§20.3).
* **Zone-number convergence, 50 / 100 / 200 zones from the IC (§21; convergence tests, production stays at 100).** G5 is
  63–86 % of the release rate G M\* Ṁ / R\* at 50 zones, 23–27 % at 100 and 6.6–7.0 % at 200: order 1.4–1.8 from 50 to 100 and
  1.8–1.9 from 100 to 200, with every set of rows. Its density-jump part in the star's outer layers converges at order 1.9
  (100 → 200), not first order as §14.6 said; its geometric part in the infall at order 1.4–1.6. **The two stars approach
  each other, towards the stack's**: ECOMPAT − stack is −0.53 / −0.26 / −0.11 dex in n_c and −41 / −27 / −12 % in T_c at 50 /
  100 / 200 zones; between 100 and 200 zones the stack's star moves by +0.05 dex and +1…4 %, ECOMPAT's by +0.22 dex and
  +22…25 %. This replaces "at 100 zones neither star is the converged one" above: the stack's 100-zone star is within 0.05
  dex, 1.2 % (T_c) and 1.3 % (R\*) of the 200-zone one, ECOMPAT's 0.33 dex, 28 % and 16 % away. Not converged at 200 zones in
  either family: the photosphere (r_ph 49 / 40 / 35 R☉), the accretion history, and the difference between dW/dt and the
  total row's gravity term (−0.8e3 L☉, 5 % of the release rate, at 100 and at 200 zones). It is located (§21.3): the
  mesh-motion part of G5 inside 0.5 R\*, where GRIDEQ moves the mesh through gas nearly at rest and the mesh velocity
  doubles with the zone number. A static test decides what it is (§21.3): not a drift of the discrete W (the W of a
  fixed smooth star changes by at most 18 L☉ under the runs' mesh motion) but the interface density of the transport. The van
  Leer value, interpolated in the zone index, is 2–9 % (200 zones) and 7–29 % (100 zones) below the density of the smooth
  mass distribution at the interfaces of the star's stretched central mesh, while the node mass of the momentum row agrees
  with it to 0.6 % and 2.5 %; the mesh moving inward therefore carries mass with it, or the gas pays for staying. Evaluated
  on the runs' states with the consistent density in place of the advected one, G5 as a whole is 0.7 % (200 zones) and 2.5 %
  (100 zones) of the release rate instead of 7 % and 26 %, and the term −62 and −88 L☉ instead of +798 and +939. That is an
  evaluation, not a run: a geometry-consistent interface density is a change of the advection of every row and was not
  implemented. As the rows stand the term limits the ledger's total to 2–5 % of the release rate for both formulations at 200 zones (stack + `ERLAB`
  −4.1 / −3.1 %, ECOMPAT +2.4 / +3.4 %), so at 200 zones the total no longer discriminates between them, while its parts
  by region do.
* **A set of energy-formulation principles (conservative total-energy budget with directly evolved components and
  compatible transfers) is assessed against these rows in §22.** They agree with the measurements and do not change the
  recommendation; the one that is not met in any variant is "the same mass fluxes in momentum, gravitational work and
  energy transport", which is G5. A residual reconciliation step is not recommended: for internal + kinetic + radiation
  energy nothing is left to reconcile, and against the binding energy the residual is not small.

## 1. Mesh, unknowns, rows

TITAN staggered mesh (Gehmeyr & Mihalas 1994, Physica D 77, 320), spherical (μ = 2): zone k between interfaces k and k+1 carries ρ, T, E
(radiation energy density), the species; interface j carries r, u, F, the exterior mass. m_k = ρ_k ΔV_k; ṁ_j = −r_j² ρᵇ_j u_rel,j is the
mass flux through interface j relative to the moving mesh (GRIDEQ), positive for inflow. Time integration: TR-BDF2 (Bank et al. 1985,
IEEE Trans. ED 32, 1992; Hosea & Shampine 1996, Appl. Numer. Math. 20, 21), γ = 2 − √2: a trapezoidal stage (θ = ½) over γΔt, then a BDF2
stage written as a backward-Euler solve (θ = 1) over (1 − γ)/(2 − γ) Δt from the reference state y* = A y_γ − B y_n,
A = 1/(γ(2 − γ)), B = (1 − γ)²/(γ(2 − γ)), A − B = 1, formed on conserved combinations (m, m e, E ΔV, the interface momenta).

Momentum row of interface j (julia/TitanResidFast.jl, "GAS MOMENTUM"), M_j = ½(m_j−1 + m_j):

    Δ(M_j u_j)/Δt − [Φ_j − Φ_j−1] + r_j^μ (p_j − p_j−1) + r_j^{−μ/2} [r̄_j^{3μ/2} q_j − r̄_j−1^{3μ/2} q_j−1] + g_j M_j − (χ̄ F/c)_j M_j = 0,
    Φ_k = ½(ṁ_k + ṁ_k+1) û_k,   û_k the van Leer value of u at the centre of zone k (upwind, slope 2ab/(a + b) for ab > 0, else 0).

Energy rows of the current stack ("dual energy"): the gas-alone first law (GASROW) where a zone is flagged `useint` (ETASW, ETASW_EOFF,
ETASWTH, ETASWTIE: 74 of 100 zones here), otherwise the conservative row of gas + kinetic + radiation energy (`etot_row`, ETOTAL with
ETOT_GRAV); the radiation energy row is the comoving-frame one everywhere.

## 2. Measurement: where the +3740 L☉ comes from (task 1)

`zone_diag` (EDIAG=2, commit 774c85b) evaluates the conservative total-energy row T_k of every zone at every converged stage and splits
it exactly (remainder `rem` at round-off) into the terms by which that row differs from the rows the solve used:

    T_k = G_k + E_k + ½(u_k M_k + u_k+1 M_k+1)        the rows of the solve (zero to Newton tolerance)
        + dΨ + dΨr                                    a different localisation of the energy flux (telescopes: no source)
        + K2a + K2b                                   kinetic energy: time discretisation of the momentum row
        + K3a                                         kinetic energy: upwinded momentum flux
        + prs + V4                                    pressure work (identically 0), viscous-work form
        + G5                                          gravity power: (mass flux) g Δr of ETOT_GRAV − ½ m_k (u_k g_k + u_k+1 g_k+1)
        + R6a + R6b + R6c                             radiation: ½ Σ_j u_j × (flux row in moment form), its sphericality form, its inertia
        + adv                                         advection limiters (gas + radiation advected together against separately)
        + old                                         old-level energy of the two rows (§2.3)

### 2.1 Table (t_en01, 60 steps, medians over the stages; BDF / TR where they differ) [L☉]

| zones  | T (imb)       | K3a    | G5            | dΨ   | V4  | adv         | R6a             | R6c        | old |
|--------|---------------|--------|---------------|------|-----|-------------|-----------------|------------|-----|
| 1–26   | 145 / 157     | +60    | +1193 / +1256 | −8   | 0.8 | −0.1 / −14  | 0               | 0          | +526 / −1271 (zones 14–20) |
| 27–50  | 473 / 496     | −224.5 | +688.4        | +8.9 | 0.5 | 0.1 / 22.5  | 0.1             | 0          | 0 |
| 51–87  | 18.6 / 8.8    | 0      | 5.3           | −0.5 | 0   | 0 / −9.8    | 13.6            | 0.2        | 0 |
| 88–100 | 3280          | 0      | 0.1           | 0    | 0   | 0           | 3273.5 / 3280.5 | 6.7 / −0.3 | 0 |

(The previous agent's numbers for the same stack, t_mg24/t_mg26: −0.4, +442, +19, +3280.)

### 2.2 What the pieces are

* **Zones 27–50, the "+442 L☉"**: G5 (+688) + K3a (−225) + dΨ (+9). G5 is not an error of either row but a *definition* mismatch: the
  internal-energy rows conserve energy with the gravity power of the momentum row, the diagnostic total row uses another discrete
  gravity power. K3a is a real loss: kinetic energy the upwinded momentum flux dissipates and the first-law row never sees (225 L☉).
* **The ± alternation** (zones 29/30/31: −166, −479, +628 L☉) is dΨ: the total row and the first-law-plus-momentum rows place the same
  energy flux at different faces. It telescopes; it is a relocalisation, not a source. With it removed no alternating residual is left.
* **Zones 88–100, +3280 L☉**: entirely R6a = ½ Σ u_j × (moment-form radiation momentum residual) at faces where FLUXTIE has replaced
  the moment equation by the closure. See §3.
* **Zones 51–87**: 13.6 L☉ of the same radiation term at partially tied faces, 5 L☉ of G5.
* **Star, zones 1–26**: total rows are used here, so T is the Newton remainder plus `old`; G5 = +1193…1256 L☉ is the amount by which
  the star's gas is driven away from the first law (§7).

### 2.3 The old-level energy (`old`)

In COUPLED zones `cpl_prepare!` re-evaluates e_g(ρ, T, x) between steps with the finished step's species and the new step's expansion
terms and writes it to the old arrays; the total row differences the solve's own value. The two differ by δe/e ≈ 1e-5…−9.4e-5 per step
in zones 14–20: +526 (BDF) / −1271 (TR) L☉ of operator-split energy jump. The compatible gas row (§4) carries the energy the previous
solve left, e_old = (e_g + E/ρ)_old − E_old/ρ_old.

## 3. Radiation: frame of the closure, frame of the energy row (task 2, H2)

### 3.1 The closure in the frame of the unknowns (TIEFRAME, commit a37f48b)

E and F of the rows are comoving-frame moments; the formal solution that provides h = H/J, f = K/J, g is static (observer at rest).
To O(u/c) (Mihalas & Mihalas 1984, *Foundations of Radiation Hydrodynamics*, §91; Mihalas & Auer 2001, JQSRT 71, 61; Krumholz, Klein &
McKee 2007, ApJ 667, 626, their eqs. 31–33; Lowrie, Morel & Hittinger 1999, ApJ 521, 432):

    F_lab = F₀ + u (E₀ + P₀),        E_lab = E₀ + 2 u F₀/c².

The closure F_lab = h c E_lab therefore reads, for the unknowns, F₀ (1 − 2hu/c) = [h c − (1 + f) u] E₀, i.e.

    F₀ = [h c − (1 + f − 2h²) u] E₀ + O(u²/c).

Limits computed from this form: h = f = 1 gives F₀ = c E₀ (a beam is a beam in both frames to this order); h = 0, f = 1/3 gives
F₀ = −(4/3) u E₀, the flux an observer moving with u sees in an isotropic field (the same number follows directly from the first
relation with F_lab = 0, P₀ = E₀/3). Outer boundary: F_lab = g c E_s − π(2g + 1) I_in becomes
F₀ = [g c E_s − π(2g + 1) I_in](1 + 2gu/c) − (1 + f) u E.

Discrete form at interface k, everything at the θ-level of the row and at the Newton iterate (nothing lagged), built so that a uniform
isotropic field is a steady state of the energy row to round-off also on the moving mesh:

    r² F₀ = r_c² c (hE)_donor + ṁ_k (E/ρ)ᵇ_k − V̇_k ½(E₋ + E₊) − r² u ½(f₋E₋ + f₊E₊) + 2 u r_c² (h²E)_donor,

ṁ_k (E/ρ)ᵇ_k the energy row's own advective flux, V̇_k = (r_new³ − r_old³)/(3Δt) the volume swept by the interface (geometric
conservation). The same shift is in SLABTIE's donor form and in the outer boundary row; 50f1609 puts it into the per-group tie and the
group boundary row of julia/MGRad.jl (§9).

### 3.2 Results

Off-path: TIEFRAME off is bit-identical (t_en07 against t_en01, every state field of the per-stage records).

t_en06, TIEFRAME alone, 500 steps = 7.26 yr:

| quantity | start | end |
|---|---|---|
| Σ T, zones 88–100 | 3280 L☉ | 3280 L☉ |
| L_F at the outer boundary (interface 101) | −1.001e5 L☉ | −6.45e4 L☉ |
| E_100/E_ref − 1 | −2.2e-6 | +3.8e-8 |
| E/E_ref − 1 at zones 88 / 92 / 96 / 99 | 0.626 / 0.278 / 0.134 / 0.034 | unchanged |
| E/E_ref − 1 at zone 80 | 1.767 | 1.841 |
| L_F at interface 100 | 6.9e3 L☉ | 1.28e5 L☉ |
| star side | photosphere zone 27, L_em(1.5 r_ph) median 7254 L☉, T_eff ≈ 4990 K, 87 Newton failures / 500 steps | |

**Negative result for H2 as stated**: the frame-consistent closure does not remove the 3280 L☉. The reason is in the energy row, not in
the closure: the comoving energy row's work term f E Δ(r²u) + (1 − 3f) E (u/r) ΔV is the discrete adjoint of the *moment* equation's
∇·P. Total (lab-frame) radiation energy is conserved by the pair {comoving energy row, moment flux row}. Where FLUXTIE replaces the
flux row by a closure, u × (∇·P + χρF/c) is no longer balanced by anything and is created: R6a.

### 3.3 The lab-frame energy row (part of ECOMPAT, commit be38295)

The radiation energy row is written as the conservation law of the lab-frame energy in the comoving unknowns, to O(u/c):

    Δ(E ΔV)/Δt + Δ(r² F₀) − Δ(ṁ (E/ρ)ᵇ) + Δ(r² u (fE)‾) + Σ_faces u_j (χ̄ F/c)_j ½ m_k = sources.

It equals the comoving row plus u × (moment equation) exactly where the moment equation is the flux row, so nothing changes there;
where a tie is used it is the statement that E_lab changes only by F_lab = F₀ + u(E₀ + P₀) through the faces and by the exchange with the
gas (the u χ F/c term is the work the radiation force does on the gas; the momentum row has its counterpart, so the pair telescopes).

Gate "zones 88–100 imbalance down ≥ 100× with nothing booked" — **met**: Σ T over zones 88–100, median / 90 % / max over 1000 stages:
t_en08 −0.000 / 0.91 / 4.8 L☉; final form t_en19 0.050 / 0.60 / 5.9 L☉ (was 3280). The residual that is left is round-off of the row itself: the
radiation energy stored in these zones is 9.8e52 erg, and one unit of round-off (2.2e-16) of it over a stage Δt of 1–2e3 s is 2.5–4.8 L☉.

Gate "E relaxes to T_CMB⁴, boundary artefact disappears" — **cannot be shown from the checkpoint**, and is not claimed. The checkpoint
carries the excess the old rows built (E/E_ref − 1 = 0.63 at zone 88, 1.77 at zone 80). Its relaxation time is r/(h c) with
h ~ 1e-6…1e-3, and r/c alone is 74 yr at the outer boundary. Transient in t_en19 (12.6 yr): the excess at zones 88 / 92 / 96 / 99
stays 0.626 / 0.278 / 0.134 / 0.034 to three digits, zone 80 goes 1.77 → 1.91, L_F(101) −1.009e5 → −4.44e4 L☉, L_F(100)
1.265e5 → 1.303e5 L☉ (the value of an isotropic lab-frame field seen by the infalling interface, −(1 + f) u E 4πr², is 1.180e5 L☉
there; u = −1.14 km/s), zone 100 itself is at E_ref to 2e-6. The test of the gate proper is the run from the initial condition,
where the excess must never build. It is **met there** (§8): E/(a T_CMB⁴) − 1 = 2.3e-6…5.4e-6 in zones 88–100 (u/c = 3.8e-6),
L_F(101) = +3e-3 L☉ before the star exists and +60…100 L☉ after, against 0.03…0.63 and −1.16e5 L☉ in z100mm at the same stage.

Star side unchanged within noise (gate): §5, gate c.

## 4. The compatible gas row (task 3, H1)

Method: Caramana, Burton, Shashkov & Whalen 1998 (JCP 146, 227): on a staggered mesh, total energy is conserved to round-off by an
internal-energy equation if the work terms of that equation are the discrete adjoints of the force terms of the momentum equation;
the kinetic-energy fix attributed to DeBar as described by Benson 1992 (CMAME 99, 235): kinetic energy lost by the momentum
advection of a staggered remap is computed explicitly and returned to the internal energy. Bryan et al. 1995 (CPC 89, 149) is the
dual-energy alternative the current stack follows. The derivation below is for this mesh: staggered, μ = 2, moving mesh with relative
mass fluxes, θ-weighted stages.

Multiply the momentum row of interface j by the θ-level velocity u_j and distribute each product to the two adjacent zones
(half the inertia and weight to each; the zone's own p, q and Φ). Per zone k, with faces j = k, k+1:

### 4.1 Time term (K2)

With m̂ = ½(m_n + m_o), δ_j = u_j^n − u_j^o and u = θ u^n + (1 − θ) u^o,

    Σ_faces u_j Δ(½ m_k u_j)/Δt = Δ(¼ m_k Σ u_j²)/Δt + (2θ − 1) ¼ m̂ Σ δ_j²/Δt − (θ² − θ + ½) ¼ (m_n − m_o) Σ δ_j²/Δt + (mass-flux terms that
    join the advective kinetic-energy flux).

The first term is the zone's kinetic energy K_k = ¼ m_k (u_k² + u_k+1²). The second is the dissipation of the implicit time
discretisation (zero for θ = ½, ¼ m̂ Σδ²/Δt for θ = 1, never negative); the third is its mass-change companion. Both are products of
differences; the row gets −(2θ − 1) ¼ m̂ Σδ²/Δt + (θ² − θ + ½) ¼ (m_n − m_o) Σδ²/Δt.

### 4.2 Pressure and viscous work (prs, V4)

u_j r_j^μ (p_k − p_k−1) summed by parts gives p_k (r_k+1^μ u_k+1 − r_k^μ u_k): the first law's p dV as it already was (prs ≡ 0). The
viscous force r_j^{−μ/2} Δ(r̄^{3μ/2} q) gives q_k r̄_k^{3μ/2} (u_k+1 r_k+1^{−μ/2} − u_k r_k^{−μ/2}); the stack's row had q (du/dr) ΔV. The
difference V4 is 0.5–1.3 L☉.

### 4.3 Gravity (G5)

No term. A row that conserves internal + kinetic energy against the momentum row's power Σ_j u_j g_j M_j needs nothing else. The
isolated-star line's ledger (julia/s5_run.jl `energies`) already uses this form.

*Corrected by §14.* The first round wrote here that Σ_j u_j g_j M_j "*is* the gravity power of this discretisation". It is the power
of the momentum row's force, and exact only where the mesh carries the mass (Lagrangian limit). Where mass crosses the
interfaces, the rate of an explicit binding energy is (mass flux) × (potential drop), which is the total row's form; the two
differ by G5 (§14.2–§14.4).

### 4.4 Momentum flux (K3a) and the transport of the sub-zone kinetic energy

The zone's part of Σ u_j (Φ_j − Φ_j−1) against the kinetic-energy flux ṁ_j ½u_j² through its faces leaves

    K3a_k = ½(ṁ_k + ṁ_k+1) (u_k+1 − u_k) (½(u_k + u_k+1) − û_k),

zero for the centred value û = ½(u_k + u_k+1), a loss for the upwind-biased one. Returning −K3a_k to zone k (the DeBar-type fix, first
implementation, be38295) conserves, but K3a is third order in the velocity differences and **not sign-definite**: with a = u_k+1 − u_k
and b the same difference of the upwind neighbour, van Leer's slope gives K3a = −½ M̄ a² (a − b)/(a + b), M̄ = ½(ṁ_k + ṁ_k+1). Measured
(t_en08, medians): heating 38 / 74 / 83 L☉ in zones 27 / 28 / 29 and *cooling* 13 / 34 / 11 L☉ in zones 24 / 25 / 26 — a heat dipole
across the shock.

The missing piece is that the staggered kinetic energy is more than ½ m ū²:

    K_k = ¼ m_k (u_k² + u_k+1²) = ½ m_k ū_k² + m_k s_k,     s_k = ⅛ (u_k+1 − u_k)²,

and the sub-zone part s is carried by the mass flux like any specific quantity. With its upwind flux G_j = ṁ_j s_up(j) (zone j for
ṁ_j ≥ 0, zone j − 1 otherwise; 0 at the boundaries) the heat returned to zone k is

    H_k = −K3a_k + (G_k+1 − G_k).

The added term telescopes, so Σ H = −Σ K3a: the same total, another place. For steady inflow (ṁ_k = ṁ_k+1 = M̄ > 0, θ = 1), monotone data,
slope 2ab/(a + b):

    H_k = ½ M̄ [½ (a − b)² + ¼ (a − b)³/(a + b)] ≥ ⅛ M̄ (a − b)² ≥ 0        (since |a − b| ≤ |a + b| for ab ≥ 0),

and at an extremum (slope 0) H_k = ⅜ M̄ a² + ⅛ M̄ b² ≥ 0. So the dissipation is non-negative, of order (Δ²u)² in smooth flow (nothing to
thermalise in smooth cold infall, whatever K/e is — no criterion "where the AV is active" is needed), ⅜–½ M̄ (Δu)² where the velocity
jumps, and it is deposited where the mass that carried the sub-zone energy arrives. Outflow is the mirror image. Not covered by the
bound: unsteady mass (the remainder ½(ṁ_k − ṁ_k+1)(s_k + s_k+1)) and the θ-level mixing of the trapezoidal stage.

Measured with the data of t_en08 (BDF stages, medians [L☉]; the two forms have the same sum 167.4):

| zone | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31 | 32 | Σ of negative parts, all zones |
|---|---|---|---|---|---|---|---|---|---|---|---|
| in place, −K3a | −0.8 | −13.0 | −34.5 | −11.0 | 38.1 | 73.9 | 82.7 | 3.7 | 10.1 | −7.7 | −68.8 |
| transported, H | 6.4 | 26.9 | 16.7 | 0.6 | 5.8 | 25.6 | 54.1 | 24.1 | 2.0 | 2.8 | −0.41 |

(most negative single value of H in the run: −15 L☉ in zone 26, where ṁ changes across the zone.) Commit e0e6b2e.

### 4.5 The kinetic energy of the TR-BDF2 reference state (KSTAR)

A stage row that differences against y* conserves E(y*). The step conserves if E(y*) = A E_γ − B E_n. That holds for everything y* is
built on linearly (m e, E ΔV) but not for the kinetic energy: y* combines m and m u, and for equal masses

    D_k = K(y*) − (A K_γ − B K_n) = A B ⅛ m_k Σ_faces (u_γ − u_n)² > 0        (A B = ¼).

Measured (t_en08): D/Δt = +0.077 L☉ time-weighted (zones 1–26 0.074, zones 27–50 0.003; per step median 0.05, max 0.48 L☉), and the
BDF stage's K2 terms return exactly that much as heat (−0.0771 L☉): every step created 1e-5 L and thermalised it. The same defect
is in the total rows of the current stack. `Trbdf2.bdf_old!` now hands D_k to the gas row, which takes it out of its old-level energy
(+D_k/Δt in the row); `restore!` empties it. Final form (t_en19): D/Δt = 0.134 L☉ time-weighted (per step median
0.050, max 0.83 L☉) is in the row and in the ledger (`kst`), and the per-zone identity closes with it.

### 4.6 The row

    Δ(m e)/Δt + (advection of e_g) + (conduction, convection) + p Δ(r^μ u) + q r̄^{3μ/2} Δ(u r^{−μ/2})
      − (2θ − 1) ¼ m̂ Σδ²/Δt + (θ² − θ + ½) ¼ (m_n − m_o) Σδ²/Δt + D_k/Δt + K3a_k − (G_k+1 − G_k)
      − (gas–radiation exchange) + Compton − sources = 0,

with e_old as in §2.3, in every zone; no `useint`, no η, no E − K. One provisional switch `ECOMPAT=1` (popiii_m2c; needs GASROW, turns
TIEFRAME on); default off.

The identity that proves it, evaluated by `zone_diag` at every converged stage (field `cres`):

    T_k − G5_k − adv_k − [rows of the solve]_k − dΨ-type relocalisations − (G_k+1 − G_k) + D_k/Δt = 0   per zone, to round-off.

## 5. Gates (task 3)

Run t_en19: tree e0e6b2e, `ECOMPAT=1 EDIAG=2`, 500 steps from the checkpoint = 12.65 yr, 1000 converged stages. Comparison runs over the
same 500 steps: z100mm itself (the current stack, 8.04 yr) and t_en06 (current rows + TIEFRAME, 7.26 yr).

**(a) Energy closure.** Sums over zones of the total-energy imbalance with nothing booked except the two definitions
(G5: which gravity power; adv: which advection limiter) [L☉], median / 90 % of |.| / max |.| over the 1000 stages:

| zones | T (imb) | T − G5 − adv + D/Δt | per-zone identity `cres` | Newton remainder of the solved rows |
|---|---|---|---|---|
| 1–26   | 1920 / 2493 / 3111 | 105.7 / 116 / 160 (flux of sub-zone energy through interface 27) | 0.000 / 0.000 / 0.000 | 0.000 / 9.6 / 103 |
| 27–50  | 595 / 792 / 877    | −104.4 / 114 / 128 (the same flux)                              | 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 |
| 51–87  | 4.7 / 4.9 / 6.4    | −0.53 / 0.62 / 0.74                                             | 0.000 / 0.000 / 0.002 | 0.000 / 0.001 / 0.005 |
| 88–100 | 0.050 / 0.60 / 5.9 | −0.009 / 0.61 / 5.8                                             | 0.000 / 0.52 / 3.96   | 0.000 / 0.38 / 4.7 |
| all    | 2655 / 3089 / 3697 | **0.000 / 9.4 / 107**                                           | 0.000 / 0.52 / 3.96   | 0.000 / 9.6 / 104 |

Time-integrated over the run (Σ row × Δt / Σ Δt): all zones +0.0009 L☉ (trapezoidal stages), +0.026 L☉ (BDF stages); mean of the
absolute value 0.011 and 0.18 L☉; i.e. 4e-6 of L_em. Ledger of the same run, time-weighted [L☉]: G5 1865…1878 (zones 1–26),
839…843 (27–50), 5.5 (51–87); upwind kinetic-energy loss returned as heat 62 (1–26) + 127 (27–50) = 189; V4 1.5; D/Δt 0.46 (BDF stages).

*The gate as written (below 1 L☉ at every stage) is not met*: |Σ| > 1 L☉ in 316 of the 1000 stages. All of it is the last column: the
residual the Newton iteration leaves in the gas and radiation rows of the exchange-stiff zones of the star (per zone the two
remainders are equal and opposite to about three digits, up to ±4e4 L☉; their sum is what is left). The identity itself —
everything that is not a row of the solve — closes per zone to 0.000 L☉ in zones 1–87 and to the round-off of the radiation row
in zones 88–100 (9.8e52 erg of stored radiation energy differenced over the stage Δt: 2.5–4.8 L☉ per unit of round-off). The solved rows of the current stack leave the same remainder
(t_en01: median −0.02, 90 % 14, max 128 L☉). Obstruction: the Newton stopping criterion (largest relative change of the unknowns
< 1e-5; `CONV=1e-3` in run.env is not read by popiii_m2c), not the formulation. Meeting 1 L☉ per stage needs a residual-based
stopping test on the two energy rows; that is a solver change and was not made here.

*No alternating per-zone residual*: `cres` is 0.000 (median and max) in every zone of 1–87. The ±100 L☉ in the range sums of column 3
is one flux, G_27 (sub-zone kinetic energy carried through interface 27), entering one range and leaving the other.

**(b) No row switching, no E − K.** By construction: the ECOMPAT branch of the gas row does not read `useint`, `ETASW*` or
`etot_row`; kinetic energy enters only through products of differences (§4). (`zone_diag` still evaluates `etot_row` as a diagnostic.)

**(c) Star side, Newton, Δt** (same 500 steps; photosphere = τ 2/3 from outside, L_em at 1.5 r_ph):

| | z100mm (current stack) | t_en06 (current rows + TIEFRAME) | t_en19 (ECOMPAT) |
|---|---|---|---|
| time covered by 500 steps | 8.04 yr | 7.26 yr | 12.65 yr |
| Newton failures per 100 steps | 11, 15, 22, 20, 20 (88) | 12, 14, 19, 20, 21 (86) | 17, 15, 14, 9, 19 (74) |
| accepted Δt, every 20th step: median / mean | 1.0e5 / 3.4e5 s | 6.8e4 / 3.3e5 s | 3.2e4 / 9.7e5 s |
| photosphere zone | 26 → 27 (between 0.23 and 4.1 yr) | 26 → 27 (before 3.0 yr) | 26 → 27 (between 0.50 and 3.9 yr) |
| r_ph [R☉] | 111.1 → 105.9 (4.1 yr) → 109.0 (8.0 yr) | 110.5 → 105.1 → 108.3 | 111.1 → 113.3 (3.9 yr) → 121.4 (10.2 yr) → 121.7 (12.6 yr) |
| L_em [L☉] | 6027 → 6632 (4.1 yr) → 7526 (8.0 yr) | 6103 → 6242 → 7381 | 6139 → 6731 (3.9 yr) → 8296 (10.2 yr) → 9145 (12.6 yr) |
| T_eff [K] | 4334 → 5012 → 4979 | 4349 → 5000 → 4987 | 4383 → 4832 → 4833 |
| M* [M☉] (perf.csv) | 0.877 → 1.186 (8.0 yr) | 0.877 → 1.161 (7.3 yr) | 0.879 → 1.031 (3.9 yr) → 1.355 (12.6 yr) |

While the photosphere is in zone 26 (first 0.5 yr of t_en19, 7 snapshots): L_em 6070–6204 L☉ (rms 0.8 %), T_eff 4338–4390 K,
r_ph 110.8–112.2 R☉. Every run makes the one zone step 26 → 27 (the photosphere is quantised by the 100 zones), and L_em rises with
M* afterwards at about the same rate (+13 % in 3.9 yr in z100mm, +36 % in 8.7 yr in t_en19). Failure rate: not worse (74 against
88). Δt: the 500 steps cover 1.6 times the time of the current stack (mean Δt 8.0e5 against 5.1e5 s), but not uniformly — the
first 300 steps of t_en19 cover only 0.22 yr at Δt 7e3–1e5 s (the current stack reaches Δt ~1e6 s after 60 steps), then Δt opens to
1e5–1.2e7 s, and the last 40 steps are at 3–6e3 s again; z100mm ends its 500 steps with 200 steps at 3–5e3 s. So the sampled
median Δt is 3 times smaller and the mean larger. "Δt not smaller" is therefore met in the mean and in the time covered, **not
in the median**. Every stack has such slow phases, at different places of the 500-step window (the booked variant t_en20 has
its fast phase first and then 300 slow steps); §6 has the anatomy and §8 the comparison without a restart transient (from the IC).

Not unchanged within noise, and reported as a difference: after the zone step the photosphere of t_en19 sits at 113–122 R☉ and
4830 K, that of the current stack at 106–109 R☉ and 4980–5010 K. This is the star's energy budget (§7), not the radiation frame:
t_en06 (TIEFRAME alone) reproduces z100mm.

**(d) Cold infall not heated.** T, ρ, u of t_en19 against the current stack (t_en01) at equal time, 3.0e6 s after the restart:

| zone | 26 | 27 | 28 | 29 | 30 | 31 | 32 | 35 | 38 | 41 | 45 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| u [km/s] | −11.2 | −22.1 | −33.7 | −43.6 | −49.8 | −51.1 | −46.0 | −31.1 | −16.9 | −9.3 | −5.3 |
| T [K] (current stack) | 5331 | 4627 | 3785 | 2068 | 1965 | 1940 | 1915 | 1775 | 1640 | 1469 | 1231 |
| ΔT/T | −1.1 % | +0.03 % | +0.2 % | +1.7 % | +1.1 % | +0.09 % | +0.08 % | +0.03 % | +0.04 % | 0.00 % | +0.02 % |

The heat returned by the transported form (medians, BDF stages): zones 24–30: 26, 18, 0.8, 5.0, 24, 52, 25 L☉ — the zones in which the
flow decelerates from 51 to 3 km/s; zones 31–50 together 9 L☉; zones 51–87 0.009 L☉; zones 88–100 0.03 L☉. The smooth supersonic
infall (K/e up to 1e4) receives nothing measurable: |ΔT/T| < 1e-3 from zone 31 outwards, although no criterion selects where to
thermalise. The largest effect, +1.7 % in zone 29 (Mach 11, first zone of the deceleration), is where the current stack's rows lose
84 L☉ of kinetic energy without booking it.

## 6. Newton iteration and time step

### 6.1 Which part costs Newton iterations: 30-step runs from the checkpoint

Debug toggles of tree be38295_dbg / e0e6b2e_dbg (not in the branch) switch the parts of the new rows separately:

| run | rows | Newton failures in 30 steps |
|---|---|---|
| t_en10 | current stack | 3 |
| t_en11 | lab-frame radiation row + TIEFRAME, gas rows of the current stack | 3 |
| t_en12 | compatible gas row, radiation rows of the current stack | 6 |
| t_en13 | compatible gas row in zones ≥ 27 only | 7 |
| t_en14 | compatible gas row in zones ≤ 26 only | 2 (Δt reaches 1.41e5 s) |
| t_en15 | compatible gas row without the K3a return | 2 (Δt reaches 1.36e5 s) |
| t_en16 | compatible gas row with the old viscous form q (du/dr) ΔV | ≥ 6 |
| t_en09 | complete, K3a returned in place | 6 |
| t_en17 | complete, K3a returned with the sub-zone transport (final form) | 6 |

In 30 steps the cost is the K3a return in the zones above the photosphere, in either form; the failing unknowns are x_H+ and x_H2 of
zones 26/27 and E of zone 27. **These 30-step tests were misleading** about the final form, as the long runs show.

### 6.2 500 steps

| | z100mm (current stack) | t_en06 (+ TIEFRAME) | t_en08 (K3a in place) | t_en20 (K3a booked, not returned) | t_en19 (final form) |
|---|---|---|---|---|---|
| time covered | 8.04 yr | 7.26 yr | 0.20 yr | 11.90 yr | 12.65 yr |
| Newton failures | 88 | 86 | 79 | 89 | 74 |
| Δt sampled every 20th step, median / mean [s] | 1.0e5 / 3.4e5 | 6.8e4 / 3.3e5 | 1.1e4 / 1.3e4 | 5.4e3 / 3.9e5 | 3.2e4 / 9.7e5 |
| Σ over all zones of T − G5 − adv + D/Δt, median [L☉] | — | — | 0.000 | −197.8 (> 1 L☉ in 1000 of 1000 stages; time-weighted −190.1 TR / −190.8 BDF) | 0.000 |
| end: r_ph [R☉] / L_em [L☉] / T_eff [K] | 109.0 / 7526 / 4979 | 108.3 / 7381 / 4987 | (0.2 yr only) | 119.9 / 8967 / 4857 | 121.7 / 9145 / 4833 |

* Only the in-place return is intrinsically slow: all 500 steps of t_en08 are at Δt ~ 1e4 s. It cools zones 24–26 by 13 / 34 / 11 L☉
  and heats zones 27–29 by 38 / 74 / 83 L☉ (§4.4): a heat dipole across the zones whose ionisation state controls the step.
* **Thermalise or book?** Booking (t_en20: the upwind loss K3a − ΔG is computed and written to the ledger but not given to the gas)
  buys nothing: fewer years per 500 steps and more failures than the transported return, a ledger that is open by 190 L☉ at every
  stage, and the same star (r_ph 120 against 122 R☉, T_eff 4860 against 4830 K). The dissipation is real (it is the kinetic energy
  of the decelerating flow in the smeared shock); it is thermalised. The branch contains no booked variant.
* The step is limited by the same unknowns in all stacks: in t_en19 by x_H2 in 283 steps, x_H+ in 215, E in 1; in t_en20 by x_H+ in
  389, x_H2 in 104, r in 3, E in 3. All of them have slow phases at Δt 2e3–1e4 s lasting 40–300 steps while the photosphere crosses
  from zone 26 to zone 27 or the ionisation front moves through a zone: t_en19 in its first 300 and last 40 steps, t_en20 in its
  last 300, z100mm in its last 200. The median Δt of a 500-step window measures where these phases fall in the window; the time
  covered and the failure count are the robust numbers, and both favour the final form.
* Because the booked variant gives the star of the final form and not that of the current stack, the difference of the star side
  in §5(c) (r_ph, T_eff) is caused by G5 — the gravity power of the star's own zones — and not by the kinetic-energy return.

### 6.3 Why "total rows everywhere" failed (t_mg25)

docs/POPIII_MULTIGROUP.md §12: with ETASW = ETASWTH = ETASWTIE = 0 (the conservative row in all zones) the ledger closes to −1.1 L☉
"but the run fails after 23 steps (the photosphere jumps from zone 26 to 31; 66 Δt rejections against 8)". In the terms of §2:

1. the total row is the first law plus u × (momentum row) plus the terms of §2's split. With ETOT_GRAV its gravity power is
   (mass flux r² u ρᵇ) × g × Δr, so in each zone the gas is given G5 less heat than the momentum row's kinetic energy requires:
   69 / 78 / 116 L☉ in zones 32 / 33 / 34 and 93 L☉ in zone 36 (medians), where the whole advected thermal power is 65–75 L☉. The
   internal energy the row implies, E − K, is then wrong by O(1) in exactly the zones with K/e ~ 1e3: the PI's caveat, but the
   error is a systematic term of the row, not round-off of the subtraction (which at K/e = 1e4 costs 4 of 16 digits);
2. it returns K3a in place (the heat dipole of §4.4), and
3. being conservative for the sum, it charges to the gas whatever the other rows create: the `old` jump of §2.3 (±0.5–1.3e3 L☉ in
   zones 14–20) and the 3.3e3 L☉ of §3 in the tied zones (the effect ETASW_TH and ETASW_TIE were introduced against, see their
   docstrings: the cloud gas of the 146-zone run cooled to 11.7 K).

So a ledger that closes, with the wrong partition between e and K. This is a quantitative argument from the measured terms, not a
controlled test: t_mg25 was not repeated with one of the three removed. The compatible row has the momentum row's own gravity
power (no G5), transports the sub-zone kinetic energy, takes e_old from the previous solve, and never forms E − K; with it the same
"one row everywhere" runs 500 steps with fewer Newton failures than the dual-energy stack (§5c).

*Corrected by §14.* Item 1 has the sign of the judgement wrong: the total row's gravity power is the one an explicit binding energy
follows (§14.2), so the row does not give the gas "less heat than required" — the momentum row gives the infall too much kinetic
energy (§14.4), and a conservative row then takes the difference out of the internal energy. The mechanism of the failure stands and
was reproduced with one term alone: ECOMPAT plus the exact mass-flux gravity term in every zone (ECGRAV=3, t_en26) conserves the
binding-energy ledger and dies after 447 steps with zone 32 at 12.8 K (§14.5).

## 7. Star interior (task 4)

### 7.1 The budget of zones 1–26 over the restart

In the current stack the star's zones are total-row zones. What that row does to the star's gas, in the terms of §2:

| term, Σ over zones 1–26 [L☉] | current stack (t_en01, medians BDF / TR) | compatible rows (t_en19, 12.65 yr) |
|---|---|---|
| G5: gravity power of the total row − that of the momentum row | +1193 / +1256, heat withheld from the gas | in no row; the diagnostic total row evaluated on the solution gives +1865 time-weighted (t_en20: +1903; t_en08, first 0.2 yr: +1.0e3) |
| `old`: e_old re-evaluated between the steps (zones 14–20) | +526 / −1271 | none: e_old is the previous solve's |
| K3a: kinetic energy lost by the upwinded momentum flux | +60, charged in place | returned with the sub-zone kinetic energy (§4.4) |
| identity residual `cres` per zone | — | 0.000 |

G5 is 16–21 % of L_em. Per zone (t_en19, medians): 88 / 126 / 335 / 64 / 127 / 41 L☉ in zones 15 / 17 / 21 / 23 / 25 / 26; above the
photosphere 69 / 78 / 116 / 93 L☉ in zones 32 / 33 / 34 / 36 (0.69–0.85e3 L☉ in zones 27–50 together). It sits where the density is
steep across an interface: the total row takes (mass flux r² u ρᵇ) × g × ½Δr on each face with the advected interface density, the
momentum row u_j g_j ½(m_j−1 + m_j). Both are consistent discretisations of ρ u g and differ by a truncation term, which is not small
against L because the gravity power of the settling star and of the infall is many times L. It is not a mesh-motion term (with
GRIDEQ the total row uses the lab-frame u, TitanResidFast.jl:150, 1659).

Which of the two is "right" is not a question of order. The kinetic energy of the solution is whatever the momentum row makes it; a
gas row built on another gravity power puts the difference into the internal energy. The total row conserves a total energy whose
potential part is defined by its own flux form; the compatible row conserves internal + kinetic energy against the momentum row's
force. At 100 zones the two stars differ (§7.2). A zone-number convergence test would show how fast each approaches the other; under
the 100-zone rule it was not run.

*Corrected by §14.* (i) "Not a mesh-motion term": nine tenths of G5 is the mass-flux term G5x of §14.3, one tenth is what the total
row itself gets wrong against the binding energy in the mesh-motion part. (ii) "Both are consistent … which is right is not a
question of order": both are consistent, but they are not equally right. The binding energy follows the flux form (§14.2), the
infall is measurably faster than free fall with the momentum row's weight (§14.4), and G5 is made of two different error terms:
a first-order one in the star's outer layers (upwinded against node-mean density) and a second-order geometric one in the wide
zones of the infall (§14.4).

### 7.2 What it does to the star

* Restart, same checkpoint (§5c, §6.2): after the photosphere's zone step the compatible rows give r_ph 113–122 R☉ and T_eff 4830 K,
  the current stack 106–109 R☉ and 4980–5010 K; L_em 9145 L☉ after 12.6 yr against 7526 L☉ after 8.0 yr. While the photosphere is
  still in zone 26 (first 0.5 yr) the two agree: L_em 6070–6204 against 6027 L☉, T_eff 4338–4390 against 4334 K.
* From the IC (§8): at equal mass T_c = 0.84e5 K against 1.19e5 K and n_c = 10^22.93 against 10^23.20 cm⁻³; Σ T over zones 1–26 of the
  diagnostic total row is 2210 L☉ in the median (G5) while L_em is 0.5–2.5e3 L☉. A star whose gas is given that heat has more
  entropy: lower central density and temperature at equal mass. This is the direction observed.
* Attribution: by exclusion, not by a controlled run. The kinetic-energy return is excluded (t_en20, booked instead of returned,
  gives the star of t_en19: 119.9 against 121.7 R☉, 4857 against 4833 K); the radiation frame is excluded (t_en06 reproduces z100mm).
  G5 and `old` are the two remaining terms of size in the star; G5 is the steady one-signed one (`old` changes sign between the two
  stages of a step). A run of the compatible rows with G5 put back as a sink would settle it; it was not made.
  *Made in the follow-up (t_en25, §14.5):* G5 put back in the star's zones removes 55–72 % of the radius difference and 55–76 % of
  the T_eff difference at equal time; 2–4 R☉ and 35–60 K remain unattributed. The comparison of this bullet list is at equal
  step count, i.e. at different times and masses (8.0 against 12.6 yr); §14.5 compares at equal time.
* "A star whose gas is given that heat has more entropy" (second bullet): the direction is as observed, and the control from the
  IC confirms that the T_c and n_c differences are the rows' (§16); but by §14.2 that heat is not matched by a decrease of the binding
  energy (ECOMPAT's ledger against W is open by it), so the higher entropy is not evidence that ECOMPAT's star is the right one.

### 7.3 The isolated-star line (S5)

julia/s5_run.jl, 100 M☉ n = 3 polytrope, R₀ = 6.5 R☉, 100 zones (NENV=40), relaxation (52 steps), then the pre-main-sequence
contraction with the energy ledger of that line (closure and hydro-step error, both in units of Σ|E_bur|; L_nuc/L ≈ 1/60 in this
phase, so 3e-5 is 5e-7 of ∫L dt). Tree e0e6b2e_s5 = e0e6b2e plus three environment toggles in the driver. The zero-age main sequence
(7.9e3 yr in the line's 140-zone logs) was **not** reached: at the Δt of 2–4 yr these settings take in the contraction that is some
thousand steps, hours at the measured rate (160 steps, relaxation and compilation included, in 11 min). What was compared is the first 90–110 evolution steps.

| run | gas energy row | integrator | step: t [yr] / Δt [yr] / Newton it. | closure / hydro error | step failures |
|---|---|---|---|---|---|
| en_s5_ref = en_s5c_ref | combined gas + radiation row (the line's) | TR-BDF2 | 100: 80.0 / 2.99 / 16; 120: 125.2 / 2.21 / 15; 140: 166.3 / 1.63 / 14; 160: 205.0 / 2.40 / 16 | 3.6e-5 / 3.56e-5; 2.7e-5 / 2.70e-5; 2.3e-5 / 2.26e-5; 2.0e-5 / 1.98e-5 | 0 |
| en_s5c_gr | GASROW, EXCAP=1e8 | TR-BDF2 | 60 … 140: every printed digit of en_s5c_ref | the same | 0 |
| en_s5c_ec | ECOMPAT, EXCAP=1e8 | TR-BDF2 | 100: 87.5 / 3.06 / 16; 120: 130.7 / 2.26 / 16; 140: 169.8 / 1.67 / 14 | 3.0e-5 / 2.96e-5; 2.3e-5 / 2.32e-5; 2.0e-5 / 1.97e-5 | 0 |
| en_s5_ec | ECOMPAT, no EXCAP | TR-BDF2 | 400: 34.3 / 0.13 / 8; 440: 38.4 / 0.14 / 2 | 9.4e-2 / 1.53e-1; 1.0e-1 / 1.63e-1 | 0 |
| en_s5b_ref | combined row | backward Euler, EFRAC=5, DXMAX=4e-3 | 140: 3.68 / 0.20 / 20; 160: 24.3 / 4.80 / 19 | 2.0e-3 / 2.04e-3; 3.0e-2 / 2.95e-2 | 0 |
| en_s5b_gr | GASROW, no EXCAP | backward Euler, same | 200: 2.59 / 0.019 / 12 | 4.0e-4 / 4.00e-4 | 11 × SingularException(6) |
| en_s5b_ec | ECOMPAT, no EXCAP | backward Euler, same | 180: 2.29 / 0.013 / 10 | 4.0e-4 / 4.03e-4 | 12 × SingularException(6) |

* **With the exchange cap the compatible rows run the isolated star like the line's own row**: same Δt, same Newton iterations, the
  ledger closing equally (2.0e-5 at 169.8 yr against 2.3e-5 at 166.3 yr). The model after the relaxation differs by 0.6 % in
  L (R = 6.5287 R☉, log L = 6.1844 against 6.5262, 6.1818): the relaxation is a dynamical phase and the rows differ in how they
  dissipate it.
* **Without the cap the gas-alone row cannot be used in an opaque interior, with or without ECOMPAT** (the control en_s5b_gr fails
  exactly like en_s5b_ec). c κ ρ Δt at the centre is 3e11 s⁻¹ × Δt: 2.8e17 at Δt = 0.03 yr, above 1/ε = 4.5e15. The gas row and the
  radiation row are then, to round-off, the same row with opposite sign in their T and E columns; the information of their sum,
  the thermal balance, is below round-off. Backward Euler meets a zero pivot; TR-BDF2 raises no exception, converges in 2–8
  iterations at Δt 0.1 yr and opens the ledger by 10 %. The combined row of the line (the sum of the two, with
  the exchange cancelled analytically) does not have the problem, and EXCAP exists for it (its docstring in TitanResidFast.jl); it is a property of GASROW, not of the compatible terms.
* GASROW + EXCAP reproducing the combined row digit for digit is the expected result (a row operation on the Jacobian), and a
  check that the cap does not change the solution at this level.
* **What would have to be run** before the star line adopts the compatible rows: the full S5 sequence at 100 zones to the ZAMS and
  along the main sequence, current rows against ECOMPAT + EXCAP=1e8 (ZAMS t, R, L, T_c, ρ_c; the ledger with burning, where
  Σ|E_bur| is no longer small); the quasi-star line the same way. Hours each.
* A general alternative to the cap, not built: replace the gas row zone by zone by (gas row + α × radiation row) with the exchange
  term cancelled analytically in the α part, α → 1 where c κ ρ Δt is large. That is a row operation: it changes neither the
  discrete solution nor what is conserved, only the conditioning, and would make EXCAP removable (§11).

## 8. From the initial condition (task 5)

Run t_en_ic1: tree e0e6b2e, the run.env and phaseB.env of z100mm (IC abn02red, 100 zones) plus `ECOMPAT=1 EDIAG=2`, nothing else
changed. Stopped by hand at step 1478 (M* = 0.33 M☉, 5.3 yr after the core forms) to free the cores; it had not failed. The
comparison is the production run z100mm itself (its own, older tree); a controlled ECOMPAT-off run from the IC with tree e0e6b2e
was not made (it needs ~6000 steps to the same mass). *Made in the follow-up: t_en_ic0, §16 (it needed 2100 steps, not 6000). The
statements below that it changes are marked.*

**How far, at what cost** (perf.csv; failures = trapezoidal + BDF Newton failures since step 0):

| M* [M☉] | t_en_ic1: step / r* [R☉] / T_c [K] / log n_c / failures | z100mm: step / r* / T_c / log n_c / failures |
|---|---|---|
| core appears | 200 (t = 2.410650e5 yr) | 200 (t = 2.378241e5 yr) |
| 0.17 | 700 / 32.5 / 7.99e4 / 22.90 / 43 | 1400 / 30.7 / 1.090e5 / 23.16 / 62 |
| 0.21 | 800 / 38.3 / 8.12e4 / 22.91 / 50 | 2100 / 34.7 / 1.132e5 / 23.18 / 96 |
| 0.26 | 1100 / 44.4 / 8.28e4 / 22.92 / 86 | 3900 / 45.8 / 1.160e5 / 23.19 / 222 |
| 0.32–0.35 | 1400 (M 0.331, 5.26 yr after the core) / 50.9 / 8.44e4 / 22.93 / 120 | 5300 (M 0.321, 4.81 yr) / 54.7 / 1.190e5 / 23.20 / 295; 6500 (M 0.349, 5.40 yr) / 59.6 / 1.186e5 / 23.19 / 335 |

* Same accretion history (0.33 M☉ 5.3 yr after the core in both), **four times fewer steps** (1200 against 5100–6300 after the
  core) and about three times fewer Newton failures (95 against 258–298 after the core). *Corrected by §16:* against the control
  on the same tree it is a third fewer steps (1200 against 1800) and 95 against 115 failures; the factor four is mostly the
  difference between this tree and settings and those of z100mm's first 13.5k steps, not the rows.
* A different star: at equal mass T_c is 0.84e5 K against 1.19e5 K and n_c 10^22.93 against 10^23.20 cm⁻³; the radius is larger at
  low mass (38.3 against 34.7 R☉ at 0.21 M☉), equal near 0.25 M☉, smaller at 0.33 M☉ (50.9 against ~56 R☉). §7 has the cause.
  *§16:* T_c and n_c at equal mass are the rows' (the control reproduces z100mm's); the radius is not settled (the control has
  45.0–46.8 R☉ at 0.31–0.34 M☉, further from z100mm than ECOMPAT is).
* The core forms at 2.41065e5 yr against 2.37824e5 yr (+1.4 %). Not attributed: the two runs differ in ECOMPAT *and* in the code
  tree, and the collapse phase was not analysed. *§16:* the control forms its core at 2.42762e5 yr (+2.08 %); the delay is the
  tree's and settings', and ECOMPAT advances the core by 0.70 % against its control.

**Ledger** (EDIAG=2, 2914 converged stages of steps 0–1430; Σ over zones [L☉], median / 90 % / max of |.|):

| | all stages | star phase only (2550 stages, steps 181–1448) |
|---|---|---|
| Σ over all zones of T − G5 − adv + D/Δt | −0.000 / 1.07 / 11.2 | −0.001 / 1.27 / 11.2 |
| per-zone identity `cres`, zones 1–87 | 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 |
| per-zone identity `cres`, zones 88–100 | 0.000 / 0.87 / 11.2 | — |
| Newton remainder of the solved rows, zones 1–26 | 0.000 / 0.005 / 1.43 | — |
| Newton remainder, zones 27–87 | ≤ 0.012 | — |
| stages with |Σ| > 1 L☉ | 304 of 2914 | 312 of 2550 |
| time-weighted Σ | −0.0000 L☉ | |

From the IC the Newton remainder of the star's rows is below 1.5 L☉ at every stage (the ±100 L☉ of the restart test belongs to
the 0.9–1.3 M☉ star of the checkpoint); all stages above 1 L☉ come from zones 88–100, where the lab-frame radiation row
differences 1e53 erg of stored CMB energy (per-zone |cres| up to 11 L☉, zero in the median). Gate (a) "below 1 L☉ per stage" is
thus met in 90 % of the stages and missed by round-off in the rest. In the star phase Σ T over zones 1–26 is 2210 L☉ in the
median (2664 at 90 %): G5, while L_em is only 0.5–2.5e3 L☉.

**Radiation in the cloud (the gate of §3.3).** E/(a T_CMB⁴) − 1 at zones 88 / 92 / 96 / 99 / 100: 5.3e-6 / 3.5e-6 / 2.3e-6 / 4.4e-6 /
4.3e-6, constant from the first record to step 1350 (u_100 = −1.14 km/s, u/c = 3.8e-6): uniform to O(u/c), **gate met**. L_F at
the outer boundary (interface 101): +3e-3…8e-3 L☉ before the core, +60…100 L☉ at steps 150–300, +390…1730 L☉ at steps 350–1350
(positive: the star's light leaving); L_F at interface 100 is 1.187–1.222e5 L☉ = the 1.158e5 L☉ an isotropic lab-frame field
shows the infalling interface, plus the star. z100mm at the same stage: E excess 0.03…0.63 in zones 88–99 and 33–44 at zone 70,
L_F(101) = −1.16e5 L☉. The −1e5 L☉ boundary artefact and the stored excess do not form.

**Photosphere** (τ = 2/3 from outside; L_em at 1.5 r_ph):

| step | 350 | 550 | 750 | 950 | 1050 | 1250 | 1350 |
|---|---|---|---|---|---|---|---|
| zone | 17 | 18 | 19 | 19 | 20 | 20 | 21 |
| r_ph [R☉] | 33.4 | 39.2 | 44.6 | 46.0 | 52.7 | 53.8 | 59.7 |
| L_em [L☉] | 506 | 836 | 1170 | 1494 | 1860 | 1900 | 2475 |
| T_eff [K] | 4766 | 4972 | 5067 | 5380 | 5236 | 5245 | 5279 |

Lowest gas temperature of the outer zones: 186.8 K (zone 83); no zone is heated or cooled anomalously.

**Open point found here** (§10): in zones 64–72, between the star-dominated and the CMB-dominated radiation field, E falls *below*
a T_CMB⁴ as the star's luminosity rises (−27 % at zones 66/67 at step 1371). *§15:* it is the flux-tie closure's, deeper in the
control (−51 %) and in z100mm (−99.98 %); and the "L_F at interface 100 = 1.187–1.222e5 L☉" above is F = 5e-6 c E at 19 pc, the
bath's tolerance level times the area, not a luminosity.

## 9. Multigroup (task 6)

### 9.1 Change and off-path identity (commit 50f1609)

julia/MGRad.jl `mg_rows!`: with TIEFRAME each group's tie row closes the comoving group flux with the discrete form of §3.1 (line 339:
the group's own h_g, f_g, E_g and advective flux; SLABTIE's donor form included), and each group's outer boundary row holds for the
lab-frame flux (line 353).

Off-path: MGROUPS=4 with TIEFRAME off on tree e0e6b2e (t_en22) against the branch point 4f1fe93 (t_en23). The 4 records both runs
have (steps 17405–17420), 27 fields each (r, ρ, T, E, F, χ, κ_E, κ_P, f, the per-group E_g, F_g, opacities and closure factors, level
populations, ionisation), are bit-identical. One field differs and is not state: the diagnostic buffer `MGHPLUS` (MGRad.jl:57, written
to mg_rec.jls by popiii_m2c.jl:1538) in its three ghost elements 1, 2, 103 of 103, which are never assigned (values such as 2e-316,
8e173). They should be initialised; a note for the multigroup branch, nothing was changed here.

### 9.2 The optical group's outflow is not a frame artefact

docs/POPIII_MULTIGROUP.md §14 reports 4.0–4.8e5 L☉ leaving through the outer boundary in the optical group with no emission to feed
it ("a draining reservoir"). Pair from the checkpoint, MGROUPS=4, 200 steps, tree e0e6b2e; optical group, zones 94–100 (optout.jl):

| | t_en21 (TIEFRAME=1) | t_en22 (TIEFRAME off) |
|---|---|---|
| time covered by 200 steps | 1.52 yr | 0.047 yr |
| stored U_opt, start → end | 7.5074e49 → 7.4802e49 erg | 7.5074e49 → 7.5073e49 erg |
| emission − absorption of the group | −1.4e-4 L☉ | −1.4e-4 L☉ |
| L_out at t − t₀ = 1.1–1.2e5 s | 1.14e4 L☉ | 1.01e4 L☉ |
| at 4.3–4.8e5 s | 3.80e4 (4.76e5 s) | 3.40e4 (4.30e5 s) |
| at 1.0e6 s | 7.8e4 (0.98e6 s) | 7.9e4 (1.02e6 s) |
| at 1.37e6 s | | 1.09e5 |
| at 2.06e6 s | 1.62e5 | |
| at 1.11e7 s | 8.5e5 | |
| at 4.42e7 s (1.40 yr) | 3.12e6 (dU_opt/dt = −2.8e6) | |

**The outflow is not gone, and the frame-consistent tie does not change it**: at equal time the two runs have the same L_out, rising
roughly in proportion to the time since the restart (≈ 1e4 L☉ per 1e5 s, still rising after 1.4 yr). The frame terms are O(u/c) ~ 1e-6
of c E and cannot matter at this level. What flows out is stored energy: the restart's partition of the grey field puts 7.5e49 erg
into the optical group in zones 94–100 (2 % of the grey field's 3.2e51 erg excess over a T_CMB⁴ in the same zones); nothing emits
or absorbs in that group out there, so the field leaves through the boundary. U_opt/L_out = 200 yr at the end of t_en21; the light
crossing time of the boundary radius is 74 yr. §14's 4.0–4.8e5 L☉ is a point on this ramp, not a steady luminosity. "Draining
reservoir" is the right description; what this work adds is where the reservoir comes from: it is part of the excess over the CMB
that the comoving energy row built in the tied zones over the checkpoint's 2.4e5 yr history (§3). From the initial condition with
the compatible rows that excess does not form (§8: E/(a T_CMB⁴) − 1 = 2–5e-6 in zones 88–100). **Not run**: multigroup from an initial
condition or checkpoint that has no excess; that is the test that would show the optical group quiet.

Total imbalance of the dual-energy ledger (EDIAG=1, median 4π·imb/Δt over the stages): 3831 (TR) / 3821 (BDF) L☉ with TIEFRAME,
3844 / 3850 L☉ without. Unchanged, as in grey (t_en06, §3.2): the creation is in the energy row.

Star side, same pair (greyana.jl, every 6th step):

| | t_en21 (TIEFRAME=1) | t_en22 (off) |
|---|---|---|
| photosphere zone, r_ph | 26 → 27 (by 0.68 yr), 111.4 → 106.3 R☉ | 26, 110.8–111.8 R☉ |
| L_em(1.5 r_ph), median (rms of ln; range) | 5986 L☉ (0.053; 5351–7025) | 5990 L☉ (0.055; 5353–7148) |
| T_eff, median | 4331 K | 4333 K |
| L_F at the outer boundary, start → end | −1.02e5 → +3.02e6 L☉ | −1.02e5 → +4.7e3 L☉ |
| median Δt | 4.3e4 s | 7.1e3 s |
| Newton failures in 200 steps | 90 | 55 |

One pair; the factor 6 in Δt is not attributed (the runs cover different times, t_en21 passes a photosphere zone change).

### 9.3 What is missing for multigroup with the compatible rows

ECOMPAT's lab-frame radiation energy row (§3.3) exists for the grey row only. The group energy rows of `mg_rows!` are comoving, with
the work term f_g E_g Δ(r²u) + (1 − 3f_g) E_g (u/r) ΔV, i.e. each group creates its share of R6a where its flux row is the closure.
The lab-frame group row was not written in this round. Its frequency-integrated sum must reduce to §3.3 (flux of u f E through the
interfaces, work u χ̄ ρ F/c on the gas); the terms that move energy between groups at O(u/c) have to be derived and checked against
the references of §12, not written from memory. popiii_m2c.jl does not refuse ECOMPAT=1 with MGROUPS > 0; that combination is
untested and cannot be conservative in the tied zones. This is the first item for the merge with the multigroup branch.

## 10. Negative results and open points

Negative results:

1. **H2 as stated.** The frame-consistent closure alone does not remove the 3280 L☉ of zones 88–100 (t_en06) and does not change the
   multigroup imbalance (3.83e3 L☉ with and without, t_en21/t_en22) or the optical group's outflow (§9). It is correct and needed;
   the creation is in the comoving energy row (§3).
2. **K3a returned where it is computed**: a heat dipole across the photosphere, 0.20 yr in 500 steps (t_en08).
3. **Booking the dissipated kinetic energy instead of thermalising it**: no gain in Δt or failures, a ledger open by 190 L☉, the
   same star (t_en20).
4. **30-step tests from the checkpoint** ranked the variants wrongly (§6.1 against §6.2).
5. **Gate (a), below 1 L☉ at every stage: not met.** The identity closes per zone to 0.000 L☉; Σ over the zones exceeds 1 L☉ in a
   third of the stages from the checkpoint (max 107 L☉: what the Newton stopping criterion leaves in the star's rows) and in 10 %
   of the stages from the IC (max 11 L☉: round-off of zones 88–100, 9.8e52 erg of stored radiation energy differenced over the stage
   Δt). Time-integrated +0.0009…+0.026 L☉ and −0.0000 L☉. The gate was not weakened. Obstructions: the stopping test is on the
   corrections, not on the energy rows' residuals; and the outer zones hold 1e53 erg in a row that has to resolve 1 L☉.
   *Follow-up (§17):* with a residual-based stopping test (ERESTOL=1) 48 of 1000 stages stay above 1 L☉ (largest 11.4 L☉), 81
   with the round-off of zones 88–100 included. Still not met.
6. **Gate (c), Δt not smaller: met in the mean and in the time covered, not in the median** (§5c).
7. **The stored excess of the checkpoint does not relax within the tests** (E/E_ref − 1 = 0.63 … 0.03 in zones 88–99 unchanged
   after 7–12 yr; r/c = 74 yr). The gate "uniform to O(u/c)" is shown from the IC (§8), not from the checkpoint.
8. **The optical group's reservoir is not removed** (§9.2).
9. **The gas-alone row, compatible or not, needs EXCAP in an opaque interior** (§7.3).

Open points:

10. **E below a T_CMB⁴ in zones 64–72 of the from-IC run.** At step 1371 E/(a T_CMB⁴) − 1 = −0.024, −0.151, −0.266, −0.265, −0.204,
    −0.161, −0.071, −0.015, −0.002 in zones 64 … 72; L_F dips from 2.48e3 to 1.69e3 L☉ at zone 67; the deficit at zone 70 grew from
    −3e-4 (step 550) to −7e-2 (step 1350); the tie weight is 1 there. Not understood. A candidate, unproven: the closure
    F = h c E_donor takes h from the static formal solution while L* grows, so the transition between the star-dominated and the
    CMB-dominated field is a closure matter (G4), which in z100mm is hidden under the excess. It is small in energy (these zones
    are optically thin, the gas does not respond: lowest T 186.8 K) but it is a radiation energy density below the bath's.
    *Follow-up (§15): cause found, the candidate confirmed by a model of the energy row alone. It is the flux tie's (the control
    has −51 %, z100mm −99.98 %), not ECOMPAT's. The fix tested (TIEX) is not general; not fixed.*
11. **The star differs from the current stack's** (§7.2); attributed to G5 by exclusion, not by a controlled run, and which gravity
    power converges faster in zone number is not known. *Follow-up (§14): attributed by a controlled run to 55–76 % (t_en25);
    an explicit binding energy follows the flux form, not the momentum row's power; no smooth force conserves it to round-off;
    the mesh weight tested against G5 gains 2 % at the strength that leaves the star alone and makes G5 larger otherwise (§14.6).
    In the star's zones the stack's total row is the one that is closed against the binding energy (§14.7). Convergence in zone
    number is still not measured.*
12. **The core forms 1.4 % later from the IC** (2.41065e5 against 2.37824e5 yr). Not attributed: ECOMPAT and the code tree both
    differ from z100mm, and no ECOMPAT-off run from the IC with tree e0e6b2e was made (~6000 steps). *Follow-up (§16): made. The
    delay is the tree's (control +2.08 %); the T_c and n_c differences are the rows'.*
13. **t_mg25's failure** is explained by the measured terms (§6.3), not by repeating it with one term removed. *Follow-up (§14.5):
    reproduced with one term alone, the exact gravity term in every zone on top of ECOMPAT (t_en26).*
14. **Multigroup**: the group energy rows are still comoving; ECOMPAT with MGROUPS > 0 is untested and not refused by the driver;
    `MGHPLUS` has three uninitialised ghost elements (§9).
15. **Isolated-star and quasi-star lines** not run through with the compatible rows (§7.3).
16. A size check, no more: z100mm's boundary inflow of 1.16e5 L☉ over 2.38e5 yr is 3.3e51 erg; the excess over a T_CMB⁴ stored in
    zones 88–100 of the checkpoint is 3.5e51 erg (estore.jl). The inflow was not followed through the run's history.
17. Stale settings met on the way: `CONV=1e-3` in z100mm's run.env is not read by popiii_m2c; the `TRESET` row of
    docs/POPIII_SWITCH_INVENTORY.md refers to code that is gone.

## 11. What becomes removable once ECOMPAT is the default

Not deleted in this round (the S4/S5 isolated-star line and the quasi-star line use the same residual). *Follow-up: whether
ECOMPAT becomes the default is open again (§14.7). The lists below hold if it does, unchanged; if ECOMPAT + ECGRAV=1 is chosen
instead, `dualswitch`, `ETASW`, `ETASWTH`, `ETASWTIE` and the `Ctx` field `useint` stay (ECGRAV=1 selects its zones with them)
and only the total-row branch of the solve goes. Line numbers in this section are those of e0e6b2e; on bb57e3e they are larger by
up to about 50. The provisional switches of the follow-up are listed at the end of this section.*

**Run settings (popiii_m2c environment):**

| Setting | Where set | Why it goes |
|---|---|---|
| `ETOTAL=1`, `ETASW=0.6` | phaseB.env; popiii_m2c.jl:1208; popiii_stack.jl:37 | no total row, no switch (gate b) |
| `ETASWTH=1e-3`, `ETASWTIE=0.5` | run.env; popiii_m2c.jl:1205–1206, 1311–1313; popiii_stack.jl:37 | selection criteria of the switch |
| `ETOTGRAVGRID` | popiii_m2c.jl:1207 | missing gravity term of the total row; §4.3 has no such term |
| `TRESET` | docs/POPIII_RUNS.md, inventory row only | the one-time T repair of the tie rows' gas-energy artefact (§3); its code is already gone from this branch, the inventory row is stale |
| `CONV=1e-3` | run.env of z100mm | already dead: popiii_m2c does not read it |
| `ECOMPAT`, `TIEFRAME` | popiii_m2c.jl:1237–1247 | the provisional switches themselves, once default |

**Code (julia/TitanResidFast.jl unless stated):**

* constants `ETOTAL`, `ETASW`, `ETASW_EOFF`, `ETASW_TH`, `ETASW_TIE` (lines 63–81), `ETOT_GRAV_GRID` (1661);
* `dualswitch` (332–357), the `Ctx` fields `etotal` and `useint` (53–54) and the `ηsw` argument of `context` (359, 390);
* the branch `elseif c.etotal && !c.useint[k]` of the gas energy row (1417). `etot_row` (1668) and `ETOT_GRAV` stay, as the
  independent diagnostic of `zone_diag` (the conservative row the compatible rows are checked against) — no longer a row of a solve;
* for this stack, the comoving radiation energy row where the closure replaces the moment equation (the `else` of 1503) and the
  static-frame forms of the tie rows: FLUXTIE donor row (the `else` of 1571), outer boundary row (1583), and the same two in
  julia/MGRad.jl (339, 353), SLABTIE included;
* `nocompat` (1664) once `zone_diag` builds its reference rows directly;
* drivers that exist only for the dual-energy experiment: julia/etotal_run.jl, julia/etot_check.jl.

**Stays, and why:**

* `GASROW` and the branch `elseif c.lrad && GASROW[]`: the quasi-star line (julia/qstar_run.jl `physics.gasrow`) and the test
  drivers stepdiag.jl, newtrace.jl, compton_box.jl, gasrow_check.jl, order_test.jl select it. ECOMPAT reuses its advection of the gas
  energy (`W.egbg`, line 1159). If the quasi-star line adopts ECOMPAT, this branch goes too.
* The combined TITAN gas + radiation row (`elseif c.lrad`): the S4/S5 isolated-star line. §7.3: with EXCAP the compatible rows
  run that star's first 170 yr with the same Δt, iterations and ledger closure; the line was not run to the ZAMS with them, so the
  combined row stays until that is done.
* `EXCAP=1e8`: not an energy-formulation switch (the same factor in both rows conserves). It bounds the round-off of the exchange
  term, and the gas-alone row cannot run an opaque interior without it (§7.3: singular Jacobian at c κ ρ Δt > 1/ε). The general
  replacement to try is the row operation of §7.3 (gas row + α × radiation row with the exchange cancelled analytically).
* julia/TitanResid.jl and julia/TitanResidRef.jl carry their own `etotal` (the reference residuals of the 2025 port tests); untouched.

docs/POPIII_SWITCH_INVENTORY.md rows EXCAP, ETOTAL, ETASW, ETASWTH, ETASWTIE, ETOTGRAVGRID, TRESET (§9 there) and its goal G5
("one conservative energy formulation") are the entries to update when the default changes.

**Provisional switches added in the follow-up, all default off.** With all of them off, tree bb57e3e reproduces e0e6b2e
bit for bit over 20 steps from the checkpoint, with the current rows (t_en33 against t_en18) and with ECOMPAT (t_en34 against
t_en19): 40 stage records, 2320 fields compared, 0 differ (`bitcmp.jl`). Line numbers of bb57e3e:

| Switch | Where | Result | For a clean code base |
|---|---|---|---|
| `ECGRAV=1` | TitanResidFast.jl:644–655, 1453, 1935; popiii_m2c.jl:1248–1251 | works; best ledger against the binding energy; needs the dual-energy flag (§14.5, §14.7) | keep only if this combination is chosen; otherwise remove |
| `ECGRAV=2`, `3` | same | fail in the cold infall (t_en32, t_en26) | remove |
| `ECGRAV=4` | same | not run | remove |
| `TIEX` (`HXEDD`) | TitanDriver.jl:141–150; TitanCore.jl:770–771; TitanResidFast.jl:1611–1613, 1942; popiii_m2c.jl:1342–1345 | the right limit once the star dominates the field; fails from the IC (§15.4) | remove once a general closure exists (§15.5); not a default |
| `ERESTOL` | Trbdf2.jl:70–106; popiii_m2c.jl:1347–1351, 1900–1901 | general; lowers the unmet energy residual, does not meet gate (a) (§17) | the only candidate to keep |
| `GEQGP` | TitanResidFast.jl:1032–1038; grideq_driver.jl:40–42, 141; popiii_m2c.jl:1193 | negative (§14.6) | remove |

**Provisional switches added in round 3, all default off.** With all of them off, tree 6642493 reproduces the trees before it
bit for bit over 20 steps from the checkpoint (`bitcmp.jl`, 40 stage records each): current rows (t_bi1 against t_en18 and
t_en24: 2320 fields, 0 differ), ECOMPAT (t_bi2 against t_en34: 2400 fields, 0 differ; against t_en19: 2320, 0 differ),
ECOMPAT + ECGRAV=1 (t_bi3 against t_en25: 2360 fields, 0 differ); and `ERLAB=1` on 6642493 reproduces d676f9b (t_bi4 against
t_en35: 2400 fields, 0 differ). Line numbers of 6642493:

| Switch | Where | Result | For a clean code base |
|---|---|---|---|
| `ERLAB` | TitanResidFast.jl:657–661, 1611, 1775, 1901; popiii_m2c.jl:1253–1260 | the lab-frame radiation row with the stack's gas rows: the stack's star, the outer cloud's 3.3e3 L☉ gone, same robustness (§18) | the candidate to make the default of the stack; the comoving branch then goes as listed above |
| `GFV` | TitanResidFast.jl:663–672, 695–711, 1402, 1463, 1936; popiii_m2c.jl:1261–1265 | free fall correct, infall's G5 −18 L☉ from the checkpoint and −8 L☉ at 200 zones from the IC; no binding energy; at 100 zones its error on the coarse central mesh gives a star 0.33–0.35 dex less dense at the centre, at 200 zones the node weight's star to 0.02–0.03 dex (§19, §21) | not a default at 100 zones (§19.6); keep only if the mesh is to resolve the centre (ratio of successive radii ≤ 1.5 there), otherwise remove |
| `GVAR` | TitanResidFast.jl:674–700, same call sites; popiii_m2c.jl:1266–1271 | the variational weight: has a binding energy, free fall to 2–4 %, best oscillation frequencies; changes the star's central balance most (§19); from the IC at 100 zones its star is 0.43–0.46 dex less dense at the centre and 36–39 % cooler than the node weight's, and the ledger is −952 / −187 L☉ against W (i) / W (iv) (§19.5) | remove unless the weight question is reopened |

## 12. References (each checked on NASA ADS; formulas checked against the limits stated in the text)

* Bank, Coughran, Fichtner, Grosse, Rose & Smith 1985, IEEE Trans. Electron Devices 32, 1992 (1985ITED...32.1992B) — TR-BDF2.
* Benson 1992, Comput. Methods Appl. Mech. Eng. 99, 235 (1992CMAME..99..235B) — review of hydrocodes; DeBar's kinetic-energy fix.
* Bryan, Norman, Stone, Cen & Ostriker 1995, Comput. Phys. Commun. 89, 149 (1995CoPhC..89..149B) — dual-energy formalism.
* Caramana, Burton, Shashkov & Whalen 1998, J. Comput. Phys. 146, 227 (1998JCoPh.146..227C) — compatible (adjoint) energy equation.
* Gehmeyr & Mihalas 1994, Physica D 77, 320 (1994PhyD...77..320G) — TITAN.
* Hosea & Shampine 1996, Appl. Numer. Math. 20, 21 (1996ApNM...20...21H) — TR-BDF2 analysis, γ = 2 − √2.
* Krumholz, Klein & McKee 2007, ApJ 667, 626 (2007ApJ...667..626K) — frame transformations of the radiation moments to O(u/c).
* Lowrie, Morel & Hittinger 1999, ApJ 521, 432 (1999ApJ...521..432L) — coupling of radiation and hydrodynamics, frames.
* Mihalas & Auer 2001, JQSRT 71, 61 (2001JQSRT..71...61M) — Lagrangian/Eulerian radiation moment equations.
* Mihalas & Mihalas 1984, Foundations of Radiation Hydrodynamics, Oxford (1984frh..book.....M) — §91.

Added in the follow-up (§14, §15), each found as refereed on NASA ADS/SciX:

* Springel 2010, MNRAS 401, 791 (2010MNRAS.401..791S; arXiv:0901.4107 read: §5.2, §5.4 eqs. 92–96) — gravity source terms on a
  moving mesh: the energy source built from the mass flux, and why the cell-centred form is not conservative.
* Jiang, Belyaev, Goodman & Stone 2013, New Astron. 19, 48 (2013NewA...19...48J) — gravitational energy flux on a fixed mesh.
  Abstract only read.
* Mullen, Hanawa & Gammie 2021, ApJS 252, 30 (2021ApJS..252...30M) — the same to round-off, with mesh refinement. Abstract only read.
* Katz, Zingale, Calder & Swesty 2016, ApJ 819, 94 (2016ApJ...819...94K) — energy non-conservation of standard gravity coupling in
  white-dwarf mergers. Abstract only read.
* Stone, Mihalas & Norman 1992, ApJS 80, 819 (1992ApJS...80..819S); Davis, Stone & Jiang 2012, ApJS 199, 9 (2012ApJS..199....9D);
  Jiang, Stone & Davis 2012, ApJS 199, 14 (2012ApJS..199...14J) — variable Eddington tensor closures from a formal solution
  (§15.5). Abstracts only read; no formula of theirs is used.

Added in round 3 (§19), each found as refereed on NASA ADS/SciX; abstracts only read, no formula of theirs is used (the weights
of §19.1 are derived there and checked against the limits computed in §19.2–§19.3):

* Zingale, Dursi, ZuHone, Calder, Fryxell, Plewa, Truran, Caceres, Olson, Ricker, Riley, Rosner, Siegel, Timmes & Vladimirova
  2002, ApJS 143, 539 (2002ApJS..143..539Z; arXiv:astro-ph/0208031) — hydrostatic equilibrium on a finite-volume mesh: the
  discrete pressure gradient and the discrete weight have to be the same quadrature or the star drifts.
* Käppeli & Mishra 2014, J. Comput. Phys. 259, 199 (2014JCoPh.259..199K) and 2016, A&A 587, A94 (2016A&A...587A..94K) —
  well-balanced finite-volume schemes with gravity: exact discrete hydrostatic balance needs a reconstruction built on the
  equilibrium, not a particular quadrature weight of g ρ.

Added in round 3 (§21.3), found as refereed on NASA ADS/SciX; not read beyond title and abstract listing, no formula of
theirs is used (the interface densities of §21.3 are computed there for power laws and measured on the runs); they are
named as the places where a reconstruction consistent with the mesh geometry is constructed:

* Mignone 2014, J. Comput. Phys. 270, 784 (2014JCoPh.270..784M) — conservative reconstruction for finite-volume methods in
  cylindrical and spherical coordinates.
* Colella & Woodward 1984, J. Comput. Phys. 54, 174 (1984JCoPh..54..174C) — interface values from the cumulative profile
  (PPM); van Leer 1977, J. Comput. Phys. 23, 276 (1977JCoPh..23..276V) — the monotonic slope the code's advection uses.

## 13. Runs, trees, scripts

All run directories under /zpool/nvme/data/tabel_scratch/titan_popiii/runs/, trees under
/zpool/nvme/data/tabel_scratch/titan_popiii/mg_trees/. Restart of `z19z100mmS5_auto0` at step 17398 unless stated; 100 zones.

| run | tree | settings beyond z100mm's | steps | used in |
|---|---|---|---|---|
| t_en01 | 774c85b | EDIAG=2 (baseline) | 60 | §2 |
| t_en06 | a37f48b | TIEFRAME=1 | 500 | §3.2 |
| t_en07 | a37f48b | TIEFRAME off | 20 | §3.2, bit identity against t_en01 |
| t_en08 | be38295 | ECOMPAT=1, K3a returned in place | 500 | §4.4, §6.2 |
| t_en09–t_en16 | be38295_dbg | one term of the row at a time | 30 each | §6.1 |
| t_en17 | e0e6b2e_dbg | transported form, per-iteration trace | 30 | §6.1 |
| t_en18 | e0e6b2e | ECOMPAT off | 20 | §5, bit identity against t_en01 |
| t_en19 | e0e6b2e | ECOMPAT=1 (the final form) | 500 | §5, §6.2, §7 |
| t_en20 | e0e6b2e_dbg | ECOMPAT=1, dissipated kinetic energy booked, not returned | 500 | §6.2, §7 |
| t_en21 | e0e6b2e | MGROUPS=4, TIEFRAME=1 | 200 | §9 |
| t_en22 | e0e6b2e | MGROUPS=4 | 200 | §9 |
| t_en23 | 4f1fe93 | MGROUPS=4 | 20 | §9.1, bit identity against t_en22 |
| t_en_ic1 | e0e6b2e | from the initial condition, ECOMPAT=1 | 1478 (stopped; last perf.csv row at step 1400) | §8, §15, §16 |
| t_en24 | e0e6b2e | EDIAG=2 (current stack with per-stage records: the control for §14) | 500 | §14 |
| t_en25 | 39923e1 | ECOMPAT=1 ECGRAV=1 | 500 | §14.2, §14.5 |
| t_en26 | bb57e3e | ECOMPAT=1 ECGRAV=3 | 447 (failed at step 17845) | §14.5 |
| t_en27 | bb57e3e | ECOMPAT=1 GEQGP=1,100 | 500 | §14.6 |
| t_en28 | bb57e3e | ECOMPAT=1 ERESTOL=1 | 500 | §17 |
| t_en29 | bb57e3e | ECOMPAT=1 TIEX=1 | 500 | §15.4 |
| t_en30 | bb57e3e | ECOMPAT=1 GEQGP=1,1000 | 500 | §14.6 |
| t_en31 | bb57e3e | ECOMPAT=1 GEQGP=1,300 | 500 | §14.6 |
| t_en32 | bb57e3e | ECOMPAT=1 ECGRAV=2 | 500 | §14.5 |
| t_en33 | bb57e3e | EDIAG=2 (new switches off) | 20 | §11, bit identity against t_en18 |
| t_en34 | bb57e3e | ECOMPAT=1 EDIAG=2 (new switches off) | 20 | §11, bit identity against t_en19 |
| t_en_ic0 | e0e6b2e | from the initial condition, ECOMPAT off (the control of §16) | 2100 (stopped) | §16 |
| t_en_ic2 | bb57e3e | from the initial condition, ECOMPAT=1 TIEX=1 | 238 (stopped) | §15.4 |
| en_s5_ref / _ec | e0e6b2e_s5 | isolated 100 M☉ star, TR-BDF2: current rows / ECOMPAT (no EXCAP) | 160 / 440 (stopped) | §7.3 |
| en_s5b_ref / _ec / _gr | e0e6b2e_s5 | the same with backward Euler, EFRAC=5, DXMAX=4e-3: current rows / ECOMPAT / GASROW (no EXCAP) | 140–180 (stopped) | §7.3 |
| en_s5c_ref / _ec / _gr | e0e6b2e_s5 | TR-BDF2: current rows / ECOMPAT + EXCAP=1e8 / GASROW + EXCAP=1e8 | 120–160 (stopped) | §7.3 |

t_en02–t_en05 were 3–20-step debugging runs of the zone_diag split (no result quoted; moved to runs/outdated/). Reference for every comparison with the
current stack: the production run z100mm (read only) and its perf.csv. All of the above are finished or were stopped by me; nothing of
mine is running.

Trees with debugging toggles that are not in the branch (`*_dbg`: one-term switches and the booked variant; `*_s5`: the three environment
toggles ECOMPAT / GASROWONLY / EXCAP at s5_run.jl:20–22) differ from the named commit only by those lines.

Scripts (all in /tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/en/): mktree.sh, mkrun.sh, mgrun.sh, s5run.sh, killtag.sh, killtag5.sh (harness);
egate.jl (gate a/b: per-zone identity and sums from ediagz.jls), eledger.jl (time-integrated ledger), kstar.jl, subk.jl, kdiss.jl
(kinetic-energy terms), greyana.jl (photosphere, L_em, T_eff, boundary flux), tprof.jl, relax.jl, outer2.jl, estore.jl (outer zones:
E against a T_CMB⁴, stored energy), bitcmp.jl, mgbit.jl, mgbit2.jl (bit identity), newt.jl, ntrace.jl, dts.jl, dtmed.py (Newton
failures, Δt), icmp.jl, icmp2.jl, rhoc.jl, perfat.py (from-IC comparison at equal mass and time), imbmed.py (EDIAG=1 medians); and
optout.jl of the multigroup harness (same scratchpad, directory mg/). They are outside the repository.

Added in the follow-up, same directory: gledger.jl (global ledger against an explicit binding energy, four definitions of W),
g5split.jl (G5 per interface, geometric and density parts; its printed time is 0.879 × the true one, the weights are consistent
between runs), imbzone.jl (conservative total row's residual and its gravity terms per zone range), uff.jl (infall velocity
against free fall), cumt.jl, starhist.jl (photosphere history at given times), ebmin.jl, icdef.jl, icmodel.jl, outerL.jl, hxdiag.jl
(E against the bath, the model of the energy row alone, flux factors), unmet.jl (unmet residual per stage), mkicrun.sh (from-IC
run directory), waitfile.sh. Trees 39923e1 (ECGRAV), cd9922c (+ TIEX, ERESTOL), bb57e3e (+ GEQGP; contains all four).

**Round 3 (§18–§21).** Trees d676f9b (ERLAB, GFV) and 6642493 (+ GVAR; contains all switches); `6642493_s5` = 6642493 with the
environment toggles GFV / GVAR at s5_run.jl and W (iv) in the S5 ledger under GVAR. 4 threads per run unless stated. "ck" =
restart of `z19z100mmS5_auto0`; "IC" = from the initial condition with the settings of t_en_ic0.

| run | tree | start, zones | settings beyond the control's | steps | used in |
|---|---|---|---|---|---|
| t_en35 | d676f9b | ck, 100 | ERLAB=1 EDIAG=2 (8 threads) | 500 | §18, §19.4 |
| t_en36 | 6642493 | ck, 100 | ECOMPAT=1 GVAR=1 EDIAG=2 | 500 | §19.4 |
| t_en37 | 6642493 | ck, 100 | ECOMPAT=1 GFV=1 EDIAG=2 | 500 | §19.4 |
| t_en38 | 6642493 | ck, 100 | ERLAB=1 GVAR=1 EDIAG=2 | 500 | §19.4 |
| t_en39 | 6642493 | ck, 100 | ERLAB=1 GFV=1 EDIAG=2 | 500 | §19.4 |
| t_bi1 … t_bi4 | 6642493 | ck, 100 | EDIAG=2 / ECOMPAT=1 / ECOMPAT=1 ECGRAV=1 / ERLAB=1 | 20 each | §11, bit identity |
| en_s5_gfv / en_s5_gvar | 6642493_s5 | isolated 100 M☉ star, 100 zones | GFV=1 / GVAR=1 | 225 / 165 (stopped) | §19.5 |
| t_en_ic3 | bb57e3e | IC, 100 | ECOMPAT=1 ECGRAV=1 EDIAG=2 (8 threads) | 2100 (stopped) | §20, §21 |
| t_ic_l | 6642493 | IC, 100 | ERLAB=1 EDIAG=2 | 2150 (stopped at M\* 0.32: Δt 2e3 s since step 1850) | §20.4, §21 |
| t_ic_lf | 6642493 | IC, 100 | ERLAB=1 GFV=1 EDIAG=2 | 1800 (stopped itself at M\* 0.34) | §19.5, §21 |
| t_ic_lv2 | 6642493 | IC, 100 | ERLAB=1 GVAR=1 EDIAG=2 | ≈ 1580 (stopped at M\* 0.25) | §19.5 (4) |
| t_cv50s | 6642493 | IC, **50** | EDIAG=2 | 1700 (stopped at M\* 0.375) | §21 |
| t_cv50e | 6642493 | IC, **50** | ECOMPAT=1 EDIAG=2 | 700 (stopped at M\* 0.351) | §21 |
| t_cv50g2 | 6642493 | IC, **50** | ECOMPAT=1 ECGRAV=1 EDIAG=2 | 1400 (stopped itself at M\* 0.34) | §21 |
| t_cv200s | 6642493 | IC, **200** | EDIAG=2 | 507 (stopped: Newton failures before the core forms) | §21.1 |
| t_cv200l | 6642493 | IC, **200** | ERLAB=1 EDIAG=2 | 1830 (stopped at M\* 0.290: Δt 5e2–2e3 s since step 1600) | §21 |
| t_cv200e | 6642493 | IC, **200** | ECOMPAT=1 EDIAG=2 | 1200 (stopped at M\* 0.226: crawling) | §21 |
| t_cv200f | 6642493 | IC, **200** | ERLAB=1 GFV=1 EDIAG=2 | ≈ 1720 (stopped at M\* 0.30: Δt 7e2–3e3 s since step 1600) | §19.5, §21 |
| t_cv200g | 6642493 | IC, **200** | ECOMPAT=1 ECGRAV=1 EDIAG=2 | ≈ 800 (stopped at M\* 0.166: Δt 3e4 s) | §21.6, §21.7 |

The 50- and 200-zone runs are convergence tests (§21). Their meshes: `ENVEUL=0,d,d,1e30 ENVCORE=0.1,d,1e30` with d = 0.46991537
(50 zones) and 0.11570781 (200 zones); the 100-zone runs have d = 0.232584376. Run directories prepared and not run: t_ic_ev
(IC, ECOMPAT=1 GVAR=1), t_cv141s, t_cv141e (141 zones); t_cv50g was stopped at its start (run as t_cv50g2) and t_ic_lv was never started (run as t_ic_lv2). All of the above
are finished or were stopped by me; nothing of mine is running.

Scripts added in round 3, same directory: gledgerw.jl (gledger.jl with windows in M\* and the constant G for the new weights:
ENV MLO, MHI, G4FIX), g5int.jl, g5reg.jl (G5 per interface range and per region of the flow, for the three weights; the splits
upw/rem and geo/den of §19.4), imbreg.jl (total-row residual per region in a window of M\*), snapstar.jl (the star at equal
M\* from the snapshots), pmass.jl (structure at equal enclosed mass), lprof.jl (luminosity profile and radiated energy),
starmesh.jl, powstat.jl, polystat.jl (§19.3), gvartest.jl (closed forms of the variational force against quadrature), newtl.jl
(newt.jl for the lab-frame radiation row), and the run harness mgrun4.sh, mkicrun4.sh, icseq.sh, seq*.sh, ledwin.sh, winana.sh,
cvana.sh, ckana.sh, starall.sh, perfat.sh, waitend.sh; and for §21.3 ledfull.sh (the parts of the gravity power against
dW/dt from gledgerw.jl), g5grid.jl (g5reg.jl with the split of G5 into its mass-flux and mesh-motion parts per region),
wdrift.jl (the static test: W of a fixed smooth distribution under the runs' mesh motion, and the interface densities against
it), g5smooth.jl (G5 per region with the smooth distribution's interface density in place of the advected one) and vlbias.jl
(the interface densities for power laws on a geometric mesh).

## 14. The gravity power (follow-up to open point 11)

ECOMPAT leaves the conservative total row open by G5 (2.7e3 L☉ time-weighted over t_en19, 2655 L☉ in the median). §4.3 and §7.1 said
that the momentum row's Σ_j u_j g_j μ_j (μ_j = ½(m_j−1 + m_j), the node mass; "M_j" in §4) *is* the gravity power of this
discretisation and that G5 is a matter of choice between two consistent forms. **That was wrong**; §14.2–§14.3 replace it. Summary of
this section in §14.7.

### 14.1 Runs

All 500 steps from the checkpoint `z19z100mmS5_auto0` (step 17398), 100 zones, EDIAG=2, settings of z100mm otherwise. "Newton
failures" is the step control's count of steps repeated because a stage did not converge (`fail_newton` of the DTERR statistics).

| run | tree | rows | time covered [yr] | Newton failures |
|---|---|---|---|---|
| t_en24 | e0e6b2e | current stack (ECOMPAT off): the control with per-stage records; z100mm itself has none | 7.42 | 88 |
| t_en19 | e0e6b2e | ECOMPAT=1 (§5–§7) | 12.65 | 74 |
| t_en25 | 39923e1 | ECOMPAT=1 ECGRAV=1: G5 put back as a sink in the zones where the current stack uses its total row | 9.65 | 76 |
| t_en32 | bb57e3e | ECOMPAT=1 ECGRAV=2: G5 as a sink in every zone | 0.118 | 87 (+ 1 error rejection) |
| t_en26 | bb57e3e | ECOMPAT=1 ECGRAV=3: the exact mass-flux gravity term (G5x, §14.3) in every zone | 0.058 (failed at step 17845, 447 steps) | 146 |
| t_en27 | bb57e3e | ECOMPAT=1 GEQGP=1,100: mesh weight on density jump × gravity power (§14.6) | 0.156 | 61 (+ 6 error rejections) |
| t_en31 | bb57e3e | ECOMPAT=1 GEQGP=1,300 | 0.217 | 76 |
| t_en30 | bb57e3e | ECOMPAT=1 GEQGP=1,1000 | 11.38 | 71 |

### 14.2 An arbiter that uses neither row (item 1b)

`gledger.jl`: from the state after every converged stage (r, u, ρ, e, E, zone masses in ediagz.jls) K = Σ_k ¼ m_k (u_k² + u_k+1²),
U = Σ m_k e_k, E_r = Σ E_k ΔV_k and an explicit binding energy W, differenced over each complete TR-BDF2 step, plus the boundary
fluxes and the non-gravitational sources (lines, ray, nuclear) of the total row's own bookkeeping:

    defect = d(K + U + E_r + W)/dt + (flux out − in) + other sources.

No gravity power of either row enters. W = W_self + W_ext; W_ext is the potential energy in the external (dark-matter) mass
EXTMASS(r) the code adds to the enclosed mass (its rate is −1.2 L☉, two quadratures agree to 0.08 L☉). Four discrete W_self, to
measure the sensitivity asked for:

    (i)   interface-based        W = −Σ_j G M_j μ_j / r_j              (M_j = gas mass inside interface j, μ_j = ½(m_j−1 + m_j))
    (ii)  zone-centred           W = −Σ_k G (M_k + ½ m_k) m_k / r_c,k    (r_c,k = ½(r_k + r_k+1))
    (iii) zone mass × mean interface potential  W = Σ_k m_k ½(ψ_k + ψ_k+1), ψ_j = −G M_j/r_j   (identical to (i) by summation by parts)
    (iv)  the exact integral −∫ G M(r) dm / r for piecewise-constant density in each zone

A ledger of z100mm itself from its snapshots is not possible (one snapshot per 50 steps, no per-step boundary fluxes); t_en24 is its
rerun with per-stage records (§14.1). Time-weighted means over the 500 steps [L☉]:

| | t_en24, current stack | t_en19, ECOMPAT | t_en25, ECOMPAT + ECGRAV=1 |
|---|---|---|---|
| d(K + U + E_r)/dt from the states | +139617.6 | +102166 | +110147.2 |
| boundary flux out − in | −101865.1 | −69797 | −76220.1 |
| other sources | +665.8 | +654 | +651.3 |
| gravity term of the total row, Σ F_j g_j ½(r_j+1 − r_j−1) | −34348.0 | −30309 | −33673.1 |
| gravity power of the momentum row, Σ u_j g_j μ_j | −37512.7 | −33024 | −36767.4 |
| G5 = the difference of the two | 3164.7 | 2715 | 3094.4 |
| G5x: exact semi-discrete dW(i)/dt − momentum row's power (§14.3) | 2946.2 | 2610 | 2886.9 |
| dW_self/dt (i) / (ii) / (iv) | −34562.0 / −34612.7 / −34354.9 | −30409.0 / −30420.2 / −30136.7 | −33876.9 / −33912.6 / −33620.4 |
| Σ of the total row's residuals on the solution | 4070.2 | 2715 | 905.2 |
| **ledger defect with W (i) / (ii) / (iv)** | **+3855.2 / +3804.5 / +4062.2** | **+2613.5 / +2602.3 / +2885.8** | **+700.3 / +664.6 / +956.8** |

(The three runs cover different times, 7.42, 12.65 and 9.65 yr, so the columns are not at equal time; the comparison is of each
run's own terms.) Reading:

* **The binding energy follows the total row's flux form, not the momentum row's power.** dW/dt − (total row's gravity term) is
  −214 / −265 / −7 L☉ (t_en24), −100 / −111 / +172 L☉ (t_en19), −204 / −240 / +52 L☉ (t_en25) for definitions (i) / (ii) / (iv):
  0.02–0.8 % of the release rate. dW/dt − (momentum row's power) is +2.6e3 … +3.2e3 L☉ in all three runs and for every
  definition: 8–9 % of the release rate and 30–45 % of L_em.
* **Sensitivity to the definition of W**: (i) and (ii) differ by 11–51 L☉, (iv) from (i) by 207–283 L☉. The spread (≤ 0.3e3 L☉) is
  a tenth of G5; it does not change which form the binding energy follows, and it is the floor below which "closes" has no meaning
  on this mesh.
* **Whose ledger closes better.** ECOMPAT's defect, 2.6e3–2.9e3 L☉, *is* G5x: it is created in the star and the infall (per zone in
  §14.4). The current stack's defect is larger, 3.8e3–4.1e3 L☉, but of different origin (imbzone.jl, Σ of the total row's residuals
  per zone range): zones 1–14 0.0, zones 15–26 154, zones 27–50 620, zones 51–100 **3296 L☉**. The last is the radiation term of the
  comoving row in the tied zones (§3), which the lab-frame row removes (t_en19: closed there; t_en25: 3.2 L☉). So: in the outer
  cloud ECOMPAT's rows close and the stack's do not; in the star the stack's total row closes against W to the definition spread
  and ECOMPAT's row is open by G5x; in the infall (zones 27–50, which the stack also runs on the gas first law) both are open.
  Neither closes. t_en25 (both repairs, §14.5) is the best of the three, 0.7e3–1.0e3 L☉.

### 14.3 Derivation (item 1c): what conserves the binding energy on this mesh

Take W(i) = Σ_j μ_j ψ_j, ψ_j = −G M_j / r_j (per steradian the code has 4πG; the external mass adds its own ψ). Interfaces move with
ṙ_j; the mass flux through interface j relative to the mesh is F_j = r_j² ρᵇ_j urel_j (ρᵇ the advected, van Leer upwinded interface
density; urel_j = u_j − ṙ_j), so ṁ_k = F_k − F_k+1, Ṁ_j = −F_j and μ̇_j = ½(F_j−1 − F_j+1). Differentiating and summing by parts,

    dW/dt = Σ_j μ_j g_j ṙ_j + Σ_j F_j D_j,      g_j = G M_j / r_j²,
    D_j = ½(ψ_j+1 − ψ_j−1) + G μ_j / r_j = ½ G [ M_j+1 (1/r_j − 1/r_j+1) + M_j−1 (1/r_j−1 − 1/r_j) ].

The first sum is the mesh carrying each node with all its mass through the field; the second is the exchanged mass times the
potential drop D_j it falls through. The term G μ_j/r_j is the change of ψ_j itself when mass crosses interface j (Ṁ_j = −F_j); with
it D_j = g_j ½(r_j+1 − r_j−1) (1 + O(Δ²)), the total row's form. Checked numerically: gledger.jl evaluates this expression on every
record and compares it with the difference of W(i) between the states; the remainder is 0.09 L☉ in the median and 2.5 L☉ in the
time-weighted mean (t_en25), a time-discretisation term (W is not linear in the unknowns). This is the structure of Springel (2010,
§5.4, eqs. 92–96): the gravitational energy change of a cell is the work on the cell moving with all its mass, plus Σ over faces of
(exchanged mass) × (potential difference); his §5.2 states that the work term built from the cell-centred momentum "may not
precisely balance the actual energy extracted from the gravitational field, which is determined by the mass fluxes", that the
violation is of first order and matters "in collapse problems that involve strong shocks … when the spatial resolution is relatively
poor". On fixed Eulerian meshes Jiang et al. (2013) and Mullen et al. (2021) obtain round-off conservation of the total energy
with a gravitational energy flux; Jiang et al. find that the errors of the non-conservative source "can affect the dynamics of
the system", and Katz et al. (2016) report that standard coupling does "not adequately conserve the total energy" in white-dwarf
mergers (of these three only the abstracts were read; their formulas are not used here).

The momentum row's power is Σ_j u_j g_j μ_j = Σ_j μ_j g_j ṙ_j + Σ_j μ_j g_j urel_j. Hence

    −dW/dt − (−Σ u g μ)  ⇒  G5x = Σ_j urel_j [ r_j² ρᵇ_j D_j − μ_j g_j ],

and G5 of §2 (the total row's form, which uses the lab velocity u in both terms) is G5x plus the mesh-motion part and the O(Δ²)
difference between D_j and g_j ½(r_j+1 − r_j−1): measured 2715 against 2610 (t_en19), 3165 against 2946 (t_en24), 3094 against 2887
(t_en25). So the statement of §7.1 that G5 "is not a mesh-motion term" has to be read with this: nine tenths of it is the mass-flux
term G5x, the rest is what the total row itself gets wrong against W (it books the mesh-motion part with ρᵇ V̂ in place of μ).

*The coordinator's reading* — Σ u g μ + d/dt Σ μ ψ telescopes only if the force on interface j is r_j² ρᵇ_j ½(ψ_j+1 − ψ_j−1), with the
same advected density as the mass flux — is right for a fixed mesh with one correction: the force must be r_j² ρᵇ_j D_j, i.e. with the
self-gravity term G μ_j / r_j (without it the ledger is open by Σ F_j G μ_j/r_j, the potential change of the matter outside interface j).
On the moving mesh this is not enough: the part of the power that belongs to ṙ_j needs the weight μ_j g_j.

*Is there a smooth force?* For any discrete W(r, m) the two velocities ṙ_j and urel_j are independent, so a force f_j that depends on
the state only and makes f_j (ṙ_j + urel_j) equal to the rate of W for every motion must satisfy both

    f_j = ∂W/∂r_j     and     f_j = r_j² ρᵇ_j (∂W/∂m_j − ∂W/∂m_j−1)

(for W(i): μ_j g_j and r_j² ρᵇ_j D_j). Together they fix the interface density of the mass flux, ρᵇ_j = μ_j g_j / (r_j² D_j) ≈
μ_j / V̂_j (V̂_j = r_j² ½(r_j+1 − r_j−1)): one centred value, independent of the direction of the flow. With an upwinded ρᵇ — which
the advection of a density profile that falls by a factor 2–4 per zone needs — **no state-only force conserves K + U + W + E_r to
round-off on a mesh that both moves and lets mass through.** The two limits are exact: Lagrangian (F = 0), where the momentum row's
μ_j g_j is the adjoint of W (the isolated-star line, §7.3, is in this limit apart from GRIDEQ's motion); and the fixed mesh, with
the force r_j² ρᵇ_j D_j, which jumps at u_j = 0 when ρᵇ is upwinded and then gives a node at rest no defined weight (a hydrostatic
star could not be posed).

What does conserve to round-off is not a force but an energy term: keep the momentum row, and put G5x into the gas row (§14.5,
ECGRAV=3). Then K + U + W + E_r is conserved, and the truncation error of the gravity power is moved from the binding-energy
ledger into the gas's internal energy — which is what every total-energy code does, and what fails where the internal energy is a
small part of the total (§6.3, §14.5). So the PI's caveat holds in this precise form: on this mesh the first law for the gas and
the conservation of the binding energy cannot both be exact; G5x measures by how much, and the only general way to have both is to
make G5x small.

### 14.4 Where G5 sits and what it is made of

`g5split.jl` (t_en19, per interface, time-weighted): with V̂_j = r_j² ½(r_j+1 − r_j−1), V̄_j = ½(ΔV_j−1 + ΔV_j) and ρ̄_j = μ_j / V̄_j,

    ρᵇ_j V̂_j − μ_j = ρ̄_j (V̂_j − V̄_j)  +  V̂_j (ρᵇ_j − ρ̄_j)
                      geometric          density

The geometric part is there for uniform density: on a mesh of ratio q = r_j+1/r_j, V̄/V̂ = (q² + 1 + q⁻²)/3 (1.124 for q = 1.35;
≈ 1 + (4/3)(ln q)²). The density part is split further with ρ_c,j, the Δr-weighted mean of the two zone densities, into
(upwinded − centred) and (centred − node mean).

| interfaces | G5 | geometric | density | of it upwinded − centred | centred − node mean |
|---|---|---|---|---|---|
| 2–14 (interior) | 189.4 | 206.9 | −17.5 | | |
| 15–27 (outer layers of the star, photosphere) | 1635.7 | 189.9 | 1445.8 | 2126.4 | −680.6 |
| 28–51 (shock, infall) | 885.4 | 1461.1 | −575.7 | −21.4 | −554.3 |
| 52–100 (outer cloud) | 5.1 | 6.8 | −1.7 | | |
| all | 2715.6 | 1864.6 | 850.9 | 2307.5 | −1456.6 |

Largest interfaces (G5 / geometric / density [L☉]; Δln r; V̄/V̂; ln(ρ_j−1/ρ_j)): 20: 175.9 / 1.0 / 174.9; 0.130; 1.0003; 0.70 —
25: 317.5 / −4.3 / 321.8; 0.044; 0.995; 1.35 — 26: 222.6 / −0.1 / 222.6; 0.040; 1.000; 1.16 — 34: 123.4 / 232.1 / −108.7; 0.262;
1.163 — 35: 145.7 / 251.8 / −106.1; 0.317; 1.195 — 36: 129.9 / 216.2 / −86.2; 0.356; 1.210.

Two different error terms, then:

* **In the star's outer layers (interfaces 19–27)** the mesh is fine (Δln r 0.04–0.15, geometric part ≤ 17 L☉ per interface) and the
  density falls by a factor 1.8–3.8 per zone. G5 there is the first-order term u g V̂ (ρᵇ − ρ̄): the mass that actually crosses the
  interface is the upwinded (outer, thin) zone's, the weight the momentum row accelerates is the node's mean.
* **In the infall (interfaces 28–51)** the density is smooth and the zones are wide (Δln r 0.2–0.37). G5 there is the second-order
  quadrature term: the momentum row weighs the node with μ_j = ρ̄ V̄, and V̄ exceeds V̂ by 8–21 %. In steady infall on a fixed mesh
  the momentum row then gives the gas per zone the Bernoulli gain g Δr × μ_j/(ρᵇ_j V̂_j): too much kinetic energy by that factor,
  because the momentum advection uses the interface mass flux F and the weight uses the node mass.

**This second term is visible in the solution without any ledger** (`uff.jl`, the checkpoint state, which is z100mm's): a
pressureless shell that falls from rest keeps the mass inside it, so u² r / (2 G M_enc) < 1 (pressure, the radiation force and the
dark-matter mass left behind only lower it). The infall is cold (c_s 3.2–3.6 km/s against |u| 17–51 km/s). Measured at interfaces
30–38: u² r/(2GM) = 0.947, **1.111, 1.035, 1.048, 1.027, 1.003**, 0.961, 0.906, 0.838; the gain of ½u² between neighbouring
interfaces over the potential drop between them is 1.12, 1.09, 1.11, 1.09, 1.07, 1.03 at interfaces 33–38, and the factor
μ_j/(ρᵇ_j V̂_j) there is 1.07, 1.09, 1.11, 1.11, 1.11, 1.10. The gas arrives at the accretion shock faster than free fall, by up to
11 % in kinetic energy, in the production run and in every run of this report: the momentum row's gravity power in the infall is
too large, and what ECOMPAT thermalises at the shock includes it. That settles "which of the two is right" for the infall: the
flux form.

The momentum row's weight μ_j g_j is nevertheless not simply "wrong": it is exact where the mesh carries the mass (Lagrangian
limit: the pointwise acceleration of interface j is g_j, and μ_j g_j is the adjoint of W). Replacing it by the finite-volume weight
g_j ρ_c,j V̂_j — the exact weight of two half zones of constant density in a point-mass field, and the more accurate hydrostatic
weight on a wide mesh (uniform density, q = 1.35: Δp from the code's weight 15 % above the exact value in a point-mass field and
10 % above for self-gravity, from the finite-volume weight 2 % above and 2 % below) — would remove the infall's geometric part
(G5 of interfaces 28–51: 885 → −21 L☉), leave the star's larger (interfaces 15–27: 1636 → 2126 L☉) and make a Lagrangian shell
fall too slowly by V̂/V̄. Not implemented: it trades one limit for the other, as §14.3 says it must. *Round 3: implemented as the
provisional switch `GFV` and measured, together with a variational weight `GVAR` (§19). The estimate for the infall holds
(−18 L☉ measured); the star's centre changes with the weight, which this paragraph did not anticipate.*

### 14.5 The controlled runs (item 1a): G5 put back as a sink

Provisional switch `ECGRAV` (default 0, with ECOMPAT; `ecgrav` in TitanResidFast.jl): a gravity term in the compatible gas row, taken
from the gas as heat, everything else of ECOMPAT unchanged (lab-frame radiation row, KSTAR, no dual-energy switch).

    1: the total row's form (G5 of §2), in the zones in which the current stack uses its total row (useint = 0: the star)
    2: the same in every zone
    3: the mass-flux form G5x of §14.3, in every zone: K + U + W(i) + E_r conserved up to the time discretisation
    4: G5x in the stack's total-row zones          (not run)

**ECGRAV=1 (t_en25)** is the run asked for: it differs from ECOMPAT only by G5 in the star, and from the current stack by everything
else. It ran the 500 steps (9.65 yr; Newton failures 76 against 74 for ECOMPAT and 88 for the stack). The term it removes is 2189 L☉
(time-weighted; zones 1–14 327, zones 15–26 1862) of a G5 of 3094 L☉; the infall's 863 L☉ stays, as it does in the stack. The star at
equal time (`starhist.jl`; photosphere = the interface of zone 27; r_ph [R☉] / T_eff [K] / L at the photosphere / L at 1.5 r_ph [L☉]):

| t − t_ck [yr] | current stack (t_en24) | ECOMPAT + ECGRAV=1 (t_en25) | ECOMPAT (t_en19) | part of the ECOMPAT − stack difference removed: r_ph, T_eff |
|---|---|---|---|---|
| 2.0 | 106.79 / 4935 / 6097 / 6101 | 109.06 / 4871 / 6031 / 6052 | 113.23 / 4784 / 6052 / 6194 | 65 %, 58 % |
| 3.0 | 105.24 / 4992 / 6199 / 6219 | 108.03 / 4913 / 6125 / 6206 | 112.87 / 4818 / 6183 / 6331 | 63 %, 55 % |
| 4.0 | 105.82 / 5005 / 6329 / 6511 | 107.91 / 4943 / 6263 / 6587 | 113.31 / 4832 / 6307 / 6727 | 72 %, 64 % |
| 5.0 | 105.95 / 5006 / 6351 / 6876 | 108.45 / 4964 / 6433 / 6992 | 113.82 / 4852 / 6467 / 7122 | 68 %, 73 % |
| 6.0 | 106.45 / 4999 / 6374 / 7039 | 109.02 / 4965 / 6507 / 7263 | 114.68 / 4860 / 6611 / 7442 | 69 %, 76 % |
| 7.0 | 107.80 / 4990 / 6493 / 7285 | 110.77 / 4940 / 6583 / 7329 | 116.62 / 4847 / 6764 / 7652 | 66 %, 65 % |
| 7.4 | 108.13 / 4992 / 6543 / 7374 | 112.32 / 4936 / 6747 / 7386 | 117.44 / 4840 / 6817 / 7721 | 55 %, 63 % |

T_c at 7.0 yr: 1.2477e5 / 1.2466e5 / 1.2429e5 K. So: **the gravity term in the star's zones is the larger part of the difference
between ECOMPAT's star and the stack's (55–72 % of the radius difference, 55–76 % of the T_eff difference), not all of it.** 2–4 R☉
and 35–60 K remain; they belong to what else differs between t_en25 and the stack (the radiation row's frame in the tied zones, the
gas row in the zones the stack runs on the gas first law with η-switching, KSTAR) and were not separated further. "By exclusion"
(§7.2) becomes: attributed by a controlled run, to two thirds.

**ECGRAV=3 (t_en26)**, the exact term in every zone, **fails.** The ledger does what the derivation says (defect against W(i) 1.4 L☉
in the median, 43 L☉ in the time-weighted mean; against W(ii) 1.7e3, W(iv) 9.6e3 L☉ — conserving one discrete W exactly says
nothing about the others). But the run dies at step 17845 ("failure to converge after 60 attempts") after 447 steps and 0.058 yr,
with 146 Newton failures: the infall's part of the term, 100–170 L☉ per zone at the start, is a heat sink in gas whose whole
internal-energy flux is of that order. Zone 32 cools from 1908 K to 12.8 K (zones 29–35 end at 666, 249, 26.7, 12.8, 270, 651,
1185 K, from 2073–1770 K), the density there responds, and the term grows to 234–851 L☉ per zone (zones 32–37). This is the failure
§6.3 described for the total row in the cold supersonic infall, reproduced with the exact term: there the gravity power's
truncation error is larger than the gas's internal energy flux, and whichever row it is put in, the gas temperature is its residue.
The kinetic energy is where the error is (§14.4: the infall is too fast), and no gas energy row can take it out again.

**ECGRAV=2 (t_en32)**, the term of ECGRAV=1 in every zone instead of in the flagged ones — the variant without a zone flag —
**fails in the same way**, more slowly. It completes the 500 steps but they cover 0.118 yr (t_en25: 9.65 yr), with 87 Newton
failures and one error rejection. The infall's part of G5 (660 L☉ in zones 27–50 over this run) is taken from gas that does not
have it: at 0.118 yr zones 27–33 are at 217, 87, 39, 75, 180, 327, 527 K (t_en19 at the same time: 4622, 3792, 2104, 1987, 1942,
1917, 1872 K), the layer has spread outwards (zone 27 at 134 against 112 R☉), and the centre has changed with it (T_c 1.182e5
against 1.215e5 K). So what makes ECGRAV=1 run is the dual-energy flag of the current stack (`dualswitch`: the term only where
the internal energy is not small against the kinetic energy), and with it the dependence of the star on which zones are flagged
(§2, §7) comes back. A gravity sink in the gas row is a row switch or it fails; there is no third variant.

### 14.6 What reduces G5 at 100 zones; the mesh weight tested (GEQGP) — negative

G5 is the error indicator of the gravity power, known per interface (§14.4). Its two parts scale differently:

* the **geometric part** (infall, interfaces 28–51: 1.46e3 L☉) is u g ρ̄ (V̂ − V̄) with V̄/V̂ − 1 = (q − 1/q)²/3 ≈ (4/3)(Δln r)²:
  second order in the zone width, 8–21 % at the infall's Δln r = 0.20–0.37;
* the **density part** (outer layers of the star, interfaces 15–27: 1.4e3 L☉, nearly all of it at 19–27) is u g V̂ (ρᵇ − ρ̄):
  first order in the density jump per zone, which is ln(ρ_j−1/ρ_j) = 0.57–1.35 at interfaces 19–27, at Δln r = 0.04–0.15.

*Round 3 (§21.2):* "first order in the density jump per zone" is a statement about one interface, not about the zone number.
The advected density ρᵇ is a second-order monotone interpolation, so where the layers are resolved the density part falls as
the square of the zone width: measured from the IC at 50, 100 and 200 zones, the density part of the outer layers falls with
order 1.4–1.7 (50 → 100) and 1.9–2.0 (100 → 200) in the zone number, and the geometric part of the infall with order 1.4–1.6,
less than the (Δln r)² of the smooth flow.

The mesh has 100 zones for ln(r_out/r_in) = ln(19 pc / 0.1 R☉) = 22.9, a mean Δln r of 0.23; the photosphere and the shock take
Δln r = 0.04–0.07 (interfaces 23–30) and are paid for with 0.3–0.67 in the interior, 0.2–0.37 in the infall and 0.2 on average
in the cloud. A geometric error of 1 % needs Δln r = 0.09 through the infall, and halving the density part needs the outer
layers' nine zones doubled. Both are zones taken from somewhere else: at 100 zones the gravity power's error can be moved, not removed, by the mesh.

**The weight asked for, tested** (`GEQGP=w_g,L_ref`, provisional, default off; `grideq_rows!`, in the implicit grid equation, of
the Newton iterate, nothing lagged): the monitor gets the term

    w_g (∂ln ρ/∂j)² × ½(|u_j| g_j + |u_j+1| g_j+1) m_j / L_ref

— the density jump across the zone times the gravity power of the zone, in units of L_ref (the gravity power per zone is
10²–10³ L☉ in the star's outer layers and the infall). Three strengths with ECOMPAT, 500 steps from the checkpoint, against t_en19
over the same time (`g5split.jl` on t_en19's records up to that time; `starhist.jl`):

| run | L_ref [L☉] | time covered [yr] | Newton failures | G5, run / t_en19 at equal time [L☉] | of it interfaces 15–27 | T_c against t_en19 |
|---|---|---|---|---|---|---|
| t_en27 | 100 | 0.156 | 61 (+ 6 error rejections) | 9926 / 1659 | 4966 / 1098 | −4.7 % |
| t_en31 | 300 | 0.217 | 76 | 4002 / 1842 | 2278 / 1234 | −2.2 % |
| t_en30 | 1000 | 11.38 | 71 | 2653 / 2706 | 1564 / 1625 | −0.8 % |

* **Strong and intermediate weight (L_ref = 100, 300): G5 is 6 and 2.2 times larger, not smaller.** The weight pulls zones into
  the star (t_en27: interfaces 15–27 at 16.6–94 R☉ instead of 27.9–111 R☉); for the whole 0.16–0.22 yr the mesh is moving through
  the star, the steps are 60–80 times shorter than in t_en19, and the star does not stay the same: T_c falls by 4.7 % and 2.2 %,
  ρ_c by 9 % and 4 %, the interior acquires velocities, and the interior interfaces 7–22, which carry at most 190 L☉ each in t_en19,
  carry 94–631 L☉ (t_en27; the excess is geometric, 9170 of 9926 L☉, although V̄/V̂ there is no worse than 1.01–1.14: G5 is
  proportional to the velocity). Neither run reaches a settled mesh in 500 steps, so they do not say what G5 would be on a mesh
  relaxed with this weight; they say that getting there from a production state disturbs the star by more than the error the
  weight was to remove.
* **Weak weight (L_ref = 1000): the mesh relaxes, the run is as robust as without it (11.38 yr, 71 failures), and G5 is 2 %
  smaller.** Interfaces 17–22 are 2–5 % closer (Δln r 0.202–0.086 against 0.211–0.088), the density jumps there 2–4 % smaller
  (0.550–0.747 against 0.567–0.763), G5 of interfaces 15–27 3.7 % smaller (1564 against 1625 L☉), the infall's unchanged (878
  against 874 L☉): the first-order scaling. The ledger defect against W is 2559 / 2546 / 2825 L☉ (definitions i / ii / iv;
  t_en19: 2613 / 2602 / 2886 over its 12.65 yr). The star's centre is 0.8 % cooler than in t_en19 at 11.38 yr (T_c 1.2428e5
  against 1.2528e5 K, ρ_c 1.6 % lower, r_ph 120.3 against 121.4 R☉), from the first 0.2 yr on.

So there is no strength at which this weight reduces G5 by an amount that matters: where it leaves the star alone it gains 2 %,
and where it acts it makes the indicator larger. A weight on the geometric part ((Δln r)² × gravity power, for the infall) was not
tested; by the zone count above it has the same budget to work with. **At 100 zones G5 remains what it is: an error of the gravity
power of 8–9 % of the release rate (30–45 % of L_em), known per interface, that the mesh does not remove.** What removes the
infall's half of it without zones is the quadrature of the weight (finite-volume weight, §14.4: 885 → −21 L☉), at the price of the
Lagrangian limit; what removes the star's half is resolution of the density gradient at interfaces 19–27. The convergence of
either star with zone number was not measured (every run here has 100 zones). *Round 3: measured at 50, 100 and 200 zones from
the IC, §21.*

### 14.7 Summary, and what it means for the default

What was measured and derived:

1. **The binding energy follows the flux form** (mass flux × potential drop), to 0.02–0.8 % of the release rate and for every
   discrete definition of W tried (spread 0.3e3 L☉). The momentum row's power Σ u g μ exceeds −dW/dt by 2.6e3–3.2e3 L☉, 8–9 % of
   the release rate, in the stack's run, in ECOMPAT's and in the controlled run alike (§14.2).
2. **No force that depends on the state alone conserves K + U + W + E_r to round-off on a mesh that both moves and passes mass
   with an upwinded interface density.** The coordinator's force is the exact one for a fixed mesh (with the self-gravity term
   G μ_j/r_j added) and discontinuous at u = 0; the momentum row's is the exact one for a Lagrangian mesh. Round-off conservation
   is possible only with an energy term (§14.3), and that fails in the cold infall, with the exact term (ECGRAV=3, t_en26) and with
   the total row's (ECGRAV=2, t_en32), as "total rows everywhere" did (§6.3, §14.5).
3. **Where the error sits** (§14.4): 1.6e3–1.9e3 L☉ in the outer layers of the star (interfaces 15–27, nearly all at 19–27;
   first order, upwinded density against the node's mean) and 0.7e3–0.9e3 L☉ in the infall (interfaces 28–51; second order,
   node volume against the flux form's). The infall's part is in the solution of every run, the production run included: the gas
   reaches the shock with up to 11 % more kinetic energy than free fall from rest allows.
4. **Attribution** (§14.5): putting G5 back as a sink in the star's zones removes 55–72 % of the radius difference and 55–76 % of
   the T_eff difference between ECOMPAT's star and the stack's; 2–4 R☉ and 35–60 K belong to the other differences and were not
   separated.
5. **The mesh does not remove it at 100 zones** (§14.6): the weight on density jump × gravity power gains 2 % where it leaves the
   star alone and makes G5 2–6 times larger where it acts.

Per region, the energy that K + U + E_r + W gains or loses without cause, time-weighted [L☉] (`imbzone.jl`, the conservative
total row's residual on each run's solution; the row follows W to the definition spread):

| zones | current stack (t_en24) | ECOMPAT (t_en19) | ECOMPAT + ECGRAV=1 (t_en25) |
|---|---|---|---|
| 1–26 (star) | 154 | 1972 | 152 |
| 27–50 (shock, infall) | 620 | 740 | 750 |
| 51–100 (outer cloud) | 3296 | 3 | 3 |
| all; ledger against W (i) / (ii) / (iv) | 4070; 3855 / 3805 / 4062 | 2715; 2614 / 2602 / 2886 | 905; 700 / 665 / 957 |

For the decision:

* **ECOMPAT as it stands is not the more conserving formulation for the star.** It closes internal + kinetic + radiation energy
  to round-off against the momentum row's forces, removes the 3.3e3 L☉ the comoving radiation row creates in the outer cloud, and
  needs no flag. But in zones 1–26 it gives the gas 2.0e3 L☉ (about 30 % of L_em) that the binding energy does not lose, where
  the current stack's total row is open by 0.15e3 L☉. Its star is larger by 6–9 R☉ and cooler by 140–175 K at 2–7 yr after the
  checkpoint, and from the IC less compact (T_c lower by 27–29 % at equal mass, §16); from the checkpoint two thirds of the
  difference is this term (the from-IC difference was not separated in the same way; *round 3: it is now, §20 — ECGRAV=1 from
  the IC restores the control's T_c to +1.4…+3 % and its radius to +3…+4 %*). §7 read the difference as heat the stack
  withholds; against W it is heat ECOMPAT adds. On this point the measurements speak against making ECOMPAT, as it stands, the
  default.
* **The combination with the best ledger is a dual formulation again.** ECOMPAT + ECGRAV=1 (compatible gas row, lab-frame
  radiation row, the gravity sink in the zones the dual-energy flag selects) has a fifth to a third of either defect and the
  stack's star to within 2–4 R☉, at the same robustness (9.65 yr, 76 failures; stack 7.42 yr, 88). It keeps `dualswitch` and its
  three parameters, and with them the dependence on which zones are flagged; of §11's list only the total-row *solve* would go.
  Without the flag the sink fails (t_en32).
* **The lab-frame radiation row and the frame-consistent closure are independent of this**: the outer
  cloud's 3.3e3 L☉ goes with them in both runs that have them. The current stack's gas rows with the lab-frame radiation row
  alone were not run (in this branch the two are switched together by `ECOMPAT`); whether the stack's total row stays closed with
  it is therefore not measured. *Round 3: run with its own switch `ERLAB` (§18). It gives the stack's star with the outer
  cloud's 3.3e3 L☉ removed, and a ledger of 479–741 L☉.*
* **The infall's error is the momentum row's, in every variant**, and no energy row repairs it (§14.5). It needs the
  finite-volume weight of §14.4 where the flow passes the mesh and the node weight where the mesh follows the flow; a weight that
  depends on urel/u is not a state-only force and was not built or tested (*round 3: the two weights were built and run, and
  §19.2 shows why no continuous form has both limits*). This is the precise form of the PI's caveat: a single
  formulation conserves K + U + E_r to round-off (shown) or tracks W to truncation in the star (the stack's total row does), not
  both, at 100 zones.
* **Neither star is the converged one.** The two gravity powers differ by 8–9 % of the release rate; the star's radius and T_eff
  move by 6–9 % and 3 % with the choice. Which is closer to the continuum can be read off the ledger against W (the stack's,
  for the star) but not confirmed without a run at higher zone number, which this task excludes. *Round 3: the zone-number
  runs were made (§21), and the ledger from the IC over 4 yr (§20.2) supports the reading for the star. One statement of this
  section needs a qualifier: "the total row's gravity term follows W" holds at the checkpoint to 0.02–0.8 % of the release rate;
  from the IC, dW/dt minus that term scatters by ±0.6e3 L☉ between windows of 0.7–1 yr and is −0.15…−0.84e3 L☉ over 4 yr
  (1–7 % of the release rate), depending on the definition of W.*

## 15. Radiation energy density below the bath (follow-up item 2; open point 10)

§10 item 10 reported E below a T_CMB⁴ by up to 27 % in zones 64–72 of the ECOMPAT run from the IC and suspected the static formal
solution. The measurements below change the finding: **the defect is not ECOMPAT's. It is the flux tie's, it is in the control and
in the production run, and there it is deeper.**

### 15.1 What "segment FAILED (z19t_en_ic1S1)" at 01:55 was

My own stop of t_en_ic1 at step 1478 (§8: "stopped by hand … to free the cores"): the run wrapper writes that line when its Julia
process ends without the end-of-run line. The run had not failed (no Newton failure sequence, no convergence abort in seg.out).

### 15.2 Where and when (`ebmin.jl`, `icdef.jl`, `outerL.jl`; E_b = a T_CMB⁴, T_CMB = 54.5 K)

min over zones of E/E_b − 1, and the zone:

| run | rows | min E/E_b − 1 (step, zone) |
|---|---|---|
| t_en_ic1, from the IC | ECOMPAT | −0.08 (150, 62); −0.49 (350, 63); −0.45 (600); −0.39 (1000); −0.28 (1200); −0.27 (1478, 67) |
| t_en_ic0, from the IC | current stack, same tree | −0.07 (250, 59); −0.36 (400, 58); −0.47 (600, 61); −0.51 (700, 63); −0.44 (1050); −0.34 (1450); −0.27 (1900); −0.07 (2100) |
| z100mm, production | current stack | −0.30 (550, 61); −0.62 (1050); −0.87 (1550); −0.97 (2050); −0.998 … −0.9998 (3050 … 13550, zones 64–72); −0.997 (17400, the checkpoint, 73); −0.89 (18050); −0.80 (32050, zone 74) |
| t_en19, from the checkpoint | ECOMPAT | −0.998 at the start (zone 73); −0.95 (5.1 yr); −0.80 (10.6 yr); −0.69 (12.65 yr, zone 74) |
| t_en24, from the checkpoint | current stack | −0.998 at the start; −0.90 (7.42 yr, zone 74) |

In z100mm 6–26 zones are below the bath by more than 1e-3 at every snapshot since step 550. At the checkpoint (r in pc; face k is
the inner face of zone k):

    zones 44–67 (1.6e-4 … 0.027 pc)   L_F = 6.03e3 … 6.25e3 L☉ (the star's), E/E_b − 1 = 2.8e3 … 0.08
    zones 68–75 (0.033 … 0.135 pc)     E/E_b − 1 = −0.02, −0.32, −0.73, −0.95, −0.995, −0.998, −0.93, −0.10;
                                       L_F falls from 6.3e3 (face 68) to 16 L☉ (face 74)
    zones 76–87 (0.16 … 1.07 pc)       E/E_b − 1 = +2.8, +5.4, +4.9, +3.1, +1.8, +1.0, 0.53, 0.34, 0.40, 0.56, 0.70, 0.71;
                                       L_F = 2.4e4, 4.1e4, 3.9e4, 2.7e4 L☉ at faces 77–80: four to seven times the star's

So the radiation field of the production run between 0.03 and 0.14 pc is 0.2–27 % of the CMB's, with a pulse of six times the
CMB's energy density outside it that carries several times the star's luminosity, 18 yr after the core formed (light takes 0.4 yr
to 0.11 pc). Eleven years later (step 32050) the minimum is still −0.80 at zone 74 and faces 78–79 carry 4.8e4 L☉ against 7.9e3 L☉
at face 45. In the earlier snapshots it is worse: at step 13550 the star's 2.0e4 L☉ are in the flux out to face 45 (2e-4 pc) only,
faces 52–73 carry 13–40 L☉, and E is 0.03–43 % of the bath in zones 51–73 (that part of z100mm was made with an older tree and
other closure settings, §16; its history was not analysed).

The control from the IC (t_en_ic0) has, besides the trough, a spurious pulse that ECOMPAT does not have: E/E_b − 1 = +19, +75,
+132, +141 in zones 64–67 at step 600 with L_F = 1.6e4–1.1e5 L☉ there (the star: 350–400 L☉), L_F at face 70 up to 2.8e5 L☉
(step 1600) against L* ≈ 2e3. That is the radiation energy the comoving row creates in the tied zones (§3); with the lab-frame
row the flux at faces 62–80 follows the star's luminosity (t_en_ic1, step 1478: 2515, 2497, 2199, 1580, 1640, 1986, 2014 L☉ at
faces 40, 56, 62, 66, 70, 74, 80). **ECOMPAT's field is the mildest of the three; the trough itself is in all of them.**

(The value 1.2e5 L☉ at face 100 that `ebmin.jl` prints for the runs from the IC is F = 5e-6 c E at 19 pc, with E within 3e-6 of the
bath: the bath's tolerance level times a large area, present from step 50, not a luminosity.)

### 15.3 Cause

In the tied zones (FLUXTIE with SLABTIE, §3: outside the photosphere, "flux tie from zone 32") the flux is not an unknown of its
own moment equation; the rows close it as F = h c E with h = H/J of the **static** formal solution for the present sources (`eddfac!`,
TitanCore.jl). The radiation energy row then is

    ∂E/∂t + r⁻² ∂(r² h c E)/∂r = sources,

a first-order advection equation with the single characteristic speed h c. Two consequences:

1. The flux is L(r) = L_fs(r) · E/E_fs: the formal solution's luminosity scaled by the ratio of the moment E to the formal J. It is
   right only where E = E_fs. When L* rises faster than the field can follow, E is still the old field (≈ E_b outside) while h is
   already that of the new source. Then E/E_fs = E_b/(E_b + L*/(4π r² c)) < 1 and *increasing outward*: more flux leaves each zone
   than enters, and E falls below E_b — the response to a brightening star has the wrong sign there.
2. The mismatch E − E_fs travels at h c, not at c. Where the bath dominates, h = L*/(4π r² c E) ≪ 1 (F/cE = 0.012 at face 73 and
   0.008 at face 74 of the checkpoint, falling as r⁻²): the trough and the pulse that holds the drained energy crawl outward ever
   more slowly. At zone 74 r/c = 0.4 yr and r/(h c) is of order 40 yr. In the static closure of the *second* moment (the variable Eddington factor K/J from a static
   formal solution: Stone, Mihalas & Norman 1992; Davis, Stone & Jiang 2012; Jiang, Stone & Davis 2012 — abstracts read) the
   flux keeps its own moment equation and with it the light travel; a static closure of the *first* moment hands the signal speed
   to the closure.

Tested directly (`icmodel.jl`, t_en_ic1): the energy row alone, zones 58–100, on the recorded mesh with the recorded h and the
recorded inner boundary, no gas terms. Model A (the closure, F = h c E) gives min E/E_b − 1 = −0.096, −0.495, −0.514, −0.361,
−0.313, −0.277, −0.304 at steps 150, 350, 600, 1000, 1200, 1300, 1478; the run has −0.082, −0.485, −0.445, −0.385, −0.282, −0.220,
−0.266. Model B (the same with the excess over the bath carried outward at c) never falls below E_b. The closure alone
reproduces the defect; no gas term, no velocity term and neither energy row is needed for it.

### 15.4 The candidate fix, tested: the tie on the excess (TIEX) — not general

Provisional switch `TIEX=1` (default off; tree cd9922c): in the tied zones F = h_x c (E − E_b), h_x = H_x/J_x the ratio of the
formal solution's moments *of the excess over the isotropic bath* (`HXEDD` from `eddfac!`). An isotropic bath then carries no flux
whatever E_fs is, and the star's light travels at h_x c ≈ c.

* **From the IC (t_en_ic2, `ECOMPAT=1 TIEX=1`): fails.** Before the core forms there is no excess in the formal solution (J_x is
  4e-6–9e-6 of the bath, the tolerance level), h_x is the ratio of two noise terms — at the clamp |h_x| = 1 in 12–30 zones, with
  9–20 sign changes between neighbours — and it carries the tolerance-level excess of the moment E inward until it piles up
  (E/E_b − 1 = 0.09, 0.10, 0.12, 0.17 at zone 43 at steps 20, 40, 80, 100). Step 200 is reached at t = 10.9 yr (t_en_ic1:
  2.4e5 yr at step 100), with 64 + 24 Newton failures and Δt = 19 s. Stopped by hand at step 238.
* **From the checkpoint, where the star dominates (t_en29, `ECOMPAT=1 TIEX=1`, tree bb57e3e, 500 steps): it does what the closure
  should do, with one unphysical transient.** h_x = 1.000 in the tied zones from zone 68 outwards (a point source seen from afar),
  no sign changes. The trough fills on the light-travel time: min E/E_b − 1 = −0.996 at the start, −0.90, −0.53, −0.29, −0.087
  after 0.02, 0.12, 0.22, 0.43 yr (zone 74: r/c = 0.37 yr), and from 3.1 yr on E ≥ E_b in every zone (minimum +7e-4, +3e-4, +1e-4
  at 3.1, 5.5, 12.6 yr; t_en19 without the switch: still −0.69 at 12.65 yr). At 12.6 yr L_F = 9131 … 8385 L☉ at faces 62 … 87
  (0.01–1.07 pc: the star's luminosity with its light-travel delay), E/E_b − 1 falls as r⁻² from 1.78 to 9.6e-5 and F/(cE) from
  0.54 to 1.2e-4. The pulse of zones 76–87 and the stored excess of zones 88–99 (§3.3, §10 item 7: "does not relax within the
  tests") leave at c: 5.9e8 L☉ through face 100 and 3.0e8 L☉ through the boundary at 12.6 yr, with E/E_b − 1 = 3.5e-4 left in zone
  88 and 0.5–9 % in zones 89–100 (from 3–63 % in zones 88–99; r/c = 74 yr at the boundary). The run is as robust as without the switch (500 steps
  = 12.64 against 12.65 yr; Newton failures 73 against 74; limiting variables x_H2 and x_H+ as before), the energy identity closes as
  before (`cres` 0.000 in zones 1–87), the star is the same (12.6 yr: r_ph 121.81 against 121.73 R☉, T_eff 4831 against 4831 K,
  T_c 1.2526e5 against 1.2527e5 K), and the gas of zones 58–92 differs by |ΔT/T| ≤ 5e-4 at equal time.
  The transient is unphysical: while E < E_b the closure makes the deficit an inward beam, F = c (E − E_b) < 0 with |F|/(cE) up
  to 7.7 (faces 72–75 at 0.02 yr, L_F = −2.3e5 … −8.3e5 L☉). |F| ≤ cE is violated until the trough is filled.

A closure for the excess is right only where the moment excess *is* the formal solution's excess. Where it is not — initial
radiation, a transient, a deficit, tolerance-level noise — the static formal solution has no angular information about it, and the
ratio of moments of a signed quantity is unbounded. So TIEX is the right limit for a dominating star (t_en29 is what the field
of §15.2 should look like) and wrong before there is one (t_en_ic2): not a general fix. It stays a provisional switch for the
record (§11); a general closure has to reduce to it where the star dominates and to "no flux" for a static bath (§15.5).

### 15.5 What a proper fix needs

The closure's angular information has to belong to the same radiation field as the moments. Taking it from the *static* transfer
equation is the standard construction of variable-Eddington-factor schemes (Stone, Mihalas & Norman 1992: "computed from a solution
of the static transfer equation"; Davis, Stone & Jiang 2012: "time-independent"; Jiang, Stone & Davis 2012; abstracts only read):
there it closes the pressure, and the flux keeps its own moment equation (way 3 below). What those schemes do where the formal
solution is ahead of the field by the light-travel time was not looked up. Three general ways, none built:

1. **A retarded formal solution.** Keep the term (1/c) ∂I/∂t in the ray sweep of `eddfac!`. Derived here, not taken from a
   source: backward differencing of (1/c) ∂I/∂t + dI/ds = χ (S − I) along a ray gives dI/ds = χ' (S' − I) with χ' = χ + 1/(cΔt)
   and S' = (χ S + I_old/(cΔt))/χ' (limits: Δt → ∞ the static equation; χ → 0 free streaming of I_old over cΔt, first order in
   time). The formal solution then is the same field the moments describe, and h is its
   H/J at any Δt (static limit for cΔt ≫ r). It needs the old intensities on every ray — on the moving mesh an interpolation in
   impact parameter and radius — and a check that the formal moments and the moment unknowns stay consistent to the tolerance.
   This is the general fix; it is a change of the transfer solver, not of the rows.
2. **Close on the mismatch.** Write I = I_fs + δI. δI obeys the transfer equation with the source −(1/c) ∂I_fs/∂t: it has the
   angular distribution of the *change* of the formal field, not of the field. A closure F = c H_fs + c h_δ (E − E_fs) with
   h_δ = Ḣ_fs/J̇_fs carries a static bath with zero flux and the star's transient at its own streaming speed (model B is its
   limit h_δ = 1). It has the same weakness as TIEX where J̇_fs vanishes and would need a bounded, sign-definite h_δ there. Not
   derived beyond this sketch and not tested in the code.
3. **No tie**: the flux's own moment equation with the static Eddington factor in the zones concerned, which is the scheme the
   tie replaced (§3 has the reasons it was introduced; they were not re-examined here).

### 15.6 Consequence for adoption

The sub-bath E is not a reason against ECOMPAT (its lab-frame row removes the part of the defect that creates energy, §15.2) and
not something ECOMPAT repairs. It is an open defect of the production stack, older than this work. What it does to the gas and
the star, measured over 12.6 yr from the checkpoint by removing it (t_en29 against t_en19, §15.4): nothing above 5e-4 in T in zones
58–92 and nothing in the star. The zones concerned are optically thin and warmer than the bath, and the light that drives H⁻
detachment and the Lyman–Werner rates is taken from the ray solution, not from E. The effect over the 2.4e5 yr of the collapse and
of z100mm's earlier, deeper state (step 13550, §15.2) was not measured.

## 16. The control from the initial condition (follow-up item 3)

Run t_en_ic0: tree e0e6b2e, 100 zones, the same run.env and phaseB.env as t_en_ic1, without `ECOMPAT=1` (`EDIAG=2` kept). Stopped
by hand at step 2100 (M* = 0.34 M☉, 5.0 yr after the core forms), the stage at which t_en_ic1 was stopped; it had not failed.

| | t_en_ic0: current stack, tree e0e6b2e | t_en_ic1: ECOMPAT, tree e0e6b2e | z100mm (production; own older tree for its first 13.5k steps) |
|---|---|---|---|
| first perf.csv row with a core: step, t [yr] | 300, 2.4276162e5 | 200, 2.4106495e5 | 200, 2.378241e5 |
| core formation time against z100mm | +2.08 % | +1.36 % | |
| M* ≈ 0.17–0.18: step / M* / r* [R☉] / T_c [K] / log n_c | 1000 / 0.176 / 28.2 / 1.092e5 / 23.17 | 700 / 0.171 / 32.5 / 7.99e4 / 22.90 | 1400 / 0.17 / 30.7 / 1.090e5 / 23.16 |
| M* ≈ 0.21 | 1300 / 0.211 / 33.3 / 1.113e5 / 23.18 | 800 / 0.21 / 38.3 / 8.12e4 / 22.91 | 2100 / 0.21 / 34.7 / 1.132e5 / 23.18 |
| M* ≈ 0.25–0.26 | 1500 / 0.256 / 38.7 / 1.142e5 / 23.19 | 1000 / 0.246 / 43.8 / 8.20e4 / 22.92 | 3900 / 0.26 / 45.8 / 1.160e5 / 23.19 |
| M* ≈ 0.31–0.34 | 1700 / 0.313 / 45.0 / 1.164e5 / 23.20; 2100 / 0.340 / 46.8 / 1.171e5 / 23.21 | 1400 / 0.331 / 50.9 / 8.44e4 / 22.93 | 5300 / 0.321 / 54.7 / 1.190e5 / 23.20 |
| steps from the core to that mass | 1800 | 1200 | 5100–6300 |
| Newton failures TR + BDF at the last row; after the core | 147 + 21; 115 | 108 + 12; 95 | 295–335; 258–298 |

What is attributed, and to what:

* **The central temperature and density at equal mass are the rows'.** The control on the same tree reproduces z100mm (T_c within
  0.2–2 %, log n_c within 0.01) and ECOMPAT differs from both: T_c lower by 27–29 % (0.80–0.84e5 against 1.09–1.17e5 K), n_c lower
  by 0.26–0.27 dex. Which term of the rows does it was not separated from the IC (no ECGRAV run from the IC was made). From the
  checkpoint, two thirds of the difference in radius and T_eff is the gravity term G5 of the star's zones (§14.5), which the
  stack's total row does not give to the gas and ECOMPAT does; against the binding energy it is the stack's row that is closed
  in the star (§14.7). If the same term is behind the from-IC difference, the lower T_c and n_c are ECOMPAT's error; that is
  the likely reading, not a measured one (from the checkpoint T_c differs by 0.4 % only, at 7 yr). *Round 3: measured. ECOMPAT
  with the G5 sink (ECGRAV=1) from the IC has the control's T_c to +1.4…+3.0 % and n_c to −0.05…−0.06 dex at equal mass (§20.1),
  so the lower T_c and n_c of ECOMPAT are this term.*
* **The later core formation is not the rows'.** The control forms the core 2.08 % later than z100mm, ECOMPAT 1.36 % later: the
  tree (and the settings z100mm acquired after its first 13.5k steps) delays it by 2 %, ECOMPAT *advances* it by 0.70 % against
  its control. The perf.csv rows are 100 steps apart and the row before the core is 1.3e3–4.4e3 yr earlier, but the times compared
  are those of the first row with a core, at which M* is 0.004–0.015 M☉ in all three runs and the next row (0.07–0.09 M☉) is
  0.6–0.7 yr later: the core's age at that row is well below a year, 4e-6 of the times compared.
* **The step count of §8 is mostly not the rows'.** After the core the control needs 1800 steps to 0.34 M☉, ECOMPAT 1200 (a third
  fewer, with 95 against 115 Newton failures), z100mm 5100–6300. "Four times fewer steps" in §8 compared ECOMPAT with a production
  run whose first 13.5k steps were made on commit 1782928 without `SLABFLUX=1 SLABTIE=1 SUBZONE=0`; corrected there.
* **The radius is not settled by this control.** At 0.31–0.34 M☉: control 45.0–46.8, ECOMPAT 50.9, z100mm 54.7 R☉. The control
  differs from z100mm (same rows, different tree and early settings) by more than ECOMPAT differs from either, so the radius at
  this stage is sensitive to something other than the rows. Not analysed. *Round 3: analysed in §20.3 — z100mm ran this phase
  with the earlier slab settings, under which its photosphere was dark (5 L☉ against 2190 L☉ at 3 R\*); the 5.4e44 erg it did not
  radiate sit in its outer layers.*
* The accretion history is the same in all three (0.33–0.34 M☉ 5.0–5.3 yr after the core).

## 17. A residual-based stopping test on the energy rows (follow-up item 4, gate a)

Provisional switch `ERESTOL` (ENV `ERESTOL=<L☉>`, default 0 = off; `newton!` in julia/Trbdf2.jl, tree cd9922c). After a stage solve
has converged on the correction test (unchanged), the measure

    4π | Σ_k (R_gas,k + R_rad,k) + Σ_j u_j R_mom,j | / L☉

is evaluated at the reached state — the energy the solved rows still create or destroy, which is what gate (a) sums (§5a) — and
further Newton iterations are made while it exceeds the tolerance. The iteration ends when the tolerance is met, when an iteration
does not halve the measure (stall), or at the iteration limit; a stalled stage is accepted as before and counted. Nothing is
lagged and the solution of a converged stage is only improved.

Run t_en28: `ECOMPAT=1 EDIAG=2 ERESTOL=1`, tree bb57e3e, 500 steps from the checkpoint, against t_en19 (the same without the test).

| | t_en19, correction test only | t_en28, ERESTOL = 1 L☉ |
|---|---|---|
| stages with \|Σ_k (R_gas + R_rad)\| > 1 L☉ (the Newton's part) | 306 of 1000 | 48 of 1000 |
| the same, 90 % of the stages below / largest [L☉] | 9.6 / 103.9 | 0.85 / 11.4 |
| stages with \|Σ net\| > 1 L☉ (gate a as in §5: the Newton's part plus the round-off of zones 88–100) | 316 | 81 |
| \|Σ net\|, 90 % of the stages below / largest [L☉] | 9.4 / 107.5 | 0.93 / 10.0 |
| stages with \|Σ cres\| > 1 L☉ (round-off of the lab-frame radiation row, zones 88–100; largest) | 39 (4.0) | 21 (4.5) |
| Δt-weighted mean of Σ_k (R_gas + R_rad) [L☉] | 0.0092 | 0.0041 |
| extra Newton iterations | | 639 in 1000 stage solves; 48 stalled above the tolerance, none at the iteration limit |
| wall time of the solver, last 400 steps [s] | 1069 | 1153 (+8 %) |
| time covered [yr]; Newton failures; median Δt of the sampled steps [s] | 12.65; 74; 3.2e4 | 12.25; 72; 3.2e4 |

* **The gate "|Σ| < 1 L☉ at every stage" is still not met**: 48 stages of 1000 stay above it, the largest at 11.4 L☉ (104 without
  the test). The stopping test removes five sixths of the cases and a factor ten of the size; it costs 0.64 Newton iterations per
  stage solve, 8 % of the solver's wall time, and nothing measurable in Δt or in the failure count.
* What is left stalls: the iteration no longer halves the measure. It is in the BDF stage (45 of the 48), in the star (zone 20 is
  the largest contributor in 8 of the first 12), and half of it (25) in the first 50 steps after the restart. The per-zone residuals
  of the gas and the radiation row there are of opposite sign and up to 4.2e3 L☉ each (4.1e4 without the test); the measure is
  their sum. Why the sum stalls at 1–11 L☉ was not found (candidates: the convergence floor of the correction test's
  linear solve in an opaque zone, where the two rows are nearly the same row, §7.3; a term of the rows that is not smooth at
  that level). A stricter stall rule was not tried.
* The round-off of the lab-frame radiation row in the outermost zones (§3.3, §5) is not the
  Newton's and no stopping test reaches it: 21 stages above 1 L☉.
* Time-integrated the Newton's part was negligible before (4e-6 of L_em) and is half of that now; the test buys a cleaner
  per-stage ledger, not a different solution (star at 12.25 yr: M* 1.3405 against 1.355 M☉ at 12.65 yr, r* 108.6 against 108.9 R☉,
  T_c 1.2529e5 against 1.2527e5 K, the difference being the 0.4 yr).

## 18. The lab-frame radiation row alone (round 3, item 1)

ECOMPAT changes two things at once: the gas energy row (the compatible first-law row in every zone, in place of the stack's
dual-energy choice with the conservative total row in the star) and the radiation energy row (lab frame in place of comoving,
§3). §14.7 named the second half alone as the minimal change to the production stack and had not run it. It now has its own
provisional switch.

**Switch.** `ERLAB` (ENV `ERLAB=1`, default off; `TitanResidFast.ERLAB`, commit d676f9b): the radiation energy row of every zone
is the lab-frame row of §3 — the row ECOMPAT uses, the same code path — while the gas rows stay the stack's (`etotal` with
`useint`: total row with the flux-form gravity where the internal energy is resolved, first-law row elsewhere). It implies
`TIEFRAME` as ECOMPAT does (the tied zones' flux is transformed with the same velocity as the row, §3.2). Nothing else changes.
With `ERLAB` off the residual is bit for bit that of the tree before it (§11: four comparisons over 20 steps, 0 fields differ).

**Run.** t_en35: tree d676f9b, the checkpoint z19z100mmS5_auto0 (step 17398), 100 zones, `ERLAB=1 EDIAG=2`, 500 steps, against
t_en24 (stack), t_en19 (ECOMPAT) and t_en25 (ECOMPAT + ECGRAV=1) of §14, all 500 steps from the same checkpoint.

### 18.1 Binding-energy ledger

The ledger of §14.2 (`gledger.jl`: K, U, E_r and an explicit W from the states after every converged stage, no gravity power of
either row), time-weighted means over the run in L☉; W (i) interface-based, (ii) zone-centred, (iv) exact for the
piecewise-constant density:

| | t_en24 stack | **t_en35 stack + ERLAB** | t_en19 ECOMPAT | t_en25 ECOMPAT + ECGRAV=1 |
|---|---|---|---|---|
| d(K + U + E_r)/dt from the states | +139617.6 | **+114559.5** | +102166 | +110147.2 |
| boundary flux out − in | −101865.1 | **−80213.0** | −69797 | −76220.1 |
| other sources | +665.8 | **+652.1** | +654 | +651.3 |
| gravity term of the total row | −34348.0 | **−34214.0** | −30309 | −33673.1 |
| gravity power of the momentum row | −37512.7 | **−37378.5** | −33024 | −36767.4 |
| G5 = the difference of the two | 3164.7 | **3164.4** | 2715 | 3094.4 |
| G5x (exact semi-discrete dW(i)/dt − momentum row's power) | 2946.2 | **2901.6** | 2610 | 2886.9 |
| dW_self/dt (i) / (ii) / (iv) | −34562.0 / −34612.7 / −34354.9 | **−34470.5 / −34518.6 / −34256.8** | −30409.0 / −30420.2 / −30136.7 | −33876.9 / −33912.6 / −33620.4 |
| Σ of the total row's residuals on the solution | 4070.2 | **784.4** | 2715 | 905.2 |
| **ledger defect with W (i) / (ii) / (iv)** | **+3855.2 / +3804.5 / +4062.2** | **+526.9 / +478.9 / +740.6** | **+2613.5 / +2602.3 / +2885.8** | **+700.3 / +664.6 / +956.8** |
| time covered [yr] | 7.42 | 7.92 | 12.65 | 9.65 |
| Newton failures | 88 | 88 | 74 | 76 |

The defect is not steady: its per-step median is −364 / −380 / −144 L☉ and the mean is carried by the first 100 steps after the
restart (cumulative mean with W (iv): 5696 L☉ after 50 steps and 0.10 yr, 1660 after 100 steps and 3.6 yr, 903 after 250 steps,
741 at the end). dW(iv)/dt agrees with the total row's gravity term to −44 L☉ (stack: −7 L☉): against
the exact binding energy of the piecewise-constant density the total row's gravity is the closed one, as found in §14.7, and the
defect is the residual of the rows (784 L☉) minus that.

### 18.2 Residual of the conservative total row per range

Time mean of the residual that the solved rows leave in the conservative total row (gas + radiation + kinetic, flux-form
gravity), in L☉, the layout of §14.7:

| zones | t_en24 stack | **t_en35 stack + ERLAB** | t_en19 ECOMPAT | t_en25 ECOMPAT + ECGRAV=1 |
|---|---|---|---|---|
| 1–26 (star) | 154 | **157** | 1972 | 152 |
| 27–50 (infall) | 620 | **626** | 740 | 750 |
| 51–100 (outer cloud) | 3296 | **1.9** | 3 | 3 |
| all | 4070 | **784** | 2715 | 905 |

Finer for t_en35: zones 1–14 −0.03; 15–26 157.0 (zone 26 alone 157: it is on the first-law row 75 % of the time, the photosphere
zone); 27–33 254.4; 34–50 371.2; 51–100 1.91. The residual that the Newton iteration leaves in the rows it actually solves (total
row where the internal energy is resolved, first-law row + lab-frame radiation row elsewhere) is −0.02 L☉ in the Δt-weighted
mean (zones 88–100: −0.0008), the same as the stack's own (+0.003): the 784 L☉ is the rows', not the iteration's.

### 18.3 The star, robustness

At equal time after the restart (photosphere zone 27; r_ph [R☉] / L_ph [L☉] / L at 1.5 r_ph [L☉] / T_eff [K] / T_c [K]):

| t [yr] | t_en24 stack | **t_en35 stack + ERLAB** | t_en25 ECOMPAT + ECGRAV=1 | t_en19 ECOMPAT |
|---|---|---|---|---|
| 2.0 | 106.79 / 6097 / 6101 / 4935 / 1.2114e5 | **107.01 / 6086 / 6089 / 4928 / 1.2109e5** | 109.06 / 6031 / 6052 / 4871 / 1.2118e5 | 113.23 / 6052 / 6194 / 4784 / 1.2139e5 |
| 4.0 | 105.82 / 6329 / 6511 / 5005 / 1.2238e5 | **105.79 / 6330 / 6516 / 5006 / 1.2237e5** | 107.91 / 6263 / 6587 / 4943 / 1.2262e5 | 113.31 / 6307 / 6727 / 4832 / 1.2280e5 |
| 6.0 | 106.45 / 6374 / 7039 / 4999 / 1.2372e5 | **106.59 / 6419 / 7052 / 5004 / 1.2366e5** | 109.02 / 6507 / 7263 / 4965 / 1.2390e5 | 114.68 / 6611 / 7442 / 4860 / 1.2370e5 |
| 7.4 | 108.13 / 6543 / 7374 / 4992 / 1.2510e5 | **108.60 / 6642 / 7341 / 5000 / 1.2508e5** | 112.32 / 6747 / 7386 / 4936 / 1.2492e5 | 117.44 / 6817 / 7721 / 4840 / 1.2447e5 |

(The sampled times differ by up to 0.1 yr between the runs; ρ_c and the mass inside the photosphere of t_en35 agree with the
stack's to 0.04 % and 0.4 %, the latter being the 0.03–0.09 yr of the sample times.)

500 steps cover 7.92 yr with 88 Newton failures (stack: 7.42 yr, 88).

### 18.4 What it shows

* **Yes: the lab-frame radiation row alone gives the stack's star with the outer cloud's 3.3e3 L☉ removed.** r_ph within 0.5 R☉
  (0.4 %), T_eff within 8 K, T_c within 0.05 %, L_ph within 1.5 % of the stack at equal time, over 7.4 yr; the residual of zones
  51–100 goes from 3296 to 1.9 L☉, zones 1–50 are unchanged (783 against 774 L☉). The steps, the failures and the time covered
  are the stack's.
* **Its ledger is at or below ECOMPAT + ECGRAV=1's**: 479–741 against 665–957 L☉ for the three definitions of W, 12 % of the
  stack's 3.8–4.1e3 L☉. What is left is the same in both and in the stack: the residual of zones 27–50 (620–750 L☉ in all four
  runs), i.e. the gravity power difference G5 of the infall, where the gas is on the first-law row with any of the choices, plus
  the photosphere zone when it is on the first-law row (157 L☉). (For the stack's rows the infall's 626 L☉ is G5 of zones 27–50,
  872 L☉, less the 274 L☉ of kinetic energy that the upwinded momentum flux dissipates at the shock and the stack's first-law
  row does not return, §19.4.) That is the subject of §19.
* **So the two halves of ECOMPAT separate cleanly.** The lab-frame radiation row removes the outer cloud's defect and does
  nothing to the star. The compatible gas row in the star changes the star (r_ph +6–9 R☉, T_eff −150 K) and opens 1.8e3 L☉ of G5
  there, which ECGRAV=1 then has to put back. Of the configurations with the code's gravity weight measured from the checkpoint, the stack with `ERLAB=1` is the
  one with the smallest ledger defect, the smallest change to the production stack, and no change to the star. (With the
  finite-volume weight of §19 on top of it the ledger closes to −0.2…+0.1e3 L☉, at the price of a different central balance.)
* The run from the IC with this configuration is t_ic_l (§20.4, §21): the control's star at equal mass (R\* within 1 %,
  T_c within 0.3 %, n_c within 0.016 dex at M\* 0.10–0.20 M☉), the outer cloud's residual 5 L☉ instead of 3.4e3 L☉, the ledger
  −375 / +313 L☉ (−3.0 / +2.5 % of the release rate) instead of +3217 / +3974 L☉.
* Not measured for this configuration: gates (a)–(d) of §5 as a set (the Newton part of gate a is given above), and the sub-bath effect of §15 (E below the equilibrium value in the infall with the
  comoving row; with the lab-frame row it is the row of §3, so §15's ECOMPAT numbers should apply, not verified here).

## 19. The gravity weight of the momentum row (round 3, item 2)

### 19.1 The three weights

Interface j's momentum row carries the gravity force F_j (per steradian). With g_j = G M_j / r_j², zone masses m_k = ρ_k ΔV_k,
ΔV_k = (r_k+1³ − r_k³)/3 and V̂_j = r_j² · ½(r_j+1 − r_j−1):

    node mass (the code's form)      F_j = g_j μ_j,          μ_j = ½(m_j−1 + m_j)
    finite volume (`GFV`)            F_j = g_j ρ_c,j V̂_j,    ρ_c,j = [ρ_j−1 (r_j − r_j−1) + ρ_j (r_j+1 − r_j)] / (r_j+1 − r_j−1)
    variational (`GVAR`)             F_j = ∂W(iv)/∂r_j at fixed zone masses

and the conservative total row uses the mass flux: its gravity power is Σ_j u_j g_j ρᵇ_j V̂_j with ρᵇ_j the advected density of
interface j, so that G5 = Σ_j u_j (g_j ρᵇ_j V̂_j − F_j) (§14.2).

*Finite volume.* ∫ρ g dV between the centres of zones j−1 and j with the density constant in each zone and r² g = G M_j across
the node: G M_j [ρ_j−1 ½(r_j − r_j−1) + ρ_j ½(r_j+1 − r_j)] = g_j ρ_c,j V̂_j. It is the weight for which a steady flow through a
fixed mesh has the right Bernoulli constant: the momentum row then reads Ṁ (u_j+½ − u_j−½) = −F_j with Ṁ = r_j² ρ u, i.e.
u du/dr = −g · F_j/(g_j ρ V̂_j), and the factor F_j/(g_j ρ V̂_j) multiplies the free-fall u² r/(2GM). The node mass has
μ_j/(ρ V̂_j) = 1.07–1.11 in the infall of the checkpoint (ρ ∝ r^−3/2 on a mesh with ratio 1.26), which is the measured
u² r/(2GM) of up to 1.111 (§14.4).

*Variational.* W(iv) = −∫ G M(r) dm/r for the density constant in each zone (§14.2) is a function of the interface radii and the
zone masses; its derivative with respect to r_j at fixed masses involves the two zones next to the interface only (the mass
inside r_j and the masses of the zones outside do not change) and has the closed form

    F_j = 4πG r_j² [ ρ_j−1 ∫ M(r) (r³ − r_j−1³) / (r² (r_j³ − r_j−1³)) dr + ρ_j ∫ M(r) (r_j+1³ − r³) / (r² (r_j+1³ − r_j³)) dr ]

(`gvar_in`, `gvar_out` in julia/TitanResidFast.jl; M(r) the enclosed mass inside the zone; integrals over the zone below and the
zone above), i.e. the zone's own ∫ρ g dV distributed to its two interfaces with the weights with which a displacement of the
interface displaces the zone's mass shells at fixed zone mass. Checks made: the closed forms against quadrature (relative
difference ≤ 5e-15, 3e-13 for a zone of relative width 1e-4, five parameter sets including a = 0); the gradient against a finite
difference of W(iv) on a record of t_en24 (1.9e-25 relative to W). The external (dark-matter) mass and a uniform field keep the
node mass (their W is not part of W(iv)).

**Switches** (provisional, default off, commits d676f9b and 6642493): `GFV=1`, `GVAR=1` (ENV; `TitanResidFast.GFV`, `.GVAR`;
mutually exclusive). They replace the weight in the momentum row, in the gravity work of the mixing-length row, and in the
diagnostic G5 of `zone_diag`; with both off the residual is bit for bit that of the tree before them (§11). The total row's
flux-form gravity term is not touched.

### 19.2 Is there a form with both limits?

Two limits are wanted. *Lagrangian* (interfaces move with the gas, u_rel = 0): the force should be μ_j g_j, or more precisely
it should be the gradient of a binding energy, so that K + U + W is conserved and a hydrostatic star is the minimum of U + W.
*Eulerian* (interface fixed, mass passes): the force should be the finite-volume weight with the density of the mass flux, so
that the gravity power of the momentum row is that of the mass flux (G5 = 0) and free fall has the right Bernoulli constant.

**No continuous, flag-free force has both limits exactly**, for two independent reasons.

1. *The two limits meet.* Write the force as F(state; ṙ_j, u_rel,j), u_j = ṙ_j + u_rel,j. "F = A(state) whenever u_rel = 0" and
   "F = B(state) whenever ṙ = 0" both apply at ṙ = u_rel = 0, gas at rest on a mesh at rest; continuity there requires A = B as
   functions of the state, and the node-mass and the finite-volume weight are different functions (ratio 0.75–1.05 on the
   checkpoint, below). Any F with both limits is therefore discontinuous at u = 0, ṙ = 0. The form that matches the *power*
   in both limits, F = (A ṙ + B u_rel)/(ṙ + u_rel), is the explicit example: it depends on the direction in which the origin is
   approached, and the origin is the hydrostatic star on a mesh that follows it, where the force has to balance the pressure
   gradient (max |u|/c_s is 9e-4 in the relaxed S5 star).
2. *G5 = 0 needs the advected density.* G5 vanishes for every flow only if F_j = g_j ρᵇ_j V̂_j, and ρᵇ_j (upwind, or any monotone
   reconstruction) changes with the sign of u_rel,j by the density contrast of the two zones. In the star's outer layers that
   contrast is a factor of several per zone: the "density-jump part" of G5, u g V̂ (ρᵇ − ρ_c), is 2.3–2.9e3 L☉ at interfaces 15–27 (§19.4),
   the larger part of all of G5. A force with that weight would jump by tens of per cent whenever the mesh velocity crosses the
   gas velocity at an interface of the hydrostatic layers. (The same discontinuity, in the gas row instead of the momentum row,
   is a likely reason why `ECGRAV=2` and `ECGRAV=3` fail, §14.5; not verified separately.)

What can be had without a flag is a state-only force. The three candidates and what each keeps:

| | node mass g μ | finite volume g ρ_c V̂ (`GFV`) | variational ∂W(iv)/∂r (`GVAR`) |
|---|---|---|---|
| gradient of a binding energy (Lagrangian limit conserves K + U + W) | yes, W (i) | **no**: ∂F_j/∂r_j+1 ≠ ∂F_j+1/∂r_j at fixed masses, relative asymmetry 0.67 at interface 3 of the checkpoint | yes, W (iv) |
| Bernoulli factor F/(g ρᵇ V̂) in the infall of the checkpoint (1 = exact) | 1.07–1.11 | 0.98–1.03 | 0.955–1.025 |
| hydrostatic row, ρ ∝ r^−1.5, mesh ratio 1.26 (the infall's mesh) | +2.9 % | −0.2 % | −2.0 % |
| hydrostatic row, ρ ∝ r^−1.8 … r^−2, mesh ratio 2 (the star's inner zones, §20.1) | +17 … +22 % | 0 … −0.3 % | −13 … −15 % |
| hydrostatic row, uniform density, mesh ratio 2 | +32 % | −25 % | 0 |

(The hydrostatic rows are §19.3. Ratios of the weights on the first record of t_en24, F_var/(gμ) | ρ_cV̂/μ | ρᵇV̂/μ:
interface 2: 0.679 | 0.864 | 0.871; 3: 0.669 | 0.750 | 0.848; 5: 0.814 | 0.933 | 0.835; 8: 0.927 | 0.990 | 0.855;
12: 0.992 | 1.017 | 0.914; 19: 1.045 | 1.046 | 0.928; 23: 1.040 | 1.034 | 0.790; 25: 1.036 | 1.027 | 0.616; 26: 1.024 | 1.018 | 0.721;
27: 1.012 | 1.009 | 0.872; 30: 0.960 | 0.969 | 0.976; 34: 0.879 | 0.901 | 0.920; 36: 0.894 | 0.910 | 0.900; 40: 0.952 | 0.958 | 0.929;
45: 0.961 | 0.975 | 0.941; 52: 0.971 | 0.990 | 0.944. The three weights differ by 25–33 % at the innermost interfaces of the star.)

Literature. That a gravity source term has to be discretised consistently with the pressure gradient for hydrostatic states to
be preserved, and that a zone-centred ρ g is not, is the subject of the well-balanced schemes: Käppeli & Mishra (2014, JCoPh
259, 199) and Käppeli & Mishra (2016, A&A 587, A94) for finite-volume Euler with gravity, Zingale et al. (2002, ApJS 143, 539)
for the mapping of hydrostatic models onto a mesh (abstracts only consulted; ADS, refereed). The energy side — a momentum
source whose power is that of the mass flux so that total energy is conserved to round-off — is Springel (2010), Jiang et al.
(2013), Katz et al. (2016) and Mullen et al. (2021) as cited in §12. The variational (gradient-of-a-discrete-energy) form of
staggered Lagrangian hydrodynamics is that of the compatible schemes (Caramana et al. 1998, §12). None of these treats a mesh
that is neither Lagrangian nor Eulerian; the statement above for that case is derived here, not taken from a source.

### 19.3 What each weight does to a hydrostatic star (computed, no simulation)

**Power-law spheres** (`powstat.jl`). Self-gravitating ρ ∝ r^−a in exact hydrostatic balance on a geometric mesh of ratio q,
zone density = zone mass / zone volume, zone pressure = the volume average of the exact pressure (what a γ-law zone holds as
internal energy). The momentum row asks r_j² (p_j−1 − p_j) = F_j; printed is F_j / [r_j² (p_j−1 − p_j)] − 1, node mass | finite
volume | variational:

| q | a = 0 | a = 1.1 | a = 1.5 | a = 1.8 | a = 2.0 |
|---|---|---|---|---|---|
| 1.05 | +0.0016 / −0.0016 / 0 | +0.0015 / −0.0003 / −0.0008 | +0.0013 / −0.0001 / −0.0009 | +0.0010 / −0.0000 / −0.0009 | +0.0008 / 0 / −0.0008 |
| 1.26 | +0.0355 / −0.0345 / 0 | +0.0350 / −0.0071 / −0.0183 | +0.0292 / −0.0022 / −0.0197 | +0.0229 / −0.0004 / −0.0188 | +0.0179 / 0 / −0.0174 |
| 1.60 | +0.146 / −0.130 / 0 | +0.150 / −0.029 / −0.072 | +0.124 / −0.009 / −0.077 | +0.097 / −0.002 / −0.073 | +0.075 / 0 / −0.068 |
| 2.00 | +0.317 / −0.247 / 0 | +0.343 / −0.060 / −0.145 | +0.283 / −0.019 / −0.153 | +0.218 / −0.003 / −0.146 | +0.167 / 0 / −0.135 |

All three are second order in q − 1. The variational force is exact for uniform density (it is the exact ∫ρ g dV of a
piecewise-constant density); the finite-volume weight is exact for ρ ∝ r^−2 (g ρ r² is then constant) and within 2 % for
a ≥ 1.5; the node mass is the worst of the three for every profile steeper than uniform. No weight is best everywhere.

**Polytrope n = 1.5** (`polystat.jl`) on the mesh of the star of t_en24 (23 zones, central mesh ratios 2.02, 1.92, 1.80, …; the
radii scaled to the polytrope's surface; zone masses from the Lane–Emden solution; rows r_j² (p_j−1 − p_j) = F_j with
p = K ρ^(5/3) of the zone's mean density).

* Residual of the row on the exact radii, relative to F_j, at interfaces 2, 3, …: node −0.234, −0.219, −0.181, −0.153, −0.131,
  −0.113, −0.097, −0.082, −0.069, −0.057, …; finite volume +0.041, +0.276, +0.209, +0.168, +0.140, +0.117, +0.098, +0.080, +0.065,
  +0.051, …; variational +0.000 (six interfaces), +0.001 … +0.029 (outwards, where the density is no longer flat).
* Discrete equilibrium (the rows solved for the radii): density of zones 1–8 relative to the exact mean density of the same
  mass shell 1.016 (node) | 1.030 (finite volume) | 0.981 (variational); surface radius R/R_exact 0.972 | 0.945 | 0.960;
  central pressure 1.027 | 1.050 | 0.968 of the exact one.
* Linear adiabatic radial oscillations about each discrete equilibrium (Γ₁ = 5/3, node-mass inertia, stiffness matrix by finite
  differences in 160-bit arithmetic), ω² in units of 4πGρ_c, three lowest modes: node 0.1480, 0.6614, 1.3677; finite volume
  0.1631, 0.7626, 1.5890; variational 0.1503, 0.7345, 1.5280. The converged values from uniform meshes of 100 and 200 zones are
  0.15056, 0.6976, 1.479 (the three weights agree there to 1e-4, 1e-3 and 1e-3 at 200 zones). On the star's mesh the fundamental
  is off by −1.7 % (node), +8.3 % (finite volume), −0.2 % (variational); the second mode by −5.2 %, +9.3 %, +5.3 %. **No weight
  has an unstable or overstable mode**: all eigenvalues are real and positive for all three, including the finite-volume force,
  which is not the gradient of an energy. (A polytrope n = 3 with Γ₁ = 4/3 is neutrally stable and the discrete equilibrium
  could not be found for any weight; not usable as a test.)

For a centrally flat polytrope the variational force is the accurate one and the finite-volume weight the least accurate; for
the r^−2-like run of density that the young star has between 0.2 and 10 R☉ (§20.1) it is the other way round. On a mesh with
ratio 1.1–1.3 the three agree to 2–3 %; the star's inner zones have ratio 1.8–2.2 at 100 zones.

### 19.4 Measured from the checkpoint (item 2a)

500 steps from the checkpoint z19z100mmS5_auto0, 100 zones, EDIAG=2, tree 6642493; the four runs with a new weight had 4 threads,
the node-weight runs of §14 and §18 had 8. The weight is changed at step 0 of the restart, so each run starts with the star out of the balance of its
new rows; how long that lasts is part of the result.

| run | rows | weight | time covered [yr] | "Newton failure" rejects | failed TR / BDF stages |
|---|---|---|---|---|---|
| t_en24 | stack | node | 7.42 | 88 | (§14) |
| t_en35 | stack + `ERLAB` | node | 7.92 | 88 | (§18) |
| t_en19 | ECOMPAT | node | 12.6 | 74 | (§14) |
| **t_en39** | stack + `ERLAB` | **finite volume** | 12.69 | 86 | 64 / 2 |
| **t_en38** | stack + `ERLAB` | **variational** | 4.07 | 58 | 38 / 0 |
| **t_en37** | ECOMPAT | **finite volume** | 0.275 | 91 | 71 / 0 |
| **t_en36** | ECOMPAT | **variational** | 0.263 | 89 | 69 / 0 |

(Failed stages: the increase of perf.csv's counters over the 500 steps; the checkpoint carries 1238 / 17.)

**The free fall.** u² r/(2GM) at interfaces 30–38 (130–1000 R☉) and the slope d(u²/2)/d(GM/r) between neighbouring interfaces
31–38 (1 = Bernoulli; `uff.jl`):

| | u² r/(2GM), interfaces 30 … 38 | slope, 31 … 38 |
|---|---|---|
| checkpoint (all runs' start) | 0.947 1.111 1.034 1.048 1.027 1.003 0.961 0.906 0.838 | 1.61 0.97 1.118 1.094 1.105 1.091 1.066 1.029 |
| node, after 0.27 yr (t_en19) | 0.962 1.111 1.036 1.050 1.030 1.006 0.965 0.910 0.844 | 1.61 0.97 1.118 1.095 1.107 1.093 1.069 1.033 |
| finite volume, ECOMPAT, after 0.27 yr (t_en37) | 0.887 1.040 0.968 0.983 0.971 0.960 0.933 0.891 0.833 | 1.51 0.90 1.026 1.000 1.025 1.031 1.029 1.010 |
| variational, ECOMPAT, after 0.26 yr (t_en36) | 0.883 1.024 0.954 0.970 0.960 0.953 0.929 0.889 0.832 | 1.49 0.88 1.005 0.980 1.011 1.021 1.024 1.007 |
| variational, stack + ERLAB, after 4.07 yr (t_en38) | 0.899 1.017 0.944 0.960 0.949 0.943 0.923 0.891 0.844 | 1.48 0.87 0.996 0.965 0.990 0.994 0.995 0.983 |
| finite volume, stack + ERLAB, after 12.7 yr (t_en39) | 0.861 1.032 0.974 0.987 0.977 0.974 0.960 0.936 0.899 | 1.53 0.90 1.024 0.988 1.008 1.012 1.017 1.009 |

With the node mass the slope is 1.09–1.12 where the node's μ/(ρᵇ V̂) is 1.07–1.11, as the Bernoulli argument of §19.1 says; with
either new weight it is 1.00 ± 0.03 after one crossing time of the infall (0.27 yr from 500 to 130 R☉) and stays there (12.7 yr
in t_en39). The maximum of u² r/(2GM) goes from 1.111 to 1.03–1.04 (finite volume) and 1.02 (variational); what is left above 1 at
interface 31 is
the first interface outside the shock, where the slope is 1.5 with every weight (the precursor, not the weight).

**G5 per interface range** (L☉, time mean over each run; `g5int.jl`; upw = u g V̂ (ρᵇ − ρ_c), the density-jump part; rem =
u (g ρ_c V̂ − F), what the weight leaves for the centred density; G5 = upw + rem):

| run | 2–14 | 15–27 | 28–51 (upw / rem) | 52–100 | all | G M\* Ṁ / R\* |
|---|---|---|---|---|---|---|
| t_en24 stack, node | 313 | 1936 | **911** (−28 / 939) | 5.0 | 3165 | 12781 |
| t_en35 stack + ERLAB, node | 307 | 1935 | **917** (−29 / 946) | 5.0 | 3164 | 12820 |
| t_en19 ECOMPAT, node | 189 | 1636 | **885** | 5.1 | 2715 | |
| t_en39 stack + ERLAB, finite volume | 258 | 2583 | **−36** (−36 / 0) | 4.5 | 2809 | 13415 |
| t_en38 stack + ERLAB, variational | 189 | 2294 | **−228** (−26 / −202) | 2.8 | 2258 | 11293 |
| t_en37 ECOMPAT, finite volume | 374 | 2932 | **−18** | 4.4 | 3291 | (transient) |
| t_en36 ECOMPAT, variational | −344 | 2350 | **−186** | 2.7 | 1823 | (transient) |

**The estimate of §14.4 is confirmed: the finite-volume weight takes the infall's G5 from 885–917 L☉ to −18 L☉ (ECOMPAT, 0.27 yr)
and −36 L☉ (stack + `ERLAB`, 12.7 yr)** (estimate −21).
The variational weight overshoots to −0.2e3 L☉ (its Bernoulli factor is 0.96–1.03, the table of §19.2). Neither touches the
larger part, interfaces 15–27 (the star's outer layers and the shock): 1.6–1.9e3 L☉ with the node mass, 2.3–2.9e3 L☉ with the new
weights. That part is the density jump (upw 2.3–2.6e3 L☉ with every weight, against which the node mass happens to leave
−0.5e3 L☉), and no state-only weight removes it (§19.2, reason 2). With the finite-volume weight G5 is the density-jump part and
nothing else (rem = 0 by construction): 2809 L☉ = 21 % of G M\* Ṁ / R\*, of which 2583 L☉ at the 13 interfaces 15–27.

**Ledger and residual of the total row** (L☉; `gledgerw.jl` with `G4FIX=1`, `imbzone.jl`):

| run | Σ residual, zones 1–26 | 27–50 | 51–100 | all | ledger defect W (i) / (ii) / (iv) |
|---|---|---|---|---|---|
| t_en24 stack, node | 154 | 620 | 3296 | 4070 | 3855 / 3805 / 4062 |
| t_en35 stack + ERLAB, node | 157 | 626 | 1.9 | 784 | 527 / 479 / 741 |
| t_en25 ECOMPAT + ECGRAV=1, node | 152 | 750 | 3 | 905 | 700 / 665 / 957 |
| **t_en39 stack + ERLAB, finite volume** | 235 | **−331** | 2.2 | **−94** | **−137 / −202 / +105** |
| **t_en38 stack + ERLAB, variational** | 178 | **−492** | −0.8 | **−314** | **+430 / +280 / +436** |
| t_en37 ECOMPAT, finite volume | 3361 | −137 | −1.8 | 3223 | +10785 / +10233 / +10783 (transient) |
| t_en36 ECOMPAT, variational | | | | | +7373 (W (iv); transient) |

**With the stack's gas rows, the lab-frame radiation row and the finite-volume weight the ledger closes: −137 / −202 / +105 L☉
over 12.7 yr (three definitions of W; ≤ 1.5 % of G M\* Ṁ / R\* = 13415 L☉), against 3855 / 3805 / 4062 L☉ for the stack.** The
cumulative defect (W (iv)) is 6.5e3 L☉ after 0.025 yr, 8.7e3 after 0.25 yr, 3.5e3 after 0.55 yr, 57 after 8.1 yr and 105 after
12.7 yr: the restart transient (the star settling into the new balance, below) is over after about a year, and the mean of the
remaining 12 yr is what the table gives. The variational run t_en38 is on the same path (10001 L☉ after 0.04 yr, 1257 after 2.3 yr,
436 at its end, 4.07 yr, still falling); the two ECOMPAT runs cover 0.27 yr and are all transient (t_en37: 33347 L☉ after 50
steps, 10783 at the end).

How the ledger closes, per zone range (`imbzone.jl`; K3a is the kinetic energy the upwinded momentum flux dissipates, which the
first-law row of the stack's infall zones does not return to the gas; it is not in ECOMPAT's residual because the compatible row
returns it):

| zones | t_en35 node: residual = G5 + rest | t_en39 finite volume | t_en38 variational |
|---|---|---|---|
| 1–14 (total row) | −0.0 (G5 351 is inside the row) | −0.0 (430) | 0.0 (404) |
| 15–26 (dual-energy switch, 6–7 % first-law) | 157 | 235 | 178 |
| 27–33 (first-law row; shock and precursor) | 254 = 210 + 282 − 238 (K3a) | −163 = 17 + 63 − 243 | −283 = −37 − 21 − 225 |
| 34–50 (first-law row; infall) | 371 = 662 − 254 − 36 (K3a) | −167 = −106 − 27 − 34 | −209 = −243 + 64 − 29 |
| 51–100 | 1.9 | 2.2 | −0.8 |
| all | 784 | −94 | −314 |

So the −94 L☉ of t_en39 is not a clean zero: it is +235 L☉ from the star's switch zones (their first-law steps carry G5 of the
density jump) against −331 L☉ from the infall, of which −0.28e3 L☉ is the K3a dissipation at the shock that the stack's first-law
row drops (the gas loses kinetic energy that no row books as heat) and −0.09e3 L☉ is the remaining G5 of the first-law zones.
Each part is ≤ 2.5 % of the release rate; that they nearly cancel is not a property of the scheme. The dissipation term is
removed by ECOMPAT's gas row, the switch-zone term by a total row with the density-jump power in it (ECGRAV=1, which is exact for
the finite-volume weight because rem = 0); neither combination was run with `GFV` to a settled state (t_en37 is transient; ECGRAV
is refused with `GFV` at values ≥ 3 only, but ECGRAV=1 + `GFV` was not run).

**The star.** At equal time (`starhist.jl`, photosphere zone 27) r_ph [R☉] / L_ph [L☉] / T_eff [K] / T_c [K] / ρ_c [g cm⁻³]:

| run | 0.26 yr | 2.0 yr | 4.0 yr |
|---|---|---|---|
| t_en35 stack + ERLAB, node | 114.1 / 6024 / 4761 / 1.215e5 / 0.2606 | 107.0 / 6086 / 4928 / 1.211e5 / 0.2589 | 105.8 / 6330 / 5006 / 1.224e5 / 0.2623 |
| t_en38 stack + ERLAB, variational | 112.0 / 5562 / 4710 / 8.670e4 / 0.1475 | 102.5 / 5719 / 4958 / 8.577e4 / 0.1445 | 99.1 / 5955 / 5093 / 8.760e4 / 0.1488 |
| t_en39 stack + ERLAB, finite volume | 115.1 / 5466 / 4626 / 1.026e5 / 0.1892 | 101.3 / 5821 / 5010 / 1.027e5 / 0.1879 | 100.1 / 6018 / 5082 / 1.044e5 / 0.1922 |
| t_en37 ECOMPAT, finite volume | 115.6 / 5647 / 4653 / 1.023e5 / 0.1821 | | |
| t_en36 ECOMPAT, variational | 113.9 / 5701 / 4700 / 8.544e4 / 0.1379 | | |

Later times, same columns (t_en39 against the two node-weight runs that cover them; the last column is the mass inside the
photosphere zone):

| run | 6.7 yr | 7.4–7.9 yr | 12.6 yr | M_ph at 12.6 yr |
|---|---|---|---|---|
| t_en35 stack + ERLAB, node | 107.5 / 6464 / 4991 / 1.245e5 / 0.2680 | 108.9 / 6568 / 4980 / 1.254e5 / 0.2707 (7.90) | | |
| t_en19 ECOMPAT, node | | 117.4 / 6817 / 4840 / 1.245e5 / 0.2678 (7.58) | 121.7 / 7275 / 4831 / 1.253e5 / 0.2698 | 1.356 |
| t_en39 stack + ERLAB, finite volume | 105.3 / 6297 / 5011 / 1.073e5 / 0.1997 | 106.8 / 6413 / 4998 / 1.086e5 / 0.2032 (7.73) | 110.1 / 6992 / 5030 / 1.116e5 / 0.2112 | 1.347 |

**The centre re-equilibrates within the first 100 steps (0.04 yr)**: ρ_c −27…−31 % and T_c −16 % with the finite-volume weight,
−43…−47 % and −28…−30 % with the variational one (p_c −42 % and −60…−63 %), with either set of gas rows, and then follows the
star's growth from there (t_en39: T_c 1.026e5 → 1.116e5 K and ρ_c 0.189 → 0.211 g cm⁻³ while M\* goes from 0.88 to 1.35 M☉; the
node runs are at 1.25e5 K and 0.27 g cm⁻³ throughout). This is the hydrostatic weight of §19.3: the three weights differ by
25–33 % at the interfaces of the star's innermost zones, where the mesh ratio is 2, and the star is in the balance of whichever
one the momentum row uses. The photosphere follows over years and by less: with the finite-volume weight r_ph is 100 R☉ against
106 R☉ at 4 yr and 107 against 109 R☉ at 7.7 yr, L_ph −2…−5 %, T_eff +20…+80 K; that is the stack's star to the accuracy at which
the stack's and ECOMPAT's stars differ from each other (117–122 R☉, 4830 K). With the variational weight at 4 yr: 99 R☉, L_ph −6 %,
T_eff +90 K.

**Robustness.** 500 steps cover 12.7 yr with the finite-volume weight on the stack + `ERLAB` rows, against 7.9 yr with the node
mass on the same rows and 12.6 yr for ECOMPAT; 86 "Newton failure" rejects against 88 and 74; 64 failed TR and 2 failed BDF
stages. The Δt-weighted mean Newton residual of the total row is −0.006 L☉. The first 400 steps (0.55 yr) are the restart
transient, so the step count says more about the transient than about the weight; no run was stopped by the solver.

### 19.5 What is lost with the Lagrangian limit (item 2b)

Three cases: an isolated hydrostatic star on a Lagrangian mesh (the S5 regression), the hydrostatic layers of the accreting
star on the GRIDEQ mesh from the checkpoint (§19.4, the centre paragraph), and the same star grown from the IC (below).

**(1) Isolated 100 M☉ star, Lagrangian mesh, 100 zones** (`s5_run.jl`, N=100, NENV=40; trees e0e6b2e_s5 for the reference and
6642493_s5 for the two weights; the script runs to the ZAMS, I stopped the runs at steps 225 and 165, past the reference's
160 steps). Here the mesh moves with the gas, G5 = 0 identically, and the only question is what the weight does to the
equilibrium and to K + U + W.

| | node mass (reference) | finite volume `GFV=1` | variational `GVAR=1` |
|---|---|---|---|
| relaxed at step (max \|u\|/c_s) | 52 | 52 (9.2e-4) | 46 (9.1e-4) |
| relaxed star: R [R☉] | 6.5262 | 6.4776 (−0.74 %) | 6.4805 (−0.70 %) |
| relaxed star: T_c [K] | 1.0473e8 | 1.0512e8 (+0.37 %) | 1.0454e8 (−0.18 %) |
| relaxed star: ρ_c [g/cm³] | 28.474 | 28.751 (+0.97 %) | 28.338 (−0.48 %) |
| step 160: t [yr] | 205.0 | 194.7 | 200.0 |
| step 160: R [R☉] / log L / log T_eff | 6.4783 / 6.2097 / 4.9080 | 6.4611 / 6.2182 / 4.9107 | 6.4643 / 6.2200 / 4.9111 |
| step 160: T_c [K] / ρ_c | 1.0632e8 / 29.74 | 1.0661e8 / 30.03 | 1.0611e8 / 29.59 |
| **energy ledger closure** (ΔE_tot − (ΣE_bur − ΣF_out)) / Σ\|E_bur\| | 2.0e-5 (2–3.6e-5 throughout) | **1.5e-1 (1.4–1.7e-1 throughout)** | 3.4e-5 (5e-8 … 7.4e-5), with W (iv) |

The structure moves by less than 1 % for both weights (the star has mesh ratios near 1 in its interior and the envelope zones
are log-spaced; this is the 2–3 % agreement of §19.3 for small mesh ratios). The energy ledger does not survive the
finite-volume weight: its force is not the gradient of a binding energy (§19.2), so there is no W whose change equals the
gravity work, and the ledger written with W (i) is open by 0.15 of the nuclear energy released, i.e. (with Σ F_out/Σ E_bur ≈ 57
during the contraction) 2.6e-3 of the energy radiated, against 3.5e-7 for the node mass and for the variational force. The
variational force keeps the closure at the reference's level with its own W (iv). **That is the measured loss: with `GFV` a
Lagrangian star no longer conserves K + U + W to round-off but to 0.3 % of what it radiates; with `GVAR` nothing is lost.**

**(2) The accreting star on the GRIDEQ mesh, from the checkpoint.** §19.4: with either new weight the star's interior
re-adjusts within the first 0.3 yr because the hydrostatic row of the innermost interfaces changes by 25–33 % (mesh ratio
1.8–2.0 there): ρ_c −27…−31 % and T_c −16 % with `GFV`, −43…−47 % and −28…−30 % with `GVAR`, against the node mass; the
photosphere stays where it was (100–110 R☉ with `GFV` against 106–109 R☉ with the node mass on the same rows). These shifts are the size §19.3 computes for the power-law spheres
(+17…+22 % for the node mass, 0 for the finite-volume weight, −13…−15 % for the variational force at mesh ratio 2 and
ρ ∝ r^−1.8…−2): the three weights are three different hydrostatic stars on a 100-zone mesh, and the node mass is the one with
the largest error in the force for the density profile this star has.

**(3) The accreting star grown from the IC, 100 zones**: t_ic_lf (stack + `ERLAB` + `GFV`) against t_ic_l (stack + `ERLAB`,
node weight), the same rows and settings, so that the weight is the only difference and there is no transient of a switched
weight. t_ic_lf ran to M\* 0.356 M☉ (1800 steps, 70 min, 88 / 12 failed TR / BDF stages, 100 Newton-failure lines; the control
t_en_ic0 to M\* 0.34: 2100 steps, 147 / 21) and was stopped there by the sequence; the core forms at 2.41394e5 yr (t_ic_l
2.41416e5).

| M\* [M☉] | weight | R\* [R☉] | r_ph [R☉] | T_c [K] | log n_c | against the 200-zone node-weight star (t_cv200l, §21): R\* / T_c / log n_c |
|---|---|---|---|---|---|---|
| 0.10 | node | 17.46 | 24.34 | 1.034e5 | 23.152 | +1.3 % / −3.9 % / −0.06 |
| 0.10 | finite volume | 14.73 | 27.18 | 8.17e4 | 22.799 | −14.5 % / −24 % / −0.41 |
| 0.15 | node | 24.68 | 32.38 | 1.076e5 | 23.174 | 0 / −2.8 % / −0.05 |
| 0.15 | finite volume | 21.38 | 28.94 | 8.54e4 | 22.828 | −13.3 % / −23 % / −0.40 |
| 0.20 | node | 30.81 | 39.78 | 1.112e5 | 23.193 | +1.3 % / −1.2 % / −0.05 |
| 0.20 | finite volume | 27.83 | 36.53 | 9.00e4 | 22.861 | −8.5 % / −20 % / −0.38 |
| 0.30 | finite volume | 38.37 | 49.21 | 9.76e4 | 22.913 | (control t_en_ic0, node: 41.62 / 52.98 / 1.156e5 / 23.200) |

At M\* 0.20 the mass shells of the finite-volume star lie at 0.535 / 1.437 / 4.064 / 10.25 R☉ (enclosed mass 0.005 / 0.02 /
0.06 / 0.13 M☉), those of the node-weight star at 0.478 / 1.412 / 4.249 / 11.07 R☉, and those of the 200-zone node-weight star
at 0.425 / 1.334 / 4.073 / 10.71 R☉: the finite-volume star is less compact inside 0.02 M☉ (+26 % at 0.005 M☉) and more
compact outside 0.06 M☉ (−4 % at 0.13 M☉, −8.5 % at the surface). Its zone 1 has less than half the density of the
node-weight star's (0.122 against 0.261 g/cm³), at log₁₀ T/ρ^(2/3) = 5.56 against 5.44.

So from the IC, at 100 zones, **the finite-volume weight makes a different young star: central density lower by 0.33–0.35
dex, T_c by 19–21 %, R\* smaller by 10–16 %**, and of the two it is the one further from the 200-zone star in every quantity
of the table. The static tests of §19.3 do not predict a difference of this size or this sign: they give 2–5 % in the central
pressure of a polytrope on the star's mesh for all three weights, and for power-law spheres the node weight is the least
accurate. The difference is there at M\* 0.05 (log n_c 22.761 against 23.084, T_c 7.71e4 against 9.32e4 K, the star 9–10
zones), i.e. it is made while the core forms and first accretes, as the differences of §20.1 are. Not traced further at 100 zones. **At 200 zones the difference is gone**
(t_cv200f against t_cv200l, §21.5): −0.033 / −0.023 / −0.016 dex in n_c and −0.4 / +0.9 / +1.8 % in T_c at M\* 0.10 / 0.15 /
0.20, R\* −6 / −8 / −1 %. The less dense centre is an error of the finite-volume weight on the coarse central mesh of the
100-zone runs (ratio of successive radii 2.2–1.6 at interfaces 2–7; 1.5–1.3 at 200 zones), larger than the node weight's
there, and it converges. Measured on the node-weight runs (end of the window M\* 0.10–0.20, §21.3): the finite-volume weight
ρ_c V̂ is below the node mass by 26 / 10 / 5.4 / 2.3 % at interfaces 3–6 of the 100-zone star (and above it by 1.6–7 % at
interfaces 8–12), and by 12.5 / 5.4 / 0.4 % at interfaces 3–5 of the 200-zone star (within ±1.3 % at interfaces 5–17), while
the node mass over V̂ is within 2.5 % (100 zones) and 0.6 % (200 zones) of the density of the smooth mass distribution there.
For ρ ∝ r⁻² the node weight is exact on any mesh and ρ_c is low by 1.6 / 5.3 / 14 % at ratios of successive radii 1.25 / 1.5 /
2 (vlbias.jl). The innermost interfaces of a `GFV` star thus carry 5–26 % less weight than the smooth distribution has there,
which has the sign of the less dense centre; the step from the weight to the central density was not computed.

Gravity power and ledger of the two runs in the window M\* 0.10–0.20 (the tables of §21.2–§21.4 have them next to the other
runs): G5 0.302 of the release rate with the finite-volume weight (all of it density jump) against 0.261 with the node weight;
of it in the infall 0.011 against 0.043; residual of the conservative total row summed over the mesh 63 against 459 L☉ (the
infall's 428 L☉ of the node weight gone); maximum of u² r/(2GM) in the infall 0.980 against 1.018; binding-energy ledger
−686 / +198 L☉ (W (i) / W (iv); −4.0 / +1.2 % of the release rate) against −375 / +313 L☉ (−3.0 / +2.5 %). Over the longer
window M\* 0.10–0.335 (3.49 yr, 1457 steps) the finite-volume run has G5 4444 L☉ = 0.271 of the release rate (infall 101 L☉),
row residual −9 L☉ and ledger −876 / −207 L☉ (−5.3 / −1.3 %). **From the IC the ledger does not show the gain it shows from
the checkpoint**: the 0.4e3 L☉ of row residual that the weight removes is smaller than the ±0.6e3 L☉ by which dW/dt and the
total row's gravity term disagree in these windows with every weight (§20.2, §21.3).

**(4) The variational weight from the IC, 100 zones**: t_ic_lv2 (stack + `ERLAB` + `GVAR`), same settings as t_ic_l and
t_ic_lf. The core forms at 2.55852e5 yr, 1.44e4 yr (6.0 %) later than with the node weight (the finite-volume weight: 22 yr
earlier); the collapse of the cloud itself is changed, which was not traced. The run reached M\* 0.25 M☉
(about 1580 steps; 55 min and 74 / 6 failed TR / BDF stages to step 1500; M\* 0.20 at step 1200, 2707 s, 60 / 6), where I
stopped it; it passed through stretches of small steps (Δt 5e3–8e3 s) at M\* 0.186–0.188 and at 0.225.

| M\* [M☉] | weight | R\* [R☉] | r_ph [R☉] | T_c [K] | log n_c |
|---|---|---|---|---|---|
| 0.05 | node / finite volume / variational | 9.19 / 7.44 / 7.57 | 47.0 / 50.4 / 53.0 | 9.32e4 / 7.71e4 / 6.01e4 | 23.084 / 22.761 / 22.659 |
| 0.10 | node / finite volume / variational | 17.46 / 14.73 / 14.81 | 24.3 / 27.2 / 31.6 | 1.034e5 / 8.17e4 / 6.34e4 | 23.152 / 22.799 / 22.696 |
| 0.15 | node / finite volume / variational | 24.68 / 21.38 / 21.75 | 32.4 / 28.9 / 30.1 | 1.076e5 / 8.54e4 / 6.57e4 | 23.174 / 22.828 / 22.719 |
| 0.20 | node / finite volume / variational | 30.81 / 27.83 / 27.35 | 39.8 / 36.5 / 35.3 | 1.112e5 / 9.00e4 / 7.06e4 | 23.193 / 22.861 / 22.768 |

The variational weight's star is the least dense at the centre: **log n_c lower than the node weight's by 0.43–0.46 dex, T_c
by 36–39 %, R\* smaller by 11–18 %** (finite volume: 0.32–0.35 dex, 17–21 %, 10–19 %). At M\* 0.20 its mass shells at
0.005 / 0.02 / 0.06 / 0.13 M☉ lie at 0.605 / 1.601 / 4.380 / 10.75 R☉ (node 0.478 / 1.412 / 4.249 / 11.07; finite volume
0.535 / 1.437 / 4.064 / 10.25): less compact inside 0.06 M☉ (+27 % at 0.005 M☉) and more compact outside, the pattern of
the finite-volume star in the inner shells, and stronger. Zone 1 has ρ 0.097 g/cm³ against the node star's 0.261, at
log₁₀ T/ρ^(2/3) = 5.52 against 5.44. This is what the checkpoint run found (§19.4: `GVAR` changes the central balance most),
and it is not a transient of a switched weight. No 200-zone `GVAR` run was made, so whether this difference converges away
as the finite-volume weight's does is not measured.

Window M\* 0.10–0.20, the three weights on the same rows (L☉; release rate G M\* Ṁ / R\*):

| | node (t_ic_l) | finite volume (t_ic_lf) | variational (t_ic_lv2) |
|---|---|---|---|
| steps / time [yr] / release rate | 989 / 1.571 / 12332 | 463 / 1.275 / 17200 | 792 / 1.267 / 17535 |
| G5 (fraction of the release rate) | 3215 (0.261) | 5194 (0.302) | 4115 (0.235) |
| of it star / outer layers / infall | 0.103 / 0.115 / 0.043 | 0.104 / 0.186 / 0.011 | 0.056 / 0.177 / 0.002 |
| residual of the conservative total row, whole mesh | 459 | 63 | −85 |
| ledger defect, W (i) / W (iv) | −375 / +313 | −686 / +198 | −952 / −187 |
| dW/dt − total row's gravity term, W (i) / W (iv) | −834 / −146 | −749 / +135 | −868 / −102 |
| max u² r/(2GM) in the infall, time mean | 1.018 | 0.980 | 0.962 |
| Newton-failure lines (to step) | 149 (2150) | 100 (1800) | 92 (1580) |

Against the energy whose gradient it is, W (iv), the variational weight's ledger is −187 L☉ (−1.1 % of the release rate; +44 L☉
in the shorter window M\* 0.10–0.18), of the size of the node weight's and the finite-volume weight's (+313, +198 L☉);
against W (i) it is the largest of the three (−952 L☉, −5.4 %). The spread of the definitions of W at 100 zones (0.7e3–0.9e3 L☉
in this window for the three weights) is unchanged, and so is dW/dt − gravity term against W (i) (−0.8e3 L☉, §21.3). Its G5
is 0.235 of the release rate, between the other two, with the infall's part gone (0.002) and the star's interior halved.
Its star is the furthest of the three from the 200-zone node-weight star (log n_c 23.239 at M\* 0.20 there).

### 19.6 Result of item 2

1. **Implemented and measured, both weights** (`GFV`, `GVAR`; provisional, default off, bit-identical when off, §11).
2. **(a) From the checkpoint the finite-volume weight does what §14.4 estimated**: the infall's G5 goes from 885–917 L☉ to
   −18…−36 L☉ (estimate −21), the slope d(u²/2)/d(GM/r) from 1.09–1.12 to 1.00 ± 0.03, the maximum of u² r/(2GM) from 1.111 to
   1.03, and with the stack's gas rows and the lab-frame radiation row the binding-energy ledger closes to −137 / −202 / +105 L☉
   (W i / ii / iv; node weight on the same rows 527 / 479 / 741 L☉, §18) over 12.7 yr at the stack's robustness (86
   Newton-failure rejects, 64 / 2 failed stages). The part of G5 that is the density jump (2.6e3 L☉ at interfaces 15–27, 21 %
   of the release rate at 100 zones) is untouched by any weight. The variational weight overshoots in the infall (−0.2e3 L☉) and
   covers a third of the time (4.07 yr).
3. **(b) What is lost with the Lagrangian limit is the binding energy.** On a Lagrangian star (S5) the finite-volume weight
   moves the structure by less than 1 % but K + U + W is no longer conserved to round-off: the ledger is open by 2.6e-3 of the
   energy radiated, against 3.5e-7 for the node weight and for the variational force. On the accreting star's hydrostatic
   layers, which GRIDEQ moves through the gas, each weight is a different hydrostatic star at 100 zones: from the checkpoint
   ρ_c −27…−31 %, T_c −16 % (finite volume), −43…−47 %, −28…−30 % (variational); from the IC −0.33…−0.35 dex in n_c, −19…−21 %
   in T_c and −10…−16 % in R\* (finite volume), −0.43…−0.46 dex, −36…−39 % and −11…−18 % (variational), both further from the 200-zone star than the node weight; at 200 zones the finite-volume star is the node-weight star to
   0.02–0.03 dex in n_c and 2 % in T_c (§21.5).
4. **(c) A form with both limits, continuous in the mesh velocity and in u, without a flag, does not exist** (§19.2): the two
   limits meet at ṙ = u_rel = 0, where a continuous form has to take one value while the two weights differ by a finite
   amount there (7–11 % in the infall, 17–33 % at the star's inner interfaces); and G5 = 0 for a flow through the mesh needs
   the advected density ρᵇ in the weight, which jumps with the sign of u_rel. A weight blended in |u_rel|/|u| is a flag with a
   smooth edge and was not built.
5. **For the default:** neither new weight at 100 zones. The node weight's error in the infall is the stack's remaining row
   residual, 3.5 % of the release rate at 100 zones from the IC and 4–5 % at the checkpoint, and it falls with the zone number
   (1.2 % at 200 zones, §21.3); the finite-volume weight removes it (row residual 0.4 % at 100 zones, −0.1 % at 200) and
   pays at 100 zones with a centre that is 0.3 dex less dense — an error that is gone at 200 zones (§21.5) — and at every zone
   number with the binding energy of the layers the mesh follows. It is a weight for a mesh that resolves the centre (ratio
   of successive radii ≤ 1.5), not for the present 100-zone mesh. `GVAR` keeps the binding energy and is the more
   accurate hydrostatic weight for a flat core, but it leaves a negative G5 in the infall and changes the centre most.

## 20. ECGRAV=1 from the IC; the radius at equal mass (round 3, item 3)

**Run.** t_en_ic3: tree bb57e3e, from the IC, 100 zones, the settings of t_en_ic0 (control, the stack's rows) and t_en_ic1
(ECOMPAT) of §16 with `ECOMPAT=1 ECGRAV=1 EDIAG=2`, 8 threads on cores 16–23. Stopped by me at step 2100 (M\* 0.3396 M☉ by the
driver's definition, the control's last value); 57 min wall, 165 failed TR stages and 7 failed BDF stages (control to step 2100:
147 and 21, 56 min; ECOMPAT failed at step 1478 with M\* 0.331, §16). The core forms at t = 2.41576e5 yr (control 2.42761e5,
ECOMPAT 2.41065e5, z100mm 2.37824e5).

### 20.1 The star at equal mass

M\* and R\* as in the driver (first interface with u < −2 km/s) but interpolated linearly between the two interfaces; values
interpolated in M\* between the two snapshots (every 50 steps) that bracket it (`snapstar.jl`). R\* [R☉] / r_ph [R☉] / T_c [K] /
log n_c:

| M\* [M☉] | control (stack) t_en_ic0 | ECOMPAT t_en_ic1 | **ECOMPAT + ECGRAV=1 t_en_ic3** | z100mm (production) |
|---|---|---|---|---|
| 0.10 | 17.31 / 25.75 / 1.0326e5 / 23.137 | 20.49 / 33.35 / 7.583e4 / 22.878 | **17.88 / 26.55 / 1.0631e5 / 23.086** | 17.69 / 25.76 / 1.0270e5 / 23.129 |
| 0.15 | 24.61 / 33.73 / 1.0737e5 / 23.159 | 28.48 / 36.58 / 7.889e4 / 22.899 | **25.37 / 33.16 / 1.1011e5 / 23.105** | 24.93 / 33.82 / 1.0752e5 / 23.155 |
| 0.20 | 30.85 / 39.81 / 1.1088e5 / 23.177 | 35.36 / 45.58 / 8.088e4 / 22.910 | **31.57 / 39.66 / 1.1316e5 / 23.120** | 32.54 / 41.99 / 1.1237e5 / 23.177 |
| 0.25 | 36.70 / 46.53 / 1.1386e5 / 23.192 | 41.47 / 52.25 / 8.225e4 / 22.917 | **37.76 / 47.14 / 1.1558e5 / 23.130** | 42.99 / 53.29 / 1.1558e5 / 23.186 |
| 0.30 | 41.62 / 52.98 / 1.1564e5 / 23.200 | 46.68 / 56.88 / 8.352e4 / 22.924 | **43.25 / 52.51 / 1.1756e5 / 23.137** | 50.70 / 63.58 / 1.1786e5 / 23.193 |
| 0.33 | 44.06 / 54.36 / 1.1703e5 / 23.206 | 48.86 / 59.59 / 8.440e4 / 22.929 | **46.01 / 55.62 / 1.1872e5 / 23.142** | 55.59 / 73.39 / 1.1801e5 / 23.192 |

Photosphere at M\* 0.33 (L_ph [L☉] / T_eff [K]): control 2236 / 5379; ECOMPAT 2525 / 5300; ECOMPAT + ECGRAV=1 2177 / 5279;
z100mm 35 / 1592 (below).

**Answer.** ECGRAV=1 restores the control's central temperature (+1.4 to +3.0 % above the control at equal mass; ECOMPAT alone
−27 to −28 %) and most of its central density (−0.05 to −0.064 dex, i.e. −11 to −14 %; ECOMPAT alone −0.26 to −0.28 dex,
−45 to −47 %), and the radius to +3 to +4 % (ECOMPAT alone +11 to +18 %). "The same G5 term" (§16) is now a measurement: the sink
of ECGRAV=1 is G5 in the zones where the stack uses the total row, and putting it into ECOMPAT gives the stack's star to the
numbers above. The remaining difference (central pressure −12 %) was not traced; it is there at M\* 0.09 already (step 300–400:
log n_c 23.078 against 23.128 at T_c 1.05e5 against 1.02e5), i.e. it is made during core formation.

**What T_c and n_c are.** They are the values of zone 1, the innermost 0.1 R☉ (`ENVCORE=0.1`, the first interface stays at
0.1000 R☉) with 2e-4 M☉. The star at M\* 0.34 has 17 zones; the mesh ratio is 2.2, 2.1, 1.9, 1.8, 1.7 at the innermost interfaces
and the density falls by a factor of three per zone from zone 2 to zone 9 (ρ ∝ r^−1.9, T within 1.07–1.31e5 K inside 4.4 R☉).
The entropy of zone 1 is the same in the runs: T_c / ρ_c^(2/3) differs by 0.005 dex between ECOMPAT and ECOMPAT + ECGRAV=1
(M\* 0.10–0.13) and by 0.04–0.05 dex between these and the control. The ECOMPAT centre is the same gas at 0.35 of the pressure: at
equal enclosed mass the pressure of the ECOMPAT star is lower by a factor 2.0–2.6 through the whole star (0.02, 0.06, 0.13, 0.30
M☉: 2.6, 2.2, 2.0, 2.0), its mass shells are 25 % further out at 0.05 M☉ and 11 % at the surface. T_c and n_c therefore measure
the compactness of the whole star, and the two stars differ in energy: thermal + gravitational energy inside R\* (zone-centre
estimate, snapshot nearest M\* 0.33) −8.02e45 erg (control) against −6.84e45 erg (ECOMPAT), a difference of 1.2e45 erg or
1.9e3 L☉ over the 5.2 yr since core formation, which is the size of G5 inside 1.25 R\* (1.8e3 L☉ in the ECOMPAT run, 2.7e3 L☉ in the control, §20.2).

### 20.2 G5 and the ledger at equal mass

Windows of M\* (the last stretch of steps inside each; `g5reg.jl`, `gledgerw.jl`, `ledwin.sh`), time means in L☉.

G5 by region over the whole stretch M\* 0.10–0.335 M☉ (4.0–4.3 yr):

| | control | ECOMPAT | ECOMPAT + ECGRAV=1 |
|---|---|---|---|
| steps / time [yr] | 1487 / 4.031 | 1169 / 4.323 | 1447 / 4.099 |
| Ṁ [M☉/yr] / G M\* Ṁ / R\* [L☉] | 5.81e-2 / 12870 | 5.29e-2 / 10275 | 5.72e-2 / 12167 |
| G5, star r < 0.5 R\* (geometric / density jump) | 1170 (1092 / 78) | 718 (643 / 75) | 1118 (1019 / 99) |
| G5, outer layers 0.5–1.25 R\* | 1488 (28 / 1460) | 1059 (20 / 1039) | 1207 (23 / 1185) |
| G5, infall 1.25–300 R\* | 543 (840 / −297) | 387 (615 / −227) | 447 (719 / −271) |
| G5, cloud | 11 | 10 | 10 |
| G5, all / as a fraction of G M\* Ṁ / R\* | 3211 / 0.250 | 2174 / 0.212 | 2783 / 0.229 |
| total row's gravity term | −27860 | −21921 | −27629 |
| maximum of u² r/(2GM) in the infall, time mean | 1.036 | 1.018 | 1.020 |

The binding-energy ledger by window, defect with W (i) / W (iv) and, in brackets, dW/dt − the total row's gravity term with
W (i) / W (iv):

| M\* window [M☉] | control | ECOMPAT | ECOMPAT + ECGRAV=1 |
|---|---|---|---|
| 0.10–0.15 | +3723 / +4581 (−216 / +642) | +2097 / +2641 (−405 / +138) | +75 / +831 (−373 / +383) |
| 0.15–0.20 | +2775 / +3444 (−1113 / −443) | +1807 / +2288 (−405 / +76) | −38 / +616 (−453 / +201) |
| 0.20–0.25 | +2822 / +3491 (−1007 / −339) | +1653 / +2081 (−475 / −48) | −543 / −28 (−949 / −434) |
| 0.25–0.30 | +3344 / +3901 (−453 / +104) | +1632 / +2008 (−336 / +40) | −503 / −6 (−911 / −414) |
| 0.30–0.335 | +2634 / +3069 (−1166 / −730) | +711 / +928 (−1429 / −1213) | −1088 / −779 (−1503 / −1194) |
| **0.10–0.335** | **+3056 / +3691** (−791 / −157) | **+1610 / +2023** (−564 / −151) | **−424 / +119** (−841 / −298) |
| Σ residual of the total row, 0.10–0.335 | 3847 | 2173 | 417 |

(W (ii) lies within 100 L☉ of W (i).) What the table says:

* Over the whole stretch **ECOMPAT + ECGRAV=1 closes the binding-energy ledger to −424 L☉ with W (i) and +119 L☉ with W (iv)**,
  −3.5 % and +1.0 % of G M\* Ṁ / R\*; **ECOMPAT alone is open by +1610 / +2023 L☉** (16–20 % of its release rate), which is its
  G5 (2174 L☉) less the 0.15–0.56e3 L☉ by which the total row's gravity term itself misses dW/dt. The compatible first-law row
  with the node-mass gravity gives the star about 1.6–2.0e3 L☉ that no binding energy pays for; this is the same sign and size as
  at the checkpoint (§14: ledger +2886 L☉ with 2715 L☉ of G5).
* The single windows scatter by ±0.6e3 L☉ about this, with all rows alike: dW/dt − (total row's gravity term) moves between
  +0.6e3 and −1.2e3 L☉ from window to window (W is evaluated on a star of 17 zones whose surface moves through the mesh; each
  window is 0.7–1.0 yr). The last window, 0.30–0.335, is the extreme one, and the first draft of this section, which had only
  that window, read it as "ECGRAV=1 over-corrects by 0.8–1.1e3 L☉ and the conservative star lies between the two". Over the
  4 yr that is not supported: the over-correction is −0.4e3 L☉ (W (i)) or absent (W (iv)), within the difference between the
  definitions of W.
* The control's +3.1…+3.7e3 L☉ is the sum of its gas rows' ledger and the outer cloud's comoving radiation row (§18: 3.3e3 L☉
  at the checkpoint). The run that separates the two from the IC is t_ic_l (§20.4).

So at 100 zones, by the binding-energy ledger over 4 yr of accretion, the stack's star (dual-energy gas rows, or ECOMPAT with
the G5 sink) is the conservative one to 1–4 % of the accretion power, and the ECOMPAT star is not: it is too large and too
cool at the centre because 17–20 % of the accretion power is not taken out of it.

### 20.3 Why the control (44–47 R☉) and z100mm (55–60 R☉) differ on the same rows

They agree until M\* ≈ 0.15 (R\* 24.6 / 24.9 R☉) and separate afterwards; the centre is the same throughout (T_c within 1 %,
log n_c within 0.015). z100mm ran this phase (to about step 13 500) on older commits (1782928, 7547bab) with `SUBZONE=1,2` and
without `SLABFLUX=1 SLABTIE=1`; the run.env of the control differs from z100mm's in exactly these settings (and EDIAG). With
those settings the photosphere of z100mm was dark: L at 3 R\* is 5 L☉ at M\* 0.34 (control 2190 L☉), T_eff 800–1600 K against
5100–5550 K, and the energy radiated since core formation is 1.07e44 erg (mean 170 L☉) against 6.45e44 erg (mean 1090 L☉)
(`lprof.jl`, trapezoid over the snapshots). The 5.4e44 erg that z100mm kept is 7 % of the star's |E_th + W| = 8.0e45 erg and sits
in its outer layers (photosphere at 73 R☉ against 54 R☉). That is the radius difference: z100mm's young star is larger because
the earlier slab treatment did not let its accretion luminosity out, not because of the energy rows. A control with z100mm's
early settings on the current tree was not run.

### 20.4 The stack with the lab-frame radiation row from the IC (t_ic_l)

Item 1 from the IC: tree 6642493, the control's settings with `ERLAB=1 EDIAG=2`, 4 threads on cores 28–31, to step 2150, M\* 0.32 M☉ (82 min; stopped by me in a stretch of small steps, Δt 2e3 s since step 1850). The
core forms at 2.41416e5 yr (control 2.42761e5: the collapse before the core is not the same to the year, as for every pair of
runs here).

The star at equal mass, t_ic_l against the control t_en_ic0 (R\* [R☉] / r_ph [R☉] / T_c [K] / log n_c):

| M\* [M☉] | control (stack) | stack + `ERLAB` | difference R\* / T_c / log n_c |
|---|---|---|---|
| 0.10 | 17.31 / 25.75 / 1.0326e5 / 23.137 | 17.46 / 24.34 / 1.0339e5 / 23.152 | +0.9 % / +0.1 % / +0.015 |
| 0.15 | 24.61 / 33.73 / 1.0737e5 / 23.159 | 24.68 / 32.38 / 1.0759e5 / 23.174 | +0.3 % / +0.2 % / +0.015 |
| 0.20 | 30.85 / 39.81 / 1.1088e5 / 23.177 | 30.81 / 39.78 / 1.1122e5 / 23.193 | −0.1 % / +0.3 % / +0.016 |
| 0.30 | 41.62 / 52.98 / 1.156e5 / 23.200 | 41.76 / 51.29 / 1.157e5 / 23.215 | +0.3 % / +0.1 % / +0.015 |

Gravity power, row residual and ledger in the window M\* 0.10–0.20 M☉ (L☉; `g5reg.jl`, `imbreg.jl`, `gledgerw.jl`):

| | control (stack) | stack + `ERLAB` |
|---|---|---|
| steps / time [yr] / G M\* Ṁ / R\* | 687 / 1.503 / 12983 | 989 / 1.571 / 12332 |
| G5, all (fraction of the release rate) | 3493 (0.269) | 3215 (0.261) |
| residual of the total row: star / outer layers / infall / **cloud** | 0 / 2 / 501 / **3407** | 0 / 26 / 428 / **4.8** |
| Σ residual of the total row | 3911 | 459 |
| **ledger defect, W (i) / (ii) / (iv)** | **+3217 / +3312 / +3974** | **−375 / −288 / +313** |
| dW/dt − the total row's gravity term, W (i) / W (iv) | −694 / +63 | −834 / −146 |

Over the longer window M\* 0.10–0.32 (t_ic_l: 1748 steps, 3.889 yr, release rate 12154 L☉): G5 2928 L☉ (0.241 of the release
rate; star 0.092, outer layers 0.109, infall 0.039), residual of the total row 0 / 84 / 319 / 5.1 L☉ (star / outer layers /
infall / cloud; Σ 408 L☉), ledger defect −584 / −501 / −7 L☉ (W i / ii / iv; −4.8 / −0.05 % of the release rate), dW/dt −
the total row's gravity term −992 / −415 L☉. In the same window the control (1441 steps) has G5 3260 L☉, Σ residual 3854 L☉
and defect +2957 L☉ with W (i); ECOMPAT (t_en_ic1, 1016 steps) has defect +1677 L☉ with W (i).

* **From the IC as from the checkpoint, the lab-frame radiation row alone gives the control's star** (R\* within 1 %, T_c
  within 0.3 %, n_c within 0.016 dex at equal mass, the difference of the core-formation time included) **and removes the outer
  cloud's 3.4e3 L☉** (4.8 L☉ left). The binding-energy ledger of the stack is then −375 / +313 L☉ (−3.0 / +2.5 % of the release
  rate) in this window, the size of ECOMPAT + ECGRAV=1's (+15 / +717 L☉ in the same window, −424 / +119 L☉ over
  M\* 0.10–0.335), without the compatible row and without a sink.
* What is left in the rows is the infall on the first-law row: 428 L☉ (3.5 % of the release rate); G5 of the same zones is
  435 L☉ in this window. The kinetic-energy terms that the stack's first-law row does not return (§4) enter the same residual
  (the upwinded momentum flux's term is −238 L☉ in the infall and +48 L☉ in the outer layers in this window, `imbreg.jl`); the
  residual was not reconciled term by term here as it was at the checkpoint (§18.4). The rest of the defect is dW/dt against
  the total row's gravity term (§21.3).

## 21. Zone-number convergence of the gravity power and of the star (round 3, item 4; convergence tests)

**These runs are convergence tests. Production stays at 100 zones**; the 50- and 200-zone runs exist only to measure how the
quantities of §14 and §20 depend on the zone number.

### 21.1 Set-up

*Mesh.* The 100-zone runs build their mesh with `ENVEUL=0,d,d,1e30 ENVCORE=0.1,d,1e30` and d = ln(1e10)/(N − 1) = 0.2326 (a
logarithmic mesh; `ENVCORE` pins the first interface at 0.1000 R☉) and then let GRIDEQ move it. The 50- and 200-zone runs use
the same construction with d = 0.46992 and d = 0.11571 and everything else unchanged: the same GRIDEQ weights, time-step
controls and switches (the last line of each run.env is in §13). Zone 1 is therefore the same innermost 0.1 R☉ at every zone
number, and T_c, n_c are the values of that zone. 4 threads per run on cores 16–31.

*Start.* From the IC, to a common stage in M\*. The alternative, a conservative remap of the 100-zone checkpoint onto 50 and
200 zones, was not used: the quantities to be compared (the central state, the star's radius) are set during core formation
and the first accretion (§20.1: the difference between the stars is there at M\* 0.09), a remapped star would carry the
100-zone star's entropy profile for many thermal times, and the remap of the tied radiation zones is a transient of its own.
From the IC there is no transient to wait out; the price is that the runs differ in the time at which the core forms and in
the accretion rate at equal mass (given with each number).

*Stage.* Windows of M\* as in §20.2; the common window that all runs reached is **M\* 0.10–0.20 M☉** (except t_cv200g, started last: 0.10–0.15, §21.6); the longer window
0.10–0.335 where a run reached it.

*Runs.*

| zones | the stack's rows | the stack's gas rows + `ERLAB=1` | `ECOMPAT=1` | `ECOMPAT=1 ECGRAV=1` | `ERLAB=1 GFV=1` |
|---|---|---|---|---|---|
| 50 | t_cv50s | — | t_cv50e | t_cv50g2 | — |
| 100 | t_en_ic0 (§16) | t_ic_l | t_en_ic1 (§16) | t_en_ic3 (§20) | t_ic_lf |
| 200 | t_cv200s: **fails before the core forms** | t_cv200l | t_cv200e | t_cv200g (to M\* 0.166; window 0.10–0.15, §21.6) | t_cv200f |

The 200-zone run with the stack's rows unchanged (t_cv200s) does not get through the collapse: Newton failures at
Δt = 2–4e6 s near log n_c = 7.9, far before the core; I stopped it at step 507 and did not patch it. With the lab-frame
radiation row (`ERLAB=1`, the same gas rows) the 200-zone run passes that phase without incident, so the 200-zone member of
the stack's family is t_cv200l. (At 100 zones t_ic_l and t_en_ic0 show what the radiation row alone changes, §20.4.) The cause
of the failure of t_cv200s was not investigated.

### 21.2 The gravity power against zone number

Window M\* 0.10–0.20 M☉, time means, `g5reg.jl` (regions by r/R\* of each record: star < 0.5 R\*, outer layers 0.5–1.25 R\*,
infall 1.25–300 R\*). G5 as a fraction of the release rate G M\* Ṁ / R\* of the same window; the accretion rate of the window
is not the same at the three zone numbers and is given.

| zones | run | rows | steps | time [yr] | Ṁ [M☉/yr] | G M\* Ṁ/R\* [L☉] | G5 [L☉] | **G5 / release rate** | star | outer layers | infall |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 50 | t_cv50s | stack | 522 | 2.416 | 4.12e-2 | 7759 | 6270 | **0.808** | 0.167 | 0.440 | 0.197 |
| 100 | t_en_ic0 | stack | 687 | 1.503 | 6.59e-2 | 12983 | 3493 | **0.269** | 0.100 | 0.122 | 0.046 |
| 100 | t_ic_l | stack + `ERLAB` | 989 | 1.571 | 6.35e-2 | 12332 | 3215 | **0.261** | 0.103 | 0.115 | 0.043 |
| 200 | t_cv200l | stack + `ERLAB` | 324 | 1.292 | 7.71e-2 | 14994 | 1049 | **0.070** | 0.025 | 0.031 | 0.013 |
| | | order 50 → 100 / 100 → 200 | | | | | | **1.63 / 1.90** | 0.70 / 2.04 | 1.94 / 1.87 | 2.21 / 1.65 |
| 50 | t_cv50e | `ECOMPAT` | 276 | 2.449 | 3.97e-2 | 5326 | 3337 | **0.627** | 0.132 | 0.336 | 0.154 |
| 100 | t_en_ic1 | `ECOMPAT` | 463 | 1.656 | 5.92e-2 | 9986 | 2344 | **0.235** | 0.083 | 0.109 | 0.042 |
| 200 | t_cv200e | `ECOMPAT` | 304 | 1.334 | 7.49e-2 | 13966 | 922 | **0.066** | 0.023 | 0.030 | 0.013 |
| | | order 50 → 100 / 100 → 200 | | | | | | **1.42 / 1.83** | 0.67 / 1.85 | 1.62 / 1.86 | 1.88 / 1.70 |
| 50 | t_cv50g2 | `ECOMPAT` + `ECGRAV=1` | 448 | 2.488 | 4.01e-2 | 7194 | 6179 | **0.859** | 0.242 | 0.410 | 0.204 |
| 100 | t_en_ic3 | `ECOMPAT` + `ECGRAV=1` | 763 | 1.510 | 6.60e-2 | 12482 | 3061 | **0.245** | 0.100 | 0.104 | 0.040 |
| | | order 50 → 100 | | | | | | **1.81** | 1.27 | 1.98 | 2.35 |
| 100 | t_ic_lf | stack + `ERLAB` + `GFV` | 463 | 1.275 | 7.80e-2 | 17200 | 5194 | **0.302** | 0.104 | 0.186 | 0.011 |
| 200 | t_cv200f | stack + `ERLAB` + `GFV` | 264 | 1.216 | 8.16e-2 | 16496 | 1240 | **0.075** | 0.032 | 0.043 | −0.0005 |
| | | order 100 → 200 | | | | | | **2.01** | 1.70 | 2.10 | — |

("Order" is log₂ of the ratio of the fractions at N and 2N; for the stack's family the 100-zone member is t_ic_l, with the same
radiation row as the 200-zone member. t_en_ic0 and t_ic_l differ by 3 % in the fraction: the radiation row does not enter G5.)

The two parts asked for, with the node-weight split of §14.4 (geometric part u g ρ̄ (V̂ − V̄), density part u g V̂ (ρᵇ − ρ̄)), as
fractions of the release rate at 50 / 100 / 200 zones and the two orders:

| | stack family (t_cv50s, t_ic_l, t_cv200l) | `ECOMPAT` (t_cv50e, t_en_ic1, t_cv200e) |
|---|---|---|
| **density-jump part of the outer layers** (0.5–1.25 R\*) | 0.393 / 0.111 / 0.0294; order **1.82 / 1.92** | 0.277 / 0.106 / 0.0282; order **1.39 / 1.90** |
| **geometric part of the infall** (1.25–300 R\*) | 0.196 / 0.0625 / 0.0242; order **1.65 / 1.37** | 0.174 / 0.0624 / 0.0232; order **1.48 / 1.43** |
| geometric part of the star (< 0.5 R\*) | 0.162 / 0.092 / 0.0214; order 0.82 / 2.10 | 0.131 / 0.070 / 0.0193; order 0.90 / 1.86 |
| density part of the infall | +0.002 / −0.020 / −0.011 | −0.019 / −0.021 / −0.010 |
| interfaces in star / outer layers / infall (end of window) | 6 / 3 / 14; 10 / 8 / 26; 21 / 18 / 47 | 6 / 4 / 14; 10 / 8 / 27; 21 / 18 / 48 |

* **G5 falls from 63–81 % of the release rate at 50 zones to 23–27 % at 100 and 6.6–7.0 % at 200**, with the stack's rows and
  with ECOMPAT alike (G5 is a property of the momentum row's weight and the advected density, not of the energy row). The order
  between 100 and 200 zones is 1.8–1.9. At that order 1 % of the release rate needs about 550 zones on this mesh construction.
* **The density-jump part of the outer layers converges at second order** between 100 and 200 zones (1.9), and at 1.4–1.8
  between 50 and 100, where the outer layers have 3–4 interfaces. This corrects the scaling given in §14.6: the part is
  proportional to the density jump at one interface, but the advected density ρᵇ is a second-order (monotone) interpolation, so
  in resolved layers ρᵇ − ρ̄ falls as the square of the zone width. It is first order only where the limiter acts.
* **The geometric part of the infall converges at order 1.4–1.6**, less than the (Δln r)² of smooth flow (§14.6), at all three
  zone numbers and in both families. Counted against the number of interfaces actually in the region (14, 26, 47) instead of
  the zone number the orders are 1.8 and 1.6. The region 1.25–300 R\* contains the upstream side of the AV-smeared shock, whose
  width in zones does not change with the zone number; I did not separate its contribution from that of the smooth infall, so
  the lower order is reported, not explained. The density part of the infall is negative and half as large (−2 % → −1 %), so
  that the infall's total (last column of the first table) falls at 1.7–2.2.
* **The star's interior is not in the asymptotic range at 50 zones** (six interfaces inside 0.5 R\*: order 0.7–0.9 from 50 to
  100, 1.9–2.1 from 100 to 200).
* With the finite-volume weight at 100 zones G5 is 0.302 of the release rate, larger than with the node weight (0.261): `GFV`
  removes the infall's part (0.011 against 0.043) and the node weight's own part in the star and the outer layers, which has
  the opposite sign to the density jump there (split of §19.4 for t_ic_l: density jump against the finite-volume weight +3892
  L☉ = 0.316, node weight's own part −677 L☉ = −0.055). With `GFV` G5 is the density jump alone, and that part of the node-weight
  runs is 1.01 / 0.316 / 0.076 of the release rate at 50 / 100 / 200 zones (orders 1.68 / 2.05). Measured with `GFV` at 200 zones (t_cv200f): 0.075 of the release rate, order 2.0 from the 100-zone
  `GFV` run, the infall's part −8 L☉.

### 21.3 The binding-energy ledger against zone number

Same window; ledger defect (Δ(K + U + E_r + W) minus what crossed the boundary, time mean) with W (i) / W (iv) (`gledgerw.jl`;
W (ii) lies between), in L☉ and as a fraction of the release rate. Its two constituents are the residual of the conservative
total row summed over the mesh (energy the rows create) and dW/dt minus the total row's gravity term (how well the discrete
W follows the flux form):

| zones | run | defect W (i) / W (iv) [L☉] | as a fraction | Σ row residual [L☉] (cloud / infall / outer layers) | dW/dt − gravity term, W (i) / W (iv) |
|---|---|---|---|---|---|
| 50 | t_cv50s, stack | +4658 / +5503 | +0.60 / +0.71 | 3607 (2896 / 711 / 0) | +1051 / +1896 |
| 100 | t_en_ic0, stack | +3217 / +3974 | +0.25 / +0.31 | 3911 (3407 / 501 / 2) | −694 / +63 |
| 100 | t_ic_l, stack + `ERLAB` | **−375 / +313** | −0.030 / +0.025 | 459 (5 / 428 / 26) | −834 / −146 |
| 200 | t_cv200l, stack + `ERLAB` | **−616 / −471** | −0.041 / −0.031 | 198 (−6 / 180 / 24) | −816 / −670 |
| 50 | t_cv50e, `ECOMPAT` | +4071 / +4585 | +0.76 / +0.86 | — | +733 / +1248 |
| 100 | t_en_ic1, `ECOMPAT` | +1939 / +2448 | +0.19 / +0.25 | — | −405 / +105 |
| 200 | t_cv200e, `ECOMPAT` | **+331 / +470** | +0.024 / +0.034 | 921 | −591 / −452 |
| 50 | t_cv50g2, `ECOMPAT` + `ECGRAV=1` | +1585 / +2460 | +0.22 / +0.34 | 723 | +862 / +1737 |
| 100 | t_en_ic3, `ECOMPAT` + `ECGRAV=1` | +15 / +717 | +0.001 / +0.057 | 430 | −415 / +286 |
| 100 | t_ic_lf, stack + `ERLAB` + `GFV` | −686 / +198 | −0.040 / +0.012 | 63 | −749 / +135 |
| 200 | t_cv200f, stack + `ERLAB` + `GFV` | −769 / −600 | −0.047 / −0.036 | −23 (−7 / −71 / 54) | −748 / −578 |

* **The outer cloud's comoving radiation row creates 2.9e3 L☉ at 50 zones and 3.4e3 L☉ at 100 zones from the IC**, as at the
  checkpoint (3.3e3 L☉, §2, §18), and with `ERLAB` it is 5 and −6 L☉: the term does not depend on the zone number, and the
  lab-frame row removes it at each (`imbreg.jl`, residual of the total row in the zones beyond 300 R\*).
* **The rows' part converges.** The stack's row residual in the infall (G5 of the first-law zones and the kinetic-energy terms
  that the stack's first-law row does not return, §18.4) is 711 / 428 / 180 L☉ = 9.2 / 3.5 / 1.2 % of the
  release rate. ECOMPAT's defect is its whole G5 plus the same W term: +0.76 / +0.19 / +0.024 of the release rate with W (i),
  and at 200 zones it differs from the stack + `ERLAB` defect by 947 L☉, which is G5 (922–1049 L☉).
* **What does not converge is dW/dt minus the total row's gravity term**: +1.05e3 / −0.83e3 / −0.82e3 L☉ with W (i) and
  +1.90e3 / −0.15e3 / −0.67e3 L☉ with W (iv) in the stack's family, i.e. −5 % of the release rate at 200 zones with either
  definition. The spread between the definitions of W does fall (845 / 688 / 145 L☉). It is the same in the ECOMPAT run
  (−591 / −452 L☉) and with the finite-volume weight (−748 / −578 L☉), and it is not a scatter between windows: the next
  window of t_cv200l, M\* 0.20–0.28 (596 steps, 1.27 yr), gives −826 / −775 L☉ (0.10–0.28: −821 / −722 L☉, row residual 143
  L☉, defect −676 / −578 L☉ = −4.6 / −3.9 % of the release rate, G5 0.063). At the checkpoint (M\* 0.88 M☉) it is −0.2e3 L☉
  with W (i) and −0.01e3 … −0.04e3 L☉ with W (iv) (§18.1).

**Where this term is (`ledfull.sh`, `g5grid.jl`).** By construction of the ledger, dW(i)/dt − (momentum row's power) = G5x +
(time discretisation), and (total row's gravity term) − (momentum row's power) = G5, so

    dW(i)/dt − total row's gravity term = G5x + (time discretisation) − G5,

and with the velocity of each interface split into the mesh velocity and the velocity of the gas against the mesh, u = ṙ +
u_rel, G5 = Σ u_j (g ρᵇ V̂ − g μ)_j is the sum of a mass-flux part (u_rel) and a mesh-motion part (ṙ). G5x (§14.3) is the
mass-flux part up to the difference between the discrete potential difference and g Δr. Window M\* 0.10–0.20, L☉:

| zones | run | G5 | = mass-flux part | + mesh-motion part | G5x | time discr. | dW(i)/dt − momentum row | dW(i)/dt − total row | mesh-motion part / release rate | G5x / release rate |
|---|---|---|---|---|---|---|---|---|---|---|
| 50 | t_cv50s, stack | 6270 | 6849 | −580 | 7320 | 0.4 | +7320 | +1051 | −0.075 | 0.943 |
| 50 | t_cv50e, `ECOMPAT` | 3337 | 3746 | −409 | 4071 | 0.3 | +4071 | +733 | −0.077 | 0.764 |
| 50 | t_cv50g2, `ECOMPAT` + `ECGRAV=1` | 6179 | 6615 | −436 | 7041 | 0.4 | +7041 | +862 | −0.061 | 0.979 |
| 100 | t_en_ic0, stack | 3493 | 2663 | +831 | 2796 | 2.7 | +2799 | −694 | +0.064 | 0.215 |
| 100 | t_ic_l, stack + `ERLAB` | 3215 | 2272 | +943 | 2379 | 2.2 | +2381 | −834 | +0.076 | 0.193 |
| 100 | t_en_ic1, `ECOMPAT` | 2344 | 1829 | +515 | 1938 | 0.9 | +1939 | −405 | +0.052 | 0.194 |
| 100 | t_en_ic3, `ECOMPAT` + `ECGRAV=1` | 3061 | 2524 | +537 | 2643 | 1.8 | +2645 | −416 | +0.043 | 0.212 |
| 200 | t_cv200l, stack + `ERLAB` | 1049 | 209 | +840 | 222 | 11.5 | +233 | −816 | +0.056 | 0.0148 |
| 200 | t_cv200e, `ECOMPAT` | 922 | 296 | +625 | 319 | 12.0 | +331 | −591 | +0.045 | 0.0228 |
| 200 | t_cv200f, stack + `ERLAB` + `GFV` | 1240 | 387 | +853 | — | — | — | −748 | +0.052 | — |

The two parts per region, stack + `ERLAB` at 100 / 200 zones (t_ic_l / t_cv200l) and ECOMPAT (t_en_ic1 / t_cv200e), L☉:

| region | mass-flux part, stack + `ERLAB` | mesh-motion part, stack + `ERLAB` | mass-flux part, `ECOMPAT` | mesh-motion part, `ECOMPAT` |
|---|---|---|---|---|
| star < 0.5 R\* | +328 / **−423** | **+939 / +798** | +298 / −274 | +528 / +594 |
| outer layers 0.5–1.25 R\* | +1471 / +462 | −58 / +6 | +1152 / +420 | −61 / −1 |
| infall 1.25–300 R\* | +496 / +183 | +28 / +19 | +400 / +164 | +17 / +16 |
| cloud | −23 / −14 | +34 / +17 | −20 / −14 | +30 / +17 |

* **The term that does not converge is the mesh-motion part of G5 in the star's interior**, and it is geometric (ṙ g ρ̄ (V̂ − V̄):
  1002 L☉ at 100 zones, 878 L☉ at 200). In the outer layers and the infall the mesh moves slowly against the flow, G5 is its
  mass-flux part, and there the total row's gravity term follows W (i) to 6–20 L☉ at 200 zones. Inside 0.5 R\* the gas is
  nearly at rest and GRIDEQ moves the mesh inward through it. At r = 0.1–0.45 R\*, end of the window: mesh velocity
  −0.002…−0.005 km/s at 50 zones, −0.02…−0.07 km/s at 100, −0.04…−0.15 km/s at 200 (gas: −0.003…−0.07 km/s at 200 zones);
  V̂/V̄ − 1 = −47…−16 %, −23…−3.5 %, −5.6…−1.3 %. The geometric factor falls at second order and the mesh velocity doubles with
  the zone number, so the product, as a fraction of the release rate, goes from 0.076 to 0.056 (ECOMPAT 0.052 to 0.045).
* **It is not a drift of the discrete W (static test, wdrift.jl).** A term weighted with the mesh velocity could be a property
  of the ledger and not of the rows: W is defined on the mesh, and the W of the same star on another mesh is another number.
  Test, with the runs' own meshes and mesh motion: for every step of the window the gas of the end state is taken as a fixed
  smooth distribution (cubic Hermite interpolation of log M against log r through the enclosed masses at the interfaces,
  three-point slopes, exact for power laws, M ∝ r³ inside the first interface); its zone masses are evaluated on the mesh of
  the start of the step, and W(end mesh) − W(start mesh) is summed over the window. A W that depended on the mass
  distribution alone would give zero. L☉, W (i) / (ii) / (iv), self-gravity of the gas, window M\* 0.10–0.20:

  | zones | run | dW/dt of the run | drift, every interface moved | drift, only the interfaces inside 0.5 R\* moved |
  |---|---|---|---|---|
  | 50 | t_cv50s, stack | −14906 / −14906 / −14061 | +93 / +151 / +396 | +45 / +12 / +146 |
  | 100 | t_ic_l, stack + `ERLAB` | −28721 / −28635 / −28033 | +7 / +107 / +129 | −9 / +61 / +46 |
  | 200 | t_cv200l, stack + `ERLAB` | −34515 / −34477 / −34370 | −14 / +25 / −38 | −18 / +18 / −53 |
  | 200 | t_cv200e, `ECOMPAT` | −31557 / −31525 / −31418 | −9 / +27 / −13 | −13 / +18 / −31 |

  For W (i), the definition the identity above is written for, the drift is at most 18 L☉ at 100 and 200 zones (2 % of the
  term); for the other two at most 129 L☉ at 100 zones and 53 L☉ at 200 (1.0 % and 0.4 % of the release rate). The −0.8e3 L☉
  are not the mesh moving through the star. Two other interpolants were tried and are not used: the piecewise-constant
  density itself (drift −6.1e3 / −6.8e3 / −9.2e3 L☉ at 50 / 100 / 200 zones: a mesh moving through a staircase, which is what
  a donor-cell remap would do, not a smooth star) and a cubic Hermite interpolant in r³ (+1.8e3…+2.0e3 L☉ at 100 and 200 zones:
  computed by hand for ρ ∝ r⁻² on a mesh with ratio of successive radii 1.4, its three-point slope gives an interface density
  17 % too high, 0.3907 against 0.3333 in units of M_j/r_j³, so its star is not the smooth one on this mesh).
* **It is the interface density of the transport.** The exact rate of W (i) is Σ ṙ_j g_j μ_j + Σ F_j Δφ_j (§14.3) and the total
  row books Σ u_j g_j ρᵇ_j V̂_j, so apart from the difference between Δφ and g Δr their difference is the mesh-motion part of
  G5, Σ ṙ_j g_j (ρᵇ_j V̂_j − μ_j), with the opposite sign: it vanishes when the advected interface density times the
  finite-volume volume is the node mass. Measured against the density of the smooth distribution at the interface,
  ρ_s = (d ln M/d ln r) M_j / r_j³ from the same interpolant (end record of the window, interfaces of the star's interior):

  | | 100 zones (t_ic_l), interfaces 3–12, r = 0.007–0.55 R\* | 200 zones (t_cv200l), interfaces 6–18, r = 0.015–0.26 R\* |
  |---|---|---|
  | ratio of successive radii | 2.28 → 1.26 | 1.39 → 1.22 |
  | advected density of the total row and of the mass flux, ρᵇ/ρ_s − 1 | −29 … −7.0 % | −8.9 … −3.7 % (−3.0 … −1.6 % at 0.31–0.57 R\*) |
  | volume-weighted mean of the two zones, ρ̄/ρ_s − 1 | −45 … −0.3 % | −10.7 … −3.0 % |
  | node mass over the finite-volume volume, (μ/V̂)/ρ_s − 1 | −4.3 … +3.5 % (−2.5 … +1.7 % from interface 5) | −1.5 … −0.2 % (within ±0.6 % for interfaces 7–24) |
  | V̂/V̄ − 1 | −47 … −2.0 % | −9.3 … −2.8 % |
  | mesh velocity / gas velocity [km/s] | −0.0004 … −0.08 / −0.0006 … −0.20 | −0.004 … −0.10 / −0.0001 … −0.014 |

  The node mass is consistent with the smooth distribution and the advected density is below it by about V̂/V̄ − 1. μ/V̂ and ρ_s
  are both centred estimates of dM/dr from the same enclosed masses, so their agreement is not an independent statement;
  the independent one is the closed case (vlbias.jl): for ρ ∝ r^(p−3) on a mesh with ratio q of successive radii and the exact
  zone averages as zone densities, the van Leer value of `advect!` (slope 2ab/(a + b) of the differences of the neighbouring
  zone values, i.e. interpolation in the zone index, no geometric weights; the same value for flow in either direction for a
  power law) is, relative to the exact density at the interface,

  | | q = 1.12 | 1.25 | 1.5 | 2.0 |
  |---|---|---|---|---|
  | van Leer in the zone index, p = 0.5 / 1 / 1.5 / 2 | −1.5 / −1.1 / −0.7 / −0.4 % | −5.5 / −4.0 / −2.8 / −1.6 % | −17 / −13 / −8.7 / −5.3 % | −40 / −31 / −23 / −14 % |
  | node mass over V̂, same p | −0.2 / 0 / +0.3 / +0.6 % | −0.6 / 0 / +1.0 / +2.5 % | −2.0 / 0 / +3.4 / +8.3 % | −5.7 / 0 / +10 / +25 % |
  | Δr-weighted mean ρ_c (the weight of `GFV`), same p | −0.3 / −0.4 / −0.5 / −0.4 % | −1.0 / −1.6 / −1.8 / −1.6 % | −3.4 / −5.3 / −5.9 / −5.3 % | −9.4 / −14 / −16 / −14 % |

  (all zero for uniform density, p = 3, except μ/V̂: +1.7 / +6.8 / +23 / +75 %, the V̄/V̂ of §14.6). The star's interior has
  d ln M/d ln r = 0.76–1.4 at these interfaces at 200 zones (0.47–2.1 at 100), where the node weight is exact (p = 1, any mesh) or nearly so, and the van
  Leer value is low by the measured amount. So in the star's interior the momentum row's weight is the accurate one and the
  density that the mass flux and the total row's gravity term use is 2–9 % (200 zones) and 7–29 % (100 zones) too low.
* **Most of G5, at every zone number, is this.** G5 evaluated on the runs' own states (stage records, weights and regions of
  the tables above) with ρ_s in place of ρᵇ, L☉, ρᵇ → ρ_s:

  | zones | run | star < 0.5 R\* | of it mesh-motion part | outer layers | infall | whole mesh | whole mesh / release rate |
  |---|---|---|---|---|---|---|---|
  | 50 | t_cv50s, stack | 1298 → −110 | −413 → +37 | 3413 → 704 | 1531 → −285 | 6269 → 307 | 0.81 → 0.040 |
  | 100 | t_ic_l, stack + `ERLAB` | 1266 → −55 | +939 → −88 | 1413 → 452 | 524 → −86 | 3215 → 311 | 0.261 → 0.025 |
  | 100 | t_en_ic1, `ECOMPAT` | 826 → −33 | +528 → −47 | 1090 → 345 | 417 → −78 | 2344 → 233 | 0.235 → 0.023 |
  | 200 | t_cv200l, stack + `ERLAB` | 375 → −8 | +798 → −62 | 469 → 120 | 202 → −7 | 1049 → 105 | 0.070 → 0.0070 |
  | 200 | t_cv200e, `ECOMPAT` | 320 → −8 | +594 → −47 | 419 → 108 | 180 → −5 | 922 → 94 | 0.066 → 0.0067 |

  With an interface density that is consistent with the enclosed masses the two rows' gravity powers would differ by a tenth
  of what they differ by now at 100 and 200 zones (a twentieth at 50), and the mesh-motion term of the star's interior would
  be −62 and −88 L☉ instead of +798 and +939. What remains is in the star's outer layers (0.8 % of the release rate at 200
  zones, 3.7 % at 100), where the profile is not smooth on the mesh and ρ_s is not a better density than ρᵇ. This is an
  evaluation on states that were computed with ρᵇ, not a run: a transport with such a density was not implemented (below).
  It also says what the "geometric part" and the "density-jump part" of §14.6 and §21.2 are in the star and the infall: one
  thing, the advected density against the node mass, split at the volume-weighted mean ρ̄, which is itself below the
  interface density by V̂/V̄.
* **What it does to the star, derived and not measured.** With an interface density low by a fraction ε the mass flux through
  an interface that moves inward through gas at rest is too small by ε: the mesh carries mass inward with it. Either the
  discrete star takes that (its binding energy grows, unbooked: the term), or it keeps its mass distribution, which needs a
  gas velocity u = −ε ṙ outward on which gravity does the work Σ u g μ, taken from U by the first-law row as p dV and by the
  total row through its gravity term (the same term to first order in ε). In both cases K + U + W + E_r loses the term, in
  both formulations, which is what the ledgers show (−816 L☉ stack + `ERLAB`, −748 with `GFV`, −591 ECOMPAT at 200 zones). Which
  of the two the star does was not separated: ε |ṙ| is at most 0.005 km/s at 200 zones, the size of the gas velocity of the
  contracting interior. The term is 9 % of the total row's gravity power inside 0.5 R\* (−9.0e3 L☉), the size of ε.
* **Why it does not fall from 100 to 200 zones:** ε falls by 3–4 (second order in the ratio of successive radii) while GRIDEQ's
  mesh velocity in the interior doubles or more, so the sum over the interior goes from 0.076 to 0.056 of the release rate
  (ECOMPAT 0.052 to 0.045). Why the mesh velocity doubles was not investigated.
* **Candidate remedy, not implemented and outside round 3:** an interface density for the transport that is consistent with
  the mesh geometry — the derivative of the cumulative mass at the interface, upwinded and limited, which is how
  reconstructions for spherical meshes are built (Mignone 2014; Colella & Woodward 1984 for the cumulative-profile
  construction; §12). It is a change of `advect!`, which every advected quantity of every row uses, so it is a change of the
  scheme and not of the energy rows, and it is not a patch of the gravity term. The dynamical test of it, and of the two
  readings above — a hydrostatic star with a prescribed mesh motion through it — was not run.
* **G5x is not a single term at 200 zones.** It falls from 0.19–0.22 of the release rate at 100 zones to 0.015–0.023 at 200,
  faster than second order, because its part in the star's interior changes sign (+2.7 % → −2.8 %; ECOMPAT +3.0 % → −2.0 %)
  and at 200 zones cancels the outer layers' (+3.1 %, order 1.95 from 100 zones) and the infall's (+1.2 %, order 1.7).
* **So the 200-zone ledgers close to 2–5 % of the release rate** (stack + `ERLAB` −4.1 / −3.1 %, `GFV` −4.7 / −3.6 %, ECOMPAT
  +2.4 / +3.4 %; W (i) / W (iv)), and they are made of different things. The stack's: +1.3 % row residual, nearly all of it in the
  infall (with `GFV` −0.1 %: +0.3 % in the outer layers, −0.4 % in the infall), and −5.4 % from dW/dt against its total row's
  gravity term, of which −798 of −816 L☉ in the star's interior. ECOMPAT's: G5x, i.e. +4.2 % from the outer layers and the
  infall and −2.0 % from the star's interior.
  **At 200 zones the ledger's total does not discriminate between the formulations**; at 100 zones it does (−3.0 / +2.5 %
  against +19 / +25 %), and the parts that are unambiguous do at both (row residual of the infall and the cloud; G5 of the
  outer layers and the infall).

### 21.4 The infall against free fall

Maximum of u² r / (2 G M) over the interfaces of the infall (1.25–300 R\*), time mean over the window:

| zones | stack family | `ECOMPAT` | stack + `ERLAB` + `GFV` |
|---|---|---|---|
| 50 | 1.236 | 1.100 | — |
| 100 | 1.018 (both t_en_ic0 and t_ic_l) | 1.009 | 0.980 |
| 200 | 0.959 | 0.954 | 0.945 |

The gas of the infall does not start from rest at infinity, so the physical value is below 1 (it is 0.95–0.96 at 200 zones
with the node weight). Taking the 200-zone value as the reference, the node weight's excess is +0.28 / +0.06 at 50 / 100 zones for the stack's
family (order 2.2) and +0.15 / +0.055 for ECOMPAT (order 1.4); the 11 % of the checkpoint (§14.4) belongs to a mesh with Δln r up to 0.37 in the
infall. With the finite-volume weight the maximum is 0.980 at 100 zones and 0.945 at 200: at 200 zones the two weights differ
by 0.014. ECOMPAT + ECGRAV=1 at 50 zones: 1.201.

### 21.5 The star at equal mass against zone number

`snapstar.jl` (interpolated in M\* between the snapshots, every 50 steps) at M\* 0.20 M☉. T_c and n_c are the values of zone 1,
the same innermost 0.1 R☉ at every zone number.

| zones | run | t − t_core [yr] | R\* [R☉] | r_ph [R☉] | L_ph [L☉] | T_eff [K] | T_c [K] | log n_c |
|---|---|---|---|---|---|---|---|---|
| 50 | t_cv50s, stack | 4.14 | 31.18 | 48.78 | 919 | 4473 | 1.007e5 | 22.705 |
| 100 | t_en_ic0, stack | 3.09 | 30.85 | 39.81 | 1375 | 5552 | 1.109e5 | 23.177 |
| 100 | t_ic_l, stack + `ERLAB` | 3.53 | 30.81 | 39.78 | 924 | 5041 | 1.112e5 | 23.193 |
| 200 | t_cv200l, stack + `ERLAB` | 3.05 | 30.43 | 35.17 | 1574 | 6132 | 1.126e5 | 23.239 |
| 50 | t_cv50e, `ECOMPAT` | 4.61 | 43.60 | 61.94 | 1938 | 4859 | 5.99e4 | 22.183 |
| 100 | t_en_ic1, `ECOMPAT` | 2.72 | 35.36 | 45.58 | 1272 | 5104 | 8.09e4 | 22.910 |
| 200 | t_cv200e, `ECOMPAT` | 2.40 | 30.84 | 37.19 | 1106 | 5458 | 9.88e4 | 23.131 |
| 50 | t_cv50g2, `ECOMPAT` + `ECGRAV=1` | 5.51 | 32.82 | 51.73 | 1331 | 4834 | 1.011e5 | 22.670 |
| 100 | t_en_ic3, `ECOMPAT` + `ECGRAV=1` | 2.57 | 31.57 | 39.66 | 1556 | 5758 | 1.132e5 | 23.120 |
| 100 | t_ic_lf, stack + `ERLAB` + `GFV` | 2.07 | 27.83 | 36.53 | 979 | 5313 | 9.00e4 | 22.861 |
| 200 | t_cv200f, stack + `ERLAB` + `GFV` | 3.01 | 30.17 | 34.35 | 1252 | 5820 | 1.146e5 | 23.223 |

The difference ECOMPAT − stack at M\* 0.10 / 0.15 / 0.20 M☉ (stack's member at 100 zones: t_en_ic0, which has all three masses):

| zones | R\* | T_c | log n_c [dex] |
|---|---|---|---|
| 50 | +54 % / +48 % / +40 % | −41 % / −41 % / −41 % | −0.53 / −0.53 / −0.52 |
| 100 | +18 % / +16 % / +15 % | −27 % / −27 % / −27 % | −0.26 / −0.26 / −0.27 |
| 200 | +4.6 % / +5.6 % / +1.3 % | −12 % / −12 % / −12 % | −0.11 / −0.11 / −0.11 |

Radius of the mass shells at M\* ≈ 0.20 M☉ (`pmass.jl`, the snapshot nearest in M\*, given in brackets), r [R☉] at enclosed
mass 0.005 / 0.02 / 0.06 / 0.13 M☉, and the number of zones inside R\*:

| zones | run (M\* of the snapshot) | r(0.005) | r(0.02) | r(0.06) | r(0.13) | zones in the star |
|---|---|---|---|---|---|---|
| 50 | t_cv50s (0.189) | 0.703 | 2.015 | 5.668 | 13.71 | 8 |
| 100 | t_en_ic0 (0.197) | 0.487 | 1.443 | 4.335 | 11.38 | 14 |
| 100 | t_ic_l (0.200) | 0.478 | 1.412 | 4.249 | 11.07 | 15 |
| 200 | t_cv200l (0.201) | 0.425 | 1.334 | 4.073 | 10.71 | 30 |
| 50 | t_cv50e (0.168) | 1.137 | 3.116 | 8.360 | 21.45 | 9 |
| 100 | t_en_ic1 (0.206) | 0.637 | 1.819 | 5.187 | 12.72 | 15 |
| 200 | t_cv200e (0.202) | 0.480 | 1.448 | 4.312 | 11.17 | 29 |
| 50 | t_cv50g2, `ECOMPAT` + `ECGRAV=1` (0.184) | 0.697 | 1.942 | 5.384 | 13.55 | 8 |
| 100 | t_ic_lf, `GFV` (0.204) | 0.535 | 1.437 | 4.064 | 10.25 | 14 |
| 200 | t_cv200f, `GFV` (0.208) | 0.421 | 1.315 | 3.970 | 10.25 | 30 |

* **The two stars approach each other, and towards the stack's.** The difference in the central density falls from 0.53 dex
  (50 zones) to 0.26 dex (100) and 0.11 dex (200), in T_c from −41 % to −27 % and −12 %, in R\* from +40…54 % to +15…18 % and
  +1…6 %. Between 100 and 200 zones the stack's star moves by +0.05…+0.06 dex in n_c, +1.3…+4 % in T_c and 0…−1.3 % in R\*
  (t_ic_l → t_cv200l, same rows; M\* 0.10, 0.15, 0.20); ECOMPAT's moves by +0.22 dex, +22…25 % and −9…−13 %. Both move the same way, to a more
  compact star, and ECOMPAT's by four to five times as much. A Richardson estimate from the three zone numbers puts the limit
  of log n_c at M\* 0.20 at 23.24–23.25 for the stack's family (increments +0.49, +0.05) and at 23.23 for ECOMPAT (increments
  +0.73, +0.22, ratio 3.3); the radius of the 0.02 M☉ shell at 1.32 R☉ and 1.30 R☉ (rough for ECOMPAT: its 50-zone snapshot is at M\* 0.168).
  **Against the 200-zone stack star at M\* 0.20, the 100-zone stack star (t_ic_l) is within 0.05 dex in n_c, 1.2 % in T_c,
  1.3 % in R\* and 3–12 % in the radii of its mass shells (the innermost 0.005 M☉ worst); the 100-zone ECOMPAT star is 0.33
  dex, 28 %, 16 % and 19–50 % away.**
* **The difference between the stars converges more slowly than G5**: 0.53 → 0.26 → 0.11 dex is order 1.0 and 1.3, and T_c's
  −41 → −27 → −12 % order 0.6 and 1.2, where G5 over the release rate has order 1.4 and 1.8. ECOMPAT's T_c is not in a
  geometric sequence yet (increments +2.1e4 and +1.8e4 K), so no limit is quoted for it.
* **The photosphere is not converged at 200 zones in either family**: r_ph 48.8 / 39.8 / 35.2 R☉ and T_eff 4470 / 5550 / 6130 K
  (stack), 61.9 / 45.6 / 37.2 R☉ and 4860 / 5100 / 5460 K (ECOMPAT). (L_ph at one mass is a snapshot value of a luminosity that
  varies by a factor of two within the window; it is given for completeness.) The photosphere lies in the infall above the
  shock, 1.15–1.55 R\*; its radius follows the density of the infall there and the shock's width in zones.
* **The accretion history is not converged either**: the core forms at 2.586e5 / 2.428e5 / 2.384e5 yr with the stack's rows
  (ECOMPAT: 2.598e5 / 2.411e5 / 2.390e5 yr) and the accretion rate in the window is 4.1e-2 / 6.4e-2 / 7.7e-2 M☉/yr (ECOMPAT
  4.0e-2 / 5.9e-2 / 7.5e-2). The comparison at equal mass is therefore a comparison of stars with accretion rates that differ
  by up to a factor 1.9 between 50 and 200 zones; within one zone number the two families differ by 3–10 %.
* **The finite-volume weight at 200 zones gives the node weight's star to 0.02–0.03 dex in n_c.** `GFV` − node at
  M\* 0.10 / 0.15 / 0.20 M☉: −0.033 / −0.023 / −0.016 dex in n_c, −0.4 / +0.9 / +1.8 % in T_c, −6.3 / −7.6 / −0.9 % in R\*
  at 200 zones (t_cv200f − t_cv200l); −0.35 / −0.35 / −0.33 dex, −21 / −21 / −19 %, −16 / −13 / −10 % at 100 zones (t_ic_lf −
  t_ic_l). So the less dense centre of the 100-zone finite-volume star (§19.5) is an error of that weight on the coarse
  central mesh (ratio of successive radii 2.2–1.6 at interfaces 2–7 at 100 zones, 1.5–1.3 at 200), and it converges: the
  finite-volume star moves by +0.36 dex in n_c between 100 and 200 zones, the node-weight star by +0.05. The outer mass
  shells go the other way: the finite-volume star has the 0.13 M☉ shell at 10.25 R☉ at 100 and at 200 zones (0.06 M☉: 4.06
  and 3.97 R☉), and the node-weight star moves towards these values (11.07 → 10.71 and 4.25 → 4.07 R☉). The inner shell
  (0.005 M☉) is at 0.42 R☉ with both weights at 200 zones.

### 21.6 Steps, failures, wall time

At the first perf.csv row (every 100 steps) with M\* ≥ 0.20 M☉ (driver's M\*), counted from the IC: step, wall time, failed TR /
BDF stages; and the Newton-failure lines of the whole run. 4 threads, except the three 100-zone runs of §16 and §20 (8
threads). The analysis jobs of this section ran on core 31, which t_ic_l was using; its wall time is an upper bound.

| zones | run | step | wall [s] | failed TR / BDF | Newton-failure lines (run to step) |
|---|---|---|---|---|---|
| 50 | t_cv50s, stack | 1300 | 2076 | 109 / 29 | 163 (1700) |
| 50 | t_cv50e, `ECOMPAT` | 600 | 932 | 37 / 1 | 52 (700) |
| 50 | t_cv50g2, `ECOMPAT` + `ECGRAV=1` | 1200 | 1533 | 83 / 4 | 124 (1400) |
| 100 | t_en_ic0, stack (8 threads) | 1300 | 2122 | 66 / 16 | 171 (2100) |
| 100 | t_ic_l, stack + `ERLAB` | 1300 | 2927 | 53 / 7 | 149 (2150) |
| 100 | t_en_ic1, `ECOMPAT` (8 threads) | 800 | 1400 | 43 / 7 | 136 (1478) |
| 100 | t_en_ic3, `ECOMPAT` + `ECGRAV=1` (8 threads) | 1200 | 1899 | 56 / 5 | 177 (2100) |
| 100 | t_ic_lf, stack + `ERLAB` + `GFV` | 800 | 2089 | 29 / 10 | 100 (1800) |
| 200 | t_cv200s, stack | — | — | — | does not pass the collapse (§21.1) |
| 200 | t_cv200l, stack + `ERLAB` | 1000 | 4901 | 53 / 10 | 155 (1830) |
| 200 | t_cv200e, `ECOMPAT` | 900 | 4535 | 88 / 9 | 133 (1200) |
| 200 | t_cv200f, stack + `ERLAB` + `GFV` | 1000 | 5066 | 99 / 20 | 180 (1720) |
| 200 | t_cv200g, `ECOMPAT` + `ECGRAV=1` | not reached (M\* 0.166 at step 800) | 4090 (to step 800) | 53 / 7 (to step 800) | 60 (800) |
| 100 | t_ic_lv2, stack + `ERLAB` + `GVAR` (§19.5) | 1200 | 2707 | 60 / 6 | 92 (1580) |

The cost per step after core formation is 1.0–1.1 s at 50 zones, 1.9–2.1 s at 100 and 4.2–4.6 s at 200 (4 threads); the number of steps to a given
mass does not grow with the zone number. Every run from the IC passes through stretches of 100–200 steps with Δt of 1e3–1e4 s
(limiter: T or E_r of a zone at or just above the photosphere) after M\* ≈ 0.20 M☉; they end by themselves in most cases
(t_ic_lf at M\* 0.25 and 0.30, t_cv200l at 0.22 and 0.25, t_ic_l at 0.20 and 0.29; t_ic_l was in one at M\* 0.32 from step 1850 to 2150, where I stopped it), and I stopped t_cv200e in one at M\*
0.226 and t_cv200l in one at M\* 0.290 (Δt 5e2–2e3 s for 230 steps). They are not specific to a row or a weight and were not investigated.

**The 200-zone run with ECOMPAT + ECGRAV=1 (t_cv200g).** Started last; it reached M\* 0.15 at step 604 and then a stretch of
small steps (Δt 3e4 s; M\* 0.162 → 0.166 in steps 700–800), where I stopped it. Its window is therefore M\* 0.10–0.15, and the
other two 200-zone runs are evaluated in the same window here (L☉):

| 200 zones, window M\* 0.10–0.15 | stack + `ERLAB` (t_cv200l) | `ECOMPAT` (t_cv200e) | `ECOMPAT` + `ECGRAV=1` (t_cv200g) |
|---|---|---|---|
| steps / time [yr] / release rate | 137 / 0.611 / 15500 | 84 / 0.628 / 14521 | 78 / 0.591 / 15847 |
| G5 (fraction of the release rate) | 1154 (0.0744) | 1021 (0.0703) | 1188 (0.0750) |
| of it star / outer layers / infall | 0.0272 / 0.0326 / 0.0144 | 0.0252 / 0.0310 / 0.0138 | 0.0279 / 0.0324 / 0.0145 |
| residual of the conservative total row, whole mesh | 256 | 1020 | 233 |
| ledger defect, W (i) / (ii) / (iv) | −396 / −349 / −193 | +411 / +450 / +585 | −306 / −262 / −94 |
| as a fraction of the release rate, W (i) / (iv) | −0.026 / −0.012 | +0.028 / +0.040 | −0.019 / −0.006 |
| dW/dt − total row's gravity term, W (i) / (iv) | −652 / −449 | −610 / −435 | −540 / −328 |
| max u² r/(2GM) in the infall, time mean | 0.960 | 0.956 | 0.956 |
| star at M\* 0.05: R\* [R☉] / T_c [K] / log n_c | 8.99 / 1.015e5 / 23.169 | 9.50 / 8.86e4 / 23.058 | 9.01 / 1.032e5 / 23.168 |
| star at M\* 0.10 | 17.23 / 1.075e5 / 23.208 | 18.02 / 9.45e4 / 23.102 | 17.14 / 1.088e5 / 23.203 |
| star at M\* 0.15 | 24.67 / 1.107e5 / 23.228 | 26.05 / 9.71e4 / 23.120 | 24.80 / 1.117e5 / 23.220 |

At 200 zones ECOMPAT with the G5 sink in the star's zones is the stack's star to 0.000 / −0.005 / −0.008 dex in n_c,
+1.6 / +1.2 / +0.9 % in T_c and +0.2 / −0.5 / +0.5 % in R\* at M\* 0.05 / 0.10 / 0.15 (ECOMPAT alone: −0.11 dex, −12…−13 %,
+5…+6 %), with the stack's row residual (233 against 256 L☉; ECOMPAT alone 1020 L☉) and a ledger of −1.9 / −0.6 % of the
release rate (W (i) / W (iv); the stack + `ERLAB` −2.6 / −1.2 %, ECOMPAT +2.8 / +4.0 %). G5 itself is the same in the three
runs (0.070–0.075): the sink does not change the gravity power, it changes where the difference goes. In this shorter
window the 200-zone values of the other two runs are those of the window 0.10–0.20 to within 0.005 in G5 / release rate
(0.0744 against 0.070 and 0.0703 against 0.066).

### 21.7 Answers to the three questions of item 4

1. **At what order does each gravity power converge?** In the zone number N, measured as log₂ of the ratio of the fractions
   of the release rate at N and 2N, window M\* 0.10–0.20 M☉:

   | | 50 → 100 | 100 → 200 |
   |---|---|---|
   | G5, whole mesh (stack family / ECOMPAT / ECOMPAT + ECGRAV=1) | 1.63 / 1.42 / 1.81 | 1.90 / 1.83 / 1.7 (200-zone value from the window M\* 0.10–0.15, §21.6) |
   | density-jump part, outer layers of the star | 1.82 / 1.39 | **1.92 / 1.90** |
   | geometric part, infall | 1.65 / 1.48 | **1.37 / 1.43** |
   | geometric part, star's interior | 0.82 / 0.90 | 2.10 / 1.86 |
   | density jump against the finite-volume weight, whole mesh (= G5 under `GFV`) | 1.68 | 2.05 (measured with `GFV`, t_ic_lf → t_cv200f: 2.01) |
   | mass-flux part of G5 (u_rel; ≈ G5x of §14.3): outer layers / infall, stack family | | 1.95 / 1.72 |
   | mesh-motion part of G5 (ṙ), star's interior (§21.3) | changes sign | **0.5 / 0.3: does not converge** |

   The density-jump part is second order once the outer layers have more than a handful of interfaces (8 at 100 zones, 18 at
   200). The geometric part of the infall is the slowest (1.4); it is 2.4 % of the release rate at 200 zones and by then a
   third of G5. The 50-zone runs are not in the asymptotic range for the star's interior. In the star's interior the split
   of G5 into its mass-flux and mesh-motion parts does not converge (each ±3–8 % of the release rate at 100 and 200 zones,
   because the mesh velocity there doubles with the zone number); their sum, which is what the rows differ by, does (order
   2). The mesh-motion part is the difference between dW/dt and the total row's gravity term (§21.3).
2. **Do the two stars approach each other, and towards which one?** Yes, towards the stack's. ECOMPAT − stack at equal mass:
   −0.53 / −0.26 / −0.11 dex in n_c, −41 / −27 / −12 % in T_c, +40…54 / +15…18 / +1…6 % in R\* at 50 / 100 / 200 zones. From
   100 to 200 zones the stack's star moves by +0.05…+0.06 dex in n_c and ECOMPAT's by +0.22 dex, in the same direction (more
   compact). The estimates of the limit from three zone numbers agree for the two families (log n_c 23.24–23.25 and 23.23 at
   M\* 0.20). ECOMPAT with the G5 sink in the star's zones (ECGRAV=1) is the stack's star at 50 zones as at 100 (log n_c 22.670
   against 22.705, T_c 1.011e5 against 1.007e5 K, R\* 32.8 against 31.2 R☉ at M\* 0.20; ECOMPAT alone 22.183, 5.99e4 K, 43.6
   R☉), and at 200 zones (t_cv200g, at M\* 0.15: log n_c 23.220 against 23.228, T_c 1.117e5 against 1.107e5 K, R\* 24.8
   against 24.7 R☉; ECOMPAT alone 23.120, 9.71e4 K, 26.1 R☉; §21.6). So the difference between the stars is G5 in the star's zones at every zone number measured, it falls as G5
   falls, and what G5 gives the ECOMPAT star is not paid for by the binding energy (§20.2): the converged star is the
   stack's, and at 100 zones the stack's rows are within 0.05 dex, 1.2 % and 1.3 % of the 200-zone values of n_c, T_c and R\*.
3. **The finite-volume weight at 100 and 200 zones.** Measured at both (t_ic_lf, t_cv200f; §19.5, §21.2–§21.6). With `GFV`
   G5 is the density jump alone: 0.302 → 0.075 of the release rate (order 2.0; node weight 0.261 → 0.070), the infall's part
   0.011 → −0.0005, the row residual summed over the mesh 63 → −23 L☉ (node weight 459 → 198 L☉), the maximum of u² r/(2GM)
   0.980 → 0.945 (node weight 1.018 → 0.959). Its star is 0.33–0.35 dex less dense at the centre than the node weight's at
   100 zones and 0.02–0.03 dex at 200 zones (T_c within 2 %): the weight works at 200 zones and not at 100. The
   binding-energy ledger is not better with it at either zone number (−4.0 / +1.2 % at 100, −4.7 / −3.6 % at 200; W (i) /
   W (iv)), because what the ledger has left is the term of the star's interior (§21.3), which no weight changes.

### 21.8 What was not run, and what is open

* The 200-zone ECOMPAT + ECGRAV=1 run (t_cv200g) beyond M\* 0.166: its numbers are for the window M\* 0.10–0.15 (§21.6), not
  for 0.10–0.20.
* The 200-zone ECOMPAT run beyond M\* 0.226 (stopped in a stretch of small steps) and the 200-zone stack + `ERLAB` run beyond
  M\* 0.29 (stopped likewise), and the 200-zone `GFV` run beyond M\* 0.30 (stopped likewise): the 200-zone numbers are for the window M\* 0.10–0.20 (and 0.20–0.28 for t_cv200l), not for 0.10–0.335.
* The 200-zone run with the stack's comoving radiation row (t_cv200s fails before the core forms; the cause was not
  investigated), and a 50-zone run with the stack's gas rows and `ERLAB` (for the 50-zone member of the stack's family the
  cloud's residual is given separately in §21.3 instead).
* 141 zones (run directories prepared, not run), which would give a second estimate of each order.
* `GVAR` from the IC beyond M\* 0.25 (t_ic_lv2, stopped by me; §19.5 (4)) and at 200 zones, ECOMPAT + `GVAR` from the IC (t_ic_ev, prepared), and ECGRAV=1 together with `GFV`.
* **Open: the interface density of the transport on the stretched mesh** (§21.3). dW/dt and the total row's gravity term
  differ by −0.8e3 L☉ at 100 and at 200 zones in every window from the IC. The static test identifies the term (not a drift
  of W; the advected interface density is below the density of the smooth mass distribution by 2–9 % at 200 zones and 7–29 %
  at 100 in the star's interior), and the same evaluation attributes nine tenths of G5 to it at 100 and 200 zones. Not done:
  a run with a geometry-consistent, upwinded and limited interface density in `advect!` (every row uses it; what it should
  be at the accretion shock is a separate question); the dynamical test, a hydrostatic star with a prescribed mesh motion
  through it, which would also show whether the star takes the mass or pays with internal energy; and why GRIDEQ's mesh
  velocity in the star's interior doubles from 100 to 200 zones. As the rows stand the ledger's total cannot test a
  formulation to better than about 5 % of the release rate during the first accretion; its parts by region can.
* Open, not explained: the order 1.4 of the geometric part in the infall (the contribution of the upstream side of the
  AV-smeared shock was not separated from the smooth infall), and why the 200-zone run needs the lab-frame radiation row to
  pass the collapse.
* Open, not investigated: the stretches of small steps after M\* ≈ 0.20 M☉ in every run from the IC (§21.6).
* The photosphere (r_ph, T_eff) and the accretion history are not converged at 200 zones in either family (§21.5).

## 22. A set of energy-formulation principles against these rows (asked during round 3)

The principles (PI, 17 Nov): a conservative total-energy budget, 𝓔 = u_g + ½ρv² + E_r − ρ G m/r with flux
v(u_g + ½ρv² + p_g + q − ρ G m/r) + F_r, together with directly evolved gas and radiation energies and compatible discrete
transfers between the components; reconcile only a small remaining residual, in a resolved reservoir; never adjust kinetic or
binding energy; large residuals mean the discretisation has to be fixed. Assessment against what is on `energy-compat`
(numbers: 100 zones from the checkpoint unless stated; W i / ii / iv).

### 22.1 Per principle

| principle | status | where the discrete form departs |
|---|---|---|
| evolve u_g and E_r directly, no small component from a subtraction | **satisfied with `ECOMPAT`** (first law in every zone, no E − K, no flag); the stack: partly (dual energy: total row in the star's zones, where K ≪ e, first law in 74 of 100 zones) | — |
| the same discrete forces and time-centred velocities in momentum and energy | **satisfied with `ECOMPAT` for pressure, viscosity, momentum advection and the time discretisation** (§4): the row is the adjoint of the momentum row with the θ-level velocities; per-zone identity 0.000 L☉ in zones 1–87 | the stack's first-law zones do not return the kinetic energy the upwinded momentum flux (190 L☉) and the implicit time discretisation dissipate |
| the same mass fluxes in momentum, gravitational work and energy transport | **not satisfied, in any variant** | the momentum row's gravity is g_j × node mass ½(m_j−1 + m_j); the mass flux and the total row's gravity term use the advected interface density ρᵇ. The two powers differ by G5: 2.7e3–3.2e3 L☉ from the checkpoint, 26 % of the release rate at 100 zones and 7 % at 200 from the IC (§21.2). The row that is consistent with the actual mass flux (`ECGRAV=3`, K + U + W + E_r to 1.4 L☉) fails in the cold infall (§14.5): it places the difference in a component that is 1e-3–1e-4 of the kinetic energy there |
| artificial viscosity as kinetic-to-thermal transfer | **satisfied with `ECOMPAT`**: the heating is q r̄^{3μ/2} Δ(u r^{−μ/2}), the discrete work of the viscous force of the momentum row (§4.2) | the stack books q (du/dr) ΔV; the difference V4 is 0.5–1.3 L☉ |
| gas–radiation exchange equal and opposite | **satisfied** (`GASROW`, in the stack and in `ECOMPAT`): the same expression, with the same cap factor, in the gas row and the radiation row | — |
| consistent radiation-force work | **satisfied with the lab-frame radiation row** (`ERLAB`, or `ECOMPAT`): Σ u_j (χ̄F/c)_j ½ m_k in the radiation row is the momentum row's term times u (§3.3) | the stack's comoving row is the adjoint of the moment equation, not of the flux tie that replaces it in zones 51–100: 3296 L☉ created there, 1.9 L☉ with `ERLAB` (§18) |
| pressure–gravity balance for near-hydrostatic flow | partly | the node weight is exact for ρ ∝ r⁻² on any mesh and within 2.5 % of the smooth distribution's weight in the star at 100 zones (§21.3); none of the three weights is well-balanced in the sense of Käppeli & Mishra; `GFV`, `GVAR` tested and not adopted at 100 zones (§19) |
| record corrections and component discrepancies | **satisfied as a diagnostic** (`EDIAG=2`: per-zone residual of the conservative total row, the identity of §4.6, and the binding-energy ledger from the states); nothing is corrected | — |
| reconcile the small residual in a resolved reservoir | not in the code | §22.3 |

Other provisional switches named in the question: `GEQGP` (mesh weight on density jump × gravity power) did not reduce G5
(§14.6); `TIEX` is not general (§15.4); `ERESTOL` is a stopping test inside the Newton iteration, not a correction (§17).

### 22.2 The binding energy of the law and of the ledger

The same budget. The ledger's W is −Σ G M dm / r with the enclosed mass, no factor ½ (§14.2: (i) at the interfaces, (ii) at the
zone centres, (iv) the exact integral for piecewise-constant density, which is the energy whose gradient `GVAR` uses), plus
the potential energy in the external mass (−1.2 L☉). The centre is regular (no central mass, no accreting inner
boundary), the outer boundary's fluxes are the total row's. The three definitions differ by ≤ 0.3e3 L☉ at 100 zones from the
checkpoint and by 38–145 L☉ at 200 zones, so a different gravitational energy is not what the defects are. Of the 700 / 665 /
957 L☉ of `ECOMPAT` + `ECGRAV=1`: 750 L☉ is the conservative row's residual in the shock and infall zones (the momentum
row's gravity where no sink is applied; `GFV` removes it, §19), 152 L☉ in the star, 3 L☉ in the cloud, and −204 / −240 / +52 L☉
the difference between dW/dt and the flux-form gravity term, which contains the spread of the definitions. The moving mesh
is where the stack's total row departs from the law: the binding energy of a zone changes by (mesh velocity) × g × (the zone's
mass) plus (mass flux) × (potential difference), and the row books the first part with the advected interface density
instead: −0.8e3 L☉, 5 % of the release rate, at 100 and at 200 zones from the IC. §21.3 shows that this is not a drift of W
but the interface density of the transport, low by 2–9 % (200 zones) and 7–29 % (100 zones) on the star's stretched mesh.

### 22.3 A reconciliation step

Not implemented, and not recommended on the present rows.

* For internal + kinetic + radiation energy there is nothing to reconcile: the identity closes per zone to round-off, and what
  is left is the Newton iteration's residual, 4e-6 of L_em time-integrated, which `ERESTOL` reduces inside the iteration (48 of
  1000 stages above 1 L☉, largest 11 L☉; §17).
* Against the binding energy the residual is not small: 0.5e3–4.1e3 L☉ at 100 zones from the checkpoint (8–60 % of L_em), 2–5 %
  of the release rate at 200 zones from the IC. By the principles' own rule this is a discretisation to fix, not a residual to
  place. Placing it in the gas has been measured, as a term inside the implicit solve: in the star's zones (`ECGRAV=1`,
  2189 L☉) it gives the stack's star and needs the dual-energy flag; in every zone (`ECGRAV=2`, `3`) it fails in the cold
  infall. A post-step correction of the same size would in addition be an input the Jacobian does not see, first order in
  time, and in the optically thin zones (51–100 here, and the infall) E_r is not thermal, so a common-temperature split has
  no meaning there.
* Cost if it were wanted for the star's zones only: one pass over the zones per step with the per-zone total-row residual
  that `EDIAG=2` already forms and δT = −δ_k / (m c_v + 4 a T³ ΔV)_k; negligible. It would be `ECGRAV=1` applied after the
  step instead of inside it.
* So: record, do not correct. The recorded quantities are the per-region residual of the conservative total row and the
  ledger against W; at 100 zones their floor is the spread of the definitions of W, 0.3e3 L☉ (0.7e3 L☉ during the first
  accretion from the IC).

### 22.4 Does any principle change the recommendation?

No. The configuration the measurements favour (§18, §20.4, §21) — the stack's gas rows with the lab-frame radiation row — is the
one in which each principle that can be met on this mesh is met where it matters: radiation-force work consistent, exchange
equal and opposite, total row where the internal energy is resolved, first law where it is not. `ECOMPAT` meets more of
the principles (no subtraction, no flag, viscous and advective dissipation returned) and is the better formulation of
internal + kinetic + radiation energy, but with the present gravity terms it gives the star 2.0e3 L☉ the binding energy
does not lose, and its star is the one further from the 200-zone star. What the principles add is the direction: "large
residuals require fixing the discretisation" applies to G5, and §21.3 names the discretisation — the interface density of the
mass flux on the stretched moving mesh. With that density consistent, the momentum row's gravity, the mass flux and the
energy transport would be the same to 0.7 % (200 zones) and 2.5 % (100 zones) of the release rate on the present states, and
the question which energy row to prefer would be nearly decoupled from the gravity power. That is to be tested by a run.

## 23. The internal-energy equation alone, compatible with the momentum row (branch `iecons`, ETOTAL=0)

Decision (T. Abel, 18 Nov): one energy variable, the gas internal energy (plus the radiation energy), no dual-energy switching and
no total-energy row; the total energy is only monitored (`ETOTMON`, `TX.etot_monitor`). For that the internal-energy row must be
the exact discrete adjoint of the momentum row, so that Σ_k (gas row + radiation row)_k + Σ_j u_j (momentum row)_j is the
conservative total-energy row sum identically. §4 derives this row (`ECOMPAT`); this section is the term-by-term check against
the code (`julia/TitanResidFast.jl`, `residual!`; momentum row at "GAS MOMENTUM", energy row at "GAS ENERGY", branch `if ECOMPAT[]`).

### 23.1 Discrete kinetic-energy equation of the momentum row, and what the gas row returns

Momentum row of interface j (mass m_j = ½(m_{j-1} + m_j) of the two adjoining zones, θ-level velocity u):

    [m_j^n u_j^n − m_j^o u_j^o]/Δt − ½[(ṁ_j + ṁ_{j+1}) û_j − (ṁ_{j-1} + ṁ_j) û_{j-1}] + r_j^μ (p_j − p_{j-1}) + g_j m_j + r_j^{−μ/2}(r̄_j^{3μ/2} q_j − r̄_{j-1}^{3μ/2} q_{j-1}) − (χ̄F/c) m_j = 0.

Multiplying by u_j and distributing half to each adjoining zone gives, per zone k (faces k, k+1), the following pieces. Every
one that is not already the first law's own term is a product of velocity *differences*, never K_new − K_old:

| piece of u·(momentum row) | term of the gas row (all functions of the Newton iterate, type T) | size on the checkpoint (docs §5, ledger, L☉ time-weighted) |
|---|---|---|
| pressure work u_j r_j^μ (p_k − p_{k-1}) | p_k Δ(r^μ u): `pg*(rmu[k+1]*u[k+1]-rmu[k]*u[k])` (already so in the old row) | 0 |
| viscous work u_j r_j^{−μ/2} Δ(r̄^{3μ/2} q) | q_k r̄_k^{3μ/2}(u_{k+1} r_{k+1}^{−μ/2} − u_k r_k^{−μ/2}): `qf*(vhi-vlo)`; the old row had q (du/dr) ΔV | 1.5 (V4) |
| time: u_j Δ(½ m u)/Δt = Δ(¼ m Σu²)/Δt + (2θ−1)¼ m̂ Σδ²/Δt − (θ²−θ+½)¼(m_n−m_o)Σδ²/Δt | `(-(2θ-1)*0.125*(mn_+mo_) + (θ^2-θ+0.5)*0.25*(mn_-mo_))*du2/dt` (δ = u^n−u^o) | θ-step and TR stages: not measured separately; TR-BDF2: below |
| BDF stage: K(y*) ≠ A K_γ − B K_n | `KSTAR[k]/dt` (Trbdf2.bdf_old!) | 0.13 (BDF stages) |
| momentum advection û_j ≠ ½(u_j + u_{j+1}) (upwinded van Leer value): u_j ṁ-flux of the momentum = ṁ ½u² flux + K3a | `0.5(ṁ_k+ṁ_k+1)(u_k+1−u_k)(½(u_k+u_k+1)−ub[k])` | 189 together with the next line (62 in zones 1–26, 127 in 27–50) |
| flux of the sub-zone kinetic energy ⅛m(Δu)², so that K3a is thermalised where the mass arrives and is sign-definite (§4.4) | `−(ksflux(k+1)−ksflux(k))`, G_j = ṁ_j (u_up+1−u_up)²/8 | (telescopes) |
| gravity u_j g_j m_j | none in the row. The total-energy row (ETOT_GRAV) books ṁ×Δψ with the advected interface density; the two differ by G5 (§14), 1.0e3–1.9e3 L☉ in the star, 0.7e3–0.85e3 in zones 27–50. `ECGRAV=2/3` is the optional term that makes the row equal to the monitor's | G5 |
| radiation force u_j (χ̄F/c)_j m_j | the lab-frame radiation row (comoving row + u × radiation momentum equation: `wk`, §3.3); gas–radiation exchange equal and opposite in both rows | 0 (zones 1–87), ≤ 4 L☉ round-off in 88–100 |
| energy of the old level | e_old = (e_g + E/ρ)_old − E_old/ρ_old of the previous solve (`c.egro − c.ero/c.do_`), `c.ego` without radiation | — |

Added in this branch: (i) the compatible row also for `lrad = 0` (pure hydrodynamics: the row used `egb`, `c.ego`; the exchange and
Compton terms are absent), selected by `ETOTAL=0` in `phase6_run.jl` and `popiii_m2c.jl` (`TX.ECOMPAT[] = TX.GASROW[] = true`);
(ii) the upwind branch of `ksflux` goes through `decide(:up, 20, …)` like every other upwind choice, so that the finite-difference
J·v replays the branch the Jacobian was built with — before this the T row of the shock zones disagreed with its band
row by 8e-3 (zone 29) in `NKCHECK` (t_ie20 / t_ie22), after it by 4e-6 (t_ie23): the flux is continuous at ṁ = 0 but has a kink there,
and the front zones have ṁ ≈ 0.

### 23.2 Heat inputs of the three-level atom (audit A5)

The chemistry heat cq (photon energy given to the gas; inside `cpL` as heating of the gas in coupled zones) and the cycle heat yq were booked once per step
from the end state of the previous step (`RAYQ`, `GASQ`). Now the hydro rows take them as functions of the Newton iterate in a coupled zone
(`W.cpQ`, `W.cpY`, from `cpl_g`, already computed in every residual) and from the booking of the last chemistry otherwise (`CQB`, `YQB`, set by
`chem_core!` in every zone, unconditionally; `cqheat`, `yqheat` choose by `cplz(k)`, the flag the residual uses for the whole attempt). Each attempt
therefore applies the heat of each zone from exactly one source, whatever the zone's status in the neighbouring steps (the first version skipped the booking
in the post pass of coupled zones, which lost the heat of a zone that was decoupled in the next step and doubled that of a zone that became coupled).
Radiation row +cq+yq, gas rows −yq, total/combined rows +cq. `RAYQ` keeps the ray source, `GASQ` the Lyman-Werner heating. The decoupling message
(`CPL_BAD`) prints the booked cq, yq the zone now uses. Old checkpoints carry cq, yq inside `RAYQ`/`GASQ` and no `CQB`, `YQB`: in the first step after such a restart the
coupled zones count them twice.

### 23.3 Gates (18 Nov; runs in `/zpool/nvme/data/tabel_scratch/titan_popiii/runs/`)

**(a) ETOTMON, `phase6_run.jl pure NAME` with `ETOTAL=0 ETOTMON=1` (`t_ie_reg/ec2_*.log`, `ec2s_SUBC.log`, per-solve CSV `ec2_*.csv`)**; G = total-energy residual / energy throughput.

| test | solves | old ETOTAL=0 (median / max) | compatible rows (median / max) | ETOTAL=1 (ETASW=0) |
|---|---|---|---|---|
| SOD | 433 | 1.4e-6 / 3.9e-4 (397 solves) | **1.6e-16 / 6.5e-13** (worst zone 99, solve 73, 2e-2 of that zone's own throughput) | 1.8e-16 / 2.6e-12 |
| SUBC | 722 | 1.6e-5 / 8.6e-2 | 2.9e-11 / 1.5e-1 (solves 54–89, 165–205: dt 1e-25…1e-11 s start-up of the benchmark, worst zones 2, 55, 71, 88) | 1.2e-14 / 4.5e-2 |
| SUBC, `stol=1e-6`, 999 solves | 999 | — | 7.6e-10 / 3.2e-9 | — |
| SUPC | 999 | 1.0e-5 / 3.7e-1 (713 solves) | 2.3e-9 / 3.1e-1 (worst zone 11, solve 153) | 4.1e-9 / 9.2e-2 |

The SOD residual is round-off. In the radiative shocks the median is 3–5 orders better than before, equal to the total-energy row's (SUPC) or
within 3 orders of it (SUBC); the maxima are the start-up solves, where the converged rows are not converged to the monitor's precision: tightening the
Newton stopping test (`stol` 1e-6) takes SUBC's maximum to 3e-9, i.e. what is left is the solver tolerance, not a term of the row.

**(b) Noh problem, spherical, γ = 5/3, 200 zones, Mach 85, e_kin/e_th = 4009 (`noh_test.jl`; `t_ie_noh/*.log`)**; T_RH = 4009 K, ρ_RH = 64 ρ₀, shock at |u₀|t/3.

| rows | shock radius / exact | T/T_RH in 0.3–0.8 r_s (median, min, max) | ρ/ρ_RH (median, min, max) |
|---|---|---|---|
| ETOTAL=0, compatible | 0.978 | 1.020, 0.914, 1.090 | 1.019, 0.906, 1.158 |
| ETOTAL=1, ETASW=0.1 | 0.853 | 0.904, 0.757, 0.995 | 1.314, 1.117, 1.699 |
| ETOTAL=1, ETASW=0.6 | 0.755 | 0.767, 0.002, 0.919 | 1.783, 0.379, 2.573 |
| ETOTAL=1, ETASW=0 | fails to converge at t = 1504 s (step 686) | | |
| TITAN's rows (unset) | 0.681 | 0.653, 0.002, 0.791 | 2.236, 0.368, 3.277 |

**(c) NKCHECK=2 on the Pop III checkpoint (`t_ie23`, with `ksflux` branch replayed, A5 included); `|Jv − J_b w|/Σ|J_b w|`, maximum over the zones (zones 28–34):** T 4.3e-6 (1.6e-6), u 1.3e-4 (3.2e-6), E_r 1.5e-5, ρ 4.6e-6, m 3.6e-6. Before the fix (t_ie20, t_ie22) T was 8.0e-3 (zone 29; the m column, 1 %).
The old ETOTAL=0 row had 4.3e-6 / 7.3e-5 (t_ie10). F_r (8.5e-2) and x_H₂ are as before.

**(d) 300 steps from the checkpoint, ETOTMON=1 (all with `NKRYLOV=1 DTERR=0.03`):**

| run | rows | t reached (timen − 7.50846e12) | rejections (Newton failures) | min dt after failure / traced | limiters | ETOTMON |G(last bdf)| median / max, worst zones |
|---|---|---|---|---|---|---|
| t_ie04 | ETOTAL=1 + ERLAB | 4.64e8 s | 12 (12) | 142 s | T@30 | — |
| t_ie05 | old ETOTAL=0 + ERLAB | 5.26e8 s | 15 (14) | 1723 s | rho@28, rho@29 | 4.9e-4 / 9.0e-4 (window max 1.3e-3); z29, z34 |
| t_ie26 | compatible, A5 off (experiment, HEAD~1 code) | 5.29e8 s | 23 (22) | 1157 s | T@30, rho@27 | 6.9e-4 / 1.6e-3; z29, z28 |
| t_ie24 | compatible + A5 (HEAD) | 0.94e8 s | 33 (29) | 1215 s | E_r@21, T@30 | 6.6e-4 / 8.9e-4; z28, z34 |
| t_ie25 | compatible + A5 + `ECGRAV=2` | −0.16e8 s (7.508299e12) | 45 (25 Newton + 20 error) | 460 s | T@25, T@30, T@35 | **1.3e-7 / 2.5e-6**; z28, z33, z23 |

(`t_ie21`, compatible without A5 and without the J·v fix, was stopped by me at step 65757 after one rejection and dt 2.9e6 s: the first 60 steps were
the best of all variants; it was not continued because the A5 edit was to go in.)

What this says. The internal-energy row is compatible with the momentum row (SOD to 1e-16, the cold-infall Noh test within 2 %), but on the Pop III
checkpoint the monitor's G stays at the old 5e-4…7e-4 because the monitor's total row books gravity as ṁ Δψ (ETOT_GRAV) and the momentum row's force does
not (G5, §14): with `ECGRAV=2` (the monitor's own gravity form put into the gas row) G falls to 1e-7, so the whole defect is the gravity term and the rest of the row is
exact. That variant, however, as §14.5 found, costs convergence (45 rejections). The 300-step runs do not meet "at least as good as ETOTAL=1 / old ETOTAL=0":
the compatible row alone has more Newton failures than both (22 against 12 and 14) at the same time covered as the old ETOTAL=0; A5 (iterate-dependent
cq, yq) makes it worse (29 failures, a fifth of the time covered). The dual-energy switches were therefore not removed.

### 23.4 Gravity of the momentum row equal to the total-energy row's work (ECOMPAT)

*Total-energy row* (`etot_row`, `ETOT_GRAV`): the work of gravity on the mass that crosses each interior face j, booked over half a zone on each side,

    W_j = r_j^μ ρᵇ_j u_j g_j (Δr_{j-1} + Δr_j) = g_j · [r_j^μ ρᵇ_j ½(r_{j+1} − r_{j-1})] · u_j,   Δr_k = ½(r_{k+1} − r_k),   g_j = 4πG·2 xm_j/r_j^μ + (1 − μ/2) g_0,

with ρᵇ the advected interface density of the continuity row (`W.db`) and u_j the θ-level velocity (the mesh velocity enters through the mass flux of the moving
mesh: `uf(j) = u_j` with `gridmoves()`, `u_rel,j` = `u_j` on a fixed mesh). The sum of the zone rows is the change of Σ m Φ (§14.2: it follows an explicit binding energy to 0.02–0.8 %).

*Momentum row* (before): force g_j m_j with the node mass m_j = ½(ρ_{j-1}ΔV_{j-1} + ρ_jΔV_j), so u_j g_j m_j ≠ W_j wherever ρᵇ_j r_j^μ ½(r_{j+1} − r_{j-1}) ≠ m_j: this is G5 (§14.4, §21.3).

*Change* (`gravcons`, `gravmass` in `TitanResidFast.jl`): with ECOMPAT the force at interface j is g_j m*_j, m*_j = r_j^μ ρᵇ_j ½(r_{j+1} − r_{j-1}), a function of the iterate (ρᵇ is
the van Leer value of the iterate's densities and of the mass flux). Then Σ_j u_j F_j equals Σ_j W_j face by face, the gas row needs no gravity term (it never had one), and
the total-energy row of the monitor is the sum of the solved rows and u × momentum rows up to Newton tolerance. Not used for a Lagrangian mass grid without grid motion
(`lgrid = 3`: no mass flux, the node mass is exact there; the dust-collapse and Bonnor–Ebert tests of `popiii_m1.jl` are unchanged). The MLT hydrostatic-gradient blend uses the same mass.
`ECGRAV` (all four variants, which put the difference G5 into the gas row) and its zone_diag term are removed; `GFV`/`GVAR` are alternatives to this force and are refused with ETOTAL=0.
The consistency does not make the force a better approximation of gravity: ρᵇ is the advected density, 7–29 % below the smooth distribution's on the star's stretched centre mesh
at 100 zones (§21.3); the star's structure changes accordingly (§23.5).

### 23.5 Gates of §23.2 and §23.4 (18 Nov)

* ETOTMON on SOD / SUBC / SUPC with ETOTAL=0 unchanged from §23.3 (`t_ie_reg/ec3_*.log`); Noh test unchanged (shock 0.978, T/T_RH 1.020, ρ/ρ_RH 1.019; `t_ie_noh/etot0_v2.log`).
* Self-gravitating test (`julia/selfgrav_energy_test.jl`, cold uniform 100 M☉ sphere on the fixed mesh, 100 zones, to 0.9 t_ff, 57 steps; `t_ie_noh/sg_etotal0.log`, `sg_etotal.log`): ETOTMON
  G with the consistent gravity force: median 2.0e-16, max 5.6e-12; with TITAN's rows (node-mass force, total row's gravity form in the monitor): median 4.9e-6, max 4.7e-5.
  On a fixed mesh the gravity gap is closed to round-off.
* NKCHECK=2 (`t_ie28`): T 4.3e-6 (zones 28–34: 1.6e-6), u 1.8e-4 (3.3e-6), E_r 4.3e-6, ρ 4.6e-6, r 4.2e-6 (was 3e-3), F_r 1.1e-1 as before.
* 300-step restart (`t_ie27`): **does not run**. 118 Newton failures, accepted Δt 6–80 s, the run stays at the start time; the failing equation is the momentum row of the innermost
  interface (corrections of u of order 1). ETOTMON G of the windows: median 6.3e-6, max 3.1e-4 (the gap G5 is closed by two orders; the rest is not yet at Newton tolerance
  because Δt is tiny). Reason: m*_j uses the advected interface density, which the van Leer scheme takes from the upwind side (the sign of the mass flux, ≈ 0 in the hydrostatic star) and which on the star's
  stretched mesh is 7–29 % below the node mass (§21.3): the force is a discontinuous function of the sign of a tiny flux and the star is no longer in hydrostatic balance
  with the pressure force. So the exact discrete compatibility of momentum row and total row is obtained at the price of a gravity force that is not a good approximation of gravity at 100 zones; it
  works where the mesh is fixed and the flow is not close to rest (Noh, free fall), not in the star. `gravcons(P) = false` restores the node-mass force.
* Heat booking with the force of the node mass (gravcons switched off by a temporary local edit, not committed), 300 steps, `t_ie29`: identical to `t_ie24` (33 rejections, 29 Newton failures, same time reached 7.508554e12 s, same G 6.6e-4): the
  CQB/YQB booking changes nothing when no zone changes status (`CPLBADDELTA=Inf`). With the default `CPLBADDELTA=100` (`t_ie30`, 150 steps) two zones were decoupled for a step; the message prints the booked heat the zone uses from then on:
  zone 32 cq 3.2e-8, yq 9.9e-8 erg cm⁻³ s⁻¹ (the heat of the previous end state, applied once, in place of the iterate's), zone 30 cq = yq = 0 (the zone had no three-level booking).
