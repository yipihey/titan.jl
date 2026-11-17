# Energy conservation of the Pop III accretion stack: measurement, derivation, compatible rows

Branch `energy-compat` (from 4f1fe93, head of `multigroup`). Commits: 774c85b (zone_diag split), a37f48b (TIEFRAME), be38295 (ECOMPAT),
50f1609 (TIEFRAME for the multigroup rows), e0e6b2e (sub-zone kinetic-energy transport, TR-BDF2 kinetic reference).
All runs: 100 zones, restart of z100mm's checkpoint `z19z100mmS5_auto0` (step 17398, M* = 0.88 M☉ in perf.csv, photosphere zone 26, 111 R☉) unless
stated. Luminosities are 4π × (per-steradian value) / L☉. Sign convention of every table: a positive row residual means the row would
create energy; a negative term added to the gas row heats.

## 0. Answer

*Is a single energy-conserving formulation possible?* Yes, and it does not need the E − K subtraction.

* The gas energy row can be the first law in **every** zone (star, shock, cold supersonic infall with K/e ~ 1e3–1e4) and still conserve
  internal + kinetic energy to round-off, if it is written as the exact discrete adjoint of the momentum row: the kinetic energy that the
  momentum row's discretisation does not keep is a sum of explicit products of velocity *differences* (never K_new − K_old), and it is
  returned to the gas as heat. On this mesh there are four such terms (§4): time discretisation, viscous-work form, upwinded momentum
  flux, and the kinetic energy of the TR-BDF2 reference state. Measured: per-zone identity residual 0.000 L☉ in zones 1–87
  (max over 1000 stages), |Σ| ≤ 1 L☉ in zones 88–100 in 90 % of the stages (round-off of the lab-frame radiation row there, §5).
* The gravity power must be booked the way the momentum row's force does it (Σ_j u_j g_j ½(m_j−1 + m_j)). The total-energy rows of the
  current stack book it as (mass flux r² u ρᵇ through the interface, ρᵇ the advected interface density) × g × Δr; the two differ per zone by
  G5 = +1.0e3 L☉ (first 0.2 yr of the restart) to +1.9e3 L☉ (time-weighted over 12.6 yr) in the star and +0.7…0.85e3 L☉ in zones
  27–50: 25–30 % of L_em in all, heating the total rows withhold from the gas. That, not the subtraction's round-off, is what makes
  "total rows everywhere" fail (t_mg25, §6.3); it is also the largest term of the dual-energy stack's ledger (§2) and the reason the
  star's structure depends on which zones are flagged (§7).
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
  (§5a). And the new rows give a *different star* than the current stack (cooler, larger; from the IC lower central temperature at
  equal mass): the current stack's result depends on a gravity term that is not the momentum row's (§7, §8).
* The isolated-star line (S5, 100 M☉, 100 zones) runs with the same rows at the same Δt, Newton cost and ledger closure as with its
  own combined row, over the first 170 yr of the contraction and only with the exchange cap EXCAP=1e8 that the accretion stack
  already uses: a gas row that carries the gas–radiation exchange explicitly is singular to round-off in an opaque interior
  (§7.3). Not run to the ZAMS.

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

No term. Σ_j u_j g_j M_j *is* the gravity power of this discretisation; a row that conserves internal + kinetic energy against it
needs nothing else. The isolated-star line's ledger (julia/s5_run.jl `energies`) already uses this form.

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
was not made (it needs ~6000 steps to the same mass).

**How far, at what cost** (perf.csv; failures = trapezoidal + BDF Newton failures since step 0):

| M* [M☉] | t_en_ic1: step / r* [R☉] / T_c [K] / log n_c / failures | z100mm: step / r* / T_c / log n_c / failures |
|---|---|---|
| core appears | 200 (t = 2.410650e5 yr) | 200 (t = 2.378241e5 yr) |
| 0.17 | 700 / 32.5 / 7.99e4 / 22.90 / 43 | 1400 / 30.7 / 1.090e5 / 23.16 / 62 |
| 0.21 | 800 / 38.3 / 8.12e4 / 22.91 / 50 | 2100 / 34.7 / 1.132e5 / 23.18 / 96 |
| 0.26 | 1100 / 44.4 / 8.28e4 / 22.92 / 86 | 3900 / 45.8 / 1.160e5 / 23.19 / 222 |
| 0.32–0.35 | 1400 (M 0.331, 5.26 yr after the core) / 50.9 / 8.44e4 / 22.93 / 120 | 5300 (M 0.321, 4.81 yr) / 54.7 / 1.190e5 / 23.20 / 295; 6500 (M 0.349, 5.40 yr) / 59.6 / 1.186e5 / 23.19 / 335 |

* Same accretion history (0.33 M☉ 5.3 yr after the core in both), **four times fewer steps** (1200 against 5100–6300 after the
  core) and about three times fewer Newton failures (95 against 258–298 after the core).
* A different star: at equal mass T_c is 0.84e5 K against 1.19e5 K and n_c 10^22.93 against 10^23.20 cm⁻³; the radius is larger at
  low mass (38.3 against 34.7 R☉ at 0.21 M☉), equal near 0.25 M☉, smaller at 0.33 M☉ (50.9 against ~56 R☉). §7 has the cause.
* The core forms at 2.41065e5 yr against 2.37824e5 yr (+1.4 %). Not attributed: the two runs differ in ECOMPAT *and* in the code
  tree, and the collapse phase was not analysed.

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
a T_CMB⁴ as the star's luminosity rises (−27 % at zones 66/67 at step 1371).

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
11. **The star differs from the current stack's** (§7.2); attributed to G5 by exclusion, not by a controlled run, and which gravity
    power converges faster in zone number is not known.
12. **The core forms 1.4 % later from the IC** (2.41065e5 against 2.37824e5 yr). Not attributed: ECOMPAT and the code tree both
    differ from z100mm, and no ECOMPAT-off run from the IC with tree e0e6b2e was made (~6000 steps).
13. **t_mg25's failure** is explained by the measured terms (§6.3), not by repeating it with one term removed.
14. **Multigroup**: the group energy rows are still comoving; ECOMPAT with MGROUPS > 0 is untested and not refused by the driver;
    `MGHPLUS` has three uninitialised ghost elements (§9).
15. **Isolated-star and quasi-star lines** not run through with the compatible rows (§7.3).
16. A size check, no more: z100mm's boundary inflow of 1.16e5 L☉ over 2.38e5 yr is 3.3e51 erg; the excess over a T_CMB⁴ stored in
    zones 88–100 of the checkpoint is 3.5e51 erg (estore.jl). The inflow was not followed through the run's history.
17. Stale settings met on the way: `CONV=1e-3` in z100mm's run.env is not read by popiii_m2c; the `TRESET` row of
    docs/POPIII_SWITCH_INVENTORY.md refers to code that is gone.

## 11. What becomes removable once ECOMPAT is the default

Not deleted in this round (the S4/S5 isolated-star line and the quasi-star line use the same residual).

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
| t_en_ic1 | e0e6b2e | from the initial condition, ECOMPAT=1 | 1478 (stopped) | §8 |
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
