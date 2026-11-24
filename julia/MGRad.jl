# ------------------------------------------------------------------ MGROUPS: implicit multigroup radiation
# Included into TitanResidFast. Design and tests: docs/POPIII_MULTIGROUP.md. Everything here is behind MGROUPS (default off: MG[] is
# false, no code path below runs and no existing path changes; bit-identical off).
"""
MGROUPS (popiii_m2c `MGROUPS=n` (default edge set) or `MGROUPS=e1,e2,…` (interior group edges in eV, ascending); default off):
the radiation field as G groups, each with its own energy density E_g and flux F_g, all Newton unknowns of the hydro system
(E_g at MGI0 + 2g − 1, F_g at MGI0 + 2g; MGI0 = 9 with COUPLED, so neqn = 9 + 2G). The grey rows of E_r, F_r become the sums
E_r = Σ_g E_g, F_r = Σ_g F_g (algebraic rows), so every term of the hydro that reads E_r, F_r and the grey Eddington factors
(the total-energy and gas+radiation rows, the pressure f E_r, the formal solution, Δt control) sees the total field unchanged, and:

* per group the moment rows of TITAN's grey system (TitanResidFast RADIATION ENERGY / MOMENTUM, comoving frame, Mihalas & Mihalas
  1984; Gehmeyr & Mihalas 1994) with the group's exchange x_cap (c κ_E,g E_g − 4π κ_P,g B_g(T)) ρ ΔV, B_g = b_g(T) σT⁴/π, and its
  interface opacity χ̄_g (harmonic, or SLABFLUX's χ_n of the group's χ_g, SUBZONE on top), the FLUXTIE closure with the group's h
  (SLABTIE: the donor h at the group's slab depth), the outer boundary F_g = c g E_g,s − π(2g+1) I_cmb b_g(T_cmb);
* the radiation force on the gas is Σ_g χ̄_g F_g/c (the grey rows' χ̄ F/c is replaced; the energy rows need nothing: their work and
  flux terms read E_r = Σ E_g and F_r = Σ F_g);
* group opacities and emissivities at the Newton iterate (stage 2b): the material returns, in the same call as the zone's grey means
  and in the residual's number type, the group record of the zone (`material_opacg!`, `mg_default!`, `mg_level`): no lagged group
  quantity enters a row. The old time level's record (MGVO) is the converged record of the previous step (TR stage only);
* the line cooling of the coupled zones (H₂/HD line transfer, the gas row's Λ) is a source of group MGLG (IR) where the gas row is
  the gas-only row, so the energy the lines carry out of the gas stays in the radiation field;
* the chemistry's energy moves (RAYQ, GASQ) come out of group MGQG (optical).
"""
const MG = Ref(false)
const MGN = Ref(0)
const MGEDGE = Ref(Float64[])                   # interior edges [eV], ascending
const MGI0 = Ref(9)
const MGEVK = Ref(1.602176634e-12 / 1.380649e-16)   # eV / k_B [K]; set by the driver from ChemNet's constants
const MGQG = Ref(2)                             # the group that carries RAYQ/GASQ (optical)
const MGH2G = Ref(0)                            # the LW group with the H₂ photon-conserving absorption (stage 2c; 0: none)
const MGH2T = Ref(0)                            # the group that receives the H₂ fluorescence (the one below the LW group)
const MGBFLOOR = 1e-200                         # floor of a group's Planck fraction (energy 1e-200 aT⁴: none; keeps ln E_g finite)
const MGEFL = 1e-3                              # floor of the group fluxes' weight in the error/Newton norms: changes of F_g below MGEFL·|F_r| do not limit Δt (the weights of E_g need none: see MGSC)
const MGCUFL = 1e-60                            # [erg cm⁻² s⁻¹] the star-light closure's flux floor: c u_g below it is not seen by the zone's photo-rates and colour temperature (MGOpac)
"""
Scale of the group unknown, per zone and group (`MGSC[g][k]`): y_g,k = asinh(E_g,k/s_g,k), with s_g,k the energy density that the zone's own emission and absorption would give,
the local emission balance s = j_g/(c κ_E,g) = 4π[κ_P,g b_g(T) + A_1 X_g(T)] σT⁴/π / (c κ_E,g) (`mg_scales!`), evaluated at the OLD state of the step (the state the step starts from) and held
fixed through all the Newton iterations of the step (and the error estimate's solves). A parametrisation, not a floor and not a model: for any s > 0 the map y ↔ E = s sinh y is a bijection of the
reals, so the discrete equations and their solution do not depend on s; s only fixes the unit in which the Newton iteration, its linear algebra and its convergence test see E_g. Above s the unknown is ln E_g, below it
linear, and the column of E_g in its own rows is s ∂R/∂E_g, the size of the group's local source terms: an empty group (E_g ≪ s) is determined by its rows at the level of the other unknowns, where a global
scale (the former constant MGES = MGCUFL/c = 3.3e-71, then 1e-150 E_r for ln E_g) left the diagonal at 1e-49 of the row's largest entry (docs §22.3, §23). `MGES` remains only as the value of s where the
balance is undefined (no emission or no absorption in the group: s ≤ 0 or not finite), which is where the choice of s is immaterial by the bijection.
"""
const MGES = MGCUFL / 2.99792458e10
const MGSTEADY = Ref(false)                     # true only inside `mg_radsolve!` (popiii_m2c): the group rows without their time derivatives (docs §26)
const MGSC = Ref(Vector{Vector{Float64}}())     # per-zone scale s_g,k of each group's unknown (see above)
const MGTCMB = Ref(0.0)
vvec(n, m) = [zeros(m) for _ in 1:n]
const MGEN = Ref(Vector{Vector{Float64}}()); const MGFN = Ref(Vector{Vector{Float64}}())     # new level (Newton iterate base)
const MGEO = Ref(Vector{Vector{Float64}}()); const MGFO = Ref(Vector{Vector{Float64}}())     # old level (start of the stage)
const MGESO = Ref(Vector{Vector{Float64}}()); const MGFSO = Ref(Vector{Vector{Float64}}())   # limiter slopes of the old E_g/ρ and face F_g
const MGFNOM = Ref(Vector{Vector{Float64}}())   # nominal F_g scale per interface (correction F_g += δ·nom)
"""
The group opacity record of a zone (stage 2b): a vector of `mg_nv()` numbers that the material fills at the Newton iterate, in the
residual's number type, in the same call that returns the zone's grey means (`material_opacg!`; PopIIIMaterial.nlo_it). Per group g
(index `mgv(g, j)`): 1 r_E = κ_E,g/κ_E, 2 r_P = κ_P,g/κ_P, 3 r_R = χ_g/χ (ratios to the grey means of the same call; the grey means carry
the run's calibration), 4 r_Ec: the continuum part of r_E (without the lines' share; enters the thin-zone flux mean, `mg_chif`),
5 K_1: H(1s) (and n ≥ 2) bound–free absorption of case-B zones [cm² g⁻¹], 6 A_1: the H(1s) Milne emission coefficient
(emission 4π A_1 X_g(T) σT⁴/π, X_g = e^{χ_H/kT} b_g(T), `mg_xfrac`), 7 κ_Ram. Per Raman pair p (`mgvp(p, j)`): 1 κ_Ram C, 2 κ_Ram ω. Last
two entries (diagnostics): the three-level n = 2 density and the colour temperature T_w of the absorption means. Table zones and zones without a composition record: `mg_default!` (r_E = r_P = r_R = 1,
rest 0: the grey field split by Planck fractions, the diffusion limit of the thick interior; the groups then sum to the grey rows).
`MGVN`, `MGVO` [mg_nv() × zones]: the Float64 records of the new level (PopIIIMaterial fill_opac!, wherever TITAN refreshes chifn) and of
the old level (copied at the start of a step, as chifo): the residual's θ-mix takes old-level group opacities from MGVO, exactly as the
grey rows take chifo, xkeo, xkpo.
"""
const MGVN = Ref(zeros(0, 0)); const MGVO = Ref(zeros(0, 0))
const MGXH = Ref(13.598434)                     # χ_H [eV] of the Milne emission (set by PopIIIMaterial to the G2 continuum's)
const MGGH = Ref(typemax(Int))                  # first ionising group (role H)
const MGLG = Ref(0)                             # group that receives the line cooling of the coupled zones (the H₂ rovibrational lines: IR)
const MGLCO = Ref(Float64[])                    # old-level line cooling per volume [erg cm⁻³ s⁻¹] (as XLAMO)
@inline mgv(g, j) = 7 * (g - 1) + j
@inline mgvp(p, j) = 7 * MGN[] + 2 * (p - 1) + j
mg_nv() = 7 * MGN[] + 2 * length(MGRAMP[]) + 2
@inline function mg_default!(o)
    fill!(o, zero(eltype(o)))
    @inbounds for g in 1:MGN[]; o[mgv(g, 1)] = one(eltype(o)); o[mgv(g, 2)] = one(eltype(o)); o[mgv(g, 3)] = one(eltype(o)); end
    o
end
"""
Flux-mean extinction of a group in a zone of column m [g cm⁻²]. χ_R: the group's Rosseland-type mean (harmonic, weight ∂B_ν/∂T), exact
where the zone is optically thick in the group (the flux has the diffusion spectrum; the groups then sum to the grey diffusion flux).
Where the zone is thin the flux has the spectrum of the energy density, and the flux mean is the absorption mean κ_a the energy row
removes (continuum κ_E,g + K_1 + κ_Ram): with a flux mean below it a streaming absorbing layer has (F/cE)² = κ_a/χ > 1 (steady moment
equations, f = 1). Weight of the thin form: e^{−χ_R m}, the transmission of the zone in the Rosseland mean. (Stage 2a took
max(χ_R, κ_E,g) everywhere with the line Planck mean in κ_E,g: 3–400 × the Rosseland mean in the stellar interior and below the
photosphere, where the groups then no longer summed to the grey diffusion flux; docs §17.)
"""
@inline mg_chif(χR, κa, m) = val(κa) > val(χR) ? χR + exp(-χR * m) * (κa - χR) : χR
"Column of zone k [g cm⁻²] from the density, the zone volume/(μ+1)-form and the face areas r^μ (0 where the faces vanish)."
@inline mg_zcol(d, dvol, rmu, k) = (a = 0.5 * (rmu[k] + rmu[k+1]); (val(a) > 0 && val(dvol[k]) > 0) ? d[k] * dvol[k] / a : zero(d[k] * dvol[k]))
"(κ_E,g, χ_g, κ_P,g, A_1) of group g (j = mgv(g, 1)) in zone k from the record gv and the grey (χ, κ_E, κ_P) of the same level; m: zone column."
@inline function mg_level(gv, j, k, χ, κE, κP, m)
    k1 = gv[j+4, k]
    (κE * gv[j, k] + k1, mg_chif(χ * gv[j+2, k], κE * gv[j+3, k] + k1 + gv[j+6, k], m), κP * gv[j+1, k], gv[j+5, k])
end
"Size MGVN/MGVO (after the Raman pairs are set) and set every zone to the default record."
function mg_valloc!(mgr)
    nv = mg_nv()
    for R_ in (MGVN, MGVO)
        size(R_[]) == (nv, mgr) || (R_[] = zeros(nv, mgr); for k in 1:mgr; mg_default!(view(R_[], :, k)); end)
    end
    length(MGLCO[]) == mgr || (MGLCO[] = zeros(mgr))
    empty!(MGWSC)                               # (workspaces are sized by mg_nv())
    nothing
end
"""
Raman scattering of H(1s) in a group gL between Lyα and the Lyman limit (docs §16; cross sections: HScatter, Kokubo 2024): H(1s) + γ(ν) ↔
H(2s) + γ(ν − ν_Lyα) destroys the LW photon (the atom returns by two-photon decay) and its inverse creates one from an optical photon and
an H(2s). `MGRAMP`: pairs (gL, g′, e1, e2) — the photons of gL between e1 and e2 [eV] have their Stokes partner in group g′. Per mass
and per pair the transfer into gL is c κ_Ram (C E_g′ − ω E_gL) and the same out of g′ (energy moves between the two groups only),
  ω = b(e1, e2; T_w)/b_gL(T_w)  (the pair's share of gL),   C = (n_2s/n_1s) e^{E_21/kT_w} b(e1, e2; T_w)/b_g′(T_w),
T_w the colour temperature of the light in the zone (NLOPACJ): by detailed balance of the Kramers–Heisenberg rates the inverse-Raman
source function is S_ν = (2hν³/c²)(n_2s/n_1s) n̄(ν − ν_Lyα) (g_1s = g_2s), in the Wien limit (n_2s/n_1s) e^{E_21/kT_w} W_g′ B_ν(T_w) for a
partner field W_g′ B_ν′(T_w). In LTE (n_2s/n_1s = e^{−E_21/kT}, E_g = (4π/c) B_g(T)) every pair is in balance exactly.
κ_Ram [cm² g⁻¹] from G2Opacity.group_fast; the zone's record (MGVN) holds κ_Ram C and κ_Ram ω per pair.
"""
const MGRAMP = Ref(Tuple{Int,Int,Float64,Float64}[])
"Set the Raman pairs (PopIIIMaterial.mg_roles!)."
mg_rampairs!(pairs) = (MGRAMP[] = pairs; nothing)
"""
Every group's own Eddington factors from a formal solution with the group's
χ_g, κ_E,g, emission and background (popiii_m2c `mg_eddfac!`, run after the grey one wherever TITAN calls it: every Newton update with VEFIMP):
f_g per zone (MGFEDD), the tie's h_g (MGHEDD) and weights (MGTIEW, from the group's own optical depth), the outer g_g (MGGEDDR). The grey
array fedd then holds Σ_g f_g E_g / Σ_g E_g (the radiation pressure of the total rows, MLT, …). Frozen in the Jacobian, as the grey factors.
"""
const MGFEDD = Ref(Vector{Vector{Float64}}()); const MGHEDD = Ref(Vector{Vector{Float64}}()); const MGTIEW = Ref(Vector{Vector{Float64}}())
const MGGEDDR = Ref(Float64[])
const MGVEFOK = Ref(false)                       # the group formal solutions have run (before the first one: the grey factors)
const MGHPLUS = Ref(Vector{Vector{Float64}}())   # (diagnostic) outgoing half moment H⁺ of each group's formal solution per interface
mgiE(g) = MGI0[] + 2g - 1; mgiF(g) = MGI0[] + 2g

"Allocate the group arrays for mgr zones (values kept where present)."
function mg_alloc!(mgr)
    G = MGN[]
    for R_ in (MGEN, MGFN, MGEO, MGFO, MGESO, MGFSO, MGFNOM, MGSC)
        (length(R_[]) == G && all(length(a) == mgr for a in R_[])) || (R_[] = vvec(G, mgr))
    end
    for a in MGSC[]; all(==(0.0), a) && fill!(a, MGES); end
    for R_ in (MGFEDD, MGHEDD, MGTIEW, MGHPLUS)
        (length(R_[]) == G && all(length(a) == mgr for a in R_[])) || (R_[] = vvec(G, mgr))
    end
    length(MGGEDDR[]) == G || (MGGEDDR[] = fill(0.5, G))
    nothing
end

# ---------------------------------------------------------------- Planck group fractions (generic in the number type)
"∫₀ˣ t³/(eᵗ − 1) dt for x < 1: Bernoulli series B₀…B₁₀ (as RadGroups of branch radgroups-3g; next term ≲ 1e-10 relative at x = 1)."
@inline function mg_plow(x)
    x2 = x * x
    x^3 * (1 / 3 - x / 8 + x2 / 60 - x2^2 / 5040 + x2^3 / 272160 - x2^4 / 13305600 + x2^5 * (5 / 66) / (13 * 3628800))
end
"∫ₓ^∞ t³/(eᵗ − 1) dt for x ≥ 1: Σₙ e^{−nx}(x³/n + 3x²/n² + 6x/n³ + 6/n⁴) (no cancellation for large x)."
@inline function mg_phigh(x)
    s = zero(x); xv = val(x)
    for n in 1:60
        e = exp(-n * x)
        s += e * (x^3 / n + 3x^2 / n^2 + 6x / n^3 + 6 / n^4)
        exp(-n * xv) < 1e-18 && break
    end
    s
end
"Fraction of σT⁴/π above x = hν/kT."
@inline mg_above(x) = val(x) >= 1 ? mg_phigh(x) * (15 / π^4) : 1 - mg_plow(x) * (15 / π^4)
"Fraction of σT⁴/π below x."
@inline mg_below(x) = val(x) >= 1 ? 1 - mg_phigh(x) * (15 / π^4) : mg_plow(x) * (15 / π^4)
"""
Planck fraction b_g(T) of group g (edges MGEDGE, eV), floored at MGBFLOOR. Interior groups whose lower edge has x ≥ 1 are the difference
of two tails (no cancellation); else of two lower integrals.
"""
@inline function mg_bfrac(g, T, edges = MGEDGE[])
    G = length(edges) + 1; s = MGEVK[] / T
    b = if g == 1
        mg_below(edges[1] * s)
    elseif g == G
        mg_above(edges[G-1] * s)
    else
        xl = edges[g-1] * s; xh = edges[g] * s
        val(xl) >= 1 ? mg_phigh(xl) * (15 / π^4) - mg_phigh(xh) * (15 / π^4) : mg_below(xh) - mg_below(xl)
    end
    val(b) < MGBFLOOR ? MGBFLOOR + zero(b) : b
end
"Σₙ e^{xc − nx}(x³/n + 3x²/n² + 6x/n³ + 6/n⁴): mg_phigh(x) e^{xc} without overflow or underflow of either factor."
@inline function mg_phs(x, xc)
    s = zero(x); xv = val(x)
    for n in 1:60
        e = exp(xc - n * x)
        s += e * (x^3 / n + 3x^2 / n^2 + 6x / n^3 + 6 / n^4)
        exp(-n * xv) < 1e-18 && break
    end
    s
end
"""
X_g(T) = e^{χ/kT} b_g(T) for a group g ≥ 2 and a threshold χ [eV] (no floor): the temperature factor of the Milne recombination emission
of a bound–free continuum with threshold χ into group g, j_g = n_e n_i λ_e³ (g_n/2g_i) ⟨σ⟩_g e^{χ/kT} b_g(T) σT⁴/π. The two exponentials
are combined term by term, so the product is exact where b_g alone underflows (cold gas): X_g → (15/π⁴)(χ/kT)³ e^{−(e_g − χ)/kT} there.
"""
@inline function mg_xfrac(g, T, χ, edges = MGEDGE[])
    G = length(edges) + 1; s = MGEVK[] / T; xc = χ * s; xl = edges[g-1] * s
    val(xl) >= 1 || return mg_bfrac(g, T, edges) * exp(xc)
    b = mg_phs(xl, xc); g < G && (b -= mg_phs(edges[g] * s, xc))
    b * (15 / π^4)
end
"Fraction of σT⁴/π between the photon energies e1 < e2 [eV] (no floor)."
@inline function mg_band(e1, e2, T)
    s = MGEVK[] / T; xl = e1 * s; xh = e2 * s
    val(xl) >= 1 ? (mg_phigh(xl) - mg_phigh(xh)) * (15 / π^4) : mg_below(xh) - mg_below(xl)
end

# ---------------------------------------------------------------- Newton hooks (chained with the COUPLED ones by popiii_m2c)
"Old-level limiter slopes of E_g/ρ (zones) and of the face value F_g·2/(ρ₋+ρ₊), as the species slopes (cpl_newstep!)."
function mg_slopes!()
    ngrs = gi(:ngrs); ngre = gi(:ngre); cadv = gf(:cadv); eps = gf(:epsadv); dO = v(:do)
    q = get!(() -> zeros(TF.MGR), BUFS, :mgq)::Vector{Float64}
    for g in 1:MGN[]
        for (kind, src, sl) in ((1, MGEO[][g], MGESO[][g]), (2, MGFO[][g], MGFSO[][g]))
            fill!(sl, 0.0)
            @inbounds for k in ngrs-2:ngre+2
                q[k] = kind == 1 ? src[k] / dO[k] : (k > 1 ? src[k] * 2 / (dO[k] + dO[k-1]) : 0.0)
            end
            @inbounds for k in ngrs:ngre+1
                a = q[k] - q[k-1]; b = q[k+1] - q[k]; sdq = a + b; den = sdq + (sdq >= 0 ? eps : -eps)
                sl[k] = a * b >= 0 ? cadv * a * b / den : 0.0
            end
        end
    end
    nothing
end
"""
The other sources of each group's energy row at the start state, per mass [erg g⁻¹ s⁻¹] (`MGSC` adds them to the emission, see `mg_scales!`): exactly the terms of `mg_rows!` that do not contain
the zone's own E_g — the chemistry and cycle heat cq + yq (booked, `CQB`/`YQB`) and RAYQ, GASQ of the group MGQG, the line cooling MGLCO of the coupled zones into the group MGLG, 
(docs §24: zone 255 of group 2 saw cq + yq = 1.7e27 erg s⁻¹ against s = 1.3e-32). The energy carried in by the mass flux is `mg_advscale!`. A scale only: magnitudes, any s > 0 is the same discrete solution.
"""
function mg_sources!(G)
    ngrs = gi(:ngrs); ngre = gi(:ngre)
    q = [zeros(ngre + 2) for _ in 1:G]
    @inbounds for k in ngrs:ngre
        g = MGQG[]
        g > 0 && (q[g][k] = (abs(isempty(CQB[]) ? 0.0 : CQB[][k]) + abs(isempty(YQB[]) ? 0.0 : YQB[][k]) + abs(isempty(RAYQ[]) ? 0.0 : RAYQ[][k]) + abs(isempty(GASQ[]) ? 0.0 : GASQ[][k])) / v(:dn)[k])
        g = MGLG[]
        (g > 0 && cplz(k)) && (q[g][k] += abs(MGLCO[][k]) / v(:dn)[k])
    end
    q
end
"""
What reaches a zone from its neighbours, in the scale. The row has the advective terms ṁ_k (E/ρ)ᵇ_k (upwind value, `advect!`) of its two faces and the flux divergence of the neighbours' fields: a zone next to one holding energy
gets E ≈ (ṁ (E/ρ)_up/m + c E_nb/(ρΔr))/(c κ_g) from them (zie3: the 1.2e-16 erg cm⁻³ of group 3 in zones 212–213 reach zone 214, which held 1.6e-84, in an outflow of 2.6e22 g s⁻¹, and zone 215 from 214, …:
a scale without them left the correction y of those zones at 1e9 – 1e12 and the Newton step at the box |δy| ≤ 4 for 25 iterations). The neighbour's value is the one its own scale allows, √(E_old² + s²) (as the Raman partner of
`mg_scales!`), a sweep up and down the grid chains it (either direction of ṁ). The free-streaming bound c E_nb/Δr is the largest inflow the faces can carry; the time term 1/(ρ c Δt) in the denominator (`mg_scales!`) is what limits it to c Δt/Δr
of the neighbour per zone. Magnitudes at the start state; a scale only.
"""
function mg_advscale!(g, ke, base)
    ngrs = gi(:ngrs); ngre = gi(:ngre); dt = gf(:dtime); xme = v(:xmen); xmo = v(:xmeo); dO = v(:do); cc = gf(:cc); rn = v(:rn)
    E = MGEO[][g]; sc = MGSC[][g]
    mdot(k) = -(xme[k] - xmo[k]) / dt                 # inward mass flux through face k (dmdt of the residual); flow = -ṁ: outward takes the inner zone
    for order in (ngrs:ngre, ngre:-1:ngrs)
        @inbounds for k in order
            m = xme[k] - xme[k+1]; qa = 0.0
            m > 0 || continue
            mdot(k) < 0 && (qa += -mdot(k) * hypot(E[k-1], sc[k-1]) / dO[k-1])                # outflow through the lower face brings in the inner zone's energy
            mdot(k + 1) > 0 && (qa += mdot(k + 1) * hypot(E[k+1], sc[k+1]) / dO[k+1])        # inflow through the upper face brings in the outer zone's
            qr = cc * (hypot(E[k-1], sc[k-1]) + hypot(E[k+1], sc[k+1])) / (dO[k] * (rn[k+1] - rn[k]))      # free-streaming inflow bound from the two neighbours, per mass
            s_ = base[k] + (qa / m + qr) / (cc * ke[k])
            (isfinite(s_) && s_ > 0) && (sc[k] = s_)
        end
    end
    nothing
end
"""
Scales of the group unknowns at the state the step starts from (`MGSC`, see there; evaluated after the predictor, `mg_predicted!`, so that the opacity records are those of that state:
after a restart from a grey checkpoint the old-level records are the defaults, which gave no Raman or Milne source at all): the energy density that the zone's own sources and sinks would give,
s_g,k = (4π em_g + q_g) / (c (κ_g + 1/(ρ Δr) + 1/(c ρ Δt))),  (q_g: `mg_sources!`; + the neighbours' inflow, `mg_advscale!`)  em_g = κ_P,g b_g(T)(σT⁴/π) + A_1 X_g(T)(σT⁴/π) + Σ_pairs κ_Ram C (c/4π) Ê_partner,  κ_g = κ_E,g + the Raman absorption κ_Ram ω,
with the group record (MGVN), the grey means, T, ρ and zone width Δr of that state: exactly the source and absorption of the group's formal solution (`mg_eddfac!`) at the old level.
Two refinements, both observed necessary (docs §23.2): (1) the term 1/(ρΔr) is the zone's escape: a zone thin in the group has E = em Δr/c, not em/(cκ) (without it s sat 1e56 above the solution in the cold outer
zones); (2) the Raman partner's energy density is Ê = √(E_old² + s_partner²) with the partner's own scale s_partner of the first pass (the emission balance without the Raman terms): a group fed by a partner
must have the scale of what the partner will be, not of what it was (an empty partner, E = 1e-72, with s = 1e-24, put the column of the partner 1e51 above the group's own column in the group's row).
Where the balance is undefined (no source or no sink: s not positive and finite) the scale is MGES, a choice that cannot enter the solution (bijection y ↔ E for every s > 0). Ghost zones take the scale of the nearest zone.
The H₂ fluorescence is not in it. (The chemistry's energy moves of group MGQG, the line cooling of MGLG and the neighbours' inflow are, docs §24: the cq + yq of zone 255 was the source of 23.6's correction.)
"""
function mg_scales!()
    ngrs = gi(:ngrs); ngre = gi(:ngre); cc = gf(:cc); to = v(:tn); do_ = v(:dn); rn = v(:rn); gv = MGVN[]; xh = MGXH[]; gH = MGGH[]
    dt = gf(:dtime); chifo = v(:chifn); xkeo = v(:xken); xkpo = v(:xkpn); plfo = v(:plfn); c4pi = cc / (4π); G = MGN[]
    em0 = [zeros(ngre + 2) for _ in 1:G]; ke0 = [zeros(ngre + 2) for _ in 1:G]; qm = mg_sources!(G)
    @inbounds for g in 1:G
        j = mgv(g, 1); sc = MGSC[][g]
        for k in ngrs:ngre
            keO, _, kpO, a1O = mg_level(gv, j, k, chifo[k], xkeo[k], xkpo[k], 0.0)
            em = kpO * plfo[k] * mg_bfrac(g, to[k])
            (g >= gH && a1O > 0) && (em += a1O * plfo[k] * mg_xfrac(g, to[k], xh))
            em0[g][k] = em; ke0[g][k] = keO + 1 / (do_[k] * (rn[k+1] - rn[k])) + 1 / (cc * do_[k] * dt)
            s_ = (4π * em + qm[g][k]) / (cc * ke0[g][k]); sc[k] = (isfinite(s_) && s_ > 0) ? s_ : MGES        # first pass: no Raman terms
        end
    end
    kefin = [copy(ke0[g]) for g in 1:G]
    isempty(MGRAMP[]) || begin
        s0 = [copy(MGSC[][g]) for g in 1:G]
        @inbounds for k in ngrs:ngre
            em = [em0[g][k] for g in 1:G]; ke = [ke0[g][k] for g in 1:G]
            for (p, pr_) in enumerate(MGRAMP[])
                kc = gv[mgvp(p, 1), k]; kw = gv[mgvp(p, 2), k]; (kc > 0 || kw > 0) || continue
                gL, gp = pr_[1], pr_[2]
                ke[gL] += kw; em[gL] += kc * c4pi * hypot(MGEN[][gp][k], s0[gp][k])
                ke[gp] += kc; em[gp] += kw * c4pi * hypot(MGEN[][gL][k], s0[gL][k])
            end
            for g in 1:G; s_ = (4π * em[g] + qm[g][k]) / (cc * ke[g]); MGSC[][g][k] = (isfinite(s_) && s_ > 0) ? s_ : MGES; kefin[g][k] = ke[g]; end
        end
    end
    for g in 1:G; mg_advscale!(g, kefin[g], copy(MGSC[][g])); end
    for (g, z) in MGROWPROBE[]
        k = ngrs + z - 1
        println("MGSCALE g ", g, " zone ", z, " s ", MGSC[][g][k], " emit0 ", em0[g][k], " ke0 ", ke0[g][k], " E ", MGEN[][g][k], " T ", to[k], " rho ", do_[k], " plf ", plfo[k], " bfrac ", mg_bfrac(g, to[k]),
                " raman ", [(p, gv[mgvp(p, 1), k], gv[mgvp(p, 2), k]) for p in eachindex(MGRAMP[])], " kappa_rec ", gv[mgv(g, 1), k], " H2g ", 0)
    end
    for g in 1:G
        sc = MGSC[][g]
        for k in ngrs-2:ngrs-1; sc[k] = sc[ngrs]; end
        for k in ngre+1:ngre+2; sc[k] = sc[ngre]; end
    end
    nothing
end
"Start of a step: old := new (fields and opacity records), slopes; provisional scales (replaced by `mg_predicted!` once the predictor has refreshed the state)."
function mg_newstep!()
    size(MGVO[]) == size(MGVN[]) && copyto!(MGVO[], MGVN[])
    for g in 1:MGN[]; copyto!(MGEO[][g], MGEN[][g]); copyto!(MGFO[][g], MGFN[][g]); end
    mg_slopes!()
    mg_predicted!()
end
"""
After the predictor (TitanCore.predict!, hook `CPLH.predicted`): the scales of the group unknowns (`mg_scales!`) and the nominal flux scale of each interface, F_g column scale
c·√(E²+s²) (the flux bound c E_g, with the unknown's own scale) at least |F_g|·1, held fixed through the step's Newton iterations.
"""
function mg_predicted!()
    ngrs = gi(:ngrs); ngre = gi(:ngre)
    mg_scales!()
    for g in 1:MGN[]
        E = MGEN[][g]; F = MGFN[][g]; nom = MGFNOM[][g]; sc = MGSC[][g]
        @inbounds for k in ngrs-2:ngre+2
            km = max(k-1, 1)
            nom[k] = max(abs(F[k]), 1e-3 * gf(:cc) * 0.5 * (hypot(E[km], sc[km]) + hypot(E[k], sc[k])))
        end
    end
    if get(ENV, "MGPIV", "0") == "1"                  # (diagnostic) the scales, fields and nominal fluxes of the zones named by MGPIVZ (default 214,245)
        for z in parse.(Int, split(get(ENV, "MGPIVZ", "214,245"), ',')), g in 1:MGN[]
            k = ngrs + z - 1; (ngrs <= k <= ngre) || continue
            println("MGPIV scales zone ", z, " group ", g, ": s ", MGSC[][g][k], " E ", MGEN[][g][k], " F ", MGFN[][g][k], " nom ", MGFNOM[][g][k], " T ", v(:tn)[k], " rho ", v(:dn)[k],
                    " kappa_E,g ", MGVN[][mgv(g, 1), k], " r_P ", MGVN[][mgv(g, 2), k], " K1 ", MGVN[][mgv(g, 5), k], " A1 ", MGVN[][mgv(g, 6), k], " ", [(MGVN[][mgvp(p, 1), k], MGVN[][mgvp(p, 2), k]) for p in eachindex(MGRAMP[])])
        end
    end
    nothing
end
"Reset after a failed iterate: new := old."
mg_reset!() = (for g in 1:MGN[]; copyto!(MGEN[][g], MGEO[][g]); copyto!(MGFN[][g], MGFO[][g]); end;
               size(MGVO[]) == size(MGVN[]) && copyto!(MGVN[], MGVO[]); nothing)
"""
The group energy density as a function of its unknown, the correction δ of y = asinh(E/s) at the Newton iterate E₀ (s the zone's scale `MGSC[g][k]`, h = √(E₀² + s²) = dE/dδ at δ = 0):
growth (δ ≥ 0): E = E₀ + h δ, linear; decay (δ < 0): E = E₀ cosh δ + h sinh δ = s sinh(asinh(E₀/s) + δ), (= E₀ e^δ for E₀ ≫ s: positive; linear below s: E = 0 and small negative
iterates representable). Both branches have the slope h at δ = 0 (the Jacobian columns are those of E itself in units of h). A Newton step that raises an empty group is a step of a linear row system and
must be taken as such: the pure sinh map applies a correction δ of 1e11 (the band step of an empty group's LW energy at zie3) as E = s sinh(1e11) = ∞, and the box and trust region (|δ| ≤ 4) then allow
only e⁴ per iteration. A parametrisation: the discrete equations and their solution do not depend on it (docs §23.2).
"""
@inline mg_eofy(E0, s, δ) = val(δ) >= 0 ? E0 + hypot(E0, s) * δ : E0 * cosh(δ) + hypot(E0, s) * sinh(δ)
"Apply the Newton correction: E_g ← `mg_eofy(E_g, s, δ)` (δ the correction of y_g: linear growth, sinh decay), F_g ← F_g + δ·nom."
function mg_update!(rhs, a, b)
    @inbounds for g in 1:MGN[]
        E = MGEN[][g]; F = MGFN[][g]; nom = MGFNOM[][g]; iE = mgiE(g); iF = mgiF(g)
        sc = MGSC[][g]
        for k in a:b; E[k] = mg_eofy(E[k], sc[k], rhs[iE, k]); end
        for k in a+1:b+1; F[k] += rhs[iF, k] * nom[k]; end
        mg_ghosts!(E, F, a, b)
    end
    if gi(:lrad) == 1               # E_r and F_r are the groups' sums (redundant unknowns, their rows are identities, see `mg_esum!`): the sums of the new group fields
        er = v(:ern); fr = v(:frn)
        @inbounds for k in a:b; s = 0.0; for g in 1:MGN[]; s += MGEN[][g][k]; end; er[k] = s; end
        @inbounds for k in a:b+1; s = 0.0; for g in 1:MGN[]; s += MGFN[][g][k]; end; fr[k] = s; end
    end
    nothing
end

"Ghost values of a group as its ghost rows: E copied, F antisymmetric at the centre (lribc = 2) and copied outside."
@inline function mg_ghosts!(E, F, a, b)
    E[a-1] = E[a]; E[b+1] = E[b]; E[b+2] = E[b+1]
    F[a] = 0.0; F[a-1] = gi(:lribc) == 2 ? -F[a+1] : F[a-1]; F[b+2] = F[b+1]
    nothing
end
"maxdel: the convergence measure is the largest |δ y_g| (the correction of the unknown itself: |δ ln E_g| for E_g ≫ s_g,k, |δE_g|/s_g,k below, s the zone's local emission balance `MGSC`; no weight, no clamp)."
function mg_dmax!(rhs, a, b)
    d = 0.0
    @inbounds for g in 1:MGN[]
        iE = mgiE(g)
        for k in a:b
            δ = rhs[iE, k]
            isfinite(δ) || (rhs[iE, k] = 0.0; d = Inf; continue)
            d = max(d, abs(δ))
        end
        for k in a:b+1; isfinite(rhs[mgiF(g), k]) || (rhs[mgiF(g), k] = 0.0; d = Inf); end
    end
    d
end

# ---------------------------------------------------------------- residual pieces
struct MGWS{T}
    egn::Vector{Vector{T}}; fgn::Vector{Vector{T}}; eg::Vector{Vector{T}}; fg::Vector{Vector{T}}
    chif::Vector{Vector{T}}; chifn::Vector{Vector{T}}; xke::Vector{Vector{T}}; emit::Vector{Vector{T}}; avchi::Vector{Vector{T}}
    erb::Vector{Vector{T}}; ersn::Vector{Vector{T}}; qe::Vector{Vector{T}}; qeo::Vector{Vector{T}}
    frb::Vector{Vector{T}}; frsn::Vector{Vector{T}}; qF::Vector{Vector{T}}; qFo::Vector{Vector{T}}
    frc::Vector{T}                                  # Σ_g χ̄_g F_g at interfaces (radiation force)
    xch::Vector{T}                                  # Σ_g x_cap,g (c κ_E,g E_g − 4π κ_P,g B_g) per zone (the GASROW gas row's exchange, per mass)
    xfer::Vector{Vector{T}}                         # energy gained from other groups per mass and time (Raman pairs; Σ_g = 0)
    gv::Matrix{T}                                   # the zones' group opacity records at the iterate (see MGVN)
    lcl::Vector{T}                                  # line cooling of the coupled zones at the iterate [erg cm⁻³ s⁻¹] (into group MGLG)
    kh2::Vector{T}; kh2g::Vector{T}                 # H₂ Lyman–Werner absorption of group MGH2G at the iterate [cm² g⁻¹]: all of it, its share that heats the gas
    qs::Matrix{T}                                   # (G+1) × zones: c u_g/c u_g(start of step) of every group, and the LW factor f_sh J(0)/J (MGSLQ, MGH2FN)
    fdd::Vector{Vector{T}}; hde::Vector{Vector{T}}; tiw::Vector{Vector{T}}; gdr::Vector{T}   # the groups' Eddington factor f_g, tie h_g, tie weight and outer g_g in the Jacobian-vector product's Duals (mg_eddfac_dual!)
    fwt::Vector{T}                                  # the grey f of the total rows there: Σ_g f_g E_g/Σ_g E_g
    sc1::Vector{T}; sc2::Vector{T}; sc3::Vector{T}  # work: χ_g, κ_E,g, source of one group
end
const MGWSC = Dict{DataType,Any}()
mgws(::Type{T}, mgr) where {T} = get!(MGWSC, T) do
    G = MGN[]; z() = [zeros(T, mgr) for _ in 1:G]; gv = zeros(T, mg_nv(), mgr)
    for k in 1:mgr; mg_default!(view(gv, :, k)); end
    MGWS{T}(z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), zeros(T, mgr), zeros(T, mgr), z(), gv, zeros(T, mgr), zeros(T, mgr), zeros(T, mgr), ones(T, G + 1, mgr),
          z(), z(), z(), zeros(T, G), zeros(T, mgr), zeros(T, mgr), zeros(T, mgr), zeros(T, mgr))
end::MGWS{T}
"""
E_r = Σ_g E_g and F_r = Σ_g F_g as functions of the iterate (MGROUPS): every use of the grey E_r, F_r in the residual (the gas and total energy rows, the radiation pressure, the grey momentum
terms) reads the sums of the groups' unknowns, so E_r and F_r are not independent unknowns: their rows are the identities (δ_E_r = 0, δ_F_r = 0: `mg_rows!`) and `mg_update!` sets the arrays to the sums.
(The former rows E_r − Σ E_g, F_r − Σ F_g left E_r a second, redundant unknown: Jacobian columns of E_r that a finite difference through `update!` never saw exactly, and a constraint that the
nonlinear group update satisfied only to second order. Docs §23.4.)
"""
function mg_esum!(W, M::MGWS{T}, k0, ngre) where {T}
    @inbounds for k in k0:ngre+2
        se = zero(T); sf = zero(T)
        for g in 1:MGN[]; se += M.egn[g][k]; sf += M.fgn[g][k]; end
        W.ern[k] = se; W.frn[k] = sf
    end
    nothing
end
"The record column of zone k in the workspace of number type T (the material sweep writes it; zones are distinct per task)."
@inline mg_gcol(::Type{T}, k) where {T} = view((MGWSC[T]::MGWS{T}).gv, :, k)
"The line cooling of coupled zone k at the iterate (TitanResidFast's coupled-zone sweep)."

"The group unknowns of the iterate (blocks k0 … ngre+2) into the workspace."
function mg_iterate!(M::MGWS{T}, x, neqn, nb, k0, ngre) where {T}
    @inbounds for g in 1:MGN[]
        en = M.egn[g]; fn = M.fgn[g]; E = MGEN[][g]; F = MGFN[][g]; nom = MGFNOM[][g]; iE = mgiE(g); iF = mgiF(g)
        for k in k0:ngre+2; en[k] = E[k]; fn[k] = F[k]; end
        for bb in 1:nb
            k = k0 + bb - 1
            en[k] = mg_eofy(E[k], MGSC[][g][k], x[(bb-1)*neqn+iE]); fn[k] = F[k] + x[(bb-1)*neqn+iF] * nom[k]
        end
    end
    nothing
end

"""
θ-level group quantities after the grey material (chifn, xken, xkpn, plfn of the iterate and their old-level copies): per level
κ_E,g = κ_E r_E + K_1, χ_g = mg_chif(χ r_R, κ_E r_Ec + K_1 + κ_Ram), κ_P,g = κ_P r_P with the zone's record of that level (iterate: M.gv, in
the residual's number type; old: MGVO), θ-mixed as the grey opacities; emission 4π[κ_P,g b_g(T) + A_1 X_g(T)] σT⁴/π with θ-weighted
coefficients and θ-weighted Planck functions (as the grey row); interface opacity χ̄_g (harmonic, SLABFLUX, SUBZONE as the grey), the
advected face values, and the Raman transfer between the pairs' groups.
"""
function mg_theta!(M::MGWS{T}, c, W, θ, ngrs, ngre, fle, flu) where {T}
    P = c.p; tn = W.tn; to = v(:to); dn = W.dn; dO = c.do_; gv = M.gv; go = MGVO[]; xh = MGXH[]; gH = MGGH[]; mix = θ < 1
    @inbounds for g in 1:MGN[]
        eg, fg, en, fn = M.eg[g], M.fg[g], M.egn[g], M.fgn[g]; EO = MGEO[][g]; FO = MGFO[][g]
        for k in ngrs-1:ngre+2
            eg[k] = θ * en[k] + (1 - θ) * EO[k]; fg[k] = θ * fn[k] + (1 - θ) * FO[k]
        end
        chi = M.chif[g]; chn = M.chifn[g]; ke = M.xke[g]; em = M.emit[g]; j = mgv(g, 1)
        for k in ngrs-1:ngre+1
            m = mg_zcol(W.d, W.dvol, W.rmu, k)
            keN, chN, kpN, a1N = mg_level(gv, j, k, W.chifn[k], W.xken[k], W.xkpn[k], m); bn = mg_bfrac(g, tn[k])
            if mix
                keO, chO, kpO, a1O = mg_level(go, j, k, c.chifo[k], c.xkeo[k], c.xkpo[k], m)
                ke[k] = θ * keN + (1 - θ) * keO; chi[k] = θ * chN + (1 - θ) * chO
                e_ = (θ * kpN + (1 - θ) * kpO) * (θ * W.plfn[k] * bn + (1 - θ) * c.plfo[k] * mg_bfrac(g, to[k]))
                if g >= gH && (val(a1N) > 0 || a1O > 0)
                    e_ += (θ * a1N + (1 - θ) * a1O) * (θ * W.plfn[k] * mg_xfrac(g, tn[k], xh) + (1 - θ) * c.plfo[k] * mg_xfrac(g, to[k], xh))
                end
            else
                ke[k] = keN; chi[k] = chN; e_ = kpN * W.plfn[k] * bn
                (g >= gH && val(a1N) > 0) && (e_ += a1N * W.plfn[k] * mg_xfrac(g, tn[k], xh))
            end
            chn[k] = chN; em[k] = 4P.cpi * e_
        end
        if g == MGH2G[]                              # H₂ Lyman–Werner absorption (new level; MGH2FN): extinction κ_H2, to the gas κ_H2,gas
            kh = M.kh2; kg = M.kh2g
            for k in ngrs:ngre; ke[k] += kg[k]; chi[k] += kh[k]; chn[k] += kh[k]; end
        end
        av = M.avchi[g]
        for k in ngrs:ngre+1; av[k] = 2 / (1 / chi[k-1] + 1 / chi[k]); end
        SLABFLUX[] && sf_avchi!(av, chi, ngrs, ngre)
        SUBZONE[] && sz_avchi!(av, chi, W.d, W.dvol, W.rmu, ngrs, ngre)
        qe, qeo, qF, qFo = M.qe[g], M.qeo[g], M.qF[g], M.qFo[g]
        for k in ngrs-1:ngre+2
            qe[k] = en[k] / dn[k]; qeo[k] = EO[k] / dO[k]
            qF[k] = fn[k] * 2 / (dn[k] + dn[k-1]); qFo[k] = FO[k] * 2 / (dO[k] + dO[k-1])
        end
        advect!(M.erb[g], M.ersn[g], qe, qeo, MGESO[][g], fle, θ, P.cadv, P.eps, ngrs, ngre, true, 30 + g)        # (own branch-decision ids per group: a shared id replays the last group's limiter/upwind decisions in the others, docs §23.2)
        advect!(M.frb[g], M.frsn[g], qF, qFo, MGFSO[][g], flu, θ, P.cadv, P.eps, ngrs, ngre, false, 40 + g)
    end
    fc = M.frc
    @inbounds for k in ngrs:ngre+1
        s = zero(T)
        for g in 1:MGN[]; s += M.avchi[g][k] * M.fg[g][k]; end
        fc[k] = s
    end
    xc = M.xch; cc = P.cc; dt = P.dt
    @inbounds for k in ngrs:ngre
        s = zero(T)
        for g in 1:MGN[]; ke = M.xke[g][k]; s += excap(cc, ke, W.d[k], dt) * (cc * ke * M.eg[g][k] - M.emit[g][k]); end
        xc[k] = s
    end
    for g in 1:MGN[]; fill!(M.xfer[g], zero(T)); end
    @inbounds for (p, pr_) in enumerate(MGRAMP[])                    # Raman pairs: group-to-group transfer (see MGRAMP)
        gL, gp = pr_[1], pr_[2]; jc = mgvp(p, 1); jw = jc + 1
        xl = M.xfer[gL]; xp = M.xfer[gp]; eL = M.eg[gL]; ep = M.eg[gp]
        for k in ngrs:ngre
            cθ = mix ? θ * gv[jc, k] + (1 - θ) * go[jc, k] : gv[jc, k]; wθ = mix ? θ * gv[jw, k] + (1 - θ) * go[jw, k] : gv[jw, k]
            (val(cθ) > 0 || val(wθ) > 0) || continue
            t = cc * (cθ * ep[k] - wθ * eL[k]); xl[k] += t; xp[k] -= t
        end
    end
    if MGH2G[] > 0                                   # H₂ fluorescence: the absorbed LW energy that does not heat the gas, into group MGH2T
        xl = M.xfer[MGH2G[]]; xt = M.xfer[MGH2T[]]; eL = M.eg[MGH2G[]]; kh = M.kh2; kg = M.kh2g
        @inbounds for k in ngrs:ngre
            t = cc * (kh[k] - kg[k]) * eL[k]; xt[k] += t; xl[k] -= t
        end
    end
    nothing
end

"""
(MGROUPS, Jacobian-vector products of NKRYLOV) The groups' Eddington factors f_g, the tie's h_g and weights, the outer g_g, and the grey f of the total rows
(Σ_g f_g E_g/Σ_g E_g) as functions of the iterate: the formal solution of each group in the residual's Duals (TitanCore.eddfac_dual), on exactly the
arrays `mg_eddfac!` (popiii_m2c) swaps into the grey formal solution — χ_g, κ_E,g and the source (emission + Raman + H₂ fluorescence)/κ_E,g of the zone's
group record at the iterate, E_g, F_g, I_cmb b_g(T_cmb). The same function the update! writes into MGFEDD, … — the product then differentiates what the
finite difference through update! sees.
"""
function mg_eddfac_dual!(M::MGWS{T}, c, P, W, a, b, rn, dn, tn, dvoln, cplz) where {T}
    G = MGN[]; V = M.gv; c4pi = P.cc / (4P.cpi); mu = c.mu
    rmun(k) = rn[k]^mu
    gH2 = MGH2G[]
    for g in 1:G
        j = mgv(g, 1); fg_ = M.fdd[g]; hg_ = M.hde[g]; tg_ = M.tiw[g]
        @inbounds for k in eachindex(fg_); fg_[k] = T(MGFEDD[][g][k]); hg_[k] = T(MGHEDD[][g][k]); tg_[k] = T(MGTIEW[][g][k]); end
        chi = M.sc1; kp = M.sc2; pl = M.sc3
        @inbounds for k in a-1:b+1
            m = (a_ = 0.5 * (rmun(k) + rmun(k + 1)); (val(a_) > 0 && val(dvoln[k]) > 0) ? dn[k] * dvoln[k] / a_ : zero(dn[k] * dvoln[k]))      # (mg_zcol on the new level)
            kE, chg, kP, a1 = mg_level(V, j, k, W.chifn[k], W.xken[k], W.xkpn[k], m)
            kE = max(kE, 1e-300); em = kP * W.plfn[k] * mg_bfrac(g, tn[k])
            (g >= MGGH[] && val(a1) > 0) && (em += a1 * W.plfn[k] * mg_xfrac(g, tn[k], MGXH[]))
            for (p_, pr_) in enumerate(MGRAMP[])
                kc = V[mgvp(p_, 1), k]; kw = V[mgvp(p_, 2), k]; (val(kc) > 0 || val(kw) > 0) || continue
                pr_[1] == g && (kE += kw; em += kc * c4pi * M.egn[pr_[2]][k])
                pr_[2] == g && (kE += kc; em += kw * c4pi * M.egn[pr_[1]][k])
            end
            if gH2 > 0 && val(M.kh2[k]) > 0
                g == gH2 && (kE += M.kh2[k]; chg += M.kh2[k])
                g == MGH2T[] && (em += (M.kh2[k] - M.kh2g[k]) * c4pi * M.egn[gH2][k])
            end
            chi[k] = chg; kp[k] = kE; pl[k] = em / kE
        end
        fe, he, hx, fx, tw, gr = EDDFAC_DUAL[](c, P, a, b, chi, kp, pl, dn, tn, M.egn[g], M.fgn[g], rn, W.xs1, W.xs2, T, DBUF, SLQBUF, nothing, false, cplz; ximr = P.ximr * MGBCMB_[g])
        @inbounds for k in a:b; fg_[k] = fe[k]; hg_[k] = he[k]; end
        @inbounds for k in a+1:b; tg_[k] = tw[k]; end
        M.gdr[g] = gr
    end
    fw = M.fwt
    @inbounds for k in a-1:b+2
        sE = zero(T); sF = zero(T)
        for g in 1:G; e_ = M.egn[g][k]; sE += e_; sF += M.fdd[g][k] * e_; end
        fw[k] = val(sE) > 0 ? sF / sE : T(c.fedd[k])
    end
    fw
end
"VEF_SPH outer-boundary energy density of a group (vef_esurf with the group's E, F, χ and background E_b,g)."
@inline function mg_esurf(rn, en, chn, dn, fn, fedd, N, xmu, cc, Eb)
    rc = 0.5 * (rn[N] + rn[N+1]); Ex = en[N] - Eb
    # f_x = (f E − E_b/3)/(E − E_b), evaluated in units of σ = max(|E|, |E_b|) (a constant of the evaluation: the same function). In Duals the quotient rule divides by (E − E_b)², which underflows to 0
    # for an empty group at the outer boundary (E ~ 1e-208: (E − E_b)² ~ 1e-416) and gave the Inf/NaN entries of the F_g row's Jacobian there (docs §23.2).
    σ = max(abs(ForwardDiff.value(en[N])), abs(ForwardDiff.value(Eb)))
    fx = (ForwardDiff.value(Ex) != 0 && σ > 0) ? clamp(((fedd[N] * en[N] - Eb / 3) / σ) / (Ex / σ), 1 / 3, 1.0) : 1 / 3        # (a Dual: the clamp's active branch carries the derivative, none at its ends)
    s = 0.5xmu * (3 - 1 / fx) * log(rn[N+1] / rc)
    Eb + exp(-s) * Ex - exp(-0.5s) * chn[N] * dn[N] * fn[N+1] * (rn[N+1] - rc) / (cc * fx)
end

"""
The radiation rows with MGROUPS: the grey rows ie/iF become identities in the (redundant) correction of E_r, F_r (interior zones / interfaces and the outer boundary
flux; E_r, F_r read as the sums of the groups in the whole residual, `mg_esum!`), and the group rows (energy per zone, momentum per interface, boundary and ghost rows as the grey ones).
"""
function mg_rows!(R, M::MGWS{T}, c, W, x, kb, P, dt, θ, xmu, cc, sph, tie) where {T}
    ngrs, ngre, neqn, k0, nb = c.ngrs, c.ngre, c.neqn, c.k0, c.nb
    ie, iF = c.ix[6], c.ix[7]
    X(i, k) = x[(kb(k) - 1) * neqn + i]
    ern, frn, rn, dn, dvoln, u, r, rmu, d, dvol, dmdt, urel = W.ern, W.frn, W.rn, W.dn, W.dvoln, W.u, W.r, W.rmu, W.d, W.dvol, W.dmdt, W.urel
    G = MGN[]
    dts = MGSTEADY[] ? Inf : dt                     # the group rows' own time derivatives (E_g, the flux inertia, the grid term of the tie): 0 in the steady solve of the initial state (`MGSTEADY`)
    # grey rows → sums
    # (identities, boundary and ghost rows included: E_r, F_r are the sums of the group unknowns, see `mg_esum!`; their grey ghost rows E_r(k±1) − E_r(k), … would read the sums too, leave the ghost
    #  unknowns of E_r, F_r in no row, and the band matrix singular)
    @inbounds for k in k0:ngre+2; R[ie, kb(k)] = X(ie, k); R[iF, kb(k)] = X(iF, k); end
    Tcmb = MGTCMB[]
    @inbounds for g in 1:G
        iE = mgiE(g); iFg = mgiF(g)
        en, fn, eg, fg = M.egn[g], M.fgn[g], M.eg[g], M.fg[g]; EO = MGEO[][g]; FO = MGFO[][g]
        ke = M.xke[g]; em = M.emit[g]; erb = M.erb[g]; frb = M.frb[g]; av = M.avchi[g]; qe_ = M.qe[g]; qeo_ = M.qeo[g]
        qg = g == MGQG[]
        vef = MGVEFOK[]
        if vef && AD_JV[] && T <: ForwardDiff.Dual      # Jacobian-vector product: the group's formal solution in Duals (mg_eddfac_dual!)
            fedd = M.fdd[g]; HE = M.hde[g]; TW = M.tiw[g]; gr = M.gdr[g]
        else
            fedd = vef ? MGFEDD[][g] : c.fedd; HE = vef ? MGHEDD[][g] : TD.HEDD; TW = vef ? MGTIEW[][g] : TD.TIEW; gr = vef ? MGGEDDR[][g] : P.geddr
        end
        # ---- energy
        for k in ngrs:ngre
            aur = 0.5 * (u[k+1] / r[k+1] + uor(u, r, k, ngrs))
            w = (en[k] * dvoln[k] - EO[k] * c.dvolo[k]) / dts + rmu[k+1] * fg[k+1] - rmu[k] * fg[k] +
                fedd[k] * eg[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) +
                (1 - 3fedd[k]) * 0.5 * xmu * eg[k] * aur * dvol[k] +
                excap(cc, ke[k], d[k], dt) * (cc * ke[k] * eg[k] - em[k]) * d[k] * dvol[k] - M.xfer[g][k] * d[k] * dvol[k]
            if qg
                isempty(RAYQ[]) || (w += RAYQ[][k] * dvol[k])
                isempty(GASQ[]) || (w += GASQ[][k] * dvol[k])
                w += (cqheat(k, W) + yqheat(k, W)) * dvol[k]
            end
            # the coupled zones' line cooling (H₂ rovibrational lines, HD: LINEIMP at the iterate) is emitted into group MGLG, with the
            # gas row's time weights: what the gas loses in lines the radiation field gains, in every zone (a gas row that holds
            # gas + radiation — the total row, the lrad row without GASROW — takes the same term back, see the gas row)
            (g == MGLG[] && cplz(k)) && (θc = max(θ, CPL_THETA[]); w -= (θc * M.lcl[k] + (1 - θc) * MGLCO[][k]) * dvol[k])
            if k == ngrs
                w += -dmdt[k+1] * erb[k+1] - urel[k] * rmu[k] * eg[k]
            elseif k == ngre
                w += dmdt[k] * erb[k]
                w += P.leobc != 3 ? urel[k+1] * rmu[k+1] * eg[k] : -dmdt[k+1] * erb[k+1]
            else
                w += -dmdt[k+1] * erb[k+1] + dmdt[k] * erb[k]
            end
            R[iE, kb(k)] = w
            if T === Float64 && !isempty(MGROWPROBE[]) && (g, k - ngrs + 1) in MGROWPROBE[]      # (diagnostic, MGROWPROBE="g,zone;g,zone") the terms of the group energy row
                tt = (en[k] * dvoln[k] - EO[k] * c.dvolo[k]) / dt; dv = rmu[k+1] * fg[k+1] - rmu[k] * fg[k]
                ab = excap(cc, ke[k], d[k], dt) * cc * ke[k] * eg[k] * d[k] * dvol[k]; ei = excap(cc, ke[k], d[k], dt) * em[k] * d[k] * dvol[k]
                println("MGROW g ", g, " zone ", k - ngrs + 1, " dt ", dt, " E ", en[k], " Eo ", EO[k], " s ", MGSC[][g][k], " time ", tt, " div ", dv, " adv ", fedd[k] * eg[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) + (1 - 3fedd[k]) * 0.5 * xmu * eg[k] * aur * dvol[k],
                        " absorb ", ab, " emit ", ei, " xfer ", M.xfer[g][k] * d[k] * dvol[k], " RAYQ ", isempty(RAYQ[]) ? 0.0 : RAYQ[][k] * dvol[k], " GASQ ", isempty(GASQ[]) ? 0.0 : GASQ[][k] * dvol[k],
                        " heat ", qg ? (cqheat(k, W) + yqheat(k, W)) * dvol[k] : 0.0, " line ", (g == MGLG[] && cplz(k)) ? (max(θ, CPL_THETA[]) * M.lcl[k] + (1 - max(θ, CPL_THETA[])) * MGLCO[][k]) * dvol[k] : 0.0,
                        " erb ", erb[k], " ke ", ke[k], " emit_density ", em[k], " T ", W.tn[k], " rho ", d[k], " dvol ", dvol[k], " kH2 ", g == MGH2G[] ? M.kh2g[k] : 0.0, " res ", w)
                println("MGROWN g ", g, " zone ", k - ngrs + 1, " dmdt ", (dmdt[k], dmdt[k+1]), " erb ", (erb[k], erb[k+1]), " ersn ", (M.ersn[g][k], M.ersn[g][k+1]), " E(k-2:k+2) ", [en[j] for j in k-2:k+2], " EO ", [EO[j] for j in k-2:k+2],
                        " rho ", [dn[j] for j in k-2:k+2], " rhoO ", [c.do_[j] for j in k-2:k+2], " slopeO ", [MGESO[][g][j] for j in k-2:k+2], " qe ", [qe_[j] for j in k-2:k+2], " qeo ", [qeo_[j] for j in k-2:k+2])
            end
        end
        R[iE, 1] = X(iE, k0)
        R[iE, kb(ngrs-1)] = en[ngrs] - en[ngrs-1]
        R[iE, kb(ngre+1)] = -en[ngre] + en[ngre+1]
        R[iE, kb(ngre+2)] = -en[ngre+1] + en[ngre+2]
        # ---- momentum
        qF, qFo = M.qF[g], M.qFo[g]; chi = M.chif[g]
        for k in ngrs+1:ngre
            inert = (qF[k] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) -
                     qFo[k] * (c.do_[k-1] * c.dvolo[k-1] + c.do_[k] * c.dvolo[k])) / (2cc^2 * dts) -
                    ((dmdt[k] + dmdt[k+1]) * frb[k] - (dmdt[k-1] + dmdt[k]) * frb[k-1]) / (2cc^2)
            drag = 0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1]) * av[k] * fg[k] / cc
            w = if !sph
                inert + rmu[k] * (fedd[k] * eg[k] - fedd[k-1] * eg[k-1]) + 0.5 * rmu[k] * (fg[k] / cc^2) * (u[k+1] - u[k-1]) +
                (0.25xmu / r[k]) * ((3fedd[k-1] - 1) * eg[k-1] * dvol[k-1] + (3fedd[k] - 1) * eg[k] * dvol[k]) + drag
            else
                qs = vefq(fedd, r, xmu, k)
                inert + rmu[k] * (qs * fedd[k] * eg[k] - fedd[k-1] * eg[k-1] / qs) + 0.5 * rmu[k] * (fg[k] / cc^2) * (u[k+1] - u[k-1]) + drag
            end
            wt = tie ? TW[k] : 0.0
            if wt > 0
                rcm = (0.5 * (r[k-1] + r[k]))^c.mu; rcp = (0.5 * (r[k] + r[k+1]))^c.mu
                hd = SLABTIE[] ? st_hdonor(max(HE[k-1], 0.0), fedd[k-1], chi[k-1], chi[k], 0.5 * d[k-1] * dvol[k-1] / rmu[k], 0.5 * d[k] * dvol[k] / rmu[k]) :
                                 max(HE[k-1], 0.0)
                B = rmu[k] * fg[k] / cc - (rcm * hd * eg[k-1] + rcp * min(HE[k], 0.0) * eg[k])
                if TIEFRAME[]                        # F_g is the comoving flux, as the grey row (TitanResidFast.TIEFRAME), per group
                    hi = min(HE[k], 0.0)
                    egup = decide(:up, 21 + g, k, val(rn[k]) >= c.ro[k]) ? eg[k] : eg[k-1]
                    B += (-dmdt[k] * erb[k] + (rn[k]^(c.mu + 1) - c.ro[k]^(c.mu + 1)) / ((c.mu + 1) * dts) * egup +
                          rmu[k] * u[k] * 0.5 * (fedd[k-1] * eg[k-1] + fedd[k] * eg[k]) -
                          2u[k] * (rcm * hd^2 * eg[k-1] + rcp * hi^2 * eg[k])) / cc
                end
                w = (1 - wt) * w + wt * B
            end
            R[iFg, kb(k)] = w
        end
        R[iFg, kb(ngrs)] = X(iFg, ngrs)
        bcmb = MGBCMB_[g]; ximg = P.ximr * bcmb
        R[iFg, kb(ngre+1)] = P.cpi * (2gr + 1) * ximg + fn[ngre+1] -
                             cc * gr * (sph ? mg_esurf(rn, en, M.chifn[g], dn, fn, fedd, ngre, xmu, cc, 4P.cpi * ximg / cc) : en[ngre])
        if TIEFRAME[]                                # the boundary condition holds for the lab-frame flux
            ub_ = W.un[ngre+1]
            R[iFg, kb(ngre+1)] = (R[iFg, kb(ngre+1)] - fn[ngre+1]) * (1 + 2gr * ub_ / cc) + fn[ngre+1] + (1 + fedd[ngre]) * ub_ * en[ngre]
        end
        R[iFg, 1] = X(iFg, k0)
        k = ngrs - 1
        R[iFg, kb(k)] = P.lribc == 2 ? fn[k+2] + fn[k] : rn[k+1]^c.mu * fn[k+1] - rn[k]^c.mu * fn[k]
        R[iFg, kb(ngre+2)] = -fn[ngre+1] + fn[ngre+2]
    end
    nothing
end
"b_g(T_cmb) of the outer background (set with the edges)."
const MGBCMB_ = Float64[]
"(diagnostic) group-row probes (group, zone) from ENV MGROWPROBE=\"g,zone;g,zone\""
const MGROWPROBE = Ref(Tuple{Int,Int}[])

"""
Set up MGROUPS: edges [eV] (interior, ascending), the eV/k_B constant, T_cmb, the unknown offset i0 (9 with COUPLED). Allocates the arrays.
"""
function mg_setup!(edges, evk, tcmb, i0, mgr)
    issorted(edges) && all(>(0), edges) || error("MGROUPS: edges must be positive and ascending (eV): $edges")
    MGROWPROBE[] = [Tuple(parse.(Int, split(q, ','))) for q in split(get(ENV, "MGROWPROBE", ""), ';') if !isempty(q)]
    MG[] = true; MGN[] = length(edges) + 1; MGEDGE[] = collect(Float64, edges); MGEVK[] = evk; MGTCMB[] = tcmb; MGI0[] = i0
    empty!(MGBCMB_); for g in 1:MGN[]; push!(MGBCMB_, tcmb > 0 ? mg_bfrac(g, tcmb) : (g == 1 ? 1.0 : MGBFLOOR)); end
    mg_alloc!(mgr); empty!(MGWSC)
    nothing
end

"""
Initial group fields from the grey field (a restart from a grey checkpoint, or the IC), the same in every zone and without a threshold: the local emission balance of the group's own row,
E_g ∝ b_g = j_g/(c κ_E,g) (j_g = κ_P,g b_g(T) σT⁴/π + the Milne term, κ_E,g: the emissivity and absorption of the energy row, `mg_rows!`/`mg_scales!`, from the zone's group record MGVN, which must be filled;
the Raman and H₂ exchange with the partner groups are not in it: they have no partner field yet), rescaled so that Σ_g E_g = E_r. The flux F_g = F_r w_g/Σ w, w_g the mean of E_g of the two zones of the interface
(each group carries the share of the field it has on that face; the groups then sum to F_r). Ghost zones take the balance of the nearest zone, with their own E_r. (Docs §25.)
"""
function mg_init_from_grey!()
    a = gi(:ngrs); b = gi(:ngre); er = v(:ern); fr = v(:frn); G = MGN[]; gv = MGVN[]; cc = gf(:cc)
    tn = v(:tn); plf = v(:plfn); chi = v(:chifn); xke = v(:xken); xkp = v(:xkpn)
    bg = [zeros(b + 2) for _ in 1:G]
    for k in a-2:b+2
        kk = clamp(k, a, b)
        for g in 1:G
            keO, _, kpO, a1O = mg_level(gv, mgv(g, 1), kk, chi[kk], xke[kk], xkp[kk], 0.0)
            em = kpO * plf[kk] * mg_bfrac(g, tn[kk])
            (g >= MGGH[] && a1O > 0) && (em += a1O * plf[kk] * mg_xfrac(g, tn[kk], MGXH[]))
            bg[g][k] = em / (cc * max(keO, 1e-300))
        end
        s = sum(bg[g][k] for g in 1:G)
        (isfinite(s) && s > 0) || error("mg_init_from_grey!: no emission balance in zone $kk (sum $s)")
        for g in 1:G; MGEN[][g][k] = er[k] * bg[g][k] / s; end
    end
    for g in 1:G; fill!(MGFN[][g], 0.0); end
    for k in a+1:b+1
        w = [0.5 * (MGEN[][g][k-1] + MGEN[][g][k]) for g in 1:G]; s = sum(w)
        for g in 1:G; MGFN[][g][k] = s > 0 ? fr[k] * w[g] / s : 0.0; end
        MGFN[][1][k] += fr[k] - sum(MGFN[][g][k] for g in 1:G)          # (round-off only)
    end
    for k in a-2:b+2                                                    # exact sums (rows E_r = Σ E_g)
        s = sum(MGEN[][g][k] for g in 1:G); s > 0 && for g in 1:G; MGEN[][g][k] *= er[k] / s; end
    end
    for g in 1:G; mg_ghosts!(MGEN[][g], MGFN[][g], a, b); end
    for g in 1:G; copyto!(MGEO[][g], MGEN[][g]); copyto!(MGFO[][g], MGFN[][g]); end
    mg_newstep!()
    nothing
end

"(diagnostic) group luminosities L_g = 4π r² F_g at interface k (values of the new level)."
mg_lum(k) = [4π * v(:rn)[k]^2 * MGFN[][g][k] for g in 1:MGN[]]

"""
Restrict the band Jacobian and the residual in COMMON to the radiation unknowns (the steady solve of the initial group fields, docs §26): every row of a non-radiation unknown (1 … MGI0: the
state variables, E_r, F_r, the coupled species) becomes the identity with a zero right-hand side (its correction is 0: the gas state is held), and the columns of those unknowns are zeroed in the group
rows. What remains is the Jacobian of the group rows (E_g, F_g of all zones, ghosts included) with respect to the group unknowns, exactly the entries the full band matrix holds.
"""
function mg_restrict_band!()
    E = (v(:em2), v(:em1), v(:e00), v(:ep1), v(:ep2)); rhs = v(:rhs); i0 = MGI0[]; ne = gi(:neqn)
    k0 = gi(:ngrs) - 2; k1 = gi(:ngre) + 2
    @inbounds for k in k0:k1
        for M in E
            for j in 1:ne, i in 1:i0; M[i, j, k] = 0.0; end
            for j in 1:i0, i in i0+1:ne; M[i, j, k] = 0.0; end
        end
        for i in 1:i0; E[3][i, i, k] = 1.0; rhs[i, k] = 0.0; end
    end
    nothing
end
