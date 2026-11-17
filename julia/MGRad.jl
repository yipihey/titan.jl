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
const MGCLAMP = Ref(3.0)                        # clamp of a correction of ln E_g per Newton iteration
const MGEMIN = 1e-150                           # floor of E_g/E_r after an update (an underflowed group makes the Newton matrix singular)
const MGCFL = Ref(1e-3)                         # convergence weight floor: |δ ln E_g| counts with E_g/(E_g + MGCFL·E_r)
const MGEFL = Ref(1e-3)                         # TR-BDF2 error-norm floor: changes of a group below MGEFL·E_r do not limit Δt
# (Tried and reverted, t_mg17/t_mg18: E_g ← E_g(1 + δ) (Newton in E_g, clamp [1e-6 − 1, 1e3]) with MGCFL = 1e-12: the ionising group then
#  flip-flops between ~1e-22 and ~1e-49 of E_r at the photosphere from step to step and Δt fell to 446–2228 s (log update ~1e4 s);
#  docs/POPIII_MULTIGROUP.md §10.)
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
MGVEF (stage 2a, popiii_m2c `MGVEF=1`, the default with MGROUPS): every group's own Eddington factors from a formal solution with the group's
χ_g, κ_E,g, emission and background (popiii_m2c `mg_eddfac!`, run after the grey one wherever TITAN calls it: every Newton update with VEFIMP):
f_g per zone (MGFEDD), the tie's h_g (MGHEDD) and weights (MGTIEW, from the group's own optical depth), the outer g_g (MGGEDDR). The grey
array fedd then holds Σ_g f_g E_g / Σ_g E_g (the radiation pressure of the total rows, MLT, …). Frozen in the Jacobian, as the grey factors.
"""
const MGVEF = Ref(false); const MGVEFOK = Ref(false)
const MGFEDD = Ref(Vector{Vector{Float64}}()); const MGHEDD = Ref(Vector{Vector{Float64}}()); const MGTIEW = Ref(Vector{Vector{Float64}}())
const MGGEDDR = Ref(Float64[])
const MGHPLUS = Ref(Vector{Vector{Float64}}())   # (diagnostic) outgoing half moment H⁺ of each group's formal solution per interface
"(diagnostic) of the last residual: Σ_g E_g/E_r − 1 worst, the group luminosities through each interface are in MGEN/MGFN"
const MGNEWTON = Ref(0)
mgiE(g) = MGI0[] + 2g - 1; mgiF(g) = MGI0[] + 2g

"Allocate the group arrays for mgr zones (values kept where present)."
function mg_alloc!(mgr)
    G = MGN[]
    for R_ in (MGEN, MGFN, MGEO, MGFO, MGESO, MGFSO, MGFNOM)
        (length(R_[]) == G && all(length(a) == mgr for a in R_[])) || (R_[] = vvec(G, mgr))
    end
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
"Start of a step: old := new (fields and opacity records), slopes, nominal flux scales."
function mg_newstep!()
    ngrs = gi(:ngrs); ngre = gi(:ngre)
    size(MGVO[]) == size(MGVN[]) && copyto!(MGVO[], MGVN[])
    for g in 1:MGN[]
        copyto!(MGEO[][g], MGEN[][g]); copyto!(MGFO[][g], MGFN[][g])
        E = MGEN[][g]; F = MGFN[][g]; nom = MGFNOM[][g]
        @inbounds for k in ngrs-2:ngre+2
            nom[k] = max(abs(F[k]), 1e-3 * gf(:cc) * 0.5 * (E[max(k-1, 1)] + E[k]), 1e-300)
        end
    end
    mg_slopes!()
    nothing
end
"Reset after a failed iterate: new := old."
mg_reset!() = (for g in 1:MGN[]; copyto!(MGEN[][g], MGEO[][g]); copyto!(MGFN[][g], MGFO[][g]); end;
               size(MGVO[]) == size(MGVN[]) && copyto!(MGVN[], MGVO[]); nothing)
"Apply the Newton correction: E_g ← max(E_g e^δ, MGEMIN·E_r) (δ the correction of ln E_g, clamped to ±MGCLAMP in mg_dmax!), F_g ← F_g + δ·nom."
function mg_update!(rhs, a, b)
    er = v(:ern)
    @inbounds for g in 1:MGN[]
        E = MGEN[][g]; F = MGFN[][g]; nom = MGFNOM[][g]; iE = mgiE(g); iF = mgiF(g)
        for k in a:b; E[k] = max(E[k] * exp(rhs[iE, k]), MGEMIN * er[k]); end
        for k in a+1:b+1; F[k] += rhs[iF, k] * nom[k]; end
        mg_ghosts!(E, F, a, b)
    end
    nothing
end
"Ghost values of a group as its ghost rows: E copied, F antisymmetric at the centre (lribc = 2) and copied outside."
@inline function mg_ghosts!(E, F, a, b)
    E[a-1] = E[a]; E[b+1] = E[b]; E[b+2] = E[b+1]
    F[a] = 0.0; F[a-1] = gi(:lribc) == 2 ? -F[a+1] : F[a-1]; F[b+2] = F[b+1]
    nothing
end
"maxdel: clamp ln E_g corrections to ±MGCLAMP and return the convergence measure (|δ ln E_g| weighted by E_g/(E_g + MGCFL·E_r))."
function mg_dmax!(rhs, a, b)
    d = 0.0; cl = MGCLAMP[]; er = v(:ern)
    @inbounds for g in 1:MGN[]
        E = MGEN[][g]; iE = mgiE(g)
        for k in a:b
            δ = rhs[iE, k]
            isfinite(δ) || (rhs[iE, k] = 0.0; d = Inf; continue)
            rhs[iE, k] = clamp(δ, -cl, cl)
            w = E[k] / (E[k] + MGCFL[] * er[k])
            d = max(d, abs(rhs[iE, k]) * w)
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
end
const MGWSC = Dict{DataType,Any}()
mgws(::Type{T}, mgr) where {T} = get!(MGWSC, T) do
    G = MGN[]; z() = [zeros(T, mgr) for _ in 1:G]; gv = zeros(T, mg_nv(), mgr)
    for k in 1:mgr; mg_default!(view(gv, :, k)); end
    MGWS{T}(z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), zeros(T, mgr), zeros(T, mgr), z(), gv, zeros(T, mgr), zeros(T, mgr), zeros(T, mgr), ones(T, G + 1, mgr))
end::MGWS{T}
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
            en[k] = E[k] * exp(x[(bb-1)*neqn+iE]); fn[k] = F[k] + x[(bb-1)*neqn+iF] * nom[k]
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
        advect!(M.erb[g], M.ersn[g], qe, qeo, MGESO[][g], fle, θ, P.cadv, P.eps, ngrs, ngre, true, 0)
        advect!(M.frb[g], M.frsn[g], qF, qFo, MGFSO[][g], flu, θ, P.cadv, P.eps, ngrs, ngre, false, 0)
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

"VEF_SPH outer-boundary energy density of a group (vef_esurf with the group's E, F, χ and background E_b,g)."
@inline function mg_esurf(rn, en, chn, dn, fn, fedd, N, xmu, cc, Eb)
    rc = 0.5 * (rn[N] + rn[N+1]); Ex = en[N] - Eb
    fx = ForwardDiff.value(Ex) != 0 ? clamp(ForwardDiff.value((fedd[N] * en[N] - Eb / 3) / Ex), 1 / 3, 1.0) : 1 / 3
    s = 0.5xmu * (3 - 1 / fx) * log(rn[N+1] / rc)
    Eb + exp(-s) * Ex - exp(-0.5s) * chn[N] * dn[N] * fn[N+1] * (rn[N+1] - rc) / (cc * fx)
end

"""
The radiation rows with MGROUPS: the grey rows ie/iF become E_r − Σ E_g, F_r − Σ F_g (interior zones / interfaces and the outer boundary
flux), and the group rows (energy per zone, momentum per interface, boundary and ghost rows as the grey ones).
"""
function mg_rows!(R, M::MGWS{T}, c, W, x, kb, P, dt, θ, xmu, cc, sph, tie) where {T}
    ngrs, ngre, neqn, k0, nb = c.ngrs, c.ngre, c.neqn, c.k0, c.nb
    ie, iF = c.ix[6], c.ix[7]
    X(i, k) = x[(kb(k) - 1) * neqn + i]
    ern, frn, rn, dn, dvoln, u, r, rmu, d, dvol, dmdt, urel = W.ern, W.frn, W.rn, W.dn, W.dvoln, W.u, W.r, W.rmu, W.d, W.dvol, W.dmdt, W.urel
    G = MGN[]
    # grey rows → sums
    @inbounds for k in ngrs:ngre
        s = zero(T); for g in 1:G; s += M.egn[g][k]; end
        R[ie, kb(k)] = ern[k] - s
    end
    @inbounds for k in ngrs+1:ngre+1
        s = zero(T); for g in 1:G; s += M.fgn[g][k]; end
        R[iF, kb(k)] = frn[k] - s
    end
    Tcmb = MGTCMB[]
    @inbounds for g in 1:G
        iE = mgiE(g); iFg = mgiF(g)
        en, fn, eg, fg = M.egn[g], M.fgn[g], M.eg[g], M.fg[g]; EO = MGEO[][g]; FO = MGFO[][g]
        ke = M.xke[g]; em = M.emit[g]; erb = M.erb[g]; frb = M.frb[g]; av = M.avchi[g]
        qg = g == MGQG[]
        vef = MGVEF[] && MGVEFOK[]
        fedd = vef ? MGFEDD[][g] : c.fedd; HE = vef ? MGHEDD[][g] : TD.HEDD; TW = vef ? MGTIEW[][g] : TD.TIEW; gr = vef ? MGGEDDR[][g] : P.geddr
        # ---- energy
        for k in ngrs:ngre
            aur = 0.5 * (u[k+1] / r[k+1] + uor(u, r, k, ngrs))
            w = (en[k] * dvoln[k] - EO[k] * c.dvolo[k]) / dt + rmu[k+1] * fg[k+1] - rmu[k] * fg[k] +
                fedd[k] * eg[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) +
                (1 - 3fedd[k]) * 0.5 * xmu * eg[k] * aur * dvol[k] +
                excap(cc, ke[k], d[k], dt) * (cc * ke[k] * eg[k] - em[k]) * d[k] * dvol[k] - M.xfer[g][k] * d[k] * dvol[k]
            if qg
                isempty(RAYQ[]) || (w += RAYQ[][k] * dvol[k])
                isempty(GASQ[]) || (w += GASQ[][k] * dvol[k])
            end
            # the coupled zones' line cooling (H₂ rovibrational lines, HD: LINEIMP at the iterate) is emitted into group MGLG, with the
            # gas row's time weights: what the gas loses in lines the radiation field gains
            (g == MGLG[] && cplz(k) && GASROW[] && !(c.etotal && !c.useint[k])) && (θc = max(θ, CPL_THETA[]); w -= (θc * M.lcl[k] + (1 - θc) * MGLCO[][k]) * dvol[k])
            if k == ngrs
                w += -dmdt[k+1] * erb[k+1] - urel[k] * rmu[k] * eg[k]
            elseif k == ngre
                w += dmdt[k] * erb[k]
                w += P.leobc != 3 ? urel[k+1] * rmu[k+1] * eg[k] : -dmdt[k+1] * erb[k+1]
            else
                w += -dmdt[k+1] * erb[k+1] + dmdt[k] * erb[k]
            end
            R[iE, kb(k)] = w
        end
        R[iE, 1] = X(iE, k0)
        R[iE, kb(ngrs-1)] = en[ngrs] - en[ngrs-1]
        R[iE, kb(ngre+1)] = -en[ngre] + en[ngre+1]
        R[iE, kb(ngre+2)] = -en[ngre+1] + en[ngre+2]
        # ---- momentum
        qF, qFo = M.qF[g], M.qFo[g]; chi = M.chif[g]
        for k in ngrs+1:ngre
            inert = (qF[k] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) -
                     qFo[k] * (c.do_[k-1] * c.dvolo[k-1] + c.do_[k] * c.dvolo[k])) / (2cc^2 * dt) -
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
                w = (1 - wt) * w + wt * B
            end
            R[iFg, kb(k)] = w
        end
        R[iFg, kb(ngrs)] = X(iFg, ngrs)
        bcmb = MGBCMB_[g]; ximg = P.ximr * bcmb
        R[iFg, kb(ngre+1)] = P.cpi * (2gr + 1) * ximg + fn[ngre+1] -
                             cc * gr * (sph ? mg_esurf(rn, en, M.chifn[g], dn, fn, fedd, ngre, xmu, cc, 4P.cpi * ximg / cc) : en[ngre])
        R[iFg, 1] = X(iFg, k0)
        k = ngrs - 1
        R[iFg, kb(k)] = P.lribc == 2 ? fn[k+2] + fn[k] : rn[k+1]^c.mu * fn[k+1] - rn[k]^c.mu * fn[k]
        R[iFg, kb(ngre+2)] = -fn[ngre+1] + fn[ngre+2]
    end
    nothing
end
"b_g(T_cmb) of the outer background (set with the edges)."
const MGBCMB_ = Float64[]

"""
Set up MGROUPS: edges [eV] (interior, ascending), the eV/k_B constant, T_cmb, the unknown offset i0 (9 with COUPLED). Allocates the arrays.
"""
function mg_setup!(edges, evk, tcmb, i0, mgr)
    issorted(edges) && all(>(0), edges) || error("MGROUPS: edges must be positive and ascending (eV): $edges")
    MG[] = true; MGN[] = length(edges) + 1; MGEDGE[] = collect(Float64, edges); MGEVK[] = evk; MGTCMB[] = tcmb; MGI0[] = i0
    empty!(MGBCMB_); for g in 1:MGN[]; push!(MGBCMB_, tcmb > 0 ? mg_bfrac(g, tcmb) : (g == 1 ? 1.0 : MGBFLOOR)); end
    mg_alloc!(mgr); empty!(MGWSC)
    nothing
end

"""
Initial group fields from the grey field (a restart from a grey checkpoint, or the IC): the CMB part partitioned at T_cmb, the excess at the
shape temperature Ts[k] (the driver passes T_r inside the photosphere and the star's T_eff outside), F_g = F_r b_g(Ts at the interface).
"""
function mg_init_from_grey!(Ts)
    a = gi(:ngrs); b = gi(:ngre); er = v(:ern); fr = v(:frn)
    arad = 4gf(:csigr) / gf(:cc); Ec = MGTCMB[] > 0 ? arad * MGTCMB[]^4 : 0.0
    for g in 1:MGN[]
        E = MGEN[][g]; F = MGFN[][g]; fill!(E, 0.0); fill!(F, 0.0)
        for k in a-2:b+2
            kk = clamp(k, a, b)
            Ex = max(er[k] - Ec, 0.0)
            E[k] = max(MGBCMB_[g] * min(er[k], Ec) + mg_bfrac(g, Ts[kk]) * Ex, MGBFLOOR * er[k])
        end
        for k in a+1:b+1; F[k] = fr[k] * mg_bfrac(g, Ts[min(k, b)]); end
        mg_ghosts!(E, F, a, b)
    end
    # exact sums (the rows E_r = Σ E_g, F_r = Σ F_g): renormalise
    for k in a-2:b+2
        s = sum(MGEN[][g][k] for g in 1:MGN[]); s > 0 && for g in 1:MGN[]; MGEN[][g][k] *= er[k] / s; end
    end
    for k in a+1:b+1
        s = sum(MGFN[][g][k] for g in 1:MGN[])
        d = fr[k] - s; MGFN[][1][k] += d          # (the remainder, round-off, to the IR group)
    end
    for g in 1:MGN[]; mg_ghosts!(MGEN[][g], MGFN[][g], a, b); end
    for g in 1:MGN[]; copyto!(MGEO[][g], MGEN[][g]); copyto!(MGFO[][g], MGFN[][g]); end
    mg_newstep!()
    nothing
end

"(diagnostic) group luminosities L_g = 4π r² F_g at interface k (values of the new level)."
mg_lum(k) = [4π * v(:rn)[k]^2 * MGFN[][g][k] for g in 1:MGN[]]
