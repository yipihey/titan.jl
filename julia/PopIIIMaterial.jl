"""
    PopIIIMaterial

Population III project, M2c: metal-free gas with non-equilibrium chemistry as a TITAN material. Composition per zone
(`y[k]` = (x_H⁺, x_H⁻, x_H₂, x_H₂⁺), PrimordialChem), RHFC02 EOS (`eos_rhfc02`: e_g includes the chemical binding energy).
The hydro step (TITAN's Newton) sees the composition frozen; after each converged step `chem!` advances chemistry,
cooling and T together in every zone at fixed (ρ, e_g) over the step just taken (`PrimordialChem.step_coupled`,
backward Euler) — operator splitting, option (A) of M0.2. Needs `TitanResidFast.ZONE` (zone index of a material call).
Hydrodynamics only for now (no opacities: neqn = 5).
"""
module PopIIIMaterial

import ..TitanDriver, ..TitanResidFast, ..TitanMaterials, ..TitanCore, ..PrimordialChem, ..PopIIIEqEOS
using ForwardDiff, Printf
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const TM = TitanMaterials; const TC = TitanCore
const PC = PrimordialChem; const EQ = PopIIIEqEOS

mutable struct PrimordialGas <: TM.Material
    y::Vector{NTuple{4,Float64}}          # per zone index k (length mgr)
    rates::Function                       # T -> named tuple of rate coefficients
    cooling::Union{Nothing,Function}      # (ρ, T, y) -> Λ [erg cm⁻³ s⁻¹]
    ρeq::Float64                          # equilibrium chemistry above this density (RHFC02: 5e-11 g cm⁻³); Inf = never
    eqtable::Bool                         # true: equilibrium EOS from PopIIIEqEOS (Saha + SCVH95, RHFC02); false: Saha only
    Λext::Vector{Float64}                 # per zone k: extra cooling [erg cm⁻³ s⁻¹] (line transfer), evaluated at T = Tref
    Tref::Vector{Float64}                 # G7: over the next step Λext × (T/Tref)^slope, implicit in T (see ext_cooling)
    slope::Vector{Float64}                #     slope = d ln Λ/d ln T (0: Λext held fixed, explicit)
end
PrimordialGas(y, rates, cooling) = PrimordialGas(y, rates, cooling, 5e-11, false)
PrimordialGas(y, rates, cooling, ρeq) = PrimordialGas(y, rates, cooling, ρeq, false)
PrimordialGas(y, rates, cooling, ρeq, eqtable) = PrimordialGas(y, rates, cooling, ρeq, eqtable, zeros(length(y)), ones(length(y)), zeros(length(y)))

"""
G7: the transfer cooling Λ_k of zone k as a function of its own temperature over the step, Λ_k (T/T_ref)^s with
s = d ln Λ_k/d ln T from the transfer (all zones perturbed together: the local response including the coupling to the
neighbours). Inside the implicit chemistry–temperature step this makes the line cooling implicit (backward Euler in T,
linearised in ln T) instead of explicit with Δt ≤ 0.1 uρ/Λ.
"""
ext_cooling(m, k) = (Λx = m.Λext[k]; Λx == 0 ? nothing :
                     m.slope[k] == 0 ? ((T) -> Λx + zero(T)) : (Tr = m.Tref[k]; s = m.slope[k]; (T) -> Λx * (T / Tr)^s))
"""
LINEIMP (popiii_m2c; always on): the H₂ line cooling of a zone is a function of the zone's current state, not a number held from the last
transfer call. Each ALI transfer call (popiii_m2c `linecool`; the only line-transfer mode) leaves per zone the radiative rates of its last population
solve (β = 1 − Λ*, J̄ = J̄_fs − Λ* S, the accelerated-Λ-iteration split) and hc, the line power the grey continuum absorbed in the zone (4π χ_c ∫J̄′dν).
Between calls those are frozen and the level populations are re-solved in statistical equilibrium at the temperature, the collider densities (n_H, n_H₂, n_He, n_H⁺ of the
iterate's composition) and the H₂ density of the iterate:
    Λ(ρ, T, x) = n_H₂ Σ_lines hν [n_u(T, n) (A β + B_ul J̄) − n_l(T, n) B_lu J̄] − hc ρ/ρ_call     (H2Levels.netrad = collcool in s.e.)
so the cooling is implicit and exact in T, ρ and the composition, and the derivatives are those of this function (ForwardDiff through the population solve: the
populations are generic in the number type). The first term is the net collisional exchange of the levels; the second returns to the gas the line photons the
continuum absorbed (thermalised by κ_P as the energy rows treat the Planck-mean absorption), so Λ is the net divergence of the line flux.
Λ is the net exchange with the zone's radiation field: it vanishes where the levels are in equilibrium with J̄ and is heating below that, so no held rate can take a zone
below the radiation temperature; J̄ ≥ β B_ν(T_CMB) per line.
"""
const LINEIMP = Ref(false)
const LIMODEL = Ref{Vector{Any}}(Any[])          # per zone k: nothing or (; LV, β, J, nH, nH2, nHe, nHp, ρ, hc, ws): the rates and the state of the call
li_model(k) = (LINEIMP[] && k <= length(LIMODEL[])) ? LIMODEL[][k] : nothing
li_invalidate!() = (fill!(LIMODEL[], nothing); nothing)
"(LINEIMP) the zone's line model from the radiative rates of the ALI transfer call (see LINEIMP); hc: the line power per volume the continuum absorbs [erg cm⁻³ s⁻¹]."
function li_set!(k, LV, β, J, nH, nH2, nHe, nHp, ρ, hc)
    length(LIMODEL[]) < k && append!(LIMODEL[], Any[nothing for _ in length(LIMODEL[])+1:k])
    LIMODEL[][k] = (nH2 > 0 && isfinite(hc)) ? (; LV, β, J, nH = max(Float64(nH), 0.0), nH2 = Float64(nH2), nHe = max(Float64(nHe), 0.0), nHp = max(Float64(nHp), 0.0), ρ = Float64(ρ), hc = Float64(hc),
                                              ws = Main.H2Levels.PopWS(LV)) : nothing
    nothing
end
const LIBAD = Ref(0)
"""
(LINEIMP) the net radiative line cooling [erg cm⁻³ s⁻¹] of the model md at temperature T and colliders (n_H, n_H₂, n_He, n_H⁺) (the continuum's absorption is added in `li_cool`) Any of the arguments may be a Dual (no work arrays then).
"""
function li_lam(md, T, nH, nH2, nHe, nHp)
    HV = Main.H2Levels
    λ = HV.li_lam_eval(md.LV, md.β::Vector{Float64}, md.J::Vector{Float64}, T, nH, nH2, nHe, nHp)
    if !isfinite(ForwardDiff.value(λ)) && LIBAD[] < 20
        LIBAD[] += 1; Tm = max(T, 1.0)
        @printf("LINEIMP: non-finite line cooling %s at T %.4g (n_H %.3e n_H2 %.3e n_He %.3e n_H+ %.3e)\n", string(λ), ForwardDiff.value(Tm), ForwardDiff.value(nH), ForwardDiff.value(nH2),
                ForwardDiff.value(nHe), ForwardDiff.value(nHp)); flush(stdout)
    end
    λ
end
"""
(LINEIMP) the line cooling of the zone with model md at (ρ, T, composition y | nothing): colliders of the iterate's composition; without one (table zone) those of the call scaled by ρ.
Net of the continuum's absorption of line photons: the gas takes back the line power the grey continuum absorbed, md.hc (its κ_P ρ J̄′ at the call, ∝ ρ at fixed J̄′), so
that the cooling is the net divergence of the line flux, −∇·F = Λ_coll − Λ_cont; the line photons absorbed by κ_P are thermalised in the gas as the Planck-mean absorption is in the energy rows.
"""
function li_cool(md, ρ, T, y)
    y === nothing && return (li_lam(md, T, md.nH, md.nH2, md.nHe, md.nHp) - md.hc) * (ρ / md.ρ)
    nHt = ρ / CNn.μH(); comp(s_) = nHt * max(y[s_], zero(nHt))
    li_lam(md, T, comp(:H), comp(:H2), comp(:He), comp(:Hp)) - md.hc * (ρ / md.ρ)
end
"Cooling function for zone k: the material's own plus the transfer term."
function zone_cooling(m, k)
    md = li_model(k)
    if md !== nothing && CPLPREP[]                   # the coupled zones' line cooling is that of the iterate: cpl_g adds it
        c = m.cooling
    elseif md !== nothing                            # LINEIMP: the line cooling of the iterate (see LINEIMP)
        g = (ρ, T, y) -> li_cool(md, ρ, T, y)
        c = m.cooling === nothing ? g : ((ρ, T, y) -> m.cooling(ρ, T, y) + g(ρ, T, y))
    else
    f = ext_cooling(m, k)
    c = f === nothing ? m.cooling : m.cooling === nothing ? ((ρ, T, y) -> f(T)) : ((ρ, T, y) -> m.cooling(ρ, T, y) + f(T))
    end
    (HDCHEM[] && !CPLPREP[] && k <= length(HDX[]) && HDX[][k] > 0) || return c      # (the coupled zones' HD cooling is that of the iterate: cpl_g)
    xhd = HDX[][k]                                   # HDCHEM=1: HD cooling with the x_HD of the start of the step (operator-split chemistry)
    c === nothing ? ((ρ, T, y) -> hd_cool(ρ, T, y, xhd)) : ((ρ, T, y) -> c(ρ, T, y) + hd_cool(ρ, T, y, xhd))
end

@inline function TX.material_gas(m::PrimordialGas, ρ, T)
    if m.eqtable && ForwardDiff.value(ρ) >= m.ρeq
        u, p = EQ.eos(ρ, T)                 # composition in equilibrium at (ρ, T): differentiated through by the Newton
        return p, u, sqrt(5.0 / 3.0 * p / ρ)
    end
    y = m.y[TX.ZONE[]]
    u, p = PC.eos_rhfc02(ρ, T, y)
    # a_s enters only the artificial viscosity: γ ≈ 5/3 estimate
    return p, u, sqrt(5.0 / 3.0 * p / ρ)
end

function TM.fill_eos!(m::PrimordialGas)
    dn, tn, pgn, egn, asn = TM.v(:dn), TM.v(:tn), TM.v(:pgn), TM.v(:egn), TM.v(:asn)
    @inbounds for k in TM.krange()
        (dn[k] > 0 && tn[k] > 0) || continue
        TX.ZONE[] = k
        pg, eg, as_ = TX.material_gas(m, dn[k], tn[k])
        pgn[k] = pg; egn[k] = eg; asn[k] = as_
    end
    nothing
end
TM.fill_opac!(::PrimordialGas) = nothing

"Copy the composition of the first/last real zone into the ghost zones."
function ghosts!(m::PrimordialGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    for k in a-2:a-1; k >= 1 && (m.y[k] = m.y[a]); end
    for k in b+1:b+2; m.y[k] = m.y[b]; end
end

"""
After a converged step: chemistry + cooling + T in every zone over the step just taken, at fixed ρ and e_g; then the
EOS fill (p, e consistent with the new T and y). Returns the largest relative change of T (for diagnostics).
"""
function chem!(m::PrimordialGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    dn, tn, egn = TM.v(:dn), TM.v(:tn), TM.v(:egn)
    dTmax = 0.0
    for k in a:b
        cool = zone_cooling(m, k)
        y1, T1 = dn[k] >= m.ρeq ? (m.eqtable ? EQ.step(dn[k], egn[k], tn[k], dt; cooling = cool) :
                                             PC.step_equilibrium(dn[k], egn[k], tn[k], dt; cooling = cool)) :
                                  PC.step_coupled(m.y[k], dn[k], egn[k], tn[k], dt, m.rates; cooling = cool)
        dTmax = max(dTmax, abs(T1 / tn[k] - 1))
        m.y[k] = y1; tn[k] = T1
    end
    ghosts!(m)
    # ghost temperatures as TITAN's boundary rules (reflecting centre, copy outside), then the EOS fill
    tn[a-1] = tn[a]; tn[b+1] = tn[b]; tn[b+2] = tn[b]
    TM.fill_eos!(m)
    dTmax
end


# ====================================================================== G5: network gas (ChemNet) with an EOS-validity switch
"""
    NetGas(x, net, photos, cooling; env, pswitch = 1e4)

Metal-free gas with the general network (ChemNet/ChemRates; composition `x[k]` = free-species abundances per zone). The
switch to the equilibrium EOS table (PopIIIEqEOS: Saha + SCVH95) is made on **EOS validity, not density**: where the
ideal-gas pressure of the zone is ≥ `pswitch` (RHFC02: non-ideal effects from 1e4 dyn cm⁻², SCVH95 from 1e6 with a
blend between) — the network (detailed balance, implicit, L-stable) already relaxes to equilibrium by itself wherever
the chemistry is fast, so no ρ_eq is needed. `env` supplies J_ν for photo-processes (`ChemNet.photorates`); default
the CMB/zero field given by the caller.
"""
mutable struct NetGas <: TM.Material
    x::Vector{Vector{Float64}}           # per zone index k
    net::Any; photos::Any
    cooling::Union{Nothing,Function}     # (ρ, T, d::Dict) -> Λ
    env::Any                             # (Trad, J)
    pswitch::Float64
    eqzone::Vector{Bool}                 # per zone: in the equilibrium (table) regime after the last chem! step
    Λext::Vector{Float64}; Tref::Vector{Float64}; slope::Vector{Float64}     # as PrimordialGas (G7)
end
NetGas(x, net, photos, cooling; env, pswitch = 1e4) =
    NetGas(x, net, photos, cooling, env, pswitch, falses(length(x)), zeros(length(x)), ones(length(x)), zeros(length(x)))

import ..ChemNet
const CNn = ChemNet

"""
S6-C1: hot branch of `NetGas` for the protostar (popiii_m2c HOT=1). `HOT[] = (; T1, T2, Δ0, gas, thermo)` with
`gas(ρ, T, k) -> (p, e)` the stellar EOS of zone k's composition, `thermo(ρ, T, k) -> (∇_ad, c_P, δ)`. Between T1 and T2 p and e are
blended with a C² weight in ln T, e_gas = e_table + w (e_star − Δ0 − e_table), Δ0 the zero-point difference of the two EOS
(neutral atoms vs PopIIIEqEOS's zero), measured in the fully ionised dilute limit. The blend is not thermodynamically consistent
inside the band, but e(ρ, T) and p(ρ, T) are continuous and differentiable, so a zone crossing it conserves energy.
`nothing` (default): no change, the runs below T1 are identical.
"""
const HOT = Ref{Any}(nothing)
@inline function hotw(T, T1, T2)
    x = log(T / T1) / log(T2 / T1)
    x * x * x * (x * (6x - 15) + 10)
end
@inline function TX.material_gas(m::NetGas, ρ, T)
    return material_gas0(m, ρ, T)
end
@inline function material_gas0(m::NetGas, ρ, T)
    k = TX.ZONE[]
    h = HOT[]
    if h !== nothing && ForwardDiff.value(T) > h.T1
        ps, es = h.gas(ρ, T, k); es -= h.Δ0
        if ForwardDiff.value(T) >= h.T2
            return ps, es, sqrt(5.0 / 3.0 * ps / ρ)
        end
        u, p = EQ.eos(ρ, T); w = hotw(T, h.T1, h.T2)
        pb = p + w * (ps - p)
        return pb, u + w * (es - u), sqrt(5.0 / 3.0 * pb / ρ)
    end
    if m.eqzone[k]
        u, p = EQ.eos(ρ, T)
    elseif TX.cplz(k)
        # COUPLED: the network EOS of the CK expansion of the Newton iterate (inside the residual) or of the current species (outside)
        w = TX.CPLW[]
        y = w === nothing ? cpl_zone(m, k, ρ, T, TX.XSN[][1][k], TX.XSN[][2][k], TX.XSN[][3][k]) : cpl_zone(m, k, ρ, T, w[1][k], w[2][k], Float64(ForwardDiff.value(w[3][k])))
        p = y[1]; u = y[2]
    else
        u, p = CNn.eos(ρ, T, m.x[k], m.net)
    end
    return p, u, sqrt(5.0 / 3.0 * p / ρ)
end
"""
S6: nuclear source of the hot branch in the implicit energy row (S5's coupling: `h.eps(ρ, T, k)` = level factor × ε(ρ, T; Y lagged), generic in the
number type, so ∂ε/∂T is in the Jacobian; 0 below T1).
"""
TX.has_source(::NetGas) = HOT[] !== nothing
TX.material_source(::NetGas, ρ, T) = (h = HOT[]; ForwardDiff.value(T) > h.T1 ? h.eps(ρ, T, TX.ZONE[]) : zero(T))
"""
S6-C1d (experiment): MLT also on the table EOS, in equilibrium zones with T > `COLDMLT[]` (0 = off, the default). `cold_thermo` is ∇_ad, c_P, δ of
PopIIIEqEOS + equilibrium radiation, from derivatives of the table, generic in the number type.
"""
const COLDMLT = Ref(0.0)
function cold_thermo(ρ, T)
    arad = 7.565723e-15
    pf(r, t) = EQ.eos(r, t)[2] + arad * t^4 / 3
    ef(r, t) = EQ.eos(r, t)[1] + arad * t^4 / r
    p = pf(ρ, T)
    pT = ForwardDiff.derivative(t -> pf(ρ, t), T); pρ = ForwardDiff.derivative(r -> pf(r, T), ρ)
    eT = ForwardDiff.derivative(t -> ef(ρ, t), T); eρ = ForwardDiff.derivative(r -> ef(r, T), ρ)
    cP = eT - (eρ - p / ρ^2) * pT / pρ
    δ = T / ρ * pT / pρ
    return p * δ / (T * ρ * cP), cP, δ
end
"""MLT thermodynamics: only the hot branch convects (and, with `COLDMLT[]` > 0, the hot equilibrium zones of the table EOS). Below T1 the huge ∇_ad switches the flux off in the chemistry gas: it must be
out of reach (it was 10, and at a T jump inside the protostar ∇ = Δln T/Δln P crossed it: an MLT flux with this placeholder c_P
of up to 1e40 erg/s between two zones, and Newton diverged)."""
function TX.material_thermo(m::NetGas, ρ, T)
    h = HOT[]
    (h !== nothing && T > h.T1) && return h.thermo(ρ, T, TX.ZONE[])
    (COLDMLT[] > 0 && T > COLDMLT[] && m.eqzone[TX.ZONE[]]) && return cold_thermo(ρ, T)
    return (1e30, 1e9, 1.0)
end
# (PERF512) the zones are independent (each writes its own entries; the COUPLED zones their own CPLCACHE entry): on MATTHREADS tasks
# through TitanResidFast.zone_sweep!, bit-identical to the serial loop. The coupled zones' CK expansion here (every update!) was ~45 % of
# the hydro step of the 512-zone fixed mesh, serial.
function TM.fill_eos!(m::NetGas)
    dn, tn, pgn, egn, asn = TM.v(:dn), TM.v(:tn), TM.v(:pgn), TM.v(:egn), TM.v(:asn)
    TX.zone_sweep!(TM.krange(), k -> @inbounds(dn[k] > 0 && tn[k] > 0)) do k
        @inbounds begin
            pg, eg, as_ = TX.material_gas(m, dn[k], tn[k])
            pgn[k] = pg; egn[k] = eg; asn[k] = as_
        end
    end
    nothing
end
"""
S6: continuum opacity of `NetGas` when TITAN's radiation is on (popiii_m2c RAD=1): `OPAC[] = (ρ, T) -> (κ_R, κ_P)`
[cm² g⁻¹], e.g. StellarOpacity (own CIA/H⁻ means in chemical equilibrium below 10^3.8 K, OPAL above). `nothing` (default):
hydrodynamics only, no opacity calls. The H₂ lines stay outside the continuum (H2Lines, as Λext).
"""
const OPAC = Ref{Any}(nothing)
"""
Non-LTE opacity of the network zones of a three-level run (SL3=1, `NLOPAC[] = true`). `OPAC[]` is a table in chemical *and* excitation
equilibrium: its electrons are Saha's at the gas T and its H levels Boltzmann's. Outside the photosphere neither holds. At 36 M☉ just outside
τ = 1 (T ≈ 5700 K, ρ ≈ 2e-12) 83 % of κ_P is H bound–free (n = 1 and n ≥ 2 from Boltzmann populations ∝ e^{−10.2 eV/kT}), so κ_P ∝ T²²
while the network zone's heat capacity at frozen composition is 9× below the table's: Newton solved an emission ∝ T²⁶ there and diverged
every few steps (Δt 1e3–5e4 s, `_sl3b`). In the three-level atom n = 2 is set by recombination and the radiation field (4 n_e q21 ≪ β):
it hardly depends on the gas T at fixed composition, and its bound–free energy exchange is already in the chemistry (tl3_heat, GASQ/RAYQ).

So in every network zone with three-level parameters, from the end of `chem!` (T0, final composition), for the next hydro step
(TITAN's exchange term is c κ_E E − 4π κ_P B):
- κ_E (absorption): continuum Planck mean of the zone's own composition: the network's H⁻ (photodetachment, associative detachment),
  H⁻ ff, H ff, CIA, no H bound–free (case B: Lyman continuum in detailed balance; n ≥ 2 in the three-level atom);
- κ_P (emission): the same with H⁻ bound–free from H⁻ in Saha with the zone's own H and electrons (Kirchhoff–Milne for radiative
  attachment with Maxwellian electrons: emission ∝ n_H n_e, whatever n_H⁻ is). With the network's E0 the gas keeps −kT per attachment;
  the table, with equilibrium (molecular, electron-poor) gas, had no such emission in the cold infall while the network released 0.754 eV
  into heat per attachment;
- κ_R: Rosseland mean of the extinction with the network's H⁻ and H bound–free from the non-LTE n = 2 of the three-level atom (n ≥ 3
  with the same departure coefficient as n = 2: an approximation, these levels matter only through the Paschen continuum), Thomson on the
  zone's electrons;
- each as a ratio to the table at (ρ0, T0) below 6000 K (where the table is the same continuum code in equilibrium composition, agreement
  ≤ 3 %), above to the continuum code in equilibrium composition (so the table's calibration, OPAL and lines above 6300 K, is kept); then
  κ = κ_table(ρ, T0) × ratio × (T/T0)^q with q the frozen-composition slope of the continuum code, so Newton sees frozen composition in
  the EOS and in the opacity alike. No threshold in T: a T switch of κ made zones bistable (sawtooth 2900/4000 K at 3000–3500 K).
Deep in (τ ≳ 1, W = 1, T_c = T) the three-level n = 2 is Boltzmann's (b2 ≈ 0.8–1.1 at 36 M☉), the electrons are Saha's, and only the
dropped H bound–free of κ_P differs from the table zones: at zone 129 (τ 0.2) the three-level bound–free heating, 2.1e-5 erg cm⁻³ s⁻¹,
and the LTE κ_P,Hbf (B − J) it replaces, 5.2e-5, agree in sign and order; further out the diluted 6000 K starlight heats by n = 2
photoionization where the gray LTE term (T > T_r) cooled.
"""
const NLOPAC = Ref(false)
const NLOPC = Ref{Vector{Any}}(Any[])          # per zone: (; T0, RE, qE, RP, qP, RR, qR) or nothing
import ..ContinuumOpacity, ..G2Opacity, ..Z0Physics
const COp = ContinuumOpacity; const G2O = G2Opacity
const NLGRID = Ref{Vector{Float64}}(Float64[])
"""
Frequency-only quantities of the nlte_means grid (midpoint ν, dν, σ = ν/c, λ in μm, ν³, σ_H⁻bf(λ), Rayleigh σ of H₂ and H), computed once:
the same functions of the same arguments as in the loop, so the means are bit-identical; they were recomputed 2–3 × 1000 times per
network zone per step.
"""
const NLFREQ = Ref{Any}(nothing)
function nlfreq()
    NLFREQ[] === nothing || return NLFREQ[]
    isempty(NLGRID[]) && (NLGRID[] = COp.grid(; n = 1000))
    νs = NLGRID[]; n = length(νs) - 1
    ν = [0.5 * (νs[i] + νs[i+1]) for i in 1:n]; dν = [νs[i+1] - νs[i] for i in 1:n]; σ = [x / COp.c for x in ν]; λμ = [1e4 / x for x in σ]
    NLFREQ[] = (; ν, dν, σ, λμ, ν3 = [x^3 for x in ν], σhm = [COp.σ_Hm_bf(x) for x in λμ], ray2 = [COp.σ_ray_H2(x * 1e3) for x in λμ],
                ray1 = [COp.σ_ray_H(x * 1e4) for x in λμ])
end
"H bound–free absorption [cm⁻¹] (Kramers, as ContinuumOpacity.α_H_kramers) from n1 and the population of n = 2; n ≥ 3 Boltzmann relative to n = 2."
function hbf_nlte(ν, T, n1, n2)
    χ = 13.598 * EVERG; σ0 = 7.906e-18; hν = COp.h * ν; α = 0.0
    hν >= χ && (α += n1 * σ0 * (χ / hν)^3)
    for n in 2:30
        hν >= χ / n^2 || continue
        nn = n2 * (n^2 / 4) * exp(-χ * (0.25 - 1 / n^2) / (CNn.kB * T))
        α += nn * σ0 * n * (χ / n^2 / hν)^3
    end
    α
end
"The same with the n = 2…30 populations precomputed (`nn[n-1]`, independent of ν): bit-identical to `hbf_nlte`, without 29 exp per frequency."
function hbf_nlte_pre(ν, n1, nn)
    χ = 13.598 * EVERG; σ0 = 7.906e-18; hν = COp.h * ν; α = 0.0
    hν >= χ && (α += n1 * σ0 * (χ / hν)^3)
    for n in 2:30
        hν >= χ / n^2 || continue
        α += nn[n-1] * σ0 * n * (χ / n^2 / hν)^3
    end
    α
end
hbf_levels(T, n2) = ntuple(i -> (n = i + 1; n2 * (n^2 / 4) * exp(-13.598 * EVERG * (0.25 - 1 / n^2) / (CNn.kB * T))), 29)
"""Planck means of absorption (κ_E, H⁻ bf from nHm) and emission (κ_P, H⁻ bf from the same nHm, Kirchhoff), H bound–free in both or neither (`hbfP`), and
Rosseland mean of the extinction (H⁻ from nHm, H bound–free always) [cm² g⁻¹]; number densities in cm⁻³."""
function nlte_means(ρ, T; n1, n2, nH2, nHe, nHp, ne, nHm, hbfP, Tw = T, nHep = 0.0, nHepp = 0.0)
    G2O.ON && return G2O.means_any(ρ, T; n1, n2, nH2, nHe, nHep, nHepp, nHp, ne, nHm, hbfP, Tw)     # G2OPAC=1 (G2Opacity)
    F = nlfreq(); t = COp.tables(); sE = 0.0; sP = 0.0; sB = 0.0; sR = 0.0; sD = 0.0; sW = 0.0; wJ = Tw != T
    cia = nH2 > 1e-3 * n1; acia = 0.0; nnl = hbf_levels(T, n2)
    @inbounds for i in eachindex(F.ν)
        ν = F.ν[i]; dν = F.dν[i]; σ = F.σ[i]; λμ = F.λμ[i]
        stim = -expm1(-COp.h * ν / (COp.kB * T))
        # CIA bands are broad (≳ 100 cm⁻¹; the grid step is 1 %): evaluated on every 8th point, held in between (cost: 30 ms per call otherwise)
        cia && i % 8 == 1 && (acia = COp.k_h2h2(t, σ, T) * nH2^2 + COp.kcia(t.H2He, σ, T) * nH2 * nHe + COp.kcia(t.H2H, σ, T) * nH2 * n1)
        a0 = acia + COp.k_Hm_ff(λμ, T) * n1 * ne * COp.kB * T + 3.692e8 * ne * nHp / (sqrt(T) * F.ν3[i]) * stim
        σhm = F.σhm[i] * stim
        ab = hbf_nlte_pre(ν, n1, nnl) * stim
        sca = F.ray2[i] * nH2 + F.ray1[i] * n1 + COp.σT * ne
        B = COp.Bν(ν, T); D = COp.dBdT(ν, T); h = hbfP ? ab : 0.0
        BJ = wJ ? COp.Bν(ν, Tw) : B
        sE += (a0 + σhm * nHm + h) * BJ * dν; sW += BJ * dν; sP += (a0 + σhm * nHm + h) * B * dν; sB += B * dν
        sR += D / max(a0 + σhm * nHm + ab + sca, 1e-300) * dν; sD += D * dν
    end
    (sE / sW / ρ, sP / sB / ρ, sD / sR / ρ)
end
"n = 2 population [cm⁻³] of the three-level atom (2s + 2p, 4 n_2s) in steady state for the zone's n1, n_e, n_p and lagged radiation p3."
function n2_tl3(T, n1, ne, np, p)
    D = p.RL + LAM2S + 4ne * q21_H(T)
    4 * (n1 * (ne * q12_H(T) + p.P12) + alphaB_recfast(T) * ne * np) / (D + p.β)
end
"(κ_E, κ_P, κ_R) of the continuum code for zone composition d (full species fractions) or, with `d = nothing`, equilibrium composition."
function nlte_kappas(ρ, T, d, p)
    nH = ρ / CNn.μH(); nHe = PC.yHe * nH
    if d === nothing
        y = EQ.composition(ρ, T); xp = y[1]; xh2 = y[3]; n1 = nH * max(1 - xp - 2xh2, 0.0); ne = nH * xp
        if CNn.HECHEM                                 # helium electrons; neutral He for the CIA and Rayleigh terms
            z1, z2 = EQ.composition_he(ρ, T); ne = nH * (xp + z1 + 2z2); nHe = nH * max(PC.yHe - z1 - z2, 0.0)
        end
        n2 = 4n1 * exp(-E21 / (CNn.kB * T)); nHm = n1 * ne / CNn.Keq([:Hm], [:H, :e], T)
        G2O.ON && CNn.HECHEM && return nlte_means(ρ, T; n1, n2 = G2O.n2_lte(T, n1, ne, n1 + nHe + nH * xh2), nH2 = nH * xh2, nHe, nHp = nH * xp, ne, nHm,
                                                  hbfP = true, nHep = nH * z1, nHepp = nH * z2)
        return nlte_means(ρ, T; n1, n2, nH2 = nH * xh2, nHe, nHp = ne, ne, nHm, hbfP = true)
    end
    n1 = nH * max(d[:H], 0.0); ne = nH * max(d[:e], 0.0); np = nH * max(d[:Hp], 0.0)
    CNn.HECHEM && (nHe = nH * max(d[:He], 0.0))
    nHmS = n1 * ne / CNn.Keq([:Hm], [:H, :e], T)
    nHm = haskey(d, :Hm) ? nH * max(d[:Hm], 0.0) : nHmS
    Tw = NLOPACJ[] && hasproperty(p, :Tc) ? p.Tc : T
    G2O.ON && return hbf_blend(hb -> nlte_means(ρ, T; n1, n2 = n2_tl3(T, n1, ne, np, p), nH2 = nH * max(d[:H2], 0.0), nHe, nHp = np, ne, nHm, hbfP = hb, Tw,
                                nHep = nH * max(get(d, :Hep, 0.0), 0.0), nHepp = nH * max(get(d, :Hepp, 0.0), 0.0)), p)
    hbf_blend(hb -> nlte_means(ρ, T; n1, n2 = n2_tl3(T, n1, ne, np, p), nH2 = nH * max(d[:H2], 0.0), nHe, nHp = np, ne, nHm, hbfP = hb, Tw), p)
end
"""
`NLOHBFDEEP=1` (with TL3DEEP=1; default off, off = bit-identical): in network zones flagged `deep` (τ ≥ 3, TL3DEEP) the H bound–free
absorption enters κ_E and κ_P (hbfP = true) as in the table, while TL3DEEP drops the three-level atom's bound–free heats there: the
n = 2 photo-exchange is carried by the grey radiation coupling c κ_E E − 4π κ_P B, implicit in the hydro Newton, instead of the lagged
chemistry heat. For network zones deep in the thick layers (G2: the network through the photosphere front, docs/POPIII_G2_MATERIAL.md §6.4).
"""
const NLOHBFDEEP = get(ENV, "NLOHBFDEEP", "0") == "1"
"The blend weight is at an end of [0, 1] and is a Float (no derivative to carry): one mean suffices. A Dual weight at an end takes the blend (its slope there is the one-sided one)."
@inline hbf_pure(w, e) = e == 0 ? (w <= 0 && !(w isa ForwardDiff.Dual && w == 0)) : (w >= 1 && !(w isa ForwardDiff.Dual && w == 1))
hbf_w(p) = NLOHBFDEEP && p !== nothing && hasproperty(p, :deep) ? p.deep : 0.0
"The means f(hbfP) with the H bound–free absorption weighted by the deep weight of p (linear blend of the two means in the transition)."
@inline function hbf_blend(f::F, p) where {F}
    w = hbf_w(p)
    hbf_pure(w, 0) && return f(false)
    hbf_pure(w, 1) && return f(true)
    map((a, b) -> exp((1 - w) * log(max(a, 1e-300)) + w * log(max(b, 1e-300))), f(false), f(true))
end
"""
`G2BLEND=Γ0,D0,P0` (or `G2BLEND=1` for 0.3,0.1,0.01; default off, off = bit-identical): the network material everywhere outside the stellar
interior, the table only where the gas is non-ideal, by a smooth weight of the state at the Newton iterate instead of the switch rules
(PSWITCH, PSWNH, PSWTAU, PSWTEX, PSWSTAY are overridden: a zone is a table zone exactly where w ≥ 0.999, staying until w < 0.99, where the
two materials are the same function). Non-ideality s = ((Γ_e/Γ0)⁴ + (D/D0)⁴ + (P/P0)⁴)^{1/4} of the LTE (table) composition at (ρ, T):
  Γ_e = e²/(a_e kT), a_e = (3/(4π n_e))^{1/3}   electron Coulomb coupling;
  D = n_e λ_e³/2, λ_e = h/√(2π m_e kT)           electron degeneracy (≪ 1: Boltzmann electrons);
  P = (4π/3) n_neutral (2a₀)³                     neutral packing: the hard-sphere exponent of HM88's ground-state occupation probability
                                                  (Hummer & Mihalas 1988, eq. 4.71), the onset of pressure dissociation/ionisation.
w = smoothstep(log₁₀(s/s₀)/log₁₀ 3), s₀ = G2BLEND's 4th number (default 1/3): network for s ≤ 1/3, table for s ≥ 1 — the blend spans a
factor 3 in s and is complete where a parameter reaches its threshold (with s₀ = 1 the network's coupled chemistry ran into zones at
Γ_e ≈ 0.3–0.9 and failed to converge there, t_g2h, docs/POPIII_G2_MATERIAL.md §7.4). EOS (ChemNet.EOSMIX, so chemistry and hydro agree):
(u, p) = (1 − w) network + w table; opacity (NLOIT network zones): ln κ = (1 − w) ln κ_network + w ln κ_table.
Thresholds from the network-vs-SCVH95 comparison (docs/POPIII_G2_MATERIAL.md §7): for s < 1 the network (HM88, H2PF) and the table
agree to |Δ ln p| ≤ 0.045, |Δ ln Γ1| ≤ 0.064, |Δu| ≤ 0.2 kT per particle (95 %); beyond s ≈ 3 they differ by up to ×10 in p.
"""
const G2B = let v = get(ENV, "G2BLEND", "0"); v == "0" ? nothing : v == "1" ? (0.3, 0.1, 0.01, 1 / 3) : (t = parse.(Float64, split(v, ",")); length(t) == 3 ? (t..., 1 / 3) : Tuple(t)) end
const QE_ = 4.80320471e-10; const A0_ = 5.29177211e-9
"(G2BLEND) non-ideality s at (ρ, T) (generic: the table composition through its bicubic interpolant)."
function g2_s(ρ, T)
    t = EQ.table(); x = EQ.gx(t, log10(ρ)); y = EQ.gy(t, log10(min(T, 1e6) + zero(T)))
    xp = clamp(EQ.bicubic(t.xHp, x, y), 0.0, 1.0); f2 = clamp(EQ.bicubic(t.fH2, x, y), 0.0, 1.0)
    z1 = zero(xp); z2 = zero(xp)
    if CNn.HECHEM && EQ.HETAB[] !== nothing
        Z1, Z2 = EQ.HETAB[]; z1 = clamp(EQ.bicubic(Z1, x, y), 0.0, PC.yHe); z2 = clamp(EQ.bicubic(Z2, x, y), 0.0, PC.yHe)
    end
    nH = ρ / CNn.μH(); ne = nH * (xp + z1 + 2z2) + 1e-300
    nneu = nH * (max(1 - xp - f2, zero(xp)) + f2 / 2 + max(PC.yHe - z1 - z2, zero(xp)))
    kT = CNn.kB * T
    Γ = QE_^2 / (kT * cbrt(3 / (4π * ne))); D = ne * (CNn.hP / sqrt(2π * CNn.me * kT))^3 / 2; P = (4π / 3) * nneu * (2A0_)^3
    ((Γ / G2B[1])^4 + (D / G2B[2])^4 + (P / G2B[3])^4)^0.25
end
"(G2BLEND) blend weight of the table at (ρ, T)."
g2_w(ρ, T) = (s = g2_s(ρ, T); x = clamp(log10(s / G2B[4]) / log10(3.0), 0.0, 1.0); x * x * (3 - 2x))
"(G2BLEND) the EOS mix of ChemNet.EOSMIX."
g2_eosmix(ρ, T, u, p) = (w = g2_w(ρ, T); ForwardDiff.value(w) == 0 ? (u, p) : (ue = EQ.eos(ρ, T); (u + w * (ue[1] - u), p + w * (ue[2] - p))))
G2B !== nothing && (CNn.EOSMIX[] = g2_eosmix)
"Lagged non-LTE opacity factors of zone k at T0 (see NLOPAC)."
function nlopc(ρ, T0, d, p)
    OPAC[] === nothing && return nothing     # no opacity table yet (the cold collapse, before RAD's opacity is set): LTE tables later
    f = 1.02
    E0, P0, R0 = nlte_kappas(ρ, T0, d, p); E1, P1, R1 = nlte_kappas(ρ, T0 * f, d, p)
    if T0 < 6000
        Rt, Pt = OPAC[](ρ, T0); Pe, Re = Pt, Rt
    else
        _, Pe, Re = nlte_kappas(ρ, T0, nothing, p)
    end
    lf = log(f)
    (; T0, RE = E0 / Pe, qE = log(E1 / E0) / lf, RP = P0 / Pe, qP = log(P1 / P0) / lf, RR = R0 / Re, qR = log(R1 / R0) / lf)
end
"""
Diagnostic switch, off by default: `CHICAP=κ[,Tmax]` caps the flux-mean opacity χ at κ [cm²/g] where T < Tmax (default 5e4 K). The energy
means κ_E, κ_P are not changed. With κ = electron scattering it removes the H/H⁻ opacity peak from the radiation force and the diffusion,
to test what drives the envelope cycle at L ≈ L_Edd. Not physics.
"""
const CHICAP = let s = split(get(ENV, "CHICAP", "0"), ","); (parse(Float64, s[1]), length(s) > 1 ? parse(Float64, s[2]) : 5e4); end
@inline chicap(κ, T) = (CHICAP[1] > 0 && T < CHICAP[2]) ? min(κ, CHICAP[1]) : κ
"""
NLOIT=1 (with NLOPAC; default off, off = bit-identical): the opacities of a network zone are those of the Newton iterate's own
populations, not the lagged ratio of the last chemistry (NLOPC). For (ρ, T, x) of the iterate
    κ_i = OPAC_i(ρ, T) · R_i,   R_E = κ_E^nl(d)/κ_P^nl(eq),  R_P = κ_P^nl(d)/κ_P^nl(eq),  R_R = κ_R^nl(d)/κ_R^nl(eq)
with κ^nl the continuum code's means (`nlte_kappas`: the network's H⁻, H ff, CIA, Rayleigh, Thomson; H bound–free from the three-level
n = 2) at the iterate's T and ρ for the zone's composition d, and κ^nl(eq) the same code in equilibrium composition at the same (ρ, T),
above NLOITW[2] (default 6000 K); below NLOITW[1] (5000 K) the base is the table itself, i.e. κ = κ^nl(d) (as NLOPAC below 6000 K: the cold,
thin table is not the continuum code there — at 1600 K, ρ 2.5e-19 the table's κ_P is 3e-17 against ≲ 1e-20 of the code in equilibrium, so
the ratio form overstated κ by 10⁴–10⁵), geometric blend smooth in ln T between (NLOPAC switches the base at 6000 K).
d: in coupled zones the CK expansion of the iterate's x_H⁺, x_H₂ (`ck_eval` with the zone's terms of the step); in split zones the
composition the chemistry left (frozen over the hydro step, as in the EOS). The three-level radiation parameters (P12, RL, β, T_c) are
those of the step (lagged with the radiation field). The ratio form keeps the table's calibration (OPAL lines above ~6300 K) and goes
to the table where the populations are LTE, so it joins the table zones at τ = PSWTAU continuously; no jump in T anywhere.
χ ≥ σ_T n_e/ρ of the iterate's electrons (extinction is never below Thomson scattering on the zone's own electrons; this is NLOFIX's floor,
and what NLOTMAX's table fallback provided in ionized gas). Values exact at every iterate; Newton's derivatives in T (central), ρ, x_H⁺ and x_H₂ (one-sided)
differences of R (ln, relative steps 1e-3 and 1e-2), reused while T stays within NLOITJ[1] and x_H⁺ within NLOITJ[2] (relative, ln;
default 0, 0: new slopes at every new iterate; reusing them, e.g. NLOITJ=3e-3,0.1, made Newton stall at ~1e-4 corrections in the warm coupled zones of the production state, and the inexact Jacobian also enters DTERR's filtered estimate); ρ within NLOITJ[1].
Retires NLOIMPT, NLODAMP, CPLNLOX, NLOFIX, NLOTMAX in network zones (they have no effect with NLOIT).
"""
const NLOIT = true      # (the method since the clean-up; the lagged-factor path NLOIT=0 and its patches are gone)
const NLOITW = let s_ = split(get(ENV, "NLOITW", "5000,6000"), ","); (parse(Float64, s_[1]), parse(Float64, s_[2])) end
const NLOITJ = let s_ = split(get(ENV, "NLOITJ", "0,0"), ","); (parse(Float64, s_[1]), parse(Float64, s_[2])) end
const NLOST = Ref{Vector{Any}}(Any[])        # per zone: (; d, p) at the end of the last chemistry, or nothing
mutable struct NloC
    v::NTuple{5,Float64}; R::NTuple{3,Float64}; ne::Float64; bPR::NTuple{2,Float64}
    vJ::NTuple{5,Float64}; qT::NTuple{3,Float64}; qx::NTuple{3,Float64}; qr::NTuple{3,Float64}; qx2::NTuple{3,Float64}; qn::NTuple{4,Float64}; hasJ::Bool
    # MGROUPS: the zone's group record at v (MGOpac.jl mg_pack; empty without MGROUPS) and its logarithmic slopes, as R and qT, qx, qr, qx2
    g::Vector{Float64}; gT::Vector{Float64}; gx::Vector{Float64}; gr::Vector{Float64}; gx2::Vector{Float64}
end
const NLOC = Ref{Vector{Union{Nothing,NloC}}}(Union{Nothing,NloC}[])          # per zone cache
const NLOITN = Threads.Atomic{Int}(0); const NLOITNS = Threads.Atomic{Int}(0); const NLOCNS = Threads.Atomic{Int}(0); const NLOITJN = Threads.Atomic{Int}(0)     # (evaluations of R, of its slopes)
nloc_reset!(n) = (length(NLOC[]) == n ? fill!(NLOC[], nothing) : (NLOC[] = Union{Nothing,NloC}[nothing for _ in 1:n]); SLIMP && (NLOQC[] = Any[nothing for _ in 1:n]); nothing)
"(SLIMP) per zone the NLOIT ratio at the star-light state: (key (ρ, T, x1, x2, q, t, g, dp), R, n_e, key of the slopes, slopes in ln q, ln t, ln g, dp[, MGROUPS: group record, its slopes])."
const NLOQC = Ref{Vector{Any}}(Any[])
"CIA kernels on the nlte_means frequencies (every 8th point), tabulated on the union of the tables' node temperatures: exact (kcia and k_h2h2 are piecewise linear in T between those nodes, held outside)."
const CIATAB = Ref{Any}(nothing)
function ciatab()
    CIATAB[] === nothing || return CIATAB[]
    F = nlfreq(); t = COp.tables()
    idx = [i for i in eachindex(F.ν) if i % 8 == 1]
    Tn = sort(unique(vcat(t.H2H2.T, t.H2He.T, t.H2H.T, t.BOR.T, [3000.0])))
    K = zeros(3, length(Tn), length(idx))
    for (s, i) in enumerate(idx), (j, Tq) in enumerate(Tn)
        K[1, j, s] = COp.k_h2h2(t, F.σ[i], Tq); K[2, j, s] = COp.kcia(t.H2He, F.σ[i], Tq); K[3, j, s] = COp.kcia(t.H2H, F.σ[i], Tq)
    end
    CIATAB[] = (; Tn, K)
end
"nlte_means with the CIA kernels from `ciatab` (NLOIT); otherwise the same arithmetic."
nlte_means_fast(ρ, T; n1, n2, nH2, nHe, nHp, ne, nHm, hbfP, Tw = T, nHep = 0.0, nHepp = 0.0) = G2O.ON ?
    G2O.means_any(ρ, T; n1, n2, nH2, nHe, nHep, nHepp, nHp, ne, nHm, hbfP, Tw) : _nlmf(nlfreq(), ciatab(), Float64(ρ), Float64(T), Float64(n1), Float64(n2), Float64(nH2),
    Float64(nHe), Float64(nHp), Float64(ne), Float64(nHm), hbfP, Float64(Tw))
function _nlmf(F, ct, ρ, T, n1, n2, nH2, nHe, nHp, ne, nHm, hbfP, Tw)       # (function barrier: F, ct concrete)
    sE = 0.0; sP = 0.0; sB = 0.0; sR = 0.0; sD = 0.0; sW = 0.0; wJ = Tw != T
    cia = nH2 > 1e-3 * n1; acia = 0.0; nnl = hbf_levels(T, n2)
    Tn = ct.Tn; K = ct.K; nT = length(Tn)
    Tc = clamp(T, Tn[1], Tn[end]); jT = clamp(searchsortedlast(Tn, Tc), 1, nT - 1); wT = (Tc - Tn[jT]) / (Tn[jT+1] - Tn[jT])
    @inbounds for i in eachindex(F.ν)
        ν = F.ν[i]; dν = F.dν[i]; λμ = F.λμ[i]
        stim = -expm1(-COp.h * ν / (COp.kB * T))
        if cia && i % 8 == 1
            s = (i - 1) ÷ 8 + 1
            k1 = (1 - wT) * K[1, jT, s] + wT * K[1, jT+1, s]; k2 = (1 - wT) * K[2, jT, s] + wT * K[2, jT+1, s]; k3 = (1 - wT) * K[3, jT, s] + wT * K[3, jT+1, s]
            acia = k1 * nH2^2 + k2 * nH2 * nHe + k3 * nH2 * n1
        end
        a0 = acia + COp.k_Hm_ff(λμ, T) * n1 * ne * COp.kB * T + 3.692e8 * ne * nHp / (sqrt(T) * F.ν3[i]) * stim
        σhm = F.σhm[i] * stim
        ab = hbf_nlte_pre(ν, n1, nnl) * stim
        sca = F.ray2[i] * nH2 + F.ray1[i] * n1 + COp.σT * ne
        B = COp.Bν(ν, T); D = COp.dBdT(ν, T); h = hbfP ? ab : 0.0
        BJ = wJ ? COp.Bν(ν, Tw) : B
        sE += (a0 + σhm * nHm + h) * BJ * dν; sW += BJ * dν; sP += (a0 + σhm * nHm + h) * B * dν; sB += B * dν
        sR += D / max(a0 + σhm * nHm + ab + sca, 1e-300) * dν; sD += D * dν
    end
    (sE / sW / ρ, sP / sB / ρ, sD / sR / ρ)
end
"""
(NLOIT, restart) NLOST from the checkpoint: network zones that had non-LTE factors (NLOPC) take their composition and the three-level
parameters of SLDIAG. MGROUPS: the zones whose three-level record is the last chemistry's (it carries `chemq`: exactly the zones that
chemistry gave an NLOST entry) — NLOPC is no longer filled under NLOIT, so the NLOPC test restores nothing and the first step after a
restart runs the network zones on the table opacity with default group records (ionising and LW groups transparent: docs §18). The grey
path keeps the NLOPC test in this round (bit-identity with popiii-s6); it has the same first-step defect.
"""
function nlost_restart!(m)
    n = length(m.x); NLOST[] = Any[nothing for _ in 1:n]; nloc_reset!(n)
    sd = SLDIAG[]; (sd === nothing || !hasproperty(sd, :tl3)) && return nothing
    for k in 1:min(n, length(NLOPC[]), length(sd.tl3))
        (m.eqzone[k] || sd.tl3[k] === nothing || (MGON ? !hasproperty(sd.tl3[k], :chemq) : NLOPC[][k] === nothing)) && continue
        NLOST[][k] = (; d = CNn.full(m.x[k], m.net), p = sd.tl3[k])
    end
    nothing
end
"(NLOIT) equilibrium reference base (bP, bR) at (ρ, T) alone."
function nlo_base(ρ, T)
    nH = ρ / CNn.μH(); nHe0 = PC.yHe * nH
    Rt, Pt = OPAC[](ρ, T); w = nlo_w(T); Pq = Rq = 1.0
    if w > 0
        # equilibrium composition (as nlte_kappas with d = nothing)
        y = EQ.composition(ρ, min(T, 1e6)); xp = y[1]; xh2 = y[3]; n1e = nH * max(1 - xp - 2xh2, 0.0); nee = nH * xp; nHee = nHe0
        if CNn.HECHEM
            z1, z2 = EQ.composition_he(ρ, min(T, 1e6)); nee = nH * (xp + z1 + 2z2); nHee = nH * max(PC.yHe - z1 - z2, 0.0)
        end
        n2e = 4n1e * exp(-E21 / (CNn.kB * T)); nHme = n1e * nee / CNn.Keq([:Hm], [:H, :e], T)
        if G2O.ON && CNn.HECHEM                   # G2OPAC: with the He ions and the LTE n = 2 of G2Opacity
            z1g, z2g = EQ.composition_he(ρ, min(T, 1e6))
            _, Pq, Rq = nlte_means_fast(ρ, T; n1 = n1e, n2 = G2O.n2_lte(T, n1e, nee, n1e + nHee + nH * xh2), nH2 = nH * xh2, nHe = nHee, nHp = nH * xp, ne = nee,
                                        nHm = nHme, hbfP = true, nHep = nH * z1g, nHepp = nH * z2g)
        else
            _, Pq, Rq = nlte_means_fast(ρ, T; n1 = n1e, n2 = n2e, nH2 = nH * xh2, nHe = nHee, nHp = nee, ne = nee, nHm = nHme, hbfP = true)
        end
    end
    bP = max(Pt, 1e-300)^(1 - w) * max(Pq, 1e-300)^w; bR = max(Rt, 1e-300)^(1 - w) * max(Rq, 1e-300)^w
    (bP, bR)
end
"(NLOIT) the ratios (R_E, R_P, R_R) and n_e [cm⁻³] of composition d (full fractions) at (ρ, T), three-level parameters p."
function nlo_ratio(ρ, T, d, p; bPR = nothing)
    t0_ = time_ns()
    nH = ρ / CNn.μH(); nHe0 = PC.yHe * nH
    bP, bR = bPR === nothing ? nlo_base(ρ, T) : bPR
    n1 = nH * max(d[:H], 0.0); ne = nH * max(d[:e], 0.0); np = nH * max(d[:Hp], 0.0)
    nHe = CNn.HECHEM ? nH * max(d[:He], 0.0) : nHe0
    nHmS = n1 * ne / CNn.Keq([:Hm], [:H, :e], T)
    nHm = haskey(d, :Hm) ? nH * max(d[:Hm], 0.0) : nHmS
    Tw = NLOPACJ[] && hasproperty(p, :Tc) ? p.Tc : T
    n2_ = n2_tl3(T, n1, ne, np, p); nH2_ = nH * max(d[:H2], 0.0)
    if MGON && !isempty(MGGIDX[])    # (the groups exist from mg_init! on; a restart evaluates the opacities once before it; Duals through `group_fast` and `mg_pack`: the slopes of the record are exact)
        nHep_ = nH * max(get(d, :Hep, 0.0), 0.0); nHepp_ = nH * max(get(d, :Hepp, 0.0), 0.0)
        y = mg_rec_blend(ρ, T, Tw, n1, ne, np, n2_, nH2_, nHe, nHm, nHep_, nHepp_, p)
        if y !== nothing
            Threads.atomic_add!(NLOITN, 1); Threads.atomic_add!(NLOITNS, Int(time_ns() - t0_))
            return (y[1] / bP, y[2] / bP, y[3] / bR), ne, y[4]
        end
    end
    E, P, R = if G2O.ON
        nHep_ = nH * max(get(d, :Hep, 0.0), 0.0); nHepp_ = nH * max(get(d, :Hepp, 0.0), 0.0)
        hbf_blend(hb -> nlte_means_fast(ρ, T; n1, n2 = n2_, nH2 = nH2_, nHe, nHp = np, ne, nHm, hbfP = hb, Tw, nHep = nHep_, nHepp = nHepp_), p)
    else
        hbf_blend(hb -> nlte_means_fast(ρ, T; n1, n2 = n2_, nH2 = nH2_, nHe, nHp = np, ne, nHm, hbfP = hb, Tw), p)
    end
    Threads.atomic_add!(NLOITN, 1); Threads.atomic_add!(NLOITNS, Int(time_ns() - t0_))
    (E / bP, P / bP, R / bR), ne, MGNOV
end
"""
MGROUPS: (κ_E, κ_P, κ_R of the grey means, group record) of one `mg_group` evaluation at the bound–free weight of p: the means as `hbf_blend`
(geometric in the transition), the record entry by entry (geometric where both ends are positive, else linear in the weight). `nothing` if the
evaluation fails.
"""
function mg_rec_blend(ρ, T, Tw, n1, ne, np, n2, nH2, nHe, nHm, nHep, nHepp, p)
    f(hb) = (y = mg_group(ρ, T; n1, n2, nH2, nHe, nHp = np, ne, nHm, hbfP = hb, Tw, nHep, nHepp); y === nothing ? nothing : (y[5], y[6], y[7], mg_pack(ρ, T, Tw, n1, ne, np, n2, hb, y)))
    w = hbf_w(p)
    hbf_pure(w, 0) && return f(false)
    hbf_pure(w, 1) && return f(true)
    a = f(false); b = f(true); (a === nothing || b === nothing) && return nothing
    gm(x, y) = exp((1 - w) * log(max(x, 1e-300)) + w * log(max(y, 1e-300)))
    (gm(a[1], b[1]), gm(a[2], b[2]), gm(a[3], b[3]), [(a[4][i] > 0 && b[4][i] > 0) ? gm(a[4][i], b[4][i]) : (1 - w) * a[4][i] + w * b[4][i] for i in eachindex(a[4])])
end
"(NLOIT) weight of the table calibration: smooth step in ln T from 0 at NLOITW[1] to 1 at NLOITW[2]."
nlo_w(T) = (x = clamp((log(T) - log(NLOITW[1])) / (log(NLOITW[2]) - log(NLOITW[1])), 0.0, 1.0); x * x * (3 - 2x))
@inline _phvec(ph_::Vector{CNn.Photo}, Tv, env) = _phvec_et(ph_, Tv, env, CNn._ph_et(Tv, env))
@inline _phvec_et(ph_::Vector{CNn.Photo}, Tv, env, ::Type{ET}) where {ET} = Tuple{CNn.Photo,ET,ET}[(p_, CNn.photorates_t(p_, Tv, env)...) for p_ in ph_]
function _phq_combine(a_::Vector{Tuple{CNn.Photo,TA,TA}}, b_::Vector{Tuple{CNn.Photo,TB,TB}}, q::TQ) where {TA,TB,TQ}
    Tuple{CNn.Photo,promote_type(TA,TB,TQ),promote_type(TA,TB,TQ)}[(a_[i][1], b_[i][2] + q * (a_[i][2] - b_[i][2]), b_[i][3] + q * (a_[i][3] - b_[i][3])) for i in eachindex(a_)]
end
"(NLOIT) composition of zone k at (ρ, T) and x_H⁺, x_H₂ (coupled zones: CK expansion with the step's terms; else the chemistry's)."
function nlo_comp(m, k, s, ρ, T, x1, x2, q = 1.0, tc = 1.0, g = nothing, dp = nothing; qs = nothing, xn = nothing)
    t = (TX.cplz(k) && k <= length(CPLTERMS[])) ? CPLTERMS[][k] : nothing
    _nlo_comp(m, k, s, t, ρ, T, x1, x2, q, tc, g, dp, qs, xn)
end
function _nlo_comp(m, k, s, t, ρ, T, x1, x2, q, tc, g, dp, qs, xn)
    t === nothing && return s.d, s.p
    t0_ = time_ns(); xnd = xn === nothing ? nothing : xneut(x1, x2, xn)
    if qs !== nothing && MGON && mg_hasrec(t)         # stage 2c: the composition in the group fields at the iterate's q_g and LW factor (qs = (q_1 … q_G, s)), cpl_tqg
        G_ = TX.MGN[]
        tq = cpl_tqg(t, ρ, T, x1, x2, ntuple(i -> qs[i], G_), qs[G_+1], tc, g, dp; need_cool = false, xn)
        x, _ = CNn.ck_eval((x1, x2), ρ, T, m.net; env = tq.env, phc = tq.phc, extra = tq.extra, xn = xnd)
        Threads.atomic_add!(NLOCNS, Int(time_ns() - t0_))
        return CNn.full_fast(x, m.net; xn = xnd), tq.env.tl3
    end
    # (not t.phc: its one-entry cache is shared with the residual's dual evaluation, which would then take this call's Float rates
    # at the same T and lose ∂(photo-rates)/∂T from the coupled zone's Jacobian)
    ph_ = t.phraw[1]::Vector{CNn.Photo}; en_ = t.phraw[2]
    phc = Tq -> _phvec(ph_, Tq, en_)
    if (q != 1 || tc != 1 || g !== nothing || dp !== nothing) && hasproperty(t, :slq)          # SLIMP: the composition at the iterate's star light (f_sh does not enter the fast species)
        e0 = t.slq.envc0; phc0 = e0 === nothing ? nothing : (Tq -> _phvec(ph_, Tq, e0))
        tq = cpl_tq(t, ρ, T, x1, x2, q, tc, 0.0, g, dp; phc, phc0, need_cool = false, xn)
        x, _ = CNn.ck_eval((x1, x2), ρ, T, m.net; env = tq.env, phc = tq.phc, extra = tq.extra, xn = xnd)
        Threads.atomic_add!(NLOCNS, Int(time_ns() - t0_))
        return CNn.full_fast(x, m.net; xn = xnd), tq.env.tl3
    end
    x, _ = CNn.ck_eval((x1, x2), ρ, T, m.net; env = t.env, phc, extra = t.extra, xn = xnd)
    Threads.atomic_add!(NLOCNS, Int(time_ns() - t0_))
    CNn.full_fast(x, m.net; xn = xnd), (hasproperty(t.env, :tl3) ? t.env.tl3 : s.p)
end
"""
(NLOIT) exact slopes of the ratios (R_E, R_P, R_R) and of n_e at (ρ, T, x_H⁺, x_H₂) and the star-light factors (q, t): ForwardDiff Duals with four
partials (ln ρ, ln T, ln x_H⁺, ln x_H₂) through the same `nlo_comp` (CK expansion) and `nlo_ratio` (means_fast) that give the values. Returns
(d ln R/d ln T, d ln R/d ln x_H⁺, d ln R/d ln x_H₂, d ln R/d ln ρ) as 3-tuples and d ln n_e/d (ln T, ln x_H⁺, ln x_H₂, ln ρ). They replaced
one-sided differences of relative step 1e-2 in x_H⁺ (1e-3 in ρ, central 1e-3 in T): the neutral fraction of the ionised gas, n_H/n = 1 − x_H⁺ ~ 1e-3,
changes by 100 % for 0.1 % in x_H⁺, so a 1 % step crossed x_H⁺ = 1 (max(n_H, 0) = 0) and returned a slope 10× too small (d ln χ/d ln x_H⁺ = −6 against
−67 in zones 15–17 of the stuck state t_nk30).
"""
function nlo_slopes(m, k, s, ρ, T, x1, x2, q, tc, g, dp, qs = nothing, xn = nothing)
    D = ForwardDiff.Dual{TX.JTag,Float64,4}
    sd(v_, i) = D(v_, ForwardDiff.Partials(ntuple(j -> j == i ? v_ : 0.0, Val(4))))      # seed in the logarithm of the variable
    ρd, Td, x1d, x2d = sd(ρ, 1), sd(T, 2), sd(x1, 3), sd(x2, 4)
    d, p = nlo_comp(m, k, s, ρd, Td, x1d, x2d, q, tc, g, dp; qs, xn)
    R, ne, rec = nlo_ratio(ρd, Td, d, p)
    sl(y, i) = ForwardDiff.partials(y, i) / ForwardDiff.value(y)
    if TX.SLIMPDBG[] > 0 && !all(i -> all(j -> isfinite(sl(R[j], i)), 1:3), 1:4)
        TX.SLIMPDBG[] -= 1; println("SLIMPDBG nlo_slopes zone $k: non-finite slope; ρ $ρ T $T x $x1 $x2 q $q tc $tc g $g dp $dp; R ", R)
    end
    gsl(i) = MGON && length(rec) > 0 ? Float64[(v_ = ForwardDiff.value(r_)) > 0 ? ForwardDiff.partials(r_, i) / v_ : 0.0 for r_ in rec] : MGNOV     # (the group record's logarithmic slopes)
    (ntuple(i -> sl(R[i], 2), 3), ntuple(i -> sl(R[i], 3), 3), ntuple(i -> sl(R[i], 4), 3), ntuple(i -> sl(R[i], 1), 3),
     (sl(ne, 2), sl(ne, 3), sl(ne, 4), sl(ne, 1)), (gsl(2), gsl(3), gsl(4), gsl(1)))
end
"""
(NLOIT, SLIMP) exact slopes of the ratios (R_E, R_P, R_R) and of the group record in the star-light variables of the iterate, ForwardDiff Duals through the same `nlo_comp` and `nlo_ratio`
that give the values (one pass; the other variables are Dual constants). Variables, in this order: with the group fields (`qs`: q_1 … q_G, LW factor s; MGROUPS) ln q_1 … ln q_G, ln s,
ln T_c/T_c,start, ln g, then the deep weight (additive, the weight of the bound–free blend); grey: ln q, ln T_c/T_c,start, ln g, deep weight. Returns
(group fields) ((slopes of R_E, R_P, R_R per variable), (record slopes per variable)) or (grey) (s_q, s_t, s_g, s_dp), each the 3-tuple of the slopes of (R_E, R_P, R_R).
They replaced central differences (relative steps 1e-2, 1e-3 in T_c; 1e-4 absolute in the deep weight).
"""
function nlo_slq(m, k, s, ρ, T, x1, x2, q, tc, g, dp, qs, bPR, xn = nothing)
    G = TX.MGN[]; mgs = qs !== nothing; nv = mgs ? G + 4 : 4
    D = ForwardDiff.Dual{TX.JTag,Float64,nv}
    sd(v_, i, lg) = D(v_, ForwardDiff.Partials(ntuple(j -> j == i ? (lg ? v_ : 1.0) : 0.0, Val(nv))))
    cs(v_) = D(v_, ForwardDiff.Partials(ntuple(_ -> 0.0, Val(nv))))
    ρd, Td, x1d, x2d = cs(ρ), cs(T), cs(x1), cs(x2)
    if mgs
        qsd = [sd(qs[i], i, true) for i in 1:G+1]
        d, p = nlo_comp(m, k, s, ρd, Td, x1d, x2d, cs(q), sd(tc, G + 2, true), sd(g, G + 3, true), sd(dp, G + 4, false); qs = qsd, xn)
    else
        d, p = nlo_comp(m, k, s, ρd, Td, x1d, x2d, sd(q, 1, true), sd(tc, 2, true), sd(g, 3, true), sd(dp, 4, false); xn)
    end
    R, _, rec = nlo_ratio(ρd, Td, d, p; bPR = bPR)
    sl(y, i) = (v_ = ForwardDiff.value(y); v_ > 0 ? ForwardDiff.partials(y, i) / v_ : 0.0)
    if TX.SLIMPDBG[] > 0 && !all(i -> all(j -> isfinite(sl(R[j], i)), 1:3), 1:nv)
        TX.SLIMPDBG[] -= 1; println("SLIMPDBG nlo_slq zone $k: non-finite slope; ρ $ρ T $T x $x1 $x2 q $q tc $tc g $g dp $dp; R ", R)
    end
    if mgs
        sE = Float64[sl(R[1], i) for i in 1:nv]; sP = Float64[sl(R[2], i) for i in 1:nv]; sR = Float64[sl(R[3], i) for i in 1:nv]
        gs = [length(rec) == 0 ? zeros(TX.mg_nv()) : Float64[sl(r_, i) for r_ in rec] for i in 1:nv]
        (sE, sP, sR), gs
    else
        ntuple(i -> ntuple(j -> sl(R[j], i), 3), 4)
    end
end
"Per zone n_e [cm⁻³] of the iterate as `nlo_it` last computed it (NaN: not computed in this fill). `fill_opac!` hands it to the formal solution as `xnen`."
const NLONED = Ref(Float64[])
const XNEDIAG = Ref(parse(Int, get(ENV, "XNEDIAG", "0"))); const XNEOLD = Ref(Float64[])
"(NLOIT) (χ, κ_E, κ_P) of network zone k at the iterate (see NLOIT)."
function nlo_it(m::NetGas, k, ρ, T, s, o = nothing)
    κR0, κP0 = OPAC[](ρ, T)
    cz = TX.cplz(k)
    w = TX.CPLW[]; xd = cz ? (w === nothing ? TX.XSN[][1][k] : w[1][k]) : 1.0; xd2 = cz ? (w === nothing ? TX.XSN[][2][k] : w[2][k]) : 1.0
    ρv = Float64(ForwardDiff.value(ρ)); Tv = Float64(ForwardDiff.value(T)); xv = Float64(ForwardDiff.value(xd)); xv2 = Float64(ForwardDiff.value(xd2))
    xnv = cz ? (w === nothing ? TX.XSN[][3][k] : Float64(ForwardDiff.value(w[3][k]))) : 0.0       # the neutral hydrogen of the iterate, direct (no 1 − x_H⁺)
    xnk = cz ? xnv : nothing
    if s === nothing                          # no three-level parameters: the table, χ floored at Thomson on the zone's electrons
        cc = CNn.compiled(m.net); ne = max(CNn.fullvec(m.x[k], cc)[cc.ie], 0.0) * ρv / CNn.μH()
        o === nothing || TX.mg_default!(o)
        return chicap(max(κR0, COp.σT * ne / ρ), T), κP0, κP0
    end
    if TX.OPAC_NOCACHE[]
        d, p = nlo_comp(m, k, s, ρv, Tv, xv, xv2; xn = xnk)
        (RE, RP, RR), ned, _ = nlo_ratio(ρv, Tv, d, p)
        o === nothing || TX.mg_default!(o)
        return chicap(max(κR0 * RR, COp.σT * ned / ρ), T), κP0 * RE, κP0 * RP
    end
    keep = k <= length(NLOC[])                # (the cache is sized by nloc_reset!; a zone beyond it is evaluated without one)
    c = keep ? NLOC[][k] : nothing; v = (ρv, Tv, xv, xv2, xnv)
    if c === nothing || c.v != v
        bPR = (c !== nothing && c.v[1] == ρv && c.v[2] == Tv) ? c.bPR : nlo_base(ρv, Tv)
        d, p = nlo_comp(m, k, s, ρv, Tv, xv, xv2; xn = xnk)
        R, ne, g = nlo_ratio(ρv, Tv, d, p; bPR = bPR)
        c === nothing ? (c = NloC(v, R, ne, bPR, (NaN, NaN, NaN, NaN, NaN), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0), false, g, MGNOV, MGNOV, MGNOV, MGNOV); keep && (NLOC[][k] = c)) :
                       (c.v = v; c.R = R; c.ne = ne; c.bPR = bPR; c.g = g)
    end
    bPR = c.bPR
    needJ = T isa ForwardDiff.Dual || xd isa ForwardDiff.Dual || ρ isa ForwardDiff.Dual
    slon = SLIMP && cz && TX.CPLSL[] !== nothing       # SLIMP: the ratio at the iterate's star-light state (below)
    qv = slon ? Float64(ForwardDiff.value(TX.CPLSL[].q[k])) : 1.0; tv = slon ? Float64(ForwardDiff.value(TX.CPLSL[].t[k])) : 1.0
    gv = slon ? Float64(ForwardDiff.value(TX.CPLSL[].g[k])) : nothing; dv = slon ? Float64(ForwardDiff.value(TX.CPLSL[].dp[k])) : nothing
    # MGROUPS: the iterate's group fields (q_1 … q_G, LW factor) as the star light of the zone's composition
    mgs = slon && MGON && TX.CPLMGQ[] !== nothing && mg_hasrec(CPLTERMS[][k])
    qsd = mgs ? ntuple(i -> TX.CPLMGQ[][i, k], TX.MGN[] + 1) : nothing
    qsv = mgs ? Float64[ForwardDiff.value(x_) for x_ in qsd] : nothing
    if needJ && !(c.hasJ && abs(log(Tv / c.vJ[2])) < NLOITJ[1] && abs(log(ρv / c.vJ[1])) < NLOITJ[1] && (!cz || (abs(log(xv / c.vJ[3])) < NLOITJ[2] && abs(log(xnv / c.vJ[5])) < NLOITJ[2] && abs(log(xv2 / c.vJ[4])) < NLOITJ[2])))
        if G2O.FASTON         # exact slopes (the Duals pass through means_fast)
            c.qT, c.qx, c.qx2, c.qr, c.qn, rs_ = nlo_slopes(m, k, s, ρv, Tv, xv, xv2, qv, tv, gv, dv, qsv, xnk)
            cz || (c.qx = (0.0, 0.0, 0.0); c.qx2 = (0.0, 0.0, 0.0))
            MGON && length(c.g) > 0 && ((c.gT, c.gx, c.gx2, c.gr) = rs_)        # (MGROUPS: the group record's slopes in (T, x_H⁺, x_H₂, ρ), exact)
        else                  # (G2OPAC ≠ 2: the means are Float64 only; differences of R, relative steps 1e-3 and 1e-2)
            εT = 1e-3; d1, p1 = nlo_comp(m, k, s, ρv, Tv * (1 + εT), xv, xv2); R1, _ = nlo_ratio(ρv, Tv * (1 + εT), d1, p1)
            d0, p0 = nlo_comp(m, k, s, ρv, Tv / (1 + εT), xv, xv2); R0, _ = nlo_ratio(ρv, Tv / (1 + εT), d0, p0)
            c.qT = ntuple(i -> log(R1[i] / R0[i]) / (2log1p(εT)), 3)          # (central: the ratio is steep in T, q ~ 10-20)
            if cz
                εx = 1e-2; d2, p2 = nlo_comp(m, k, s, ρv, Tv, xv * (1 + εx), xv2); R2, _ = nlo_ratio(ρv, Tv, d2, p2; bPR = bPR)
                c.qx = ntuple(i -> log(R2[i] / c.R[i]) / log1p(εx), 3)
                d4, p4 = nlo_comp(m, k, s, ρv, Tv, xv, xv2 * (1 + εx)); R4, _ = nlo_ratio(ρv, Tv, d4, p4; bPR = bPR)
                c.qx2 = ntuple(i -> log(R4[i] / c.R[i]) / log1p(εx), 3)
            end
            εr = 1e-3; d3, p3 = nlo_comp(m, k, s, ρv * (1 + εr), Tv, xv, xv2); R3, _ = nlo_ratio(ρv * (1 + εr), Tv, d3, p3)
            c.qr = ntuple(i -> log(R3[i] / c.R[i]) / log1p(εr), 3)
            c.qn = (0.0, cz ? 1.0 : 0.0, 0.0, 0.0)          # (n_e ∝ x_H⁺)
        end
        c.vJ = v; c.hasJ = true; Threads.atomic_add!(NLOITJN, 1)
    end
    RE, RP, RR = c.R
    nebase = c.ne; gb = c.g; gsl = nothing; dualq = false; lsl = nothing
    if slon      # SLIMP: the ratio at the iterate's star-light state (value exact, slopes in ln q, ln t, ln g by differences; those in ρ, T, x_H⁺, x_H₂ above, at q, t, g)
        qd = TX.CPLSL[].q[k]; td = TX.CPLSL[].t[k]; gd = TX.CPLSL[].g[k]; dd = TX.CPLSL[].dp[k]
        length(NLOQC[]) < k && (NLOQC[] = vcat(NLOQC[], Any[nothing for _ in length(NLOQC[])+1:k]))
        e_ = NLOQC[][k]; key = (v..., qv, tv, gv, dv); mgs && (key = (key..., qsv...))
        if e_ === nothing || e_[1] != key
            Rq, neq, gq = nlo_ratio(ρv, Tv, nlo_comp(m, k, s, ρv, Tv, xv, xv2, qv, tv, gv, dv; qs = qsv, xn = xnk)...; bPR = bPR)
            e_ = (key, Rq, neq, e_ === nothing ? nothing : e_[4], e_ === nothing ? ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)) : e_[5], gq, e_ === nothing ? nothing : e_[7])
        end
        dual = qd isa ForwardDiff.Dual || td isa ForwardDiff.Dual || gd isa ForwardDiff.Dual || dd isa ForwardDiff.Dual || (mgs && any(x_ -> x_ isa ForwardDiff.Dual, qsd))
        if dual && !(e_[4] !== nothing && (e_[4] == key || all((i == 8 ? abs(key[i] - e_[4][i]) : abs(log(key[i] / e_[4][i]))) < NLOITJ[i <= 2 ? 1 : 2] for i in 1:8)))
            if mgs
                sl_, gsl_ = nlo_slq(m, k, s, ρv, Tv, xv, xv2, qv, tv, gv, dv, qsv, bPR, xnk)
                e_ = (e_[1], e_[2], e_[3], key, sl_, e_[6], gsl_)
            else
                e_ = (e_[1], e_[2], e_[3], key, nlo_slq(m, k, s, ρv, Tv, xv, xv2, qv, tv, gv, dv, nothing, bPR, xnk), e_[6], e_[7])
            end
        end
        NLOQC[][k] = e_
        RE, RP, RR = e_[2]; nebase = e_[3]; gb = e_[6]
        if dual
            if mgs
                # logarithms of the star-light variables (q_1 … q_G, s, t, g) against their values at the slope point, and the deep weight's difference
                G_ = TX.MGN[]
                lsl = (ntuple(i -> log(qsd[i] / qsv[i]), G_ + 1)..., log(td / tv), log(gd / gv), dd - dv)
                sl_ = e_[5]                                                # (per variable (E, P, R) slopes)
                RE = RE * exp(sum(sl_[1][i] * lsl[i] for i in eachindex(lsl))); RP = RP * exp(sum(sl_[2][i] * lsl[i] for i in eachindex(lsl))); RR = RR * exp(sum(sl_[3][i] * lsl[i] for i in eachindex(lsl)))
                gsl = e_[7]; dualq = true
            else
            lq = log(qd / qv); lτ = log(td / tv); lg = log(gd / gv); ld = dd - dv; sq, st, sg, sdp = e_[5]
            RE = RE * exp(sq[1] * lq + st[1] * lτ + sg[1] * lg + sdp[1] * ld); RP = RP * exp(sq[2] * lq + st[2] * lτ + sg[2] * lg + sdp[2] * ld); RR = RR * exp(sq[3] * lq + st[3] * lτ + sg[3] * lg + sdp[3] * ld)
            end
        end
    end
    lt = zero(T); lx = lt; lx2 = lt; lr = lt
    if needJ
        lt = log(T / Tv); lx = cz ? log(xd / xv) : zero(lt); lr = log(ρ / ρv); lx2 = cz ? log(xd2 / xv2) : zero(lt)
        RE = RE * exp(c.qT[1] * lt + c.qx[1] * lx + c.qr[1] * lr + c.qx2[1] * lx2); RP = RP * exp(c.qT[2] * lt + c.qx[2] * lx + c.qr[2] * lr + c.qx2[2] * lx2)
        RR = RR * exp(c.qT[3] * lt + c.qx[3] * lx + c.qr[3] * lr + c.qx2[3] * lx2)
    end
    ned = needJ ? nebase * exp(c.qn[1] * lt + c.qn[2] * lx + c.qn[3] * lx2 + c.qn[4] * lr) : nebase
    if o !== nothing                     # MGROUPS: the zone's group record in the number type of the iterate
        n = length(gb)
        if n == 0
            TX.mg_default!(o)
        else
            sT = needJ && length(c.gT) == n; sr = needJ && length(c.gr) == n; sx = needJ && cz && length(c.gx) == n && length(c.gx2) == n
            if sT || sr || sx || dualq
                @inbounds for i in 1:n
                    b = gb[i]
                    if b == 0
                        o[i] = zero(eltype(o))
                    else
                        e = zero(eltype(o))
                        sT && (e += c.gT[i] * lt); sr && (e += c.gr[i] * lr); sx && (e += c.gx[i] * lx + c.gx2[i] * lx2)
                        if dualq
                            @inbounds for j in eachindex(lsl); e += gsl[j][i] * lsl[j]; end
                        end
                        o[i] = b * exp(e)
                    end
                end
            else
                copyto!(o, gb)
            end
        end
    end
    k <= length(NLONED[]) && (NLONED[][k] = Float64(ForwardDiff.value(ned)))     # n_e of the iterate, for fill_opac!'s xnen (the scattering share of the formal solution)
    if TX.AD_JV[] && !slon && typeof(ned) <: ForwardDiff.Dual
        neL = get(TX.SLQBUF, (:locop_ne, typeof(ned)), nothing)
        neL !== nothing && k <= length(neL) && (@inbounds (neL::Vector{typeof(ned)})[k] = ned)
    end
    return chicap(max(κR0 * RR, COp.σT * ned / ρ), T), κP0 * RE, κP0 * RP
end
"(diagnostic, NLOITDIAG=n: n sweeps of fill_opac!) per network zone the NLOIT opacities against those the lagged NLOPAC ratio would give at the same state."
const NLOITDIAG = Ref(parse(Int, get(ENV, "NLOITDIAG", "0"))); const NLOITDK = Ref(0); const NLOITSW = Ref(false)
function nloit_diag(m, k, ρ, T, new)
    js = Int(TF.getscalar(:jstep))
    if js != NLOITDK[]
        k == Int(TF.getscalar(:ngrs)) || return nothing          # (one full sweep per step, the first one of the step)
        NLOITDK[] = js; NLOITDIAG[] -= 1; NLOITSW[] = true
    end
    k >= Int(TF.getscalar(:ngre)) && (NLOITSW[] = false)
    NLOITSW[] || return nothing
    (m.eqzone[k] || k > length(NLOST[]) || (s = NLOST[][k]) === nothing) && return nothing
    o = nlopc(ρ, T, s.d, s.p); o === nothing && return nothing
    κR0, κP0 = OPAC[](ρ, T)
    old = (κR0 * o.RR, κP0 * o.RE, κP0 * o.RP)
    @printf("NLOITDIAG step %d zone %d T %.1f rho %.3e x_e %.3e | chi new %.4e old %.4e (table %.4e) | kE new %.4e old %.4e | kP new %.4e old %.4e (table %.4e)\n",
            TF.getscalar(:jstep), k - Int(TF.getscalar(:ngrs)) + 1, T, ρ, s.d[:e], new[1], old[1], κR0, new[2], old[2], new[3], old[3], κP0); flush(stdout)
    nothing
end
# (PERF512) the Float64 opacity of NetGas in convcoef's parallel pre-sweep uses OPAC_NOCACHE so NLOC/NLONED are updated only by the final zone
TX.opac_pure(::NetGas) = true
@inline function TX.material_opac(m::NetGas, ρ, T)
    if NLOIT && NLOPAC[]
        k = TX.ZONE[]
        if k <= length(m.eqzone) && !m.eqzone[k]
            G2B === nothing && return nlo_it(m, k, ρ, T, k <= length(NLOST[]) ? NLOST[][k] : nothing)
            χn, κEn, κPn = nlo_it(m, k, ρ, T, k <= length(NLOST[]) ? NLOST[][k] : nothing)
            w = g2_w(ρ, T); ForwardDiff.value(w) == 0 && return χn, κEn, κPn
            κRt, κPt = OPAC[](ρ, T)                  # G2BLEND: geometric blend with the table opacity at the iterate
            return exp((1 - w) * log(χn) + w * log(chicap(κRt, T))), exp((1 - w) * log(κEn) + w * log(κPt)), exp((1 - w) * log(κPn) + w * log(κPt))
        end
        κR, κP = OPAC[](ρ, T)
        return chicap(κR, T), κP, κP
    end
    κR, κP = OPAC[](ρ, T)
    begin        # (NLOPAC off: before the radiation set-up) a network zone without the non-LTE factors (no star-light rates, so no three-level parameters) took the LTE table, whose
        # equilibrium electrons vanish at ~1000 K: χ 3e-9 cm²/g at 21 pc in c128 (13 Nov) against electron scattering on the zone's own
        # x_e ≈ 9e-5 (3e-5 cm²/g) in its corrected neighbours. Floor χ at σ_T n_e/ρ of the zone's composition.
        k = TX.ZONE[]
        if k <= length(m.eqzone) && !m.eqzone[k]
            cc = CNn.compiled(m.net); ne = max(CNn.fullvec(m.x[k], cc)[cc.ie], 0.0) / CNn.μH()
            κR = max(κR, COp.σT * ne)
        end
    end
    return chicap(κR, T), κP, κP
end
"MGROUPS: as `material_opac`, and the zone's group record into `o` (network zones: at the iterate, `nlo_it`; table zones: the default record)."
@inline function TX.material_opacg!(o, m::NetGas, ρ, T)
    k = TX.ZONE[]
    (NLOIT && NLOPAC[] && k <= length(m.eqzone) && !m.eqzone[k]) && return nlo_it(m, k, ρ, T, k <= length(NLOST[]) ? NLOST[][k] : nothing, o)
    TX.mg_default!(o)
    TX.material_opac(m, ρ, T)
end
"`SLWFLUX=1` (default off): star-light dilution W from max(c(E_r − E_cmb), F_r) instead of c(E_r − E_cmb) alone (see starlight_setup)."
const SLWFLUX = get(ENV, "SLWFLUX", "0") == "1"
function TM.fill_opac!(m::NetGas)
    OPAC[] === nothing && return nothing
    dn, tn = TM.v(:dn), TM.v(:tn)
    chifn, xken, xkpn, xnen = TM.v(:chifn), TM.v(:xken), TM.v(:xkpn), TM.v(:xnen)
    if TX.COUPLED[]      # per-zone n_e of this fill, written by nlo_it (read below only for coupled network zones)
        length(NLONED[]) < length(dn) && (NLONED[] = zeros(length(dn)))
        fill!(NLONED[], NaN)
        XNEDIAG[] > 0 && (length(XNEOLD[]) < length(dn) && (XNEOLD[] = zeros(length(dn))); fill!(XNEOLD[], NaN))
    end
    # (PERF512) zones independent: threaded as fill_eos! (NLOITDIAG prints per zone, so it keeps the serial order)
    sweep! = (NLOIT && NLOITDIAG[] > 0) ? ((f, ks, ok) -> (for k in ks; ok(k) || continue; TX.ZONE[] = k; f(k); end)) : TX.zone_sweep!
    sweep!(TM.krange(), k -> @inbounds(dn[k] > 0 && tn[k] > 0)) do k
      @inbounds begin
        χ, κE, κP = (MGON && k <= size(TX.MGVN[], 2)) ? TX.material_opacg!(view(TX.MGVN[], :, k), m, dn[k], tn[k]) : TX.material_opac(m, dn[k], tn[k])
        chifn[k] = χ; xken[k] = κE; xkpn[k] = κP
        NLOIT && NLOITDIAG[] > 0 && nloit_diag(m, k, dn[k], tn[k], (χ, κE, κP))
        # electrons for the scattering term of the formal solution (the template state's xnen is not ours)
        xnen[k] = (CNn.HECHEM ? (m.eqzone[k] ? (Tq = min(tn[k], 1e6); z = EQ.composition_he(dn[k], Tq); EQ.composition(dn[k], Tq)[1] + z[1] + 2z[2]) :
                                 CNn.full_fast(m.x[k], m.net)[:e]) :
                   (m.eqzone[k] ? EQ.composition(dn[k], min(tn[k], 1e6))[1] : CNn.full_fast(m.x[k], m.net)[:Hp])) * dn[k] / CNn.μH()
        # coupled network zone: the electrons of the iterate (x_H⁺, x_H₂ of the Newton vector, via nlo_it's CK expansion), the same n_e that
        # chifn's Thomson floor uses, not the composition of the last chemistry
        if TX.COUPLED[] && !m.eqzone[k] && TX.cplz(k) && k <= length(NLONED[]) && isfinite(NLONED[][k])
            XNEDIAG[] > 0 && k <= length(XNEOLD[]) && (XNEOLD[][k] = xnen[k])
            xnen[k] = NLONED[][k]
        end
      end
    end
    if XNEDIAG[] > 0 && TX.COUPLED[]            # (XNEDIAG=n, diagnostic: the first n fills print the iterate's n_e against the last chemistry's)
        XNEDIAG[] -= 1; rs = Tuple{Float64,Int}[]
        for k in TM.krange()
            (k <= length(XNEOLD[]) && isfinite(XNEOLD[][k]) && XNEOLD[][k] > 0 && xnen[k] > 0 && !m.eqzone[k] && TX.cplz(k)) && push!(rs, (log(xnen[k] / XNEOLD[][k]), k))
        end
        if !isempty(rs)
            lo, hi = extrema(rs); a_ = Int(TF.getscalar(:ngrs))
            @printf("XNEDIAG jstep %d: %d coupled zones, n_e(iterate)/n_e(last chemistry): min %.4g (zone %d) max %.4g (zone %d); zones with |ln| > 0.1: %d, > 1: %d\n", TF.getscalar(:jstep),
                    length(rs), exp(lo[1]), lo[2] - a_ + 1, exp(hi[1]), hi[2] - a_ + 1, count(r -> abs(r[1]) > 0.1, rs), count(r -> abs(r[1]) > 1, rs))
        end
    end
    nothing
end

"Free-species vector from the equilibrium table's composition (H₂, H⁺; H⁻ and Li⁺ in Saha with the table's electrons)."
function x_from_table(m::NetGas, ρ, T)
    y = EQ.composition(ρ, T); nH = ρ / CNn.μH()
    xH2 = y[3]; xHp = y[1]; xH = 1 - xHp - 2xH2
    ne = xHp * nH
    z1, z2 = CNn.HECHEM ? EQ.composition_he(ρ, T) : (0.0, 0.0)
    CNn.HECHEM && (ne = (xHp + z1 + 2z2) * nH)
    xHm = xH * ne / CNn.Keq([:Hm], [:H, :e], T)
    out = zeros(length(m.net.free))
    for (i, s) in enumerate(m.net.free)
        out[i] = s === :Hp ? xHp : s === :H2 ? xH2 : s === :Hm ? xHm : s === :Hep ? z1 : s === :Hepp ? z2 :
                 s === :Lip ? (K = CNn.Keq([:Li], [:Lip, :e], T); CNn.yLi * K / (K + max(ne, 1e-300))) : 0.0
    end
    out
end

"""
S6: hysteresis of the EOS switch: a zone in the table regime stays there until p < PSWHYS × pswitch (1 = the sharp switch of G5). Warm gas
of the accretion layer sitting on pswitch otherwise flips every step between table and network, which differ by ~10 % in T at fixed e.
"""
const PSWHYS = Ref(1.0)
const PSWTMIN = Ref(0.0)
"""
Table to network, falling pressure: PSWTEX[] > 0 keeps a table zone on the table until T < PSWTEX[], whatever its pressure. The network has no photoionization by the local radiation field, so
at 8500 K it holds x(H+) well below Saha at fixed e, and the recombination energy went into heat: T 8509 -> 13938 K in one step at the surface of a 28 M_sun protostar, then Newton failed at every dt.
E_ion/E_th at rho = 1e-11..1e-10 is <= 4 % at 5000 K, 15-48 % at 6000 K. Zones the ray reaches go to the network before this test.
"""
const PSWTEX = Ref(0.0)
"""
PSWSTAY[] = false drops the PSWTEX "stay" rule above, so a table zone below the pressure switch stays on the table only by the τ rule (PSWTAU).
With the three-level network (SL3) the network carries the star-light ionization the rule was written to stand in for; the rule otherwise holds table
zones out to τ ≈ 0.4 while T ≥ PSWTEX, so the boundary becomes a temperature switch in the zones where table and network disagree most.
"""
const PSWSTAY = Ref(true)
"""
PSWTAU[] > 0: a zone with T >= PSWTEX[] and Rosseland depth from the surface >= PSWTAU[] (half that to stay) uses the table whatever its pressure. Where radiation is thermalised, J = B holds H at Saha.
The network has no radiation field: at the outer edge of an opaque layer (p ~ 6, T ~ 1e4 K, 39 M_sun) it sits at coronal x(H+) ~ 0.05 against 0.99, Newton heats the zone by a factor of 2 at frozen
composition and chem! brings it back each step, until Newton fails. On the table the ionization heat capacity is inside Newton. The switch at fixed e lowers T once; radiation reheats the zone.
"""
const PSWTAU = Ref(0.0)
"""
PSWRAMP=n (default 0 = off): a gradual hand-over from the table to the network. From the step the run (re)starts, a zone still on the table
stays there, except the outermost (j − j0) ÷ n of them, which follow the switch rules (counted from the outermost table zone below the hot branch at the start; with the table rules off, PSWITCH=1e30
PSWTAU=0 PSWSTAY=0: they leave the table, one zone every n steps from the outside in). Network zones are never put back on the table by it.
"""
const PSWRAMP = parse(Int, get(ENV, "PSWRAMP", "0")); const RAMPJ0 = Ref(-1); const RAMPK0 = Ref(0)
"""
PSWNH=n (default 0 = off; cm⁻³): the table also where n_H ≥ n and T ≥ PSWTMIN (hysteresis as the pressure switch: stay while n_H ≥ PSWHYS n,
T ≥ 0.95 PSWTMIN), for the non-ideal regime. The network's EOS is ideal (H partition sum to n = 15, no occupation probabilities, no
pressure ionization); against the table's SCVH it departs in partly ionized gas from n_H ~ 1e17 (x_e 0.82 vs 0.71, u −59 % at 1e18, 2e4 K)
and everywhere at n_H ≳ 1e20 (10 M☉ centre: x_e 0.85 vs 0.04, p −39 %). With PSWITCH=1e30 PSWTAU=0 PSWSTAY=0 this is the minimal table.
"""
const PSWNH = parse(Float64, get(ENV, "PSWNH", "0"))
"SLBLEND[]: star-light colour temperature and dilution blended smoothly in τ (see starlight_setup); false: switch at the τ = 2/3 zone."
const SLBLEND = Ref(false)
"""
SLEB[] (with SLBLEND): the star-light colour temperature of zone k is local, T_r at τ_k + 2/3 (Eddington–Barbier: the outward light at depth τ_k
comes from about Δτ = 2/3 further in), interpolated in log τ, in place of the photosphere's T_eff. In the extended precursor τ is nearly flat
(0.8 → 0.4 over 400–1400 R_sun at 36 M_sun), the τ = 2/3 zone hopped between 122 and 130 from step to step, and T_eff = (L/4πr²σ)^¼ there
swung 5580 ↔ 7100 K at constant L: every lit zone's spectrum jumped by 1500 K, faster than its ionization relaxes (C α_B n_e ~ 1e-5 s⁻¹).
"""
const SLEB = Ref(false)
"""
NLOPACJ[] (with NLOPAC): κ_E of a network zone, the absorption mean of the gray E_r, weighted with the star light's spectrum B_ν(T_c) (T_c of
the three-level parameters) in place of B_ν(T_gas). In the precursor T_gas ≈ 4800 K under T_c ≈ 6000 K light; H⁻ bound–free stops at 1.64 μm,
where B(4800 K) has more of its energy than B(6000 K).
"""
const NLOPACJ = Ref(false)
const NSWIN = Ref(0); const NSWOUT = Ref(0)
"T at fixed (ρ, u) for network composition x (a vector, or a function of T), bisection in ln T, NaN without a bracket."
function T_at_u(ρ, u, x, net, T0)
    F(T) = CNn.eos(ρ, T, x isa Function ? x(T) : x, net)[1] - u
    lo, hi = T0 / 1.5, T0 * 1.5
    while F(lo) > 0 && lo > 20; lo /= 1.5; end
    while F(hi) < 0 && hi < 1e7; hi *= 1.5; end
    (F(lo) > 0 || F(hi) < 0) && return NaN
    for _ in 1:100; mid = sqrt(lo * hi); F(mid) > 0 ? (hi = mid) : (lo = mid); hi / lo - 1 < 1e-9 && break; end
    sqrt(lo * hi)
end
"After a converged hydro step: network step per zone in the ideal regime, equilibrium-table step where p ≥ pswitch."
const CHEMFALLBACK = Ref(false); const NCHEMFB = Ref(0)
"""
CHEMFBDENSE=1 (default): a failed network step of a zone whose pressure is above the table switch (p ≥ pswitch; only its T < PSWTMIN keeps
it on the network) takes the table-equilibrium fallback too, in every phase, counted in NCHEMFB. Such gas is in chemical equilibrium to a
good approximation (n ≳ 1e15 cm⁻³), and it is where the backward-Euler Newton has failed in the collapse runs: the H₂ dissociation layer
just below PSWTMIN at Δt of hundreds of seconds (full3, core growth stage, step 11935: n_H 4.5e18, T 4346 K, Δt 515 s, twelve halvings
without convergence), which without the fallback ends the run. 0: the fallback only where CHEMFALLBACK is set (the accretion phase).
"""
const CHEMFBDENSE = Ref(get(ENV, "CHEMFBDENSE", "1") == "1")
"""
M3: the star's ionizing light (IonRay), `nothing` = off. `RAY[] = (; start, source, env, Ec, Q, xmin)`: `start(a, b)` -> (photosphere zone, photons per group) or `nothing`; `source` as IonRay.source;
`env`: ChemNet environment for lit zones (case B); `Q`: energy taken per zone, shared with TitanResidFast.RAYQ. Zones are done outward from the photosphere, each with the photons the last one left. A zone
where the photons of one step could ionize more than `xmin` of its hydrogen stays on the network whatever its pressure, because the table assumes equilibrium.
"""
const RAY = Ref{Any}(nothing)
"""
M3a: the star's light below the Lyman edge in the network zones (`STARLIGHT[] = (; pHm, pH, shield)`, `nothing` = off). Each zone sees
J_ν = W B_ν(T_col): outside the photosphere (τ_R = 2/3, from the outside) T_col = T_eff of the star and W = c(E_r − aT_CMB⁴)/(4σT_eff⁴), so the
dilution and attenuation are those of TITAN's own transfer; inside it T_col = (E_r/a)^¼ and W = 1. Rates, lagged over the step (start-of-step T):
- H⁻ + γ → H + e (`pHm`, John 1988), and optionally (`pH`, `nothing` = off, the default) H(n ≥ 2) + γ → H⁺ + e (the network's H levels without
  n = 1: Lyman-continuum photons are IonRay's, photon-conserving; a dilute Wien tail applied locally would ionize every zone). No extra heat for
  these: in RAD runs the network has no radiative cooling either, all radiative energy exchange is the gray κ_P term (LTE opacities), and ChemNet
  pays the binding energy out of u. That is harmless for H⁻ (x ~ 1e-8, 0.75 eV) but not for H(n ≥ 2): near the photosphere k_H ~ 1 s⁻¹ ionizes
  in ≪ dt what u cannot pay for (13.6 eV per H against ~1 eV of thermal energy), T collapses and ChemNet's Newton fails (zone at 7700 K, 2026-11-05).
  Doing it right needs the H bound-free heating and recombination losses in the network and the H bound-free part removed from κ_P there.
- Lyman-Werner H₂ dissociation (two-step Solomon), k = 1.1e8 F_ν(12.87 eV) s⁻¹ (AAZN97 k27, F_ν = 4πJ_ν), times the H₂ self-shielding
  f_sh(N_H2, T) of Wolcott-Green, Haiman & Bryan (2011, eq. 10, α = 1.1) and the shielding by H I Lyman lines, their eq. 13 (β = 1.6, γ = 0.15,
  N_HI/2.85e23 cm⁻²); fits valid to a factor 2 for 1e13 < N_H2 < 1e21, 1e22 < N_HI < 1e24, 500 < T < 5000 K, static gas. Columns from the
  photosphere out to the zone centre (start-of-step composition). `shield = :static` (that), `:sobolev` (the H₂ column of one Sobolev length
  n v_th/|du/dr| when smaller: the infall's Doppler shifts move the H₂ lines out of the shielding; N_HI stays static, since the H I Lyman lines
  shield through damping wings thousands of km/s wide at N_HI ~ 1e24), `:none` (no shielding), `:off` (no LW). Heat
  4.48 + 0.4 eV per dissociation (binding + fragment kinetic energy; UV-pumping heat in dense gas is left out), moved from E_r to the gas in the
  hydro step (TitanResidFast.GASQ), not in the chemistry: there, without the gas's radiative cooling, the H⁻ → H₂ → LW cycle heated a zone by
  ×2.5 in one step.
`SLDIAG[]`: per-step diagnostics (photosphere zone, T_eff, and per zone k_H⁻, k_H, k_LW, f_sh).
"""
const STARLIGHT = Ref{Any}(nothing); const SLDIAG = Ref{Any}(nothing)
const EVERG = 1.602176634e-12
"a/expm1(x); for a Dual x as a e^{-x}/(1 − e^{-x}): the derivative of 1/expm1 divides by expm1², which overflows for x > 355 (the Wien tail of a cold T_c)."
@inline _div_expm1(a, x) = a / expm1(x)
@inline _div_expm1(a, x::ForwardDiff.Dual) = a * exp(-x) / (-expm1(-x))
_bnu(ν, T) = (x = CNn.hP * ν / (CNn.kB * T); x > 700 ? 0.0 : _div_expm1(2CNn.hP * ν^3 / CNn.c^2, x))
"H₂ self-shielding, Wolcott-Green, Haiman & Bryan (2011) eq. (10) with α = 1.1; b = (2kT/m_H₂)^½."
fsh_h2(N, T) = (x = N / 5e14; b5 = sqrt(CNn.kB * T / 1.6735575e-24) / 1e5; 0.965 / (1 + x / b5)^1.1 + 0.035 / sqrt(1 + x) * exp(-8.5e-4 * sqrt(1 + x)))
"Shielding of H₂ by H I Lyman-series lines, Wolcott-Green & Haiman (2011) eq. (13)."
fsh_hi(N) = (x = N / 2.85e23; (1 + x)^-1.6 * exp(-0.15x))
"""
Three-level hydrogen (1s, n = 2, continuum) for the network zones under STARLIGHT (`STARLIGHT[].tl3`): the Peebles (1968) C factor as coded in
ChemistryKernels.jl (recombination.jl, `_peebles_k2_from_beta`, RECFAST constants), with the Hubble rate replaced by the Sobolev velocity
gradient of the flow and electron collisions added. Per 2s atom (2s and 2p in statistical equilibrium, n_2p = 3 n_2s):
  C = (R_Lyα + Λ_2s + 4 n_e q21) / (R_Lyα + Λ_2s + 4 n_e q21 + β),
R_Lyα = 3A_21 (1 − e^{−τ_S})/τ_S, τ_S = 3A_21 λ_α³ n_1s / (8π (dv/dr)) (→ 8π (dv/dr)/(λ_α³ n_1s), CK's 1/(K n_1s) with H → dv/dr), dv/dr the
angle mean (|∂u/∂r| + 2|u|/r)/3; Λ_2s = 8.2245809 s⁻¹; β = W α_B(T_c) (C_R T_c)^{3/2} e^{−χ_2/kT_c} (RECFAST's Milne β of a blackbody, diluted by W);
α_B the RECFAST fit (CK `recfast_alpha`); q12 = Λ_ceHI/E_21 from Black (1981, via Cen 1992, as EmissionKernels ceHI: 7.5e-19 e^{−118348/T}
/(1 + (T/1e5)^½) erg cm³ s⁻¹), q21 = q12 e^{E_21/kT}/4. Net rate per volume
  dn_p/dt = (1 − C) n_1s (n_e q12 + P12) − C α_B n_e n_p,   P12 = (R_Lyα + Λ_2s) n̄_α,  n̄_α = W/(e^{E_21/kT_c} − 1)
(pumping by the local field at Lyα). These replace the Boltzmann (LTE) excited levels of the H photo pair and, on the spot, the 1s recombination
(Lyman-continuum photons: IonRay). With T_c = T, W = 1 the pair gives Saha exactly (ground g = 2, the network's Z_H = 2 to 1e-6 below 1e4 K);
with collisions dominant (n_e q21 ≫ R_Lyα, Λ, β) n_2 is Boltzmann at T and the old LTE treatment is recovered.
Energy: ChemNet pays χ_H out of u per net ionization (E0). The chemistry also gets, per ionization through n = 2, the photon energy ⟨hν⟩ (Balmer
continuum mean, σ ∝ ν⁻³ in W B(T_c)) + E_21 P12/(n_e q12 + P12) (the pumped share of the excitation), and per net recombination loses
χ_H + 0.75 kT to radiation; that energy leaves (or enters) the gray field through RAYQ. Inside the chemistry so that the collisional-excitation
cost (10.2 eV per ionization, ≫ kT) limits the rate implicitly in T. The recombine/re-ionize cycle's heat (1 − C) α_B n_e n_p (⟨hν⟩ − χ_2 −
0.75 kT) goes to the gas in the hydro step (GASQ), where it meets the gas's radiative cooling. Left out: H–H collisional excitation (Drawin),
n ≥ 3 except through α_B, and the Lyα cooling of excitations that do not ionize (in the gray κ_P's lines).
"""
const TL3DBG = Ref(parse(Int, get(ENV, "TL3DBG", "0")))
const LYA_LAM = 1.215668e-5; const LAM2S = 8.2245809; const A21LYA = 6.2648e8
const E21 = 10.19881 * EVERG; const CHI2 = CNn.χH - E21; const CR_REC = 1.799920e10      # 2π m_e k/h² [cm⁻² K⁻¹]
alphaB_recfast(T) = (t = T / 1e4; 4.309e-13 * t^-0.6166 / (1 + 0.6703 * t^0.53))          # cm³ s⁻¹ (= CK recfast_alpha × 1e6)
beta2s_bb(Tc) = alphaB_recfast(Tc) * (CR_REC * Tc)^1.5 * exp(-CHI2 / (CNn.kB * Tc))         # s⁻¹ per 2s atom, J = B(T_c)
q12_H(T) = 7.5e-19 * exp(-118348.0 / T) / (1 + sqrt(T / 1e5)) / E21
q21_H(T) = 7.5e-19 / (1 + sqrt(T / 1e5)) / E21 / 4
tl3_C(T, p) = (D = p.RL + LAM2S + 4p.ne * q21_H(T); D / (D + p.β))
"Mean photon energy of n = 2 photoionization (σ ∝ ν⁻³) in J ∝ B(T_c)."
function balmer_mean(Tc)
    ν0 = CHI2 / CNn.hP; s0 = 0.0; s1 = 0.0
    xe = log(1 + 40CNn.kB * Tc / CHI2); xs = Tc isa ForwardDiff.Dual ? [xe * (i / 100) for i in 0:100] : range(0.0, xe, length = 101); νa = ν0 * exp(xs[1])
    for i in 1:100
        νb = ν0 * exp(xs[i+1]); ν = 0.5(νa + νb); w = _bnu(ν, Tc) / ν^4 * (νb - νa); s0 += w; s1 += w * CNn.hP * ν; νa = νb
    end
    s0 > 0 ? G2O._qd(s1, s0) : CHI2
end
"The three reactions of the three-level H atom (zero unless the zone's env carries `tl3`)."
function tl3_reactions()
    on(env) = hasproperty(env, :tl3) && env.tl3 !== nothing
    [CNn.Reaction("H + e → H(n=2) + e, then n = 2 photoionization (three-level, (1−C) q12)", [:H, :e], [:Hp, :e, :e],
                  (T, env) -> on(env) ? (1 - tl3_C(T, env.tl3)) * q12_H(T) : zero(T)),
     CNn.Reaction("H + γ(Lyα) → H(n=2), then n = 2 photoionization (three-level, (1−C) P12)", [:H], [:Hp, :e],
                  (T, env) -> on(env) ? (1 - tl3_C(T, env.tl3)) * env.tl3.P12 + zero(T) : zero(T)),
     CNn.Reaction("H+ + e → H(n≥2) → H(1s) (three-level, C α_B)", [:Hp, :e], [:H],
                  (T, env) -> on(env) ? tl3_C(T, env.tl3) * alphaB_recfast(T) : zero(T))]
end
"""
The network with the three-level H atom: `tl3_reactions` added, and mutual neutralization H⁻ + H⁺ → 2H (GA08 k5) split by the same C: it
leaves one H in n ≥ 2 (mostly n = 3), which the star light ionizes again with probability 1 − C, so in lit zones only C of it recombines
(the rest is H⁻ + H⁺ → H + H⁺ + e). Without the split, at x_e ~ 0.4 it recombined ~1 % of the H per step and put 12.9 eV per event into
the gas (the network has no radiative losses in RAD runs). Its detailed-balance reverse is unchanged.
"""
function tl3_network(net)
    on(env) = hasproperty(env, :tl3) && env.tl3 !== nothing
    rx = CNn.Reaction[]
    for r in net.reactions
        if startswith(r.name, "H- + H+ → 2H") && !occursin("reverse", r.name)
            k0 = r.k
            push!(rx, CNn.Reaction(r.name * " × C (three-level)", r.R, r.P, (T, env) -> on(env) ? tl3_C(T, env.tl3) * k0(T, env) : k0(T, env)))
            push!(rx, CNn.Reaction("H- + H+ → H(n≥2) + H+, re-ionized: H + H+ + e ((1−C) × GA08 k5)", [:Hm, :Hp], [:H, :Hp, :e],
                                   (T, env) -> on(env) ? (1 - tl3_C(T, env.tl3)) * k0(T, env) : zero(T)))
        else
            push!(rx, r)
        end
    end
    CNn.Network(net.species, net.free, MGON ? vcat(rx, tl3_reactions(), mg_h1s_reactions()) : vcat(rx, tl3_reactions()))
end
_fv(x::ForwardDiff.Dual) = _fv(ForwardDiff.value(x)); _fv(x::Real) = Float64(x)
"""
LOCALJ=1 (default off): in the zones with three-level parameters, the photo-processes of the network (He I, He II, H⁻, Li, from the levels with
their Milne inverses) see the diluted blackbody of the star-light rates, J = W B(T_c) (T_c from T_eff in thin gas to the local T_r at τ ≥ 3,
W = 1 there), instead of the static field of `m.env` (the CMB). Without it He in the thermalised zones that leave the table is ionized only
collisionally against radiative recombination (coronal balance): x_e 1.10 instead of 1.167 (full He²⁺) at 7e4 K, 1e-11 g/cm³ (t_nt4).
W is recovered from the parameters' β = W β_2s(T_c). The rates are then computed directly (no table in T).
"""
const LOCALJ = get(ENV, "LOCALJ", "0") == "1"
# LOCALJTAB (with LOCALJ=1, default 1): the rates from ChemNet.photorates_lj's tables (W a factor, T_c and T tabulated); 0: direct integrals
const LOCALJTAB = get(ENV, "LOCALJTAB", "1") == "1"
# The field is the static one (CMB) plus the diluted blackbody, J = J_static + W B(T_c) (`ljs` = the static env: its rates are added to the table's, ChemNet._photorates_t).
# The former LOCALJWMIN switch (static env alone where W < w, the far thin cloud of t_tau3) is gone: there W B(T_c) ≪ J_static and the sum is the static env.
# MGROUPS: the photo list's rates come from the static field plus the group fields (`mgq`, ChemNet.photorates_mgq).
ljon(p) = LOCALJ && p !== nothing && hasproperty(p, :Tc) && p.Tc > 0
tl3env(m, p) = (MGON && p !== nothing && hasproperty(p, :mgq)) ? merge(m.env, (; tabkey = m.env, mgq = p.mgq, tl3 = p)) : ljon(p) ?
    (W_ = p.β / beta2s_bb(p.Tc); Tc_ = p.Tc; Js = m.env.J;
     LOCALJTAB ? merge(m.env, (; Trad = Tc_, J = ν -> W_ * _bnu(ν, Tc_), static = false, ljW = W_, ljs = m.env, tl3 = p)) :
                 merge(m.env, (; Trad = Tc_, J = ν -> Js(ν) + W_ * _bnu(ν, Tc_), static = false, ljs = m.env, tl3 = p))) :
    merge(m.env, (; tabkey = m.env, tl3 = p))
"""
CHEMFBNET=n (default 0 = off): the chemistry fallback of a network zone below the equilibrium table's density (PopIIIEqEOS.rho_min, where
EQ.step would throw: the far cloud, n_H ~ 0.3, t_tau3 at 455 M☉) stays on the network: the step again in n equal sub-steps (u carried from
sub-step to sub-step through the network EOS); if that fails too, the composition is held and T solved from u at fixed composition.
"""
const CHEMFBNET = parse(Int, get(ENV, "CHEMFBNET", "0")); const NCHEMFBNET = Ref(0)
function chem_fallback_net(m, k, ρ, u, T0, dt, env, ph, cools, exs)
    x = copy(m.x[k]); T = T0; ua = u
    try
        h = dt / CHEMFBNET
        for _ in 1:CHEMFBNET
            x, T = CNn.step(x, ρ, ua, T, h, m.net; env, photos = ph, cooling = cools, extra = exs)
            (all(isfinite, x) && isfinite(T) && T > 0) || error("non-finite")
            ua = CNn.eos(ρ, T, x, m.net)[1]
        end
        return (x, T)
    catch
    end
    x = m.x[k]; f(t) = CNn.eos(ρ, t, x, m.net)[1] - u
    lo, hi = 1.0, 1e7
    (f(lo) <= 0 <= f(hi)) || return (x, T0)
    for _ in 1:200
        mid = sqrt(lo * hi); f(mid) > 0 ? (hi = mid) : (lo = mid); hi / lo < 1 + 1e-12 && break
    end
    (x, sqrt(lo * hi))
end
"""
TL3DEEP=1 (default off): no three-level heats (chemistry heat, cycle heat, LYABB bound–bound) in zones at τ ≥ 3, where the star-light blend
is thermalised (T_c = T_r, W = 1). There the gas–radiation exchange is the implicit κ_P term; the three-level terms do not cancel in LTE (the
cycle heat leaves (⟨hν⟩(T_c) − χ_2 − 0.75 kT) ≈ a few % of kT per cycle) and, multiplied by α_B n_e n_p of hot ionized gas, cooled the deep He II
zone that left the table by ~2000 K per step (zone 73, 2.6e4 K, cycle heat −2.65e-3 erg cm⁻³ s⁻¹ lagged over 1.3e8 s against ρe ≈ 6 erg cm⁻³).
"""
const TL3DEEP = get(ENV, "TL3DEEP", "0") == "1"
"Weight of the thermalised-zone treatment (TL3DEEP) at optical depth τ: smooth step in ln τ from 0 at τ = 2 to 1 at τ = 4.5 (centre 3, the former discrete flag τ ≥ 3)."
sl_deepw(τ) = (x = clamp(log(max(τ, 1e-30) / 2) / log(2.25), 0, 1); x * x * (3 - 2x))
"Energy terms of the three-level atom per volume: (chemistry heat, cycle heat for GASQ) at (ρ, T, full composition d); × (1 − `p.deep`) with TL3DEEP."
function tl3_heat(ρ, T, d, p)
    wd = TL3DEEP && hasproperty(p, :deep) ? p.deep : 0.0
    wd >= 1 && return (zero(T), zero(T))
    nH = ρ / CNn.μH(); C = tl3_C(T, p); xH = max(d[:H], zero(T)); xe = max(d[:e], zero(T)); xp = max(d[:Hp], zero(T))
    Q12 = nH * xe * q12_H(T); I = (1 - C) * xH * (Q12 + p.P12); rec = alphaB_recfast(T) * nH * xe * xp
    if iszero(p.P12)
        I_Ein = I * p.hν
    elseif _val(Q12 + p.P12) < 1e-150
        I_Ein = (1 - C) * xH * ((Q12 + p.P12) * p.hν + E21 * p.P12)
    else
        Ein = p.hν + E21 * p.P12 / (Q12 + p.P12 + 1e-300)
        I_Ein = I * Ein
    end
    chemq = nH * (I_Ein - C * rec * (CNn.χH + 0.75CNn.kB * T))
    cycq = nH * (1 - C) * rec * (p.hν - CHI2 - 0.75CNn.kB * T)
    LYABB[] && (chemq += lya_bb(ρ, T, d, p))
    (1 - wd) * chemq, (1 - wd) * cycq
end
"Bound–bound n = 1 ↔ 2 heat of the three-level atom per volume [erg cm⁻³ s⁻¹] (see LYABB)."
function lya_bb(ρ, T, d, p)
    nH = ρ / CNn.μH(); xH = max(d[:H], zero(T)); xe = max(d[:e], zero(T)); xp = max(d[:Hp], zero(T))
    # branching of n = 2 as in the reactions and in chemq: C and D with the lagged p.ne, so the de-excited share 4 p.ne q21/(D + β) ≤ 1 of
    # every n = 2 atom is the one the network's C assumes (with the current n_e in the collisional term and p.ne in D, a zone ionizing
    # within the step de-excited 41× its n = 2 atoms: energy from nothing, 3654 K → 5e5 K, l2 step 99005); collisional excitation with
    # the current n_e (as the network's H + e → H(n=2) rate), its radiative share (R + Λ)/(D + β) ≤ 1
    rec = alphaB_recfast(T) * nH * xe * xp; D = p.RL + LAM2S + 4p.ne * q21_H(T); C = D / (D + p.β)
    nH * E21 * C * ((xH * p.P12 + rec) * 4p.ne * q21_H(T) - xH * nH * xe * q12_H(T) * (p.RL + LAM2S)) / D
end
"""
LYABB=1: the bound–bound energy exchange of the three-level atom with the gas, missing from tl3_heat's ionization bookkeeping. The atom's
exact radiative balance with the gas is E21 [n1 P12 − n2s (R_Lyα + Λ_2s)] + β n2s ⟨hν⟩ − α_B n_e n_p (χ_2 + 0.75 kT) (n2s in steady state);
tl3_heat's chemq + cycq differ from it by
  E21 C [(n1 P12 + α_B n_e n_p) 4 n_e q21 − n1 n_e q12 (R_Lyα + Λ_2s)] / D,   D = R_Lyα + Λ_2s + 4 n_e q21:
n = 2 atoms made by Lyα/two-photon pumping (excitation temperature of the diluted starlight T_c) or by recombination that are collisionally
de-excited heat the gas; collisional excitations that decay radiatively (Lyα escape, two-photon continuum) cool it. With W = 1, T_c = T the
sum vanishes (detailed balance; without the term, LTE gas lost E21 β n2s 4 n_e q21/(D+β)). The energy goes to/from the gray field via RAYQ
with the rest of chemq, implicit in T in the chemistry. Left out: H–H (Drawin) collisions.
"""
const LYABB = Ref(get(ENV, "LYABB", "0") == "1")
"""
`SLIMP>0` (with COUPLED, SLBLEND=1; default off, off = bit-identical): the star light of a coupled zone at the Newton iterate. Every
star-light quantity of the zone's chemistry follows the iterate: the dilution W (the zone's own E_r and its two interface fluxes, SLWFLUX:
c u = max(c(E_r − E_cmb), ½(F_k + F_k+1)), W = c u/(4σT_c⁴) capped at 1), the colour temperature T_c = T_src^(1−wq) T_r^wq, the Lyman-Werner
factor and Lyα escape. Every star-light rate is linear in W at fixed T_c: the photo-rates of J = W B_ν(T_c) (R ∝ W; the Milne inverse
α = α_spont + W α_stim), kHm, kH, kLW (∝ W), the three-level β = W β_2s(T_c) and P12 = (R_Lyα + Λ_2s) W/(e^{E21/kT_c} − 1). So with the state
S = (q, t, fs, g) of `sl_sweep!` (q = W(iterate)/W(start of step), t = T_c(iterate)/T_c(start), fs, g) the prepared terms are rebuilt exactly
(`cpl_tq`), and S is a local variable of the coupled zone (`cpl_zone_q`, local 7×8 Jacobian), of the NLOIT opacity composition (`nlo_comp`,
slopes in ln q, ln t, ln g) and of the post-step species (`cpl_finish!`).
"""
const SLIMP = parse(Int, get(ENV, "SLIMP", "0")) > 0
"(SLIMP) per zone the star-light record of the prepare pass (W0, T_c0, wq, T_src, mode, E_cmb, ρ, T, fsh, LW on); nothing: no star light."
const SLW0 = Ref(Any[])
"(SLIMP) a star exists at the last prepare pass (a zone at τ ≥ 2/3 and T_eff > 0)."
const SLSTAR = Ref(false)
const SDPH = get(ENV, "SDPH", "0") == "1"; const SDPHJ = Ref(0); a_sd() = Int(TF.getscalar(:ngrs))   # (stalldiag) photosphere index log
"The angle-mean velocity gradient of the Lyα escape (|∂u/∂r| + 2|u|/r)/3 of the zone between r1 and r2 with edge velocities u1, u2."
sl_grad(r1, r2, u1, u2) = (abs(u2 - u1) / (r2 - r1) + 2abs(0.5(u1 + u2)) / (0.5(r1 + r2))) / 3
"""
Window of the band Jacobian: the colour-coded Duals of the residual carry the derivatives of blocks k-2..k+2 only (TitanResidFast.assemble!;
blocks 3 or more apart share a seed slot). A non-local quantity of zone k is therefore built from the Duals of zones j with |j − k| ≤ 2 and the
values of all others; the rest of its dependence on the iterate is seen by the Jacobian-vector products of the Newton–Krylov solve.
"""
@inline winv(x::ForwardDiff.Dual{<:Any,<:Any,1}, j, k) = x
@inline winv(x::ForwardDiff.Dual, j, k) = abs(j - k) <= 2 ? x : ForwardDiff.value(x)
@inline winv(x, j, k) = x
@inline _val(x) = x isa ForwardDiff.Dual ? ForwardDiff.value(x) : x
"On-weight of the Lyman-Werner shielding column of a zone at optical depth τ: 1 outside, 0 inside the photosphere, smooth in ln τ between τ = 0.5 and 0.9 (the smooth form of the zone index k > k_ph, τ_k < 2/3)."
sl_lwon(τ) = (x = clamp(log(max(τ, 1e-30) / 0.5) / log(0.9 / 0.5), 0, 1); 1 - x * x * (3 - 2x))
"""
(SLIMP) the star-light state S = (q, t, fs, g, ts, wq) of every coupled zone a..b (`cplz(k)`) at the iterate (arrays of any number type), all
from the iterate's r, ρ, T, u, E_r, F_r and species (the opacity χ of the iterate's last `update!` enters the optical depth as a Float):
- optical depth τ_j = Σ_{i≥j} χ_i ρ_i Δr_i from the outer boundary to the inner interface of zone j (τ_{b+1} = 0);
- colour temperature (SLBLEND) T_c = T_src^(1−wq) T_r^wq, T_r = (max(E, E_cmb)/a)^¼, wq = smoothstep(log10(τ/0.3)), T_src = T_r at optical
  depth τ + 2/3 (SLEB: log E_r–log τ interpolation between the bracketing interfaces) or T_eff (not SLEB), and dilution W = 1 (wq ≥ 1) or
  min(c u/(4σT_c⁴), 1), c u = max(c(E − E_cmb), ½(F_k + F_k+1)) (SLWFLUX); q = W/W0, t = T_c/T_c0 against the prepare pass;
- T_eff = (F/σ)^¼ at τ = 2/3, linear in τ between the two interfaces that bracket it (the smooth form of the photosphere zone k_ph);
- Lyman-Werner factor fs = on-weight(τ_k) f_sh(N_H2, T) f_sh(N_HI), the columns N = Σ_{j<k} w_j n_j Δr_j + ½ w_k n_k Δr_k over the zones
  outside the photosphere (weight w_j = sl_lwon(τ_j)), n_H2 = n_H x_H₂, n_HI = n_H (1 − x_H⁺ − 2x_H₂); `shield = :sobolev`: N_H2 ≤ n_H2 v_th/(√2 |∂u/∂r|);
- g, the velocity gradient of Lyα escape (`sl_grad`).
Non-local sums keep the Duals of the zones in the band window of the consumer zone (`winv`); the sums over all other zones enter as values.
`mgqt` (MGROUPS): `mgqt(k, T_c)` fills zone k's group fields q_g (the star light of the zone's chemistry, MGOpac.jl `sl_qt_mg`) and returns the optical group's; T_c, T_src, w_q are the iterate's as in the grey path; f_sh is not used (the H₂ closure).
"""
function sl_sweep!(S, a, b, rn, dn, tn, un, ern, frn, xs1, xs2, cplz; mgqt = nothing)
    sl = STARLIGHT[]; arad = 7.565723e-15; σ = 5.670374419e-5; chi = TM.v(:chifn); isD = eltype(dn) <: ForwardDiff.Dual; CMBE = arad * sl.tcmb^4
    T = eltype(dn); isD1 = T <: ForwardDiff.Dual{<:Any,<:Any,1}; inv_mu = 1 / CNn.μH()
    # optical depth: τ[j] by values, and the Dual corrections (zero value) that a consumer zone k picks up from the zones i that depend on
    # blocks inside its window (zone i = k−3 reaches block k−2 through r_{i+1}; the partials of block k−3 land in a seed slot the band never reads)
    τv = get!(() -> zeros(b + 2), TX.SLQBUF, (:sl_tauv, b + 2))::Vector{Float64}; fill!(τv, 0.0)
    for j in b:-1:a; τv[j] = τv[j+1] + chi[j] * _val(dn[j]) * _val(rn[j+1] - rn[j]); end
    chiD = if isD && TX.MATERIAL[] !== nothing
        z0_ = TX.ZONE0[]; sl0_ = TX.CPLSL[]; TX.CPLSL[] = nothing
        loc = get!(() -> (Vector{T}(undef, length(xs1)), Vector{T}(undef, length(xs1))), TX.SLQBUF, (:locop, T, length(xs1)))::Tuple{Vector{T},Vector{T}}
        chiL, kpL = loc
        if isD1
            neL = get!(() -> zeros(T, length(xs1)), TX.SLQBUF, (:locop_ne, T))::Vector{T}
            length(neL) < length(xs1) && resize!(neL, length(xs1))
        end
        cdv = get!(() -> Vector{T}(undef, b), TX.SLQBUF, (:sl_cdv, T, b))::Vector{T}
        try
            TX.zone_sweep!(a:b, _ -> true) do i
                if isD1
                    m_ = TX.MATERIAL[]
                    if m_ isa NetGas && m_.eqzone[i]
                        Tq = min(tn[i], 1e6)
                        @inbounds neL[i] = (CNn.HECHEM ? (z_ = EQ.composition_he(dn[i], Tq); EQ.composition(dn[i], Tq)[1] + z_[1] + 2z_[2]) : EQ.composition(dn[i], Tq)[1]) * dn[i] * inv_mu
                    else
                        @inbounds neL[i] = TM.v(:xnen)[i] * (dn[i] / _val(dn[i]))
                    end
                end
                ch_, _, kp_ = TX.material_opac(TX.MATERIAL[], dn[i], tn[i])
                @inbounds chiL[i] = ch_; @inbounds kpL[i] = kp_
                vc_ = _val(ch_); @inbounds cdv[i] = (isD1 || !(vc_ > 0)) ? ch_ : chi[i] * (ch_ / vc_)
            end
            for i in 1:a-1; @inbounds cdv[i] = zero(T); end
            cdv
        finally TX.ZONE0[] = z0_; TX.CPLSL[] = sl0_ end
    else nothing end
    corτ(i, k) = (c_ = (chiD === nothing ? chi[i] : winv(chiD[i], i, k)) * winv(dn[i], i, k) * (winv(rn[i+1], i + 1, k) - winv(rn[i], i, k)); c_ - _val(c_))
    τwin(j, k) = j > b ? zero(T) : isD1 ? τv[j] + sum((corτ(i, k) for i in j:b); init = zero(T)) : isD ? τv[j] + sum((corτ(i, k) for i in max(j, k - 3):min(b, k + 2)); init = zero(T)) : τv[j]
    ph = τv[a] >= 2 / 3                                               # (a photosphere inside the grid)
    jph = ph ? findlast(j -> τv[j] >= 2 / 3, a:b) + a - 1 : a         # bracketing interfaces jph, jph+1: τ ≥ 2/3 > τ
    SDPH && (SDPHJ[] != jph && println("SDPH sweep jstep ", Int(TF.getscalar(:jstep)), " thet ", TF.getscalar(:thet), ": iterate jph ", jph - a + 1, " (was ", SDPHJ[] - a + 1, ")"); SDPHJ[] = jph)
    # Lyman-Werner columns: contributions of the zones outside the photosphere, with the Duals of the window and of the transition zones
    lw = sl !== nothing && sl.shield !== :off
    lwcols = lw && sl.shield !== :none
    trans = get!(() -> Int[], TX.SLQBUF, :sl_trans)::Vector{Int}; empty!(trans)
    for j in a:b; 0 < sl_lwon(τv[j]) < 1 && push!(trans, j); end      # zones whose weight depends on τ
    nH(j) = dn[j] * inv_mu
    c2(j, k, w) = w * winv(nH(j), j, k) * max(winv(xs2[j], j, k), 0) * (winv(rn[j+1], j + 1, k) - winv(rn[j], j, k))
    c1(j, k, w) = w * winv(nH(j), j, k) * max(1 - winv(xs1[j], j, k) - 2winv(xs2[j], j, k), 0) * (winv(rn[j+1], j + 1, k) - winv(rn[j], j, k))
    P2 = get!(() -> zeros(b + 2), TX.SLQBUF, (:sl_P2, b + 2))::Vector{Float64}; fill!(P2, 0.0)
    P1 = get!(() -> zeros(b + 2), TX.SLQBUF, (:sl_P1, b + 2))::Vector{Float64}; fill!(P1, 0.0) # value prefix sums of the weighted column integrands over the zones below k
    if lwcols
        for j in a:b
            w = sl_lwon(τv[j]); P2[j+1] = P2[j] + _val(w * nH(j) * max(xs2[j], 0) * (rn[j+1] - rn[j])); P1[j+1] = P1[j] + _val(w * nH(j) * max(1 - xs1[j] - 2xs2[j], 0) * (rn[j+1] - rn[j]))
        end
    end
    function interp_tr(τp, k)                                         # T_r at optical depth τp: log E_r–log τ between the bracketing interfaces
        τ_a = τwin(a, k)
        τp >= τ_a && return (max(winv(ern[a], a, k), CMBE) / arad)^0.25
        τpv = _val(τp)
        for j in a:b-1
            if τv[j] >= τpv > τv[j+1]
                t1, t2 = log(max(τwin(j, k), 1e-30)), log(max(τwin(j + 1, k), 1e-30)); w = t1 > t2 ? (log(τp) - t2) / (t1 - t2) : zero(t1)
                l1, l2 = log(max(winv(ern[j], j, k), CMBE)), log(max(winv(ern[j+1], j + 1, k), CMBE))
                return (exp(l2 + w * (l1 - l2)) / arad)^0.25
            end
        end
        (max(winv(ern[b], b, k), CMBE) / arad)^0.25
    end
    teff(i, k) = (max(winv(frn[i], i, k), zero(frn[i])) / σ)^0.25
    for k in a:b
        cplz(k) || continue
        S.g[k] = sl_grad(rn[k], rn[k+1], un[k], un[k+1])
        TL3DEEP && (S.dp[k] = sl_deepw(τwin(k, k)))
        (k <= length(SLW0[]) && (r = SLW0[][k]) !== nothing) || continue
        E = ern[k]; Fk = frn[k]; Fk1 = frn[k+1]; τk = τwin(k, k)
        Tsrc = if SLEB[]
            interp_tr(τk + 2 / 3, k)
        else                                                              # T_eff at τ = 2/3
            w = (τwin(jph, k) - 2 / 3) / (τwin(jph, k) - τwin(jph + 1, k)); (1 - w) * teff(jph, k) + w * teff(jph + 1, k)
        end
        wq = clamp((log10(max(τk, 1e-30)) + log10(1 / 0.3)) / 1.0, 0, 1); wq = wq * wq * (3 - 2wq)
        cu = SLWFLUX ? max(CNn.c * (E - CMBE), 0.5 * (max(Fk, zero(Fk)) + max(Fk1, zero(Fk1)))) : CNn.c * (E - CMBE)
        Tr = (max(E, CMBE + zero(E)) / arad)^0.25
        Tc = Tsrc^(1 - wq) * Tr^wq
        if mgqt !== nothing              # MGROUPS: the star light is the group fields (mgqt fills every group's q_g; the LW factor is the H₂ closure's, mg_h2lw!)
            S.q[k] = mgqt(k, Tc); S.t[k] = Tc / r.Tc0; S.ts[k] = Tsrc; S.wq[k] = wq; S.fs[k] = zero(Tc)
            continue
        end
        W = wq >= 1 ? one(cu) : min(max(cu, zero(cu)) / (4σ * Tc^4), one(cu))
        S.q[k] = max(W / r.W0, 1e-6); S.t[k] = Tc / r.Tc0; S.ts[k] = Tsrc; S.wq[k] = wq
        if lw
            on = sl_lwon(τk)
            f = one(on) + zero(on)
            if lwcols
                N2 = P2[k] + 0.5 * c2(k, k, on); N1 = P1[k] + 0.5 * c1(k, k, on)
                if isD                                                   # the Duals inside the window and of the transition zones
                    @inbounds for j in (isD1 ? a : max(a, k - 3)):k - 1
                        wj = sl_lwon(τwin(j, k))
                        c2j = c2(j, k, wj); c1j = c1(j, k, wj)
                        N2 += c2j - _val(c2j); N1 += c1j - _val(c1j)
                    end
                    if !isD1
                        @inbounds for j in trans
                            j < k - 3 || continue
                            wj = sl_lwon(τwin(j, k))
                            c2j = c2(j, k, wj); c1j = c1(j, k, wj)
                            N2 += c2j - _val(c2j); N1 += c1j - _val(c1j)
                        end
                    end
                end
                if sl.shield === :sobolev
                    dudr = abs(un[k+1] - un[k]) / (rn[k+1] - rn[k])
                    dudr > 0 && (N2 = min(N2, nH(k) * max(xs2[k], 0) * sqrt(2CNn.kB * tn[k] / 1.6735575e-24) / sqrt(2) / dudr))
                end
                f = fsh_h2(N2, tn[k]) * fsh_hi(N1)
            end
            S.fs[k] = on * f
        end
    end
    S
end
"Per-step star-light setup: photosphere, T_eff, and the LW columns; returns a function k -> (kHm, kH, kLW, fsh) for a network zone k."
function starlight_setup(m::NetGas, a, b, tauz_)
    sl = STARLIGHT[]; dn, tn, rn = TM.v(:dn), TM.v(:tn), TM.v(:rn); er, fr, un = TM.v(:ern), TM.v(:frn), TM.v(:un)
    kph = a; for k in b:-1:a; tauz_[k] >= 2 / 3 && (kph = k; break); end
    L = 4π * rn[kph+1]^2 * fr[kph+1]; Te = L > 0 ? (L / (4π * rn[kph+1]^2 * 5.670374419e-5))^0.25 : 0.0
    SLIMP && CPLPREP[] && (SLSTAR[] = tauz_[kph] >= 2 / 3 && Te > 0)
    SDPH && CPLPREP[] && println("SDPH prep jstep ", Int(TF.getscalar(:jstep)), ": kph ", kph - a_sd() + 1, " SLSTAR ", SLSTAR[], " Te ", round(Te, digits = 1), " tau(z28-34) ", join((@sprintf("%.3f", tauz_[a_sd() + z - 1]) for z in 28:34), " "))
    MGON && (MGKPH[] = kph; length(MGX2[]) == length(dn) || (MGX2[] = zeros(length(dn))))
    arad = 7.565723e-15; Ecmb = arad * sl.tcmb^4
    nz = length(dn); NH2 = zeros(nz); NHI = zeros(nz)
    s2 = 0.0; s1 = 0.0
    for k in kph+1:b
        d = CNn.full(m.x[k], m.net); nH = dn[k] / CNn.μH(); dr = rn[k+1] - rn[k]
        nh2 = nH * max(d[:H2], 0.0); nhi = nH * max(d[:H], 0.0)
        NH2[k] = s2 + 0.5nh2 * dr; NHI[k] = s1 + 0.5nhi * dr
        if MGON
            MGX2[][k] = max(d[:H2], 0.0)
            NH2[k] = mg_ncol(NH2[k], nh2, tn[k], abs(un[k+1] - un[k]) / dr)   # (stage 2c: the same column as mg_h2lw!)
        elseif sl.shield === :sobolev
            dudr = abs(un[k+1] - un[k]) / dr
            if dudr > 0
                vth = sqrt(2CNn.kB * tn[k] / 1.6735575e-24)
                NH2[k] = min(NH2[k], nh2 * vth / sqrt(2) / dudr)     # (N_HI stays static: its Lyman lines shield through Lorentz wings ≫ the infall's Δv)
            end
        end
        s2 += nh2 * dr; s1 += nhi * dr
    end
    νLW = 12.87 * EVERG / CNn.hP
    diag = (; kph, Te, L, kHm = zeros(nz), kH = zeros(nz), kLW = zeros(nz), fsh = ones(nz), tl3 = Vector{Any}(nothing, nz))
    tl3on = hasproperty(sl, :tl3) && sl.tl3
    SLDIAG[] = diag
    # T_r at optical depth τ (log–log interpolation in the zones a..b; τ beyond the innermost zone: its T_r), for SLEB
    function tr_at_tau(τ)
        τ >= tauz_[a] && return (max(er[a], Ecmb) / arad)^0.25
        for j in a:b-1
            if tauz_[j] >= τ > tauz_[j+1]
                t1, t2 = log(max(tauz_[j], 1e-30)), log(max(tauz_[j+1], 1e-30)); w = t1 > t2 ? (log(τ) - t2) / (t1 - t2) : 0.0
                l1, l2 = log(max(er[j], Ecmb, 1e-300)), log(max(er[j+1], Ecmb, 1e-300))
                return (exp(l2 + w * (l1 - l2)) / arad)^0.25
            end
        end
        (max(er[b], Ecmb) / arad)^0.25
    end
    # c × (star-light energy density): E_r above the CMB; with SLWFLUX also at least the zone's outward flux (free streaming: c u = F).
    # Far out (21 pc in c128, 13 Nov) the star's u is ~1e-3 of the CMB's, E_r − E_cmb went negative in one zone and W = 0 there.
    slcu(k) = SLWFLUX ? max(CNn.c * (er[k] - Ecmb), 0.5 * (max(fr[k], 0.0) + max(fr[k+1], 0.0))) : CNn.c * (er[k] - Ecmb)
    function rates(k)
        wq = 0.0
        if SLBLEND[] && Te > 0
            # smooth in τ (no jump at the moving τ = 2/3 zone): colour temperature from T_eff (τ ≤ 0.3) to the local T_r (τ ≥ 3), log blend,
            # dilution W from the local energy density (W = 1 at T_c = T_r). The sharp switch changed x_ss of the zone just outside τ = 2/3
            # by ×2 at equal energy density (diluted 6500 K vs a 4750 K blackbody).
            Tr = (max(er[k], Ecmb) / arad)^0.25
            Tsrc = SLEB[] ? tr_at_tau(tauz_[k] + 2 / 3) : Te
            wq = clamp((log10(max(tauz_[k], 1e-30)) + log10(1 / 0.3)) / 1.0, 0.0, 1.0); wq = wq * wq * (3 - 2wq)
            Tc = Tsrc^(1 - wq) * Tr^wq
            W = wq >= 1 ? 1.0 : min(max(slcu(k), 0.0) / (4 * 5.670374419e-5 * Tc^4), 1.0)
        elseif k > kph && Te > 0
            Tc = Te; W = max(slcu(k), 0.0) / (4 * 5.670374419e-5 * Te^4)
        else
            Tc = (max(er[k], Ecmb) / arad)^0.25; W = 1.0
        end
        ((MGON || W > 0) && Tc > 0) || return (0.0, 0.0, 0.0, 1.0)   # (MGROUPS: the grey W gates nothing)
        Wg = nothing; kHon = false; cu0 = nothing      # MGROUPS: the star light is the group fields (MGOpac.jl mg_cu, mg_wg)
        if MGON
            therm = (SLBLEND[] && Te > 0) ? wq >= 1 : !(k > kph && Te > 0)
            cu0 = mg_cu(k, Ecmb); Wg = mg_wg(cu0, Tc); kHon = !therm && !(TL3DEEP && tauz_[k] >= 3.0)
            CPLPREP[] || mg_wrec!(k, W, Wg, Tc)        # (diagnostic: the grey W next to the group W_g, for mg_record)
        end
        if Wg === nothing
            env = (Trad = Tc, J = ν -> W * _bnu(ν, Tc))
            kHm = sl.pHm === nothing ? 0.0 : CNn.photorate(sl.pHm, tn[k], env); kH = sl.pH === nothing ? 0.0 : CNn.photorate(sl.pH, tn[k], env)
        else
            kHm = sl.pHm === nothing ? 0.0 : CNn.photorates_g(sl.pHm, tn[k], Tc, cu0, MGEHZ[]); kH = kHon ? mg_kion(cu0, tn[k], Tc) : 0.0   # (kH: the record's; the network's own reaction acts, mg_h1s_reactions)
        end
        f = 1.0; kLW = 0.0
        if Wg !== nothing                               # MGROUPS: the LW group's field with the photon-conserving closure (mg_fshJ)
            if k > kph
                fs_, J_, J0_ = mg_fshJ(NH2[k], tn[k]); f = fs_ * J0_ / J_; kLW = mg_klw(cu0[MGROLE[].LW], Tc) * f
            end
        elseif k > kph && sl.shield !== :off
            f = sl.shield === :none ? 1.0 : fsh_h2(NH2[k], tn[k]) * fsh_hi(NHI[k])
            kLW = 1.1e8 * 4π * W * _bnu(νLW, Tc) * f
        end
        diag.kHm[k] = kHm; diag.kH[k] = kH; diag.kLW[k] = kLW; diag.fsh[k] = f
        if SLIMP && CPLPREP[]                            # SLIMP: the zone's star light at the prepare pass (reference of sl_sweep!)
            md = (SLBLEND[] && Te > 0) ? :blend : (k > kph && Te > 0) ? :te : :local
            # no star yet — no zone at τ ≥ 2/3 (the cloud before the first core: kph defaults to the centre zone and T_eff is that of the
            # centre's flux, ~T_CMB) or T_eff = 0: no record, q = t = 1 exactly and the zone takes the SLIMP=0 path (the ICEUL start failed
            # with NaN: Dual zeros in β, P12 at T_c ~ T_CMB, see cpl_tq)
            # (SLW0 is sized by cpl_prepare! before the chemistry: the zones run on CHEMTHREADS tasks, a resize here raced and lost records)
            (SLSTAR[] && k <= length(SLW0[])) && (SLW0[][k] = (; W0 = W, Tc0 = Tc, wq, Tsrc = (md === :blend ? (SLEB[] ? tr_at_tau(tauz_[k] + 2 / 3) : Te) : Tc), mode = md, Ecmb,
                         ρ = dn[k], T = tn[k], fsh = f, lw = k > kph && sl.shield !== :off))
            # MGROUPS: the record carries every group's c u_g at the start of the step and the LW factor (sl_qt_mg, cpl_tqg)
            (MGON && SLSTAR[] && k <= length(SLW0[]) && SLW0[][k] !== nothing) &&
                (SLW0[][k] = merge(SLW0[][k], (; cu0, s0 = f, lw = k > kph, kHon)))
        end
        if tl3on                                        # three-level H parameters (start-of-step state), see tl3_reactions
            d = CNn.full(m.x[k], m.net); nH = dn[k] / CNn.μH(); n1 = nH * max(d[:H], 1e-30)
            g = sl_grad(rn[k], rn[k+1], un[k], un[k+1])
            τS = 3A21LYA * LYA_LAM^3 * n1 / (8π * max(g, 1e-30))
            RL = 3A21LYA * (τS < 1e-6 ? 1.0 : -expm1(-τS) / τS)
            nb = (Wg === nothing ? W : Wg[MGROLE[].Lya]) / expm1(E21 / (CNn.kB * Tc))
            diag.tl3[k] = (; β = (Wg === nothing ? W : Wg[MGROLE[].Ba]) * beta2s_bb(Tc), RL, ne = nH * max(d[:e], 0.0), P12 = (RL + LAM2S) * nb, hν = balmer_mean(Tc), Tc)
            TL3DEEP && (diag.tl3[k] = merge(diag.tl3[k], (; deep = sl_deepw(tauz_[k]))))      # (the start-of-step weight: the operator-split chemistry; the Newton's is the iterate's, sl_sweep!)
            ov_ = LYAOVR[]; (k <= length(ov_) && ov_[k] !== nothing) && (diag.tl3[k] = merge(diag.tl3[k], ov_[k]))   # LYART: transferred Lyα
            Wg === nothing || (diag.tl3[k] = merge(diag.tl3[k], (; mgq = mg_env(cu0, Wg, Tc, kHon))))   # MGROUPS: the photo list's rates from the group fields (tl3env)
        end
        # (MGROUPS with the three-level atom: H⁻ photodetachment is in the network's photo list with the group rates, not in the extra terms,
        # and H(1s) photoionisation with its ground-state recombination is the network's reaction pair mg_h1s_reactions)
        ((Wg !== nothing && tl3on) ? 0.0 : kHm, Wg !== nothing ? 0.0 : kH, kLW, f)
    end
    tl3on && LYART[] && lya_transfer!(m, a, b, tauz_, kph, rates, diag)
    rates
end
include(joinpath(@__DIR__, "LyaTransfer.jl"))
include(joinpath(@__DIR__, "MGOpac.jl"))                # MGROUPS (default off)
"dx/dt of the star-light reactions (rates lagged) for ChemNet.step's `extra`. No heat here: the LW energy goes to TitanResidFast.GASQ (see `lwheat`)."
function starlight_terms(m::NetGas, ρ, kHm, kH, kLW)
    fr_ = m.net.free; iHp = findfirst(==(:Hp), fr_); iHm = findfirst(==(:Hm), fr_); iH2 = findfirst(==(:H2), fr_)
    nH = ρ / CNn.μH()
    extra = (x, xn = nothing) -> begin
        f = zeros(promote_type(eltype(x), typeof(kH), typeof(kLW)), length(x)); d = CNn.full_fast(x, m.net; xn)
        iHm === nothing || (f[iHm] -= kHm * x[iHm])
        f[iHp] += kH * max(d[:H], zero(d[:H]))
        iH2 === nothing || (f[iH2] -= kLW * x[iH2])
        f
    end
    extra
end
"CHEMTHREADS=1: the zones of `chem!(::NetGas)` in parallel (Threads.@threads) when no ray is set (the ray's photons pass the zones in order). Default off."
const CHEMTHREADS = Ref(get(ENV, "CHEMTHREADS", "0") == "1")
const CHEMLOCK = ReentrantLock()
"""
CHEMSCHED=zone (default): the chemistry tasks take zones one at a time from a shared counter. `chunk`: `Threads.@threads`,
which gives every task one contiguous block of zones; the expensive zones (sub-cycled, fallback) are neighbours, so a few
tasks then carry most of the work. The zones are independent: the results are the same.
"""
const CHEMSCHED = Ref(Symbol(get(ENV, "CHEMSCHED", "zone")))
"""
CHEMSUB=fT,fx (default off): a three-level (tl3) network zone whose single backward-Euler chemistry step changes T by more than the factor
fT, or n_e by more than the factor fx (with x_e > 1e-4 at either end), is integrated again in halved sub-steps, each held to the same test
(to depth CHEMSUBMAX), with n_e in the Peebles C refreshed from the current composition between sub-steps.
Why: C (and with it the share 1 − C of collisionally excited n = 2 atoms that the star light ionizes) is frozen at the start-of-step n_e.
A dense zone just out of the table (n 9e13, T 7000 K, x_e 3e-5, β 7e4 /s: C 2e-3) then has an electron avalanche with an e-folding time
of 6 s inside a step of 1e4 s; backward Euler on a growing mode with h λ ≫ 1 has spurious roots, and Newton ended on one at 2.6e6 K, fully
ionized, the recombine/re-ionize loop feeding ⟨hν⟩ per cycle: 5.5e14 erg/g from nothing (full3 step 34704, zone 157, the photosphere zone
in the accretion shock; the radiation spike that followed stopped the run). With the sub-steps the same zones
(n 5e13–2.4e14, T 5300–5500 K at the start of the step) end at x_e 2e-4–1.6e-3 with T 4–5 % lower, and no step heats (test t_cs1, steps
34700–34900); where a single step still converges on the right branch it overshoots x_e by 2–5×.
"""
const CHEMSUB = Ref{Union{Nothing,NTuple{2,Float64}}}(let s = get(ENV, "CHEMSUB", "0"); s == "0" ? nothing : (v = parse.(Float64, split(s, ",")); (v[1], length(v) > 1 ? v[2] : 3.0)) end)
const CHEMSUBMAX = Ref(parse(Int, get(ENV, "CHEMSUBMAX", "16"))); const NCHEMSUB = Ref(0)
"""
`CHEMTMIN=T` (K, default 0 = off): a chemistry step (single or sub-stepped) that returns T below this floor is treated as failed
and the zone takes the table fallback for that step (CHEMFALLBACK). Added 12 Nov after f128 stopped at 46.6 M☉: zone 69 (795 R☉,
n 8.8e11) went 1516 K -> 0.3 K in the sub-stepped chemistry (single step 57 K) and the next hydro solve was NaN. Below the CMB
(55 K at z = 19) such a T is a spurious root of the implicit step, not physics.
"""
const CHEMTMIN = Ref(parse(Float64, get(ENV, "CHEMTMIN", "0"))); const NCHEMTMIN = Ref(0)
"""
CHEMZPROF=1 (diagnostic, default off): wall time of the chemistry per zone; every 500 steps the twelve zones that took most of it.
"""
const CHEMZPROF = get(ENV, "CHEMZPROF", "0") == "1"
const ZPT = Ref(Float64[]); const ZPN = Ref(0)
"The chemistry–temperature step of the H network (`chem!` adds the deuterium trace network with HDCHEM=1)."
function chem_core!(m::NetGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    dn, tn, egn, pgn = TM.v(:dn), TM.v(:tn), TM.v(:egn), TM.v(:pgn)
    dTmax = 0.0
    h = HOT[]
    if PSWRAMP > 0 && RAMPJ0[] < 0 && !CPLPREP[]   # the ramp starts at the outermost table zone below the hot branch
        RAMPJ0[] = Int(TF.getscalar(:jstep)); RAMPK0[] = a - 1
        for k in a:b; m.eqzone[k] && !(h !== nothing && tn[k] > h.T1) && (RAMPK0[] = k); end
        @printf("PSWRAMP: from step %d, table zones from %d inwards leave one every %d steps\n", RAMPJ0[], RAMPK0[] - a + 1, PSWRAMP)
    end
    ray = RAY[]; kph = typemax(Int); Nin = nothing
    if ray !== nothing
        fill!(ray.Q, 0.0)
        st = ray.start(a, b); st === nothing || ((kph, Nin) = st)
    end
    rn = TM.v(:rn)
    tauz = zeros(length(dn))
    if PSWTAU[] > 0 || STARLIGHT[] !== nothing
        chi_ = TM.v(:chifn); tt = 0.0
        for k in b:-1:a; tt += chi_[k] * dn[k] * (rn[k+1] - rn[k]); tauz[k] = tt; end
    end
    slrates = STARLIGHT[] === nothing ? nothing : starlight_setup(m, a, b, tauz)
    # the star-light terms of zone k for ChemNet.step: (extra, cooling) composed with the ray's / the zone's own
    function slterms(k, extra0, cool0)
        slrates === nothing && return extra0, cool0
        kHm, kH, kLW, _ = slrates(k)
        # LOCALJ: the network's own H⁻ photo-process (with its stimulated attachment) sees the same J = W B(T_c); the extra would count it twice
        LOCALJ && ljon(SLDIAG[].tl3[k]) && (kHm = 0.0)
        all(isfinite, (kHm, kH, kLW)) || error("starlight: non-finite rates in zone $k: kHm $kHm kH $kH kLW $kLW, T $(tn[k]), ρ $(dn[k]), E_r $(TM.v(:ern)[k])")
        (kHm == 0 && kH == 0 && kLW == 0) && return extra0, cool0
        ex = starlight_terms(m, dn[k], kHm, kH, kLW)
        ex2 = extra0 === nothing ? ex : (x -> extra0(x) .+ ex(x))
        ex2, cool0
    end
    # LW energy absorbed in zone k over the step, from the end-of-step H₂ (steady under LW): 4.48 eV binding + 0.4 eV fragment kinetic energy
    # per dissociation (UV pumping heat left out); the chemistry is energy-neutral (it takes the binding energy from u and gives it back on re-formation)
    lwheat!(k, x1) = (MGON || slrates === nothing || isempty(TX.GASQ[])) ? nothing :
        (TX.GASQ[][k] = SLDIAG[].kLW[k] * max(CNn.full(x1, m.net)[:H2], 0.0) * dn[k] / CNn.μH() * (4.48 + 0.4) * EVERG; nothing)
    # (COUPLED prepare pass: the hydro step about to run reads these from the last chemistry; they must not be cleared here)
    CPLPREP[] || slrates === nothing || isempty(TX.GASQ[]) || fill!(TX.GASQ[], 0.0)
    CPLPREP[] || isempty(TX.CQB[]) || (fill!(TX.CQB[], 0.0); fill!(TX.YQB[], 0.0))
    CPLPREP[] || ray !== nothing || isempty(TX.RAYQ[]) || fill!(TX.RAYQ[], 0.0)            # (RAYQ also carries the three-level H's photon energy)
    CPLPREP[] || isempty(TX.CHEMB[]) || fill!(TX.CHEMBON[], false)
    ph3 = filter(p -> p.R !== :H, m.photos)                                 # three-level zones: no Boltzmann-level H photo pair
    if NLOPAC[] && !CPLPREP[]                 # (COUPLED prepare pass: keep the non-LTE opacities of the last chemistry for the hydro step)
        length(NLOPC[]) == length(dn) || (NLOPC[] = Any[nothing for _ in dn])
        fill!(NLOPC[], nothing)
        if NLOIT                              # NLOIT: the zones' compositions and three-level parameters of this chemistry instead
            length(NLOST[]) == length(dn) || (NLOST[] = Any[nothing for _ in dn])
            fill!(NLOST[], nothing); nloc_reset!(length(dn))
        end
    end
    dTz = zeros(length(tn))
    function zone!(k)
        CPLPREP[] && (m.eqzone[k] || (h !== nothing && tn[k] > h.T1)) && return nothing     # COUPLED prepare: network zones only, no change of state
        if h !== nothing && tn[k] > h.T1             # hot branch: no chemistry, no re-solve of T with the table (driver: burn and mix)
            m.x[k] = x_from_table(m, dn[k], min(tn[k], 1e6)); m.eqzone[k] = true
            return nothing
        end
        cool = zone_cooling(m, k)
        if Nin !== nothing && k > kph
            V = 4π / 3 * (rn[k+1]^3 - rn[k]^3)
            if sum(Nin) * dt / (dn[k] / CNn.μH() * V) > ray.xmin
                src = ray.source(dn[k], rn[k], rn[k+1], Nin)
                cool2 = (ρ, T, d) -> (cool === nothing ? zero(T) : cool(ρ, T, d)) - src.heat(d)
                m.eqzone[k] = false
                ex3, cool3 = slterms(k, src.extra, cool2)
                x1, T1 = try
                    CNn.step(m.x[k], dn[k], egn[k], tn[k], dt, m.net; env = ray.env, photos = m.photos, cooling = cool3, extra = ex3, accepted = src.accepted)
                catch e
                    println("ray zone failed: k ", k, " x ", m.x[k], " free ", m.net.free, " ρ ", dn[k], " u ", egn[k], " T ", tn[k], " dt ", dt, " Nin ", sum(Nin), " r ", rn[k], " ", rn[k+1])
                    rethrow()
                end
                lwheat!(k, x1)
                out = src.out(x1, dt)
                ray.Q[k] = sum((Nin .- out) .* ray.Ec) / V
                m.x[k] = x1; Nin = out
                dTz[k] = abs(T1 / tn[k] - 1); tn[k] = T1
                return nothing
            end
            Nin = nothing                                # the rest of the light cannot change any zone further out
        end
        # PSWTMIN[] > 0: also require T above it (5 % hysteresis). Below ~5000 K ionization holds negligible energy, so the network and the table agree there
        # to ~1 % in T at fixed e; at 7000-9500 K the network's H ionization lags the table's Saha by ×5, and switching there drops T by up to 16 %.
        ptab = (pgn[k] >= m.pswitch && tn[k] >= PSWTMIN[]) || (m.eqzone[k] && pgn[k] >= PSWHYS[] * m.pswitch && tn[k] >= 0.95 * PSWTMIN[])
        PSWNH > 0 && !ptab && (ptab = dn[k] / CNn.μH() >= (m.eqzone[k] ? PSWHYS[] : 1.0) * PSWNH && tn[k] >= (m.eqzone[k] ? 0.95 : 1.0) * PSWTMIN[])
        tab = ptab || (PSWSTAY[] && m.eqzone[k] && PSWTEX[] > 0 && tn[k] >= PSWTEX[]) || (PSWTAU[] > 0 && tn[k] >= PSWTEX[] && tauz[k] >= (m.eqzone[k] ? 0.5 : 1.0) * PSWTAU[])
        G2B !== nothing && (tab = g2_w(dn[k], tn[k]) >= (m.eqzone[k] ? 0.99 : 0.999))     # G2BLEND: the table only where the blend is complete
        CPLPREP[] && (tab = false)
        PSWRAMP > 0 && !CPLPREP[] && m.eqzone[k] && k <= RAMPK0[] - (Int(TF.getscalar(:jstep)) - RAMPJ0[]) ÷ PSWRAMP && (tab = true)
        if tab
            coolt = cool === nothing ? nothing : ((ρ, T, y) -> cool(ρ, T, nothing))
            (CPLPOST[] && TX.cplz(k)) && (coolt = nothing)      # COUPLED: the zone's cooling over this step was in the hydro solve
            _, T1 = EQ.step(dn[k], egn[k], tn[k], dt; cooling = coolt)
            m.x[k] = x_from_table(m, dn[k], T1); m.eqzone[k] = true
        else
            eqprev = m.eqzone[k]
            # (back to the ideal regime: the network continues from the table composition). With the three-level atom, from the table composition
            # at the zone's current T, not m.x: m.x of a table zone is still advected, and took in the ionized network gas upstream (x_e 0.3–0.4)
            # while the table's e had radiated that energy away; the network then started 3e12 erg/g below u(10 K, x) (zones at 1800–1900 R_sun).
            eqprev && STARLIGHT[] !== nothing && hasproperty(STARLIGHT[], :tl3) && STARLIGHT[].tl3 && (m.x[k] = x_from_table(m, dn[k], tn[k]))
            m.eqzone[k] && (m.eqzone[k] = false)
            exs, cools = slterms(k, nothing, cool)
            p3 = (slrates === nothing || m.eqzone[k]) ? nothing : SLDIAG[].tl3[k]
            envk, phk = m.env, m.photos; fb = false
            cools0 = cools                           # (SLIMP) the zone's cooling without the three-level heat
            p3e = p3                                 # the three-level parameters at the end of the chemistry (n_e refreshed if CHEMSUB sub-stepped)
            if p3 !== nothing
                envk = tl3env(m, p3); phk = ph3
                c0 = cools
                cools = (ρ, T, d) -> (c0 === nothing ? zero(T) : c0(ρ, T, d)) - tl3_heat(ρ, T, d, p3)[1]
                # CHEMCKHSS: n_e of the three-level heat from the current composition, as in ChemNet.step_ck's C
                if (CNn.CK[] && CNn.CKHSS[] > 0) || TX.COUPLED[]
                    cools = (ρ, T, d) -> (c0 === nothing ? zero(T) : c0(ρ, T, d)) - tl3_heat(ρ, T, d, merge(p3, (; ne = ρ / CNn.μH() * max(d[:Hp], zero(d[:Hp])))))[1]
                end
            end
            if CNn.HM88                              # HM88: perturber densities of the zone (composition at the start of the step) for φ
                dh_ = CNn.full(m.x[k], m.net); nh_ = dn[k] / CNn.μH()
                hm_ = (nh_ * (max(dh_[:H], 0.0) + max(dh_[:He], 0.0) + max(dh_[:H2], 0.0)), nh_ * max(dh_[:e], 0.0))
                envk = merge(envk, (; hm = hm_, tabkey = hasproperty(envk, :tabkey) ? envk.tabkey : envk))
            end
            if !isempty(TX.CHEMB[])                  # BALSPLIT: the hydro step already gave the zone c̄ (TitanResidFast.CHEMB)
                cb = eqprev ? 0.0 : TX.CHEMB[][k]; TX.CHEMBON[][k] = !eqprev
                if cb != 0
                    c1 = cools; cools = (ρ, T, d) -> (c1 === nothing ? zero(T) : c1(ρ, T, d)) + ρ * cb
                end
            end
            if CPLDBGZ[] > 0 && k - a + 1 == CPLDBGZ[] && CPLDBGN[] > 0      # (diagnostic) the zone's cooling terms at its current state
                CPLDBGN[] -= 1; d0 = CNn.full(m.x[k], m.net)
                c1 = cool === nothing ? 0.0 : cool(dn[k], tn[k], d0); tl = p3 === nothing ? (0.0, 0.0) : tl3_heat(dn[k], tn[k], d0, p3)
                ex = exs === nothing ? nothing : exs(m.x[k])
                @printf("CPLDBGZ %s step %d zone %d r %.1f T %.1f ρ %.3e x_e %.3e: zone_cooling %.4e  tl3_heat %.4e %.4e  cools %.4e  extra %s  dt %.3e  u %.6e\n",
                        CPLPREP[] ? "prep" : CPLPOST[] ? "post" : "split", TF.getscalar(:jstep), k - a + 1, rn[k] / 6.957e10, tn[k], dn[k], d0[:e], c1, tl[1], tl[2],
                        cools === nothing ? 0.0 : cools(dn[k], tn[k], d0), ex === nothing ? "-" : string(round.(ex, sigdigits = 3)), dt, egn[k])
            end
            if CPLPREP[]                             # COUPLED prepare: the zone's terms for the hydro solve
                envc = envk; phc_ = phk; cache = Ref((NaN, Tuple{CNn.Photo,Float64,Float64}[]))
                phc = photo_cache(cache, phc_, envc)
                if SLIMP && p3 !== nothing && isempty(TX.CHEMB[])   # SLIMP: what scales with W (see sl_q, cpl_tq)
                    lj = !MGON && (hasproperty(envk, :ljW) || (ljon(p3) && hasproperty(envk, :J)))
                    envc0 = !lj ? nothing : hasproperty(envk, :ljW) ? merge(envk, (; J = ν -> 0.0, ljW = 0.0)) : merge(envk, (; J = envk.ljs.J))      # (W = 0: the static field alone)
                    cache0 = Ref((NaN, Tuple{CNn.Photo,Float64,Float64}[]))
                    phc0 = envc0 === nothing ? nothing :
                        photo_cache(cache0, phc_, envc0)
                    hm0 = if CNn.HM88                        # He and the He-ion electrons of the start of the step: the part of the perturber densities the slow species do not carry
                        dh0 = CNn.full(m.x[k], m.net); (max(dh0[:He], 0.0), max(dh0[:e] - dh0[:Hp], 0.0))
                    else nothing end
                    slq = (; p3, lj, envc0, phc0, hm0, c0 = cools0, ph = phc_, net = m.net, mm = m,
                           rec = (k <= length(SLW0[]) ? SLW0[][k] : nothing))
                    CPLTERMS[][k] = (; env = envk, phc, cooling = cools, extra = exs, phraw = (phc_, envc), lcool = cool, slq)
                    return nothing
                end
                CPLTERMS[][k] = (; env = envk, phc, cooling = cools, extra = exs, phraw = (phc_, envc), lcool = cool)   # (phraw: NLOIT's uncached photo-rates; lcool: the line cooling alone)
                return nothing
            end
            if CPLPOST[] && TX.cplz(k) && p3 !== nothing
                p3e = merge(p3, (; ne = dn[k] / CNn.μH() * max(CNn.full(m.x[k], m.net)[:Hp], 0.0)))
            end
            x1, T1 = (CPLPOST[] && TX.cplz(k)) ? (m.x[k], tn[k]) : try
                # u below the zone's own composition at 10 K (x and e advected by different schemes; with x_e ~ 0.1-0.5 from the three-level
                # atom the mismatch reaches ~1 eV per H): no solution for this x, take the fallback below
                p3 !== nothing && egn[k] < CNn.eos(dn[k], 10.0, m.x[k], m.net)[1] && error("tl3: u below u(10 K, x) in zone $k")
                xs1, Ts1 = CNn.step(m.x[k], dn[k], egn[k], tn[k], dt, m.net; env = envk, photos = phk, cooling = cools, extra = exs)
                ckhss = CNn.CK[] && CNn.CKHSS[] > 0
                ckhss && p3 !== nothing && (p3e = merge(p3, (; ne = dn[k] / CNn.μH() * max(CNn.full(xs1, m.net)[:Hp], 0.0))))
                if p3 !== nothing && CHEMSUB[] !== nothing && isempty(TX.CHEMB[]) && !ckhss
                    fT, fx = CHEMSUB[]; nHk = dn[k] / CNn.μH(); xel(x) = max(CNn.full(x, m.net)[:e], 0.0)
                    bad(xa, Ta, xb, Tb) = (ea = xel(xa); eb = xel(xb); !(isfinite(Tb) && all(isfinite, xb)) || abs(log(Tb / Ta)) > log(fT) ||
                                           (max(ea, eb) > 1e-4 && max(ea, eb) > fx * min(ea, eb)))
                    if bad(m.x[k], tn[k], xs1, Ts1)
                        dmax = Ref(0)
                        function sub(xa, Ta, ua, hh, depth, pa)
                            dmax[] = max(dmax[], depth)
                            ca = (ρ, T, d) -> (c0 === nothing ? zero(T) : c0(ρ, T, d)) - tl3_heat(ρ, T, d, pa)[1]
                            xb, Tb = CNn.step(xa, dn[k], ua, Ta, hh, m.net; env = tl3env(m, pa), photos = phk, cooling = ca, extra = exs)
                            (depth >= CHEMSUBMAX[] || !bad(xa, Ta, xb, Tb)) && return (xb, Tb, merge(pa, (; ne = nHk * xel(xb))))
                            xc, Tc_, pc = sub(xa, Ta, ua, hh / 2, depth + 1, pa)
                            sub(xc, Tc_, Float64(CNn.eos(dn[k], Tc_, xc, m.net)[1]), hh / 2, depth + 1, pc)
                        end
                        xs0, Ts0 = xs1, Ts1
                        xm, Tm, pm = sub(m.x[k], tn[k], Float64(egn[k]), dt / 2, 1, p3)
                        xs1, Ts1, p3e = sub(xm, Tm, Float64(CNn.eos(dn[k], Tm, xm, m.net)[1]), dt / 2, 1, pm)
                        lock(() -> (NCHEMSUB[] += 1), CHEMLOCK)
                        (dmax[] >= 4 || TL3DBG[] > 0) && @printf("chemsub k %d step %d r %.3e R_sun n %.2e dt %.2e depth %d: T %.1f -> %.1f (one step: %.1f), x_e %.2e -> %.2e (one step: %.2e)\n", k, TF.getscalar(:jstep), rn[k] / 6.957e10, nHk, dt, dmax[], tn[k], Ts1, Ts0, xel(m.x[k]), xel(xs1), xel(xs0))
                    end
                end
                CHEMTMIN[] > 0 && Ts1 < CHEMTMIN[] && (lock(() -> (NCHEMTMIN[] += 1), CHEMLOCK); error("chem: T $(Ts1) below CHEMTMIN in zone $k (from $(tn[k]))"))
                (STARLIGHT[] !== nothing && !(all(isfinite, xs1) && isfinite(Ts1))) && error("starlight: non-finite chemistry in zone $k: x $xs1 T $Ts1 from x $(m.x[k]) T $(tn[k]) ρ $(dn[k]) rates $(slrates(k))")
                (xs1, Ts1)
            catch e
                # Fixed mesh: e and x are advected by different schemes, and e may then have no solution for the zone's x. Hot dense gas: table equilibrium, same e.
                STARLIGHT[] !== nothing && !(e isa ErrorException) && println("starlight zone failed: k ", k, " x ", m.x[k], " ρ ", dn[k], " u ", egn[k], " T ", tn[k], " dt ", dt, " rates ", slrates(k), " eq ", m.eqzone[k], " E_r ", TM.v(:ern)[k])
                ((CHEMFALLBACK[] || (CHEMFBDENSE[] && pgn[k] >= m.pswitch)) && (e isa ErrorException || (p3 !== nothing && e isa ArgumentError))) || rethrow()
                lock(() -> (NCHEMFB[] += 1), CHEMLOCK); fb = true
                TL3DBG[] > 0 && @printf("tl3fb k %d r %.3e T %.0f tauz %.2e eqprev %s xe %.2e u %.4e u(10K,x) %.4e u(T,x) %.4e: %s\n", k, rn[k] / 6.957e10, tn[k], tauz[k], eqprev, CNn.full(m.x[k], m.net)[:e], egn[k], CNn.eos(dn[k], 10.0, m.x[k], m.net)[1], CNn.eos(dn[k], tn[k], m.x[k], m.net)[1], first(sprint(showerror, e), 60))
                if CHEMFBNET > 0 && dn[k] < EQ.rho_min()
                    # below the equilibrium table: keep the network (CHEMFBNET sub-steps; failing that, the composition held and T from u)
                    lock(() -> (NCHEMFBNET[] += 1), CHEMLOCK)
                    chem_fallback_net(m, k, dn[k], Float64(egn[k]), tn[k], dt, envk, phk, cools, exs)
                else
                    coolt = cool === nothing ? nothing : ((ρ, T, y) -> cool(ρ, T, nothing))
                    _, T2 = EQ.step(dn[k], egn[k], tn[k], dt; cooling = coolt)
                    (x_from_table(m, dn[k], T2), T2)
                end
            end
            fb && !isempty(TX.CHEMB[]) && (TX.CHEMBON[][k] = false)
            m.eqzone[k] || lwheat!(k, x1)
            if p3 !== nothing && !m.eqzone[k] && !fb
                cq, yq = tl3_heat(dn[k], T1, CNn.full(x1, m.net), p3e)
                isempty(TX.CQB[]) || (TX.CQB[][k] = cq)          # booked for the next step's rows of a zone that is not coupled there (TX.CQB)
                isempty(TX.YQB[]) || (TX.YQB[][k] = yq)
                SLDIAG[].tl3[k] = merge(p3e, (; chemq = cq, cycq = yq, qbb = lya_bb(dn[k], T1, CNn.full(x1, m.net), p3e)))
                NLOPAC[] && NLOIT && (NLOST[][k] = (; d = CNn.full(x1, m.net), p = p3e))
                if TL3DBG[] > 0 && TF.getscalar(:jstep) % TL3DBG[] == 0 && k % 6 == 0
                    d0 = CNn.full(m.x[k], m.net); d1 = CNn.full(x1, m.net)
                    @printf("tl3dbg k %d r %.3e tauz %.2e T %.0f -> %.0f  xe %.2e -> %.2e  C %.2e β %.2e RL %.2e  chemq %.2e cycq %.2e  u %.3e  u(T0,x0) %.3e  eqprev %s\n", k, rn[k] / 6.957e10, tauz[k], tn[k], T1, d0[:e], d1[:e], tl3_C(T1, p3), p3.β, p3.RL, cq, yq, egn[k], CNn.eos(dn[k], tn[k], m.x[k], m.net)[1], eqprev)
                end
            end
            m.x[k] = x1
        end
        dTz[k] = abs(T1 / tn[k] - 1); tn[k] = T1
        nothing
    end
    CHEMZPROF && length(ZPT[]) != length(tn) && (ZPT[] = zeros(length(tn)))
    zrun!(k) = CHEMZPROF ? (t0_ = time_ns(); zone!(k); ZPT[][k] += (time_ns() - t0_) * 1e-9; nothing) : zone!(k)   # CHEMZPROF: time per zone (one task per zone, no race)
    if CHEMTHREADS[] && ray === nothing && Threads.nthreads() > 1    # zones are independent (start-of-step state, per-zone outputs)
        if CHEMSCHED[] === :zone
            next = Threads.Atomic{Int}(a)
            @sync for _ in 1:Threads.nthreads()
                Threads.@spawn while true
                    kz = Threads.atomic_add!(next, 1)
                    kz > b && break
                    zrun!(kz)
                end
            end
        else
            Threads.@threads :dynamic for k in a:b; zrun!(k); end
        end
    else
        for k in a:b; zrun!(k); end
    end
    dTmax = maximum(view(dTz, a:b))
    if CHEMZPROF && (ZPN[] += 1) % 500 == 0           # the zones that took most of the network time over the last 500 steps
        zt = ZPT[]; tot = sum(view(zt, a:b)); js = sort(collect(a:b); by = k -> -zt[k])[1:min(12, b - a + 1)]
        @printf("CHEMZPROF step %d: %.1f s in zones over 500 steps (%d table, %d network now); top zones:\n", TF.getscalar(:jstep), tot, count(k -> m.eqzone[k], a:b), count(k -> !m.eqzone[k], a:b))
        for k in js
            @printf("    zone %3d  %5.1f %%  r %.3e R_sun  n %.3e  T %8.1f  p/pswitch %.3g  %s\n", k - a + 1, 100 * zt[k] / tot, TM.v(:rn)[k] / 6.957e10, dn[k] / CNn.μH(), tn[k], pgn[k] / m.pswitch, m.eqzone[k] ? "table" : "network")
        end
        fill!(zt, 0.0)
    end
    for k in a-2:a-1; k >= 1 && (m.x[k] = m.x[a]; m.eqzone[k] = m.eqzone[a]); end
    for k in b+1:b+2; m.x[k] = m.x[b]; m.eqzone[k] = m.eqzone[b]; end
    tn[a-1] = tn[a]; tn[b+1] = tn[b]; tn[b+2] = tn[b]
    TM.fill_eos!(m)
    dTmax
end

# ------------------------------------------------------------------------------------------ deuterium and HD (HDCHEM=1)
"""
HDCHEM=1 (default off): D, D⁺ and HD as a trace network carried with the H network, and optically thin HD cooling.

State per zone: (x_D⁺, x_HD) per H nucleus, appended to the zone's `m.x[k]` behind the free species of the network (the advection
and every remap treat `m.x[k]` as a vector; `ChemNet.full` and `fullvec` read only the free species). x_D = y_D − x_D⁺ − x_HD with
y_D = D/H by number (ENV `DH`, the variable the nuclear network uses; 2.5e-5).
After the chemistry–temperature step of a zone the three-state linear system takes one backward-Euler step at the zone's new T and
composition: rates k50–k56 of ChemistryKernels.jl (Savin 2002; Galli & Palla 2002; Clark et al. 2011; Kreckel et al. 2010 form for
D + H⁻), D ionisation and recombination with the H rates, HD photodissociation with the zone's Lyman–Werner rate. The step is written
as the stationary state of the chain with a restart at rate 1/Δt to the old state: sums of positive products only, exact for any
rate × Δt, and y_D is conserved. Table and hot zones: the stationary state.
Cooling: n_HD [W(T) − W(T_CMB)], W = Λ_LTE/(1 + Λ_LTE/(n_HI Λ_low)) (Coppola et al. 2011; Wrathmall, Gusdorf & Flower 2007: HD–H
collisions only), inside the implicit chemistry–temperature step, with the x_HD of the start of the step. The coupled zones of the hydro Newton
have the x_HD of the iterate instead (`cpl_g`): the backward-Euler step of the chain from the start-of-step state over the time elapsed in the stage, at the
iterate's T, ρ, composition and Lyman–Werner rate (local, exact derivatives).
Zones without the two entries (restart from a run without HD, a rebuilt mesh): table zones get the stationary state; network zones
the stationary state with `HDINIT=eq` (default), or with `HDINIT=<x_HD>` that HD abundance and D⁺/D = H⁺/H.
Not included: HD in the EOS (binding energy, particle number: 2.5e-5 of H), HD–H₂ and HD–He collisions, HD line optical depths,
pumping of HD by the star, the nuclear network's deuterium (the star's zones are fully ionised, so their D is D⁺ here either way).
"""
const HDCHEM = Ref(get(ENV, "HDCHEM", "0") == "1")
const HDYD = Ref(parse(Float64, get(ENV, "DH", "2.5e-5")))
const HDINIT = Ref(get(ENV, "HDINIT", "eq"))
const HDTCMB = Ref(parse(Float64, get(ENV, "TCMB", "0")))
const HDX = Ref(Float64[])                 # per zone: x_HD at the start of the step (the cooling term of zone_cooling)
const HDLAM = Ref(Float64[])               # per zone: HD cooling at the end of the step [erg cm⁻³ s⁻¹] (diagnostic)
const HDNPAD = Ref(0)
hd_k50(T) = T <= 2e5 ? max(0.0, 2.0e-10 * T^0.402 * exp(-37.1 / T) - 3.31e-17 * T^1.48) : 2.5e-8 * (T / 2e5)^0.402      # H⁺ + D → H + D⁺
hd_k51(T) = 2.06e-10 * T^0.396 * exp(-33.0 / T) + 2.03e-9 * T^-0.332                                                    # H + D⁺ → H⁺ + D
hd_k52(T) = T <= 1e4 ? (l = log10(T); max(0.0, 1.0e-9 * (0.417 + 0.846l - 0.137l^2))) : 1.609e-9                        # H₂ + D⁺ → HD + H⁺
hd_k53(T) = 1.1e-9 * exp(-488.0 / T)                                                                                    # HD + H⁺ → H₂ + D⁺
hd_k54(T) = T <= 2e3 ? (l = log10(T); exp10(-56.4737 + l * (5.88886 + l * (7.19692 + l * (2.25069 + l * (-2.16903 + l * 0.317887)))))) :
                       3.17e-10 * exp(-5207.0 / T)                                                                      # H₂ + D → HD + H
# (the published fit turns upward below 78.5 K: below 200 K an Arrhenius tail matched in value and slope, as ChemistryKernels.jl)
hd_k55(T) = T <= 200 ? exp(log(5.25e-11) - 4430.0 / 200 + 173900.0 / 200^2 - (4430.0 - 2 * 173900.0 / 200) * (1 / T - 1 / 200)) :
                       5.25e-11 * exp(-4430.0 / T + 173900.0 / T^2)                                                     # HD + H → H₂ + D
hd_k56(T) = 1.35e-9 * (T^0.098493 + 0.32852 * T^0.5561 + 2.771e-7 * T^2.1826) /
            (1 + 6.191e-3 * T^1.0461 + 8.9712e-11 * T^3.0424 + 3.2576e-14 * T^3.7741)                                   # D + H⁻ → HD + e
"HD cooling per molecule in LTE [erg/s] (Coppola et al. 2011) and the low-density coefficient for HD–H [erg cm³/s] (Wrathmall et al. 2007)."
hd_lte(T) = (t = min(max(T, 10.0), 3.0e4); l = log10(t); exp10(min(-55.5725 + l * (56.649 + l * (-37.9102 + l * (12.698 + l * (-2.02424 + l * 0.122393)))), 0.0)))
hd_low(T) = (t = min(max(T, 10.0), 6.0e3); l = log10(t / 1e3); exp10(-23.175780 + l * (1.5035261 + l * (0.40871403 + l * (0.17849311 + l * (-0.077291388 + l * 0.10031326))))))
hd_w(T, nHI) = (l = hd_lte(T); lo = nHI * hd_low(T); (l > 0 && lo > 0) ? l / (1 + l / lo) : zero(l * lo))
"HD cooling [erg cm⁻³ s⁻¹] of gas with composition d (ChemNet.full; nothing: table zone, no HD) and x_HD; < 0 below T_CMB."
function hd_cool(ρ, T, d, xHD)
    d === nothing && return zero(T)
    nH = ρ / CNn.μH(); nHI = nH * max(d[:H], 0.0); Tc = HDTCMB[]
    nH * xHD * (hd_w(T, nHI) - (Tc > 0 ? hd_w(Tc, nHI) : 0.0))
end
"Transition rates [1/s] between the states (D, D⁺, HD): (D→D⁺, D⁺→D, D⁺→HD, HD→D⁺, D→HD, HD→D)."
function hd_rates(T, nH, d, kLW)
    t = min(max(T, 10.0), 1e6)
    g(s) = max(get(d, s, zero(t)), 0.0)
    xH, xp, xe, xH2, xm = g(:H), g(:Hp), g(:e), g(:H2), g(:Hm)
    (nH * (PC.k1_aazn(t) * xe + hd_k50(t) * xp), nH * (PC.k2_aazn(t) * xe + hd_k51(t) * xH), nH * hd_k52(t) * xH2,
     nH * hd_k53(t) * xp, nH * (hd_k54(t) * xH2 + hd_k56(t) * xm), nH * hd_k55(t) * xH + kLW)
end
"""
One backward-Euler step of the chain from p0 = fractions (D, D⁺, HD) of y_D, idt = 1/Δt (0: the stationary state). p(I − Δt Q) = p0 is
the stationary state of the chain with every state restarting at rate 1/Δt into p0; for three states that is a sum over spanning trees.
"""
function hd_step(q, p0, idt)
    q12, q21, q23, q32, q13, q31 = q
    if idt > 0
        isfinite(idt) || return p0
        q21 += idt * p0[1]; q31 += idt * p0[1]; q12 += idt * p0[2]; q32 += idt * p0[2]; q13 += idt * p0[3]; q23 += idt * p0[3]
    end
    w1 = q21 * q31 + q21 * q32 + q23 * q31
    w2 = q12 * q32 + q12 * q31 + q32 * q13
    w3 = q13 * q23 + q13 * q21 + q23 * q12
    s = w1 + w2 + w3
    (isfinite(s) && s > 0) || return p0
    (w1 / s, w2 / s, w3 / s)
end
"""
The (D, D⁺, HD) fractions of the coupled zone at the iterate (HDCHEM): the backward-Euler step of the three-state chain (`hd_step`) from the start-of-step
fractions p0 over the time elapsed since the start of the step (stage: γΔt, Δt) at the iterate's T, ρ, composition d and Lyman–Werner rate kLW; at the end of the
step it is the step the operator-split chemistry (`chem!`) takes from the same state. Generic in the number type.
"""
hd_x(ρ, T, d, kLW, p0) = (te = TF.getscalar(:timen) - TF.getscalar(:timeo); hd_step(hd_rates(T, ρ / CNn.μH(), d, kLW), p0, te > 0 ? 1 / te : Inf))
hd_fracs(xDp, xHD) = (yD = HDYD[]; p2 = clamp(xDp / yD, 0.0, 1.0); p3 = clamp(xHD / yD, 0.0, 1.0 - p2); (1.0 - p2 - p3, p2, p3))
function hd_default(m, k, x)
    dn, tn = TM.v(:dn), TM.v(:tn); yD = HDYD[]
    d = CNn.full(x, m.net); nH = dn[k] / CNn.μH()
    if m.eqzone[k] || HDINIT[] == "eq"
        p = hd_step(hd_rates(tn[k], nH, d, 0.0), (1.0, 0.0, 0.0), 0.0)
        return (yD * p[2], yD * p[3])
    end
    xhd = min(parse(Float64, HDINIT[]), yD); xH = max(Float64(d[:H]), 0.0); xp = max(Float64(d[:Hp]), 0.0)
    ((yD - xhd) * (xH + xp > 0 ? xp / (xH + xp) : 0.0), xhd)
end
"Give every zone its (x_D⁺, x_HD) entries (no-op where they are present)."
function hd_pad!(m::NetGas)
    nf = length(m.net.free); n = 0
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    for k in a:b
        length(m.x[k]) == nf || continue
        p = hd_default(m, k, m.x[k]); m.x[k] = vcat(m.x[k], p[1], p[2]); n += 1
    end
    for k in eachindex(m.x)                                        # ghosts and the unused tail of the arrays: as the nearest zone
        (a <= k <= b || length(m.x[k]) != nf) && continue
        j = k < a ? a : b; m.x[k] = vcat(m.x[k], m.x[j][nf+1], m.x[j][nf+2])
    end
    n > 0 && (HDNPAD[] += n; @printf("HDCHEM: (x_D+, x_HD) set in %d zones (HDINIT=%s, D/H = %.3g)\n", n, HDINIT[], HDYD[]))
    m
end
"""
Restart from a checkpoint written without `HECHEM`: insert the He⁺ and He²⁺ fractions (Saha at the zone's ρ, T from the
equilibrium table) in front of the HD tail of every zone's vector (no-op where they are present). Ghost entries use the nearest zone.
"""
function he_pad!(m::NetGas)
    CNn.HECHEM || return m
    i1 = findfirst(==(:Hep), m.net.free); i1 === nothing && return m
    nf = length(m.net.free); want = nf + (HDCHEM[] ? 2 : 0); a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dn, tn = TM.v(:dn), TM.v(:tn); n = 0
    for k in eachindex(m.x)
        length(m.x[k]) == want - 2 || continue
        j = clamp(k, a, b); z = EQ.composition_he(dn[j], min(tn[j], 1e6))
        m.x[k] = vcat(m.x[k][1:i1-1], z[1], z[2], m.x[k][i1:end]); a <= k <= b && (n += 1)
    end
    n > 0 && @printf("HECHEM: (x_He+, x_He++) set from Saha in %d zones\n", n)
    m
end
hd_strip!(m::NetGas) = (nf = length(m.net.free); for k in eachindex(m.x); length(m.x[k]) > nf && (m.x[k] = m.x[k][1:nf]); end; m)
"Per-zone HD record for the snapshots (zones a:b): x_D⁺, x_HD and the HD cooling."
function hd_snap(m::NetGas, a, b)
    hd_pad!(m); nf = length(m.net.free); lam = HDLAM[]
    (; Dp = [m.x[k][nf+1] for k in a:b], HD = [m.x[k][nf+2] for k in a:b], Λ = length(lam) >= b ? lam[a:b] : zeros(b - a + 1))
end
# ------------------------------------------------------------------------------------------ COUPLED (TitanResidFast.COUPLED)
"""
COUPLED=1 (popiii_m2c; default off): x_H⁺ and x_H₂ of the network zones are unknowns of the hydro Newton (TitanResidFast.COUPLED).
Per step: `cpl_prepare!` (before the hydro step) runs `chem_core!` in a prepare mode that changes no state and records per network zone
the terms the zone's chemistry step would use (env with the three-level parameters, photo-rates, cooling closure including the
three-level heat with self-consistent n_e, star-light `extra`), loads the iterate/old species and the old rates and cooling.
`cpl_zone` gives (p, e, S_H⁺, S_H₂, Λ) of a zone at (ρ, T, x_H⁺, x_H₂) with a 4-wide local ForwardDiff Jacobian, lifted onto the
residual's dual numbers (cached per zone on the Float values). After the step `cpl_chem!` sets the coupled zones' free species to the
CK expansion of the solved species and runs `chem!` in a post mode: those zones skip the network step (their T and x are the hydro
solution, their cooling was in the energy row), the rest is unchanged (HD step, three-level bookkeeping RAYQ/GASQ/SLDIAG/NLOPC, the
table switch — with no second cooling for a zone that leaves the network —, table and hot zones as before).
Not coupled: the IonRay path (a ray zone is not prepared), CHEMB (BALSPLIT) — both off in the runs this was built for.
"""
const CPLPREP = Ref(false); const CPLPOST = Ref(false)
const CPLDBGZ = Ref(parse(Int, get(ENV, "CPLDBGZ", "0"))); const CPLDBGN = Ref(parse(Int, get(ENV, "CPLDBGN", "20")))
const CPLTERMS = Ref(Any[]); const CPLCACHE = Ref(Any[])
"(HDCHEM) per zone the start-of-step fractions (D, D⁺, HD) of the coupled zones' HD chain (hd_fracs), the old state of `hd_x`."
const CPLHD0 = Ref(Any[])
"""
The zone's photo-rates as a function of T: cached on the last Float64 T (`cache = Ref((T, rates))`); a Dual T (the local Jacobian sweeps of
`cpl_vj`, `nlo_slopes`) is evaluated directly and carries ∂(rates)/∂T (ChemNet.ck_eval used to take the rates at the value of T: the Jacobian then lacked
the T dependence of the photo-ionisation/recombination coefficients, e.g. 0.03 instead of 0.12 for d ln n_e/d ln T of the He ions in zones 17–19 of the stuck state).
"""
photo_cache(cache, ph_, env) = Tv -> (Tv isa ForwardDiff.Dual || (CNn.HM88 && hasproperty(env, :hm) && env.hm[1] isa ForwardDiff.Dual)) ? _phvec(ph_, Tv, env) :
    (Tv == cache[][1] || (cache[] = (Tv, _phvec(ph_, Tv, env))); cache[][2])
# The neutral hydrogen of the slow species, 1 − x_H⁺ − 2x_H₂, with the value `xn` given directly (the state carries it: x_H⁺ = 1 − 1e-9 cannot hold a
# neutral fraction to better than 1e-7) and the derivatives of 1 − x₁ − 2x₂ (−1, −2): the value is replaced, not recomputed. `xn = nothing`: 1 − x₁ − 2x₂.
@inline xneut(x1, x2, ::Nothing) = 1 - x1 - 2x2
@inline function xneut(x1, x2, xn)
    r = 1 - x1 - 2x2
    r + (xn - ForwardDiff.value(r))
end
"""
Lyα escape rate R_Lyα [s⁻¹] per 2p atom of the three-level model at the iterate: τ_S = 3A_21 λ_α³ n_1s / (8π g), n_1s = n_H x_H with
x_H = 1 − x_H⁺ − 2x_H₂ (the slow species of the iterate; the fast ones are ≤ 1e-6), g the angle-mean velocity gradient
(|∂u/∂r| + 2|u|/r)/3 of the iterate (`sl_grad`).
"""
function lya_escape(ρ, x1, x2, g; xn = nothing)
    n1 = ρ / CNn.μH() * max(xneut(x1, x2, xn), 1e-30)
    τS = 3A21LYA * LYA_LAM^3 * n1 / (8π * max(g, 1e-30))
    3A21LYA * (τS < 1e-6 ? one(τS) : -expm1(-τS) / τS)
end
"""
(SLIMP) the prepared terms of a zone at the star light of the iterate: (q, tc, fs, g) = (W/W_start, T_c/T_c,start, the LW factor
f_sh × (smooth on-weight of τ) of `sl_sweep!`, the velocity gradient g of the zone). Everything that follows from them is recomputed here
exactly: the photo-rates of J = W B(T_c) (R ∝ W; the Milne inverse α = α_spont + W α_stim), kLW = 1.1e8 · 4π W B_ν(T_c) fs (kHm, kH at the
value of T of the iterate), the three-level β = W β_2s(T_c), R_Lyα(ρ, x_H⁺, x_H₂, g), P12 = (R_Lyα + Λ_2s) W/(e^{E21/kT_c} − 1) and ⟨hν⟩(T_c). The
colour temperature, q, fs, g and (ρ, x) are Duals (the photo-rate tables are interpolated in ln T_c, the quadrature of the direct rates carries T_c in its nodes). A zone without a star-light record (no star at the prepare pass) has q = tc = 1 and keeps the prepared W and T_c; its
R_Lyα is still the iterate's. `phc`, `phc0` may be given uncached (NLOIT's path).
"""
function cpl_tq(t, ρ, T, x1, x2, q, tc, fs, g, dp = nothing; phc = t.phc, phc0 = t.slq.phc0, need_cool::Bool = true, xn = nothing)
    s_ = t.slq; p3_0 = s_.p3; extrab0 = t.extra; r = s_.rec
    tenv0 = t.env
    tenv = if CNn.HM88 && hasproperty(tenv0, :hm) && s_.hm0 !== nothing       # HM88: the perturber densities of the iterate (H, H₂, H⁺ from the slow species; the He part of the start)
        nHt = ρ / CNn.μH(); merge(tenv0, (; hm = (nHt * (max(xneut(x1, x2, xn), 0) + s_.hm0[1] + x2), nHt * (x1 + s_.hm0[2]))))
    else
        tenv0
    end
    p3_1 = dp === nothing ? p3_0 : merge(p3_0, (; deep = dp))          # TL3DEEP: the iterate's deep weight
    RL = g === nothing ? p3_1.RL : lya_escape(ρ, x1, x2, g; xn)            # (g = nothing: the prepared R_Lyα, for the opacity's base value)
    if r !== nothing                             # the terms at the iterate's colour temperature (Float) and dilution (q: below)
        Tc = r.Tc0 * tc; W0 = r.W0; sl = STARLIGHT[]
        nb = _div_expm1(W0, E21 / (CNn.kB * Tc))
        p3b = merge(p3_1, (; β = W0 * beta2s_bb(Tc), P12 = (RL + LAM2S) * nb, hν = balmer_mean(Tc), Tc, RL))
        Jb = s_.lj && !hasproperty(tenv, :ljW) ? (Js = tenv.ljs.J; ν -> Js(ν) + W0 * _bnu(ν, Tc)) : (ν -> W0 * _bnu(ν, Tc))          # (J = J_static + W B(T_c))
        envb = s_.lj ? (hasproperty(tenv, :ljW) ? merge(tenv, (; Trad = Tc, J = Jb, ljW = W0, tl3 = p3b)) :
                                                   merge(tenv, (; Trad = Tc, J = Jb, tl3 = p3b))) : merge(tenv, (; tl3 = p3b))
        ph_ = s_.ph::Vector{CNn.Photo}; phc1 = Tv -> _phvec(ph_, Tv, envb)
        if t.extra !== nothing
            e1 = (Trad = Tc, J = Jb); Tv = Float64(ForwardDiff.value(T))
            kHm = (sl.pHm === nothing || (LOCALJ && ljon(s_.p3))) ? 0.0 : CNn.photorate(sl.pHm, Tv, e1)
            kH = sl.pH === nothing ? 0.0 : CNn.photorate(sl.pH, Tv, e1)
            kLW = 1.1e8 * 4π * W0 * _bnu(12.87 * EVERG / CNn.hP, Tc) * fs
            extrab = starlight_terms(s_.mm, r.ρ, kHm, kH, kLW)
        else
            kLW = zero(fs); extrab = extrab0
        end
    else
        p3b = merge(p3_1, (; RL, P12 = iszero(p3_1.P12) ? p3_1.P12 : (RL + LAM2S) * p3_1.P12 / (s_.p3.RL + LAM2S)))
        envb = merge(tenv, (; tl3 = p3b))
        phc1 = phc; kLW = 0.0; extrab = extrab0
    end
    # (a rate that is exactly 0 stays a Float 0: q·0 is a Dual zero, and tl3_heat's E21 P12/(Q12 + P12 + 1e-300) then has 0/0 partials
    #  where Q12 underflows (cold gas, T_c ~ T_CMB); the derivative of q·0 is 0 anyway)
    p3q = merge(p3b, (; β = iszero(_val(p3b.β)) ? _val(p3b.β) : q * p3b.β, P12 = iszero(_val(p3b.P12)) ? _val(p3b.P12) : q * p3b.P12))
    env = merge(envb, (; tl3 = p3q))
    phq = (s_.lj && phc0 !== nothing) ? (Tv -> _phq_combine(phc1(Tv), phc0(Tv), q)) : phc1      # (b: W = 0, the static field)
    ex = extrab === nothing ? nothing : ((x, xn_ = nothing) -> q .* (xn_ === nothing ? extrab(x) : extrab(x, xn_)))
    c0 = s_.c0
    cooling = need_cool ? ((ρ_, T_, d) -> (c0 === nothing ? zero(T_) : c0(ρ_, T_, d)) - tl3_heat(ρ_, T_, d, merge(p3q, (; ne = ρ_ / CNn.μH() * max(d[:Hp], zero(d[:Hp])))))[1]) : nothing
    (; env, phc = phq, cooling, extra = ex, kLW)
end
cpl_g(m::NetGas, k, z; xn = nothing) = _cpl_g(m, k, CPLTERMS[][k], z, xn)
function _cpl_g(m::NetGas, k, t0, z, xn)
    ρ = z[1]; T = z[2]
    t = if MGON && length(z) > 9                     # stage 2c: z = (ρ, T, x1, x2, q_1 … q_G, s, g, dp, t), see cpl_tqg
        G_ = length(z) - 8
        mg_hasrec(t0) ? cpl_tqg(t0, ρ, T, z[3], z[4], ntuple(i -> z[4+i], G_), z[5+G_], z[end], z[6+G_], z[7+G_]; xn) : t0
    else
        (length(z) >= 9 && hasproperty(t0, :slq)) ? cpl_tq(t0, ρ, T, z[3], z[4], z[5], z[6], z[7], z[8], z[9]; xn) : t0
    end
    xnd = xn === nothing ? nothing : xneut(z[3], z[4], xn)         # the neutral hydrogen, direct value, derivative of 1 − x₁ − 2x₂
    x, fs = CNn.ck_eval((z[3], z[4]), ρ, T, m.net; env = t.env, phc = t.phc, extra = t.extra, xn = xnd)
    u, p = CNn.eos(ρ, T, x, m.net; xn = xnd)
    dfull = CNn.full_fast(x, m.net; xn = xnd)
    Λ = t.cooling === nothing ? zero(u) : t.cooling(ρ, T, dfull)
    md = li_model(k); Λl = zero(u)                         # (Λl: the line cooling alone, which the MGROUPS hydro solve emits into the IR group)
    md === nothing || (λl_ = li_cool(md, ρ, T, dfull); Λ += λl_; Λl += λl_)      # H2 line cooling of the iterate (LINEIMP)
    hd0 = HDCHEM[] && k <= length(CPLHD0[]) ? CPLHD0[][k] : nothing
    hd0 === nothing || (λh_ = hd_cool(ρ, T, dfull, HDYD[] * hd_x(ρ, T, dfull, hasproperty(t, :kLW) ? t.kLW : 0.0, hd0)[3]); Λ += λh_; Λl += λh_)      # HD cooling of the iterate (HDCHEM)
    # three-level heats of the iterate (chemistry heat cq, cycle heat yq; n_e self-consistent as in the cooling closure): the hydro solve
    # takes them from the gray field implicitly (TitanResidFast cpQ, cpY) instead of the lagged RAYQ/GASQ bookkeeping
    p3_ = hasproperty(t.env, :tl3) ? t.env.tl3 : nothing
    cq, yq = p3_ === nothing ? (zero(u), zero(u)) : tl3_heat(ρ, T, dfull, merge(p3_, (; ne = ρ / CNn.μH() * max(dfull[:Hp], zero(dfull[:Hp])))))
    if CPLDBGZ[] > 0 && z isa Vector{Float64} && k - Int(TF.getscalar(:ngrs)) + 1 == CPLDBGZ[] && CPLDBGN[] > 0
        CPLDBGN[] -= 1; d = CNn.full(x, m.net); p3 = hasproperty(t.env, :tl3) ? t.env.tl3 : nothing
        nHs = ρ / CNn.μH(); p3s = p3 === nothing ? nothing : merge(p3, (; ne = nHs * max(d[:Hp], 0.0)))
        @printf("CPLDBGZ cpl_g zone %d T %.1f ρ %.3e x(H+,H2) %.3e %.3e x_e %.3e x_H- %.3e: Λ %.4e  tl3(lagged ne) %.4e  tl3(self) %.4e  lya_bb(self) %.4e  S %.3e %.3e  u %.6e\n",
                CPLDBGZ[], T, ρ, z[3], z[4], d[:e], d[:Hm], Λ, p3 === nothing ? 0.0 : tl3_heat(ρ, T, d, p3)[1], p3s === nothing ? 0.0 : tl3_heat(ρ, T, d, p3s)[1],
                p3s === nothing ? 0.0 : lya_bb(ρ, T, d, p3s), fs[1], fs[2], u)
    end
    [p, u, fs[1], fs[2], Λ, cq, yq, Λl]      # (Λl: MGROUPS gives the lines' share of Λ to the IR group, mg_rows!)
end
"(diagnostic, with CPLJCHECK) zone k's local (p, e, S₁, S₂, Λ, cq, yq, Λ_lines): ForwardDiff Jacobian against central differences in (ρ, T, x_H⁺, x_H₂)."
function cpl_gcheck(m::NetGas, k)
    z0 = [TM.v(:dn)[k], TM.v(:tn)[k], TX.XSN[][1][k], TX.XSN[][2][k]]
    J = ForwardDiff.jacobian(z -> cpl_g(m, k, z), z0); names = ("p", "e", "S1", "S2", "Λ", "cq", "yq", "Λl"); vn = ("ρ", "T", "x1", "x2")
    for j in 1:4, rh in (1e-6, 1e-4, 1e-2)
        h = rh * z0[j]; zp = copy(z0); zm = copy(z0); zp[j] += h; zm[j] -= h
        fd = (cpl_g(m, k, zp) .- cpl_g(m, k, zm)) ./ 2h
        @printf("CPLGCHECK zone %d d/d%s (FD step %.0e): %s\n", k - Int(TF.getscalar(:ngrs)) + 1, vn[j], rh,
                join((@sprintf("%s AD %.5e FD %.5e", names[i], J[i, j], fd[i]) for i in 1:8 if max(abs(J[i, j]), abs(fd[i])) > 0), " | "))
    end
    y0 = cpl_g(m, k, z0); d0 = CNn.full(cpl_ckx(m, k, z0), m.net)
    @printf("CPLGCHECK zone %d values: Λ %.6e cq %.6e yq %.6e | x_e %.6e x_H %.6e x_H+ %.6e\n", k - Int(TF.getscalar(:ngrs)) + 1, y0[5], y0[6], y0[7], d0[:e], d0[:H], d0[:Hp])
    nothing
end
cpl_ckx(m, k, z) = (t = CPLTERMS[][k]; CNn.ck_eval([z[3], z[4]], z[1], z[2], m.net; env = t.env, phc = t.phc, extra = t.extra)[1])
cpl_lift(y, J, ρ::Float64, T::Float64, a::Float64, b::Float64) = Tuple(y)
function cpl_lift(y, J, ρ::D, T::D, a::D, b::D) where {D<:ForwardDiff.Dual}
    P1, P2, P3, P4 = ForwardDiff.partials(ρ), ForwardDiff.partials(T), ForwardDiff.partials(a), ForwardDiff.partials(b)
    ntuple(i -> D(y[i], J[i, 1] * P1 + J[i, 2] * P2 + J[i, 3] * P3 + J[i, 4] * P4), 8)
end
cpl_lift(y, J, ρ, T, a, b) = (pr = promote(ρ, T, a, b); cpl_lift(y, J, pr...))
"(p, e, S_H⁺, S_H₂, Λ) of coupled zone k at (ρ, T, x_H⁺, x_H₂), with derivatives (see COUPLED above)."
const CPLPN = Threads.Atomic{Int}(0); const CPLPH = Threads.Atomic{Int}(0); const CPLPT = Threads.Atomic{Int}(0)   # (CPLPROF) misses, hits, ns
"""
Jacobian reuse (CPLJLAG, default "3e-3,0.1"): a zone's local Jacobian is recomputed only when ρ or T (or, SLIMP, the star-light state q, t, fs, g) moved by
more than the first number or x_H⁺, x_H₂ by more than the second (relative, ln) since it was taken; the values are always exact. Float evaluations (the EOS refresh in
update!) take values only. The reuse serves the band Jacobian, the preconditioner (the Newton becomes inexact, the converged state is unchanged): the Jacobian-vector
product of the Krylov solve (`TX.AD_JV`, NewtonKrylov.adjv!) is the exact derivative at the iterate, so there a Jacobian is reused only at the very same point.
`CPLJLAG=0`: a new Jacobian at every new point (also in the band).
"""
const CPLJLAG = let s_ = split(get(ENV, "CPLJLAG", "3e-3,0.1"), ","); s_[1] == "0" ? (0.0, 0.0) : (parse(Float64, s_[1]), parse(Float64, s_[2])) end
@inline cpl_close(vJ, v) = vJ !== nothing && (TX.AD_JV[] ? vJ == v : CPLJLAG[1] > 0 && abs(log(v[1] / vJ[1])) < CPLJLAG[1] && abs(log(v[2] / vJ[2])) < CPLJLAG[1] &&
                           abs(log(v[3] / vJ[3])) < CPLJLAG[2] && abs(log(v[4] / vJ[4])) < CPLJLAG[2])
"A local Jacobian is reused (cpl_close and its star-light variants) only at the same neutral hydrogen xn when the product is exact (TX.AD_JV)."
@inline cpl_closex(f, c, v, xn) = f(c[3], v) && (!TX.AD_JV[] || c[6] == xn)
function cpl_vj(m, k, v, xn = nothing)
    t0 = time_ns(); z0 = [v...]
    res = ForwardDiff.DiffResults.JacobianResult(zeros(8), z0)
    ForwardDiff.jacobian!(res, z -> cpl_g(m, k, z; xn), z0)
    Threads.atomic_add!(CPLPN, 1); Threads.atomic_add!(CPLPT, Int(time_ns() - t0))
    copy(ForwardDiff.DiffResults.value(res)), copy(ForwardDiff.DiffResults.jacobian(res))
end
function cpl_zone(m::NetGas, k, ρ, T, a, b, xn)
    v = (Float64(ForwardDiff.value(ρ)), Float64(ForwardDiff.value(T)), Float64(ForwardDiff.value(a)), Float64(ForwardDiff.value(b)))
    needJ = !(ρ isa Float64 && T isa Float64 && a isa Float64 && b isa Float64)
    c = CPLCACHE[][k]                                         # (v, y(v), v_J, J(v_J), xn): xn, the neutral hydrogen of the point, is part of the value's key
    if c === nothing || c[1] != v || c[5] != xn
        if needJ && !(c !== nothing && cpl_closex(cpl_close, c, v, xn))
            y, J = cpl_vj(m, k, v, xn); c = (v, y, v, J, xn, xn)
        else
            y = cpl_g(m, k, [v...]; xn); c = c === nothing ? (v, y, nothing, nothing, xn, nothing) : (v, y, c[3], c[4], xn, c[6])
            Threads.atomic_add!(CPLPH, 1)
        end
        CPLCACHE[][k] = c
    elseif needJ && !cpl_closex(cpl_close, c, v, xn)
        y, J = cpl_vj(m, k, v, xn); c = (v, y, v, J, xn, xn); CPLCACHE[][k] = c
    else
        Threads.atomic_add!(CPLPH, 1)
    end
    needJ ? cpl_lift(c[2], c[4], ρ, T, a, b) : Tuple(c[2])
end
const CPLCACHEQ = Ref(Any[])
@inline cpl_closeq(vJ, v) = cpl_close(vJ, v) && (TX.AD_JV[] || abs(log(v[5] / vJ[5])) < CPLJLAG[1] && abs(log(v[6] / vJ[6])) < CPLJLAG[1] &&     # (the star-light state moves the rates steeply: the
                            abs(v[7] - vJ[7]) <= CPLJLAG[1] * max(v[7], vJ[7], 1e-3) && abs(log(v[8] / vJ[8])) < CPLJLAG[1] && abs(v[9] - vJ[9]) <= CPLJLAG[1])       #  tolerance of ρ, T, not that of the species)
"(SLIMP) local values and Jacobian at v = (ρ, T, x1, x2, q, t, fs, g, dp): AD in all variables (T_c enters the photo-rate tables as a Dual)."
function cpl_vjq(m, k, v, xn = nothing)
    y, J = cpl_vj(m, k, v, xn)
    if TX.SLIMPDBG[] > 0 && !all(isfinite, J)
        TX.SLIMPDBG[] -= 1; bad = [(i, j) for i in 1:8, j in 1:size(J, 2) if !isfinite(J[i, j])]
        println("SLIMPDBG cpl_vjq zone $k v $v: non-finite J at $bad; y $y; rec ", CPLTERMS[][k].slq.rec)
    end
    y, J
end
function cpl_liftq(y, J, x1::D, x2::D, x3::D, x4::D, x5::D, x6::D, x7::D, x8::D, x9::D) where {D<:ForwardDiff.Dual}
    P1, P2, P3 = ForwardDiff.partials(x1), ForwardDiff.partials(x2), ForwardDiff.partials(x3)
    P4, P5, P6 = ForwardDiff.partials(x4), ForwardDiff.partials(x5), ForwardDiff.partials(x6)
    P7, P8, P9 = ForwardDiff.partials(x7), ForwardDiff.partials(x8), ForwardDiff.partials(x9)
    @inbounds ntuple(i -> D(y[i], J[i, 1] * P1 + J[i, 2] * P2 + J[i, 3] * P3 + J[i, 4] * P4 + J[i, 5] * P5 + J[i, 6] * P6 + J[i, 7] * P7 + J[i, 8] * P8 + J[i, 9] * P9), 8)
end
cpl_liftq(y, J, x1::Float64, x2::Float64, x3::Float64, x4::Float64, x5::Float64, x6::Float64, x7::Float64, x8::Float64, x9::Float64) = Tuple(y)
cpl_liftq(y, J, x1, x2, x3, x4, x5, x6, x7, x8, x9) = (pr = promote(x1, x2, x3, x4, x5, x6, x7, x8, x9); cpl_liftq(y, J, pr...))
function cpl_liftq(y, J, xs...)
    D = promote_type(map(typeof, xs)...); D <: ForwardDiff.Dual || return Tuple(y)
    Ps = map(x -> ForwardDiff.partials(D(x)), xs)
    ntuple(i -> D(y[i], sum(J[i, j] * Ps[j] for j in eachindex(Ps))), 8)
end
"(SLIMP) cpl_zone with the star-light state (q, t, fs, g, dp) of the iterate as local variables (cached per zone on the Float values)."
function cpl_zone_q(m::NetGas, k, ρ, T, a, b, q, tc, fs, g, dp, xn)
    V(x) = Float64(ForwardDiff.value(x))
    v = (V(ρ), V(T), V(a), V(b), V(q), V(tc), V(fs), V(g), V(dp))
    needJ = any(x -> x isa ForwardDiff.Dual, (ρ, T, a, b, q, tc, fs, g, dp))
    length(CPLCACHEQ[]) < k && (CPLCACHEQ[] = vcat(CPLCACHEQ[], Any[nothing for _ in length(CPLCACHEQ[])+1:k]))
    c = CPLCACHEQ[][k]
    if c === nothing || c[1] != v || c[5] != xn
        if needJ && !(c !== nothing && cpl_closex(cpl_closeq, c, v, xn))
            y, J = cpl_vjq(m, k, v, xn); c = (v, y, v, J, xn, xn)
        else
            y = cpl_g(m, k, [v...]; xn); c = c === nothing ? (v, y, nothing, nothing, xn, nothing) : (v, y, c[3], c[4], xn, c[6])
        end
        CPLCACHEQ[][k] = c
    elseif needJ && !cpl_closex(cpl_closeq, c, v, xn)
        y, J = cpl_vjq(m, k, v, xn); c = (v, y, v, J, xn, xn); CPLCACHEQ[][k] = c
    end
    needJ ? cpl_liftq(c[2], c[4], ρ, T, a, b, q, tc, fs, g, dp) : Tuple(c[2])
end
"Per zone (x_H⁺, x_H₂, x_HI) the last `cpl_finish!` wrote into the free vector m.x (whose entries cannot carry x_HI = 1e-9 to more than 1e-7): `cpl_prepare!` takes x_HI from it while m.x is unchanged."
const XNKEEP = Ref(Any[])
function cpl_prepare!(m::NetGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)); nf = length(m.net.free)
    dn, tn = TM.v(:dn), TM.v(:tn); mgr = length(dn); TX.cpl_alloc!(mgr)
    CPLHD0[] = HDCHEM[] ? Any[length(x) > nf ? hd_fracs(x[nf+1], x[nf+2]) : nothing for x in m.x] : Any[]
    length(CPLTERMS[]) == mgr || (CPLTERMS[] = Any[nothing for _ in 1:mgr]); fill!(CPLTERMS[], nothing)
    length(CPLCACHE[]) == mgr || (CPLCACHE[] = Any[nothing for _ in 1:mgr]); fill!(CPLCACHE[], nothing)
    SLIMP && (length(CPLCACHEQ[]) == mgr || (CPLCACHEQ[] = Any[nothing for _ in 1:mgr]); fill!(CPLCACHEQ[], nothing);
             length(SLW0[]) == mgr || (SLW0[] = Any[nothing for _ in 1:mgr]); fill!(SLW0[], nothing))
    fill!(TX.CPLZ[], false)                       # chem_core!'s closing fill_eos! uses m.x
    CPLPREP[] = true; sd0 = SLDIAG[]               # (the prepare pass's star-light setup must not replace the last chemistry's diagnostics)
    try chem_core!(m) finally CPLPREP[] = false; SLDIAG[] = sd0 end
    SLIMP && TX.SLIMPDBG[] > 0 && println("SLIMPDBG prepare jstep $(Int(TF.getscalar(:jstep))): $(count(!isnothing, SLW0[])) zones with a star-light record",
                                          isempty(SLW0[]) || all(isnothing, SLW0[]) ? "" : ", e.g. " * string(SLW0[][findfirst(!isnothing, SLW0[])]))
    sl = CNn.ckinfo(m.net).slow
    for k in 1:mgr
        x = m.x[min(k, length(m.x))]
        for i in 1:2; TX.XSN[][i][k] = max(x[sl[i]], 1e-40); end
        # the neutral hydrogen: the one the last cpl_finish! took from the solved state when m.x still holds its x_H⁺, x_H₂ (the free vector cannot carry it to better than 1e-16 of 1),
        # else 1 − 2x_H₂ − x_H⁺ of the stored values (a restart; at least the 1e-40 of the species floor: the logit z of a zone at x_H⁺ = 1 exactly must stay finite)
        kp = k <= length(XNKEEP[]) ? XNKEEP[][k] : nothing
        TX.XSN[][3][k] = (kp !== nothing && kp[1] == TX.XSN[][1][k] && kp[2] == TX.XSN[][2][k]) ? kp[3] : max(1 - 2TX.XSN[][2][k] - TX.XSN[][1][k], 1e-40)
    end
    for k in a:b; TX.CPLZ[][k] = CPLTERMS[][k] !== nothing; end
    NLOIT && nloc_reset!(mgr)                     # (NLOIT: coupled zones' opacity from the terms just prepared)
    # (PERF512) the zones' values first, on MATTHREADS tasks (independent; each fills its own CPLCACHE entry): the loop below then reads
    # them from the cache, the same numbers as computing them there
    let z0 = TX.ZONE0[]
        try TX.zone_sweep!(a:b, TX.cplz) do k
            try cpl_zone(m, k, dn[k], tn[k], TX.XSN[][1][k], TX.XSN[][2][k], TX.XSN[][3][k]) catch end    # (a failing zone fails again, in order, below)
            nothing
        end
        finally TX.ZONE0[] = z0 end
    end
    for k in a:b
        TX.cplz(k) || continue
        y = cpl_zone(m, k, dn[k], tn[k], TX.XSN[][1][k], TX.XSN[][2][k], TX.XSN[][3][k])
        TX.XSRCO[][1][k] = y[3]; TX.XSRCO[][2][k] = y[4]; TX.XLAMO[][k] = y[5]; MGON && (TX.MGLCO[][k] = y[8])
    end
    nothing
end
"(diagnostic, SLSIZE=z1:z2) per zone of the range, the star-light inputs at the start of the step (prepare pass) and at the solved state: T_src, w_q, T_c, W, f_sh × on-weight, R_Lyα, k_LW."
const SLSIZE = Ref{Union{Nothing,UnitRange{Int}}}(haskey(ENV, "SLSIZE") ? (z = parse.(Int, split(ENV["SLSIZE"], ":")); z[1]:z[2]) : nothing)
function sl_size_report(S, dn)
    a = Int(TF.getscalar(:ngrs)); νLW = 12.87 * EVERG / CNn.hP; kLW(W, Tc, fs) = 1.1e8 * 4π * W * _bnu(νLW, Tc) * fs
    for z in SLSIZE[]
        k = a + z - 1; t = CPLTERMS[][k]
        (k <= length(SLW0[]) && (r = SLW0[][k]) !== nothing && t !== nothing) || continue
        RL0 = t.slq.p3.RL; fs0 = r.lw ? r.fsh : 0.0; RL1 = lya_escape(dn[k], TX.XSN[][1][k], TX.XSN[][2][k], S.g[k]; xn = TX.XSN[][3][k]); Tc1 = r.Tc0 * S.t[k]; W1 = r.W0 * S.q[k]
        @printf("SLSIZE step %d zone %d: T_src %.1f -> %.1f | w_q %.4f -> %.4f | T_c %.1f -> %.1f | W %.4e -> %.4e | fs %.4e -> %.4e | R_Lya %.4e -> %.4e | k_LW %.4e -> %.4e\n",
                Int(TF.getscalar(:jstep)), z, r.Tsrc, S.ts[k], r.wq, S.wq[k], r.Tc0, Tc1, r.W0, W1, fs0, S.fs[k], RL0, RL1, kLW(r.W0, r.Tc0, fs0), kLW(W1, Tc1, S.fs[k]))
    end
end
"After the hydro step (and advect_species!, which moved the HD tail): the coupled zones' free species := CK expansion of the solution."
function cpl_finish!(m::NetGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)); nf = length(m.net.free); dn, tn = TM.v(:dn), TM.v(:tn); sl_ = CNn.ckinfo(m.net).slow
    n = length(dn); S = nothing; mgqs = nothing
    if SLIMP                                       # SLIMP: the fast species at the star light of the solved state
        S = (; q = ones(n), t = ones(n), fs = zeros(n), g = ones(n), ts = zeros(n), wq = zeros(n), dp = zeros(n)); er_, fr_ = TM.v(:ern), TM.v(:frn)
        mgqt = nothing
        if MGON                                    # stage 2c: the groups' fields of the solved state (Float), then the H₂ LW closure at its colour temperature
            mgqs = ones(TX.MGN[] + 1, n); Mq_ = (; egn = TX.MGEN[], qs = mgqs); mgqt = (k, Tc) -> sl_qt_mg(Mq_, k, Tc)
        end
        sl_sweep!(S, a, b, TM.v(:rn), dn, tn, TM.v(:un), er_, fr_, TX.XSN[][1], TX.XSN[][2], TX.cplz; mgqt)
        MGON && SLSTAR[] && mg_h2fill!(S.t, mgqs, zeros(n), zeros(n))
        SLSIZE[] === nothing || sl_size_report(S, dn)
    end
    for k in a:b
        TX.cplz(k) || continue
        t = CPLTERMS[][k]
        if SLIMP && hasproperty(t, :slq)
            if MGON && mg_hasrec(t)
                G_ = TX.MGN[]
                t = cpl_tqg(t, dn[k], tn[k], TX.XSN[][1][k], TX.XSN[][2][k], ntuple(g -> mgqs[g, k], G_), mgqs[G_+1, k], S.t[k], S.g[k], S.dp[k]; xn = TX.XSN[][3][k])
            else
                t = cpl_tq(t, dn[k], tn[k], TX.XSN[][1][k], TX.XSN[][2][k], S.q[k], S.t[k], S.fs[k], S.g[k], S.dp[k]; xn = TX.XSN[][3][k])
            end
        end
        x, _ = CNn.ck_eval((TX.XSN[][1][k], TX.XSN[][2][k]), dn[k], tn[k], m.net; env = t.env, phc = t.phc, extra = t.extra, xn = TX.XSN[][3][k])
        old = m.x[k]; m.x[k] = length(old) > nf ? vcat(x, old[nf+1:end]) : copy(x)
        length(XNKEEP[]) < length(dn) && (XNKEEP[] = Any[nothing for _ in dn])
        XNKEEP[][k] = (m.x[k][sl_[1]], m.x[k][sl_[2]], TX.XSN[][3][k])
    end
    nothing
end
cpl_chem!(m::NetGas) = (cpl_finish!(m); CPLPOST[] = true; try chem!(m) finally CPLPOST[] = false end)
"TitanResidFast.CPL_LOCAL: zone k's species rows alone at its current (ρ, T) over dt (see TitanResidFast.cpl_local_all!); nothing if not converged."
function cpl_local!(m::NetGas, k, dt)
    dn, tn, dvn, do_, dvo, dmdt = TM.v(:dn), TM.v(:tn), TM.v(:dvoln), TM.v(:do), TM.v(:dvolo), TM.v(:dmdt)
    ρ = dn[k]; T = tn[k]; mk = ρ * dvn[k]; mo = do_[k] * dvo[k]
    xo = (TX.XSO[][1][k], TX.XSO[][2][k]); xu = (TX.XSN[][1][k+1], TX.XSN[][2][k+1]); xd = (TX.XSN[][1][k-1], TX.XSN[][2][k-1])
    fo = dmdt[k+1]; fi = dmdt[k]                 # mass flux terms of the row: −ṁ_{k+1} x_b(k+1) + ṁ_k x_b(k); fle = −ṁ ≥ 0: outward
    # unknowns (z, ln x_H₂), z the logit ln(x_H⁺/x_HI) of the hydrogen ionisation (TX.hlogit): the neutral fraction of an ionised zone is resolved
    function R(z)
        x2 = exp(z[2]); x1, xn_ = TX.hlogit(z[1], x2); S_ = cpl_g(m, k, [ρ + zero(x1), T + zero(x1), x1, x2]; xn = ForwardDiff.value(xn_))
        xs = (x1, x2)
        ntuple(i -> (mk * xs[i] - mo * xo[i]) / dt - mk * S_[2+i] - fo * (-fo >= 0 ? xs[i] : xu[i]) + fi * (-fi >= 0 ? xd[i] : xs[i]), 2)
    end
    z = [log(max(TX.XSN[][1][k], 1e-40) / TX.XSN[][3][k]), log(max(TX.XSN[][2][k], 1e-40))]
    for it in 1:12
        r0 = collect(R(z)); J = ForwardDiff.jacobian(zz -> collect(R(zz)), z)
        δ = -(J \ r0); all(isfinite, δ) || return nothing
        δ .= clamp.(δ, -3.0, 3.0); z .+= δ
        if maximum(abs, δ) < 1e-4           # (only Newton's initial guess)
            x2 = exp(z[2]); x1, xn_ = TX.hlogit(z[1], x2)
            return (x1, x2, xn_)
        end
    end
    nothing
end
function chem!(m::NetGas)
    nf = length(m.net.free)
    if !HDCHEM[]
        any(x -> length(x) != nf, m.x) && hd_strip!(m)             # (restart from a run with HDCHEM=1)
        return chem_core!(m)
    end
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    dn, tn = TM.v(:dn), TM.v(:tn)
    hd_pad!(m)
    pas = [(x[nf+1], x[nf+2]) for x in m.x]
    hd_strip!(m)                                                   # the network step sees its own free species only
    HDX[] = [max(p[2], 0.0) for p in pas]
    dT = chem_core!(m)
    yD = HDYD[]; lam = zeros(length(m.x)); sd = STARLIGHT[] === nothing ? nothing : SLDIAG[]
    for k in a:b
        d = CNn.full(m.x[k], m.net); nH = dn[k] / CNn.μH()
        kLW = (sd === nothing || k > length(sd.kLW)) ? 0.0 : sd.kLW[k]
        p0 = hd_fracs(pas[k]...)
        p = hd_step(hd_rates(tn[k], nH, d, kLW), p0, m.eqzone[k] ? 0.0 : (dt > 0 ? 1 / dt : Inf))
        m.x[k] = vcat(m.x[k], yD * p[2], yD * p[3])
        lam[k] = Float64(hd_cool(dn[k], tn[k], d, yD * p[3]))
    end
    for k in 1:a-1; m.x[k] = copy(m.x[a]); end
    for k in b+1:length(m.x); m.x[k] = copy(m.x[b]); end
    HDLAM[] = lam
    dT
end

end # module
