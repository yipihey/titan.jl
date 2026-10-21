"""
    StarMaterial

S4: `StarGas`, the TITAN material for a star whose composition is evolved: per-zone NucNet abundance vectors Y[i,k]
(19 species, mol/g; `NucNet.SP`), burning, and convective mixing.

* Composition seen by the physics (`xyz`): X = Σ_{p,d} A Y (hydrogen), Y_He = Σ_{³He,⁴He,⁷Li,⁷Be,⁸B} A Y, Z = Σ_{rest} A Y
  (CNO and F; split C:N:O by `zmix` inside `StellarEOS`). ΣA Y = 1 holds to 1e-15 in the network. Opacity: StellarOpacity κ_R(ρ,T,X)
  (Y = 1−X, Z_CNO = 0 in the tables) incl. conduction, κ_P; `opac = :z0` restores the Z0Physics table of `StellarGas`.
* Gas only (E_r is TITAN's). `has_source`/`material_source(m, ρ, T)`: ε_nuc(ρ, T; Y_lagged) net of neutrinos, generic in the number
  type (ForwardDiff through `NucNet.rates_k`; the screening plasma uses Float64 ρ and the lagged Y, so the ρ-dependence of the screening
  enhancement, a few % of ∂ε/∂ρ, is not in the Jacobian). The Y of a zone is fixed during the Newton iteration of a step.
* `burn!(m)` after every converged step (as `PrimordialGas`' `chem!`): per zone `NucNet.burn` at fixed ρ and the temperature of the
  step (`burnT = :end`: T_n; `:theta`: θT_n + (1−θ)T_o) over Δt → new Y. Energy: the hydro step has already put
  Einj = Δt (θ ε(ρ_n,T_n) + (1−θ) ε(ρ_o,T_o)) (ρΔV)_θ into the zone (read back from the residual: `TX.step_diag`);
  the network released Ebur = ε_burn Δt (time average over the substeps, = Σ_r N_A (Q_r − ⟨E_ν⟩_r) nev_r exactly). The difference
  Ebur − Einj (per gram of the zone) is added to the gas energy e at fixed ρ, and T is re-solved from the EOS with the NEW Y
  (Newton, 1e-13): total energy (gas + cumulative nuclear release) is conserved to the Newton tolerance; `ledger` records Einj, Ebur
  (`dabs/Eabs` = Σ|Ebur − Einj|/Σ|Ebur| is the size of the operator-splitting correction, `resid` the energy error of the T re-solve).
* `mix!(m)` after `burn!`: implicit conservative diffusion of all species on the Lagrangian mass grid,
  ∂Y/∂t = ∂/∂m[(4πr²ρ)² D ∂Y/∂m], D = v ℓ/3 at interfaces with F_conv > 0 and MLT weight > 0, where
  v = ℓ (gδ/(8H_P))^½ (∇−∇')^½ (Kippenhahn–Weigert eq. for the mean convective velocity; ζ² = ∇−∇' of the MLT block, same
  ℓ = αH_P, same g, δ, ∇ as F_conv, so v ℓ²... consistent with F_conv = ρ c_P T (gδ)^½ ℓ² ζ³/(4√2 H_P^{3/2})). Tridiagonal, one
  factorisation for all species. Σ_k dm_k Y_k is conserved to round-off. Composition change at fixed e: T is re-solved at fixed
  ρ and the zone's gas energy e, so mixing changes no zone's energy (energy exactly conserved; nuclear binding energy moves with
  the material).
* S5 (consistent coupling, TR-BDF2): `install_stage_hook!` books, at every converged implicit solve (θ-step, TR stage, BDF stage, BE fallback), the
  energy injected by the residual's source and lost through the boundary with the weights of the scheme; `burn!` integrates the network along
  the T/ρ history of the step (`burnT = :path`; E_inj of every stage is known, so TR-BDF2 needs no special case); `solve_step!` iterates a per-zone
  level factor f_k on ε in the residual (state restored from a snapshot between iterations) until Σ|E_bur − E_inj| ≲ 3e-3 ΣE_bur and no zone's mismatch exceeds `ztol` of its own energy (default 2e-4, the star runs 3e-3, i.e. ΔT/T ≲ 0.2 %; the central zones release 1e-4 of the total, a sum criterion alone left 5 % T_c excursions after the post-step correction). Without a stage
  hook `burn!` behaves as in S4 (E_inj from `step_diag`, constant T_n, f = 1), which the S4 tests rely on.
* Time-step control: after burn!, `m.dtcap` = Δt min(dX_max/max_k|ΔX_H|, f_e e/max_{i,k}|ΔY_i|Q_i N_A MeV, ...), consulted by
  `TitanCore.timstp!` through `TitanCore.DT_HOOKS` and also applied at once to `dtime` (timstp! runs before burn!).
"""
module StarMaterial

import ..TitanDriver, ..TitanCore, ..TitanResidFast, ..TitanMaterials, ..Z0Physics, ..StellarEOS, ..StellarOpacity, ..NucNet
using ForwardDiff
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const TM = TitanMaterials; const TC = TitanCore
const SE = StellarEOS; const Z0 = Z0Physics; const SO = StellarOpacity; const NN = NucNet
const NS = NN.NS; const MeV = NN.MeV; const NA = NN.NA

const IH = NN.IDX[:p]; const ID = NN.IDX[:d]
const HEL = [NN.IDX[s] for s in (:he3, :he4, :li7, :be7, :b8)]
const METALS = [i for i in 1:NS if !(i == IH || i == ID || i in HEL)]
# Energy-weighted time-step criterion: Qspec[i] = MeV released when species i is consumed by its main channel (p: 4p→⁴He per proton,
# D: d(p,γ), ³He: ³He(³He,2p), ⁴He: 3α per α, ⁷Li: 7Li(p,α)). Catalysts and short-lived intermediates (CNO, ⁷Be, ⁸B, ...) have Q = 0:
# they sit in quasi-equilibrium within the implicit burn step and their (large) binding energies are not released by changing them.
const Qspec = let q = zeros(NN.NS)
    q[NN.IDX[:p]] = 6.68; q[NN.IDX[:d]] = 5.494; q[NN.IDX[:he3]] = 12.86; q[NN.IDX[:he4]] = 7.275 / 3; q[NN.IDX[:li7]] = 17.35; q
end
# binding energy per nucleus [MeV] from the AME2020 atomic mass excesses (keV): B = Z Δ(¹H) + (A−Z) Δ(n) − Δ
const BIND = [(NN.ZZ[i] * NN.MEXC[:p] + (NN.AA[i] - NN.ZZ[i]) * NN.MEXC[:n] - NN.MEX[i]) * 1e-3 for i in 1:NS]

"Energy booked at every converged implicit solve of the step (TitanDriver.STAGE_HOOK), with the weights of the scheme."
mutable struct StageAcc
    Einj::Vector{Float64}                 # per zone: erg injected by the residual's nuclear source (4π included)
    Fout::Float64; Fin::Float64           # energy through the outer / inner boundary, erg
    Tg::Vector{Float64}; dg::Vector{Float64}   # zone T, ρ at the end of the TR stage
    last::Symbol                          # kind of the last stage booked (:single, :tr, :bdf, :be)
    ngam::Int
    nstage::Int
end
StageAcc(n) = StageAcc(zeros(n + 2), 0.0, 0.0, zeros(n + 2), zeros(n + 2), :none, 0, 0)

mutable struct StarGas{T} <: TM.Material
    Y::Matrix{Float64}                    # Y[i, k] mol/g
    zmix::NamedTuple{(:C, :N, :O),NTuple{3,Float64}}
    tab::T
    opac::Symbol                          # :stellar (StellarOpacity) | :z0 (Z0Physics table, as StellarGas)
    grey::Bool
    coulomb::Bool; pairs::Bool
    kcap::Float64                         # soft upper limit of κ_R/κ_P (Γ cap, S5 hand-over; Inf = off)
    nuc::Bool                             # nuclear source in the residual and burn! on
    screen::Symbol
    burnT::Symbol
    dXmax::Float64; efrac::Float64
    dtcap::Float64
    capE::Float64                         # energy part of the post-step cap (burn!)
    Yprev::Matrix{Float64}                # abundances at the start of burn! (for the net, mixing-aware ΔX_H)
    mixing::Bool
    ledger::NamedTuple{(:Einj, :Ebur, :dabs, :Eabs, :resid, :nsteps),NTuple{6,Float64}}   # cumulative totals (erg, 4π incl.); corr = dabs/Eabs
    last::NamedTuple
    Dmix::Vector{Float64}                 # diffusion coefficient at interfaces of the last mix!
    lambda_mix::Float64                   # fraction of the composition change mixed in the last mix! (1 unless tmax limited it)
    fscale::Vector{Float64}               # S5: last measured E_bur/E_inj⁰ per zone (diagnostic)
    fslope::Vector{Float64}               # S5: d(E_bur/E_inj⁰)/dΔ [1/s]: predictor of the level factor, f = 1 + fslope·Δ, see solve_step!
    fcur::Vector{Float64}                 # S5: level factor on ε in the residual of the step being solved (iterated by solve_step!)
    acc::StageAcc                         # S5: energy booked by the stage hook during the step
    fcoupling::Bool                       # use fscale (true) or the plain ε of S4 (false)
end
function StarGas(n::Integer; X = 0.76, Y = 0.24, Z = 0.0, zmix = SE.ZMIX, tab = Z0.table(), opac = :stellar, grey = false,
                 coulomb = true, pairs = true, kcap = Inf, nuc = false, screen = :chugunov, burnT = :path, dXmax = 1e-3, efrac = 0.05, mixing = true, fcoupling = true)
    Yv = zeros(NS, n)
    s = zmix.C + zmix.N + zmix.O
    for k in 1:n
        Yv[IH, k] = X; Yv[NN.IDX[:he4], k] = Y / 4
        Yv[NN.IDX[:c12], k] = Z * zmix.C / s / 12; Yv[NN.IDX[:n14], k] = Z * zmix.N / s / 14; Yv[NN.IDX[:o16], k] = Z * zmix.O / s / 16
    end
    StarGas(Yv, zmix, tab, opac, grey, coulomb, pairs, kcap, nuc, screen, burnT, dXmax, efrac, Inf, Inf, copy(Yv), mixing,
            (Einj = 0.0, Ebur = 0.0, dabs = 0.0, Eabs = 0.0, resid = 0.0, nsteps = 0.0), (;), zeros(n + 2), 1.0, ones(n + 2), zeros(n + 2), ones(n + 2), StageAcc(n), fcoupling)
end

"Metal mass fractions below this are set to zero in the EOS/opacity composition (their effect is far below 1e-18)."
const ZCUT = 1e-18
"(X, Y_He, Z) mass fractions of an abundance vector."
@inline function xyz(Y)
    X = NN.AA[IH] * Y[IH] + NN.AA[ID] * Y[ID]
    Yh = 0.0; for i in HEL; Yh += NN.AA[i] * Y[i]; end
    Z = 0.0; for i in METALS; Z += NN.AA[i] * Y[i]; end
    Z < ZCUT && (Z = 0.0)             # trace metals (3α makes 1e-240 of C): the Saha/EOS code gives NaN derivatives for such values
    X, Yh, Z
end
@inline zonecomp(m::StarGas, k) = (t = xyz(view(m.Y, :, k)); SE.composition(t[1], t[2], t[3], m.zmix))
@inline zonecomp(m::StarGas) = zonecomp(m, TX.ZONE[])
@inline hmass(m::StarGas, k) = NN.AA[IH] * m.Y[IH, k] + NN.AA[ID] * m.Y[ID, k]

"Set the abundances of zone k (a vector of 19 mol/g)."
setY!(m::StarGas, k, Y) = (m.Y[:, k] .= Y; m)
function ghosts!(m::StarGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    for k in max(a - 2, 1):a-1; m.Y[:, k] .= m.Y[:, a]; end
    for k in b+1:min(b + 2, size(m.Y, 2)); m.Y[:, k] .= m.Y[:, b]; end
end

"""
S6a: `TFLOOR[]` (default 0 = off): the gas EOS/opacity/thermodynamics are evaluated at T for T ≥ TFLOOR and at TFLOOR + 0.05 (T − TFLOOR) below. The Saha branch of `StellarEOS` and the opacity tables end near 3e3 K;
the optically thin top layers of an accreting protostar (gas decoupled from the 4600 K radiation field: κ_P(T_gas) → 0) cooled below it and Newton fell off the tables (NaN).
"""
const TFLOOR = Ref(0.0)
@inline tfl(T) = TFLOOR[] > 0 && T < TFLOOR[] ? TFLOOR[] + 0.05 * (T - TFLOOR[]) : T       # continuous, slope 0.05 below the floor (a flat floor makes the energy row singular in T)
@inline function TX.material_gas(m::StarGas, ρ, T)
    T = tfl(T)
    g = SE.gas_state(ρ, T, zonecomp(m); pairs = m.pairs, coulomb = m.coulomb)
    return g.p, g.e, sqrt(5.0 / 3.0 * g.p / ρ)
end
@inline function TX.material_gasne(m::StarGas, ρ, T)
    T = tfl(T)
    g = SE.gas_state(ρ, T, zonecomp(m); pairs = m.pairs, coulomb = m.coulomb)
    return g.p, g.e, sqrt(5.0 / 3.0 * g.p / ρ), g.ne
end
"C¹ soft cap: identity below 0.7c, saturating to c above."
@inline softcap(κ, c) = isinf(c) ? κ : (κ <= 0.7c ? κ : c - 0.3c * exp(-(κ - 0.7c) / (0.3c)))
"S6a: opacities are evaluated at max(T, TFLOOR_OPAC[]) (default 0 = off): the tables end at 3e3 K and extrapolate to garbage below (cold gas in the thin atmosphere of an accreting protostar)."
const TFLOOR_OPAC = Ref(0.0)
@inline function TX.material_opac(m::StarGas, ρ, T)
    T = tfl(T)
    TFLOOR_OPAC[] > 0 && T < TFLOOR_OPAC[] && (T = TFLOOR_OPAC[] + 0 * T)
    if m.opac === :z0
        κR = Z0.kap(m.tab, ρ, T)
        m.grey && return κR, κR, κR
        κP = Z0.kapP(m.tab, ρ, T)
        return κR, κP, κP
    end
    X = hmass(m, TX.ZONE[])
    κR = softcap(SO.kappa_R(ρ, T, X), m.kcap)
    m.grey && return κR, κR, κR
    κP = softcap(SO.kappa_P(ρ, T, X), m.kcap)
    return κR, κP, κP
end
function TX.material_thermo(m::StarGas, ρ, T)
    T = tfl(T)
    st = SE.eos(ρ, T; comp = zonecomp(m), pairs = m.pairs, coulomb = m.coulomb, radiation = true)
    return st.∇ad, st.cP, st.δ
end

# ----------------------------------------------------------------------------------------------- nuclear source
"""
    eps_nuc(ρ, T, Y; screen) -> ε [erg g⁻¹ s⁻¹]

Net of neutrino losses, generic in the number type of ρ, T (screening plasma from Float64 ρ and the abundances Y).
"""
function eps_nuc(ρ, T, Y; screen = :chugunov)
    pl = NN.Plasma(Float64(ForwardDiff.value(ρ)), Y)
    k = NN.rates_k(T, ρ, pl, NN.ye_of(Y); screen)
    e = zero(T) * 1.0
    for (ir, r) in enumerate(NN.NET)
        e += (NA * MeV * (r.Q - r.Enu)) * NN.event(r, k[ir], Y)
    end
    e
end
TX.has_source(m::StarGas) = m.nuc
TX.material_source(m::StarGas, ρ, T) = m.nuc ? (m.fcoupling ? m.fcur[TX.ZONE[]] : 1.0) * eps_nuc(ρ, T, view(m.Y, :, TX.ZONE[]); screen = m.screen) : zero(T)

# ----------------------------------------------------------------------------------------------- state fills
function TM.fill_eos!(m::StarGas)
    dn, tn = TM.v(:dn), TM.v(:tn)
    pgn, egn, asn, xnen = TM.v(:pgn), TM.v(:egn), TM.v(:asn), TM.v(:xnen)
    @inbounds for k in TM.krange()
        ρ, T = dn[k], tn[k]
        (ρ > 0 && T > 0) || continue
        TX.ZONE[] = k
        pg, eg, as_, ne = TX.material_gasne(m, ρ, T)
        pgn[k] = pg; egn[k] = eg; asn[k] = as_; xnen[k] = ne
    end
    nothing
end
function TM.fill_opac!(m::StarGas)
    dn, tn = TM.v(:dn), TM.v(:tn)
    chifn, xken, xkpn, xnen = TM.v(:chifn), TM.v(:xken), TM.v(:xkpn), TM.v(:xnen)
    @inbounds for k in TM.krange()
        (dn[k] > 0 && tn[k] > 0) || continue
        TX.ZONE[] = k
        χ, κE, κP = TX.material_opac(m, dn[k], tn[k])
        chifn[k] = χ; xken[k] = κE; xkpn[k] = κP
        xnen[k] = TX.material_gasne(m, dn[k], tn[k])[4]
    end
    nothing
end

"Refresh ghost zones, EOS, opacity, E_r + e_g sums after the state arrays were changed by burn!/mix!."
refresh!() = TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())

# ----------------------------------------------------------------------------------------------- T from e
"T such that e_gas(ρ, T; comp) = e (Newton in ln T with c_V, safeguarded). Returns (T, relative error of e)."
function T_of_e(ρ, e, T0, comp; pairs = true, coulomb = true)
    TFLOOR[] > 0 && T0 < TFLOOR[] && (T0 = TFLOOR[] + 0.05 * (T0 - TFLOOR[]))      # (TFLOOR: the EOS argument is T_eff = tfl(T); solved for T_eff, mapped back below)
    T, err = _T_of_e(ρ, e, T0, comp; pairs, coulomb)
    TFLOOR[] > 0 && T < TFLOOR[] && (T = TFLOOR[] + (T - TFLOOR[]) / 0.05)
    T, err
end
function _T_of_e(ρ, e, T0, comp; pairs = true, coulomb = true)
    T = T0; lo = 0.0; hi = Inf
    err = Inf
    for it in 1:40
        st = SE.eos(ρ, T; comp, pairs, coulomb)
        f = st.e - e; err = abs(f) / abs(e)
        err < 1e-12 && return T, err
        f > 0 ? (hi = T) : (lo = T)
        Tn = T - f / st.cV
        Tn = Tn <= lo || Tn >= hi ? (isfinite(hi) && lo > 0 ? sqrt(lo * hi) : (isfinite(hi) ? 0.5 * hi : 2T)) : Tn
        Tn = clamp(Tn, 0.3T, 3T)
        T = Tn
    end
    st = SE.eos(ρ, T; comp, pairs, coulomb)
    T, abs(st.e - e) / abs(e)
end

# ----------------------------------------------------------------------------------------------- stage hook
const γTR = 2 - sqrt(2.0); const ATR = 1 / (γTR * (2 - γTR))      # TR-BDF2: y_{n+1} = A y_γ − B y_n + (1−γ)/(2−γ) Δ f(y_{n+1})
"""
    install_stage_hook!(m)

Registers `TitanDriver.STAGE_HOOK`: after every converged implicit solve (θ-step `:single`, TR-BDF2 stage 1 `:tr`, stage 2 `:bdf`,
backward-Euler fallback `:be`) one Float64 residual evaluation books, per zone, the energy injected by the nuclear source
(srcw·dt_stage) and the energy that left through the outer boundary (fout), weighted as the scheme combines the stages
(y_{n+1} = A y_γ − B y_n + C₂Δ f(y_{n+1}) ⇒ weight A on stage 1, 1 on stage 2; A γ + C₂ = 1). The state T, ρ at the end of
stage 1 is kept as the middle knot of the temperature history that `burn!` integrates.
"""
function install_stage_hook!(m::StarGas)
    TD.STAGE_HOOK[] = function (kind, dt)
        a = m.acc
        dn, tn = TM.v(:dn), TM.v(:tn)
        sd = TX.stage_diag(dt)
        w = 1.0
        if kind === :tr || kind === :single || kind === :be
            fill!(a.Einj, 0.0); a.Fout = 0.0; a.Fin = 0.0; a.nstage = 0
        end
        kind === :tr && (w = ATR)
        for k in eachindex(a.Einj); k <= length(sd.srcw) && (a.Einj[k] += w * 4π * sd.srcw[k] * dt); end
        a.Fout += w * 4π * sd.fout; a.Fin += w * 4π * sd.fin
        if kind === :tr; copyto!(a.Tg, tn); copyto!(a.dg, dn); end
        a.last = kind; a.nstage += 1
        nothing
    end
    nothing
end

# ----------------------------------------------------------------------------------------------- burn!
"ln-quadratic (3 knots at s = 0, γ, 1) or ln-linear (2 knots at 0, 1) interpolation."
@inline function pathval(s, v0, vg, v1, tr::Bool)
    if tr
        γ = γTR
        return v0 * (s - γ) * (s - 1) / γ + vg * s * (s - 1) / (γ * (γ - 1)) + v1 * s * (s - γ) / (1 - γ)
    end
    (1 - s) * v0 + s * v1
end

"""
    burn_zone(m, k, dt, tr, legacy) -> (Y_new, e_released [erg/g])

Network integration of zone k over the step along the T/ρ history booked by the stage hook (see `burn!`). Does not change `m`.
"""
function burn_zone(m::StarGas, k, dt, tr::Bool, legacy::Bool)
    dn, tn, to, do_ = TM.v(:dn), TM.v(:tn), TM.v(:to), TM.v(:do)
    acc = m.acc; ρ = dn[k]; Y0 = m.Y[:, k]
    if m.burnT === :end || legacy
        res = NN.burn(Y0, ρ, tn[k], dt; screen = m.screen); return res.Y, res.ε * dt
    end
    l0, l1 = log(to[k]), log(tn[k]); lg = tr ? log(acc.Tg[k]) : 0.0
    r0, r1 = log(do_[k]), log(ρ); rg = tr ? log(acc.dg[k]) : 0.0
    dl = max(abs(l1 - l0), tr ? abs(lg - l0) : 0.0)
    ns = tn[k] < 2e6 ? 1 : clamp(ceil(Int, 16dl / 0.15), 1, 6)
    Yk = Y0; eb = 0.0
    for j in 1:ns
        sj = (j - 0.5) / ns
        Tj = exp(pathval(sj, l0, lg, l1, tr)); ρj = exp(pathval(sj, r0, rg, r1, tr))
        res = NN.burn(Yk, ρj, Tj, dt / ns; screen = m.screen)
        eb += res.ε * dt / ns; Yk = res.Y
    end
    Yk, eb
end

"""
    burn!(m::StarGas; dt = timen − timeo) -> NamedTuple

After a converged step. Per zone the network is integrated over the step along the temperature/density history of the step
(`burnT = :path`: ln T and ln ρ interpolated through the knots (t_n, T_o), [(t_γ, T_γ) for TR-BDF2], (t_{n+1}, T_n), in `nseg` equal
segments with `nseg` ≈ 16 Δln T / 0.15, 1..6; `:end`: constant T_n, S4 behaviour). The energy booked by the stage hook as injected
by the residual E_inj is compared with the released E_bur; their difference is added to the zone's gas energy at fixed ρ (total
energy exactly conserved). The ratio E_bur/E_inj⁰ (E_inj⁰ = E_inj/f used, f = 1 + fslope Δ) gives the new slope fslope = (E_bur/E_inj⁰ − 1)/Δ, so the residual of the next
step (any Δ) injects ε·(1 + fslope Δ) (`fcoupling`): the implicit T-dependence of the residual's heating is kept, its level follows the network
(the dominant mismatch, the growth of ε through the self-produced CNO and the T history, is linear in Δ).
"""
function burn!(m::StarGas; dt = TF.getscalar(:timen) - TF.getscalar(:timeo), sd = nothing)
    m.nuc || return (; Einj = 0.0, Ebur = 0.0, corr = 0.0, resid = 0.0, dXH = 0.0, dtcap = Inf)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dn, tn, egn, to, do_ = TM.v(:dn), TM.v(:tn), TM.v(:egn), TM.v(:to), TM.v(:do)
    dvoln = TM.v(:dvoln)
    acc = m.acc
    m.Yprev .= m.Y
    tr = acc.last === :bdf
    legacy = acc.nstage == 0                      # no stage hook booked this step: S4 behaviour (E_inj from step_diag, constant T_n, f = 1)
    if legacy
        sd === nothing && (sd = TX.step_diag(dt))
        fill!(acc.Einj, 0.0); for k in eachindex(sd.srcw); acc.Einj[k] = 4π * sd.srcw[k] * dt; end
        tr = false
    end
    Einj = 0.0; Ebur = 0.0; Eabs = 0.0; dcorr = 0.0; resid = 0.0; dXH = 0.0; dEfrac = 0.0; dEsp = 0.0
    fnew = ones(length(m.fscale)); slope_new = zeros(length(m.fscale))
    for k in a:b
        ρ = dn[k]; mk = 4π * dn[k] * dvoln[k]
        Y0 = m.Y[:, k]
        Yk, eb = burn_zone(m, k, dt, tr, legacy)
        ei = acc.Einj[k]                         # erg injected in the zone by the residual (stage weights included)
        Einj += ei; Ebur += eb * mk; Eabs += abs(eb * mk)
        dcorr += abs(eb * mk - ei)
        # next step's level factor: network release over this step / what the residual would have injected with f = 1
        ei0 = ei / m.fcur[k]                         # what the residual would have injected with the level factor 1
        fnew[k] = (m.fcoupling && !legacy && ei0 > 0 && eb * mk > 1e-6 * ei0) ? clamp(eb * mk / ei0, FLIM[1], FLIM[2]) : 1.0
        slope_new[k] = fnew[k] == 1.0 ? 0.0 : (fnew[k] - 1.0) / dt        # predictor for the next step's f (solve_step!)
        et = egn[k] + (eb * mk - ei) / mk
        m.Y[:, k] .= Yk
        comp = zonecomp(m, k)
        Tn, er = T_of_e(ρ, et, tn[k], comp; pairs = m.pairs, coulomb = m.coulomb)
        resid = max(resid, er); tn[k] = Tn
        dXH = max(dXH, abs(NN.AA[IH] * (Yk[IH] - Y0[IH]) + NN.AA[ID] * (Yk[ID] - Y0[ID])))
        for i in 1:NS
            dEsp = max(dEsp, abs(Yk[i] - Y0[i]) * Qspec[i] * NA * MeV / abs(et))
        end
        dEfrac = max(dEfrac, abs(eb) / abs(et))
    end
    m.fscale .= fnew; m.fslope .= slope_new; acc.nstage = 0
    ghosts!(m)
    tn[a-1] = tn[a]; tn[b+1] = tn[b]; tn[b+2] = tn[b]
    refresh!()
    cap = Inf
    dXH > 0 && !REGION_CAP[] && (cap = min(cap, dt * m.dXmax / dXH))
    dE = max(dEsp, dEfrac)
    capE = dE > 0 ? dt * m.efrac / dE : Inf
    m.capE = capE; cap = min(cap, capE)
    m.dtcap = cap
    TF.setscalar!(:dtime, min(TF.getscalar(:dtime), cap))
    corr = Ebur == 0 ? 0.0 : dcorr / abs(Ebur)
    L = m.ledger
    m.ledger = (Einj = L.Einj + Einj, Ebur = L.Ebur + Ebur, dabs = L.dabs + dcorr, Eabs = L.Eabs + Eabs, resid = max(L.resid, resid), nsteps = L.nsteps + 1)
    m.last = (; Einj, Ebur, corr, resid, dXH, dtcap = cap, dEsp, dEfrac, Fout = acc.Fout, kind = acc.last)
end

"""
    install_state_predictor!(; maxratio = 1.5)

Newton's initial guess of a step (TitanCore.PEQ_HOOK): T and E_r of every zone are extrapolated in time, ln y(t_n + Δ) = ln y_n + (Δ/Δ_prev)
(ln y_n − ln y_{n−1}) (ratio capped at `maxratio`; the grid is extrapolated by `peq!` as before). TITAN starts from the old state,
which in a star evolving on t_KH ≫ t_dyn is off by the change of the last step in every zone; for a convective core the error of ∇ − ∇_ad
(∼ 1e-3 against the solution's 1e-7) puts Newton on the (∇ − ∇_ad)^{3/2} law of the MLT flux far from its root, where it converges only
linearly, by the factor 1 − 1/p = ⅓ per iteration (8–14 iterations per step, and dt limited by the iteration count).
"""
function install_state_predictor!(; maxratio = 1.5)
    hist = Ref{Any}(nothing)
    TC.PEQ_HOOK[] = function ()
        a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
        tn, ern = TM.v(:tn), TM.v(:ern)
        t_old = TF.getscalar(:timeo); dt = TF.getscalar(:dtime)
        cur = (; t = t_old, T = copy(tn), E = copy(ern))
        h = hist[]
        if h !== nothing && t_old > h.t
            r = clamp(dt / (t_old - h.t), 0.0, maxratio)
            for k in a:b
                (h.T[k] > 0 && h.E[k] > 0 && tn[k] > 0 && ern[k] > 0) || continue
                tn[k] *= clamp((cur.T[k] / h.T[k])^r, 0.8, 1.25); ern[k] *= clamp((cur.E[k] / h.E[k])^r, 0.5, 2.0)
            end
        end
        hist[] = cur
        nothing
    end
    nothing
end

"""
    solve_step!(m, stepfn; tol = 2e-3, maxit = 4) -> (mismatch, iterations)

One implicit step with a nuclear source that is consistent with the network (`stepfn()` = `TitanDriver.STEPFN[]`, which starts from
`peq`). The residual's source is ε(ρ,T;Y_lag)·f_k with a per-zone level factor f_k = `fcur`. Starting from the predictor
f_k = 1 + fslope_k Δ, the step is solved, the network is integrated along the step's T/ρ history (`burn_zone`, no commit) and
f_k ← E_bur,k/E_inj⁰,k (E_inj⁰ = E_inj/f used) until Σ|E_bur − E_inj| ≤ tol Σ E_bur; between iterations the state is restored from a
snapshot. The implicit T-dependence of the heating stays in the Jacobian, so the thermostat of the burning core is implicit, and what
`burn!` has to correct afterwards (E_bur − E_inj) is ≲ tol of the energy released: no operator-splitting heat pulse at the end of a step.
"""
const VERBOSE_SPLIT = Ref(false)
"Range of the level factor f_k on ε in the residual: zones where ³He, D or the CNO catalysts relax within a step (cool edge of the core) release far less than ε(Y_lag)·Δ."
const FLIM = (1e-3, 1e2)
function solve_step!(m::StarGas, stepfn; tol = 2e-3, maxit = 4, ztol = 2e-4)
    if !(m.nuc && m.fcoupling); stepfn(); return (NaN, 1, NaN); end
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dt_plan = TF.getscalar(:dtime)
    m.fcur .= clamp.(1 .+ m.fslope .* dt_plan, FLIM[1], FLIM[2])
    snap = TD.snapshot()
    mism = Inf; zmis = Inf; it = 0
    for outer it in 1:maxit
        stepfn()
        dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
        tn, dn, dvoln, egn = TM.v(:tn), TM.v(:dn), TM.v(:dvoln), TM.v(:egn); acc = m.acc
        tr = acc.last === :bdf
        sumb = 0.0; sumd = 0.0; zmis = 0.0; fnewv = copy(m.fcur)
        for k in a:b
            (tn[k] < 1e6 || acc.Einj[k] <= 0) && (fnewv[k] = 1.0; continue)
            Yk, eb = burn_zone(m, k, dt, tr, false)
            mk = 4π * dn[k] * dvoln[k]; ebk = eb * mk; sumb += ebk; sumd += abs(ebk - acc.Einj[k])
            zmis = max(zmis, abs(ebk - acc.Einj[k]) / (mk * abs(egn[k])))      # the heat the post-step correction would add, relative to the zone's own energy
            ei0 = acc.Einj[k] / m.fcur[k]
            fnewv[k] = (ebk > 1e-6 * ei0) ? clamp(ebk / ei0, FLIM[1], FLIM[2]) : 1.0
            VERBOSE_SPLIT[] && abs(ebk - acc.Einj[k]) > 1e-3 * 1e47 && println("         zone ", k - a + 1, " Ebur ", ebk, " Einj ", acc.Einj[k], " f ", m.fcur[k], " -> ", fnewv[k], " T ", tn[k])
        end
        mism = sumb > 0 ? sumd / sumb : 0.0
        VERBOSE_SPLIT[] && println("      split it ", it, " mism ", round(mism; sigdigits = 3), " zone-wise ", round(zmis; sigdigits = 3))
        (mism <= tol && zmis <= ztol || it == maxit) && break
        TD.restore!(snap); m.fcur .= fnewv
    end
    (mism, it, zmis)
end

"S5: the X_H part of the time-step cap comes from `predict_cap` only (rates averaged over the mixed regions of the last `mix!`), not from the per-zone burn of `burn!`: a convective core is homogenised within the step, so the central zone's own ΔX_H (10× the region's mean) is not what the lagged composition sees."
const REGION_CAP = Ref(false)
const WHICH = Ref{Any}(nothing)
"""
    predict_cap(m) -> dt [s]

Composition/energy time-step criterion on the converged state of the step just taken, from the instantaneous rates of the
network with the current abundances (what the next hydro step will lag): per zone dt ≤ dX_max/|dX_H/dt|, dt ≤ f_e e/ε,
dt ≤ f_e e/(max_i |dY_i/dt| Q_i N_A MeV). Called through `TC.DT_HOOKS` inside `timstp!`.
"""
function predict_cap(m::StarGas; verbose = false)
    m.nuc || return Inf
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dn, tn, egn = TM.v(:dn), TM.v(:tn), TM.v(:egn)
    cap = Inf
    fHs = zeros(b - a + 1); wts = [dn[k] * TM.v(:dvoln)[k] for k in a:b]
    for k in a:b
        Yk = view(m.Y, :, k)
        kr = NN.rates_k(tn[k], dn[k], NN.Plasma(dn[k], Yk), NN.ye_of(Yk); screen = m.screen)
        f = NN.ydot(Yk, kr)
        ε = 0.0
        for (ir, r) in enumerate(NN.NET); ε += NA * MeV * (r.Q - r.Enu) * NN.event(r, kr[ir], Yk); end
        fHs[k-a+1] = NN.AA[IH] * f[IH] + NN.AA[ID] * f[ID]
        c1 = ε != 0 ? m.efrac * abs(egn[k]) / abs(ε) : Inf
        fbs = [abs(f[i]) * Qspec[i] for i in 1:NS] .* (NA * MeV); fb = maximum(fbs)
        c2 = fb > 0 ? m.efrac * abs(egn[k]) / fb : Inf
        verbose && min(c1, c2) <= cap && (WHICH[] = (zone = k, eps = ε, species = NN.SP[argmax(fbs)], cap_eps = c1, cap_bind = c2, T = tn[k], rho = dn[k]))
        cap = min(cap, c1, c2)
    end
    # X_H rate: averaged over the connected mixed (convective) regions of the last mix! (mixing redistributes the change)
    D = m.Dmix; j = 1
    while j <= length(fHs)
        j2 = j
        while j2 < length(fHs) && j2 + a < length(D) && D[a+j2] > 0; j2 += 1; end
        fm = sum(fHs[j:j2] .* wts[j:j2]) / sum(wts[j:j2])
        fm != 0 && (cap = min(cap, m.dXmax / abs(fm)))
        j = j2 + 1
    end
    cap
end

"Install the composition criterion in timstp! (replaces an earlier one)."
function install_dt_hook!(m::StarGas)
    filter!(h -> !(hasproperty(h, :tag) && false), TC.DT_HOOKS)
    empty!(TC.DT_HOOKS)
    push!(TC.DT_HOOKS, () -> min(m.dtcap, predict_cap(m)))
    nothing
end

# ----------------------------------------------------------------------------------------------- mixing
"""
    mix_species!(Y, D, dm, area, dt) -> Y

Conservative implicit diffusion of the columns of `Y[:, 1:n]` (zones) with D[k] at the interface between zones k−1 and k (D[1] and
D[n+1] unused), zone masses dm[k] (g), interface mass-flux area densities area[k] = 4πr_k²ρ_k (g/cm). Returns Y updated in place.
"""
function mix_species!(Y::AbstractMatrix, D::AbstractVector, dm::AbstractVector, area::AbstractVector, dt)
    n = size(Y, 2)
    c = zeros(n + 1)                       # flux coefficient at interfaces 2..n
    for k in 2:n
        D[k] > 0 || continue
        c[k] = area[k]^2 * D[k] / (0.5 * (dm[k-1] + dm[k]))
    end
    lo = zeros(n); di = zeros(n); up = zeros(n)
    for k in 1:n
        di[k] = dm[k] / dt + c[k] + c[k+1]
        k > 1 && (lo[k] = -c[k])
        k < n && (up[k] = -c[k+1])
    end
    # Thomas, one factorisation for all species; solve for the change δ = Y' − Y (A δ = −L Y, flux form) and refine once,
    # so that the round-off of the conserved sum scales with |δ|, not with |Y|
    cp = zeros(n); w = zeros(n)
    w[1] = di[1]; cp[1] = up[1] / di[1]
    for k in 2:n; w[k] = di[k] - lo[k] * cp[k-1]; cp[k] = up[k] / w[k]; end
    ns = size(Y, 1); dp = zeros(n); dp2 = zeros(n); δ = zeros(n); rr = zeros(n); F = zeros(n + 1)
    function solve!(out, rhs)
        dp[1] = rhs[1] / w[1]
        for k in 2:n; dp[k] = (rhs[k] - lo[k] * dp[k-1]) / w[k]; end
        out[n] = dp[n]
        for k in n-1:-1:1; out[k] = dp[k] - cp[k] * out[k+1]; end
    end
    for i in 1:ns
        for k in 2:n; F[k] = c[k] * (Y[i, k] - Y[i, k-1]); end          # F[1] = F[n+1] = 0 (no flux through the ends)
        for k in 1:n; rr[k] = F[k+1] - F[k]; end                          # −(L Y)_k·(−1): dm dY/dt = F_{k+1} − F_k
        solve!(δ, rr)
        for rep in 1:2                                                    # iterative refinement with the exact flux-form residual
            for k in 2:n; F[k] = c[k] * (Y[i, k] + δ[k] - Y[i, k-1] - δ[k-1]); end
            for k in 1:n; rr[k] = F[k+1] - F[k] - dm[k] / dt * δ[k]; end
            solve!(dp2, rr); δ .+= dp2
        end
        for k in 1:n; Y[i, k] += δ[k]; end
    end
    Y
end

"""
    mix!(m; dt, sd = TX.step_diag(dt)) -> NamedTuple

Convective mixing over the step just taken (see the module docstring). `tmax` (S5): partial mixing, at most the fraction of the composition change that shifts T of any zone by ≤ tmax in the fixed-e re-solve (the rest follows next step). `bridge = nb` (S5): convective regions are contiguous: gaps of ≤ nb interfaces (∇ ≤ ∇_ad at a single noisy interface) are filled with D = min of the neighbours, and a core is extended to the centre; without it the intermittently unmixed core drifted in composition and the next mixing re-solved T at fixed e (ΔT/T up to 24 % in the central zone). `capmix = false` (S5): mixing does not limit dt. The net (burn + mixing) change of X_H per zone limited dt to dXmax per step, but a redistribution of composition conserves the mass of every species: when the convective core switched off and on between steps (its flux is zero for ∇ ≤ ∇_ad), the step that mixed a core burnt unmixed for one step moved X_H by 0.1 in a zone and cut dt to 1e-6 yr, which kept the core unmixed: a death spiral (15 and 100 M☉ runs). Returns (cons = max_i |Δ(Σ dm Y_i)|/Σ dm Y_i,
nconv = number of interfaces with D > 0, Dmax, t_turn = min over convective interfaces of ℓ/v, dEmax = max zone energy change).
"""
function mix!(m::StarGas; dt = TF.getscalar(:timen) - TF.getscalar(:timeo), sd = TX.step_diag(dt), D = nothing, capmix = true, bridge = 0, tmax = Inf)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dn, tn, egn, rn, dvoln = TM.v(:dn), TM.v(:tn), TM.v(:egn), TM.v(:rn), TM.v(:dvoln)
    n = b - a + 1
    Di = zeros(n + 1)                       # Di[j] at the interface between zones a+j-2 and a+j-1 (j = 2..n)
    for j in 2:n
        k = a + j - 1
        Di[j] = D === nothing ? (sd.fcv[k] > 0 && sd.vcv[k] > 0 ? sd.vcv[k] * sd.lcv[k] / 3 : 0.0) : D[k]
    end
    if bridge > 0 && D === nothing      # convective regions are contiguous: fill gaps of ≤ `bridge` interfaces and extend a core to the centre
        jf = findfirst(>(0), view(Di, 2:n))
        if jf !== nothing
            jf += 1
            jf <= 2 + 2bridge && (for j in 2:jf-1; Di[j] = Di[jf]; end)
            j = jf
            while j <= n
                if Di[j] > 0; j += 1; continue; end
                j2 = j; while j2 <= n && Di[j2] <= 0; j2 += 1; end
                (j2 <= n && j2 - j <= bridge) && (for jj in j:j2-1; Di[jj] = min(Di[j-1], Di[j2]); end)
                j = j2
            end
        end
    end
    m.Dmix = zeros(length(dn)); for j in 2:n; m.Dmix[a+j-1] = Di[j]; end
    any(>(0), Di) || return (; cons = 0.0, nconv = 0, Dmax = 0.0, t_turn = Inf, dEmax = 0.0)    # (burn!'s local cap stands)
    dm = [4π * dn[a+j-1] * dvoln[a+j-1] for j in 1:n]
    area = zeros(n + 1); for j in 2:n; k = a + j - 1; area[j] = 4π * rn[k]^2 * 0.5 * (dn[k-1] + dn[k]); end
    Y = m.Y[:, a:b]
    tot0 = [sum(Y[i, j] * dm[j] for j in 1:n) for i in 1:NS]
    mix_species!(Y, Di, dm, area, dt)
    tot1 = [sum(Y[i, j] * dm[j] for j in 1:n) for i in 1:NS]
    cons = maximum(abs(tot1[i] - tot0[i]) / tot0[i] for i in 1:NS if tot0[i] > 0; init = 0.0)
    dEmax = 0.0
    if tmax < Inf       # S5: partial mixing. Re-solving T at fixed e and ρ for a large composition change (the boundary of a growing core moves
                        # over a zone: ΔX_H = 0.2) shifts T of that zone by 10–20 %, which the next hydro step must undo (rejected / diverging steps).
                        # Mix only the fraction λ of the change that keeps max|Δ ln T| ≤ tmax; the rest follows in the next steps (Σ dm Y conserved: linear in λ).
        Y0 = m.Y[:, a:b]; ΔY = Y .- Y0
        function maxdT(λ)
            Ys = copy(m.Y); mx = 0.0
            for j in 1:n
                any(ΔY[:, j] .!= 0) || continue
                k = a + j - 1; Ys[:, k] .= Y0[:, j] .+ λ .* ΔY[:, j]
                comp = (t = xyz(view(Ys, :, k)); SE.composition(t[1], t[2], t[3], m.zmix))
                Tn, _ = T_of_e(dn[k], egn[k], tn[k], comp; pairs = m.pairs, coulomb = m.coulomb)
                mx = max(mx, abs(log(Tn / tn[k])))
            end
            mx
        end
        mx1 = maxdT(1.0)
        if mx1 > tmax
            λ = clamp(tmax / mx1, 0.02, 1.0)
            for _ in 1:3; mxl = maxdT(λ); mxl <= 1.05tmax && break; λ *= 0.7tmax / mxl; end
            Y .= Y0 .+ λ .* ΔY; m.lambda_mix = λ
        else
            m.lambda_mix = 1.0
        end
    end
    for j in 1:n
        k = a + j - 1
        any(Y[:, j] .!= m.Y[:, k]) || continue
        m.Y[:, k] .= Y[:, j]
        Tn, _ = T_of_e(dn[k], egn[k], tn[k], zonecomp(m, k); pairs = m.pairs, coulomb = m.coulomb)
        dEmax = max(dEmax, abs(Tn / tn[k] - 1)); tn[k] = Tn
    end
    ghosts!(m)
    tn[a-1] = tn[a]; tn[b+1] = tn[b]; tn[b+2] = tn[b]
    refresh!()
    if m.nuc && capmix             # net (burn + mixing) change of X_H per zone → time-step cap (S5: off, see below)
        dXn = maximum(abs(hmass(m, k) - (NN.AA[IH] * m.Yprev[IH, k] + NN.AA[ID] * m.Yprev[ID, k])) for k in a:b)
        m.dtcap = min(m.capE, dXn > 0 ? dt * m.dXmax / dXn : Inf)
        TF.setscalar!(:dtime, min(TF.getscalar(:dtime), m.dtcap))
    end
    tt = Inf
    for j in 2:n; k = a + j - 1; Di[j] > 0 && (tt = min(tt, sd.lcv[k] / sd.vcv[k])); end
    (; cons, nconv = count(>(0), Di), Dmax = maximum(Di), t_turn = tt, dEmax)
end

# ----------------------------------------------------------------------------------------------- S6a: composition advection
"""
    advect_composition!(Y, M, δ, Yin; order = 2) -> (; M = M′, cfl)   (Y is updated in place)

Conservative transport of the abundances `Y[i, j]` (mol/g, zones j = 1..n of old mass `M[j]`) through a mass grid whose interfaces move
relative to the matter (the adaptive grid, lgrid = 1; or the Lagrangian grid of a star whose mass interfaces are re-set). Interface j
(j = 1..n+1) lies below zone j; `δ[j]` is the change of the mass interior to it over the step. δ[j] > 0: mass δ[j] crosses from zone j into
zone j−1 (δ[1] > 0 leaves through the inner boundary: not allowed when it would need material, δ[1] < 0 is an error); δ[n+1] > 0 is mass
entering from outside with abundances `Yin`. New zone masses M′[j] = M[j] + δ[j+1] − δ[j].
Upwind (donor-cell, `order = 1`) or second order (`order = 2`): the donated slab of a zone has the mean abundance of a linear
reconstruction whose slope is the central difference of the neighbouring zones, scaled by ONE factor θ per zone (Barth–Jespersen over all species,
so that Σ A_i Y_i = 1 is preserved exactly, a common factor times a linear combination) such that the reconstruction stays within the neighbouring
zones' values at both zone edges (monotone, Y ≥ 0). Σ_j M_j Y_ij changes only by the inflow: conservation to round-off. Requires |δ| ≤ the donor's mass
(CFL ≤ 1); returns M′ and the largest ratio |δ|/M_donor (`cfl`) in the NamedTuple.
"""
function advect_composition1!(Y::AbstractMatrix, M::AbstractVector, δ::AbstractVector, Yin::AbstractVector; order = 2)
    ns, n = size(Y); length(M) == n && length(δ) == n + 1 || error("advect_composition!: sizes")
    δ[1] < -1e-14 * M[1] && error("advect_composition!: material would have to enter through the inner boundary (δ[1] = $(δ[1]))")
    # slopes d Y/d m per zone, with a common limiter θ_j
    σ = zeros(ns, n); θ = ones(n)
    if order >= 2 && n >= 3
        for j in 2:n-1
            dm = 0.5 * (M[j-1] + M[j]) + 0.5 * (M[j] + M[j+1])          # distance (in mass) between the centres of the neighbours
            th = 1.0
            for i in 1:ns
                s = (Y[i, j+1] - Y[i, j-1]) / dm
                s == 0 && continue
                lo = min(Y[i, j-1], Y[i, j], Y[i, j+1]); hi = max(Y[i, j-1], Y[i, j], Y[i, j+1])
                d = 0.5 * M[j] * s                                       # change from the centre to either edge
                d > 0 ? (th = min(th, (hi - Y[i, j]) / d, (Y[i, j] - lo) / d)) : (th = min(th, (lo - Y[i, j]) / d, (Y[i, j] - hi) / d))
            end
            th = clamp(th, 0.0, 1.0); θ[j] = th
            for i in 1:ns; σ[i, j] = th * (Y[i, j+1] - Y[i, j-1]) / dm; end
        end
    end
    Yn = zeros(ns, n); Mn = zeros(n); cfl = 0.0
    F = zeros(ns, n + 1)                   # amount of species crossing interface j downward (mol): δ[j] × slab abundance
    for jf in 1:n+1
        d = δ[jf]; d == 0 && continue
        if d > 0
            if jf == n + 1; F[:, jf] .= d .* Yin
            else
                M[jf] > 0 || error("advect_composition!: empty donor")
                cfl = max(cfl, d / M[jf])
                for i in 1:ns; F[i, jf] = d * (Y[i, jf] - σ[i, jf] * (M[jf] - d) / 2); end
            end
        else
            jd = jf - 1                    # donor below the interface (for jf = n+1: outflow through the outer boundary)
            cfl = max(cfl, -d / M[jd])
            for i in 1:ns; F[i, jf] = d * (Y[i, jd] + σ[i, jd] * (M[jd] + d) / 2); end     # d < 0: slab of mass |d| at the top of zone jd
        end
    end
    for j in 1:n
        Mn[j] = M[j] + δ[j+1] - δ[j]
        Mn[j] > 0 || error("advect_composition!: zone $j empties (M = $(M[j]), δ_below = $(δ[j]), δ_above = $(δ[j+1]))")
        for i in 1:ns; Yn[i, j] = (M[j] * Y[i, j] + F[i, j+1] - F[i, j]) / Mn[j]; end
    end
    Y .= Yn
    (; M = Mn, cfl)
end


"""
    advect_composition_down!(Y, M, δ, Yin; order = 2, cmax = 0.5) -> (; M, cfl)

Same transport as `advect_composition1!` for a grid whose interior masses only GROW (δ[j] ≥ 0: matter flows downward through every interface, as in the self-similar
accretion grid), for ANY CFL number: interfaces whose donor would give away more than `cmax` of its mass use the new abundances of the donor (backward-Euler upwind),
the others the explicit second-order flux of `advect_composition1!`. Sweeping from the top zone down, zone j has M′_j = M_j + δ_{j+1} − δ_j and
M′_j Y′_j = M_j Y_j + F_{j+1} − F_j with F_j = δ_j Y′_j (implicit) or δ_j Y_face,j (explicit), i.e. an implicit interface gives Y′_j = (M_j Y_j + F_{j+1})/(M_j + δ_{j+1}).
Unconditionally positive; Σ M Y changes only by the inflow δ[n+1] Yin.
"""
function advect_composition_down!(Y::AbstractMatrix, M::AbstractVector, δ::AbstractVector, Yin::AbstractVector; order = 2, cmax = 0.5)
    ns, n = size(Y); length(M) == n && length(δ) == n + 1 || error("advect_composition_down!: sizes")
    all(>=(-1e-14 * maximum(M)), δ) || error("advect_composition_down!: needs δ ≥ 0")
    δ[1] > 1e-14 * M[1] && error("advect_composition_down!: δ[1] must be 0 (regular centre)")
    σ = zeros(ns, n)
    if order >= 2 && n >= 3
        for j in 2:n-1
            dm = 0.5 * (M[j-1] + M[j]) + 0.5 * (M[j] + M[j+1]); th = 1.0
            for i in 1:ns
                s = (Y[i, j+1] - Y[i, j-1]) / dm; s == 0 && continue
                lo = min(Y[i, j-1], Y[i, j], Y[i, j+1]); hi = max(Y[i, j-1], Y[i, j], Y[i, j+1]); d = 0.5 * M[j] * s
                d > 0 ? (th = min(th, (hi - Y[i, j]) / d, (Y[i, j] - lo) / d)) : (th = min(th, (lo - Y[i, j]) / d, (Y[i, j] - hi) / d))
            end
            th = clamp(th, 0.0, 1.0)
            for i in 1:ns; σ[i, j] = th * (Y[i, j+1] - Y[i, j-1]) / dm; end
        end
    end
    Yn = zeros(ns, n); Mn = zeros(n); cfl = 0.0
    Fin = δ[n+1] .* Yin                      # flux through the top interface (inflow)
    for j in n:-1:1
        Mn[j] = M[j] + δ[j+1] - δ[j]
        Mn[j] > 0 || error("advect_composition_down!: zone $j empties")
        d = max(δ[j], 0.0); c = d / M[j]; cfl = max(cfl, c)
        if j == 1 || d == 0
            for i in 1:ns; Yn[i, j] = (M[j] * Y[i, j] + Fin[i]) / Mn[j]; end
            j > 1 || (Fin = zeros(ns))
            d == 0 && j > 1 && (Fin = zeros(ns))
        elseif c > cmax          # implicit upwind interface
            for i in 1:ns; Yn[i, j] = (M[j] * Y[i, j] + Fin[i]) / (M[j] + δ[j+1]); end
            Fin = d .* Yn[:, j]
        else                      # explicit second-order flux
            Fex = [d * (Y[i, j] - σ[i, j] * (M[j] - d) / 2) for i in 1:ns]
            for i in 1:ns; Yn[i, j] = (M[j] * Y[i, j] + Fin[i] - Fex[i]) / Mn[j]; end
            Fin = Fex
        end
    end
    Y .= Yn
    (; M = Mn, cfl)
end

"""
    advect_composition!(Y, M, δ, Yin; order = 2, cflmax = 0.5) -> (; M, cfl, nsub)

`advect_composition1!` in `nsub = ceil(cfl/cflmax)` equal substeps (interface masses interpolated linearly in time), so that no donor gives away more than
`cflmax` of its mass at a time (the adaptive grid can move an interface through several thin zones in one implicit step).
"""
function advect_composition!(Y::AbstractMatrix, M::AbstractVector, δ::AbstractVector, Yin::AbstractVector; order = 2, cflmax = 0.5)
    n = size(Y, 2)
    if all(>=(-1e-14 * maximum(M)), δ) && abs(δ[1]) <= 1e-14 * M[1]
        r = advect_composition_down!(Y, M, δ, Yin; order, cmax = cflmax)
        return (; M = r.M, cfl = r.cfl, nsub = 1)
    end
    c = 0.0
    for jf in 1:n+1
        d = δ[jf]; d == 0 && continue
        jd = d > 0 ? min(jf, n) : max(jf - 1, 1)
        c = max(c, abs(d) / M[jd])
    end
    # (the donor masses change during the substeps; the estimate is refined by checking the last substep)
    nsub = max(1, ceil(Int, c / cflmax))
    Mc = copy(M); cfl = 0.0
    for isub in 1:nsub
        r = advect_composition1!(Y, Mc, δ ./ nsub, Yin; order); Mc = r.M; cfl = max(cfl, r.cfl)
    end
    (; M = Mc, cfl, nsub)
end

"""
    advect!(m::StarGas, xme_start; xin = nothing, order = 2, fac = nothing)

Transport `m.Y` of the zones ngrs..ngre through the grid motion of the step just taken (lgrid = 1): the interior mass of interface k is
xmtot − xmen[k] (per steradian); `xme_start` are the exterior masses `xmen[ngrs:ngre+1]` and `M_start` the zone masses 4π ρ ΔV (per steradian × 4π)
at the START of the step (saved by the caller). Inflow through the outer interface (leobc = 2) carries `xin` (default: the initial abundances in `m.Yin`).
Then T of every changed zone is re-solved at fixed ρ and gas energy e (as `mix!`), and the EOS/opacity refreshed. Returns (cons, cfl, dTmax):
cons = max_i |Δ(Σ M Y_i) − inflow| / Σ M Y_i.
"""
function advect!(m::StarGas, xme_start::AbstractVector, M_start::AbstractVector; xin = nothing, order = 2, fac = nothing)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dn, tn, egn, dvoln, xmen = TM.v(:dn), TM.v(:tn), TM.v(:egn), TM.v(:dvoln), TM.v(:xmen)
    n = b - a + 1
    δ = [-(xmen[k] - xme_start[k-a+1]) * 4π for k in a:b+1]            # mass (g) gained by the interior of interface k
    if fac !== nothing        # prescribed self-similar grid (ALE_MDOT): delta_k = f m_k exactly (the difference of two x_me of size x_mtot has no significant places left for f below 1e-12)
        xmtot = TF.getscalar(:xmtot)
        δ = [k == a ? 0.0 : (xmtot - xme_start[k-a+1]) * 4π * fac * (isempty(TD.ALE_PSI[]) || k > b ? 1.0 : TD.ALE_PSI[][k]) for k in a:b+1]
    end
    # zone masses from the same interface masses (Δ x_me), so that M′_j = M_j + δ_{j+1} − δ_j holds identically: ρΔV of the Newton solution differs from Δx_me by the
    # converged residual of the mass row (≈1e-10 relative to x_mtot), which for a surface zone of 1e-7 M is a visible fraction of its mass
    M_start = [4π * (xme_start[j] - xme_start[j+1]) for j in 1:n]
    Ysub = m.Y[:, a:b]
    Yin = xin === nothing ? Ysub[:, n] : xin
    tot0 = [sum(M_start[j] * Ysub[i, j] for j in 1:n) for i in 1:NS]
    infl = δ[n+1] .* Yin
    r = advect_composition!(Ysub, copy(M_start), δ, Yin; order)
    tot1 = [sum(r.M[j] * Ysub[i, j] for j in 1:n) for i in 1:NS]
    ref = maximum(tot0 .* NN.AA)       # species whose mass is below 1e-9 of the largest are excluded (trace species: the round-off of the largest zones dominates their relative error)
    cons = maximum(abs(tot1[i] - tot0[i] - infl[i]) / max(tot0[i], 1e-300) for i in 1:NS if tot0[i] * NN.AA[i] > 1e-9 * ref; init = 0.0)
    dT = 0.0
    for j in 1:n
        k = a + j - 1
        any(Ysub[:, j] .!= m.Y[:, k]) || continue
        big = maximum(abs(Ysub[i, j] - m.Y[i, k]) / max(abs(m.Y[i, k]), 1e-30) for i in 1:NS if m.Y[i, k] > 1e-12; init = 0.0) > 1e-10      # (a round-off change of Y needs no re-solve of T)
        m.Y[:, k] .= Ysub[:, j]
        big || continue
        Tn, _ = T_of_e(dn[k], egn[k], tn[k], zonecomp(m, k); pairs = m.pairs, coulomb = m.coulomb)
        dT = max(dT, abs(Tn / tn[k] - 1)); tn[k] = Tn
    end
    ghosts!(m)
    tn[a-1] = tn[a]; tn[b+1] = tn[b]; tn[b+2] = tn[b]
    refresh!()
    (; cons, cfl = r.cfl, dTmax = dT, nsub = r.nsub)
end

# ----------------------------------------------------------------------------------------------- diagnostics
"Total nuclear power ∫ε dm [erg/s] on the current state (total mass, 4π included)."
function nuclear_power(m::StarGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dn, tn, dvoln = TM.v(:dn), TM.v(:tn), TM.v(:dvoln)
    sum(4π * dn[k] * dvoln[k] * eps_nuc(dn[k], tn[k], view(m.Y, :, k); screen = m.screen) for k in a:b)
end
"Total mass Σ dm and the mass of each species (g) in the grid."
function species_mass(m::StarGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dn, dvoln = TM.v(:dn), TM.v(:dvoln)
    [sum(4π * dn[k] * dvoln[k] * NN.AA[i] * m.Y[i, k] for k in a:b) for i in 1:NS]
end

end # module
