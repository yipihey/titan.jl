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
end
function StarGas(n::Integer; X = 0.76, Y = 0.24, Z = 0.0, zmix = SE.ZMIX, tab = Z0.table(), opac = :stellar, grey = false,
                 coulomb = true, pairs = true, kcap = Inf, nuc = false, screen = :chugunov, burnT = :end, dXmax = 1e-3, efrac = 0.05, mixing = true)
    Yv = zeros(NS, n)
    s = zmix.C + zmix.N + zmix.O
    for k in 1:n
        Yv[IH, k] = X; Yv[NN.IDX[:he4], k] = Y / 4
        Yv[NN.IDX[:c12], k] = Z * zmix.C / s / 12; Yv[NN.IDX[:n14], k] = Z * zmix.N / s / 14; Yv[NN.IDX[:o16], k] = Z * zmix.O / s / 16
    end
    StarGas(Yv, zmix, tab, opac, grey, coulomb, pairs, kcap, nuc, screen, burnT, dXmax, efrac, Inf, Inf, copy(Yv), mixing,
            (Einj = 0.0, Ebur = 0.0, dabs = 0.0, Eabs = 0.0, resid = 0.0, nsteps = 0.0), (;), zeros(n + 2))
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

@inline function TX.material_gas(m::StarGas, ρ, T)
    g = SE.gas_state(ρ, T, zonecomp(m); pairs = m.pairs, coulomb = m.coulomb)
    return g.p, g.e, sqrt(5.0 / 3.0 * g.p / ρ)
end
@inline function TX.material_gasne(m::StarGas, ρ, T)
    g = SE.gas_state(ρ, T, zonecomp(m); pairs = m.pairs, coulomb = m.coulomb)
    return g.p, g.e, sqrt(5.0 / 3.0 * g.p / ρ), g.ne
end
"C¹ soft cap: identity below 0.7c, saturating to c above."
@inline softcap(κ, c) = isinf(c) ? κ : (κ <= 0.7c ? κ : c - 0.3c * exp(-(κ - 0.7c) / (0.3c)))
@inline function TX.material_opac(m::StarGas, ρ, T)
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
TX.material_source(m::StarGas, ρ, T) = m.nuc ? eps_nuc(ρ, T, view(m.Y, :, TX.ZONE[]); screen = m.screen) : zero(T)

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

# ----------------------------------------------------------------------------------------------- burn!
"""
    burn!(m::StarGas; dt = timen − timeo) -> NamedTuple

After a converged step. See the module docstring. Returns the step's ledger entry (Einj, Ebur, corr, resid, max ΔX_H, dtcap).
"""
function burn!(m::StarGas; dt = TF.getscalar(:timen) - TF.getscalar(:timeo), sd = nothing)
    m.nuc || return (; Einj = 0.0, Ebur = 0.0, corr = 0.0, resid = 0.0, dXH = 0.0, dtcap = Inf)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dn, tn, egn, to = TM.v(:dn), TM.v(:tn), TM.v(:egn), TM.v(:to)
    dvoln = TM.v(:dvoln)
    θ = TF.getscalar(:thet)
    sd === nothing && (sd = TX.step_diag(dt))    # srcw: what the residual injected (erg/s per steradian)
    m.Yprev .= m.Y
    Einj = 0.0; Ebur = 0.0; Eabs = 0.0; dcorr = 0.0; resid = 0.0; dXH = 0.0; dEfrac = 0.0; dEsp = 0.0
    for k in a:b
        ρ = dn[k]; mk = 4π * dn[k] * dvoln[k]
        Tb = m.burnT === :end ? tn[k] : θ * tn[k] + (1 - θ) * to[k]
        Y0 = m.Y[:, k]
        res = NN.burn(Y0, ρ, Tb, dt; screen = m.screen)
        eb = res.ε * dt                          # erg/g released (net of ν), time-averaged over substeps
        ei = 4π * sd.srcw[k] * dt                # erg injected in the zone by the residual
        Einj += ei; Ebur += eb * mk; Eabs += abs(eb * mk)
        dcorr += abs(eb * mk - ei)
        # energy: e_new = e_n + (Ebur − Einj)/m at fixed ρ; T re-solved with the new composition
        et = egn[k] + (eb * mk - ei) / mk
        m.Y[:, k] .= res.Y
        comp = zonecomp(m, k)
        Tn, er = T_of_e(ρ, et, tn[k], comp; pairs = m.pairs, coulomb = m.coulomb)
        resid = max(resid, er); tn[k] = Tn
        dXH = max(dXH, abs(NN.AA[IH] * (res.Y[IH] - Y0[IH]) + NN.AA[ID] * (res.Y[ID] - Y0[ID])))
        for i in 1:NS
            dEsp = max(dEsp, abs(res.Y[i] - Y0[i]) * Qspec[i] * NA * MeV / abs(et))
        end
        dEfrac = max(dEfrac, abs(eb) / abs(et))
    end
    ghosts!(m)
    tn[a-1] = tn[a]; tn[b+1] = tn[b]; tn[b+2] = tn[b]
    refresh!()
    cap = Inf
    dXH > 0 && (cap = min(cap, dt * m.dXmax / dXH))
    dE = max(dEsp, dEfrac)
    capE = dE > 0 ? dt * m.efrac / dE : Inf
    m.capE = capE; cap = min(cap, capE)
    m.dtcap = cap
    TF.setscalar!(:dtime, min(TF.getscalar(:dtime), cap))
    corr = Ebur == 0 ? 0.0 : dcorr / abs(Ebur)
    L = m.ledger
    m.ledger = (Einj = L.Einj + Einj, Ebur = L.Ebur + Ebur, dabs = L.dabs + dcorr, Eabs = L.Eabs + Eabs, resid = max(L.resid, resid), nsteps = L.nsteps + 1)
    m.last = (; Einj, Ebur, corr, resid, dXH, dtcap = cap, dEsp, dEfrac)
end

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

Convective mixing over the step just taken (see the module docstring). Returns (cons = max_i |Δ(Σ dm Y_i)|/Σ dm Y_i,
nconv = number of interfaces with D > 0, Dmax, t_turn = min over convective interfaces of ℓ/v, dEmax = max zone energy change).
"""
function mix!(m::StarGas; dt = TF.getscalar(:timen) - TF.getscalar(:timeo), sd = TX.step_diag(dt), D = nothing)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dn, tn, egn, rn, dvoln = TM.v(:dn), TM.v(:tn), TM.v(:egn), TM.v(:rn), TM.v(:dvoln)
    n = b - a + 1
    Di = zeros(n + 1)                       # Di[j] at the interface between zones a+j-2 and a+j-1 (j = 2..n)
    for j in 2:n
        k = a + j - 1
        Di[j] = D === nothing ? (sd.fcv[k] > 0 && sd.vcv[k] > 0 ? sd.vcv[k] * sd.lcv[k] / 3 : 0.0) : D[k]
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
    if m.nuc                       # net (burn + mixing) change of X_H per zone → time-step cap
        dXn = maximum(abs(hmass(m, k) - (NN.AA[IH] * m.Yprev[IH, k] + NN.AA[ID] * m.Yprev[ID, k])) for k in a:b)
        m.dtcap = min(m.capE, dXn > 0 ? dt * m.dXmax / dXn : Inf)
        TF.setscalar!(:dtime, min(TF.getscalar(:dtime), m.dtcap))
    end
    tt = Inf
    for j in 2:n; k = a + j - 1; Di[j] > 0 && (tt = min(tt, sd.lcv[k] / sd.vcv[k])); end
    (; cons, nconv = count(>(0), Di), Dmax = maximum(Di), t_turn = tt, dEmax)
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
