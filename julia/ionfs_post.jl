# IonFS diagnostics as a post-process of a run's snapshots: does ionizing or Lyman–Werner light get ahead of the photosphere?
#
#   julia -t 4 --project=<env> ionfs_post.jl <run dir> [snaps_*.jls ...]      (default: every <run dir>/snaps_*.jls)
#
# Appends one line per snapshot not yet done to <run dir>/ionfs_post.csv (incremental: steps already in the file are skipped).
# Per snapshot, the same construction as the in-run IONFSIMP diagnostics (PopIIIMaterial IonFS.jl / IonFSImp.jl, ionfs_snap.csv),
# without the hydro coupling:
#  1. grey τ from the outer boundary (snapshot chi ρ Δr); r_ph = the first interface outside τ = 2/3 (and outside the table zones),
#     L_run = 4π r_ph² F_r(r_ph), T_eff from it.
#  2. source: the static spherical formal solution (Spectrum.monochromatic: tangent rays, Feautrier, ALI with Thomson scattering, thermal
#     inner boundary) from the innermost zone with grey τ > 3e3 (at least 6 zones) to r_ph. Opacity: ContinuumOpacity.alpha (CIA, H⁻,
#     Rayleigh, Thomson, H bound–free with Milne emission of the zone's n_e n_p) + H⁺/He⁺/He²⁺ free–free + He I / He II ground bound–free
#     (Verner, ChemRates.σ_verner), from the snapshot's composition. IONFSSUB: the photosphere inside the front zone as sub-points
#     (grey Eddington T⁴ = ¾ T_eff⁴ (τ + 2/3) on the zone's grey τ scale, source opacity scaled to it). Groups H I 13.598–24.587,
#     He I –54.418, He II –300 eV, LW 11.2–13.598 eV (8 frequencies each), and 24 below 11.2 eV for L_fs; photon numbers scaled by
#     fnorm = L_run/L_fs (the formal solution lacks line opacity in the optical windows, see docs/POPIII_IONFS_IMPLICIT.md).
#  3. ray: photon-conserving attenuation through the snapshot's absorbers zone by zone outward (H, He, He⁺ with the group cross sections
#     photon-weighted over the source spectrum); LW: the photon-conserving H₂ self-shielding of IonFSImp (k = σ_d φ f_sh(N) f_HI, the
#     carrier's used photons giving the equivalent static column; Draine & Bertoldi 1996 ⟨p_diss⟩ = 0.15; Wolcott-Green et al. 2011
#     f_sh with Sobolev, Wolcott-Green & Haiman 2011 f_HI), and STARLIGHT's SLSHIELD rate at the same zone for comparison
#     (1.1e8·4π W B_ν(12.87 eV, T_c) f_sh f_HI with SLBLEND/SLEB/SLWFLUX as the runs; TCMB from ENV, default 54.5 K).
#  4. break-out alarm (see columns alarm_*): radius r_a = max(r_ph + 10 R☉, 2 r_ph).
# Columns (radii R☉, rates s⁻¹, heating erg g⁻¹ s⁻¹):
#  step, t_yr, Mstar_Msun (enclosed mass at r_ph), r_ph_Rsun, T_eff, T_ph, L_run_Lsun, L_fs_Lsun, fnorm,
#  Qem_LW, Qem_H, Qem_HeI, Qem_HeII (leaving r_ph), Qesc_* (through the outer boundary),
#  rIF_H_Rsun, rIF_HeI_Rsun, rIF_HeII_Rsun (first interface beyond which x_H⁺ < ½, (x_He⁺+x_He²⁺)/y_He < ½, x_He²⁺/y_He < ½; = r_ph: none),
#  rH2_01_Rsun, rH2_05_Rsun (first zone centre outside r_ph with f_H₂ = 2x_H₂ ≥ 0.1, 0.5; NaN: none),
#  Dray, Dsl (H₂ dissociations per s outside r_ph: ray, SLSHIELD), Qion_heat, LW_heat (erg/s outside r_ph),
#  r_a_Rsun, QH_a, QLW_a (photons crossing r_a), kLWray_a, tflow_a (s), alarm_H (Q_H(r_a) > ALARMQH, default 1e42 s⁻¹, or rIF_H > r_a),
#  alarm_LW (k_LW,ray(r_a) t_flow(r_a) > ALARMLW, default 0.1: LW dissociates the infall faster than 10 flow times), alarm (either),
#  then per R in 1e2, 1e3, 1e4, 1e5: QLW_R, QH_R, QHeI_R, QHeII_R, T_R, u_kms_R, xe_R, fH2_R, kLWray_R, kLWsl_R, qray_R, tflow_theat_R,
#  and tsolve_s (wall time of the snapshot).
using Serialization, Printf, LinearAlgebra, ForwardDiff
isdefined(Main, :Z0Physics) || include(joinpath(@__DIR__, "Z0Physics.jl"))
isdefined(Main, :Spectrum) || include(joinpath(@__DIR__, "Spectrum.jl"))
isdefined(Main, :ContinuumOpacity) || include(joinpath(@__DIR__, "ContinuumOpacity.jl"))
ENV["HECHEM"] = "1"
isdefined(Main, :PrimordialChem) || include(joinpath(@__DIR__, "PrimordialChem.jl"))
isdefined(Main, :ChemNet) || include(joinpath(@__DIR__, "ChemNet.jl"))
isdefined(Main, :ChemRates) || include(joinpath(@__DIR__, "ChemRates.jl"))
module IFPost
import ..Spectrum, ..ContinuumOpacity, ..ChemNet, ..ChemRates
using Printf, Serialization
const SP = Spectrum; const CO = ContinuumOpacity; const CN = ChemNet; const CR = ChemRates
const hP = CN.hP; const kB = CN.kB; const cc = CN.c; const me = CN.me; const EV = 1.602176634e-12
const RS = 6.957e10; const LS = 3.83e33; const MS = 1.989e33; const σSB = 5.670374419e-5; const arad = 7.565723e-15
const EDGES = [13.598, 24.587, 54.418, 300.0]; const LWB = (11.2, 13.598); const NSUB = 8
const PD = 0.15; const K27 = 1.1e8; const E12 = 12.87
const TCMB = parse(Float64, get(ENV, "TCMB", "54.5"))
const ALARMQH = parse(Float64, get(ENV, "ALARMQH", "1e42")); const ALARMLW = parse(Float64, get(ENV, "ALARMLW", "0.1"))
const RADII = (1e2, 1e3, 1e4, 1e5)
σs(s, ν) = s === :H ? CR.σ_verner(ν, 13.598, CR.VERNER_HI) : s === :He ? CR.σ_verner(ν, 24.587, CR.VERNER_HeI) :
           s === :Hep ? CR.σ_verner(ν, 54.418, CR.VERNER_HeII) : 0.0
const SPC = (:H, :He, :Hep)
"Source opacity at the formal solution's points (as PopIIIMaterial.ifs_coeffs3), with a per-point factor sc."
struct POpac
    T::Vector{Float64}; nH0::Vector{Float64}; nHp::Vector{Float64}; ne::Vector{Float64}; nHm::Vector{Float64}; nH2::Vector{Float64}
    nHe0::Vector{Float64}; nHe1::Vector{Float64}; nHe2::Vector{Float64}; sc::Vector{Float64}
end
const IF_E0HM = 0.754 * EV; const IF_χH = 13.598 * EV; const IF_σ0 = 7.906e-18; const IF_χHe1 = 24.587387 * EV; const IF_χHe2 = 54.417763 * EV
function coeffs3(o::POpac, k, ν)
    T = o.T[k]; kT = kB * T; hν = hP * ν; c2 = 2hP / cc^2
    a, s, p = CO.alpha(ν, T; nH = o.nH0[k], nH2 = o.nH2[k], nHe = o.nHe0[k], nHp = o.nHp[k], ne = o.ne[k], nHm = o.nHm[k])
    B = CO.Bν(ν, T); stim = -expm1(-hν / kT)
    ab = a - p.hm_bf - p.h + 3.692e8 * o.ne[k] * (o.nHp[k] + o.nHe1[k] + 4o.nHe2[k]) / (sqrt(T) * ν^3) * stim
    em = ab * B
    λe3 = (hP^2 / (2π * me * kT))^1.5
    if hν > IF_E0HM
        σ = CO.σ_Hm_bf(1e4 * cc / ν)
        ab += o.nHm[k] * σ * stim; em += o.nH0[k] * o.ne[k] * λe3 / 4 * c2 * σ * ν^3 * exp(-(hν - IF_E0HM) / kT)
    end
    nep = o.ne[k] * o.nHp[k]; lnep = nep > 0 ? log(nep * λe3) : -Inf; ln1 = o.nH0[k] > 0 ? log(o.nH0[k]) : -Inf
    for n in 1:30
        χn = IF_χH / n^2; hν >= χn || continue
        σ = n == 1 ? σs(:H, ν) : IF_σ0 * n * (χn / hν)^3
        ls = log(n^2) + lnep + (χn - hν) / kT; nst = ls > -700 ? exp(ls) : 0.0
        nn = n == 1 ? o.nH0[k] : n^2 * exp(min(lnep + χn / kT, ln1 - (IF_χH - χn) / kT))
        ab += σ * max(nn - nst, 0.0); em += σ * nst * c2 * ν^3
    end
    for (n0, ni, χ, gr, sp) in ((o.nHe0[k], o.nHe1[k], IF_χHe1, 1 / 4, :He), (o.nHe1[k], o.nHe2[k], IF_χHe2, 1.0, :Hep))
        hν >= χ || continue
        σ = σs(sp, ν); ls = o.ne[k] * ni > 0 ? log(o.ne[k] * ni * gr * λe3) + (χ - hν) / kT : -Inf; nst = ls > -700 ? exp(ls) : 0.0
        ab += σ * max(n0 - nst, 0.0); em += σ * nst * c2 * ν^3
    end
    ab = max(ab, 1e-300)
    (ab * o.sc[k], s * o.sc[k], em / ab)
end
function freqs()
    νs = Float64[]; dν = Float64[]; grp = Int[]
    Elo = exp.(range(log(0.05), log(LWB[1]), length = 25))
    for i in 1:24; push!(νs, sqrt(Elo[i] * Elo[i+1]) * EV / hP); push!(dν, (Elo[i+1] - Elo[i]) * EV / hP); push!(grp, 0); end
    for (b_, gid) in zip(vcat([LWB], [(EDGES[g], EDGES[g+1]) for g in 1:3]), (4, 1, 2, 3))
        e1, e2 = b_; q = (e2 / e1)^(1 / NSUB)
        for i in 1:NSUB; push!(νs, e1 * q^(i - 0.5) * EV / hP); push!(dν, e1 * (q^i - q^(i - 1)) * EV / hP); push!(grp, gid); end
    end
    push!(νs, E12 * EV / hP); push!(dν, 0.0); push!(grp, -1)
    νs, dν, grp
end
const FREQ = freqs()
# ---------------------------------------------------------------- H₂ self-shielding (as PopIIIMaterial / IonFSImp)
fsh_h2(N, T) = (x = N / 5e14; b5 = sqrt(kB * T / 1.6735575e-24) / 1e5; 0.965 / (1 + x / b5)^1.1 + 0.035 / sqrt(1 + x) * exp(-8.5e-4 * sqrt(1 + x)))
fsh_hi(N) = (x = N / 2.85e23; (1 + x)^-1.6 * exp(-0.15x))
function Gsh(N, T)
    x = N / 5e14; b5 = sqrt(kB * T / 1.6735575e-24) / 1e5; a = 8.5e-4
    5e14 * (9.65 * b5 * (1 - (1 + x / b5)^(-0.1)) + (0.07 / a) * (exp(-a) - exp(-a * sqrt(1 + x))))
end
function ginv(g, T)
    g <= 0 && return 0.0
    lo = g; hi = 1e30; Gsh(hi, T) <= g && return hi
    for _ in 1:90; mid = sqrt(lo * hi); Gsh(mid, T) < g ? (lo = mid) : (hi = mid); hi / lo < 1 + 1e-13 && break; end
    sqrt(lo * hi)
end
"LW in a zone: (Δτ_LW of the carrier, dissociations per s) from the carrier N_LW entering, Φ, σ_d, the zone's x_H₂ n_H Δr, T, sob, f_HI."
function lwzone(N, Φ, σd, ΔN, T, Nsob, fhi)
    (Φ > 0 && σd > 0 && N > 0 && ΔN > 0) || return (0.0, 0.0)
    q = clamp(1 - N / Φ, 0.0, 1.0); N1 = ginv(q * PD / σd, T)
    fst = ΔN < 1e-4 * max(N1, 5e14) ? fsh_h2(N1 + ΔN / 2, T) : (Gsh(N1 + ΔN, T) - Gsh(N1, T)) / ΔN
    fe = isfinite(Nsob) ? max(fst, fsh_h2(Nsob, T)) : fst
    A = σd * Φ / PD * ΔN * fe; τ = A / N
    (τ, PD * N * (-expm1(-τ)) * fhi)
end
_bnu(ν, T) = (x = hP * ν / (kB * T); x > 700 ? 0.0 : 2hP * ν^3 / cc^2 / expm1(x))

# ---------------------------------------------------------------- one snapshot
function process(s)
    t0 = time()
    n = length(s.ρ); r = s.r; ρ = s.ρ; T = s.T; sp = s.spn; μH = CN.μH(); nH = ρ ./ μH; yHe = CN.yHe
    Y(q) = (i = findfirst(==(q), sp); i === nothing ? zeros(n) : max.(s.Y[i, :], 0.0))
    xH, xH2, xHe, xHep, xHepp, xHp, xHm, xe = Y(:H), Y(:H2), Y(:He), Y(:Hep), Y(:Hepp), Y(:Hp), Y(:Hm), Y(:e)
    τz = zeros(n + 1); for k in n:-1:1; c = s.chi[k]; τz[k] = τz[k+1] + (isfinite(c) && c > 0 ? c : 0.0) * ρ[k] * (r[k+1] - r[k]); end
    kph = 1; for k in n:-1:1; τz[k] >= 2 / 3 && (kph = k; break); end
    k0 = kph + 1; while k0 <= n && s.eq[k0]; k0 += 1; end            # (the first zone outside the photosphere and the table zones)
    k0 > n - 1 && return nothing
    rb = r[k0]; Lrun = 4π * rb^2 * s.fr[k0]
    kc = 1; for k in 1:k0-1; τz[k] > 3e3 && (kc = k); end
    kc = max(1, min(kc, k0 - 6))
    atm0 = SP.atmosphere((r = r[kc:k0], d = ρ[kc:k0-1], T = T[kc:k0-1]); split = 3.0)
    zone = copy(atm0.zone .+ (kc - 1)); rr = copy(atm0.r); jump = copy(atm0.jump); Tp = copy(atm0.T)
    kf = k0 - 1; χf = s.chi[kf] * ρ[kf]; τzc = χf * (rb - rr[end]); nsub = 0
    Te = Lrun > 0 ? (Lrun / (4π * rb^2 * σSB))^0.25 : 0.0
    if Lrun > 0 && τzc > 1
        rcf = rr[end]
        for τ in (300.0, 100.0, 30.0, 10.0, 3.0, 1.0, 0.3, 0.1, 0.03, 0.01, 0.001, 0.0)
            τ < 0.9τzc || continue
            push!(rr, τ == 0 ? rb * (1 - 1e-9) : rb - (rb - rcf) * τ / τzc); push!(zone, kf); push!(jump, false)
            push!(Tp, min(T[kf], Te * (0.75 * (τ + 2 / 3))^0.25)); nsub += 1
        end
    else
        push!(zone, kf); push!(rr, rb * (1 - 1e-9)); push!(jump, false); push!(Tp, T[kf])
    end
    g(x) = [x[z] * nH[z] for z in zone]
    op = POpac(Tp, g(xH), g(xHp), g(xe), g(xHm), g(xH2), g(xHe), g(xHep), g(xHepp), ones(length(rr)))
    νs, dν, grp = FREQ
    if nsub > 0
        for j in length(rr)-nsub+1:length(rr)
            sR = 0.0; wR = 0.0
            for i in eachindex(νs)
                grp[i] >= 0 || continue
                ab, sc_, _ = coeffs3(op, j, νs[i]); x = hP * νs[i] / (kB * Tp[j])
                dB = x > 600 ? 0.0 : CO.Bν(νs[i], Tp[j]) * x * exp(x) / expm1(x) / Tp[j]
                sR += dB * dν[i] / (ab + sc_); wR += dB * dν[i]
            end
            sR > 0 && (op.sc[j] = χf / (wR / sR))
        end
    end
    atm = (r = rr, ρ = ρ[zone], T = Tp, zone = zone, jump = jump)
    _, R = SP.rays(atm.r; ncore = 8); R = SP.weights!(atm, R)
    L = zeros(length(νs))
    Threads.@threads for i in eachindex(νs)
        L[i] = max(SP.monochromatic(atm, op, R, νs[i]; maxit = 60, maxit_direct = true)[1], 0.0)
    end
    Lfs = sum(L[i] * dν[i] for i in eachindex(νs) if grp[i] >= 0)
    N = zeros(4); Ec = zeros(4); σ = zeros(3, 3)
    for i in eachindex(νs)
        gi = grp[i]; gi <= 0 && continue
        n_ = L[i] * dν[i] / (hP * νs[i]); N[gi] += n_; Ec[gi] += n_ * hP * νs[i]
        gi <= 3 && for (j, q) in enumerate(SPC); σ[gi, j] += n_ * σs(q, νs[i]); end
    end
    for gg in 1:4
        if N[gg] > 0; Ec[gg] /= N[gg]; gg <= 3 && (σ[gg, :] ./= N[gg])
        else
            e1, e2 = gg <= 3 ? (EDGES[gg], EDGES[gg+1]) : LWB; Ec[gg] = sqrt(e1 * e2) * EV
            gg <= 3 && for (j, q) in enumerate(SPC); σ[gg, j] = σs(q, Ec[gg] / hP); end
        end
    end
    f = (Lfs > 0 && Lrun > 0) ? Lrun / Lfs : 1.0
    σd = N[4] > 0 ? K27 * L[end] / N[4] : 0.0; Φ = f * N[4]; Nem = f .* N
    # star light for the SLSHIELD comparison (STARLIGHT's starlight_setup: SLBLEND, SLEB, SLWFLUX)
    Ecmb = arad * TCMB^4
    trtau(τ) = begin
        τ >= τz[1] && return (max(s.er[1], Ecmb) / arad)^0.25
        for j in 1:n-1
            if τz[j] >= τ > τz[j+1]
                t1, t2 = log(max(τz[j], 1e-30)), log(max(τz[j+1], 1e-30)); w = t1 > t2 ? (log(τ) - t2) / (t1 - t2) : 0.0
                l1, l2 = log(max(s.er[j], Ecmb, 1e-300)), log(max(s.er[j+1], Ecmb, 1e-300))
                return (exp(l2 + w * (l1 - l2)) / arad)^0.25
            end
        end
        (max(s.er[n], Ecmb) / arad)^0.25
    end
    # ray outward from k0
    Nk = copy(Nem); Qint = fill(NaN, 4, n + 1); Dz = zeros(n); kray = zeros(n); ksl = zeros(n); qion = zeros(n)
    sH2 = 0.0; sHI = 0.0
    for k in k0:n
        for g_ in 1:4; Qint[g_, k] = Nk[g_]; end
        dr = r[k+1] - r[k]; V = 4π / 3 * (r[k+1]^3 - r[k]^3); nHk = nH[k]
        τi = [(xH[k] * σ[g_, 1] + xHe[k] * σ[g_, 2] + xHep[k] * σ[g_, 3]) * nHk * dr for g_ in 1:3]
        fhi = fsh_hi(sHI + 0.5 * nHk * xH[k] * dr)
        dudr = abs(s.u[k+1] - s.u[k]) / dr; vth = sqrt(2kB * T[k] / 1.6735575e-24)
        sobL = dudr > 0 ? vth / sqrt(2) / dudr : Inf
        τL, D = lwzone(Nk[4], Φ, σd, xH2[k] * nHk * dr, T[k], xH2[k] * nHk * sobL, fhi)
        Dz[k] = D; kray[k] = xH2[k] > 0 ? D / (xH2[k] * nHk * V) : 0.0
        for g_ in 1:3; qion[k] += Nk[g_] * (-expm1(-τi[g_])) * Ec[g_] / V; end
        # SLSHIELD rate at this zone (STARLIGHT: column to the zone centre, static or Sobolev)
        NH2c = min(sH2 + 0.5 * nHk * xH2[k] * dr, xH2[k] * nHk * sobL)
        Tr = (max(s.er[k], Ecmb) / arad)^0.25; Tsrc = trtau(τz[k] + 2 / 3)
        wq = clamp((log10(max(τz[k], 1e-30)) + log10(1 / 0.3)) / 1.0, 0.0, 1.0); wq = wq * wq * (3 - 2wq)
        Tc = Tsrc^(1 - wq) * Tr^wq
        cu = max(cc * (s.er[k] - Ecmb), 0.5 * (max(s.fr[k], 0.0) + max(s.fr[k+1], 0.0)))
        W = wq >= 1 ? 1.0 : min(max(cu, 0.0) / (4σSB * Tc^4), 1.0)
        ksl[k] = K27 * 4π * W * _bnu(E12 * EV / hP, Tc) * fsh_h2(NH2c, T[k]) * fhi
        sH2 += nHk * xH2[k] * dr; sHI += nHk * xH[k] * dr
        for g_ in 1:3; Nk[g_] *= exp(-min(τi[g_], 700.0)); end
        Nk[4] *= exp(-min(τL, 700.0))
    end
    for g_ in 1:4; Qint[g_, n+1] = Nk[g_]; end
    front(c) = (for k in k0:n; c(k) || return r[k] / RS; end; r[n+1] / RS)
    rIH = front(k -> xHp[k] >= 0.5); rIHe1 = front(k -> (xHep[k] + xHepp[k]) / yHe >= 0.5); rIHe2 = front(k -> xHepp[k] / yHe >= 0.5)
    h2f(th) = (for k in k0:n; 2xH2[k] >= th && return 0.5 * (r[k] + r[k+1]) / RS; end; NaN)
    Dr = sum(Dz[k0:n]); Ds = sum(ksl[k] * xH2[k] * nH[k] * 4π / 3 * (r[k+1]^3 - r[k]^3) for k in k0:n)
    Pion = sum(qion[k] * 4π / 3 * (r[k+1]^3 - r[k]^3) for k in k0:n); Plw = Dr * (4.48 + 0.4) * EV
    Qat(g_, ki) = (q = Qint[g_, ki]; isnan(q) ? 0.0 : q)
    ifc(Rcm) = (i = findfirst(x -> x >= Rcm, r[k0:n+1]); i === nothing ? n + 1 : i + k0 - 1)
    ra = max(rb + 10RS, 2rb); ka = ifc(ra); kza = clamp(ka - 1, k0, n)
    tfa = 0.5 * (r[kza] + r[kza+1]) / max(abs(0.5 * (s.u[kza] + s.u[kza+1])), 1e-30)
    alH = Qat(2, ka) > ALARMQH || rIH * RS > ra; alLW = kray[kza] * tfa > ALARMLW
    Mst = (hasproperty(s, :m) ? s.m[k0] : NaN) / MS
    cols = Any[s.step, s.t / 3.156e7, Mst, rb / RS, Te, T[kf], Lrun / LS, Lfs / LS, f,
               (Nem[g_] for g_ in (4, 1, 2, 3))..., (Qat(g_, n + 1) for g_ in (4, 1, 2, 3))..., rIH, rIHe1, rIHe2, h2f(0.1), h2f(0.5),
               Dr, Ds, Pion, Plw, ra / RS, Qat(1, ka), Qat(4, ka), kray[kza], tfa, Int(alH), Int(alLW), Int(alH || alLW)]
    for Rr in RADII
        ki = ifc(Rr * RS); kz = clamp(ki - 1, k0, n)
        q = (qion[kz] + Dz[kz] * (4.48 + 0.4) * EV / (4π / 3 * (r[kz+1]^3 - r[kz]^3))) / ρ[kz]
        uR = 0.5 * (s.u[kz] + s.u[kz+1]); tfl = 0.5 * (r[kz] + r[kz+1]) / max(abs(uR), 1e-30)
        append!(cols, [(Qat(g_, ki) for g_ in (4, 1, 2, 3))..., T[kz], uR / 1e5, xe[kz], 2xH2[kz], kray[kz], ksl[kz], q, tfl / (s.e[kz] / max(q, 1e-300))])
    end
    push!(cols, time() - t0)
    cols
end
function header()
    h = ["step", "t_yr", "Mstar_Msun", "r_ph_Rsun", "T_eff", "T_ph", "L_run_Lsun", "L_fs_Lsun", "fnorm",
         "Qem_LW", "Qem_H", "Qem_HeI", "Qem_HeII", "Qesc_LW", "Qesc_H", "Qesc_HeI", "Qesc_HeII",
         "rIF_H_Rsun", "rIF_HeI_Rsun", "rIF_HeII_Rsun", "rH2_01_Rsun", "rH2_05_Rsun", "Dray", "Dsl", "Qion_heat", "LW_heat",
         "r_a_Rsun", "QH_a", "QLW_a", "kLWray_a", "tflow_a", "alarm_H", "alarm_LW", "alarm"]
    for Rr in RADII
        t = string(Int(Rr)); append!(h, ["QLW_$t", "QH_$t", "QHeI_$t", "QHeII_$t", "T_$t", "u_kms_$t", "xe_$t", "fH2_$t",
                                        "kLWray_$t", "kLWsl_$t", "qray_$t", "tflow_theat_$t"])
    end
    push!(h, "tsolve_s"); h
end
function run(dir, files)
    out = joinpath(dir, "ionfs_post.csv")
    done = Set{Int}()
    if isfile(out)
        for l in eachline(out); (isempty(l) || startswith(l, "step")) && continue; push!(done, parse(Int, split(l, ",")[1])); end
    end
    CO.tables()
    nnew = 0; nskip = 0
    for f in files
        open(f) do io
            while !eof(io)
                s = try deserialize(io) catch; break end
                s.step in done && (nskip += 1; continue)
                c = try process(s) catch e
                    @printf("ionfs_post: step %d failed: %s\n", s.step, first(sprint(showerror, e), 200)); nothing
                end
                c === nothing && continue
                new = !isfile(out)
                open(out, "a") do o
                    new && println(o, join(header(), ","))
                    println(o, join((x isa Integer ? string(x) : @sprintf("%.6e", x) for x in c), ","))
                end
                push!(done, s.step); nnew += 1
            end
        end
    end
    @printf("ionfs_post: %s: %d snapshots added, %d already done\n", dir, nnew, nskip)
end
end # module
@eval Spectrum coeffs3(o::$(IFPost.POpac), k, ν) = $(IFPost.coeffs3)(o, k, ν)
let dir = ARGS[1]
    files = length(ARGS) > 1 ? ARGS[2:end] : sort(filter(f -> startswith(basename(f), "snaps_") && endswith(f, ".jls"), readdir(dir; join = true)))
    IFPost.run(dir, files)
end
