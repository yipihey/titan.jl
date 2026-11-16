# Ionizing-continuum precursor of an accreting Pop III protostar snapshot: where the photons above 13.6 / 24.6 / 54.4 eV come
# from, how far they get, and what H(1s) / He photo-ionization rates they would give the network zones. Post-processing only.
#
#   julia -t 8 --project=<env> ionprec.jl <run dir> <snaps_*.jls> <step|last> [out.jls]
#
# 1. The explorer's own continuum (snapspec cache <run>/spec/<key>.jls, if present): Q above 13.6 eV at the edge and through r_ph.
# 2. A formal solution on an EUV grid (13.6–300 eV) with two opacity models, both from the code:
#    (a) SnapSpec.Opac (the explorer's: H bound–free of n ≤ 30 with recombination (Milne) emission of the zone's own n_e n_p, H⁻,
#        H⁺ free–free, Thomson; no He opacity);
#    (b) SpecFull.FOpac pass-1 populations (adds He I / He II ground bound–free with the snapshot's He ion fractions, H levels
#        n ≤ 20 LTE relative to the zone's H⁺, Gaunt ff on H⁺/He⁺/He²⁺) on SpecFull.points (refined thick photosphere zones).
#    Spectrum.monochromatic (tangent rays, Feautrier ALI, coherent scattering): J_ν and H_ν at every point.
#    Q(r) = ∫ 16π² r² H_ν/(hν) dν and Γ_X = 4π ∫ σ_X J_ν/(hν) dν with the network's own cross sections (ChemRates.σ_verner,
#    Verner et al. 1996: H I 1s, He I 1s², He II 1s).
# 3. The network's own H ionization in each zone ahead of the photosphere, from the snapshot with the star-light rules of
#    PopIIIMaterial.starlight_setup (SLBLEND, SLEB, SLWFLUX as the run): three-level (1 − C)(n_e q12 + P12), collisional k1 n_e,
#    recombination C α_B n_e (case B, on the spot); and the equilibrium x_e with Γ_1s added.
using Serialization, Printf, LinearAlgebra
ENV["HECHEM"] = "1"
include(joinpath(@__DIR__, "snapspec.jl"))
isdefined(Main, :ModelAtom) || include(joinpath(@__DIR__, "ModelAtom.jl"))
isdefined(Main, :SpecFull) || include(joinpath(@__DIR__, "SpecFull.jl"))
include(joinpath(@__DIR__, "PrimordialChem.jl")); include(joinpath(@__DIR__, "ChemNet.jl")); include(joinpath(@__DIR__, "ChemRates.jl"))
const SS = SnapSpec; const SP = Spectrum; const CO = ContinuumOpacity; const SF = SpecFull; const CR = ChemRates; const CN = ChemNet
const hP = CO.h; const kB = CO.kB; const cl_ = CO.c; const eVe = 1.602176634e-12
const EH = 13.598; const EHe1 = 24.587; const EHe2 = 54.418
σH(ν) = CR.σ_verner(ν, EH, CR.VERNER_HI); σHe1(ν) = CR.σ_verner(ν, EHe1, CR.VERNER_HeI); σHe2(ν) = CR.σ_verner(ν, EHe2, CR.VERNER_HeII)

rundir = ARGS[1]; v = readstream(ARGS[2])
x = ARGS[3] == "last" ? v[end] : v[findfirst(y -> y.step == parse(Int, ARGS[3]), v)]
n = length(x.ρ); τ = SS.greytau(x); kph = findlast(>(2 / 3), τ); d = SS.densities(x)
Lrun = 4π * x.r[kph]^2 * x.fr[kph]
@printf("step %d  t %.6e yr  nz %d  kph %d (r_ph %.1f R☉)  T_eff(run) %.0f K  L_run(r_ph) %.3e L☉  L_run(edge) %.3e L☉\n",
        x.step, x.t / yr, n, kph, x.r[kph] / Rs, x.Teff, Lrun / Ls, 4π * x.r[end-2]^2 * x.fr[end-2] / Ls)
"Photons/s above E_th of a blackbody sphere (radius R, T) — Planck integral on a fine grid."
function Qbb(R, T, Eth)
    s = 0.0; νs = exp.(range(log(Eth * eVe / hP), log(Eth * eVe / hP) + log(1 + 60kB * T / (Eth * eVe)), length = 4001))
    for i in 1:length(νs)-1; ν = 0.5(νs[i] + νs[i+1]); s += π * CO.Bν(ν, T) / (hP * ν) * (νs[i+1] - νs[i]); end
    4π * R^2 * s
end
Teff_ph = (Lrun / (4π * x.r[kph]^2 * 5.670374e-5))^0.25
@printf("blackbody check at r_ph, T_eff %.0f K: Q(>13.6) %.3e  Q(>24.6) %.3e  Q(>54.4) %.3e s⁻¹\n", Teff_ph,
        Qbb(x.r[kph], Teff_ph, EH), Qbb(x.r[kph], Teff_ph, EHe1), Qbb(x.r[kph], Teff_ph, EHe2))

# ---------------- 1. the explorer's spectrum
f = joinpath(rundir, "spec", SS.key(x) * ".jls")
if isfile(f)
    s = deserialize(f); ν = SS.NU; q = SS.QB - 1 / SS.QB
    Qs(L, E) = sum((isfinite(L[j]) ? L[j] : 0.0) * q / hP for j in eachindex(ν) if hP * ν[j] * SS.QB >= E * eVe; init = 0.0)
    Ls_(L, E) = sum((isfinite(L[j]) ? L[j] : 0.0) * ν[j] * q for j in eachindex(ν) if hP * ν[j] * SS.QB >= E * eVe; init = 0.0)
    @printf("explorer spectrum (snapspec cache %s): grid ends at %.1f Å\n", basename(f), 1e8 * cl_ / ν[end])
    @printf("  Q(>13.6 eV): edge %.3e  r_ph %.3e s⁻¹;  L(>13.6 eV): edge %.3e  r_ph %.3e L☉;  Q(>24.6, last band only) edge %.3e r_ph %.3e\n",
            Qs(s.cont.Le, EH), Qs(s.cont.Lp, EH), Ls_(s.cont.Le, EH) / Ls, Ls_(s.cont.Lp, EH) / Ls, Qs(s.cont.Le, EHe1), Qs(s.cont.Lp, EHe1))
    println("  λ[Å]   νLν edge [L☉]  νLν r_ph [L☉]  r(τν=2/3) [R☉]")
    for j in eachindex(ν)
        λ = 1e8 * cl_ / ν[j]; λ < 1300 || continue
        @printf("  %7.1f  %11.3e  %11.3e  %11.3e\n", λ, s.cont.Le[j] * ν[j] / Ls, s.cont.Lp[j] * ν[j] / Ls, s.cont.rτ[j] / Rs)
    end
else
    println("no snapspec cache for this step: ", f)
end

# ---------------- 2. formal solution on an EUV grid
# (every point at least 0.05 % above an edge: SpecFull's model-atom thresholds (H 13.5984 eV) sit a hair above ChemRates' 13.598, and a
#  grid point exactly at 13.598 eV is transparent in FOpac: it leaked 1e47 photons/s to the edge in a first version of this script)
Egrid = sort(vcat(exp.(range(log(13.6057 * 1.0005), log(300.0), length = 120)), [EHe1 * 0.999, EHe1 * 1.0005, EHe2 * 0.999, EHe2 * 1.0005]))
νE = Egrid .* eVe ./ hP; nν = length(νE)
CO.tables()
function solve(atm, op, R)
    np = length(atm.r); J = zeros(np, nν); H = zeros(np, nν); κ = zeros(np, nν); χ = zeros(np, nν); B = zeros(np, nν); Le = zeros(nν)
    Threads.@threads for i in 1:nν
        L, Ji, _, m = SP.monochromatic(atm, op, R, νE[i]; maxit = 80, moments = true, maxit_direct = true)
        Le[i] = L; J[:, i] .= max.(Ji, 0.0); H[:, i] .= m.H; κ[:, i] .= m.κ; χ[:, i] .= m.χ; B[:, i] .= m.B
    end
    (; J, H, κ, χ, B, Le)
end
"∫ f(ν) dν over ν ≥ E_th on the grid (trapezoid)."
function ∫ν(fv, Eth)
    s = 0.0
    for i in 1:nν-1
        Egrid[i] >= Eth * 0.9999 || continue
        s += 0.5 * (fv[i] + fv[i+1]) * (νE[i+1] - νE[i])
    end
    s
end
kc = findlast(>(SS.TAUCUT), τ); kc = kc === nothing ? 1 : min(kc, n - 8)
rc = [0.5 * (x.r[k] + x.r[k+1]) for k in 1:n]
atmA = SP.atmosphere((r = x.r, d = x.ρ, T = x.T); rmin = kc > 1 ? rc[kc] : nothing, split = 3.0)
opA = SS.Opac(atmA, d); _, RA = SP.rays(atmA.r; ncore = kc > 1 ? 8 : 1); RA = SP.weights!(atmA, RA)
t0 = time(); solA = solve(atmA, opA, RA); @printf("\n(a) SnapSpec.Opac: %d points from zone %d, %.1f s\n", length(atmA.r), kc, time() - t0)
A = SF.points(x, d, τ); pts = A.P
atmB = (r = [p.r for p in pts], ρ = [p.ρ for p in pts], T = [p.T for p in pts], zone = [p.zone for p in pts], jump = A.jump)
_, RB = SP.rays(atmB.r; ncore = A.kc > 1 ? 8 : 1); RB = SP.weights!(atmB, RB)
opB = SF.FOpac(pts, nothing)
t0 = time(); solB = solve(atmB, opB, RB); @printf("(b) SpecFull.FOpac (pass-1 populations, He bound–free): %d points (%d refined zones), %.1f s\n", length(atmB.r), A.nref, time() - t0)

function photon_profile(lab, atm, sol, comp)
    np = length(atm.r)
    QH = [∫ν([16π^2 * atm.r[k]^2 * sol.H[k, i] / (hP * νE[i]) for i in 1:nν], EH) for k in 1:np]
    QHe1 = [∫ν([16π^2 * atm.r[k]^2 * sol.H[k, i] / (hP * νE[i]) for i in 1:nν], EHe1) for k in 1:np]
    QHe2 = [∫ν([16π^2 * atm.r[k]^2 * sol.H[k, i] / (hP * νE[i]) for i in 1:nν], EHe2) for k in 1:np]
    LEUV = [∫ν([16π^2 * atm.r[k]^2 * sol.H[k, i] for i in 1:nν], EH) for k in 1:np]
    ΓH = [∫ν([4π * σH(νE[i]) * sol.J[k, i] / (hP * νE[i]) for i in 1:nν], EH) for k in 1:np]
    ΓHe1 = [∫ν([4π * σHe1(νE[i]) * sol.J[k, i] / (hP * νE[i]) for i in 1:nν], EHe1) for k in 1:np]
    ΓHe2 = [∫ν([4π * σHe2(νE[i]) * sol.J[k, i] / (hP * νE[i]) for i in 1:nν], EHe2) for k in 1:np]
    # net H(1s) rate: 4π∫σ (J − S_th)/(hν): what the diffuse field does beyond the on-the-spot balance (S_th = the point's thermal
    # source, H 1s bound–free dominated above 13.6 eV in the H I zones; J = S_th there is case B's cancellation)
    ΓHnet = [∫ν([4π * σH(νE[i]) * (sol.J[k, i] - sol.B[k, i]) / (hP * νE[i]) for i in 1:nν], EH) for k in 1:np]
    # radial optical depth from the edge just above each edge: absorption only, and extinction
    iH = findfirst(E -> E > EH, Egrid); i1 = findfirst(E -> E > EHe1, Egrid); i2 = findfirst(E -> E > EHe2, Egrid)
    τa = zeros(np, 3); τx = zeros(np, 3); dτz = Dict{Int,Float64}()        # dτz: absorption at 912 Å inside each zone (its points' segments)
    for (c, i) in enumerate((iH, i1, i2)), k in np-1:-1:1
        dr = atm.r[k+1] - atm.r[k]
        τa[k, c] = τa[k+1, c] + SP.dtau(sol.κ[k, i], sol.κ[k+1, i], dr, atm.jump[k]); τx[k, c] = τx[k+1, c] + SP.dtau(sol.χ[k, i], sol.χ[k+1, i], dr, atm.jump[k])
    end
    # 912 Å absorption per zone: κ at the zone's centre point(s) × zone width (the run's zone, so independent of the point layout)
    for z in unique(atm.zone)
        ip = findall(==(z), atm.zone); dτz[z] = sum(sol.κ[ip, iH]) / length(ip) * (x.r[z+1] - x.r[z])
    end
    rph = x.r[kph]; τout = Dict{Int,Float64}(); acc = 0.0
    for z in sort(collect(keys(dτz))); z >= kph || continue; acc += dτz[z]; τout[z] = acc; end   # τ_912 from r_ph out to the zone's outer interface
    Le = sol.Le
    @printf("\n%s: at the edge Q(>13.6) %.3e  Q(>24.6) %.3e  Q(>54.4) %.3e s⁻¹, L(>13.6 eV) %.3e L☉\n", lab, ∫ν(Le ./ (hP .* νE), EH), ∫ν(Le ./ (hP .* νE), EHe1), ∫ν(Le ./ (hP .* νE), EHe2), ∫ν(Le, EH) / Ls)
    println("  pt zone  r[R☉]      T[K]     n_H     x_H0     x_e    τgrey | Q_H(r)    Q_HeI(r)  Q_HeII(r) L_EUV[L☉] | τabs912  τext912 τabs504 τabs228 | J912/B(T)  Γ_H1s     Γnet_1s   Γ_HeI     Γ_HeII  | Δτ912zone τ912(rph→)")
    for k in 1:np
        z = atm.zone[k]; i = iH
        @printf("  %3d %3d %10.4e %9.0f %9.2e %8.2e %6.3f %8.2e | %9.2e %9.2e %9.2e %9.2e | %8.2e %8.2e %8.2e %8.2e | %9.2e %9.2e %10.2e %9.2e %9.2e | %9.2e %9.2e\n",
                k, z, atm.r[k] / Rs, atm.T[k], comp.nH[k], comp.xH0[k], comp.xe[k], 0.5(τ[z] + τ[z+1]), QH[k], QHe1[k], QHe2[k], LEUV[k] / Ls,
                τa[k, 1], τx[k, 1], τa[k, 2], τa[k, 3], sol.J[k, i] / max(CO.Bν(νE[i], atm.T[k]), 1e-300), ΓH[k], ΓHnet[k], ΓHe1[k], ΓHe2[k], dτz[z], get(τout, z, NaN))
    end
    (; QH, QHe1, QHe2, LEUV, ΓH, ΓHnet, ΓHe1, ΓHe2, τa, τx, dτz, τout)
end
μH = SS.μH
compA = (nH = x.ρ[atmA.zone] ./ μH, xH0 = d.nH[atmA.zone] ./ (x.ρ[atmA.zone] ./ μH), xe = d.ne[atmA.zone] ./ (x.ρ[atmA.zone] ./ μH))
compB = (nH = [p.ρ / μH for p in pts], xH0 = [p.nH0 / (p.ρ / μH) for p in pts], xe = [p.ne / (p.ρ / μH) for p in pts])
PA = photon_profile("(a) explorer opacity", atmA, solA, compA)
PB = photon_profile("(b) SpecFull opacity with He", atmB, solB, compB)

# ---------------- 3. the network's own ionization in the zones ahead of the photosphere
const arad_ = 7.565723e-15; const σSB = 5.670374419e-5; const Tcmb = 54.5
const LYA_LAM = 1.215668e-5; const LAM2S = 8.2245809; const A21LYA = 6.2648e8                  # PopIIIMaterial
const E21 = 10.19881 * eVe; const χH_ = CN.χH; const CHI2 = χH_ - E21; const CR_REC = 1.799920e10
alphaB(T) = (t = T / 1e4; 4.309e-13 * t^-0.6166 / (1 + 0.6703 * t^0.53))
beta2s_bb(Tc) = alphaB(Tc) * (CR_REC * Tc)^1.5 * exp(-CHI2 / (kB * Tc))
q12_H(T) = 7.5e-19 * exp(-118348.0 / T) / (1 + sqrt(T / 1e5)) / E21
q21_H(T) = 7.5e-19 / (1 + sqrt(T / 1e5)) / E21 / 4
er = x.er; fr = x.fr; un = x.u; rn = x.r
const kphs = findlast(>=(2 / 3), τ[1:n])                              # starlight_setup's kph
Lst = 4π * rn[kphs+1]^2 * fr[kphs+1]; Te = (Lst / (4π * rn[kphs+1]^2 * σSB))^0.25; Ecmb = arad_ * Tcmb^4
function tr_at_tau(t)
    t >= τ[1] && return (max(er[1], Ecmb) / arad_)^0.25
    for j in 1:n-1
        if τ[j] >= t > τ[j+1]
            t1, t2 = log(max(τ[j], 1e-30)), log(max(τ[j+1], 1e-30)); w = t1 > t2 ? (log(t) - t2) / (t1 - t2) : 0.0
            l1, l2 = log(max(er[j], Ecmb, 1e-300)), log(max(er[j+1], Ecmb, 1e-300)); return (exp(l2 + w * (l1 - l2)) / arad_)^0.25
        end
    end
    (max(er[n], Ecmb) / arad_)^0.25
end
slcu(k) = max(cl_ * (er[k] - Ecmb), 0.5 * (max(fr[k], 0.0) + max(fr[k+1], 0.0)))
iy(s) = findfirst(==(s), x.spn)
photos = CR.modern(; li = true, he = true)[2]; pHe = photos[findfirst(p -> p.R === :He, photos)]; pHep = photos[findfirst(p -> p.R === :Hep, photos)]
envcmb = (Trad = Tcmb, J = ν -> CO.Bν(ν, Tcmb), static = true)
"Zone k of the network: three-level parameters and the H, He rates; Γs from the formal solution (points of zone k, mean)."
function zone_rates(k, Γ1s, ΓHe1)
    nH = x.ρ[k] / μH; T = x.T[k]; xH = x.Y[iy(:H), k]; xe = x.Y[iy(:e), k]; xp = x.Y[iy(:Hp), k]; xHe = x.Y[iy(:He), k]; xHep = x.Y[iy(:Hep), k]
    Tr = (max(er[k], Ecmb) / arad_)^0.25; Tsrc = tr_at_tau(τ[k] + 2 / 3)
    wq = clamp((log10(max(τ[k], 1e-30)) + log10(1 / 0.3)) / 1.0, 0.0, 1.0); wq = wq * wq * (3 - 2wq)
    Tc = Tsrc^(1 - wq) * Tr^wq; W = wq >= 1 ? 1.0 : min(slcu(k) / (4σSB * Tc^4), 1.0)
    n1 = nH * max(xH, 1e-30); g = (abs(un[k+1] - un[k]) / (rn[k+1] - rn[k]) + 2abs(0.5(un[k] + un[k+1])) / rc[k]) / 3
    τS = 3A21LYA * LYA_LAM^3 * n1 / (8π * max(g, 1e-30)); RL = 3A21LYA * (τS < 1e-6 ? 1.0 : -expm1(-τS) / τS)
    nb = W / expm1(E21 / (kB * Tc)); β = W * beta2s_bb(Tc); ne = nH * xe; P12 = (RL + LAM2S) * nb
    Cf(ne_) = (D = RL + LAM2S + 4ne_ * q21_H(T); D / (D + β))
    C = Cf(ne)
    I3 = (1 - C) * (ne * q12_H(T) + P12); Icoll = PrimordialChem.k1_aazn(T) * ne; rec = C * alphaB(T) * ne
    # equilibrium x_p with Γ_1s added (n_e = n_H (x_p + He electrons), C and the n_e terms recomputed), from the snapshot's x_p
    xHe_e = xe - xp
    function xeq(Γ)
        y = xp
        for _ in 1:200
            ne_ = nH * (y + xHe_e); C_ = Cf(ne_)
            Rr = max((1 - C_) * (ne_ * q12_H(T) + P12) + PrimordialChem.k1_aazn(T) * ne_ + Γ, 0.0) / (C_ * alphaB(T) * ne_)
            yn = Rr / (1 + Rr) * (1 - 2x.Y[iy(:H2), k]); abs(yn - y) < 1e-10 * max(y, 1e-30) && (y = yn; break); y = 0.5(y + yn)
        end
        y
    end
    αHe = CN.photorates(pHe, T, envcmb)[2]                                        # the network's He I recombination (Milne, all levels)
    (; k, r = rc[k] / Rs, T, nH, xH, xe, xp, Tc, W, τS, C, β, P12, I3, Icoll, rec, Γ1s, ratio = Γ1s / (I3 + Icoll),
       xp0 = xeq(0.0), xp1 = xeq(Γ1s), xHe, xHep, αHe, ΓHe1, xHep_eq = ΓHe1 / (ΓHe1 + αHe * ne) * yHe_)
end
const yHe_ = 1 / 12
function network_table(lab, atm, P)
    println("\n$lab: network zones (eq = 0) ahead of the photosphere; rates per neutral atom [s⁻¹]. T_eff(starlight) $(round(Int, Te)) K, kph $kphs")
    println("  zone r[R☉]     T     n_H      x_H     x_e    | T_c     W        C        | I_3lev    I_coll    Γnet_1s   Γnet/I    Γgross_1s      C α_B n_e | x_p snap  x_p eq(0)  x_p eq(Γ) Δx_e     | x_He+     Γ_HeI     α_He n_e  x_He+ eq(Γ)")
    rows = []
    for k in kphs+1:n
        x.eq[k] && continue
        ip = findall(==(k), atm.zone); isempty(ip) && continue
        Γ = sum(P.ΓHnet[ip]) / length(ip); Γg = sum(P.ΓH[ip]) / length(ip); Γ1 = sum(P.ΓHe1[ip]) / length(ip)
        z = zone_rates(k, Γ, Γ1); push!(rows, z)
        @printf("  %3d %9.3e %6.0f %8.2e %8.2e %7.4f | %6.0f %8.2e %8.2e | %9.2e %9.2e %10.2e %10.2e %9.2e %9.2e | %9.3e %9.3e %9.3e %9.2e | %9.2e %9.2e %9.2e %9.2e\n",
                z.k, z.r, z.T, z.nH, z.xH, z.xe, z.Tc, z.W, z.C, z.I3, z.Icoll, z.Γ1s, z.ratio, Γg, z.rec, z.xp, z.xp0, z.xp1, z.xp1 - z.xp0,
                z.xHep, z.ΓHe1, z.αHe * z.nH * z.xe, z.xHep_eq)
    end
    rows
end
NA = network_table("(a)", atmA, PA); NB = network_table("(b)", atmB, PB)
length(ARGS) >= 4 && serialize(ARGS[4], (; step = x.step, Egrid, atmA, atmB, solA, solB, PA, PB, NA, NB))

# ---------------- 4. Strömgren budget at the H I front: photons arriving vs case-B recombinations per zone
"Q(r) of a solution interpolated (linear in r) to radius r0."
function Qat(atm, P, r0)
    i = searchsortedlast(atm.r, r0); (1 <= i < length(atm.r)) || return (NaN, NaN, NaN)
    w = (r0 - atm.r[i]) / (atm.r[i+1] - atm.r[i]); f(v) = (1 - w) * v[i] + w * v[i+1]
    (f(P.QH), f(P.QHe1), f(P.QHe2))
end
println("\nQ through spheres [s⁻¹] (H, He I, He II), (a) explorer opacity | (b) SpecFull opacity:")
knet = findfirst(k -> !x.eq[k] && k > kph, 1:n)
for (lab, r0) in (("r_ph (kph inner interface)", x.r[kph]), ("kph+1 interface", x.r[kph+1]), ("first network zone inner edge", x.r[knet]))
    qa = Qat(atmA, PA, r0); qb = Qat(atmB, PB, r0)
    @printf("  %-30s r %8.1f R☉: %9.2e %9.2e %9.2e | %9.2e %9.2e %9.2e\n", lab, r0 / Rs, qa..., qb...)
end
@printf("  max Q_H inside: (a) %.3e at r = %.1f R☉;  (b) %.3e at %.1f R☉\n", maximum(PA.QH), atmA.r[argmax(PA.QH)] / Rs, maximum(PB.QH), atmB.r[argmax(PB.QH)] / Rs)
println("\nzone  r_in–r_out [R☉]       T     eq  x_H        α_B n_e n_p V [s⁻¹]  Q_H(r_in) (a)  J912/B(T) (a)")
for k in kph-4:min(kph+5, n)
    V = 4π / 3 * (x.r[k+1]^3 - x.r[k]^3); nH = x.ρ[k] / μH
    rec = alphaB(x.T[k]) * nH^2 * x.Y[iy(:e), k] * x.Y[iy(:Hp), k] * V
    ip = findall(==(k), atmA.zone); jb = isempty(ip) ? NaN : sum(solA.J[ip, 1] ./ CO.Bν.(νE[1], atmA.T[ip])) / length(ip)
    @printf("  %3d %8.1f–%8.1f %7.0f %d %9.2e   %9.2e            %9.2e      %9.2e\n", k, x.r[k] / Rs, x.r[k+1] / Rs, x.T[k], x.eq[k], x.Y[iy(:H), k], rec, Qat(atmA, PA, x.r[k])[1], jb)
end
