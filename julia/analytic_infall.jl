# The simple model of the infall's hydrogen: no radiative transfer. Each point of a steady free fall at Ṁ sees the
# star's light, diluted geometrically and nothing else:  J_ν(r) = W(r) I_ν,  W = ½[1 − √(1 − (R/r)²)],
# with I_ν the star's surface intensity (a blackbody at T_eff, or the emergent spectrum of the TITAN star without
# infall, star_t0, scaled to luminosity L). The statistical equilibrium of NLTEH.solve_point (12 levels, Sobolev lines
# with this J as the local continuum, the same collision rates) gives the populations; no X-rays, no diffuse field,
# no attenuation, helium neutral. Temperatures: T⁴ = T_eff⁴ W (dilute radiative equilibrium) or, for the direct
# comparison, the temperatures of a non-LTE solution (temp=NAME: lib/spectra/nlte_NAME.jls, which also sets the radii).
#   Z0_TABLE=ext TITAN_NOLIB=1 julia --project=env analytic_infall.jl [mdot=3] [teff=7310] [L=1.4204e44] [rmax=120]
#        [star=bb|t0] [temp=dilute|NAME] [out=NAME]
# Prints the ionized fraction, b₂, the Hα source function and Sobolev depth along the infall; the neutral column, the
# Thomson depth, and the net line power (Σ hν β [n_u A − (n_l B_lu − n_u B_ul) J]) of Hα, Hβ, Paα, Lyα and all lines.
# Writes lib/spectra/analytic_OUT.jls.
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf
const SP = Spectrum; const Z0 = Z0Physics; const NH = NLTEH
opt = Dict(Symbol(split(a, "=")[1]) => split(a, "=")[2] for a in ARGS)
g(k, d) = parse(Float64, get(opt, k, d))
Mdot = g(:mdot, "3") * Z0.Msun / Z0.yr; Teff = g(:teff, "7310"); L = g(:L, "1.4204e44"); rmax = g(:rmax, "120")
M = 1e6 * Z0.Msun; R0 = 8.4352e15
Rst = sqrt(L / (4π * Z0.σSB * Teff^4))                    # photospheric radius for this L and T_eff
X = 0.76; mH = Z0.mH
ν = SP.grid(; λmax_μm = 30.0, n = 260, nedge = 12)
wν = [i == 1 ? 0.5(ν[2] - ν[1]) : i == length(ν) ? 0.5(ν[end] - ν[end-1]) : 0.5(ν[i+1] - ν[i-1]) for i in eachindex(ν)]
star = get(opt, :star, "bb")
Iν = if star == "bb"
    [SP.Bν(x, Teff) for x in ν]
else                                                      # emergent spectrum of the star alone, as a surface intensity
    d0 = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "star_t0.jls"))
    Lt = [exp(SP_interp) for SP_interp in (begin
        i = clamp(searchsortedlast(d0.ν, x), 1, length(d0.ν) - 1); t = (log(x) - log(d0.ν[i])) / (log(d0.ν[i+1]) - log(d0.ν[i]))
        log(max(d0.L[i], 1e-300)) * (1 - t) + log(max(d0.L[i+1], 1e-300)) * t end for x in ν)]
    Ltot = sum(0.5 * (Lt[i] + Lt[i+1]) * (ν[i+1] - ν[i]) for i in 1:length(ν)-1)
    Lt .* (L / Ltot) ./ (4π^2 * Rst^2)
end
# fuv=f: multiply the star's intensity below 1300 Å (the far UV that pumps the Lyman lines) by f, a sensitivity test
# of the one input everything rests on (the TITAN star's far-UV excess over a blackbody, from a photosphere without lines)
fuv = g(:fuv, "1"); Iν = [Z0.c / ν[q] * 1e8 < 1300 ? fuv * Iν[q] : Iν[q] for q in eachindex(ν)]
W(r) = 0.5 * (1 - sqrt(max(1 - (Rst / r)^2, 0.0)))
tmode = get(opt, :temp, "dilute")
rs, Ts = if tmode == "dilute"
    r = exp.(range(log(1.0015Rst), log(rmax * R0), length = 160)); (r, [Teff * W(x)^0.25 for x in r])
else
    dn = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$tmode.jls")); (dn.r, dn.T)
end
vff(r) = sqrt(2Z0.G * M / r)
# trap=1: the ionized infall scatters the star's light (Thomson, coherent), which raises the mean intensity above the
# geometric value. Spherical Eddington approximation (K = J/3) in a pure scattering envelope carrying the star's
# flux, dJ_ν/dr = −3 n_e σ_T H_ν with H_ν = L_ν/(16π² r²), integrated inwards from J → W I_ν at large r:
#   J_ν(r) = W(r) I_ν + 3 H_ν(r) τ*(r),   τ*(r) = r² ∫_r^∞ n_e σ_T r'⁻² dr',   H_ν = I_ν R*²/(4r²).
# For n_e ∝ r^(-3/2) (an ionized free fall), τ* = τ_es/5, so J/J_geom ≈ 1 + 0.3 τ_es near the star. τ* depends on the
# ionization, so the whole profile is iterated (from the geometric solution) until it settles.
trap = get(opt, :trap, "0") == "1"; verbose = get(opt, :verbose, "0") == "1"
σT = 6.6524e-25
# caseB=1: the Lyman continuum in detailed balance (no net photoionization from n = 1), as the full solution finds
# (its Lyman-continuum J is recombination radiation reabsorbed on the spot; nlte_rates.jl).
# atten=1: the starlight is absorbed on its way out by bound–free absorption from the excited levels,
#   J_ν → J_ν exp(−τ_ν(r)),  τ_ν(r) = ∫_{R*}^r Σ_{i≥2} n_i σ_i(ν) dr'  (radial, from the photosphere),
# accumulated in the outward sweep from the populations already solved further in.
# ots2=1: on-the-spot Balmer continuum, the n = 2 analogue of case B: recombinations to n = 2 emit Balmer-continuum
# photons that are absorbed nearby by another n = 2 atom, so neither counts (J holds only the attenuated starlight).
Jfull = haskey(opt, :jfull) ? deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlterates_$(opt[:jfull]).jls")) : nothing
band = haskey(opt, :band) ? parse.(Float64, split(opt[:band], ",")) : [0.0, 0.0]
Jfull === nothing || (length(Jfull.r) == length(rs) && all(abs.(Jfull.r ./ rs .- 1) .< 1e-9)) || error("jfull: radii differ from temp=")
caseB = get(opt, :caseB, "0") == "1"; atten = get(opt, :atten, "0") == "1"; ots2 = get(opt, :ots2, "0") == "1"; albedo = get(opt, :albedo, "0") == "1"
function sweep(τes)
    out = []; τν = zeros(length(ν))
    for (k, (r, T)) in enumerate(zip(rs, Ts))
        ρ = Mdot / (4π * r^2 * vff(r)); nH = X * ρ / mH
        if atten && k > 1
            p = out[end]; dr = r - p.r
            # albedo=1: a photoionization from n = 2 is followed by a recombination to n = 2 (re-emitting a Balmer-
            # continuum photon) with probability a = α₂/Σ_{i≥2} α_i, so only the fraction 1 − a of the n = 2 bound–free
            # opacity removes Balmer-continuum light (the rest re-emits it, approximated as forward scattering)
            fa = 1.0
            if albedo
                α = [NH.Φ(i, p.T) * sum(4π * wν[q] * NH.σbf(i, ν[q]) * 2Z0.h * ν[q]^3 / Z0.c^2 * exp(-Z0.h * ν[q] / (Z0.kB * p.T)) / (Z0.h * ν[q])
                     for q in eachindex(ν) if Z0.h * ν[q] < 100Z0.eV) for i in 2:NH.NL]
                fa = 1 - α[1] / sum(α)
            end
            for q in eachindex(ν)
                κ = fa * p.n[2] * NH.σbf(2, ν[q]) + sum(p.n[i] * NH.σbf(i, ν[q]) for i in 3:NH.NL)
                τν[q] += κ * dr
            end
        end
        J = (W(r) .* Iν .+ (trap ? 3 * τes[k] * Rst^2 / (4r^2) : 0.0) .* Iν) .* exp.(-τν)       # τes holds τ*
        # jfull=NAME (with temp=NAME): use the full solution's own J_ν (lib/spectra/nlterates_NAME.jls), except in the
        # wavelength band band=λ1,λ2 [Å], where the simple-model J above is kept: which part of the radiation field matters
        if Jfull !== nothing
            Jf = Jfull.J[k, :]
            J = [band[1] <= Z0.c / ν[q] * 1e8 <= band[2] ? J[q] : Jf[q] for q in eachindex(ν)]
        end
        n0 = isempty(out) ? fill(1e-6nH, NH.NL) : out[end].n ./ out[end].nH .* nH
        np0 = isempty(out) ? nH : out[end].np / out[end].nH * nH
        nn, npp, ne, β, P = NH.solve_point(T, nH, 0.0, vff(r), r, -0.5, ν, wν, J, 0.0, n0, np0; caseB, ots2)
        push!(out, (r = r, T = T, nH = nH, n = nn, np = npp, ne = ne, β = β, P = P, J = J,
                    att = [exp(-τν[clamp(searchsortedlast(ν, Z0.c / (λ * 1e-8)), 1, length(ν))]) for λ in (3000.0, 6000.0)]))
    end
    out
end
function τprofile(out; star = true)        # star: τ* (weighted by (r/r')²); else the plain Thomson depth
    τ = zeros(length(out)); I = 0.0
    for k in length(out)-1:-1:1
        I += 0.5 * (out[k].ne / (star ? out[k].r^2 : 1.0) + out[k+1].ne / (star ? out[k+1].r^2 : 1.0)) * σT * (out[k+1].r - out[k].r)
        τ[k] = (star ? out[k].r^2 : 1.0) * I
    end
    τ
end
τes = zeros(length(rs))
# seed=NAME: start the trapping iteration from the Thomson-depth profile of an earlier result (analytic_NAME.jls), to
# test whether the answer depends on where the iteration starts
if haskey(opt, :seed)
    s0 = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "analytic_$(opt[:seed]).jls"))
    xs = [s0.xH[clamp(searchsortedlast(s0.r, r), 1, length(s0.r))] for r in rs]
    ne0 = [x * X * Mdot / (4π * r^2 * vff(r)) / mH for (x, r) in zip(xs, rs)]
    acc = 0.0
    for k in length(rs)-1:-1:1; acc += 0.5 * (ne0[k] / rs[k]^2 + ne0[k+1] / rs[k+1]^2) * σT * (rs[k+1] - rs[k]); τes[k] = rs[k]^2 * acc; end
    @printf("seeded from analytic_%s: τ* at the base %.3f\n", opt[:seed], τes[1])
end
sw = sweep(τes)
if trap
    for it in 1:40
        global τes, sw
        τn = τprofile(sw)
        d = maximum(abs.(τn .- τes)) ; τes = 0.5 .* (τes .+ τn)          # damped
        sw = sweep(τes)
        @printf("trap iteration %2d: Thomson depth %.3f, τ* at the base %.3f (change %.3e)\n", it, τprofile(sw; star = false)[1], τprofile(sw)[1], d)
        d < 1e-3 && break
    end
end
res = []
for (k, q) in enumerate(sw)
    r, T, nH, nn, npp, ne, β, P, J = q.r, q.T, q.nH, q.n, q.np, q.ne, q.β, q.P, q.J
    verbose && k % 10 == 0 && @printf("   r/R0 %7.2f attenuation at 3000 Å %.3e, 6000 Å %.3e; J(3000)/[W B] %.3f\n", r / R0, q.att..., J[clamp(searchsortedlast(ν, Z0.c / 3000e-8), 1, length(ν))] / (W(r) * SP.Bν(Z0.c / 3000e-8, Teff)))
    # net power of selected lines
    Pl = Dict{Tuple{Int,Int},Float64}()
    for (l, u) in ((2, 3), (2, 4), (3, 4), (1, 2))
        νl = (NH.En[u] - NH.En[l]) / Z0.h; f = NH.fosc(l, u)
        Blu = 4π * NH.πe2mc * f / (Z0.h * νl); Bul = NH.gn[l] / NH.gn[u] * Blu; A = 2Z0.h * νl^3 / Z0.c^2 * Bul
        i = clamp(searchsortedlast(ν, νl), 1, length(ν) - 1); Jl = J[i] + (J[i+1] - J[i]) * (νl - ν[i]) / (ν[i+1] - ν[i])
        Pl[(l, u)] = Z0.h * νl * β[l, u] * (nn[u] * A - (nn[l] * Blu - nn[u] * Bul) * Jl)
    end
    ν32 = (NH.En[3] - NH.En[2]) / Z0.h
    SL = 2Z0.h * ν32^3 / Z0.c^2 / (nn[2] * NH.gn[3] / (nn[3] * NH.gn[2]) - 1)
    i = clamp(searchsortedlast(ν, ν32), 1, length(ν) - 1); Jc = J[i] + (J[i+1] - J[i]) * (ν32 - ν[i]) / (ν[i+1] - ν[i])
    K = NH.πe2mc * NH.fosc(2, 3) * Z0.c / ν32 * (nn[2] - NH.gn[2] / NH.gn[3] * nn[3])
    push!(res, (r = r, T = T, nH = nH, n = nn, np = npp, ne = ne, P = P, Pl = Pl, SJ = SL / Jc, τ = K / (vff(r) / r),
                b2 = nn[2] / (npp * ne * NH.Φ(2, T))))
end
# radial integrals (trapezoid in r)
I(f) = sum(0.5 * (f(res[i]) + f(res[i+1])) * (res[i+1].r - res[i].r) for i in 1:length(res)-1)
NHI = I(q -> q.nH - q.np); τT = I(q -> q.ne * 6.6524e-25)
Lof(f) = I(q -> 4π * q.r^2 * f(q))
Ltot = Lof(q -> q.P); LHa = Lof(q -> q.Pl[(2, 3)]); LHb = Lof(q -> q.Pl[(2, 4)]); LPa = Lof(q -> q.Pl[(3, 4)]); LLy = Lof(q -> q.Pl[(1, 2)])
@printf("Ṁ = %.2f Msun/yr, T_eff = %.0f K, L = %.3e, R* = %.4f R0, star = %s, T = %s, r to %.0f R0\n", Mdot * Z0.yr / Z0.Msun, Teff, L, Rst / R0, star, tmode, rmax)
@printf("%8s %6s %9s %9s %9s %9s\n", "r/R0", "T", "x_H", "b2", "S/J(Hα)", "τ(Hα)")
for (k, q) in enumerate(res)
    (k % max(1, length(res) ÷ 20) == 0 || k == 1) || continue
    @printf("%8.3f %6.0f %9.4f %9.3e %9.3e %9.3e\n", q.r / R0, q.T, q.np / q.nH, q.b2, q.SJ, q.τ)
end
@printf("neutral column %.3e cm⁻², Thomson depth %.3f\n", NHI, τT)
@printf("net line power: all %.3e, Hα %.3e, Hβ %.3e, Paα %.3e, Lyα %.3e erg/s; Hα/Hβ %.2f\n", Ltot, LHa, LHb, LPa, LLy, LHa / LHb)
for x in (5.0, 12.0, 30.0)
    @printf("  fraction of the line power beyond %4.0f R0: %.2f\n", x, Lof(q -> q.r > x * R0 ? q.P : 0.0) / Ltot)
end
out = get(opt, :out, "")
isempty(out) || serialize(joinpath(@__DIR__, "..", "lib", "spectra", "analytic_$out.jls"),
    (opt = opt, r = [q.r for q in res], T = [q.T for q in res], xH = [q.np / q.nH for q in res], b2 = [q.b2 for q in res],
     SJ = [q.SJ for q in res], τ = [q.τ for q in res], P = [q.P for q in res], PHa = [q.Pl[(2, 3)] for q in res],
     NHI = NHI, τT = τT, Ltot = Ltot, LHa = LHa, LHb = LHb, LPa = LPa, LLy = LLy))
