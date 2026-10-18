# Emergent hydrogen line spectrum of a non-LTE solution (nlte_run.jl): profiles of the Lyman, Balmer, Paschen and
# Brackett lines between levels ≤ 12, computed as in line_profile.jl (ray integration through the infall with Gaussian
# profiles, thermal + ξ microturbulence) with the non-LTE populations and source functions and the non-LTE continuum,
# then put on the continuum L_ν of the same solution to give the full spectrum.
#   Z0_TABLE=ext TITAN_NOLIB=1 julia -t 2 --project=env nlte_spectrum.jl NAME [OUT]
# NAME: spectrum of spectra_run.jl (lib/spectra/NAME.jls); OUT: non-LTE solution lib/spectra/nlte_OUT.jls (default NAME).
# Writes lib/spectra/linespec_OUT.jls: per line (series, lower, upper, λ_vac, v, F/F_c, EW, L_line) and the composite
# spectrum on a wavelength grid. Wavelengths are vacuum (Rydberg constant for hydrogen). Lines are computed
# independently and multiplied onto the continuum; overlaps (near series limits) are not treated.
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf
const SP = Spectrum; const Z0 = Z0Physics; const NH = NLTEH
name = ARGS[1]; on = length(ARGS) >= 2 ? ARGS[2] : name
d = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "$name.jls"))
dn = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$on.jls"))
S = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(d.run)_relax_snaps.jls"))
s = S[argmin([abs(x.t - d.t) for x in S])]
atm = haskey(dn, :Tstar) ? merge(d.atm, (T = dn.Tstar,)) : d.atm
M = 1e6 * Z0.Msun; R0 = 8.4352e15
rsh = s.r[d.ksh+1]
inf = dn.sel; rI = atm.r[inf]; TI = atm.T[inf]; rmax = rI[end]; uabs = dn.vabs
function interp(x, xs, ys)
    x <= xs[1] && return ys[1]; x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x); t = (log(x) - log(xs[i])) / (log(xs[i+1]) - log(xs[i]))
    ys[i] + t * (ys[i+1] - ys[i])
end
const ξturb = parse(Float64, get(ENV, "LP_XI", "10")) * 1e5
nlfull = Any[nothing for _ in eachindex(atm.r)]; for (q, j) in enumerate(inf); nlfull[j] = dn.nl[q]; end
op = NH.NOpac(SP.Opac(atm), nlfull)
_, R = SP.rays(atm.r; ncore = 12); R = SP.weights!(atm, R)
RH = 109677.58                                                    # Rydberg constant for H [cm⁻¹]
SERIES = Dict(1 => "Ly", 2 => "H", 3 => "Pa", 4 => "Br")
const GREEK = collect("αβγδεζηθικλμ")
LINES = [(l, u) for l in 1:4 for u in l+1:(l == 1 ? 6 : l == 2 ? 10 : l == 3 ? 9 : 8)]
vgrid = collect(range(-3000e5, 3000e5, length = 301))
n = 3000; gq = [(exp(5 * (i - 1) / (n - 1)) - 1) / (exp(5) - 1) for i in 1:n]
# rays: samples along z, far side to near side (far side only if the ray misses the star), crowded towards the star
RAYZ = map(1:length(R)) do j
    p = R[j].p
    znear = sqrt(max(rmax^2 - p^2, 0.0)); zfar = p < rsh ? sqrt(rsh^2 - p^2) : -znear
    (znear > zfar && p < rmax) || return nothing
    zs = p < rsh ? zfar .+ (znear - zfar) .* gq : (zp = znear .* gq; vcat(-reverse(zp), zp[2:end]))
    rz = sqrt.(p^2 .+ zs .^ 2)
    (z = zs, r = rz, w = [r < rsh ? NaN : interp(r, rI, uabs) * z / r for (z, r) in zip(zs, rz)],
     Δv = [sqrt(2Z0.kB * interp(r, rI, TI) / Z0.mH + ξturb^2) for r in rz])
end
out = []
for (l, u) in LINES
    λ = 1e8 / (RH * (1 / l^2 - 1 / u^2)); ν0 = Z0.c / (λ * 1e-8)
    label = SERIES[l] * string(GREEK[u-l])
    gl, gu = NH.gn[l], NH.gn[u]; f = NH.fosc(l, u)
    _, _, _, ci = SP.monochromatic(atm, op, R, ν0; moments = :I)
    ps, Ic = ci.p, ci.I
    Fc = sum(0.5 * (Ic[j] * ps[j] + Ic[j+1] * ps[j+1]) * (ps[j+1] - ps[j]) for j in 1:length(ps)-1)
    # inverted points (K < 0: only the few points right above the shock where the X-ray upper bound empties n = 1, see
    # nlte_run.jl) are treated as transparent in this line instead of amplifying it
    Kr = [max(NH.πe2mc * f * λ * 1e-8 * (st.n[l] - gl / gu * st.n[u]), 0.0) for st in dn.nl]
    SLr = [K > 0 ? 2Z0.h * ν0^3 / Z0.c^2 / (st.n[l] * gu / (st.n[u] * gl) - 1) : 0.0 for (K, st) in zip(Kr, dn.nl)]
    ninv = count(==(0.0), Kr); ninv > 0 && println("  $label: $ninv inverted point(s) set transparent")
    KZ = [ry === nothing ? nothing : [r < rI[1] ? 0.0 : interp(r, rI, Kr) for r in ry.r] for ry in RAYZ]
    SZ = [ry === nothing ? nothing : [interp(r, rI, SLr) for r in ry.r] for ry in RAYZ]
    F = zeros(length(vgrid))
    Threads.@threads for iv in eachindex(vgrid)
        v = vgrid[iv]; Iv = similar(Ic)
        for j in eachindex(ps)
            I = Ic[j]; ry = RAYZ[j]
            if ry !== nothing
                zs, wz, Δvz, Kz, Sz = ry.z, ry.w, ry.Δv, KZ[j], SZ[j]
                for i in 1:length(zs)-1
                    (isfinite(wz[i]) && isfinite(wz[i+1])) || continue
                    lo, hi = minmax(wz[i], wz[i+1]); dv = max(Δvz[i], Δvz[i+1])
                    (v < lo - 5dv || v > hi + 5dv) && continue
                    m = max(1, ceil(Int, (hi - lo) / (0.25 * min(Δvz[i], Δvz[i+1]))))
                    dz = (zs[i+1] - zs[i]) / m
                    for q in 1:m
                        t = (q - 0.5) / m
                        w = wz[i] + t * (wz[i+1] - wz[i]); D = Δvz[i] + t * (Δvz[i+1] - Δvz[i])
                        x = (v - w) / D; abs(x) > 6 && continue
                        τ = (Kz[i] + t * (Kz[i+1] - Kz[i])) * exp(-x^2) / (sqrt(π) * D) * dz
                        S_ = Sz[i] + t * (Sz[i+1] - Sz[i])
                        I = I * exp(-τ) + S_ * (-expm1(-τ))
                    end
                end
            end
            Iv[j] = I
        end
        F[iv] = sum(0.5 * (Iv[j] * ps[j] + Iv[j+1] * ps[j+1]) * (ps[j+1] - ps[j]) for j in 1:length(ps)-1) / Fc
    end
    EW = sum(0.5 * ((1 - F[i]) + (1 - F[i+1])) * (vgrid[i+1] - vgrid[i]) for i in 1:length(vgrid)-1) / Z0.c * λ   # Å, >0 absorption
    Lc_λ = 8π^2 * Fc * ν0 / λ            # continuum L_λ [erg s⁻¹ Å⁻¹] at the line (L_ν = 8π² ∫ I p dp)
    Lline = -EW * Lc_λ
    k = argmax(F)
    flush(stdout)
    @printf("%-5s %2d–%2d λ %9.2f Å: peak %8.3f at %+5.0f km/s, min %.3f, EW %+9.2f Å, L_line %+.3e erg/s\n",
            label, l, u, λ, F[k], vgrid[k] / 1e5, minimum(F), EW, Lline)
    push!(out, (name = label, l = l, u = u, λ = λ, v = vgrid ./ 1e5, F = F, EW = EW, Lline = Lline, Lcλ = Lc_λ))
end
# composite spectrum: continuum L_ν of the non-LTE solution (its own grid), times every line profile, on a λ grid
λgrid = sort(unique(vcat(exp.(range(log(900.0), log(30000.0), length = 4000)),
                         [x.λ * (1 + v / 2.99792458e5) for x in out for v in range(-3000, 3000, length = 121)])))
νc = dn.ν; Lc = dn.L
Lν = map(λgrid) do λ
    ν = Z0.c / (λ * 1e-8); i = clamp(searchsortedlast(νc, ν), 1, length(νc) - 1)
    t = (log(ν) - log(νc[i])) / (log(νc[i+1]) - log(νc[i]))
    a, b = max(Lc[i], 1e-300), max(Lc[i+1], 1e-300); exp(log(a) * (1 - t) + log(b) * t)
end
prof = ones(length(λgrid))
for x in out
    for (i, λ) in enumerate(λgrid)
        v = (λ / x.λ - 1) * 2.99792458e5
        abs(v) < 3000 || continue
        j = clamp(searchsortedlast(x.v, v), 1, length(x.v) - 1); t = (v - x.v[j]) / (x.v[j+1] - x.v[j])
        prof[i] *= x.F[j] * (1 - t) + x.F[j+1] * t
    end
end
νLν = [Z0.c / (λ * 1e-8) * Lν[i] * prof[i] for (i, λ) in enumerate(λgrid)]
νLνc = [Z0.c / (λ * 1e-8) * Lν[i] for (i, λ) in enumerate(λgrid)]
serialize(joinpath(@__DIR__, "..", "lib", "spectra", "linespec_$on.jls"), (name = name, nlte = on, lines = out, λ = λgrid, νLν = νLν, νLνc = νLνc))
L(x) = out[findfirst(o -> o.name == x, out)].Lline
@printf("Balmer decrement L(Hα)/L(Hβ) = %.2f (case B 2.86); L(Paα)/L(Hα) = %.3f; L(Lyα)/L(Hα) = %.3f; total lines %.3e erg/s\n",
        L("Hα") / L("Hβ"), L("Paα") / L("Hα"), L("Lyα") / L("Hα"), sum(o.Lline for o in out))
println("saved lib/spectra/linespec_$on.jls")
