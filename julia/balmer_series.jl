# The Balmer series up to n = NMAX near the Balmer limit, for the Balmer break as observers measure it: one ray
# integration through the infall with all Balmer lines at once (their opacities and source functions added at every
# point, so that overlapping lines are treated correctly), on a fine wavelength grid from 3550 to 4150 Å, on top of the
# non-LTE continuum (which already contains the Balmer-continuum absorption and the jump at 3646 Å).
#   Z0_TABLE=ext TITAN_NOLIB=1 julia -t 2 --project=env balmer_series.jl NAME [nmax=40] [bhigh=b12|lte] [dv=15]
# Populations: n ≤ 12 from the non-LTE solution; for 12 < n ≤ nmax, n_n = b_n n_p n_e Φ_n(T) with b_n = b_12 (bhigh=b12)
# or b_n = 1 (bhigh=lte), which brackets the unknown upper levels (they matter only for the line source functions; the
# opacity is set by n = 2). Line profiles: Gaussian, thermal + ξ = 10 km/s, in the free-fall velocity field.
# Outputs: L_ν(λ); the break ⟨f_ν⟩(4000–4100 Å)/⟨f_ν⟩(3620–3720 Å) at full resolution and after smoothing to R = 100
# and R = 30 (Gaussian in ln λ, as a NIRSpec PRISM spectrum); the same without lines. Writes lib/spectra/balmer_NAME_TAG.jls.
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf
const SP = Spectrum; const Z0 = Z0Physics; const NH = NLTEH
name = ARGS[1]
opt = Dict(Symbol(split(a, "=")[1]) => split(a, "=")[2] for a in ARGS[2:end])
nmax = parse(Int, get(opt, :nmax, "40")); bhigh = get(opt, :bhigh, "b12"); dv = parse(Float64, get(opt, :dv, "15")) * 1e5
tag = "n$(nmax)_$bhigh"
d = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "$name.jls"))
dn = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$name.jls"))
S = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(d.run)_relax_snaps.jls"))
s = S[argmin([abs(x.t - d.t) for x in S])]
atm = haskey(dn, :Tstar) ? merge(d.atm, (T = dn.Tstar,)) : d.atm
M = 1e6 * Z0.Msun; c = Z0.c; h = Z0.h; kB = Z0.kB
rsh = s.r[d.ksh+1]; inf = dn.sel; rI = atm.r[inf]; TI = atm.T[inf]; rmax = rI[end]
function interp(x, xs, ys)
    x <= xs[1] && return ys[1]; x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x); t = (log(x) - log(xs[i])) / (log(xs[i+1]) - log(xs[i])); ys[i] + t * (ys[i+1] - ys[i])
end
ξturb = 10e5
nlfull = Any[nothing for _ in eachindex(atm.r)]; for (q, j) in enumerate(inf); nlfull[j] = dn.nl[q]; end
op = NH.NOpac(SP.Opac(atm), nlfull)
_, R = SP.rays(atm.r; ncore = 12); R = SP.weights!(atm, R)
# level populations up to nmax at every non-LTE point
λe3(T) = Z0.λe3(T)
Φn(n, T) = n^2 * exp(Z0.χH / (n^2 * kB * T)) / λe3(T)          # n_n*/(n_p n_e), g_n = 2n², g_p = 1, g_e = 2
pop = map(eachindex(dn.nl)) do q
    st = dn.nl[q]; T = TI[q]
    b12 = st.n[12] / (st.np * st.ne * Φn(12, T))
    [n <= 12 ? st.n[n] : (bhigh == "lte" ? 1.0 : b12) * st.np * st.ne * Φn(n, T) for n in 1:nmax]
end
νline(u) = Z0.χH * (1 / 4 - 1 / u^2) / h
lines = 3:nmax
# per line: opacity factor K = (πe²/mc) f λ (n_2 − g_2 n_u/g_u) and source function at every point; inverted points
# (the X-ray layer) transparent
Kr = Dict{Int,Vector{Float64}}(); Sr = Dict{Int,Vector{Float64}}()
for u in lines
    ν0 = νline(u); λ0 = c / ν0; f = NH.fosc(2, u); gl, gu = 8.0, 2.0u^2
    Kr[u] = [max(NH.πe2mc * f * λ0 * (p[2] - gl / gu * p[u]), 0.0) for p in pop]
    Sr[u] = [Kr[u][q] > 0 ? 2h * ν0^3 / c^2 / (pop[q][2] * gu / (pop[q][u] * gl) - 1) : 0.0 for q in eachindex(pop)]
end
uabs = dn.vabs
n = 1500; gq = [(exp(5 * (i - 1) / (n - 1)) - 1) / (exp(5) - 1) for i in 1:n]
RAYZ = map(1:length(R)) do j
    p = R[j].p
    znear = sqrt(max(rmax^2 - p^2, 0.0)); zfar = p < rsh ? sqrt(rsh^2 - p^2) : -znear
    (znear > zfar && p < rmax) || return nothing
    zs = p < rsh ? zfar .+ (znear - zfar) .* gq : (zp = znear .* gq; vcat(-reverse(zp), zp[2:end]))
    rz = sqrt.(p^2 .+ zs .^ 2)
    (z = zs, r = rz, w = [r < rsh ? NaN : interp(r, rI, uabs) * z / r for (z, r) in zip(zs, rz)],
     Δv = [sqrt(2kB * interp(r, rI, TI) / Z0.mH + ξturb^2) for r in rz],
     K = Dict(u => [r < rI[1] ? 0.0 : interp(r, rI, Kr[u]) for r in rz] for u in lines),
     S = Dict(u => [interp(r, rI, Sr[u]) for r in rz] for u in lines))
end
# wavelength grid (vacuum), uniform in ln λ
λs = exp.(range(log(3550.0), log(4150.0), step = dv / c))
# continuum: emergent intensity I(p) at a set of frequencies across the range (a pair straddling the edge), interpolated
λedge = 1e8 / (109677.58 / 4)
λc = sort(vcat(collect(range(3550.0, λedge - 1.0, length = 12)), collect(range(λedge + 1.0, 4150.0, length = 24))))
Icgrid = Vector{Vector{Float64}}(undef, length(λc)); ps = Float64[]
for (k, λ) in enumerate(λc)
    _, _, _, ci = SP.monochromatic(atm, op, R, c / (λ * 1e-8); moments = :I)
    Icgrid[k] = ci.I; global ps = ci.p
end
function Ic_at(λ)
    # interpolate on the same side of the edge only
    sel = λ < λedge ? findall(<(λedge), λc) : findall(>(λedge), λc)
    xs = λc[sel]; i = clamp(searchsortedlast(xs, λ), 1, length(xs) - 1); t = (λ - xs[i]) / (xs[i+1] - xs[i])
    Icgrid[sel[i]] .* (1 - t) .+ Icgrid[sel[i+1]] .* t
end
Lnu = zeros(length(λs)); Lc = zeros(length(λs))
Threads.@threads for iλ in eachindex(λs)
    λ = λs[iλ]; Ic = Ic_at(λ); Iv = similar(Ic)
    # lines that can matter at this wavelength (within 3000 km/s + 6 thermal widths)
    act = [u for u in lines if abs(c / νline(u) * 1e8 / λ - 1) < 3.2e8 / c]
    for j in eachindex(ps)
        I = Ic[j]; ry = RAYZ[j]
        if ry !== nothing && !isempty(act)
            for i in 1:length(ry.z)-1
                (isfinite(ry.w[i]) && isfinite(ry.w[i+1])) || continue
                D = max(ry.Δv[i], ry.Δv[i+1]); lo, hi = minmax(ry.w[i], ry.w[i+1])
                m = max(1, ceil(Int, (hi - lo) / (0.25 * min(ry.Δv[i], ry.Δv[i+1])))); dz = (ry.z[i+1] - ry.z[i]) / m
                for qq in 1:m
                    t = (qq - 0.5) / m; w = ry.w[i] + t * (ry.w[i+1] - ry.w[i]); Dt = ry.Δv[i] + t * (ry.Δv[i+1] - ry.Δv[i])
                    τt = 0.0; Sτ = 0.0
                    for u in act
                        v = (λ * 1e-8 * νline(u) / c - 1) * c            # this wavelength's velocity offset from line u
                        x = (v - w) / Dt; abs(x) > 6 && continue
                        Ku = ry.K[u][i] + t * (ry.K[u][i+1] - ry.K[u][i])
                        Ku <= 0 && continue
                        τ = Ku * exp(-x^2) / (sqrt(π) * Dt) * dz
                        τt += τ; Sτ += τ * (ry.S[u][i] + t * (ry.S[u][i+1] - ry.S[u][i]))
                    end
                    τt > 0 && (I = I * exp(-τt) + (Sτ / τt) * (-expm1(-τt)))
                end
            end
        end
        Iv[j] = I
    end
    Lnu[iλ] = 8π^2 * sum(0.5 * (Iv[j] * ps[j] + Iv[j+1] * ps[j+1]) * (ps[j+1] - ps[j]) for j in 1:length(ps)-1)
    Lc[iλ] = 8π^2 * sum(0.5 * (Ic[j] * ps[j] + Ic[j+1] * ps[j+1]) * (ps[j+1] - ps[j]) for j in 1:length(ps)-1)
end
# smoothing to resolving power Rp (Gaussian in ln λ with FWHM 1/Rp)
function smooth(y, Rp)
    σ = 1 / (Rp * 2sqrt(2log(2))); lnλ = log.(λs)
    [sum(y[k] * exp(-0.5 * ((lnλ[k] - lnλ[i]) / σ)^2) for k in eachindex(y)) / sum(exp(-0.5 * ((lnλ[k] - lnλ[i]) / σ)^2) for k in eachindex(y)) for i in eachindex(y)]
end
win(y, a, b) = (sel = (λs .>= a) .& (λs .<= b); sum(y[sel]) / count(sel))
brk(y) = win(y, 4000, 4100) / win(y, 3620, 3720)
L100 = smooth(Lnu, 100); L30 = smooth(Lnu, 30); Lc100 = smooth(Lc, 100)
@printf("%s, Balmer lines to n = %d, upper levels above 12: %s\n", name, nmax, bhigh)
@printf("break ⟨f_ν⟩(4000–4100)/⟨f_ν⟩(3620–3720): continuum only %.2f (R=100: %.2f); with the Balmer series %.2f, R=100 %.2f, R=30 %.2f\n",
        brk(Lc), brk(Lc100), brk(Lnu), brk(L100), brk(L30))
@printf("jump at the edge (continuum, 3640 vs 3652 Å): %.2f\n", Lc[argmin(abs.(λs .- 3652))] / Lc[argmin(abs.(λs .- 3640))])
for λ in (3600, 3640, 3660, 3680, 3700, 3720, 3750, 3800, 3900, 4050)
    k = argmin(abs.(λs .- λ)); @printf("  λ %5d Å: L_ν/L_ν,cont %.3f (R=100: %.3f)\n", λ, Lnu[k] / Lc[k], L100[k] / Lc100[k])
end
serialize(joinpath(@__DIR__, "..", "lib", "spectra", "balmer_$(name)_$tag.jls"), (λ = λs, Lnu = Lnu, Lc = Lc, L100 = L100, L30 = L30,
          brk = (cont = brk(Lc), cont100 = brk(Lc100), full = brk(Lnu), R100 = brk(L100), R30 = brk(L30)), nmax = nmax, bhigh = bhigh))
println("saved lib/spectra/balmer_$(name)_$tag.jls")
