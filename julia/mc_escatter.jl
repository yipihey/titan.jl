# Monte Carlo transfer of a hydrogen line through the infall of a non-LTE solution (nlte_run.jl), with electron
# scattering: what the Thomson-thick, partly ionized infall does to the line profile.
#   Z0_TABLE=ext TITAN_NOLIB=1 julia -t 2 --project=env mc_escatter.jl NAME [line=Halpha|Hbeta] [N=400000] [es=1] [res=1]
#        [Te=…] [seed=1] [out=TAG]
# Photons are emitted where the solution's net Sobolev line power is positive, in proportion to it, with directions
# weighted by the angle-dependent Sobolev escape probability (1 − e^{−τ(μ)})/τ(μ), τ(μ) = K/Q(μ), and at the comoving
# line centre plus a thermal (+ ξ = 10 km/s) Doppler offset. They are then followed in straight lines through spherical
# shells of the solution (log-interpolated in r):
# * es=1: Thomson scattering on electrons (n_e of the solution, σ_T, dipole phase function); each electron moves with
#   the infall (−v_ff r̂) plus a Maxwellian thermal velocity at the gas temperature (or Te=…, a test); first-order
#   Doppler shift ν' = ν [1 + (n̂' − n̂)·v_e/c].
# * res=1: Sobolev resonances: where the photon's comoving frequency crosses line centre it interacts with probability
#   1 − e^{−τ(μ)} and is re-emitted at comoving line centre in a new escape-weighted direction (pure line scattering,
#   no destruction: an upper limit on line trapping).
# * Photons entering the star (r < r_min, the inner edge of the non-LTE region) are absorbed; at r > r_max they escape.
# Output: the escaping luminosity per velocity bin, V = −(lab-frame Doppler) in km/s (positive = redshift), normalized
# to the line luminosity, and fits of exponential wings on both sides. Writes lib/spectra/mc_NAME_TAG.jls.
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf, Random
const Z0 = Z0Physics; const NH = NLTEH
name = ARGS[1]
opt = Dict(Symbol(split(a, "=")[1]) => split(a, "=")[2] for a in ARGS[2:end])
line = get(opt, :line, "Halpha"); (l, u) = line == "Halpha" ? (2, 3) : line == "Hbeta" ? (2, 4) : line == "Paalpha" ? (3, 4) : error("line")
Nph = parse(Int, get(opt, :N, "400000")); es = get(opt, :es, "1") == "1"; res = get(opt, :res, "1") == "1"
Te_fix = haskey(opt, :Te) ? parse(Float64, opt[:Te]) : 0.0; seed = parse(Int, get(opt, :seed, "1"))
# theta=θ_o: bipolar outflow cones |cos θ| > cos θ_o (axis z) at +v_ff, infall elsewhere (as axisym_lines.jl); escaping
# photons are then also binned by inclination (|n_z|): 0–30°, 30–60°, 60–90° from the axis
θo = parse(Float64, get(opt, :theta, "0")); cosθo = cosd(θo)
tag = get(opt, :out, "$(line)_es$(Int(es))_res$(Int(res))" * (Te_fix > 0 ? "_Te$(opt[:Te])" : "") * (θo > 0 ? "_th$(Int(θo))" : ""))
dn = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$name.jls"))
M = 1e6 * Z0.Msun; c = Z0.c; σT = 6.6524e-25; ξturb = 10e5
r = dn.r; nr = length(r); lr = log.(r)
ν0 = (NH.En[u] - NH.En[l]) / Z0.h; λ0 = c / ν0; f = NH.fosc(l, u)
Blu = 4π * NH.πe2mc * f / (Z0.h * ν0); Bul = NH.gn[l] / NH.gn[u] * Blu; A = 2Z0.h * ν0^3 / c^2 * Bul
K = [NH.πe2mc * f * λ0 * (st.n[l] - NH.gn[l] / NH.gn[u] * st.n[u]) for st in dn.nl]      # Sobolev opacity factor [s⁻¹]
ne = [st.ne for st in dn.nl]; T = dn.T
# net line power per volume (Sobolev, with the solution's escape probability and the local continuum at the line)
Jl = dn.Jhα    # the continuum J at Hα; for other lines only the emission weights change slightly (net ≈ gross here)
P = [Z0.h * ν0 * st.β[l, u] * (st.n[u] * A - (line == "Halpha" ? (st.n[l] * Blu - st.n[u] * Bul) * Jl[q] : 0.0)) for (q, st) in enumerate(dn.nl)]
w = [max(P[q], 0.0) * dn.dV[q] for q in 1:nr]; w[K .<= 0] .= 0.0          # no emission from the inverted X-ray layer
cw = cumsum(w) ./ sum(w)
vff(x) = sqrt(2Z0.G * M / x)
rmin, rmax = r[1], r[end]
# log-interpolation of a radial table
function at(tab, x)
    x <= r[1] && return tab[1]; x >= r[end] && return tab[end]
    i = searchsortedlast(r, x); t = (log(x) - lr[i]) / (lr[i+1] - lr[i]); tab[i] + t * (tab[i+1] - tab[i])
end
lnne = log.(max.(ne, 1e-30)); neat(x) = exp(at(lnne, x))
lnK = log.(max.(K, 1e-300)); Kat(x) = (k = at(K, x); k <= 0 ? 0.0 : exp(at(lnK, x)))
Q(x, μ) = vff(x) / x * abs((1 - μ^2) - 0.5 * μ^2)                    # |dv_n/dn| for free fall (d ln v/d ln r = −1/2)
esc(τ) = τ < 1e-6 ? 1.0 : -expm1(-τ) / τ
# direction from the local escape probability: sample μ = r̂·n̂ by rejection, azimuth uniform
function emit_dir(rng, x)
    k = Kat(x)
    while true
        μ = 2rand(rng) - 1
        rand(rng) <= esc(k / max(Q(x, μ), 1e-30)) && return μ
    end
end
# unit vectors: position p (3-vector), direction n
function dir_from_mu(rng, p, μ)
    rh = p ./ sqrt(sum(abs2, p))
    a = abs(rh[1]) < 0.9 ? [1.0, 0, 0] : [0, 1.0, 0]
    e1 = a .- sum(a .* rh) .* rh; e1 ./= sqrt(sum(abs2, e1)); e2 = [rh[2]*e1[3]-rh[3]*e1[2], rh[3]*e1[1]-rh[1]*e1[3], rh[1]*e1[2]-rh[2]*e1[1]]
    φ = 2π * rand(rng); s = sqrt(max(1 - μ^2, 0.0))
    μ .* rh .+ s * cos(φ) .* e1 .+ s * sin(φ) .* e2
end
# gas velocity (infall) at p; photon frequency carried as V_lab = −(Doppler), in cm/s: ν = ν0 (1 − V/c)
vgas(p) = (x = sqrt(sum(abs2, p)); (θo > 0 && abs(p[3]) / x > cosθo ? 1.0 : -1.0) * vff(x) .* p ./ x)
thermal(x) = sqrt(2Z0.kB * at(T, x) / Z0.mH + ξturb^2)
Vbins = collect(range(-4000e5, 4000e5, length = 401)); H = [zeros(length(Vbins) - 1) for _ in 1:Threads.nthreads()]
Hi = [zeros(length(Vbins) - 1, 3) for _ in 1:Threads.nthreads()]      # by inclination bin
stats = [zeros(Int, 4) for _ in 1:Threads.nthreads()]    # escaped, absorbed by the star, electron scatterings, resonances
Threads.@threads for tid in 1:Threads.nthreads()
    rng = MersenneTwister(seed * 1000 + tid)
    for _ in 1:cld(Nph, Threads.nthreads())
        # emission point: radius from the cumulative weight (uniform within the shell in log r), random direction
        q = searchsortedfirst(cw, rand(rng)); q = clamp(q, 1, nr)
        x = q == 1 ? r[1] : exp(lr[q-1] + rand(rng) * (lr[q] - lr[q-1]))
        p = dir_from_mu(rng, [0.0, 0.0, 1.0], 2rand(rng) - 1) .* x
        n = dir_from_mu(rng, p, emit_dir(rng, x))
        V = -sum(vgas(p) .* n) + thermal(x) * randn(rng) / sqrt(2)          # lab-frame: comoving centre + thermal
        τ = -log(rand(rng)); alive = true; nsc = 0
        δold = NaN; armed = false      # after a line emission, resonance checks start once the photon has left its own
                                        # resonance zone (|δ| > 3 thermal widths): the local escape is already in β(μ)
        while alive
            x = sqrt(sum(abs2, p))
            if x > rmax
                k = searchsortedlast(Vbins, V)
                if 1 <= k < length(Vbins)
                    H[tid][k] += 1; ib = clamp(1 + floor(Int, acosd(min(abs(n[3]), 1.0)) / 30), 1, 3); Hi[tid][k, ib] += 1
                end
                stats[tid][1] += 1; break
            end
            if x < rmin * 0.9999
                stats[tid][2] += 1; break
            end
            ds = 0.02 * x
            # resonance: comoving offset δ = V − (−vgas·n) crosses zero within the step
            δ = V + sum(vgas(p) .* n)
            armed = armed || abs(δ) > 3thermal(x)
            if res && armed && isfinite(δold) && sign(δ) != sign(δold)
                μ = sum(p .* n) / x
                k = Kat(x)
                if k > 0 && rand(rng) < -expm1(-k / max(Q(x, μ), 1e-30))
                    stats[tid][4] += 1
                    n = dir_from_mu(rng, p, emit_dir(rng, x))
                    V = -sum(vgas(p) .* n) + thermal(x) * randn(rng) / sqrt(2)
                    δold = NaN; armed = false; nsc += 1; nsc > 10^5 && break
                    continue
                end
            end
            δold = δ
            if es
                dτ = neat(x) * σT * ds
                if τ < dτ                                                      # electron scattering inside this step
                    p = p .+ (τ / dτ * ds) .* n
                    x = sqrt(sum(abs2, p))
                    Tx = Te_fix > 0 ? Te_fix : at(T, x)
                    ve = vgas(p) .+ sqrt(Z0.kB * Tx / 9.1093837e-28) .* randn(rng, 3)
                    # dipole phase function by rejection
                    n2 = n
                    while true
                        n2 = dir_from_mu(rng, [0.0, 0.0, 1.0], 2rand(rng) - 1)
                        cθ = sum(n .* n2); rand(rng) * 2 <= 1 + cθ^2 && break
                    end
                    V = V - sum((n2 .- n) .* ve)          # ν' = ν[1 + (n'−n)·v_e/c] and V = −c δν/ν
                    n = n2; τ = -log(rand(rng)); stats[tid][3] += 1; δold = NaN; armed = true
                    nsc += 1; nsc > 10^5 && break
                    continue
                end
                τ -= dτ
            end
            p = p .+ ds .* n
        end
    end
end
h = sum(H); hi = sum(Hi); st = sum(stats); Vc = 0.5 .* (Vbins[1:end-1] .+ Vbins[2:end]) ./ 1e5; dV = (Vbins[2] - Vbins[1]) / 1e5
prof = h ./ (sum(h) * dV)                     # per km/s, normalized to the escaping line luminosity
Lline = sum(w)                               # emitted (net) line power [erg/s]
fesc = st[1] / Nph
k = argmax(prof); m = prof[k]; hh = findall(>(m / 2), prof); fwhm = Vc[hh[end]] - Vc[hh[1]]
# exponential wings: fit ln prof vs V over 600–2500 km/s on each side (where prof > 1e-3 of peak)
function wing(sel)
    xs = Vc[sel]; ys = log.(prof[sel]); ok = isfinite.(ys); xs, ys = xs[ok], ys[ok]
    length(xs) < 5 && return NaN
    xm = sum(xs) / length(xs); ym = sum(ys) / length(ys)
    s = sum((xs .- xm) .* (ys .- ym)) / sum((xs .- xm) .^ 2); -1 / s
end
Wr = wing(findall(v -> 600 <= v <= 2500, Vc) ∩ findall(>(1e-3m), prof)); Wb = -wing(findall(v -> -2500 <= v <= -600, Vc) ∩ findall(>(1e-3m), prof))
cen = sum(Vc .* prof) * dV
fr(a, b) = sum(prof[(Vc .>= a) .& (Vc .< b)]) * dV
@printf("%s %s: es=%d res=%d%s; %d photons: escaped %.3f, absorbed by the star %.3f; %.2f electron scatterings and %.2f line resonances per photon\n",
        name, line, es, res, Te_fix > 0 ? " Te=$(Te_fix)" : "", Nph, fesc, st[2] / Nph, st[3] / Nph, st[4] / Nph)
@printf("  emitted (net) line power %.3e erg/s, escaping %.3e; profile: peak at %+.0f km/s, FWHM %.0f km/s, centroid %+.0f km/s\n", Lline, Lline * fesc, Vc[k], fwhm, cen)
@printf("  exponential wing scale (600–2500 km/s): blue %.0f, red %.0f km/s (exponential FWHM 2 ln2 W: %.0f, %.0f); fraction beyond ±1000 km/s: blue %.3f, red %.3f\n",
        Wb, Wr, 2log(2) * Wb, 2log(2) * Wr, fr(-1e9, -1000), fr(1000, 1e9))
profi = [hi[:, b] ./ (sum(hi[:, b]) * dV) for b in 1:3]
if θo > 0
    for b in 1:3
        pb = profi[b]; mb = maximum(pb); hb = findall(>(mb / 2), pb)
        @printf("  inclination %2d–%2d°: FWHM %.0f km/s, centroid %+.0f km/s, beyond ±1000 km/s: blue %.3f, red %.3f\n", 30(b-1), 30b,
                Vc[hb[end]] - Vc[hb[1]], sum(Vc .* pb) * dV, sum(pb[Vc .< -1000]) * dV, sum(pb[Vc .>= 1000]) * dV)
    end
end
serialize(joinpath(@__DIR__, "..", "lib", "spectra", "mc_$(name)_$tag.jls"),
          (name = name, line = line, es = es, res = res, Te = Te_fix, V = Vc, prof = prof, fesc = fesc, Lline = Lline, stats = st ./ Nph,
           fwhm = fwhm, Wb = Wb, Wr = Wr, centroid = cen, θo = θo, profi = profi))
println("saved lib/spectra/mc_$(name)_$tag.jls")
