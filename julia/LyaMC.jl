# S8a: Monte Carlo Lyα radiative transfer, 1D spherical, moving medium (radial velocities at the zone interfaces, linear inside a zone).
#
# Physics (docs/POPIII_NOTES.md, "S8a"): resonant scattering with the exact Voigt profile (Humlíček W4 of LymanWings), partial frequency redistribution in the atom frame with exact
# parallel-velocity sampling (Zheng & Miralda-Escudé 2002: split Cauchy/Gaussian rejection, exact for any split point), recoil, dipole angular redistribution (core: 2/3 isotropic + 1/3 Rayleigh =
# Hamilton's 11/24 + 3μ²/24, within the fine-structure splitting of the atom-frame frequency; wings: Rayleigh), the Kramers–Heisenberg wing cross section of `LymanWings` (ratio to the Lorentz
# wing, tabulated), Doppler shifts from the bulk flow along the ray (sub-stepped), destruction by collisional 2p → 2s mixing followed by two-photon decay, continuous absorption (H⁻ photodetachment
# etc.) and H₂ Lyman-band lines near 1216 Å. Core skipping (Smith et al. 2015 idea, exact variant): a scattering with |x| < x_c in a zone whose line mean free path at x_c is much shorter than
# the zone is replaced by a precomputed frequency walk (the real redistribution kernel, tables of exit frequency, number of core scatterings N and visits per frequency bin), the destruction acts on all N
# scatterings: (1 − ε)^N. The skipped photon does not move, so skipping is applied only to a photon deeper than D = safe · x_c · l(x_c) (l = mean free path at x_c) from both zone
# interfaces and from r = 0, in zones with D ≤ 0.2 Δr (the core walk is otherwise explicit); inside the layer D the photons are exact: a skipped walk next to a boundary changes the
# reflection of the layer and the force on it (test e: N_H = 1e19 shell, force 7.5 L/c with skipping everywhere against 40.6 exact; 40.8 with the layer).
# Force: Σ over the momentum balance of the radiation field, F_k = A_k Π_k − A_{k+1} Π_{k+1} + ∫(u − Π_rr)/r dV (crossing estimator of the radial momentum flux Π_rr at the interfaces and
# Δμ along the tracks, ∫(1 − μ²)/r ds = Δμ), which has no sign cancellation per event (the event estimator Σ(hν/c)μ_in has variance ∝ N_sc).
#
# Units: frequencies as the lab-frame velocity offset w = c (ν − ν0)/ν0 [cm/s]; x = (w − u μ)/b in the fluid frame (b = √(2kT/m_H)). Rates: each source has a photon rate [s⁻¹];
# every output is a rate (s⁻¹, erg s⁻¹, dyn).
module LyaMC
using Random, Printf
include(joinpath(@__DIR__, "LymanWings.jl"))
const LW = LymanWings

const c = LW.c; const hP = LW.hP; const kB = LW.kB; const mH = LW.mH; const eV = LW.eV
const ν0 = LW.ν_line(2); const E0 = hP * ν0
const A21 = LW.A_tot(2)                  # 6.265e8 s⁻¹
const σint = LW.πe2mc * LW.FN[2]         # ∫σ dν
const A2G = 8.2293                       # 2s → 1s two-photon rate [s⁻¹]
const DFS = 1.0969e10                    # 2p fine-structure splitting [Hz]

bth(T) = sqrt(2kB * T / mH)
dopp(T) = ν0 * bth(T) / c
damp(T) = A21 / (4π * dopp(T))
sig0(T) = σint / (sqrt(π) * dopp(T))
recoil_g(T) = E0 / (mH * c * bth(T))
voigt(a, x) = LW.voigt_H(a, x)

# ------------------------------------------------------------------ atomic data: destruction per line scattering, sources
"""
Probability per Lyα scattering of destruction by 2p → 2s mixing followed by two-photon decay (the chain 2p → 2s → 2p … is summed): q(2p → 2s) = 4.7e-4 (T/1e4)^{1/2} cm³ s⁻¹ for protons,
half of it for electrons (Seaton 1955; Osterbrock & Ferland 2006), 2s → 2p three times larger (statistical weights), A(2s → 2γ) = 8.23 s⁻¹. Saturates at A_2γ/(3A_21) = 4.4e-9 for n ≫ n_cr.
"""
function eps_collisional(T, ne, np)
    q = 4.7e-4 * sqrt(T / 1e4)
    Cps = q * (np + 0.5 * ne); Csp = 3 * Cps
    Pps = Cps / (A21 + Cps); Pd = A2G / (A2G + Csp)
    Pps * Pd / (1 - Pps * (1 - Pd))
end
"Fraction of case-B recombinations that end as Lyα (Cantalupo et al. 2008, Dijkstra 2014)."
f_lya_caseB(T) = (T4 = T / 1e4; clamp(0.686 - 0.106 * log10(T4) - 0.009 * T4^-0.44, 0.3, 0.75))
"Case-B recombination coefficient [cm³ s⁻¹] (Hui & Gnedin 1997)."
alphaB(T) = (λ = 315614.0 / T; 2.753e-14 * λ^1.5 / (1 + (λ / 2.74)^0.407)^2.242)
"Collisional 1s → 2p excitation rate coefficient by electrons [cm³ s⁻¹] (Callaway et al. 1985 fit as in Dijkstra 2014)."
gamma_1s2p(T) = 2.41e-6 / sqrt(T) * (T / 1e4)^0.22 * exp(-E0 / (kB * T))
"H⁻ photodetachment cross section [cm²] at wavelength λ [μm] (John 1988 / Wishart 1979, as in ContinuumOpacity)."
function sigma_Hm(λ)
    (λ >= 1.6419 || λ < 0.125) && return 0.0
    C = (152.519, 49.534, -118.858, 92.536, -34.194, 4.982)
    x = 1 / λ - 1 / 1.6419; s = 0.0
    for n in 6:-1:1; s = C[n] + s * sqrt(x); end
    1e-18 * λ^3 * x^1.5 * s
end

# H₂ B–X Lyman-band lines within ±3 Å of Lyα that start in v″ ≤ 3 (Herczeg et al. 2002, ApJ 572, 310, Table 6 pumping transitions; λ in Å (vacuum), v″, J″, E″ [eV])
const H2LINES = [
    (1212.425, 1, 11, 1.36), (1213.356, 1, 14, 1.79), (1213.677, 2, 12, 1.93), (1214.465, 1, 15, 1.94), (1214.781, 3, 5, 1.65), (1214.995, 3, 6, 1.72),
    (1215.726, 2, 6, 1.27), (1216.070, 2, 5, 1.20), (1217.031, 3, 2, 1.50), (1217.038, 3, 1, 1.48), (1217.205, 2, 0, 1.00), (1217.410, 0, 19, 2.20),
    (1217.643, 2, 1, 1.02), (1217.904, 1, 13, 1.64), (1217.982, 0, 18, 2.02), (1218.521, 1, 14, 1.79), (1219.089, 2, 2, 1.04), (1219.368, 2, 1, 1.02)]
const H2_F = Ref(0.01)         # oscillator strength of the lines (order of magnitude of the Lyman-band lines; the f of these lines are not in the local data)
const H2_GAMMA = 1.7e9         # B-state decay rate [s⁻¹] (radiative + predissociation: the damping of the line)
"LTE populations n(v″, J″) [cm⁻³] of the lines' lower levels for H₂ density `nH2` at T (rigid-rotor + vibration partition function, ortho/para weights 3:1)."
function h2_pops(nH2, T)
    kT = T * kB / eV
    Z = 0.0
    for J in 0:40
        g = (iseven(J) ? 1 : 3) * (2J + 1); Z += g * exp(-0.00737 * J * (J + 1) / kT)        # rotational ladder of v = 0
    end
    Z *= 1 + exp(-0.516 / kT) + exp(-1.0 / kT)
    [nH2 * (iseven(l[3]) ? 1 : 3) * (2l[3] + 1) * exp(-l[4] / kT) / Z for l in H2LINES]
end

# ------------------------------------------------------------------ KH wing correction  C(Δ) = σ_KH(ν) / σ_Lorentz(ν),  Δ = (ν − ν0)/ν0
const LG0 = -5.0; const LG1 = log10(0.3); const NLG = 200
const WING_P = Ref(Float64[]); const WING_N = Ref(Float64[])
function init_wings!()
    isempty(WING_P[]) || return
    lg = range(LG0, LG1, length = NLG)
    cp = Float64[]; cn = Float64[]
    for l in lg, (sg, arr) in ((1, cp), (-1, cn))
        ν = ν0 * (1 + sg * 10.0^l); push!(arr, LW.σ_coh(ν) / (σint * LW.φ_L(ν, 2)))
    end
    WING_P[] = cp; WING_N[] = cn
end
@inline function Cwing(Δ)
    aΔ = abs(Δ); aΔ < 1e-5 && return 1.0
    t = (log10(min(aΔ, 0.3)) - LG0) / (LG1 - LG0) * (NLG - 1) + 1
    i = clamp(floor(Int, t), 1, NLG - 1); f = t - i
    arr = Δ > 0 ? WING_P[] : WING_N[]
    (1 - f) * arr[i] + f * arr[i+1]
end

# ------------------------------------------------------------------ sampling
@inline function ray_mu(rng)
    cc = 2 * (2 * rand(rng) - 1)
    s3 = cc + sqrt(cc * cc + 1)
    s = cbrt(s3)
    s - 1 / s
end
const U0X = (0.9, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0); const U0Y = (0.8, 0.9, 1.35, 1.8, 2.55, 3.0, 3.25, 3.5, 4.0)
@inline function u0_of(x)
    x <= 0.9 && return 0.0
    x >= 10 && return 4.0
    i = 1
    @inbounds while U0X[i+1] < x; i += 1; end
    @inbounds U0Y[i] + (U0Y[i+1] - U0Y[i]) * (x - U0X[i]) / (U0X[i+1] - U0X[i])
end
"Atom velocity component along the photon direction [units of b] for a photon at fluid-frame x: p(u) ∝ e^{-u²}/((u − x)² + a²), exact (two-region rejection, Zheng & Miralda-Escudé 2002)."
@inline function sample_upar(rng, x, a)
    s = x < 0 ? -1.0 : 1.0; xa = abs(x)
    if xa < 0.9
        while true
            u = xa + a * tan(π * (rand(rng) - 0.5))
            rand(rng) < exp(-u * u) && return s * u
        end
    end
    u0 = u0_of(xa)
    LA = 1 / ((u0 - xa)^2 + a * a)
    MA = LA * 0.5 * sqrt(π) * (1 + erf1(u0))
    th0 = atan((u0 - xa) / a)
    MB = exp(-u0 * u0) / a * (π / 2 - th0)
    pA = MA / (MA + MB)
    while true
        if rand(rng) < pA
            # Gaussian N(0, 1/2) truncated to u < u0, accept with L(u)/L(u0)
            u = randn(rng) * 0.7071067811865476
            u < u0 || continue
            rand(rng) < ((u0 - xa)^2 + a * a) / ((u - xa)^2 + a * a) && return s * u
        else
            # Cauchy truncated to u > u0, accept with e^{-u²}/e^{-u0²}
            θ = th0 + rand(rng) * (π / 2 - th0)
            u = xa + a * tan(θ)
            rand(rng) < exp(-(u * u - u0 * u0)) && return s * u
        end
    end
end
# erf(z) for the envelope mass (Abramowitz–Stegun 7.1.26, |error| < 1.5e-7)
@inline function erf1(z)
    sg = z < 0 ? -1.0 : 1.0; z = abs(z); t = 1 / (1 + 0.3275911 * z)
    sg * (1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * exp(-z * z))
end

"""
One resonant scattering of a photon of fluid-frame frequency x (Doppler units): returns (x_out, cosΘ). Atom velocity: parallel component from the exact conditional distribution, perpendicular
Gaussian; atom-frame coherence with recoil `g`; angular redistribution isotropic (prob 2/3) or Rayleigh inside the fine-structure core |x_atom| < xfs, Rayleigh outside.
"""
@inline function scatter_x(rng, x, a, g, xfs)
    u = sample_upar(rng, x, a)
    xa = x - u
    μs = (abs(xa) < xfs && rand(rng) < 2 / 3) ? 2 * rand(rng) - 1 : ray_mu(rng)
    s2 = sqrt(max(1 - μs * μs, 0.0))
    xo = xa - g * (1 - μs) + u * μs + randn(rng) * 0.7071067811865476 * s2
    (xo, μs)
end

# ------------------------------------------------------------------ core-skipping tables
const XB = 0.1; const XBM = 12.0; const NXB = 240          # visit / J(x) bins: x ∈ [−12, 12]
@inline xbin(x) = (abs(x) >= XBM ? 0 : floor(Int, (x + XBM) / XB) + 1)
const TABCACHE = Dict{Any,Any}()          # skip tables per (ln T bucket, x_c, M): built once per process
struct SkipTab
    xc::Float64; nodes::Vector{Float64}; xs::Matrix{Float64}; ns::Matrix{Int32}; vis::Matrix{Float64}        # vis[node, bin]: mean visits to the bin before the exit
end
"Frequency-walk table for the core |x| < xc: for entry nodes |x| (the scattering at the entry frequency is the first of the walk) the exit frequency (signed) and number of scatterings of M walks."
function skiptab(a, g, xfs, xc; M = 384, Nn = 9, seed = 7)
    nodes = collect(range(0.0, xc * 0.98, length = Nn))
    xs = zeros(M, Nn); ns = zeros(Int32, M, Nn); vis = zeros(Nn, NXB)
    Threads.@threads for j in 1:Nn
        rng = Xoshiro(seed + 101 * j); x0 = nodes[j]
        for m in 1:M
            x = x0; n = 0
            while true
                x, _ = scatter_x(rng, x, a, g, xfs); n += 1
                abs(x) >= xc && break
                ib = xbin(x); ib > 0 && (vis[j, ib] += 1)
                n > 10^7 && break
            end
            xs[m, j] = x; ns[m, j] = n
        end
    end
    vis ./= M
    SkipTab(xc, nodes, xs, ns, vis)
end

# ------------------------------------------------------------------ medium
"""
Zones of the transfer. `r` interfaces (n + 1), `nHI` [cm⁻³] neutral H, `T` [K], `u` radial velocity at the interfaces [cm/s] (> 0 outwards), `ne`, `np` for the destruction probability,
`kabs` continuous absorption [cm⁻¹] at Lyα (H⁻ photodetachment, …), `nH2` for the Lyman-band lines. `eps_extra` multiplies the destruction probability (scaling of `inner_column`).
`xcs` = candidate core-skipping frequencies (largest admissible one per zone: line mean free path at x_c × `safe` × x_c ≤ zone width); `xcs = ()` switches skipping off.
"""
struct Zones
    n::Int; r::Vector{Float64}; u::Vector{Float64}; du::Vector{Float64}
    nHI::Vector{Float64}; T::Vector{Float64}; b::Vector{Float64}; a::Vector{Float64}; nsig::Vector{Float64}; g::Vector{Float64}; xfs::Vector{Float64}
    eps::Vector{Float64}; kabs::Vector{Float64}; epsA::Vector{Float64}
    xc::Vector{Float64}; lskip::Vector{Float64}; tab::Vector{Int}; tabs::Vector{SkipTab}
    h2nu::Vector{Vector{Float64}}; h2amp::Vector{Vector{Float64}}; h2D::Vector{Vector{Float64}}; h2a::Vector{Vector{Float64}}
    vol::Vector{Float64}; wscale::Float64; dlayer::Float64
end
function build_zones(r, nHI, T, u; ne = zeros(length(T)), np = zeros(length(T)), kabs = zeros(length(T)), nH2 = zeros(length(T)), eps_extra = 1.0, g_extra = 1.0, wscale = 1.0, dlayer = 0.5, xcs = (6.0, 3.0, 2.0, 1.2), safe = 4.0,
                     tabM = 192, epsfloor = 0.0)
    init_wings!()
    n = length(T); @assert length(r) == n + 1 && length(u) == n + 1
    du = [(u[k+1] - u[k]) / (r[k+1] - r[k]) for k in 1:n]
    b = bth.(T); a = damp.(T); D = dopp.(T); ns = nHI .* sig0.(T); g = recoil_g.(T) .* g_extra; xfs = DFS ./ D
    eps = [min(eps_collisional(T[k], ne[k], np[k]) * eps_extra + epsfloor, 1.0) for k in 1:n]
    xc = zeros(n); lskip = zeros(n); tabid = zeros(Int, n); tabs = SkipTab[]; cache = Dict{Any,Int}()          # (tables of the run are shared through TABCACHE)
    h2nu = [Float64[] for _ in 1:n]; h2amp = [Float64[] for _ in 1:n]; h2D = [Float64[] for _ in 1:n]; h2a = [Float64[] for _ in 1:n]
    epsA = zeros(n)
    for k in 1:n
        if nH2[k] > 0
            pops = h2_pops(nH2[k], T[k]); bH2 = sqrt(2kB * T[k] / (2mH))
            for (l, (λ, v, J, E)) in enumerate(H2LINES)
                ν = c / (λ * 1e-8); Dl = ν * bH2 / c
                push!(h2nu[k], ν); push!(h2D[k], Dl); push!(h2a[k], H2_GAMMA / (4π * Dl)); push!(h2amp[k], pops[l] * σint * H2_F[] / (sqrt(π) * Dl))
            end
        end
        epsA[k] = kabs[k] > 0 ? kabs[k] / (ns[k] * 0.6 + kabs[k]) : 0.0
        if !isempty(h2nu[k]) && ns[k] > 0
            # H₂ lines inside the Doppler core: probability per core scattering that the interaction is a line absorption, averaged over the core profile (weights H(a, x))
            sa = 0.0; sw = 0.0
            for x in range(-6.0, 6.0, length = 121)
                ν = ν0 * (1 + x * b[k] / c); kA = 0.0
                for l in eachindex(h2nu[k]); xl = (ν - h2nu[k][l]) / h2D[k][l]; abs(xl) < 40 && (kA += h2amp[k][l] * voigt(h2a[k][l], xl)); end
                kL = ns[k] * voigt(a[k], x); hw = voigt(a[k], x); sa += hw * (kA + kabs[k]) / (kL + kA + kabs[k]); sw += hw
            end
            epsA[k] = max(epsA[k], sa / sw)
        end
        # core skipping: largest candidate for which the photon at x_c has many scatterings per zone crossing
        dr = r[k+1] - r[k]
        if ns[k] > 0
            for xcand in sort(collect(xcs); rev = true)
                lmfp = 1 / (ns[k] * voigt(a[k], xcand))
                if lmfp * max(xcand, 1.0) * safe <= 0.2 * dr
                    xc[k] = xcand; lskip[k] = lmfp * max(xcand, 1.0) * safe
                    key = (round(log(T[k]) / 0.25), xcand, tabM)
                    id = get!(cache, key) do
                        push!(tabs, get!(TABCACHE, key) do; skiptab(a[k], g[k], xfs[k], xcand; M = tabM); end); length(tabs)
                    end
                    tabid[k] = id; break
                end
            end
        end
    end
    vol = [4π / 3 * (r[k+1]^3 - r[k]^3) for k in 1:n]
    Zones(n, collect(float.(r)), collect(float.(u)), du, collect(float.(nHI)), collect(float.(T)), b, a, ns, g, xfs, eps, collect(float.(kabs)), epsA, xc, lskip, tabid, tabs, h2nu, h2amp, h2D, h2a, vol, float(wscale), float(dlayer))
end

# ------------------------------------------------------------------ sources
"""
Photon sources. `kind`: `:zones` (rates per zone, born thermally at line centre, uniform in volume, isotropic), `:shell` (photons from radius `r0` outward with the frequency spectrum
`wedges`/`wts` — velocity-offset bin edges [cm/s] and weights — and Lambert directions), `:point` (at r = 0, isotropic, line centre), `:line` (at r0, outward, line centre, thermal).
"""
struct Source
    kind::Symbol; name::Symbol; rate::Float64; npk::Int; r0::Float64; zrates::Vector{Float64}; wedges::Vector{Float64}; wcdf::Vector{Float64}
end
src_zones(rates; npk, name = :zones) = (s = sum(rates); Source(:zones, name, s, npk, 0.0, cumsum(rates) ./ max(s, 1e-300), Float64[], Float64[]))
src_point(rate; npk, name = :point) = Source(:point, name, rate, npk, 0.0, Float64[], Float64[], Float64[])
src_sphere(r0, rate; npk, name = :sphere) = Source(:sphere, name, rate, npk, r0, Float64[], Float64[], Float64[])
src_line(r0, rate; npk, name = :line) = Source(:line, name, rate, npk, r0, Float64[], Float64[], Float64[])
function src_shell(r0, rate, wedges, wts; npk, name = :shell)
    cdf = cumsum(wts) ./ sum(wts)
    Source(:shell, name, rate, npk, r0, Float64[], collect(float.(wedges)), cdf)
end

# ------------------------------------------------------------------ tallies
mutable struct Tally
    nz::Int
    nsc::Vector{Float64}; nsk::Vector{Float64}; mom::Vector{Float64}; hexch::Vector{Float64}; hrec::Vector{Float64}; Egas::Vector{Float64}
    nHm::Vector{Float64}; nH2::Vector{Float64}; n2g::Vector{Float64}; Jpath::Vector{Float64}; Jwin::Vector{Float64}; Prr::Vector{Float64}; Gth::Vector{Float64}; Xi::Vector{Float64}; Jskip::Vector{Float64}
    Jx::Matrix{Float64}; momsrc::Matrix{Float64}
    out::Matrix{Float64}; inn::Vector{Float64}
    emit::Float64; nout::Float64; nin::Float64; nlost::Float64; Eemit::Float64; Eout::Float64; Ein::Float64; E2g::Float64; Eabs::Float64; Eskipmiss::Float64
    nscphot::Float64
    nbw::Int
end
Tally(nz, nsrc, nbw) = Tally(nz, zeros(nz), zeros(nz), zeros(nz), zeros(nz), zeros(nz), zeros(nz), zeros(nz), zeros(nz), zeros(nz), zeros(nz), zeros(nz), zeros(nz), zeros(nz), zeros(nz + 1), zeros(nz), zeros(nz, NXB), zeros(nsrc, nz),
                              zeros(nsrc, nbw), zeros(nbw), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nbw)
function merge!(A::Tally, B::Tally)
    for f in (:nsc, :nsk, :mom, :hexch, :hrec, :Egas, :nHm, :nH2, :n2g, :Jpath, :Jwin, :Prr, :Gth, :Xi, :Jskip, :Jx, :momsrc, :out, :inn)
        getfield(A, f) .+= getfield(B, f)
    end
    for f in (:emit, :nout, :nin, :nlost, :Eemit, :Eout, :Ein, :E2g, :Eabs, :Eskipmiss, :nscphot); setfield!(A, f, getfield(A, f) + getfield(B, f)); end
    A
end

struct Opts
    inner::Symbol            # :absorb (star / inner boundary absorbs), :escape (counted as an inner escape: slab test), :reflect (albedo `albedo`)
    albedo::Function         # w -> reflection probability at the inner boundary (inner = :reflect)
    maxsc::Int               # packet cut (explicit scatterings)
    wmax::Float64; nbw::Int  # escape-spectrum bins: w ∈ [−wmax, wmax]
    dxmax::Float64           # max frequency change per ray sub-step (Doppler units, or 10 % of |x|)
    xwin::Float64            # J_α window |x| < xwin
    xsafe_skip::Bool
end
Opts(; inner = :absorb, albedo = w -> 0.0, maxsc = 20_000_000, wmax = 2e9, nbw = 2000, dxmax = 0.25, xwin = 1.0, xsafe_skip = true) = Opts(inner, albedo, maxsc, wmax, nbw, dxmax, xwin, xsafe_skip)

# ------------------------------------------------------------------ transport
@inline energy(w) = E0 * (1 + w / c)
@inline function kappas(Z::Zones, k, x, wf)
    κL = Z.nsig[k] * voigt(Z.a[k], x) * Cwing(Z.wscale * wf / c)
    κA = Z.kabs[k]
    nl = length(Z.h2nu[k])
    if nl > 0
        ν = ν0 * (1 + wf / c)
        for l in 1:nl
            xl = (ν - Z.h2nu[k][l]) / Z.h2D[k][l]
            abs(xl) < 40 && (κA += Z.h2amp[k][l] * voigt(Z.h2a[k][l], xl))
        end
    end
    (κL, κA)
end
@inline function wbin(w, wmax, nbw)
    i = floor(Int, (w + wmax) / (2wmax) * nbw) + 1
    clamp(i, 1, nbw)
end

"Trace one packet (rate `wt` photons/s) to its end. Position r, radial direction cosine μ, lab-frame velocity offset w, zone k."
function trace!(rng, Z::Zones, Tl::Tally, o::Opts, isrc, r, μ, w, wt, k)
    nexp = 0
    cinv = 1 / c
    (isfinite(w) && isfinite(μ) && isfinite(r) && isfinite(wt)) || (Tl.nlost += isfinite(wt) ? wt : 0.0; return)      # FULL-1: non-finite packet (source sampling)
    @inbounds while true
        rk = Z.r[k]; rk1 = Z.r[k+1]; rμ = r * μ
        sout = -rμ + sqrt(max(rμ * rμ + rk1 * rk1 - r * r, 0.0))
        sin_ = Inf
        if μ < 0
            d = rμ * rμ - (r * r - rk * rk)
            if d > 0; sin_ = max(-rμ - sqrt(d), 0.0); end
        end
        crossin = sin_ < sout
        smax = crossin ? sin_ : sout
        b = Z.b[k]; uk = Z.u[k]; duk = Z.du[k]
        τt = -log(1 - rand(rng))
        s = 0.0; hit = false; sint = 0.0
        lossy = Z.nsig[k] > 0 || Z.kabs[k] > 0 || !isempty(Z.h2nu[k])
        if lossy
            rs = r; μs = μ
            us = (uk + duk * (rs - rk)); wf = w - us * μs; xs = wf / b
            κL, κA = kappas(Z, k, xs, wf); κs = κL + κA
            while s < smax
                # sub-step
                dxds = abs(duk * μs * μs + us * (1 - μs * μs) / max(rs, 1e-300)) / b
                dsx = dxds > 0 ? max(o.dxmax, 0.1 * abs(xs)) / dxds : Inf
                ds = min(smax - s, dsx, κs > 0 ? 1.5 * τt / κs : Inf)
                ds <= 0 && (ds = smax - s)
                s2 = s + ds
                r2 = sqrt(max(r * r + 2 * rμ * s2 + s2 * s2, 0.0)); μ2 = r2 > 0 ? clamp((rμ + s2) / r2, -1.0, 1.0) : μs
                u2 = uk + duk * (r2 - rk); wf2 = w - u2 * μ2; x2 = wf2 / b
                κL2, κA2 = kappas(Z, k, x2, wf2); κe = κL2 + κA2
                Δτ = 0.5 * (κs + κe) * ds
                if Δτ >= τt
                    fr = τt / Δτ; sint = s + fr * ds; hit = true
                    dpath = fr * ds; xm = xs + 0.5 * fr * (x2 - xs)
                    Tl.Jpath[k] += wt * dpath; Tl.Prr[k] += wt * dpath * (0.5 * (μs + μs + fr * (μ2 - μs)))^2
                    rh = sqrt(max(r * r + 2 * rμ * sint + sint * sint, 0.0)); Tl.Gth[k] += wt * ((rh > 0 ? clamp((rμ + sint) / rh, -1.0, 1.0) : μ) - μs)
                    ib = xbin(xm); ib > 0 && (Tl.Jx[k, ib] += wt * dpath)
                    abs(xm) < o.xwin && (Tl.Jwin[k] += wt * dpath)
                    break
                end
                τt -= Δτ
                xm = 0.5 * (xs + x2)
                Tl.Jpath[k] += wt * ds; Tl.Prr[k] += wt * ds * (0.5 * (μs + μ2))^2
                Tl.Gth[k] += wt * (μ2 - μs)
                ib = xbin(xm); ib > 0 && (Tl.Jx[k, ib] += wt * ds)
                abs(xm) < o.xwin && (Tl.Jwin[k] += wt * ds)
                s = s2; rs = r2; μs = μ2; us = u2; wf = wf2; xs = x2; κs = κe
            end
        else
            Tl.Jpath[k] += wt * smax; Tl.Prr[k] += wt * smax * (μ * μ); Tl.Gth[k] += wt * ((let re = sqrt(max(r * r + 2 * rμ * smax + smax * smax, 0.0)); re > 0 ? clamp((rμ + smax) / re, -1.0, 1.0) : μ end) - μ)
            s = smax
        end
        if !hit
            # boundary
            rn = crossin ? rk : rk1
            μn = rn > 0 ? clamp((rμ + smax) / rn, -1.0, 1.0) : μ
            r = rn; μ = μn
            Tl.Xi[crossin ? k : k + 1] += wt * abs(μn)            # radial momentum flux through the interface (crossing estimator of Π_rr)
            if crossin
                if k == 1
                    if o.inner === :escape
                        Tl.nin += wt; Tl.Ein += wt * energy(w); Tl.inn[wbin(w, o.wmax, o.nbw)] += wt; return
                    elseif o.inner === :reflect && rand(rng) < o.albedo(w)
                        μ = sqrt(rand(rng)); Tl.Xi[1] += wt * μ; continue
                    else
                        Tl.nin += wt; Tl.Ein += wt * energy(w); Tl.inn[wbin(w, o.wmax, o.nbw)] += wt; return
                    end
                end
                k -= 1
            else
                if k == Z.n
                    Tl.nout += wt; Tl.Eout += wt * energy(w); Tl.out[isrc, wbin(w, o.wmax, o.nbw)] += wt; return
                end
                k += 1
            end
            continue
        end
        # interaction at s = sint
        r2 = sqrt(max(r * r + 2 * rμ * sint + sint * sint, 0.0)); μ = r2 > 0 ? clamp((rμ + sint) / r2, -1.0, 1.0) : μ; r = r2
        us = uk + duk * (r - rk); wf = w - us * μ; x = wf / b
        κL, κA = kappas(Z, k, x, wf)
        if rand(rng) * (κL + κA) > κL
            # continuous / H₂ absorption
            Ew = energy(w)
            if rand(rng) * κA < Z.kabs[k]; Tl.nHm[k] += wt; else Tl.nH2[k] += wt; end
            Tl.Egas[k] += wt * Ew; Tl.Eabs += wt * Ew
            Tl.mom[k] += wt * Ew / c * μ; Tl.momsrc[isrc, k] += wt * Ew / c * μ
            return
        end
        # ---- scattering
        Ew = energy(w)
        Tl.nsc[k] += wt; Tl.mom[k] += wt * Ew / c * μ; Tl.momsrc[isrc, k] += wt * Ew / c * μ; Tl.nscphot += wt
        nexp += 1
        nexp > o.maxsc && (Tl.nlost += wt; return)
        eps = Z.eps[k]; epsA = Z.epsA[k]
        g = Z.g[k]; bc = b * cinv
        tid = Z.tab[k]
        if tid > 0 && abs(x) < Z.xc[k] && r > Z.lskip[k] && min(r - rk, rk1 - r) > Z.dlayer * Z.lskip[k]
            tb = Z.tabs[tid]
            ax = abs(x); j = clamp(floor(Int, ax / tb.xc / 0.98 * (length(tb.nodes) - 1)) + 1, 1, length(tb.nodes) - 1)
            jj = (rand(rng) < (ax - tb.nodes[j]) / (tb.nodes[j+1] - tb.nodes[j])) ? j + 1 : j
            m = rand(rng, 1:size(tb.xs, 1)); xo = tb.xs[m, jj] * (x < 0 ? -1 : 1); N = Int(tb.ns[m, jj])
            Tl.nsk[k] += wt * N; Tl.nsc[k] += wt * (N - 1); Tl.nscphot += wt * (N - 1)
            # visits in the core: J(x) from the table (expected path 1/κ per visit)
            pskip = 0.0
            for ib in 1:NXB
                v = tb.vis[jj, ib]; v == 0 && continue
                xb = -XBM + (ib - 0.5) * XB; xb = x < 0 ? -xb : xb; ibb = xbin(xb)
                κLv, κAv = kappas(Z, k, xb, xb * b); pth = wt * v / (κLv + κAv)
                Tl.Jpath[k] += pth; Tl.Prr[k] += pth / 3; Tl.Jskip[k] += pth; r > 0 && (Tl.Gth[k] += pth * (2 / 3) / r); ibb > 0 && (Tl.Jx[k, ibb] += pth); abs(xb) < o.xwin && (Tl.Jwin[k] += pth)
                pskip += pth
            end
            # momentum of the skipped absorptions: the radial cosine of a photon re-emitted isotropically at r and absorbed after a flight l is (r μ + l)/r', mean (2/3) l/r (the rotation of the radial
            # direction along the ray); the skipped chain does not move the photon, so this expectation replaces the deposit of its N − 1 interior absorptions
            pskip > 0 && r > 0 && (dep = (2 / 3) * (Ew / c) * pskip / r; Tl.mom[k] += dep; Tl.momsrc[isrc, k] += dep)
            pk = 1 - (1 - eps - epsA)^N
            if rand(rng) < pk
                # destroyed during the core walk
                if rand(rng) * (eps + epsA) < eps
                    Tl.n2g[k] += wt; Tl.E2g += wt * Ew
                else
                    Tl.nHm[k] += wt; Tl.Egas[k] += wt * Ew; Tl.Eabs += wt * Ew
                end
                return
            end
            μn = 2 * rand(rng) - 1
            wn = xo * b + us * μn
            isfinite(wn) || (Tl.nlost += wt; return)         # FULL-1: a non-finite photon frequency (found once, 1e4 yr into the production run) ends the packet as lost instead of the run
            Tl.hexch[k] += wt * E0 * bc * (x - xo); Tl.hrec[k] += wt * E0 * bc * g * N
            Tl.Egas[k] += wt * (Ew - energy(wn))
            μ = μn; w = wn
            continue
        end
        if eps > 0 && rand(rng) < eps
            Tl.n2g[k] += wt; Tl.E2g += wt * Ew; return
        end
        xo, μsc = scatter_x(rng, x, Z.a[k], g, Z.xfs[k])
        s2 = sqrt(max(1 - μsc * μsc, 0.0))
        μn = clamp(μsc * μ + s2 * sqrt(max(1 - μ * μ, 0.0)) * cos(2π * rand(rng)), -1.0, 1.0)
        wn = xo * b + us * μn
        isfinite(wn) || (Tl.nlost += wt; return)
        Tl.hexch[k] += wt * E0 * bc * (x - xo); Tl.hrec[k] += wt * E0 * bc * g
        Tl.Egas[k] += wt * (Ew - energy(wn))
        μ = μn; w = wn
    end
end

# ------------------------------------------------------------------ emission
function emit_packet(rng, Z::Zones, s::Source)
    if s.kind === :zones
        ξ = rand(rng); k = searchsortedfirst(s.zrates, ξ); k = clamp(k, 1, Z.n)
        rk = Z.r[k]; rk1 = Z.r[k+1]
        r = cbrt(rk^3 + rand(rng) * (rk1^3 - rk^3)); μ = 2 * rand(rng) - 1
        x = randn(rng) * 0.7071067811865476 + Z.a[k] * tan(π * (rand(rng) - 0.5))
        u = Z.u[k] + Z.du[k] * (r - rk)
        return (r, μ, x * Z.b[k] + u * μ, k)
    elseif s.kind === :point
        k = 1; μ = 2 * rand(rng) - 1
        x = randn(rng) * 0.7071067811865476 + Z.a[k] * tan(π * (rand(rng) - 0.5))
        return (0.0, μ, x * Z.b[k] + 0.0, k)
    elseif s.kind === :sphere
        r = s.r0; k = clamp(searchsortedlast(Z.r, r), 1, Z.n); μ = 2 * rand(rng) - 1
        x = randn(rng) * 0.7071067811865476 + Z.a[k] * tan(π * (rand(rng) - 0.5))
        u = Z.u[k] + Z.du[k] * (r - Z.r[k])
        return (r, μ, x * Z.b[k] + u * μ, k)
    elseif s.kind === :line
        r = s.r0; k = clamp(searchsortedlast(Z.r, r), 1, Z.n); μ = sqrt(rand(rng))
        x = randn(rng) * 0.7071067811865476 + Z.a[k] * tan(π * (rand(rng) - 0.5))
        u = Z.u[k] + Z.du[k] * (r - Z.r[k])
        return (r, μ, x * Z.b[k] + u * μ, k)
    else  # :shell
        r = s.r0; k = clamp(searchsortedlast(Z.r, r), 1, Z.n); μ = sqrt(rand(rng))
        i = clamp(searchsortedfirst(s.wcdf, rand(rng)), 1, length(s.wcdf))
        w = s.wedges[i] + rand(rng) * (s.wedges[i+1] - s.wedges[i])
        return (r, μ, w, k)
    end
end

# ------------------------------------------------------------------ driver
"""
    run(Z, sources; nbatch = 8, seed = 1, opts = Opts()) -> Result

Monte Carlo of all sources; the packets of each source are split over `nbatch` independent batches (threads) whose scatter gives the error estimates. Fills `LyaMC.result`.
"""
function run(Z::Zones, sources::Vector{Source}; nbatch = 8, seed = 1, opts::Opts = Opts(), store = true)
    nsrc = length(sources); n = Z.n
    tl = [Tally(n, nsrc, opts.nbw) for _ in 1:nbatch]
    Threads.@threads :dynamic for ib in 1:nbatch
        rng = Xoshiro(seed * 1000003 + ib * 7919)
        T = tl[ib]
        for (is, s) in enumerate(sources)
            nb = s.npk ÷ nbatch; nb < 1 && (nb = 1)
            wt = s.rate / (nb * nbatch)
            for _ in 1:nb
                r, μ, w, k = emit_packet(rng, Z, s)
                T.emit += wt; T.Eemit += wt * energy(w)
                s.kind === :shell && r == Z.r[k] && (T.Xi[k] += wt * μ)           # injected radial momentum flux at an interface
                trace!(rng, Z, T, opts, is, r, μ, w, wt, k)
            end
        end
    end
    tot = Tally(n, nsrc, opts.nbw); for T in tl; merge!(tot, T); end
    res = finalize(Z, sources, tl, tot, opts)
    store && (RESULT.r = res)
    res
end

# per-zone physics from the tallies
function finalize(Z::Zones, sources, tl, tot::Tally, o::Opts)
    n = Z.n; nb = length(tl)
    Dnu = dopp.(Z.T)
    V = Z.vol
    Jint = tot.Jpath ./ (4π .* V)                                 # ∫ J^phot dν [photons cm⁻² s⁻¹ sr⁻¹]
    Jalpha = tot.Jwin ./ (4π .* V .* (2 * o.xwin .* Dnu))        # J^phot_ν averaged over |x| < xwin  [photons cm⁻² s⁻¹ Hz⁻¹ sr⁻¹]
    Rsc = [Z.nHI[k] > 0 ? tot.nsc[k] / (Z.nHI[k] * V[k]) : 0.0 for k in 1:n]       # scatterings per H atom per second (the WF coupling P_α/…)
    force = tot.mom                                                 # radial force on the gas of the zone [dyn]
    # colour temperature from the odd part of ln J(x) (J(x) tilted by recoil: d ln J/dx = −2 g T_k/T_c, Hirata 2006)
    Tc = copy(Z.T); Tcerr = fill(Inf, n); Tvalid = falses(n)
    function slope(Jx, k)
        sw = 0.0; sx = 0.0
        for ib in NXB÷2+1:NXB-8
            xp = -XBM + (ib - 0.5) * XB; im = NXB + 1 - ib
            ap = Jx[k, ib]; am = Jx[k, im]
            (ap > 0 && am > 0 && xp > 0.3 && xp < 6) || continue
            y = log(ap / am); wgt = 1 / (1 / ap + 1 / am)
            sx += wgt * xp * (-y / 2) ; sw += wgt * xp^2
        end
        sw > 0 ? sx / sw : NaN
    end
    for k in 1:n
        Z.nHI[k] > 0 || continue
        s = slope(tot.Jx, k); s0 = 2 * Z.g[k]
        sb = [slope(tl[ib].Jx, k) for ib in 1:nb]; sb = filter(isfinite, sb)
        if isfinite(s) && s > 0 && length(sb) >= 3
            m = sum(sb) / length(sb); sd = sqrt(sum((sb .- m) .^ 2) / (length(sb) - 1) / length(sb))
            tcc = Z.T[k] * s0 / s; Tcerr[k] = tcc * sd / max(s, 1e-300)
            if sd < 0.15 * s && 0.3 < tcc / Z.T[k] < 3.0; Tc[k] = tcc; Tvalid[k] = true; end        # (a tilt measured with > 15 % scatter or implausible stays T_k: the true departures are ≲ 1e-2)
        end
    end
    # noise of the force: scatter of the batches
    fb = hcat([t.mom for t in tl]...)
    fm = vec(sum(fb, dims = 2)) ./ nb
    ferr = [nb > 1 ? sqrt(sum((fb[k, :] .- fm[k]) .^ 2) / (nb - 1) / nb) : 0.0 for k in 1:n] .* nb      # (error of the sum over batches, in the units of `force`)
    # pressure-tensor estimator (track length, low variance): f = −V [dP_rr/dr + (3 P_rr − u)/r]
    uE = E0 .* tot.Jpath ./ (c .* V); Prr = E0 .* tot.Prr ./ (c .* V)
    rcen = [0.5 * (Z.r[k] + Z.r[k+1]) for k in 1:n]
    forceP = zeros(n)
    for k in 1:n
        k0 = max(k - 1, 1); k1 = min(k + 1, n)
        dP = k1 > k0 ? (Prr[k1] - Prr[k0]) / (rcen[k1] - rcen[k0]) : 0.0
        forceP[k] = -V[k] * (dP + (3 * Prr[k] - uE[k]) / rcen[k])
    end
    # momentum balance of the radiation field (exact in steady state): force on zone k = A_k Π_k − A_{k+1} Π_{k+1} + ∫ (u − Π_rr)/r dV   (crossing + track-length estimators, no sign cancellation per event)
    # skipped core scatterings never cross an interface: their isotropic stress Π_c = u_c/3 at the interfaces is interpolated from the zone values of u_c (visit tallies)
    function core_if(Jsk)
        uc = E0 .* Jsk ./ (c .* V); gas = [Z.nHI[k] > 0 for k in 1:n]; Pif = zeros(n + 1)
        for i in 1:n+1
            kl = i - 1; kr = i
            if kl >= 1 && kr <= n && gas[kl] && gas[kr]
                w = (Z.r[i] - rcen[kl]) / (rcen[kr] - rcen[kl]); Pif[i] = ((1 - w) * uc[kl] + w * uc[kr]) / 3
            else
                k1 = (kr <= n && gas[kr]) ? kr : kl; k1 >= 1 && k1 <= n && gas[k1] || continue
                k2 = k1 == kr ? k1 + 1 : k1 - 1
                if k2 >= 1 && k2 <= n && gas[k2]
                    slope = (uc[k1] - uc[k2]) / (rcen[k1] - rcen[k2]); Pif[i] = max(uc[k1] + slope * (Z.r[i] - rcen[k1]), 0.0) / 3
                else
                    Pif[i] = uc[k1] / 3
                end
            end
        end
        [4π * (Z.r[k]^2 * Pif[k] - Z.r[k+1]^2 * Pif[k+1]) for k in 1:n]
    end
    Xi = tot.Xi; cif(t) = Z.dlayer > 0 ? zeros(n) : core_if(t.Jskip)        # (with a boundary layer the interface stress is carried by explicit photons: no interpolation of the skipped core stress)
    forcePi = [E0 / c * (Xi[k] - Xi[k+1] + tot.Gth[k]) for k in 1:n] .+ cif(tot)
    fb2 = hcat([[E0 / c * (t.Xi[k] - t.Xi[k+1] + t.Gth[k]) for k in 1:n] .+ cif(t) for t in tl]...)
    fm2 = vec(sum(fb2, dims = 2)) ./ nb
    ferrPi = [nb > 1 ? sqrt(sum((fb2[k, :] .- fm2[k]) .^ 2) / (nb - 1) / nb) : 0.0 for k in 1:n] .* nb
    emit = tot.emit
    (; forcePi, forcePi_err = ferrPi, forceP, uE, Prr, r = Z.r, vol = V, T = Z.T, nHI = Z.nHI, Jalpha, Jint, Rsc, Tcolor = Tc, Tcolor_err = Tcerr, Tcolor_valid = Tvalid, force = forcePi, force_err = ferrPi, force_event = force, force_event_err = ferr, nsc = tot.nsc, nskip = tot.nsk,
       heat_recoil = tot.hrec, heat_exch = tot.hexch, heat_Hm = tot.nHm .* (E0 - 0.755 * eV), E_2g_zone = tot.n2g .* E0, Egas_lab = tot.Egas, n_Hm = tot.nHm, n_H2 = tot.nH2, n_2g = tot.n2g, momsrc = tot.momsrc, Jx = tot.Jx, sources = [s.name for s in sources],
       emit, f_out = tot.nout / emit, f_in = tot.nin / emit, f_2g = sum(tot.n2g) / emit, f_Hm = sum(tot.nHm) / emit, f_H2 = sum(tot.nH2) / emit, f_lost = tot.nlost / emit,
       nsc_per_photon = tot.nscphot / emit, spec_out = vec(sum(tot.out, dims = 1)), spec_out_src = copy(tot.out), spec_in = copy(tot.inn), wmax = o.wmax, nbw = o.nbw,
       Eemit = tot.Eemit, Eout = tot.Eout, Ein = tot.Ein, E2g = tot.E2g, Eabs = tot.Eabs, Egas_total = sum(tot.Egas), batches = tl, xwin = o.xwin)
end

"The latest result: `LyaMC.result.Jalpha`, `.Tcolor`, `.force`, … (a NamedTuple behind a property forwarder; `nothing` until `run` has been called)."
mutable struct Holder; r::Any; end
Base.getproperty(h::Holder, s::Symbol) = s === :r ? getfield(h, :r) : getproperty(getfield(h, :r), s)
Base.propertynames(h::Holder) = propertynames(getfield(h, :r))
const RESULT = Holder(nothing)
const result = RESULT

"Velocity-bin centres [cm/s] of an escape spectrum."
wcenters(res) = [(-res.wmax + (i - 0.5) * 2 * res.wmax / res.nbw) for i in 1:res.nbw]


"""
    hi21_adapter(dir, tag) -> f

Adapter for `HI21.use_lyamc!(f)` (S8b): `f(snapshot)` returns `(; J, Tc)` per zone of the snapshot — J = hν_α × J_α^phot [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹], the mean intensity of the line radiation averaged
over the Doppler core |x| < 1 of the zone; Tc [K], T_k where the tilt of the line profile is not measurable — read from the `s8a_<tag>_J_<j>.csv` files that the coupled run writes at the barrier
of every snapshot (the file whose time is closest to `snapshot.t`; the zones agree: same grid). `LyaMC.result.Jalpha` / `.Tcolor` hold the live values (photon units; multiply Jalpha by `LyaMC.E0`).
"""
function hi21_adapter(dir, tag)
    function f(s)
        files = filter(x -> occursin(Regex("^s8a_$(tag)_J_\\d+\\.csv\$"), x), readdir(dir))
        isempty(files) && error("no s8a_$(tag)_J_*.csv in $dir")
        ts = [parse(Float64, strip(split(readline(joinpath(dir, x)), "=")[2])) for x in files]
        fn = files[argmin(abs.(ts .- s.t))]
        rows = [split(l, ",") for l in readlines(joinpath(dir, fn))[3:end]]
        J = [parse(Float64, r[5]) for r in rows]; Tc = [parse(Float64, r[6]) for r in rows]
        n = length(s.T); length(J) == n || error("snapshot has $n zones, $fn has $(length(J))")
        (; J, Tc)
    end
    f
end

end # module
