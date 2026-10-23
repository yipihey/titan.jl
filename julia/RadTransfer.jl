"""
    RadTransfer

S7a: frequency-dependent continuum transfer through the spherical envelope; photo-rates, photo-heating, Lyman–Werner dissociation and radiation
force per zone. Opt-in (nothing in the collapse drivers calls it). Method (notes: POPIII_NOTES.md, S7a):

- **Grid** `FreqGrid`: bins in ν (IR 0.1 eV … 100 eV), bin edges at the H⁻ (0.754 eV), Paschen/Balmer (1.51, 3.40), LW (11.2–13.598) and H I,
  He I, He II (13.598, 24.587, 54.418 eV) edges; cross sections are bin averages (3-point Gauss–Legendre), the bin-integrated source power `Source.L` [erg/s].
- **Opacity** per zone and bin from (T, n_H, abundances per H nucleus): bound–free of every `ChemNet.Photo` (level lists, Boltzmann populations of
  the Milne pair), net of stimulated emission when `stim`; recombination emission by the Milne relation (η = κ_stim 2hν³/c², `emit`); everything
  else from `ContinuumOpacity.alpha` (CIA, H⁻ ff, H/He⁺/He²⁺/He⁻/H₂⁻ ff, H₂⁺, Rayleigh, Thomson) with LTE emission κB; the bound–free terms that
  `alpha` carries for H (Kramers) and H⁻ are removed when the network's `Photo` supplies them.
- **Direct field** (central source, finite star inside the first interface): photon-conserving radial attenuation per zone and bin,
  L_out = L_in e^{−τ}; the absorbed power L_in(1 − e^{−τ}) is shared among the absorbers by their net κ; the effective mean intensity of the zone is
  J = absorbed / (4π κ V) (= L_in (1 − e^{−τ})/τ Δr/(4πV)), so that 4π∫σJ/(hν) n V is exactly the number of photons absorbed (Abel, Norman & Madau 1999;
  Mellema et al. 2006).
- **Diffuse field**: tangent-ray formal solution (zone-constant source, exact exponential across each chord segment; rays p = r_j plus core rays through the
  inner void), J and F at the interfaces from the μ-moments with piecewise-linear I(μ). Sources: thermal and recombination emission, the scattered direct
  beam; electron (and Rayleigh) scattering is coherent and isotropic and is iterated (Λ iteration, `nscat` sweeps) — not true absorption.
  The recombination continua are therefore transported (diffuse field, on-the-spot is not imposed; `AtomicHHe.ONSPOT[] = false`).
- **Outputs** (`RadField`): J, F, per-absorber photo-rates Γ, heating G (photon energy − ionisation energy), radiation acceleration (1/ρc)∫χF dν,
  k_LW with self-shielding, ledgers of the photon and energy budgets.
- **ChemNet bridge** `ZoneEnv`: per zone `env.J`, `env.kLW` and fast T-tabulated Milne rates (`PhotoQuad`) for `ChemNet.photorates_t` / `photoenergy_t`.
"""
module RadTransfer

import ..ChemNet, ..ContinuumOpacity
const CN = ChemNet; const CO = ContinuumOpacity
const hP = CN.hP; const c = CN.c; const kB = CN.kB; const eV = CN.eV
const σSB = 5.670374419e-5

# S7c-2 switches (read at construction of `Opts` / `default_grid`; env RT_LYMAN=1, RT_LYA=1 for the drivers that cannot pass them)
"Lyman-series scattering of atomic H (Kramers–Heisenberg with damping wings and thermal cores, LymanWings.jl) instead of the Rayleigh polynomial. Default off (= S7a/S7b)."
const LYMAN = Ref(get(ENV, "RT_LYMAN", "0") == "1")
"Trapped-Lyα radiation force (McKee & Tan 2008 eq. 22) as a separate component of `ARAD`. Default off."
const LYA = Ref(get(ENV, "RT_LYA", "0") in ("1", "mc"))
"S8a: RT_LYA=mc replaces the MT08 force of Lyα itself (window around the line, the injected and grid line sources) by the Monte Carlo of LyaMC.jl; the other Lyman lines keep the MT08 approximation. Default off."
const LYAMC = Ref(get(ENV, "RT_LYA", "0") == "mc")
"Half-width |ν/ν_α − 1| of the window around Lyα that the Monte Carlo transports (RT_LYA_MCWIN)."
const LYA_MC_WIN = Ref(parse(Float64, get(ENV, "RT_LYA_MCWIN", "0.12")))
"Lines Ly-2 … Ly-n resolved by `default_grid` when `lyman` (extra bin edges, see LymanWings.line_edges)."
const LYMAN_LINES = Ref(parse(Int, get(ENV, "RT_LYMAN_LINES", "12")))
"Line-centre Lyα luminosity [erg/s] injected at the inner boundary (the H II region inside r_s: 0.68 of the absorbed ionising photons); set by the driver."
const LYA_INJECT = Ref(0.0)
"Stellar radius [cm] for the blackbody limit J ≤ B of the trapped radiation (0 = no limit; used when `Source.Rstar` is 0)."
const LYA_RSTAR = Ref(0.0)

# ------------------------------------------------------------------ frequency grid
struct FreqGrid
    edges::Vector{Float64}; ν::Vector{Float64}; dν::Vector{Float64}; E::Vector{Float64}       # E: eV
end
nfreq(g::FreqGrid) = length(g.ν)
function grid_from_edges(edges)
    ν = sqrt.(edges[1:end-1] .* edges[2:end])
    FreqGrid(collect(Float64, edges), ν, diff(edges), ν .* hP ./ eV)
end
"""
    default_grid(; Emin = 0.1, Emax = 100.0, dIR = 0.1, dopt = 0.08, dLW = 0.01, dUV = 0.02)

Bins of constant Δln ν per segment between the break energies [eV] (Emin, 0.754, 1.511, 3.4, 11.2, χ_H, χ_He, χ_He⁺, Emax); ≈ 250 bins.
"""
function default_grid(; Emin = 0.1, Emax = 100.0, dIR = 0.1, dopt = 0.08, dLW = 0.01, dUV = 0.02, lyman = (LYMAN[] || LYA[]), nlines = LYMAN_LINES[])
    bp = [Emin, 0.754, 1.511, 3.4, 11.2, CN.χH / eV, CN.χHe1 / eV, CN.χHe2 / eV, Emax]
    dl = [dIR, dopt, dopt, dopt, dLW, dUV, dUV, dUV]
    E = Float64[]
    for s in 1:8
        n = max(1, ceil(Int, log(bp[s+1] / bp[s]) / dl[s])); append!(E, bp[s] .* (bp[s+1] / bp[s]) .^ ((0:n-1) ./ n))
    end
    push!(E, Emax)
    ed = E .* eV ./ hP
    lyman && (ed = merge_edges(ed, LymanWings.line_edges(2, nlines)))
    grid_from_edges(ed)
end
"Union of two ascending edge lists; an edge of `extra` closer than 2e-6 (relative) to an existing edge is dropped."
function merge_edges(ed, extra)
    out = collect(Float64, ed)
    for e in extra
        (e <= ed[1] || e >= ed[end]) && continue
        j = searchsortedfirst(out, e)
        (j <= length(out) && abs(out[j] / e - 1) < 2e-6) && continue
        (j > 1 && abs(out[j-1] / e - 1) < 2e-6) && continue
        insert!(out, j, e)
    end
    out
end
"Grid from explicit edge energies [eV]."
edge_grid(Eedges) = grid_from_edges(collect(Float64, Eedges) .* eV ./ hP)
"Index of the bin containing E [eV] (0 outside)."
function ibin(g::FreqGrid, E)
    ν = E * eV / hP; (ν < g.edges[1] || ν >= g.edges[end]) && return 0
    searchsortedlast(g.edges, ν)
end
const GL3 = ((-sqrt(3 / 5), 5 / 18), (0.0, 8 / 18), (sqrt(3 / 5), 5 / 18))     # nodes on [-1,1] (weights sum to 1/2·2 → normalised below)
"Bin average of f(ν) over bin i (3-point Gauss–Legendre, weights normalised to 1)."
function binavg(f, g::FreqGrid, i)
    a, b = g.edges[i], g.edges[i+1]; m = 0.5 * (a + b); hw = 0.5 * (b - a); s = 0.0
    for (x, w) in GL3; s += w * f(m + x * hw); end
    s / sum(w for (_, w) in GL3)
end
Bν(ν, T) = 2hP * ν^3 / c^2 / expm1(hP * ν / (kB * T))

# ------------------------------------------------------------------ sources
"Source: power per bin L [erg/s] (bin-integrated), the stellar radius [cm] (inside the first interface) and the total luminosity represented."
struct Source
    L::Vector{Float64}; Rstar::Float64; Ltot::Float64
end
"Blackbody of effective temperature Teff, luminosity L (power outside the grid is dropped; `Ltot` keeps the full L)."
function blackbody_source(g::FreqGrid, Teff, L; Rstar = sqrt(L / (4π * σSB * Teff^4)))
    Lb = zeros(nfreq(g))
    for i in 1:nfreq(g)
        Lb[i] = L * binavg(ν -> π * Bν(ν, Teff), g, i) * g.dν[i] / (σSB * Teff^4)
    end
    Source(Lb, Rstar, L)
end
"Blackbody normalised to emit Ndot photons/s above the H ionisation edge (Iliev et al. tests 2, 5)."
function blackbody_photons(g::FreqGrid, Teff, Ndot; Rstar = 0.0)
    s = blackbody_source(g, Teff, 1.0)
    N1 = sum(s.L[i] / (hP * g.ν[i]) for i in 1:nfreq(g) if g.edges[i] >= CN.χH / hP * (1 - 1e-9))
    Source(s.L .* (Ndot / N1), Rstar, sum(s.L) * Ndot / N1)
end
"Monochromatic source: Ndot photons/s of energy E [eV], placed in the bin containing E."
function mono_source(g::FreqGrid, E, Ndot; Rstar = 0.0)
    L = zeros(nfreq(g)); i = ibin(g, E); L[i] = Ndot * hP * g.ν[i]
    Source(L, Rstar, sum(L))
end

# ------------------------------------------------------------------ absorbers (ChemNet photo-processes on the grid)
struct Absorber
    p::CN.Photo
    R::Symbol; ion::Symbol; Eion::Float64          # absorber species, first product, ionisation energy E0(P) − E0(R) [erg]
    Ex::Vector{Float64}; g::Vector{Float64}
    σ::Matrix{Float64}                               # nlev × nf bin-averaged cross sections [cm²]
    Eth::Float64
end
function absorbers(photos, g::FreqGrid)
    out = Absorber[]
    for p in photos
        nl = length(p.levels); σm = zeros(nl, nfreq(g))
        for (l, lev) in enumerate(p.levels), i in 1:nfreq(g)
            g.edges[i+1] <= lev[3] / hP && continue
            σm[l, i] = binavg(lev[4], g, i)
        end
        Eion = sum(CN.SPECIES[s].E0 for s in p.P) - CN.SPECIES[p.R].E0
        push!(out, Absorber(p, p.R, p.P[1], Eion, [lev[1] for lev in p.levels], [lev[2] for lev in p.levels], σm, minimum(lev[3] for lev in p.levels)))
    end
    out
end

# ------------------------------------------------------------------ the medium and the options
"Options: `stim` stimulated emission as negative absorption, `emit` Milne recombination emission, `co` the other continua of ContinuumOpacity,
`diffuse` solve the diffuse field, `twophoton` H two-photon continuum, `nscat` scattering sweeps."
Base.@kwdef struct Opts
    stim::Bool = true; emit::Bool = true; co::Bool = true; diffuse::Bool = true; twophoton::Bool = false
    nscat::Int = 6; scat_tol::Float64 = 1e-4; conserve::Bool = true
    lyman::Bool = LYMAN[]          # S7c-2: Kramers–Heisenberg scattering of H I (Lyman wings and cores) replaces Dalgarno's Rayleigh polynomial
    lya::Bool = LYA[]              # S7c-2: trapped-Lyα radiation force (MT08 eq. 22) added to `arad` (and kept in `ARAD_LYA`)
    lya_vt::Float64 = 0.0          # extra (turbulent) speed in the Lyα Doppler width [cm/s]
    lya_mc_win::Float64 = (LYAMC[] ? LYA_MC_WIN[] : 0.0)    # S8a: > 0: bins within |ν/ν_α − 1| < lya_mc_win and the line sources are left to the Monte Carlo (LyaMC)
    lya_xh::Float64 = 3.0          # half-width of the stellar-continuum window of the trapped force, in units of x̂ = x/(aτ̄)^{1/3}
end
"Zones: interfaces r (n+1), T, nH (H nuclei cm⁻³) and X (species → abundance per H nucleus; missing = 0)."
struct Medium
    r::Vector{Float64}; T::Vector{Float64}; nH::Vector{Float64}; X::Dict{Symbol,Vector{Float64}}
end
Base.length(m::Medium) = length(m.T)
volume(r, k) = 4π / 3 * (r[k+1]^3 - r[k]^3)
xs(m::Medium, s::Symbol, k) = haskey(m.X, s) ? m.X[s][k] : 0.0

mutable struct Opacity
    kabs::Matrix{Float64}      # net absorption (≥ 0), cm⁻¹
    ssc::Matrix{Float64}       # scattering
    eta::Matrix{Float64}       # emissivity erg cm⁻³ s⁻¹ Hz⁻¹ sr⁻¹
    kbf::Vector{Matrix{Float64}}     # gross bound–free of each absorber
    knet::Vector{Matrix{Float64}}    # net (gross − stimulated)
    kother::Matrix{Float64}    # net absorption of the non-network continua
    sHI::Matrix{Float64}       # S7c-2: the part of `ssc` that is scattering by neutral H (Dalgarno's Rayleigh, or the Kramers–Heisenberg sum when `lyman`), only filled when `lyman` or `lya`
end
function new_opacity(n, nf, na)
    Opacity(zeros(n, nf), zeros(n, nf), zeros(n, nf), [zeros(n, nf) for _ in 1:na], [zeros(n, nf) for _ in 1:na], zeros(n, nf), zeros(n, nf))
end

"Lowest temperature at which ContinuumOpacity.alpha is called. 0 (default): alpha is finite and continuous from 10 K (ck/alpha-lowT: H⁻ ff held at its 1400 K validity edge, H₂⁺ LTE estimate in logs); the former 800 K floor also mis-evaluated the H₂ CIA of cold gas."
const TCO_MIN = Ref(0.0)
"ContinuumOpacity.alpha for the other continua; the bf terms the network supplies are removed (returns absorption, scattering)."
function co_other(ν, T, nH, X, abs_R)
    n(s) = nH * X(s)
    nH1 = n(:H); nHp = n(:Hp)
    a, sc, parts = CO.alpha(ν, T; nH = nH1, nH2 = n(:H2), nHe = n(:He), nHp = nHp, ne = n(:e), nHm = n(:Hm), nHep = n(:Hep), nHepp = n(:Hepp),
                             nH2p = X(:H2p) > 0 ? n(:H2p) : nothing)
    :H in abs_R && (a -= CO.α_H_kramers(ν, T, nH1, 0.0, 0.0))      # hydrogenic bf (the ff of H⁺ stays)
    :Hm in abs_R && (a -= parts.hm_bf)
    a, sc
end

"Two-photon continuum of H: ψ(y) = 12 y (1 − y) (Spitzer & Greenstein 1951), A_2γ = 8.23 s⁻¹; fraction f of the n ≥ 2 recombinations that decay through 2s."
two_photon_eta(ν, rate_B; f2γ = 0.32) = (ν0 = 10.2 * eV / hP; ν >= ν0 ? 0.0 : (y = ν / ν0; hP * ν * 2 * 12y * (1 - y) / ν0 * f2γ * rate_B / (4π)))

"Bin-averaged Lyman-series scattering table of a grid (LymanWings.table; cached on the bin edges)."
const LYTAB = Dict{UInt,Any}()
function lytable(g::FreqGrid)
    key = hash(g.edges)
    get!(LYTAB, key) do
        LymanWings.table(g.edges)
    end
end

function fill_opacity!(op::Opacity, g::FreqGrid, ab::Vector{Absorber}, med::Medium, o::Opts)
    n = length(med); nf = nfreq(g); abs_R = Symbol[a.R for a in ab]
    ltab = (o.lyman || o.lya) ? lytable(g) : nothing
    Threads.@threads for k in 1:n
        T = med.T[k]; nH = med.nH[k]; kT = kB * T
        Xk(s) = xs(med, s, k)
        for (ia, a) in enumerate(ab)
            nR = nH * Xk(a.R); nion = nH * Xk(a.ion); ne = nH * Xk(:e)
            Z = a.p.Z(T); nl = length(a.Ex)
            pop = [a.g[l] * exp(-a.Ex[l] / kT) / Z for l in 1:nl]
            lnK = CN.lnKeq([a.R], a.p.P, T)
            nn = nion * ne                                   # (two-body products; all photo pairs here have P = [ion, e])
            lnrec = nn > 0 ? log(nn) - lnK : -Inf
            for f in 1:nf
                s = 0.0
                for l in 1:nl; s += pop[l] * a.σ[l, f]; end
                s == 0 && continue
                kg = nR * s
                ks = 0.0
                if (o.stim || o.emit) && isfinite(lnrec)       # Σ_l (g_l/Z) σ_l exp(ln(n_ion n_e/K) − (E_x + hν)/kT), the exponent combined per level (no over/underflow)
                    for l in 1:nl; a.σ[l, f] > 0 && (ks += a.g[l] / Z * a.σ[l, f] * exp(lnrec - (a.Ex[l] + hP * g.ν[f]) / kT)); end
                end
                op.kbf[ia][k, f] = kg
                op.knet[ia][k, f] = kg - (o.stim ? ks : 0.0)
                o.emit && (op.eta[k, f] += ks * 2hP * g.ν[f]^3 / c^2)
            end
        end
        if o.twophoton
            # n ≥ 2 recombination rate of H (excited photo pair, Milne at J = 0): from the emission integral is not needed — use the caller's rate via X[:alpha_B]
            haskey(med.X, :aB) && for f in 1:nf; op.eta[k, f] += nH^2 * Xk(:Hp) * Xk(:e) * two_photon_eta(g.ν[f], med.X[:aB][k]); end
        end
        if o.co
            for f in 1:nf
                a_, sc = co_other(g.ν[f], max(T, TCO_MIN[]), nH, Xk, abs_R); a_ = max(a_, 0.0)
                op.kother[k, f] = a_; op.ssc[k, f] = sc
                op.eta[k, f] += a_ * Bν(g.ν[f], T)
            end
        end
        if ltab !== nothing
            # S7c-2: scattering of the neutral H atoms from the Kramers–Heisenberg sum (damping wings, thermal cores) replaces Dalgarno's polynomial of ContinuumOpacity (still
            # inside `sc` when o.co): σ̄_bin(T) is the bin average (Voigt cores resolved), not the value at the bin centre. `sHI` keeps the neutral-H part of the scattering.
            nHa = nH * Xk(:H)
            if nHa > 0
                for f in 1:nf
                    sc0 = o.co ? CO.σ_ray_H(c / g.ν[f] * 1e8) * nHa : 0.0
                    if o.lyman
                        s1 = nHa * LymanWings.sigma_bin(ltab, f, T)
                        op.ssc[k, f] = max(op.ssc[k, f] - sc0 + s1, 0.0); op.sHI[k, f] = s1
                    else
                        op.sHI[k, f] = sc0
                    end
                end
            end
        end
        for f in 1:nf
            tot = op.kother[k, f]
            for ia in eachindex(ab); tot += op.knet[ia][k, f]; end
            op.kabs[k, f] = max(tot, 0.0)
        end
    end
    op
end

# ------------------------------------------------------------------ the field
"Radiation acceleration [cm s⁻²] of the last `transfer` call (per zone): the hook for the momentum equation (S7b adds ρ a to the force of each zone)."
const ARAD = Ref(Float64[])
"Trapped-Lyα component of `ARAD` (S7c-2; zeros when `Opts.lya` is off). `ARAD[]` is the sum of the direct beam, the diffuse field and this component."
const ARAD_LYA = Ref(Float64[])
mutable struct RadField
    g::FreqGrid; r::Vector{Float64}
    Jdir::Matrix{Float64}; Jdiff::Matrix{Float64}; J::Matrix{Float64}        # n × nf zone mean intensity [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹]
    Ldir::Matrix{Float64}                                                   # (n+1) × nf direct luminosity at the interfaces [erg/s per bin]
    Fdiff::Matrix{Float64}                                                  # (n+1) × nf diffuse flux at the interfaces [erg s⁻¹ cm⁻² Hz⁻¹]
    Gamma::Vector{Vector{Float64}}                                          # per absorber: photo-rate per atom [s⁻¹]
    heat::Vector{Float64}                                                   # photo-heating of the network absorbers + others [erg cm⁻³ s⁻¹]
    heat_other::Vector{Float64}                                             # of the non-network continua
    ionrate::Vector{Float64}                                                # ionisation-energy sink rate of the network absorbers [erg cm⁻³ s⁻¹]
    arad::Vector{Float64}                                                   # radiation acceleration [cm s⁻²] (radial)
    kLW::Vector{Float64}; JLW::Vector{Float64}; fsh::Vector{Float64}
    ledger::NamedTuple
    iters::Int
end

"""
    transfer(g, ab, src, med, o::Opts; Jext = nothing, rho = nothing) -> RadField

One steady-state transfer: direct (photon-conserving) + diffuse field, photo-rates Γ, heating, radiation acceleration (needs `rho`, g cm⁻³),
LW rate (needs :H2 in `med.X`). `Jext` (length nf): isotropic intensity incident at the outer boundary.
"""
function transfer(g::FreqGrid, ab::Vector{Absorber}, src::Source, med::Medium, o::Opts = Opts(); Jext = nothing, rho = nothing, op = nothing)
    n = length(med); nf = nfreq(g); r = med.r; na = length(ab)
    op === nothing && (op = new_opacity(n, nf, na))
    fill_opacity!(op, g, ab, med, o)
    V = [volume(r, k) for k in 1:n]; Δr = diff(r)
    Ld = zeros(n + 1, nf); Jdir = zeros(n, nf); Pabs = zeros(n, nf); Psc = zeros(n, nf)
    arad_hi = zeros(n)                                    # S7c-2: the beam's force from the scattering by neutral H (single scattering), subtracted from the trapped-Lyα force
    La = o.lya ? zeros(n + 1, nf) : zeros(0, 0)          # luminosity attenuated by true absorption only (the flux of the scattered, trapped photons: ∇·F = 0)
    # ---- direct beam
    for f in 1:nf
        L = src.L[f]; Ld[1, f] = L
        if o.lya
            Lp = L; La[1, f] = L
            for k in 1:n; Lp *= exp(-op.kabs[k, f] * Δr[k]); La[k+1, f] = Lp; end
        end
        for k in 1:n
            ext = op.kabs[k, f] + op.ssc[k, f]; τ = ext * Δr[k]; t = -expm1(-τ)
            Ld[k+1, f] = L * (1 - t)
            if ext > 0
                Pabs[k, f] = L * t * op.kabs[k, f] / ext; Psc[k, f] = L * t * op.ssc[k, f] / ext
                o.lya && rho !== nothing && (arad_hi[k] += L * t * op.sHI[k, f] / ext / (c * rho[k] * V[k]))
            end
            # effective mean intensity of the zone: absorbed/(4π κ V) with absorbed = L t κ_abs/ext  →  L (t/τ) Δr/(4π V) · (κ_abs/ext)(ext/κ) = L (t/τ) Δr/(4πV) (extinction beam)
            Jdir[k, f] = (τ > 1e-300 ? L * t / τ : L) * Δr[k] / (4π * V[k]) / (g.dν[f])
            L = L * (1 - t)
        end
    end
    # ---- diffuse field
    Jdiff = zeros(n, nf); Fdiff = zeros(n + 1, nf); iters = 0
    if o.diffuse
        iters = solve_diffuse!(Jdiff, Fdiff, g, op, Psc, V, r, med, o, Jext)
    end
    Jtot = Jdir .+ Jdiff
    # ---- rates, heating
    Gam = [zeros(n) for _ in 1:na]; heat = zeros(n); heat_o = zeros(n); ionr = zeros(n)
    for k in 1:n
        T = med.T[k]; kT = kB * T
        for (ia, a) in enumerate(ab)
            nl = length(a.Ex); Z = a.p.Z(T); nR = med.nH[k] * xs(med, a.R, k)
            pop = [a.g[l] * exp(-a.Ex[l] / kT) / Z for l in 1:nl]
            γ = 0.0; Hh = 0.0
            for f in 1:nf
                s = 0.0
                for l in 1:nl; s += pop[l] * a.σ[l, f]; end
                s == 0 && continue
                x = 4π * s * Jtot[k, f] * g.dν[f] / (hP * g.ν[f])
                γ += x; Hh += x * (hP * g.ν[f] - a.Eion)
            end
            Gam[ia][k] = γ; heat[k] += nR * Hh; ionr[k] += nR * γ * a.Eion
        end
        for f in 1:nf; heat_o[k] += 4π * max(op.kother[k, f], 0.0) * Jtot[k, f] * g.dν[f]; end
        heat[k] += heat_o[k]
    end
    # ---- radiation acceleration. Direct beam: momentum removed from the beam, (L_in − L_out)/c per zone (exactly conserving); diffuse field: (1/ρc) ∫ χ F dν with
    # F at the zone centre (mean of the interface values). The radial component only (spherical symmetry).
    arad = zeros(n); arad_dir = zeros(n)
    if rho !== nothing
        for k in 1:n, f in 1:nf
            arad_dir[k] += (Ld[k, f] - Ld[k+1, f]) / (c * rho[k] * V[k])
            arad[k] += (op.kabs[k, f] + op.ssc[k, f]) * 0.5 * (Fdiff[k, f] + Fdiff[k+1, f]) * g.dν[f] / (rho[k] * c)
        end
        arad .+= arad_dir
    end
    # ---- trapped Lyα (MT08 eq. 22), a separate component
    alya = zeros(n)
    if o.lya && rho !== nothing
        GH = zeros(n); for (ia, a) in enumerate(ab); (a.R === :H && occursin("ground", a.p.name)) && (GH .+= Gam[ia]); end        # ionising-photon (> 13.6 eV, from the ground level) photoionisation rate per H atom: the stellar driver of the H II gas
        mcw = o.lya_mc_win > 0
        winbins = mcw ? lya_window(g, o.lya_mc_win) : nothing
        alya = lya_accel(g, lytable(g), La, med, rho, o; Linject = mcw ? 0.0 : LYA_INJECT[], Rstar = (src.Rstar > 0 ? src.Rstar : LYA_RSTAR[]), Lsrc = src.L,
                         GammaH = mcw ? nothing : GH, skipbins = winbins, nolines = mcw) .- arad_hi
        arad .+= alya
    end
    ARAD_LYA[] = alya
    # ---- ledger (photons above χ_H and energy)
    iH = findfirst(i -> g.edges[i] >= CN.χH / hP * (1 - 1e-9), 1:nf)
    iH === nothing && (iH = nf + 1)          # (grids without ionising bins)
    Nsrc = sum(src.L[f] / (hP * g.ν[f]) for f in iH:nf; init = 0.0); Ndir_out = sum(Ld[n+1, f] / (hP * g.ν[f]) for f in iH:nf; init = 0.0)
    Ndif_out = sum(4π * r[n+1]^2 * Fdiff[n+1, f] * g.dν[f] / (hP * g.ν[f]) for f in iH:nf; init = 0.0)
    Nabs = 0.0; Nem = 0.0
    for k in 1:n
        Nem += sum(4π * op.eta[k, f] * V[k] * g.dν[f] / (hP * g.ν[f]) for f in iH:nf; init = 0.0)         # 4π ∫ η dV
        Nabs += sum(4π * op.kabs[k, f] * Jtot[k, f] * V[k] * g.dν[f] / (hP * g.ν[f]) for f in iH:nf; init = 0.0)
    end
    Eabs = sum(4π * op.kabs[k, f] * Jtot[k, f] * V[k] * g.dν[f] for k in 1:n, f in 1:nf)
    ledger = (Nsrc = Nsrc, Ndir_out = Ndir_out, Ndif_out = Ndif_out, Nabs = Nabs, Nem = Nem,
              Eabs = Eabs, Eheat = sum(heat .* V) + sum(ionr .* V), Esrc = sum(src.L), Edir_out = sum(Ld[n+1, :]))
    # ---- LW
    kLW = zeros(n); JLW = zeros(n); fsh = zeros(n)
    if haskey(med.X, :H2)
        kLW, JLW, fsh = lw_rates(g, med, Jdir, Jdiff, Jext)
    end
    ARAD[] = arad
    RadField(g, r, Jdir, Jdiff, Jtot, Ld, Fdiff, Gam, heat, heat_o, ionr, arad, kLW, JLW, fsh, ledger, iters)
end

# ------------------------------------------------------------------ diffuse field: tangent rays
"""
Solve the diffuse field. Fills Jdiff (zone means), Fdiff (interfaces). Sources per zone and bin: η (thermal/recombination) + the scattered direct beam
(Psc/(4π V dν)) + σ_s J_diff (iterated).
"""
function solve_diffuse!(Jdiff, Fdiff, g, op, Psc, V, r, med, o, Jext)
    n = size(Jdiff, 1); nf = size(Jdiff, 2)
    # ray impact parameters: core rays through the inner void and p = r_j (j = 1 … n+1)
    core = r[1] > 0 ? r[1] .* [0.0, 0.35, 0.6, 0.8, 0.93] : Float64[]
    p = vcat(core, r[1:n+1]); np = length(p)
    kmin = [min(p_ < r[1] ? 1 : findlast(k -> r[k] <= p_ * (1 + 1e-12), 1:n+1), n) for p_ in p]
    ds = zeros(np, n)
    for j in 1:np, k in kmin[j]:n
        ds[j, k] = sqrt(max(r[k+1]^2 - p[j]^2, 0.0)) - sqrt(max(r[k]^2 - p[j]^2, 0.0))
    end
    nj_of = [count(pj -> pj <= r[k] * (1 + 1e-12), p) for k in 1:n+1]
    sca_any = any(>(0), op.ssc)
    Jprev = zeros(n, nf); iters = 0
    jem0 = [op.eta[k, f] + Psc[k, f] / (4π * V[k] * g.dν[f]) for k in 1:n, f in 1:nf]       # isotropic emissivity per sr and Hz
    for sweep in 1:(sca_any ? o.nscat : 1)
        JI = zeros(n + 1, nf)
        Threads.@threads for f in 1:nf
            Iin = zeros(np, n + 1); Iout = zeros(np, n + 1)
            χ = [op.kabs[k, f] + op.ssc[k, f] for k in 1:n]
            S = [χ[k] > 0 ? (jem0[k, f] + op.ssc[k, f] * Jprev[k, f]) / χ[k] : 0.0 for k in 1:n]
            Ib = Jext === nothing ? 0.0 : Jext[f]
            for j in 1:np
                I = Ib; Iin[j, n+1] = I
                for k in n:-1:kmin[j]
                    d = ds[j, k]; τ = χ[k] * d
                    I = χ[k] > 0 ? I * exp(-τ) + S[k] * (-expm1(-τ)) : I + jem0[k, f] * d
                    Iin[j, k] = I
                end
                Iout[j, kmin[j]] = Iin[j, kmin[j]]
                for k in kmin[j]:n
                    d = ds[j, k]; τ = χ[k] * d
                    I = χ[k] > 0 ? I * exp(-τ) + S[k] * (-expm1(-τ)) : I + jem0[k, f] * d
                    Iout[j, k+1] = I
                end
            end
            for k in 1:n+1
                rk = r[k]; rk <= 0 && continue
                nj = nj_of[k]; Jk = 0.0; Fk = 0.0
                for j in 1:nj-1
                    μa = sqrt(max(1 - (p[j+1] / rk)^2, 0.0)); μb = sqrt(max(1 - (p[j] / rk)^2, 0.0)); dμ = μb - μa
                    dμ <= 0 && continue
                    ib, ia = Iin[j, k], Iin[j+1, k]; ob, oa = Iout[j, k], Iout[j+1, k]
                    Jk += 0.25 * dμ * ((ib + ob) + (ia + oa))
                    Fk += dμ / 6 * ((ob - ib) * (2μb + μa) + (oa - ia) * (μb + 2μa))
                end
                Fdiff[k, f] = 2π * Fk; JI[k, f] = Jk
            end
        end
        Jnew = zeros(n, nf)
        for k in 1:n, f in 1:nf
            j1 = (k == 1 && r[1] <= 0) ? JI[2, f] : JI[k, f]
            Jnew[k, f] = 0.5 * (j1 + JI[k+1, f])
        end
        # Zone-mean intensity for the absorption. The mean of the two interface values is accurate in thin zones (τ ≲ 0.5) but not in thick ones (J falls
        # across the zone); there the absorbed power follows from the flux divergence, A = 4π V η_tot + L_k − L_{k+1}, J_eff = A/(4π κ V). The thin zones are
        # then scaled per bin so that the total absorbed power closes the budget exactly: Σ A = Σ 4πVη_tot − L_out + L_in(1) (conservation to round-off).
        if o.conserve
            for f in 1:nf
                Aspec = 0.0; Athin = 0.0; Athick = 0.0; thin = Int[]
                for k in 1:n
                    Ek = 4π * V[k] * jem0[k, f]; Aspec += Ek
                    Ak = Ek + 4π * r[k]^2 * Fdiff[k, f] - 4π * r[k+1]^2 * Fdiff[k+1, f]
                    kk = op.kabs[k, f]
                    if kk > 0 && (kk + op.ssc[k, f]) * (r[k+1] - r[k]) > 0.5
                        Jnew[k, f] = max(Ak, 0.0) / (4π * kk * V[k]); Athick += 4π * kk * Jnew[k, f] * V[k]
                    elseif kk > 0
                        push!(thin, k); Athin += 4π * kk * Jnew[k, f] * V[k]
                    end
                end
                target = Aspec - 4π * r[n+1]^2 * Fdiff[n+1, f] + 4π * r[1]^2 * Fdiff[1, f] - Athick
                if Athin > 0 && target > 0
                    sc = target / Athin
                    0.2 < sc < 5.0 && for k in thin; Jnew[k, f] *= sc; end
                end
            end
        end
        iters = sweep
        dmax = 0.0
        for i in eachindex(Jnew); Jnew[i] > 0 && (dmax = max(dmax, abs(Jnew[i] - Jprev[i]) / Jnew[i])); end
        Jprev = Jnew
        sca_any && dmax < o.scat_tol && break
    end
    Jdiff .= Jprev
    iters
end

# ------------------------------------------------------------------ Lyman–Werner
"Doppler parameter of H₂ [cm/s] (thermal)."
bdopp(T) = sqrt(2kB * T / (2 * CN.mH))
"Self-shielding factor, Draine & Bertoldi (1996) eq. 37 (α = 2) or Wolcott-Green, Haiman & Bryan (2011) eq. 10 (α = 1.1, default)."
function fsh_lw(NH2, T; α = 1.1)
    x = NH2 / 5e14; b5 = bdopp(T) / 1e5
    0.965 / (1 + x / b5)^α + 0.035 / sqrt(1 + x) * exp(-8.5e-4 * sqrt(1 + x))
end
const KLW_COEF = 1.38e9        # s⁻¹ per (erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹) of J_LW (Abel et al. 1997; Wolcott-Green et al. 2011)
const ELW = 12.87              # eV
const LWSHIELD = Ref(:wghb11)
"""
    lw_rates(g, med, Jdir, Jdiff, Jext) -> (k_LW, J_LW, f_sh)

k_LW = 1.38e9 [J_dir(12.87 eV) f_sh(N_in) + (J_diff + J_ext)(12.87 eV) f_sh(N_out)] per H₂ molecule, N_in the H₂ column from the source to the zone centre
(the direct beam), N_out from the zone centre to the surface (diffuse and external field); returns the combined unshielded J_LW and the effective f_sh.
"""
function lw_rates(g, med, Jdir, Jdiff, Jext; α = LWSHIELD[] === :db96 ? 2.0 : 1.1)
    n = length(med); r = med.r; i = ibin(g, ELW)
    nH2 = [med.nH[k] * xs(med, :H2, k) for k in 1:n]; dr = diff(r)
    Nin = zeros(n); Nout = zeros(n); acc = 0.0
    for k in 1:n; Nin[k] = acc + 0.5 * nH2[k] * dr[k]; acc += nH2[k] * dr[k]; end
    tot = acc; acc = 0.0
    for k in n:-1:1; Nout[k] = acc + 0.5 * nH2[k] * dr[k]; acc += nH2[k] * dr[k]; end
    kLW = zeros(n); JLW = zeros(n); fs = zeros(n)
    for k in 1:n
        Jd = Jdir[k, i]; Jo = Jdiff[k, i] + (Jext === nothing ? 0.0 : Jext[i])
        fi = fsh_lw(Nin[k], med.T[k]; α); fo = fsh_lw(Nout[k], med.T[k]; α)
        kLW[k] = KLW_COEF * (Jd * fi + Jo * fo); JLW[k] = Jd + Jo
        fs[k] = JLW[k] > 0 ? (Jd * fi + Jo * fo) / JLW[k] : 1.0
    end
    kLW, JLW, fs
end
"The LW photodissociation reaction H₂ + γ → 2H for a ChemNet network (rate from `env.kLW`, 0 if the env has none). Opt-in: `with_lw(net)`."
lw_reaction() = CN.Reaction("H2 + γ(LW) → 2H (k = 1.38e9 J_LW f_sh)", [:H2], [:H, :H], (T, env) -> hasproperty(env, :kLW) ? env.kLW : 0.0)
with_lw(net::CN.Network) = CN.Network(net.species, net.free, vcat(net.reactions, [lw_reaction()]))

# ------------------------------------------------------------------ ChemNet bridge
"""
Quadrature of a Photo on fixed ν-nodes (as `ChemNet.photorates`, with the node range for T up to `Tmax`) with the cross sections tabulated once;
per zone `prepare!` fills the J values at the nodes, `eval` returns (R, A, RE, AE, AX) at T (the Milne sums of ChemNet.photorates / photoenergy).
"""
struct PhotoQuad
    p::CN.Photo
    ν::Vector{Vector{Float64}}; q::Vector{Vector{Float64}}        # per level: node frequencies and 4πσ dν/(hν)
    Ex::Vector{Float64}; g::Vector{Float64}
    J::Vector{Vector{Float64}}                                      # J at the nodes (filled by prepare!)
    Rl::Vector{Float64}; REl::Vector{Float64}                       # per level: 4π∫σJ/(hν) dν and 4π∫σJ dν (T independent)
end
function PhotoQuad(p::CN.Photo; n = 240, Tmax = 1e6)
    νs = Vector{Float64}[]; qs = Vector{Float64}[]
    for (Ex, gl, Eth, σf) in p.levels
        ν0 = Eth / hP; edges = [ν0 * exp(x) for x in range(0.0, log(1 + 60kB * Tmax / Eth + 30.0), length = n + 1)]
        ν = 0.5 .* (edges[1:end-1] .+ edges[2:end]); dν = diff(edges)
        push!(νs, ν); push!(qs, [4π * σf(ν[i]) / (hP * ν[i]) * dν[i] for i in 1:n])
    end
    nl = length(νs)
    PhotoQuad(p, νs, qs, [l[1] for l in p.levels], [l[2] for l in p.levels], [zeros(length(v)) for v in νs], zeros(nl), zeros(nl))
end
"""
Per-zone env of ChemNet (`rates_fn` for `photorates_t`, `energy_fn` for `photoenergy_t`, `J`, `kLW`). Photo-rates R and absorbed energies RE come from the
transfer's own bins (bin-average cross sections × J_f dν_f: identical to the Γ of `transfer`, hence photon conserving) in `:binned` mode, or from the
quadrature nodes in `:thermal` mode (J = B analytic, where R and the Milne coefficient must share their nodes for exact detailed balance). The Milne
coefficient and the emitted energies (A, AE, AX) are quadratures on fixed nodes, tabulated in log T (0.02 dex, Catmull–Rom, filled lazily per zone) unless `exact`.
"""
mutable struct ZoneEnv
    Trad::Float64
    J::Function
    static::Bool
    kLW::Float64
    rates_fn::Function
    energy_fn::Function
    edges::Vector{Float64}; tab::Vector{Float64}
    quads::IdDict{Any,PhotoQuad}
    cache::Dict{Tuple{UInt,Int},NTuple{3,Float64}}
    exact::Bool
    mode::Symbol
    absb::IdDict{Any,Absorber}
    g::FreqGrid
end
function ZoneEnv(g::FreqGrid, photos; exact = false, Trad = 0.0)
    ze = ZoneEnv(Trad, ν -> 0.0, false, 0.0, (p, T) -> error("unset"), (p, T) -> error("unset"), g.edges, zeros(nfreq(g)),
                 IdDict{Any,PhotoQuad}(), Dict{Tuple{UInt,Int},NTuple{3,Float64}}(), exact, :binned, IdDict{Any,Absorber}(), g)
    for (p, a) in zip(photos, absorbers(photos, g)); ze.quads[p] = PhotoQuad(p); ze.absb[p] = a; end
    ze.rates_fn = (p, T) -> (v = quad_eval(ze, p, T); (v[1], v[2]))
    ze.energy_fn = (p, T) -> (v = quad_eval(ze, p, T); (v[3], v[4], v[5]))
    ze.J = ν -> 0.0
    ze
end
"Refresh the per-level photo-integrals from the current J (after `ze.J` or the table changed) and clear the T cache."
function prepare!(ze::ZoneEnv)
    empty!(ze.cache)
    g = ze.g
    for (p, pq) in ze.quads
        a = ze.absb[p]
        for l in eachindex(pq.ν)
            if ze.mode === :thermal
                @inbounds for i in eachindex(pq.ν[l]); pq.J[l][i] = ze.J(pq.ν[l][i]); end
                pq.Rl[l] = sum(pq.q[l][i] * pq.J[l][i] for i in eachindex(pq.ν[l]))
                pq.REl[l] = sum(pq.q[l][i] * pq.J[l][i] * hP * pq.ν[l][i] for i in eachindex(pq.ν[l]))
            else
                @inbounds for i in eachindex(pq.ν[l]); pq.J[l][i] = ze.J(pq.ν[l][i]); end
                r = 0.0; re = 0.0
                for f in 1:nfreq(g)
                    s = a.σ[l, f]; s == 0 && continue
                    r += 4π * s * ze.tab[f] * g.dν[f] / (hP * g.ν[f]); re += 4π * s * ze.tab[f] * g.dν[f]
                end
                pq.Rl[l] = r; pq.REl[l] = re
            end
        end
    end
    ze
end
"Set the zone: bin-wise intensity (piecewise constant) and the LW rate."
function set_zone!(ze::ZoneEnv, Jbins::AbstractVector, kLW = 0.0)
    ze.tab .= Jbins; ze.kLW = kLW; ze.mode = :binned
    edges = ze.edges; tab = ze.tab; nf = length(tab)
    ze.J = ν -> (i = searchsortedlast(edges, ν); (i >= 1 && i <= nf) ? tab[i] : 0.0)
    prepare!(ze)
end
"""
    set_zone_light!(ze, Jbins, kLW = 0.0)

S7b: as `set_zone!` but only the photo-integrals R and RE (bin sums, T independent) are refreshed; the nodes of the Milne quadratures keep their J (set by a previous `set_zone!`, e.g. J = 0 or the CMB)
and the T-tabulated Milne coefficients (the expensive part) are kept. Valid where J_ν ≪ 2hν³/c² (stimulated recombination negligible), which holds everywhere outside the immediate neighbourhood of a hot star
(J/(2hν³/c²) ~ 1e-9 at 1e16 cm from 1e6 L☉). Cost: one pass over the bins per level.
"""
function set_zone_light!(ze::ZoneEnv, Jbins::AbstractVector, kLW = 0.0)
    ze.tab .= Jbins; ze.kLW = kLW; ze.mode = :binned
    edges = ze.edges; tab = ze.tab; nf = length(tab)
    ze.J = ν -> (i = searchsortedlast(edges, ν); (i >= 1 && i <= nf) ? tab[i] : 0.0)
    g = ze.g
    for (p, pq) in ze.quads
        a = ze.absb[p]
        for l in eachindex(pq.ν)
            r = 0.0; re = 0.0
            @inbounds for f in 1:nfreq(g)
                s = a.σ[l, f]; s == 0 && continue
                r += 4π * s * tab[f] * g.dν[f] / (hP * g.ν[f]); re += 4π * s * tab[f] * g.dν[f]
            end
            pq.Rl[l] = r; pq.REl[l] = re
        end
    end
    ze
end

"Set J = B_ν(T) exactly (analytic, not binned): the detailed-balance state."
function set_thermal!(ze::ZoneEnv, T)
    ze.J = ν -> Bν(ν, T); ze.kLW = 0.0; ze.mode = :thermal
    prepare!(ze)
end
"A, AE, AX of the Milne sums at T (quadrature nodes)."
function quad_milne(pq::PhotoQuad, T)
    p = pq.p; lnK = CN.lnKeq([p.R], p.P, T); Z = p.Z(T); kT = kB * T
    A = 0.0; AE = 0.0; AX = 0.0
    for l in eachindex(pq.ν)
        w = pq.g[l] / Z; Ex = pq.Ex[l]; ν = pq.ν[l]; q = pq.q[l]; J = pq.J[l]
        Al = 0.0; Ae = 0.0; pre = -(Ex) / kT - lnK
        for i in eachindex(ν)
            a = q[i] * (2hP * ν[i]^3 / c^2 + J[i]) * exp(pre - hP * ν[i] / kT)
            Al += a; Ae += a * hP * ν[i]
        end
        A += w * Al; AE += w * Ae; AX += w * Ex * Al
    end
    (A, AE, AX)
end
const DLT = 0.02
function quad_eval(ze::ZoneEnv, p, T)
    pq = ze.quads[p]; Tv = float(T); Z = p.Z(Tv); kT = kB * Tv
    R = 0.0; RE = 0.0
    for l in eachindex(pq.ν)
        wb = pq.g[l] / Z * exp(-pq.Ex[l] / kT); R += wb * pq.Rl[l]; RE += wb * pq.REl[l]
    end
    if ze.exact
        m = quad_milne(pq, Tv); return (R, m[1], RE, m[2], m[3])
    end
    x = clamp(log10(Tv), 0.0, 6.0); i = clamp(floor(Int, x / DLT) + 1, 2, round(Int, 6.0 / DLT) - 2); w = (x - (i - 1) * DLT) / DLT
    node(j) = get!(ze.cache, (objectid(p), j)) do; v = quad_milne(pq, 10^((j - 1) * DLT)); ntuple(m -> log(max(v[m], 1e-300)), 3); end
    nd = (node(i - 1), node(i), node(i + 1), node(i + 2))
    out = ntuple(m -> (p0 = nd[1][m]; p1 = nd[2][m]; p2 = nd[3][m]; p3 = nd[4][m]; exp(p1 + 0.5w * (p2 - p0 + w * (2p0 - 5p1 + 4p2 - p3 + w * (3(p1 - p2) + p3 - p0))))), 3)
    o = ntuple(m -> out[m] < 1e-290 ? 0.0 : out[m], 3)
    (R, o[1], RE, o[2], o[3])
end

# ------------------------------------------------------------------ thin helpers for drivers
"""
Zone hook for `PopIIIMaterial.NetGas.zonehook`-style drivers: `zonehook(ze, rf)` returns `(k, d) -> …` that loads zone k's J_ν and k_LW of the field into the env
(call it before the zone's `ChemNet.step`; with it the network, the photo-heating in `AtomicHHe.cooling(photos, ze)` and the LW reaction all see this zone's field).
Index offset: `k0` = index of the first real zone of the driver minus 1.
"""
zonehook(ze::ZoneEnv, rf; k0 = 0) = (k, d) -> (set_zone!(ze, view(rf.J, k - k0, :), isempty(rf.kLW) ? 0.0 : rf.kLW[k - k0]); nothing)

"Abundance dictionary (per H nucleus) of every zone from a ChemNet network state: x[k] free vector → Dict of vectors."
function medium_from_net(r, T, nH, xs_, net)
    n = length(T); X = Dict{Symbol,Vector{Float64}}()
    for k in 1:n
        d = CN.full(xs_[k], net)
        for (s, v) in d; get!(X, s, zeros(n))[k] = v; end
    end
    Medium(r, T, nH, X)
end

include(joinpath(@__DIR__, "RadSweep.jl"))      # S7b: ion_sweep, interface_accel
include(joinpath(@__DIR__, "LymanWings.jl"))    # S7c-2: H I Lyman-series scattering cross section, MT08 trapped-Lyα pressure
LymanWings.LYB_CONV[] = LYAMC[]
include(joinpath(@__DIR__, "RadLya.jl"))        # S7c-2: trapped-Lyα radiation force on the grid

end # module
