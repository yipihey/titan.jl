"""
    ChemNet

Population III project, G3: a general reaction-network engine for metal-free gas (replaces the hand-coded network of
PrimordialChem for new work; PrimordialChem's RHFC02 network stays as the reference it must reproduce).

- **Species** carry mass, charge, element content, an energy zero E0 (relative to free nuclei and electrons at rest —
  only differences within an element matter) and an internal partition function Z_int(T) (electronic × nuclear-spin ×
  rotational × vibrational; nuclear spins kept consistently so that equilibrium constants are exact).
- **Equilibrium constants** for any reaction from the species: K = Π_P (Z_tr Z_int e^{−E0/kT}) / Π_R (…), Z_tr =
  (2π m k T/h²)^{3/2}.
- **Reactions** as data: reactants, products, rate function k(T, env) (env: n_H, J_ν access, …). `db = true` adds the
  reverse reaction with k_rev = k_fwd / K (detailed balance), so the network's own steady state is the Saha/K equilibrium.
- **Photo-processes** from cross sections σ(ν) and the mean intensity J_ν supplied per zone (`env.J(ν)`; default: a
  blackbody at `env.Trad`): photo-rate 4π∫σJ/(hν)dν; the inverse (radiative attachment/recombination) from the Milne
  relation with stimulated emission, α = (n_A/(n_B n_e))_eq 4π∫σ/(hν)(2hν³/c² + J_ν)e^{−hν/kT}dν — so J = B gives the
  equilibrium exactly.
- **Unknowns**: number fractions x_i = n_i/n_H of all species except the eliminated ones — one per element (H, He, Li:
  from nuclei conservation) and e⁻ (charge). Conservation holds to round-off by construction.
- **Energy**: u = Σ x_i (c_v,i T + E0_i)/μ_H per gram (c_v: 3/2 k for atoms/ions/e⁻, RHFC02 A9 for H₂, 5/2 k for H₂⁺) —
  the generalisation of PrimordialChem.eos_rhfc02.
- **Step**: backward Euler in x with the rates at the new T plus u(ρ, T, x) = u − dt Λ/ρ (chemistry and T together),
  Newton with the exact Jacobian (ForwardDiff), damping that keeps every abundance (incl. the eliminated ones) ≥ 0,
  recursive halving.
"""
module ChemNet

using ForwardDiff, Printf

const kB = 1.380649e-16; const hP = 6.62607015e-27; const c = 2.99792458e10; const eV = 1.602176634e-12
const mH = 1.6735575e-24; const me = 9.1093837e-28; const mp = 1.67262192e-24
const yHe = 1 / 12; const yLi = 4.15e-10

# ------------------------------------------------------------------ species
struct Species
    name::Symbol
    m::Float64            # mass [g]
    q::Int                # charge
    el::NTuple{3,Int}     # (H, He, Li) nuclei
    E0::Float64           # energy zero [erg]
    Zint::Function        # T -> internal partition function (incl. nuclear spin)
    cv::Function          # T -> heat capacity per particle in units of k (thermal energy = ∫ c_v dT ≈ c_v T, see u_th)
end
const θrot_H2 = 85.3; const θvib_H2 = 6332.0
function Zrot_H2(T)
    Zr = zero(T); for J in 0:60; Zr += (iseven(J) ? 1 : 3) * (2J + 1) * exp(-J * (J + 1) * θrot_H2 / T); end; Zr
end
Zvib(T, θ) = 1 / (1 - exp(-θ / T))
"RHFC02 A9 (ON98): 1/(γ_H₂ − 1) = 5/2 + x² eˣ/(eˣ − 1)², x = 6100 K/T — energy per H₂ (thermal) = this × kT."
u_H2(T) = (x = 6100 / T; 5 / 2 + x^2 * exp(x) / expm1(x)^2)
const χH = 13.598eV; const EA_H = 0.754eV; const D0_H2 = 4.478eV; const D0_H2p = 2.65eV
const χLi = 5.3917eV
"Levels of atomic H kept in its partition function and photo-ionisation (n ≤ NH_LEVELS; the sum diverges without a
cut-off — at the densities where excited levels matter they are pressure-dissolved and SCVH95 takes over)."
const NH_LEVELS = 8
"Electronic partition function of H, levels n ≤ NH_LEVELS (g = 2n², E_n = χ(1 − 1/n²))."
ZH_el(T) = (z = zero(T); for n in 1:NH_LEVELS; z += 2n^2 * exp(-χH * (1 - 1 / n^2) / (kB * T)); end; z)
"H₂⁺ rotation (B = 29.9 cm⁻¹ → θ = 43 K) and vibration (2297 cm⁻¹ → 3305 K), ortho/para weights 3/1, electron spin 2 [approximate]."
Zint_H2p(T) = (Zr = zero(T); for J in 0:60; Zr += (iseven(J) ? 1 : 3) * (2J + 1) * exp(-J * (J + 1) * 43.0 / T); end; 2Zr * Zvib(T, 3305.0))

const SPECIES = Dict{Symbol,Species}(
    # nuclear spin: H nucleus 2 states; Li-7 nucleus 4; He-4 1 (only ratios within an element matter)
    :H   => Species(:H, mH, 0, (1, 0, 0), -χH, T -> 2 * ZH_el(T), T -> 1.5),     # × nuclear spin 2
    :Hp  => Species(:Hp, mp, 1, (1, 0, 0), 0.0, T -> 2 + zero(T), T -> 1.5),
    :Hm  => Species(:Hm, mH + me, -1, (1, 0, 0), -(χH + EA_H), T -> 2 + zero(T), T -> 1.5),
    :H2  => Species(:H2, 2mH, 0, (2, 0, 0), -(2χH + D0_H2), T -> Zrot_H2(T) * Zvib(T, θvib_H2), T -> NaN),   # c_v via u_H2
    :H2p => Species(:H2p, 2mH - me, 1, (2, 0, 0), -(χH + D0_H2p), Zint_H2p, T -> 2.5),
    :e   => Species(:e, me, -1, (0, 0, 0), 0.0, T -> 2 + zero(T), T -> 1.5),
    :He  => Species(:He, 4.0026 * 1.66053907e-24, 0, (0, 1, 0), 0.0, T -> 1 + zero(T), T -> 1.5),
    :Li  => Species(:Li, 7.016 * 1.66053907e-24, 0, (0, 0, 1), -χLi, T -> 2 * 4 + zero(T), T -> 1.5),
    :Lip => Species(:Lip, 7.016 * 1.66053907e-24 - me, 1, (0, 0, 1), 0.0, T -> 4 + zero(T), T -> 1.5),
)
Ztr(s::Species, T) = (2π * s.m * kB * T / hP^2)^1.5
"Equilibrium constant Π_P n / Π_R n [cm^{-3 Δν}] of R → P (vectors of species names, repeats allowed)."
function Keq(R, P, T)
    lnK = zero(T)
    for p in P; s = SPECIES[p]; lnK += log(Ztr(s, T) * s.Zint(T)) - s.E0 / (kB * T); end
    for r in R; s = SPECIES[r]; lnK -= log(Ztr(s, T) * s.Zint(T)) - s.E0 / (kB * T); end
    exp(lnK)
end

# ------------------------------------------------------------------ reactions
struct Reaction
    name::String
    R::Vector{Symbol}; P::Vector{Symbol}
    k::Function           # (T, env) -> rate coefficient [cm^{3(n_R−1)} s⁻¹]
end
"Photo-process X + γ → products with cross section σ(ν, T) (T: the species' level populations, e.g. atomic H's excited
levels in Boltzmann ratio) above the threshold energy Eth."
struct Photo
    name::String
    R::Symbol; P::Vector{Symbol}; σ::Function; Eth::Float64
    milne::Bool           # add the inverse (radiative attachment/recombination) from the Milne relation
end
struct Network
    species::Vector{Symbol}          # all species in the network (incl. eliminated; He as inert spectator)
    free::Vector{Symbol}             # the unknowns
    reactions::Vector{Reaction}
end

"Reaction list with detailed-balance reverses added for the entries flagged `db`."
function with_db(list)
    out = Reaction[]
    for (r, db) in list
        push!(out, r)
        db && push!(out, Reaction(r.name * " (reverse, detailed balance)", r.P, r.R, (T, env) -> r.k(T, env) / Keq(r.R, r.P, T)))
    end
    out
end

"Photo-rate [s⁻¹] 4π∫σJ/(hν)dν and the Milne inverse coefficient α for A + e → X + γ (stimulated emission included)."
function photorates(p::Photo, T, env; n = 400)
    ν0 = p.Eth / hP; νs = [ν0 * exp(x) for x in range(0.0, log(1 + 60kB * max(T, env.Trad) / p.Eth + 30.0), length = n + 1)]
    R = 0.0; A = 0.0
    for i in 1:n
        ν = 0.5 * (νs[i] + νs[i+1]); dν = νs[i+1] - νs[i]; σ = p.σ(ν, T); J = env.J(ν)
        R += 4π * σ * J / (hP * ν) * dν
        A += 4π * σ / (hP * ν) * (2hP * ν^3 / c^2 + J) * exp(-hP * ν / (kB * T)) * dν
    end
    # α = (n_X/(n_products))_eq × A, (n_X / Π n_P)_eq = 1/K(X → P)
    R, A / Keq([p.R], p.P, T)
end

# ------------------------------------------------------------------ composition
"Eliminated species: H (H nuclei), Li (Li nuclei), e (charge); He fixed at y_He."
function full(x, net)
    d = Dict{Symbol,eltype(x)}()
    for (i, s) in enumerate(net.free); d[s] = x[i]; end
    get!(d, :He, yHe + zero(eltype(x)))
    hasLi = :Li in net.species
    H = 1 - sum((SPECIES[s].el[1] * d[s] for s in net.free if s != :H); init = zero(eltype(x)))
    d[:H] = H
    hasLi && (d[:Li] = yLi - sum((SPECIES[s].el[3] * d[s] for s in net.free if SPECIES[s].el[3] > 0); init = zero(eltype(x))))
    d[:e] = sum((SPECIES[s].q * d[s] for s in keys(d) if s != :e); init = zero(eltype(x)))
    d
end
μH() = mH * (1 + 4yHe)    # mass per H nucleus (Li negligible)

"Internal energy per gram and pressure."
function eos(ρ, T, x, net)
    d = full(x, net); nH = ρ / μH()
    u = zero(T); N = zero(T)
    for (s, xs) in d
        sp = SPECIES[s]
        u += xs * ((s == :H2 ? u_H2(T) : sp.cv(T)) * kB * T + sp.E0)
        N += xs
    end
    (u / μH(), N * nH * kB * T)
end

"dx/dt for the free species at H-nucleus density nH."
function dxdt(x, nH, T, net, env, photo)
    d = full(x, net); dd = Dict(s => zero(eltype(x)) for s in keys(d))
    for r in net.reactions
        k = r.k(T, env); rate = k * nH^(length(r.R) - 1)
        for s in r.R; rate *= d[s]; end
        for s in r.R; dd[s] -= rate; end
        for s in r.P; dd[s] += rate; end
    end
    for (p, R, α) in photo                          # precomputed per zone: photo-rate R, inverse coefficient α
        rate = R * d[p.R]; for s in p.P; dd[s] += rate; end; dd[p.R] -= rate
        inv_ = α * nH^(length(p.P) - 1); for s in p.P; inv_ *= d[s]; end    # per H nucleus: α n_H^(k−1) Π x
        for s in p.P; dd[s] -= inv_; end; dd[p.R] += inv_
    end
    [dd[s] for s in net.free]
end

"""
    step(x0, ρ, u, T0, dt, net; env, photos = Photo[], cooling = nothing) -> (x, T)

Chemistry and temperature together over dt at fixed ρ and u (u − dt Λ/ρ with `cooling(ρ, T, d)`, d the full composition).
"""
function step(x0, ρ, u, T0, dt, net; env, photos = Photo[], cooling = nothing, tol = 1e-10, maxit = 60)
    nH = ρ / μH(); nf = length(net.free)
    function resid(z, x0, h)
        x = z[1:nf]; T = exp(z[nf+1])
        ph = [(p, photorates(p, ForwardDiff.value(T), env)...) for p in photos]   # photo-rates lagged in T (Float)
        f = dxdt(x, nH, T, net, env, ph)
        r = x .- x0 .- h .* f
        ue = eos(ρ, T, x, net)[1] - u + (cooling === nothing ? zero(T) : h * cooling(ρ, T, full(x, net)) / ρ)
        vcat(r, ue / (kB * T / μH()))
    end
    function attempt(x0, T0, h)
        z = vcat(collect(Float64, x0), log(T0))
        for it in 1:maxit
            R = resid(z, x0, h)
            Jm = ForwardDiff.jacobian(zz -> resid(zz, x0, h), z)
            δ = -(Jm \ R); all(isfinite, δ) || return (z, false)
            λ = 1.0
            # keep every abundance ≥ 0, incl. the eliminated H, Li, e (linear in x)
            for i in 1:nf; z[i] + λ * δ[i] < 0 && (λ = min(λ, 0.9 * z[i] / -δ[i])); end
            d0 = full(z[1:nf], net); d1 = full(z[1:nf] .+ δ[1:nf], net)
            for s in (:H, :Li, :e)
                haskey(d0, s) || continue
                dv = d1[s] - d0[s]; d0[s] + λ * dv < 0 && (λ = min(λ, 0.9 * d0[s] / -dv))
            end
            abs(δ[nf+1]) * λ > 0.5 && (λ = 0.5 / abs(δ[nf+1]))
            z .+= λ .* δ
            λ == 1.0 && maximum(abs.(δ[1:nf]) ./ max.(abs.(z[1:nf]), 1e-20)) < tol && abs(δ[nf+1]) < tol && return (z, true)   # 1e-20: abundance floor
        end
        (z, false)
    end
    function adv(x0, T0, h, depth)
        z, ok = attempt(x0, T0, h)
        ok && return (z[1:nf], exp(z[nf+1]))
        depth >= 12 && error("ChemNet.step: Newton failed (nH = $nH, T0 = $T0, dt = $dt)")
        x1, T1 = adv(x0, T0, h / 2, depth + 1)
        adv(x1, T1, h / 2, depth + 1)
    end
    adv(x0, T0, dt, 0)
end

"""
Steady state of the network at fixed (ρ, T) (for tests): pseudo-transient continuation — backward-Euler steps at fixed
T with dt growing from 1e-8 s by ×10 to 1e30 s (robust from far-from-equilibrium starts).
"""
function steady(x0, ρ, T, net; env, photos = Photo[])
    nH = ρ / μH(); nf = length(net.free)
    ph = [(p, photorates(p, T, env)...) for p in photos]
    x = collect(Float64, x0); dt = 1e-8
    while dt <= 1e30
        xo = copy(x)
        for it in 1:100
            F(xx) = xx .- xo .- dt .* dxdt(xx, nH, T, net, env, ph)
            R = F(x); Jm = ForwardDiff.jacobian(F, x)
            δ = -(Jm \ R); all(isfinite, δ) || break
            λ = 1.0
            for i in 1:nf; x[i] + λ * δ[i] < 0 && (λ = min(λ, 0.9 * x[i] / -δ[i])); end
            d0 = full(x, net); d1 = full(x .+ δ, net)
            for s in (:H, :Li, :e); haskey(d0, s) || continue; dv = d1[s] - d0[s]; d0[s] + λ * dv < 0 && (λ = min(λ, 0.9 * d0[s] / -dv)); end
            x .+= λ .* δ
            λ == 1.0 && maximum(abs.(δ) ./ max.(abs.(x), 1e-30)) < 1e-13 && break
        end
        dt *= 10
    end
    x
end

end # module
