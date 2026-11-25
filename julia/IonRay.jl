"""
    IonRay

M3: the star's light as a radial, multi-group, photon-conserving ray (Abel, Norman & Madau 1999; zone by zone outward, as C²-Ray, Mellema et al. 2006),
coupled to ChemNet. It supplies what one gray radiation field cannot: the spectrum above and below the Lyman edge.

For group g with photon rate N_g [s⁻¹] entering a shell (r1, r2), absorbers s with abundance x_s = n_s/n_H give Δτ_g = Σ_s x_s n_H σ_{s,g} (r2 − r1), and
    absorbed: N_g (1 − e^{−Δτ_g}),   leaving: N_g e^{−Δτ_g}.
The rate per H nucleus for species s is Σ_g N_g (1 − e^{−Δτ_g}) (x_s n_H σ_{s,g} Δr/Δτ_g) / (n_H V), written as a function of the new abundances so that
ChemNet's Newton solves it implicitly: a zone can never use more photons than enter it, however long the step, and what it leaves goes to the next zone.
The absorbed photon energy goes into u; ChemNet's u holds the binding energies, so the ionization energy is paid from it and the rest is heat.
"""
module IonRay

import ..ChemNet, ..ChemRates
const CN = ChemNet; const CR = ChemRates
const hP = CN.hP; const c = CN.c; const kB = CN.kB; const eV = CN.eV
const σSB = 5.670374419e-5

σ_H(ν) = CR.σ_verner(ν, 13.598, CR.VERNER_HI)
"H₂ + γ → H₂⁺ + e: AAZN97 eq. (σ₂₄ in their numbering), the fit with threshold 15.42 eV."
function σ_H2(ν)
    E = hP * ν / eV
    E < 15.42 ? 0.0 : E < 16.5 ? 6.2e-18 * E - 9.4e-17 : E < 17.7 ? 1.4e-18 * E - 1.48e-17 : 2.5e-14 * E^-2.71
end
"Groups: `Ec` centre [erg], `dE` width, and for each absorber (species, product-side rule in `apply!`) its cross section per group."
struct Groups
    Ec::Vector{Float64}; dE::Vector{Float64}
    species::Vector{Symbol}          # absorbers
    σ::Matrix{Float64}               # group × absorber [cm²]
end

"Log-spaced groups in [E1, E2] eV (n of them) plus the Lyman band up to Emax [eV]. σ evaluated at the group centre of 8 sub-points, weighted by 1/E (flat in photon number)."
function groups(; Elo = 13.598, Emax = 300.0, n = 32, species = [:H], σs = Dict(:H => σ_H, :H2 => σ_H2))
    edges = exp.(range(log(Elo), log(Emax), length = n + 1)) .* eV
    Ec = sqrt.(edges[1:end-1] .* edges[2:end]); dE = diff(edges)
    σ = zeros(n, length(species))
    for (g, s) in Iterators.product(1:n, eachindex(species))
        ps = range(edges[g], edges[g+1], length = 9); pc = 0.5 .* (ps[1:end-1] .+ ps[2:end])
        σ[g, s] = sum(σs[species[s]](e / hP) for e in pc) / 8
    end
    Groups(Ec, dE, species, σ)
end

"Photons per second in each group for a blackbody of luminosity L and temperature Teff (L = 4πR²σT⁴): N_g = L π∫B_ν/(hν)dν / (σT⁴)."
function blackbody(G::Groups, L, Teff)
    N = zeros(length(G.Ec))
    for g in eachindex(N)
        E1 = G.Ec[g] - 0.0; dE = G.dE[g]
        lo = E1 - dE / 2; ps = range(lo, lo + dE, length = 9); s = 0.0
        for i in 1:8
            E = 0.5 * (ps[i] + ps[i+1]); ν = E / hP
            B = 2hP * ν^3 / c^2 / expm1(E / (kB * Teff))
            s += π * B / E * (ps[i+1] - ps[i]) / hP
        end
        N[g] = L * s / (σSB * Teff^4)
    end
    N
end

"""
    source(ρ, r1, r2, Nin, G, net) -> (; extra, heat, accepted, out)

Closures for one zone over one step: `extra(x)` the dx/dt added to ChemNet's free species [s⁻¹ per H], `heat(d)` the absorbed energy [erg cm⁻³ s⁻¹] from the full composition `d`
(for ChemNet's `cooling`, with a minus sign), `accepted` for ChemNet.step, `out(x, dt)` the photons per group leaving. Products: H → H⁺ + e, H₂ → H₂⁺ + e.
"""
function source(ρ, r1, r2, Nin, G::Groups, net)
    nH = ρ / CN.μH(); dr = r2 - r1; V = 4π / 3 * (r2^3 - r1^3)
    iH = findfirst(==(:Hp), net.free); iH2 = findfirst(==(:H2), net.free); iH2p = findfirst(==(:H2p), net.free)
    sp = G.species
    function dτ(g, xs)
        s = zero(eltype(xs)); for j in eachindex(sp); s += xs[j] * nH * G.σ[g, j] * dr; end; s
    end
    # abundances of the absorbers in the order of G.species, from the free vector
    xabs(x) = (d = CN.full(x, net); [max(d[s], zero(d[s])) for s in sp])
    function extra(x)
        xs = xabs(x); T = eltype(xs)
        f = zeros(T, length(net.free))
        for g in eachindex(Nin)
            Nin[g] > 0 || continue
            τ = dτ(g, xs); τ > 0 || continue
            ab = Nin[g] * (-expm1(-τ))
            for j in eachindex(sp)
                dx = ab * (xs[j] * nH * G.σ[g, j] * dr / τ) / (nH * V)
                if sp[j] === :H
                    f[iH] += dx                          # H → H⁺ (H is the eliminated species)
                elseif sp[j] === :H2
                    f[iH2] -= dx; f[iH2p] += dx          # H₂ → H₂⁺
                end
            end
        end
        f
    end
    function heat(d)
        xs = [max(d[s], zero(d[s])) for s in sp]
        q = zero(eltype(xs))
        for g in eachindex(Nin)
            Nin[g] > 0 || continue
            q += Nin[g] * (-expm1(-dτ(g, xs))) * G.Ec[g]
        end
        q / V
    end
    # photons used, summed over the sub-steps ChemNet.step actually took: ending on a transparent zone must not hand the whole step's photons to the next one
    used = zeros(length(Nin)); tused = Ref(0.0)
    function accepted(x, h)
        xs = xabs(x); tused[] += h
        for g in eachindex(Nin); used[g] += Nin[g] * (-expm1(-dτ(g, xs))) * h; end
    end
    # photons per group leaving over the whole step, from the count above (dt: the step). No sub-steps taken yet (a zone not stepped): transmission at x.
    function out(x, dt)
        tused[] > 0 && return [max(Nin[g] - used[g] / dt, 0.0) for g in eachindex(Nin)]
        xs = xabs(x); [Nin[g] * exp(-dτ(g, xs)) for g in eachindex(Nin)]
    end
    (; extra, heat, accepted, out)
end

# =====================================================================================================================
# Option A Extension: Smooth C^∞ Multigroup Photon-Conserving Ray (H⁻ + LW with H₂ self-shielding + H I / H₂ / He I / He II)
# =====================================================================================================================

"Smooth C^∞ (1 − e^{−τ})/τ with Taylor series below |τ| = 1e-4: exact to 1e-16 in value & all Dual derivatives."
@inline function phi1(τ)
    if abs(τ) < 1e-4
        return 1 - τ * (0.5 - τ * (1/6 - τ * (1/24 - τ / 120)))
    else
        return -expm1(-τ) / τ
    end
end

"Smooth C^∞ non-negative projection ½(x + √(x² + ε²)): avoids derivative corners at x = 0."
@inline spos(x, ε = 1e-30) = 0.5 * (x + sqrt(x * x + ε * ε))

"C² quintic smootherstep on [0, 1] with zero 1st and 2nd derivatives at both endpoints."
@inline function smoother5(u)
    u <= 0 && return zero(u)
    u >= 1 && return one(u)
    u^3 * (10 - 15u + 6u^2)
end

σ_HeI(ν)  = CR.σ_verner(ν, 24.587, CR.VERNER_HeI)
σ_HeII(ν) = CR.σ_verner(ν, 54.418, CR.VERNER_HeII)
σ_Hm(ν)   = CR.σ_Hm(ν)

"""
    f_shield_H2(NH2, T) -> f_sh ∈ (0, 1]

Smooth C^∞ H₂ Lyman–Werner self-shielding factor (Wolcott-Green, Haiman & Bryan 2011; Draine & Bertoldi 1996)
using α = 1.1 for warm primordial protostellar envelopes. Differentiable in both NH2 and T (no clamp corners).
"""
@inline function f_shield_H2(NH2, T)
    x = spos(NH2, 1e-30) / 5.0e14
    Tm = spos(T - 10.0, 1.0) + 10.0
    b5 = sqrt(2 * kB * Tm / (2.01565006 * 1.66053907e-24)) / 1.0e5
    sx = sqrt(1 + x)
    0.965 / (1 + x / b5)^1.1 + 0.035 / sx * exp(-8.5e-4 * sx)
end

"Extended photon-conserving ray groups: H⁻ (0.754–11.2 eV), LW (11.2–13.598 eV), and EUV (13.598–300 eV, H/H₂/He I/He II)."
struct FullGroups
    Ec::Vector{Float64}             # group centres [erg]
    dE::Vector{Float64}             # group widths [erg]
    Elo::Vector{Float64}            # lower edges [eV]
    Ehi::Vector{Float64}            # upper edges [eV]
    species::Vector{Symbol}         # ionizing/detaching absorbers [:Hm, :H, :H2, :He, :Hep]
    χ::Vector{Float64}              # ionization/detachment thresholds [erg]
    σ::Matrix{Float64}              # (ngroup × nspecies) cross sections [cm²]
    g_lw::Int                       # index of the Lyman–Werner (11.2–13.598 eV) group
    σ_lw_eff::Float64               # effective Solomon dissociation cross section [cm²] (k_diss = σ_lw_eff * F_LW * f_sh)
end

"""
    groups_full(; n_opt = 6, n_euv = 28, Emax = 300.0) -> FullGroups

Multi-band group set with exact edge alignment at H⁻ (0.754 eV), LW (11.2 eV), H I (13.598 eV),
H₂ (15.42 eV), He I (24.587 eV), and He II (54.418 eV).
"""
function groups_full(; n_opt = 6, n_euv = 28, Emax = 300.0)
    # Exact threshold-aligned sub-grids (eV)
    e_opt = collect(range(0.7542, 11.2; length = n_opt + 1))
    e_lw  = [11.2, 13.5984]
    n1 = max(3, round(Int, 0.15 * n_euv))
    n2 = max(4, round(Int, 0.25 * n_euv))
    n3 = max(5, round(Int, 0.30 * n_euv))
    n4 = max(5, n_euv - n1 - n2 - n3)
    e_h   = exp.(range(log(13.5984), log(15.42);   length = n1 + 1))
    e_h2  = exp.(range(log(15.42),   log(24.5874); length = n2 + 1))
    e_he1 = exp.(range(log(24.5874), log(54.4178); length = n3 + 1))
    e_he2 = exp.(range(log(54.4178), log(Emax);    length = n4 + 1))
    edges_eV = vcat(e_opt[1:end-1], e_lw[1:end-1], e_h[1:end-1], e_h2[1:end-1], e_he1[1:end-1], e_he2)
    ng = length(edges_eV) - 1
    Elo = edges_eV[1:end-1]; Ehi = edges_eV[2:end]
    g_lw = n_opt + 1
    Ec = sqrt.(Elo .* Ehi) .* eV
    dE = (Ehi .- Elo) .* eV
    species = [:Hm, :H, :H2, :He, :Hep]
    χ = [0.7542, 13.5984, 15.42, 24.5874, 54.4178] .* eV
    σfns = (σ_Hm, σ_H, σ_H2, σ_HeI, σ_HeII)
    σ = zeros(ng, length(species))
    for g in 1:ng, (j, fn) in enumerate(σfns)
        e1 = Elo[g] * eV; e2 = Ehi[g] * eV
        ps = range(e1, e2; length = 9); pc = 0.5 .* (ps[1:end-1] .+ ps[2:end])
        σ[g, j] = sum(fn(e / hP) for e in pc) / 8
    end
    # Effective LW Solomon dissociation cross section: Abel+97 k_diss = 1.38e9 J_LW [s⁻¹]
    # For isotropic/radial flux F_LW = N_LW / (4π r²), J_LW = (F_LW / 4π) * (hν_LW / Δν_LW) => σ_lw_eff ≈ 3.7e-18 cm²
    Δν_lw = (13.5984 - 11.2) * eV / hP
    ν_lw  = 12.4 * eV / hP
    σ_lw_eff = 1.38e9 * (hP * ν_lw / (4π * Δν_lw)) / 1.0 # ≈ 3.73e-18 cm²
    FullGroups(Ec, dE, Elo, Ehi, species, χ, σ, g_lw, σ_lw_eff)
end

"Blackbody photon rate N_g [s⁻¹] per group for FullGroups."
function blackbody_full(G::FullGroups, L, Teff)
    ng = length(G.Ec); N = zeros(ng)
    (L <= 0 || Teff <= 50) && return N
    R2 = L / (4π * σSB * Teff^4)
    for g in 1:ng
        e1 = G.Elo[g] * eV; e2 = G.Ehi[g] * eV
        ps = range(e1, e2; length = 9); pc = 0.5 .* (ps[1:end-1] .+ ps[2:end])
        de = (e2 - e1) / 8
        s = 0.0
        for e in pc
            ν = e / hP; x = e / (kB * Teff)
            x > 500 && continue
            B = 2hP * ν^3 / c^2 / expm1(x)
            s += (π * B / e) * (de / hP)
        end
        N[g] = 4π * R2 * s
    end
    N
end

"""
    source_full(ρ, r1, r2, Nin, NH2_in, T_gas, G::FullGroups, net; binding_in_u = true)
        -> (; extra, heat, accepted, out, NH2_out)

Smooth C^∞ photon-conserving radial ray closures for H⁻ photodetachment, LW H₂ Solomon dissociation
(with Wolcott-Green/Draine-Bertoldi self-shielding), and H I / H₂ / He I / He II EUV photoionization.
All terms use `phi1(τ)` and `spos(x)` so ForwardDiff Dual Jacobians are strictly smooth with zero hard edges.
"""
function source_full(ρ, r1, r2, Nin, NH2_in, T_gas, G::FullGroups, net; binding_in_u = true)
    nH = ρ / CN.μH(); dr = r2 - r1; rmid = 0.5 * (r1 + r2)
    V = 4π / 3 * (r2^3 - r1^3); area = 4π * rmid^2
    iHp   = findfirst(==(:Hp),   net.free)
    iHm   = findfirst(==(:Hm),   net.free)
    iH2   = findfirst(==(:H2),   net.free)
    iH2p  = findfirst(==(:H2p),  net.free)
    iHep  = findfirst(==(:Hep),  net.free)
    iHepp = findfirst(==(:Hepp), net.free)
    sp = G.species; ns = length(sp); ng = length(Nin)

    function xabs(x)
        d = CN.full(x, net)
        T = eltype(x)
        xs = Vector{T}(undef, ns)
        @inbounds for j in 1:ns
            s = sp[j]
            xs[j] = haskey(d, s) ? spos(d[s]) : zero(T)
        end
        xs
    end

    @inline function dτ_g(g, xs)
        s = zero(eltype(xs))
        @inbounds for j in 1:ns
            s += xs[j] * nH * G.σ[g, j] * dr
        end
        s
    end

    function extra(x)
        xs = xabs(x); Tx = eltype(xs)
        f = zeros(Tx, length(net.free))
        # 1. Photon-conserving continuum absorption (H⁻, H I, H₂, He I, He II) using smooth phi1(τ)
        @inbounds for g in 1:ng
            Nin[g] > 0 || continue
            τ = dτ_g(g, xs)
            ph = phi1(τ)                 # (1 - exp(-τ)) / τ, strictly C^∞ at τ -> 0!
            fac = Nin[g] * ph * dr / V   # [cm³ s⁻¹] * nH cancels in dx = Nin * phi1(τ) * xs[j]*σ*dr / V
            for j in 1:ns
                σgj = G.σ[g, j]; σgj > 0 || continue
                dx = fac * xs[j] * σgj
                sj = sp[j]
                if sj === :H && iHp !== nothing
                    f[iHp] += dx
                elseif sj === :Hm && iHm !== nothing
                    f[iHm] -= dx
                elseif sj === :H2 && iH2 !== nothing && iH2p !== nothing
                    f[iH2] -= dx; f[iH2p] += dx
                elseif sj === :He && iHep !== nothing
                    f[iHep] += dx
                elseif sj === :Hep && iHep !== nothing && iHepp !== nothing
                    f[iHep] -= dx; f[iHepp] += dx
                end
            end
        end
        # 2. Smooth self-shielded Lyman–Werner H₂ Solomon dissociation (H₂ + γ_LW -> 2H)
        if iH2 !== nothing && Nin[G.g_lw] > 0
            xH2 = xs[3]
            dNH2 = xH2 * nH * dr
            fsh = f_shield_H2(NH2_in + 0.5 * dNH2, T_gas)
            # Photon-conserving cap via phi1 on effective LW optical depth τ_lw = xH2 * nH * σ_lw * fsh * dr
            τlw = xH2 * nH * G.σ_lw_eff * fsh * dr
            # Photodissociation rate per H nucleus: k_diss * xH2 = (Nin_lw / area) * σ_lw_eff * fsh * phi1(τlw) * xH2
            # Notice (Nin_lw / V) * (1 - exp(-τlw)) / nH = (Nin_lw / (nH * V)) * τlw * phi1(τlw)
            dx_diss = (Nin[G.g_lw] / area) * G.σ_lw_eff * fsh * phi1(τlw) * xH2
            f[iH2] -= dx_diss
        end
        f
    end

    function heat(d)
        Tx = promote_type(typeof(T_gas), eltype(values(d)))
        xs = Vector{Tx}(undef, ns)
        @inbounds for j in 1:ns
            s = sp[j]
            xs[j] = haskey(d, s) ? spos(Tx(d[s])) : zero(Tx)
        end
        q = zero(Tx)
        @inbounds for g in 1:ng
            Nin[g] > 0 || continue
            τ = dτ_g(g, xs)
            ph = phi1(τ)
            if binding_in_u
                q += Nin[g] * τ * ph * G.Ec[g]
            else
                for j in 1:ns
                    σgj = G.σ[g, j]; σgj > 0 || continue
                    ex = max(G.Ec[g] - G.χ[j], 0.0)
                    q += Nin[g] * ph * (xs[j] * nH * σgj * dr) * ex
                end
            end
        end
        # LW Solomon dissociation pumping heat (~0.4 eV per dissociated H₂)
        if Nin[G.g_lw] > 0
            xH2 = xs[3]
            fsh = f_shield_H2(NH2_in + 0.5 * xH2 * nH * dr, T_gas)
            τlw = xH2 * nH * G.σ_lw_eff * fsh * dr
            q += Nin[G.g_lw] * (4π * rmid^2 / area) * τlw * phi1(τlw) * (0.4 * eV)
        end
        q / V
    end

    used = zeros(ng); tused = Ref(0.0); dNH2_acc = Ref(0.0)
    function accepted(x, h)
        xs = xabs(x); tused[] += h
        @inbounds for g in 1:ng
            τ = dτ_g(g, xs)
            used[g] += Nin[g] * τ * phi1(τ) * h
        end
        dNH2_acc[] += xs[3] * nH * dr * h
    end

    function out(x, dt)
        xs = xabs(x)
        if tused[] > 0
            Nout = [max(Nin[g] - used[g] / dt, 0.0) for g in 1:ng]
            NH2_o = NH2_in + dNH2_acc[] / dt
        else
            Nout = [Nin[g] * exp(-dτ_g(g, xs)) for g in 1:ng]
            NH2_o = NH2_in + xs[3] * nH * dr
        end
        # Attenuate LW band smoothly by ratio of self-shielding factors across the shell (never exceeds 1)
        if Nin[G.g_lw] > 0
            ratio = f_shield_H2(NH2_o, T_gas) / max(f_shield_H2(NH2_in, T_gas), 1e-30)
            Nout[G.g_lw] = min(Nout[G.g_lw], Nin[G.g_lw] * clamp(ratio, 0.0, 1.0))
        end
        Nout, NH2_o
    end

    (; extra, heat, accepted, out)
end

end # module
