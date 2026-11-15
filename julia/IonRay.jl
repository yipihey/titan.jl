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

end # module
