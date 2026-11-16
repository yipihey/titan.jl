"""
    G2Opacity

G2 ("one non-LTE material everywhere", docs/POPIII_G2_MATERIAL.md): the continuum means of a network zone
(PopIIIMaterial.nlte_means) with the continuum ingredients that the equilibrium-table side (Z0Physics, used for the table's κ_P
above 10^3.8 K; OPAL for its κ_R) has and the network's code lacked. Same arguments and same conventions as nlte_means
(number densities from the zone's own composition; `hbfP` = H bound–free in κ_E/κ_P or not; κ_E weighted with B_ν(Tw)), plus the
helium ions. Switch `G2OPAC=1` (default off; off = no call, runs bit-identical): PopIIIMaterial.nlte_means / nlte_means_fast return
this instead.

Added against nlte_means (sources: docs/Z0_EOS_OPACITY.md and julia/z0/; each is Z0Physics' own function, evaluated on the network's
populations instead of Saha's):
- frequency grid to 25 kT/hc at 3e5 K (nlte_means' grid stops at 1e5 cm⁻¹ = 12.4 eV: above ~3e4 K the Rosseland and Planck weights
  lie beyond it and its means were the low-frequency free–free alone, too high by up to 10³);
- level bookkeeping: `n1` is the zone's neutral H (all levels, ChemNet's x_H), n ≥ 2 follow from `n2`, the ground state is the rest,
  n1 − Σ_{n≥2} p_n (nlte_means took all neutral H as ground state and added the excited levels on top: above ~3e4 K, where the
  partition function is ≫ 2, that overstated the bound–free opacity by the partition function, ×250 at 1e5 K); H⁻ free–free and
  Rayleigh on the ground state; `n2_lte` gives the LTE n = 2 for an equilibrium composition;
- H bound–free with Gaunt factors (n = 1 Stobbe exact, n ≥ 2 Seaton fit: Z0Physics.gaunt_bf) instead of Kramers g = 1, levels
  n = 2 … 16 from the zone's n = 2 by Boltzmann × Hummer & Mihalas (1988) occupation probabilities (Z0Physics.occupation), and the
  dissolved-level pseudo-continua below each edge (Hubeny, Hummer & Lanz 1994), as Z0Physics.kappa_nu;
- He I ground bound–free (Verner et al. 1996, Z0Physics.σHeI_ground) on the zone's neutral He; He II bound–free (hydrogenic Z = 2,
  levels in Boltzmann with the zone's He⁺, occupation probabilities, Gaunt factors);
- free–free of H⁺ and He⁺ (Z = 1) and He²⁺ (Z = 2) with the van Hoof et al. (2014) Gaunt factors (Z0Physics.gaunt_ff) instead of
  H⁺ alone with g = 1;
- Rayleigh cross sections of H and H₂ held at their 1300 Å values at shorter wavelengths (Z0Physics.σRay; the fits diverge towards Lyα).
Not added (see the doc): bound–bound lines of H/He (in OPAL, and Z0Physics' line lists are LTE and ~0.05–4 s per point), He⁻ and H₂⁻
free–free, H₂ lines, pressure effects on H⁻.
"""
module G2Opacity

import ..Z0Physics, ..ContinuumOpacity
const Z0 = Z0Physics; const COp = ContinuumOpacity

const ON = get(ENV, "G2OPAC", "0") == "1"
const NLEV = 16
const TTOP = 3e5            # the grid reaches 25 kT/hc at this T (the HOT branch's T1)

const FREQ = Ref{Any}(nothing)
"Frequency-only quantities (computed once): midpoints, widths, λ, ν³, σ(H⁻ bf), capped Rayleigh σ of H₂ and H, σ(He I ground)."
function freq()
    FREQ[] === nothing || return FREQ[]
    νs = COp.grid(; n = 1500, T = TTOP); n = length(νs) - 1
    ν = [0.5 * (νs[i] + νs[i+1]) for i in 1:n]; dν = [νs[i+1] - νs[i] for i in 1:n]; σ = ν ./ COp.c; λμ = 1e4 ./ σ
    λc = max.(λμ, 0.13)                                           # 1300 Å
    FREQ[] = (; ν, dν, σ, λμ, ν3 = ν .^ 3, σhm = COp.σ_Hm_bf.(λμ), ray2 = COp.σ_ray_H2.(λc .* 1e3), ray1 = COp.σ_ray_H.(λc .* 1e4),
                σHeI = Z0.σHeI_ground.(ν))
end

"""
    means(ρ, T; n1, n2, nH2, nHe, nHep, nHepp, nHp, ne, nHm, nHmS, hbfP, Tw = T) -> (κ_E, κ_P, κ_R) [cm² g⁻¹]

As PopIIIMaterial.nlte_means (n1: H atoms, taken as the ground state; n2: population of n = 2 [cm⁻³]; nHe: neutral He) with the
additions listed in the module docstring.
"""
function means(ρ, T; n1, n2, nH2, nHe, nHep = 0.0, nHepp = 0.0, nHp, ne, nHm, nHmS, hbfP, Tw = T)
    ρ = Float64(ρ); T = Float64(T); F = freq(); t = COp.tables(); kT = COp.kB * T
    sE = 0.0; sP = 0.0; sB = 0.0; sR = 0.0; sD = 0.0; sW = 0.0; wJ = Tw != T
    Nion = nHp + nHep + nHepp; Nneu = n1 + nHe + nH2
    χH = Z0.χH; χ2 = Z0.χHe2; KR = Z0.KRAMERS
    # level populations: H n ≥ 2 from n2 (Boltzmann relative to n = 2, occupation probabilities w_n/w_2); He II from He⁺ in LTE
    wH = ntuple(n -> Z0.occupation(n, 1, Nion, Nneu), NLEV); wHe = ntuple(n -> Z0.occupation(n, 2, Nion, Nneu), NLEV)
    pH0 = ntuple(n -> n == 1 ? 0.0 : n2 * (n^2 / 4) * exp(-χH * (0.25 - 1 / n^2) / kT) * wH[n] / max(wH[2], 1e-300), NLEV)
    ng = max(n1 - sum(pH0), 0.0)                                  # ground state
    pH = ntuple(n -> n == 1 ? ng : pH0[n], NLEV)
    UHe = sum(n^2 * exp(-χ2 * (1 - 1 / n^2) / kT) * wHe[n] for n in 1:NLEV)
    pHe = ntuple(n -> nHep * n^2 * exp(-χ2 * (1 - 1 / n^2) / kT) * wHe[n] / max(UHe, 1e-300), NLEV)
    cia = nH2 > 1e-3 * n1; acia = 0.0
    νH = χH / COp.h; νHe = χ2 / COp.h
    @inbounds for i in eachindex(F.ν)
        ν = F.ν[i]; dν = F.dν[i]; σ = F.σ[i]; λμ = F.λμ[i]; ν3 = F.ν3[i]
        stim = -expm1(-COp.h * ν / kT)
        cia && i % 8 == 1 && (acia = COp.k_h2h2(t, σ, T) * nH2^2 + COp.kcia(t.H2He, σ, T) * nH2 * nHe + COp.kcia(t.H2H, σ, T) * nH2 * n1)
        # free–free: H⁻ (John 1988), H⁺ + He⁺ (Z = 1) and He²⁺ (Z = 2) with van Hoof Gaunt factors
        aff = COp.k_Hm_ff(λμ, T) * ng * ne * kT
        aff += 3.692e8 * ne / (sqrt(T) * ν3) * stim * ((nHp + nHep) * Z0.gaunt_ff(1, T, ν) + (nHepp > 0 ? 4nHepp * Z0.gaunt_ff(2, T, ν) : 0.0))
        # H bound–free (Gaunt factors, occupation probabilities, dissolved levels)
        ab = 0.0
        for n in 1:NLEV
            pH[n] <= 0 && continue
            νn = νH / n^2
            if ν >= νn
                ab += pH[n] * KR / (n^5 * ν3) * Z0.gaunt_bf(n, ν / νn)
            elseif n < NLEV
                q = 1 / n^2 - ν / νH
                if q <= 1 / 1600
                    d = 1 - Z0.occupation(1 / sqrt(max(q, 1e-12)), 1, Nion, Nneu) / max(wH[n], 1e-300)
                    d > 0 && (ab += d * pH[n] * KR / (n^5 * ν3) * Z0.gaunt_bf(n, 1.0))
                end
            end
        end
        ab *= stim
        # helium bound–free: He I ground (Verner), He II hydrogenic Z = 2
        ahe = nHe * F.σHeI[i]
        if nHep > 0 && ν >= νHe / NLEV^2
            for n in 1:NLEV
                νn = νHe / n^2
                ν >= νn && (ahe += pHe[n] * KR * 16 / (n^5 * ν3) * Z0.gaunt_bf(n, ν / νn))
            end
        end
        ahe *= stim
        a0 = acia + aff + ahe
        σhm = F.σhm[i] * stim
        sca = F.ray2[i] * nH2 + F.ray1[i] * ng + COp.σT * ne
        B = COp.Bν(ν, T); D = COp.dBdT(ν, T); h = hbfP ? ab : 0.0
        BJ = wJ ? COp.Bν(ν, Tw) : B
        sE += (a0 + σhm * nHm + h) * BJ * dν; sW += BJ * dν; sP += (a0 + σhm * nHmS + h) * B * dν; sB += B * dν
        sR += D / max(a0 + σhm * nHm + ab + sca, 1e-300) * dν; sD += D * dν
    end
    (sE / sW / ρ, sP / sB / ρ, sD / sR / ρ)
end

"""
LTE population of H n = 2 [cm⁻³] for `nH0` neutral H atoms (all levels): Boltzmann with occupation probabilities over n ≤ $(NLEV)
(Nion, Nneu as in `means`).
"""
function n2_lte(T, nH0, Nion, Nneu)
    kT = COp.kB * Float64(T); χH = Z0.χH
    w = ntuple(n -> Z0.occupation(n, 1, Nion, Nneu), NLEV)
    U = sum(n^2 * exp(-χH * (1 - 1 / n^2) / kT) * w[n] / w[1] for n in 1:NLEV)
    nH0 * 4 * exp(-χH * 0.75 / kT) * w[2] / w[1] / U
end

end # module
