"""
    StellarOpacityLow

S3, low-temperature side of `StellarOpacity` (used by the table builder and the tests, not at run time):

- `equilibrium(ρ, T, X)`: chemical/ionisation equilibrium of H/He gas with the primordial lithium electrons, for any X
  (Y = 1 − X): H₂ ⇌ 2H (PrimordialChem.K_H2), H ⇌ H⁺ + e⁻, He ⇌ He⁺ + e⁻ (He⁺⁺ neglected: used to log T ≈ 4.2), Li ⇌ Li⁺ + e⁻
  (⁷Li/H = 4.15e-10, as popiii_regimes.jl and PrimordialChem.ne_with_Li), H⁻ by Saha. Ideal gas, no pressure effects.
  (PopIIIEqEOS.composition is a table for X ≈ 0.75 with SCVH95 at high pressure; this solver is the same Saha chemistry
  for general X.)
- `means_at(ρ, T, X)`: Rosseland (extinction) and Planck (absorption) means of the continuum of ContinuumOpacity (H₂–H₂,
  H₂–He, H₂–H CIA; H⁻ bound–free and free–free; hydrogenic H; Rayleigh on H₂ and H; Thomson). Differences from
  `ContinuumOpacity.means`: He Rayleigh added (σ_He = 5.484e-14 λ⁻⁴ [1 + (2.44e5 + 5.94e10/(λ²−2.9e5))/λ²]², λ in Å
  [from memory of ATLAS SIGHE, unverified; matters only at X ≲ 0.05]) and both Rayleigh formulae capped at λ = 1300 Å
  (as Z0Physics). No lines of H₂ or atoms (H₂ lines are `H2Lines`' job).
  Many (ρ, X) at one T share the frequency-dependent coefficients (`prep`), so a table column is cheap.
- `z0_means(ρ, T, X)`: the Z0Physics continuum (Saha H/He, H bf n ≤ 12, He⁺ bf, ff with g = 1, H⁻, Rayleigh H, Thomson) for
  general X: Rosseland and Planck mean (no molecules, no lines) — the κ_P above 10⁴ K.
"""
module StellarOpacityLow

import ..PrimordialChem as PC
import ..ContinuumOpacity as CO
import ..Z0Physics as Z0

const kB = CO.kB; const h = CO.h; const c = CO.c; const me = CO.me; const eV = CO.eV
const mH = PC.mH; const mHe = Z0.mHe
const χHe = 24.587389 * eV

"Equilibrium composition: number densities [cm⁻³] (nHn = H nuclei, nH atoms, nH2 molecules, nHe neutral He, nHp, ne, nHm)."
function equilibrium(ρ, T, X; li = true)
    NH = max(X, 1e-12) * ρ / mH; NHe = (1 - X) * ρ / (4 * mH)
    L = (2π * me * kB * T / h^2)^1.5
    KH2 = PC.K_H2(T); Kion = PC.K_ion(T); KHe = 4L * exp(-χHe / (kB * T)); KLi = PC.K_Li(T)
    cc = KH2 > 0 ? 2NH / KH2 : Inf
    function parts(ne)
        b = 1 + Kion / ne
        x = 2 / (b + sqrt(b^2 + 4cc)); fp = (b - 1) * x
        (x, fp, NHe * KHe / (KHe + ne), li ? PC.yLi * NH * KLi / (KLi + ne) : 0.0)
    end
    g(lne) = (ne = exp(lne); (x, fp, nhe, nli) = parts(ne); fp * NH + nhe + nli - ne)
    lo = -460.0; hi = log(NH + 2NHe + 1e-300) + 1
    if g(lo) <= 0
        lne = lo
    else
        for _ in 1:200
            mid = 0.5(lo + hi); g(mid) > 0 ? (lo = mid) : (hi = mid)
            hi - lo < 1e-13 && break
        end
        lne = 0.5(lo + hi)
    end
    ne = exp(lne); (x, fp, nhe, nli) = parts(ne)
    nH = x * NH; nHp = fp * NH; nH2 = max(NH - nH - nHp, 0.0) / 2
    nHm = nH * ne * (h^2 / (2π * me * kB * T))^1.5 / 4 * exp(0.754195 * eV / (kB * T))
    (; nHn = NH, nH, nH2, nHe = NHe - nhe, nHep = nhe, nHp, ne, nHm)
end

σ_ray_He(λÅ) = (λ = max(λÅ, 1300.0); 5.484e-14 / λ^4 * (1 + (2.44e5 + 5.94e10 / (λ^2 - 2.9e5)) / λ^2)^2)

"Frequency grid [Hz]: log-uniform in wavenumber, edges."
function grid(T; n = 4000)
    σlo = min(10.0, 0.03 * kB * T / (h * c)); σhi = max(1e5, 25kB * T / (h * c))
    [c * 10^x for x in range(log10(σlo), log10(σhi), length = n + 1)]
end

"Composition-independent coefficients at temperature T on the grid."
function prep(T; n = 4000)
    t = CO.tables(); νe = grid(T; n); νm = 0.5 .* (νe[1:end-1] .+ νe[2:end]); dν = diff(νe)
    σ = νm ./ c; λμ = 1e4 ./ σ
    stim = -expm1.(-h .* νm ./ (kB * T))
    (; νm, dν, T,
       k22 = [CO.k_h2h2(t, s, T) for s in σ], k2He = [CO.kcia(t.H2He, s, T) for s in σ], k2H = [CO.kcia(t.H2H, s, T) for s in σ],
       sHm = [CO.σ_Hm_bf(l) for l in λμ] .* stim, kHmff = [CO.k_Hm_ff(l, T) for l in λμ] .* (kB * T),
       cbf = [CO.α_H_kramers(ν, T, 1.0, 0.0, 0.0) for ν in νm], cff = [CO.α_H_kramers(ν, T, 0.0, 1.0, 1.0) for ν in νm],
       rH2 = [CO.σ_ray_H2(max(l * 1e3, 130.0)) for l in λμ], rH = [CO.σ_ray_H(max(l * 1e4, 1300.0)) for l in λμ],
       rHe = [σ_ray_He(l * 1e4) for l in λμ],
       B = [CO.Bν(ν, T) for ν in νm], D = [CO.dBdT(ν, T) for ν in νm])
end

"Means (κR, κP, κRabs) [cm² g⁻¹] from the prepared coefficients and a composition."
function means(p, ρ, s)
    sP = 0.0; sB = 0.0; sR = 0.0; sD = 0.0; sRa = 0.0
    @inbounds for i in eachindex(p.νm)
        a = p.k22[i] * s.nH2^2 + p.k2He[i] * s.nH2 * s.nHe + p.k2H[i] * s.nH2 * s.nH + p.sHm[i] * s.nHm + p.kHmff[i] * s.nH * s.ne +
            p.cbf[i] * s.nH + p.cff[i] * s.nHp * s.ne
        sc = p.rH2[i] * s.nH2 + p.rH[i] * s.nH + p.rHe[i] * s.nHe + CO.σT * s.ne
        w = p.dν[i]
        sP += a * p.B[i] * w; sB += p.B[i] * w
        sR += p.D[i] / max(a + sc, 1e-300) * w; sRa += p.D[i] / max(a, 1e-300) * w; sD += p.D[i] * w
    end
    (; κR = sD / sR / ρ, κP = sP / sB / ρ, κRabs = sD / sRa / ρ)
end
means_at(ρ, T, X; li = true, n = 4000) = means(prep(T; n), ρ, equilibrium(ρ, T, X; li))

"Z0Physics continuum Rosseland and Planck means for general X (Saha, no molecules)."
function z0_means(ρ, T, X)
    s = Z0.saha(ρ, T; comp = Z0.Comp(X, 1 - X)); pH, pHe = Z0.level_pops(T, s)
    numR = 0.0; denR = 0.0; numP = 0.0; denP = 0.0
    ug = Z0.ugrid_with_edges(T)
    for i in 1:length(ug)-1
        u = 0.5(ug[i] + ug[i+1]); du = ug[i+1] - ug[i]; ν = u * kB * T / h
        wR = u^4 * exp(u) / expm1(u)^2; wP = u^3 / expm1(u)
        ab, sc = Z0.kappa_nu(ν, T, s, pH, pHe)
        numR += wR / (ab + sc) * du; denR += wR * du
        numP += wP * ab * du; denP += wP * du
    end
    (; κR = denR / numR / ρ, κP = numP / denP / ρ)
end

end # module
