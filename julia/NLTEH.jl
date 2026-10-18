"""
    NLTEH

Non-LTE hydrogen in the infall of a TITAN snapshot, on top of Spectrum.jl (requires Z0Physics and Spectrum loaded).

Model atom: H I levels n = 1…12 (l-collapsed, g = 2n²) and H II. Temperatures and densities are fixed (TITAN's, with
the Unsöld–Lucy correction below the shock as in spectra_run.jl); only the hydrogen level populations and the electron
density are solved for, in the points flagged non-LTE (the infall). Helium stays in LTE.

Rates
* Bound–bound radiative: Sobolev escape probabilities in the infall's velocity field (moving gas, v ≫ v_thermal),
  with the local continuum: net upward rate = β[(n_l B_lu − n_u B_ul) J_c − n_u A_ul], β the angle-averaged escape
  probability with τ(μ) = (πe²/m_e c) f λ (n_l − g_l n_u/g_u) / |μ² dv/dr + (1−μ²) v/r|.
  Oscillator strengths: Kramers × the bound–bound gaunt factor of Johnson (1972).
* Bound–free radiative: photoionisation and (spontaneous + induced) recombination with the continuum J_ν from the full
  spherical transfer, which itself uses the non-LTE opacities and emissivities. Cross sections: Verner et al. (1996) for
  n = 1, Kramers × the Johnson (1972) bound–free gaunt factor for n ≥ 2. Only photons below 100 eV count here.
* X-rays: the energy absorbed above 100 eV (by H and He bound–free) produces secondary ionisations of H at 36 eV per
  ion pair (the high-yield limit for a mostly neutral gas, an upper bound on the yield).
* Collisions with electrons: excitation by the van Regemorter (1962) formula, ionisation by Seaton's approximation
  (both as given in Mihalas 1978), multiplied by `cscale` for sensitivity tests; inverse rates by detailed balance.

`solve_point` is the statistical equilibrium of one point; `NOpac` gives the transfer the non-LTE opacities.
"""
module NLTEH

import ..Z0Physics, ..Spectrum
using LinearAlgebra
const Z0 = Z0Physics; const SP = Spectrum
const h = Z0.h; const kB = Z0.kB; const c = Z0.c; const eV = Z0.eV
const NL = 12
const χn = [Z0.χH / n^2 for n in 1:NL]                  # ionisation energy of level n [erg]
const gn = [2.0 * n^2 for n in 1:NL]
const En = [Z0.χH * (1 - 1 / n^2) for n in 1:NL]        # excitation energy [erg]
const πe2mc = 0.026540                                   # πe²/(m_e c) [cm² s⁻¹]

"Johnson (1972) gaunt-factor coefficients (g0, g1, g2) for lower level n."
function gcoef(n)
    n == 1 && return (1.1330, -0.4059, 0.07014)
    n == 2 && return (1.0785, -0.2319, 0.02947)
    (0.9935 + 0.2328 / n - 0.1296 / n^2, -(0.6282 - 0.5598 / n + 0.5299 / n^2) / n, (0.3887 - 1.181 / n + 1.470 / n^2) / n^2)
end
"Absorption oscillator strength n → m (Kramers × Johnson's g(x), x = 1 − (n/m)²); Lyα 0.4165, Hα 0.6404."
function fosc(n, m)
    g0, g1, g2 = gcoef(n); x = 1 - (n / m)^2
    32 / (3π * sqrt(3)) / (n^5 * m^3) * (1 / n^2 - 1 / m^2)^(-3) * (g0 + g1 / x + g2 / x^2)
end
"Photoionisation cross section of level n at ν [cm²]; 0 below threshold."
function σbf(n, ν)
    ν < χn[n] / h && return 0.0
    n == 1 && return SP.σH1(max(h * ν / eV, 13.6))
    g0, g1, g2 = gcoef(n); x = ν * h / χn[n]
    2.815e29 / (n^5 * ν^3) * (g0 + g1 / x + g2 / x^2)
end
"n_i*/(n_p n_e): Saha–Boltzmann factor of level n [cm³]."
Φ(n, T) = 0.5gn[n] * exp(χn[n] / (kB * T)) / Z0.λe3(T)

"e^y E₁(y) (Abramowitz & Stegun 5.1.53 for y ≤ 1, 5.1.56 above)."
function expE1(y)
    if y <= 1
        E1 = -log(y) - 0.57721566 + 0.99999193y - 0.24991055y^2 + 0.05519968y^3 - 0.00976004y^4 + 0.00107857y^5
        return exp(y) * E1
    end
    (y^2 + 2.334733y + 0.250621) / (y^2 + 3.330657y + 1.681534) / y
end
"van Regemorter excitation rate coefficient l → u [cm³ s⁻¹] (neutral atom, P(y) = 0.276 e^y E₁(y))."
function qexc(l, u, T)
    ΔE = En[u] - En[l]; y = ΔE / (kB * T)
    1.2521e-4 / sqrt(T) * (Z0.χH / ΔE) * fosc(l, u) * 0.276 * expE1(y) * exp(-y)
end
"Seaton collisional ionisation rate coefficient of level n [cm³ s⁻¹] (ḡ = 0.1, 0.2, 0.3 for n = 1, 2, ≥ 3)."
function qion(n, T)
    u0 = χn[n] / (kB * T); gbar = n == 1 ? 0.1 : n == 2 ? 0.2 : 0.3
    1.55e13 / sqrt(T) * gbar * σbf(n, χn[n] / h * 1.000001) * exp(-u0) / u0
end

"""
Angle-averaged Sobolev escape probability for line opacity factor K = (πe²/mc) f λ (n_l − g_l n_u/g_u) [s⁻¹] at a
point falling at speed v [cm/s] at radius r with d ln v / d ln r = s.
"""
function βsob(K, v, r, s; nμ = 64)
    K <= 0 && return 1.0
    acc = 0.0
    for i in 1:nμ
        μ = (i - 0.5) / nμ                       # symmetric in μ: average over [0, 1]
        Q = abs(v / r * ((1 - μ^2) + μ^2 * s))
        τ = K / max(Q, 1e-30)
        acc += τ < 1e-6 ? 1.0 : -expm1(-τ) / τ
    end
    acc / nμ
end

"""
Statistical equilibrium of one point. T, nH (total hydrogen), neHe (electrons from helium), v, r, s (velocity field),
ν and wν (frequency grid and trapezoid weights), J (continuum mean intensity on the grid), Γx (X-ray secondary
ionisations per ground-state atom per s), start populations n0 (length NL) and np0. Returns (n, np, ne, β, P) with β
the escape probabilities of the lines (matrix, l < u) and P the net power escaping in the lines [erg cm⁻³ s⁻¹].
"""
function solve_point(T, nH, neHe, v, r, s, ν, wν, J, Γx, n0, np0; cscale = 1.0, caseB = false, ots2 = false, nit = 200, tol = 1e-10)
    n = copy(n0); np = np0; ne = np + neHe
    x = h .* ν ./ (kB * T)
    # bound–free integrals (only below 100 eV; the X-rays enter through Γx)
    Rph = zeros(NL); Rrc = zeros(NL)
    for i in 1:NL, q in eachindex(ν)
        (h * ν[q] > 100eV) && continue
        σ = σbf(i, ν[q]); σ == 0 && continue
        f = 4π * wν[q] * σ / (h * ν[q])
        Rph[i] += f * J[q]
        Rrc[i] += f * (2h * ν[q]^3 / c^2 + J[q]) * exp(-x[q])
    end
    if caseB      # Lyman continuum in detailed balance (no net radiative rate)
        Rph[1] = 0.0; Rrc[1] = 0.0
    end
    if ots2       # on-the-spot Balmer continuum (simple model only, analytic_infall.jl): radiative recombinations to n = 2
        Rrc[2] = 0.0      # are re-absorbed nearby and do not count; J must then hold only the (attenuated) external light
    end
    Jline = zeros(NL, NL)
    lν = log.(ν)
    for l in 1:NL-1, u in l+1:NL
        νl = (En[u] - En[l]) / h
        i = clamp(searchsortedlast(ν, νl), 1, length(ν) - 1)
        t = (log(νl) - lν[i]) / (lν[i+1] - lν[i])
        # interpolate J/B_ν(T), which is smooth across the Wien tail, then multiply by B at the line
        ra, rb = J[i] / max(SP.Bν(ν[i], T), 1e-300), J[i+1] / max(SP.Bν(ν[i+1], T), 1e-300)
        Jline[l, u] = (ra * (1 - t) + rb * t) * SP.Bν(νl, T)
    end
    β = ones(NL, NL)
    for it in 1:nit
        M = zeros(NL + 1, NL + 1)      # M[i, j]: rate (per particle in j) from j into i; column sums vanish
        add!(from, to, rate) = (M[to, from] += rate; M[from, from] -= rate)
        for l in 1:NL-1, u in l+1:NL
            νl = (En[u] - En[l]) / h; λ = c / νl; f = fosc(l, u)
            K = πe2mc * f * λ * (n[l] - gn[l] / gn[u] * n[u])
            β[l, u] = βsob(K, v, r, s)
            Blu = 4π * πe2mc * f / (h * νl); Bul = gn[l] / gn[u] * Blu
            A = 2h * νl^3 / c^2 * Bul
            Cup = cscale * ne * qexc(l, u, T); Cdn = Cup * gn[l] / gn[u] * exp((En[u] - En[l]) / (kB * T))
            add!(l, u, β[l, u] * Blu * Jline[l, u] + Cup)
            add!(u, l, β[l, u] * (A + Bul * Jline[l, u]) + Cdn)
        end
        for i in 1:NL
            Ci = cscale * ne * qion(i, T)
            up = Rph[i] + Ci + (i == 1 ? Γx : 0.0)
            dn = ne * Φ(i, T) * (Rrc[i] + Ci)        # per proton: n_i* = n_p n_e Φ_i
            add!(i, NL + 1, up)
            add!(NL + 1, i, dn)
        end
        b = zeros(NL + 1); M[1, :] .= 1.0; b[1] = nH
        sol = M \ b
        sol = max.(sol, 1e-300 * nH)
        dn_ = maximum(abs.(sol[1:NL] ./ n .- 1)); dp = abs(sol[end] / np - 1)
        # n_e enters the recombination rates: n_p ∝ 1/n_e would make a plain update oscillate, so the geometric mean
        n .= sol[1:NL]; np = sol[end]; ne_new = sqrt(ne * (np + neHe)); dp = max(dp, abs(ne_new / ne - 1)); ne = ne_new
        max(dn_, dp) < tol && break
    end
    # net power escaping in the lines (erg cm⁻³ s⁻¹): Σ hν β [n_u A − (n_l B_lu − n_u B_ul) J_c]
    P = 0.0
    for l in 1:NL-1, u in l+1:NL
        νl = (En[u] - En[l]) / h; f = fosc(l, u)
        Blu = 4π * πe2mc * f / (h * νl); Bul = gn[l] / gn[u] * Blu; A = 2h * νl^3 / c^2 * Bul
        P += h * νl * β[l, u] * (n[u] * A - (n[l] * Blu - n[u] * Bul) * Jline[l, u])
    end
    n, np, np + neHe, β, P
end

"Non-LTE opacity for Spectrum.monochromatic: points with `nl[k] !== nothing` use the non-LTE hydrogen state."
struct NOpac
    base::SP.Opac
    nl::Vector{Any}          # nothing, or (n = level populations, np, ne)
end

function SP.coeffs3(o::NOpac, k, ν)
    st = o.nl[k]
    st === nothing && return SP.coeffs3(o.base, k, ν)
    T = o.base.T[k]; s = o.base.s[k]; pHe = o.base.pHe[k]
    kT = kB * T; xν = h * ν / kT; stim = -expm1(-xν); E = h * ν / eV
    n, np, ne = st.n, st.np, st.ne
    λμ = c / ν * 1e4; λÅ = c / ν * 1e8
    # hydrogen bound–free (non-LTE): κ = Σ σ (n_i − n_i* e^{−x}), η = (2hν³/c²) Σ σ n_i* e^{−x}
    κH = 0.0; ηH = 0.0
    for i in 1:NL
        σ = σbf(i, ν); σ == 0 && continue
        nst = np * ne * Φ(i, T)
        κH += σ * (n[i] - nst * exp(-xν)); ηH += σ * nst * exp(-xν)
    end
    ηH *= 2h * ν^3 / c^2
    # LTE (thermal) parts: helium bound–free, free–free, H⁻ (in LTE relative to n₁ and n_e)
    κth = 0.0
    for i in 1:Z0.NLEV_HE2
        ν >= Z0.χHe2 / (h * i^2) && (κth += pHe[i] * (i == 1 ? SP.σHe2(E) : 2.815e29 * 16 / (i^5 * ν^3)))
    end
    ν >= Z0.χHe1 / h && (κth += s.nHe0 * SP.σHe1(E))
    κth *= stim
    κth += 3.69e8 * ne * (np + s.nHe1 + 4s.nHe2) / sqrt(T) / ν^3 * stim
    nHm = n[1] * ne * exp(Z0.χHm / kT) / (4Z0.λe3(T))
    κth += nHm * Z0.σHm_bf(λμ) * stim + n[1] * Z0.kHm_ff(λμ, T) * ne * kT
    κ = max(κH + κth, 1e-300)
    sc = ne * Z0.σT + n[1] * Z0.σRay(λÅ)
    κ, sc, (ηH + κth * SP.Bν(ν, T)) / κ
end

end # module
