# Z = 0 Rosseland and Planck means in LTE from frequency-dependent sources (atomic data and profiles: atomic.jl):
#   electron scattering (Thomson), H Rayleigh scattering,
#   H I and He II bound-free (levels n ≤ NLEV_*, hydrogenic with bound-free Gaunt factors, dissolved-level pseudo-continua),
#   He I bound-free (ground state Verner et al. 1996, excited states hydrogenic in n*),
#   free-free on H+, He+ and He++ with the van Hoof et al. 2014 Gaunt factors,
#   H− bound-free and free-free (John 1988 fits),
#   bound–bound lines of H I, He II (exact hydrogenic f, Stark + Doppler + radiative profiles) and He I (NIST multiplets).
# NOT included: H2/H2+/CIA molecules (below ~3000–3500 K), non-LTE.
# Populations: Saha gives the ground-state density (nH0, nHe1 of `saha` carry g₁ only), excited levels follow from Boltzmann
# with occupation probabilities; `opacity_state` makes the ionization balance consistent with the resulting partition functions.

const NLEV_H   = 16
const NLEV_HE2 = 16

# --- H− cross sections (John 1988, A&A 193, 189) -----------------------------------------
"H− bound-free cross section (cm²), λ in μm."
function σHm_bf(λμ)
    λ0 = 1.6419
    (λμ >= λ0 || λμ < 0.125) && return 0.0
    C = (152.519, 49.534, -118.858, 92.536, -34.194, 4.982)
    x = 1 / λμ - 1 / λ0
    f = sum(C[n] * x^((n - 1) / 2) for n in 1:6)
    return 1e-18 * λμ^3 * x^1.5 * f
end

const _HmA1 = ((0.0,2483.346,-3449.889,2200.040,-696.271,88.283),
               (0.0,285.827,-1158.382,2427.719,-1841.400,444.517),
               (0.0,-2054.291,8746.523,-13651.105,8624.970,-1863.864),
               (0.0,2827.776,-11485.632,16755.524,-10051.530,2095.288),
               (0.0,-1341.537,5303.609,-7510.494,4400.067,-901.788),
               (0.0,208.952,-812.939,1132.738,-655.020,132.985))
const _HmA2 = ((518.1021,473.2636,-482.2089,115.5291,0.0,0.0),
               (-734.8666,1443.4137,-737.1616,169.6374,0.0,0.0),
               (1021.1775,-1977.3395,1096.8827,-245.6490,0.0,0.0),
               (-479.0721,922.3575,-521.1341,114.2430,0.0,0.0),
               (93.1373,-178.9275,101.7963,-21.9972,0.0,0.0),
               (-6.4285,12.3600,-7.0920,1.5097,0.0,0.0))
"H− free-free absorption per neutral H atom per unit electron pressure (cm⁴/dyn), incl. stimulated emission."
function kHm_ff(λμ, T)
    λμ < 0.1823 && return 0.0
    A = λμ > 0.3645 ? _HmA1 : _HmA2
    θ = 5040.0 / T
    s = 0.0
    for n in 1:6
        s += θ^((n + 1) / 2) * (λμ^2 * A[1][n] + A[2][n] + A[3][n] / λμ + A[4][n] / λμ^2 +
                                A[5][n] / λμ^3 + A[6][n] / λμ^4)
    end
    return 1e-29 * s
end

"Rayleigh scattering on H I (Dalgarno), λ in Å; capped below 1300 Å (near Lyα)."
function σRay(λÅ)
    λ = max(λÅ, 1300.0)
    return 5.799e-13 / λ^4 + 1.422e-6 / λ^6 + 2.784 / λ^8
end

# --- monochromatic opacity (cm⁻¹) --------------------------------------------------------
"""
Level populations (cm⁻³) of H I and He II for a Saha state `s` (`saha` or `opacity_state`): nH0/nHe1 of `saha` are the ground-state
densities, n_n = n_1 w_n n² exp(−E_n/kT) / w_1.  With `opacity_state` the arrays it computed are returned.
"""
function level_pops(T, s)
    hasproperty(s, :pH) && return s.pH, s.pHe
    kT = kB * T
    Nion = s.nHp + s.nHe1 + s.nHe2; Nneu = s.nH0 + s.nHe0
    w1H = occupation(1, 1, Nion, Nneu); w1He = occupation(1, 2, Nion, Nneu)
    pH  = [s.nH0  * n^2 * exp(-χH  * (1 - 1 / n^2) / kT) * occupation(n, 1, Nion, Nneu) / w1H  for n in 1:NLEV_H]
    pHe = [s.nHe1 * n^2 * exp(-χHe2 * (1 - 1 / n^2) / kT) * occupation(n, 2, Nion, Nneu) / w1He for n in 1:NLEV_HE2]
    return pH, pHe
end

const VERNER_HEI = (E0 = 1.361e1, σ0 = 9.492e2, ya = 1.469, P = 3.188, yw = 2.039, y0 = 4.434e-1, y1 = 2.136)   # verner/phfit2.f, ph2(·, 2, 2)
"He I ground-state photoionization cross section [cm²], Verner et al. 1996."
function σHeI_ground(ν)
    E = h * ν / eV; E < χHe1 / eV && return 0.0
    v = VERNER_HEI; x = E / v.E0 - v.y0; y = sqrt(x^2 + v.y1^2)
    1e-18 * v.σ0 * ((x - 1)^2 + v.yw^2) * y^(0.5 * v.P - 5.5) * (1 + sqrt(y / v.ya))^(-v.P)
end

"""
Ionization balance with the occupation-probability partition functions of H I, He I and He II (used for the opacities only; the EOS
keeps `saha`). The returned state extends `saha`'s: nH0 and nHe1 are total neutral H and He⁺ (all levels), pH, pHe the level populations,
wH, wHe occupation probabilities, heI the excited He I groups (ν_th, n*, population).
"""
function opacity_state(ρ, T; comp = COMP)
    s0 = saha(ρ, T; comp); kT = kB * T
    nH = s0.nH; nHe = s0.nHe; SH, S1, S2 = _saha_rates(T)
    E_H = [χH * (1 - 1 / n^2) for n in 1:NLEV_H]; E_He = [χHe2 * (1 - 1 / n^2) for n in 1:NLEV_HE2]
    lev = hei_data().levels
    nsHeI = [sqrt(Ry∞ / max(χHe1 - E * eV, 1e-3 * eV)) for (E, _) in lev]
    bH = [n^2 * exp(-E_H[n] / kT) for n in 1:NLEV_H]; bHe = [n^2 * exp(-E_He[n] / kT) for n in 1:NLEV_HE2]
    bI = [g * exp(-E * eV / kT) for (E, g) in lev]
    ne = s0.ne; Nion = s0.nHp + s0.nHe1 + s0.nHe2; Nneu = s0.nH0 + s0.nHe0
    wH = ones(NLEV_H); wHe = ones(NLEV_HE2); wI = ones(length(lev)); UH = UHe = UI = 1.0
    xH = s0.xH; f0 = 0.0; y1 = 0.0; y2 = 0.0
    for pass in 1:12
        for n in 1:NLEV_H;   wH[n]  = occupation(n, 1, Nion, Nneu); end
        for n in 1:NLEV_HE2; wHe[n] = occupation(n, 2, Nion, Nneu); end
        for i in eachindex(lev); wI[i] = occupation(nsHeI[i], 1, Nion, Nneu); end
        UH = max(sum(wH .* bH), 1e-60); UHe = max(sum(wHe .* bHe), 1e-60); UI = max(sum(wI .* bI), 1e-60)   # floors: fully pressure-ionized
        # n₊ n_e / N₀ = S / U (ground-state-normalised S)
        sH = SH / UH; s1 = S1 * UHe / UI; s2 = S2 / UHe
        lo = log(max(1e-30, 1e-12 * (nH + nHe))); hi = log(nH + 2nHe)
        for _ in 1:70
            m = 0.5(lo + hi); e = exp(m)
            _fractions(e, nH, nHe, sH, s1, s2)[5] > e ? (lo = m) : (hi = m)
        end
        nen = exp(0.5(lo + hi)); ne = pass == 1 ? nen : sqrt(ne * nen)
        xH, f0, y1, y2, _ = _fractions(ne, nH, nHe, sH, s1, s2)
        Nion = nH * xH + nHe * (y1 + y2); Nneu = nH * (1 - xH) + nHe * f0
    end
    nHI = nH * (1 - xH); nHe0 = nHe * f0; nHeII = nHe * y1
    pH = [nHI * wH[n] * bH[n] / UH for n in 1:NLEV_H]; pHe = [nHeII * wHe[n] * bHe[n] / UHe for n in 1:NLEV_HE2]
    # He I excited: lump by integer n*
    grp = Dict{Int,Vector{Float64}}()
    for (i, (E, _)) in enumerate(lev)
        (E < 15.0 || nsHeI[i] > 12) && continue
        p = nHe0 * wI[i] * bI[i] / UI; k = round(Int, nsHeI[i])
        v = get!(grp, k, [0.0, 0.0]); v[1] += p; v[2] += p * (χHe1 - E * eV)
    end
    heI = [(ν = g[2] / max(g[1], 1e-300) / h, ns = Float64(k), p = g[1]) for (k, g) in grp if g[1] > 0]
    (; ne, nH0 = nHI, nHp = nH * xH, nHe0, nHe1 = nHeII, nHe2 = nHe * y2, nH, nHe, xH, pH, pHe, wH, wHe, heI,
       pHe0 = lev, wI, bI, UI, n1He = nHe0)
end

"""
Absorption (incl. stimulated emission) and scattering coefficients (cm⁻¹) at frequency ν, continuum only.
`pH`, `pHe` are level populations from `level_pops`.
"""
function kappa_nu(ν, T, s, pH, pHe)
    kT  = kB * T
    u   = h * ν / kT
    stim = -expm1(-u)
    λcm = c / ν; λμ = λcm * 1e4; λÅ = λcm * 1e8
    ab = 0.0; ν3 = ν^3
    dis = hasproperty(s, :wH)
    νH = χH / h
    for n in 1:NLEV_H
        pH[n] <= 0 && continue
        νn = νH / n^2
        if ν >= νn
            ab += pH[n] * KRAMERS / (n^5 * ν3) * gaunt_bf(n, ν / νn)
        elseif dis && n < NLEV_H       # dissolved levels n* > n below the edge (Hubeny, Hummer & Lanz 1994)
            q = 1 / n^2 - ν / νH; q <= 1 / 1600 && begin
                d = 1 - occupation(1 / sqrt(max(q, 1e-12)), 1, s.nHp + s.nHe1 + s.nHe2, s.nH0 + s.nHe0) / max(s.wH[n], 1e-300)
                d > 0 && (ab += d * pH[n] * KRAMERS / (n^5 * ν3) * gaunt_bf(n, 1.0))
            end
        end
    end
    νHe = χHe2 / h
    for n in 1:NLEV_HE2
        pHe[n] <= 0 && continue
        νn = νHe / n^2
        if ν >= νn
            ab += pHe[n] * KRAMERS * 16 / (n^5 * ν3) * gaunt_bf(n, ν / νn)
        elseif dis && n < NLEV_HE2
            q = 1 / n^2 - ν / νHe; q <= 1 / 1600 && begin
                d = 1 - occupation(1 / sqrt(max(q, 1e-12)), 2, s.nHp + s.nHe1 + s.nHe2, s.nH0 + s.nHe0) / max(s.wHe[n], 1e-300)
                d > 0 && (ab += d * pHe[n] * KRAMERS * 16 / (n^5 * ν3) * gaunt_bf(n, 1.0))
            end
        end
    end
    if dis
        n1 = s.nHe0 / s.UI * s.wI[1] * s.bI[1]
        ab += n1 * σHeI_ground(ν)
        for g in s.heI
            ν >= g.ν && (ab += g.p * KRAMERS / (g.ns^5 * ν3))
        end
    else
        ab += s.nHe0 * σHeI_ground(ν)
    end
    ab *= stim
    ab += 3.69e8 * s.ne / sqrt(T) / ν3 * stim * ((s.nHp + s.nHe1) * gaunt_ff(1, T, ν) + 4s.nHe2 * gaunt_ff(2, T, ν))
    n1H = pH[1]
    nHm = n1H * s.ne * exp(χHm / kT) / (4λe3(T))
    ab += nHm * σHm_bf(λμ) * stim
    ab += n1H * kHm_ff(λμ, T) * s.ne * kT
    sc = s.ne * σT + n1H * σRay(λÅ)
    return ab, sc
end

"Σ over lines of (ν₀, S): the integrated strengths (profile-free), for the Planck mean."
function line_strengths(s, T; nmaxU = 30, nmaxL = 8)
    kT = kB * T; out = Tuple{Float64,Float64}[]
    Nion = s.nHp + s.nHe1 + s.nHe2; Nneu = s.nH0 + s.nHe0
    for (Z, pops, ws, χ, nl) in ((1, s.pH, s.wH, χH, NLEV_H), (2, s.pHe, s.wHe, χHe2, NLEV_HE2))
        for n in 1:min(nmaxL, nl)
            pops[n] <= 0 && continue
            wn = max(ws[n], 1e-300)
            for np in n+1:nmaxU
                r = occupation(np, Z, Nion, Nneu) / wn; r < 1e-3 && break
                push!(out, (χ / h * (1 / n^2 - 1 / np^2), πe2mc * hyd_f(n, np) * pops[n] * min(r, 1.0)))
            end
        end
    end
    UI = s.UI
    for l in hei_data().lines
        pl = s.nHe0 * l.gi * exp(-l.Ei * eV / kT) / UI; pl <= 0 && continue
        ns = sqrt(Ry∞ / max(χHe1 - l.Ek * eV, 1e-3 * eV)); wu = occupation(ns, 1, Nion, Nneu); wu < 1e-3 && continue
        push!(out, ((l.Ek - l.Ei) * eV / h, πe2mc * l.f * pl * wu))
    end
    out
end

"""
Profiles of the lines that matter at (ρ, T) (peak opacity above `tol` × the local continuum estimate κ_floor [cm⁻¹]); hydrogenic lines
with exact Stark patterns, He I multiplets Doppler + radiative + electron impact (n* from the levels).
"""
function line_list(s, T; κ_floor = 0.0, tol = 1e-3, nmaxU = 30, nmaxL = 8)
    kT = kB * T; ne = s.ne
    Fi = 2.6031 * ecgs * (s.nHp + s.nHe1 + 2.828 * s.nHe2)^(2 / 3); Fe = 2.6031 * ecgs * ne^(2 / 3)
    Ls = LineProf[]
    for (Z, pops, ws, m, χ, nl) in ((1, s.pH, s.wH, mH, χH, NLEV_H), (2, s.pHe, s.wHe, mHe, χHe2, NLEV_HE2))
        Ci = ecgs * a0cm * Fi / (Z * h); Ce = ecgs * a0cm * Fe / (Z * h)
        Nion = s.nHp + s.nHe1 + s.nHe2; Nneu = s.nH0 + s.nHe0
        for n in 1:min(nmaxL, nl)
            pops[n] <= 0 && continue
            wn = max(ws[n], 1e-300)
            for np in n+1:nmaxU
                wu = occupation(np, Z, Nion, Nneu); r = wu / wn; r < 1e-3 && break
                ν0 = χ / h * (1 / n^2 - 1 / np^2); D = ν0 * sqrt(2kT / (m * c^2))
                f = hyd_f(n, np); S = πe2mc * f * pops[n] * min(r, 1.0)
                γe = impact_width(n, np, Z, ne, T); pat = stark_pattern(n, np)
                W = max(D, γe, maximum(pat.X) * Ci)
                S / (π * W) * (-expm1(-h * ν0 / kT)) < tol * κ_floor && continue
                γr = 0.6670 * ν0^2 / c^2 * f * (n^2 / np^2) / (4π) * 2
                push!(Ls, line_profile(ν0, S, D, γr, γe, pat.X, pat.I, Ci, Ce))
            end
        end
    end
    hd = hei_data(); UI = s.UI
    for l in hd.lines
        ν0 = (l.Ek - l.Ei) * eV / h; D = ν0 * sqrt(2kT / (mHe * c^2))
        pl = s.nHe0 * l.gi * exp(-l.Ei * eV / kT) / UI       # w_l ≈ 1 for the low levels that carry the lines
        pl <= 0 && continue
        ns = sqrt(Ry∞ / max(χHe1 - l.Ek * eV, 1e-3 * eV)); nl = sqrt(Ry∞ / max(χHe1 - l.Ei * eV, 1e-3 * eV))
        wu = occupation(ns, 1, s.nHp + s.nHe1 + s.nHe2, s.nH0 + s.nHe0); wu < 1e-3 && continue
        S = πe2mc * l.f * pl * wu; γe = impact_width(nl, ns, 1, ne, T)
        S / (π * max(D, γe)) * (-expm1(-h * ν0 / kT)) < tol * κ_floor && continue
        γr = 0.6670 * ν0^2 / c^2 * l.f * (l.gi / l.gk) / (4π) * 2
        push!(Ls, line_profile(ν0, S, D, γr, γe, Float64[], Float64[], 0.0, 0.0))
    end
    Ls
end

const UGRID = exp.(range(log(1e-3), log(60.0), length = 4000))

"""
Frequency grid (u = hν/kT) with every opacity discontinuity inserted as a breakpoint, so that
midpoint quadrature never straddles an edge: H I edges n = 1..NLEV_H, He II edges n = 1..NLEV_HE2,
He I ground edge, H− bound-free threshold (1.6419 μm) and cutoff (0.125 μm), the H− free-free fit
switches (0.3645, 0.1823 μm) and the Rayleigh cap (1300 Å). Without this, edges sliding through
the fixed grid as T changes gave ~1 % quadrature noise in κ_P (κ_R was barely affected).
"""
function ugrid_with_edges(T)
    kT = kB * T
    edges = Float64[]
    for n in 1:NLEV_H;   push!(edges, χH / (n^2 * kT)); end
    for n in 1:NLEV_HE2; push!(edges, χHe2 / (n^2 * kT)); end
    push!(edges, χHe1 / kT)
    for λμ in (1.6419, 0.125, 0.3645, 0.1823, 0.13)
        push!(edges, h * c / (λμ * 1e-4) / kT)
    end
    filter!(u -> UGRID[1] < u < UGRID[end], edges)
    return sort!(unique!(vcat(UGRID, edges)))
end

"""
Rosseland and Planck means (cm²/g) at (ρ, T). Planck: absorption only, continuum + Σ lines. Rosseland: continuum + scattering, with the lines
resolved on an added per-line frequency grid when `lines_R` (slower). `lines = false` drops lines from both.
"""
function kappa_means(ρ, T; comp = COMP, lines = true, lines_R = true)
    s = opacity_state(ρ, T; comp)
    pH, pHe = s.pH, s.pHe
    kT = kB * T
    ug = ugrid_with_edges(T); N = length(ug) - 1
    abc = zeros(N); scc = zeros(N)
    numR = 0.0; denR = 0.0; numP = 0.0; denP = 0.0
    for i in 1:N
        u = 0.5(ug[i] + ug[i+1]); du = ug[i+1] - ug[i]
        ν = u * kT / h
        wR = u^4 * exp(u) / expm1(u)^2; wP = u^3 / expm1(u)
        ab, sc = kappa_nu(ν, T, s, pH, pHe); abc[i] = ab; scc[i] = sc
        numR += wR / (ab + sc) * du; denR += wR * du
        numP += wP * ab * du;        denP += wP * du
    end
    κc = denR / numR       # cm⁻¹
    lines || return denR / numR / ρ, numP / denP / ρ
    for (ν0, S) in line_strengths(s, T)
        u0 = h * ν0 / kT; u0 > 300 && continue
        numP += u0^3 / expm1(u0) * (-expm1(-u0)) * S * (h / kT)
    end
    lines_R || return denR / numR / ρ, numP / denP / ρ
    Ls = line_list(s, T; κ_floor = κc, tol = 1e-3)
    isempty(Ls) && return denR / numR / ρ, numP / denP / ρ
    # Rosseland with lines: common frequency nodes around each line, continuum piecewise constant on the cells of ug
    g = copy(ug)
    for l in Ls
        uc = h * l.ν0 / kT; (uc < ug[1] || uc > ug[end]) && continue
        w0 = h * exp(l.lx0) / kT * 1e3
        push!(g, uc)
        for x in exp.(range(log(0.02w0), log(max(30.0, 2w0)), length = 90))
            uc + x < ug[end] && push!(g, uc + x); uc - x > ug[1] && push!(g, uc - x)
        end
    end
    sort!(g); unique!(g); M = length(g) - 1
    κ = zeros(M); um = 0.5 .* (g[1:end-1] .+ g[2:end]); ix = 1
    for i in 1:M
        while ix < N && ug[ix+1] <= um[i]; ix += 1; end
        κ[i] = abc[ix] + scc[ix]
    end
    for l in Ls
        uc = h * l.ν0 / kT; i0 = clamp(searchsortedfirst(um, uc), 1, M)
        for dir in (1, -1)
            i = dir == 1 ? i0 : i0 - 1
            while 1 <= i <= M
                u = um[i]; ν = u * kT / h
                v = l.S * (-expm1(-u)) * prof_val(l, ν)
                v < 1e-4 * κ[i] && break
                κ[i] += v; i += dir
            end
        end
    end
    numR = 0.0; denR = 0.0
    for i in 1:M
        u = um[i]; du = g[i+1] - g[i]; wR = u^4 * exp(u) / expm1(u)^2
        numR += wR / κ[i] * du; denR += wR * du
    end
    return denR / numR / ρ, numP / denP / ρ
end
kappa_ross(ρ, T) = kappa_means(ρ, T)[1]
kappa_planck(ρ, T) = kappa_means(ρ, T)[2]

# --- table on (log T, log R), R = ρ / T6³ -------------------------------------------------
struct KapTable
    lT::Vector{Float64}
    lR::Vector{Float64}
    lk::Matrix{Float64}        # log10 κ_Rosseland
    lkP::Matrix{Float64}       # log10 κ_Planck (absorption)
end
KapTable(lT, lR, lk) = KapTable(lT, lR, lk, fill(NaN, size(lk)))   # v1 tables (no Planck)

"""
Table on (log T, log R): κ_P and the continuum κ_R on every node; the line correction Δ = log₁₀(κ_R,full/κ_R,cont), slow to compute, on a coarse
grid (Δlog T 0.05, Δlog R 0.25; zero outside log T 3.4–5.8), Catmull–Rom interpolated.
"""
function build_table(; lT = 3.3:0.01:8.0, lR = -12.0:0.05:1.0)
    lTv = collect(lT); lRv = collect(lR)
    cT = collect(3.4:0.05:5.8); cR = collect(floor(lRv[1]):0.25:ceil(lRv[end]))
    Δ = zeros(length(cT), length(cR))
    Threads.@threads :dynamic for k in 1:length(cT)*length(cR)
        i, j = Tuple(CartesianIndices((length(cT), length(cR)))[k])
        T = 10^cT[i]; ρ = 10^cR[j] * (T / 1e6)^3
        Δ[i, j] = log10(kappa_means(ρ, T; lines = true, lines_R = true)[1] / kappa_means(ρ, T; lines = false)[1])
    end
    @info "line correction grid done"
    lk = zeros(length(lTv), length(lRv)); lkP = similar(lk)
    Threads.@threads :dynamic for i in eachindex(lTv)
        T = 10^lTv[i]
        for j in eachindex(lRv)
            ρ = 10^lRv[j] * (T / 1e6)^3
            kR, kP = kappa_means(ρ, T; lines = true, lines_R = false)
            x = (lTv[i] - cT[1]) / (cT[2] - cT[1]); y = (lRv[j] - cR[1]) / (cR[2] - cR[1])     # 0-based node coordinates
            d = (0 <= x < length(cT) - 1 && 0 <= y < length(cR) - 1) ? _cr2(Δ, x, y, length(cT), length(cR)) : 0.0
            lk[i, j] = log10(kR) + max(d, 0.0); lkP[i, j] = log10(kP)
        end
    end
    return KapTable(lTv, lRv, lk, lkP)
end

# Catmull–Rom cubic weights: interpolating and C¹ across cell edges.
@inline function _crw(t)
    t2 = t * t; t3 = t2 * t
    return ((-t3 + 2t2 - t) / 2, (3t3 - 5t2 + 2) / 2, (-3t3 + 4t2 + t) / 2, (t3 - t2) / 2)
end

"""
C¹ bicubic (Catmull–Rom) interpolation of a log10 table in (log T, log R); generic in the
number type, so ForwardDiff derivatives are continuous across cell edges. Outside the table
the arguments are clamped (derivative then zero) — `in_table` reports whether that happened.
"""
function _interp(tab::KapTable, M::Matrix{Float64}, ρ, T)
    lT = log10(T); lR = log10(ρ) - 3(lT - 6)
    nT, nR = length(tab.lT), length(tab.lR)
    x = (lT - tab.lT[1]) / (tab.lT[2] - tab.lT[1])
    y = (lR - tab.lR[1]) / (tab.lR[2] - tab.lR[1])
    x = clamp(x, 0, nT - 1 - 1e-9)
    if EXTRAP_LOWR[] && pval(y) < 0
        # below the lowest log R: linear extrapolation in log R with the slope between the first two columns
        # (fully ionized, thin gas: free-free/bound-free ∝ ρ, slope 1 for κ_P; electron scattering, slope 0 for κ_R)
        v0 = _cr2(M, x, zero(y), nT, nR); v1 = _cr2(M, x, one(y), nT, nR)
        return 10^(v0 + y * (v1 - v0))
    end
    y = clamp(y, 0, nR - 1 - 1e-9)
    return 10^_cr2(M, x, y, nT, nR)
end
"Catmull–Rom value of table M at fractional node coordinates (x, y) (both inside the table)."
function _cr2(M, x, y, nT, nR)
    i = floor(Int, pval(x)) + 1; j = floor(Int, pval(y)) + 1
    j = min(j, nR - 1); i = min(i, nT - 1)
    wx = _crw(x - (i - 1)); wy = _crw(y - (j - 1))
    v = zero(x * y)
    for a in 1:4, b in 1:4
        ii = clamp(i - 2 + a, 1, nT); jj = clamp(j - 2 + b, 1, nR)
        v += wx[a] * wy[b] * M[ii, jj]
    end
    v
end
"""
EXTRAP_LOWR[]: below the table's lowest log R (= log ρ − 3 log T₆) extrapolate log κ linearly in log R instead of
clamping. Clamping evaluates the opacity of gas up to 1e3 times denser: in a 1e7–1e8 K, 1e-17 g cm⁻³ accretion
corona it overstates the Planck mean, hence free-free cooling, by that factor. Default true.
"""
const EXTRAP_LOWR = Ref(true)
kap(tab::KapTable, ρ, T) = _interp(tab, tab.lk, ρ, T)
kapP(tab::KapTable, ρ, T) = _interp(tab, tab.lkP, ρ, T)

"Old bilinear lookup (v1), kept for regression comparison."
function kap_bilinear(tab::KapTable, ρ, T)
    lT = log10(T); lR = log10(ρ) - 3(lT - 6)
    x = clamp((lT - tab.lT[1]) / (tab.lT[2] - tab.lT[1]), 0, length(tab.lT) - 1.000001)
    y = clamp((lR - tab.lR[1]) / (tab.lR[2] - tab.lR[1]), 0, length(tab.lR) - 1.000001)
    i = floor(Int, pval(x)) + 1; j = floor(Int, pval(y)) + 1
    fx = x - (i - 1); fy = y - (j - 1)
    v = (1 - fx) * (1 - fy) * tab.lk[i, j] + fx * (1 - fy) * tab.lk[i+1, j] +
        (1 - fx) * fy * tab.lk[i, j+1] + fx * fy * tab.lk[i+1, j+1]
    return 10^v
end

in_table(tab::KapTable, ρ, T) = (lT = log10(pval(T)); lR = log10(pval(ρ)) - 3(lT - 6);
    tab.lT[1] <= lT <= tab.lT[end] && tab.lR[1] <= lR <= tab.lR[end])
