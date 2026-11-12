"""
    H2Lines

Population III project, M3: transfer of the H₂ ro-vibrational (quadrupole) lines in spherical symmetry, as RHFC02
App. A1.4/A2.4: levels v = 0–2, J = 0–20 (all quadrupole lines among them), LTE level populations, absorption
coefficient Eq. A11, Gaussian thermal profile, comoving-frame Doppler shifts from the bulk velocity, tangent rays
(impact parameters = interface radii), formal solution of each frequency bin along each ray; the cooling rate of shell
i is Λ_i = [L(r_{i+1}) − L(r_i)]/ΔV_i (Eq. A10).

Line data: ExoMol RACPPK (Roueff et al. 2019; `~/Projects/popiii_lit/H2lines/`): level energies (cm⁻¹), total
degeneracies g = (2J+1) g_ns, A-values (s⁻¹). RHFC02 used Turner, Kirby-Docken & Dalgarno 1977 (not available here):
the modern A-values are the "modern" switch, not a reproduction of RHFC02's numbers. The partition function sums all
302 levels of the list.

Formal solution along a ray, per segment of optical depth Δτ with the source function linear in τ (short
characteristics):  I_out = I_in e^{−Δτ} + S_out (w0 − w1/Δτ) + S_in (w1/Δτ),  w0 = 1 − e^{−Δτ},  w1 = w0 − Δτ e^{−Δτ}
(Δτ ≪ 1: series), exact for S linear in τ, accurate in the diffusion limit.
"""
module H2Lines

using Printf

const h = 6.62607015e-27; const c = 2.99792458e10; const kB = 1.380649e-16; const mH2 = 2.01565006 * 1.66053907e-24

struct LineList
    E::Vector{Float64}        # level energies [erg], all levels (partition function)
    g::Vector{Float64}
    # selected lines (RHFC02 subset)
    iu::Vector{Int}; il::Vector{Int}; A::Vector{Float64}; ν::Vector{Float64}
end

function load(dir = joinpath(get(ENV, "POPIII_LIT", joinpath(homedir(), "Projects", "popiii_lit")), "H2lines"); vmax = 2, Jmax = 20)
    st = [split(l) for l in readlines(joinpath(dir, "1H2__RACPPK.states")) if !isempty(strip(l))]
    n = length(st); E = zeros(n); g = zeros(n); J = zeros(Int, n); v = zeros(Int, n)
    for s in st
        k = parse(Int, s[1]); E[k] = parse(Float64, s[2]) * h * c; g[k] = parse(Float64, s[3])
        J[k] = parse(Int, s[4]); v[k] = parse(Int, s[5])
    end
    iu = Int[]; il = Int[]; A = Float64[]; ν = Float64[]
    for l in readlines(joinpath(dir, "1H2__RACPPK.trans"))
        s = split(l); isempty(s) && continue
        u, lo = parse(Int, s[1]), parse(Int, s[2])
        (v[u] <= vmax && v[lo] <= vmax && J[u] <= Jmax && J[lo] <= Jmax) || continue
        push!(iu, u); push!(il, lo); push!(A, parse(Float64, s[3])); push!(ν, (E[u] - E[lo]) / h)
    end
    LineList(E, g, iu, il, A, ν)
end

"Partition function Z(T) and the LTE fraction of H₂ molecules in each level."
function populations(L::LineList, T)
    w = L.g .* exp.(-L.E ./ (kB * T)); w ./ sum(w)
end

"Optically thin LTE line cooling per H₂ molecule [erg s⁻¹]: Σ_lines f_u A hν."
function thin_per_molecule(L::LineList, T)
    f = populations(L, T)
    sum(f[L.iu[i]] * L.A[i] * h * L.ν[i] for i in eachindex(L.A))
end

Bν(ν, T) = 2h * ν^3 / c^2 / expm1(h * ν / (kB * T))
"Largest number of bins used by a line and number of lines solved by transfer in the last call."
const LAST = Ref{Any}((nbins = 0, nthick = 0))

"""
    cooling(L, r, u, ρ, T, nH2; nbin = 30, vmax_extra = 0.0) -> (Λ, Lr)

Line cooling rate Λ_k [erg cm⁻³ s⁻¹] of each zone and the line luminosity L(r_i) [erg s⁻¹] at each interface.
r: interfaces (n+1, r[1] = 0 allowed), u: interface velocities (n+1), ρ, T, nH2: zones (n). Zone values are uniform
within the zone; the velocity along a ray inside a zone is interpolated linearly in r between the interfaces.
Frequency bins per line (G6a): span ±(3 v_th + |u|) and bin width ≤ ½ v_th,min over the zones that carry a fraction
`cover` of that line's emission or opacity column (nbin ≤ bins ≤ nbinmax); comoving shift per segment. Lines with a
radial line-centre depth < 2 tol use the thin emission (relative error ≲ tol).
"""
function cooling(L::LineList, r, u, ρ, T, nH2; nbin = 30, fnlte = nothing, tol = 0.01, nbinmax = 200, cover = 0.999,
                 nlte = nothing, nlte_mode = :ali, nlte_tol = 1e-3, nlte_maxit = 60, nlte_init = :ep, nlte_lambda = :diag, inc = nothing)
    if nlte !== nothing
        return cooling_nlte(L, r, u, T, nH2, nlte; nbin, tol, nbinmax, cover, mode = nlte_mode, ntol = nlte_tol,
                            maxit = nlte_maxit, init = nlte_init, lam = nlte_lambda, inc)
    end
    n = length(T); @assert length(r) == n + 1
    vth = [sqrt(2kB * T[k] / mH2) for k in 1:n]
    # thin-line shortcut by tolerance: the thin approximation overestimates a uniform sphere's line luminosity by a
    # relative ≈ (3/8)τ_c for τ_c ≪ 1 (escape from the centre): τ_c < 2 tol keeps the error ≲ tol (h2lines_test §3)
    τthin = 2tol
    nbmax_used = Threads.Atomic{Int}(0)
    fn = fnlte === nothing ? ones(n) : fnlte          # non-LTE factor on the source function (upper-level depletion)
    # LTE populations per zone
    pops = [populations(L, T[k]) for k in 1:n]
    nl = length(L.A)
    Lr = zeros(n + 1)
    # rays: impact parameters p_j = r[j], j = 1..n (r[1] = 0: the central ray); intersections with interfaces k ≥ j
    Lline = zeros(nl, n + 1); nthick = Threads.Atomic{Int}(0)
    Threads.@threads for i in 1:nl
        ν0 = L.ν[i]; gu = L.g[L.iu[i]]; gl = L.g[L.il[i]]
        # line-integrated opacity coefficient per zone: κ0_k = (hν/4π)(c²/2hν³) A (g_u/g_l n_l − n_u)   [cm⁻¹ × Hz]
        κ0 = [h * ν0 / (4π) * c^2 / (2h * ν0^3) * L.A[i] * nH2[k] * (gu / gl * pops[k][L.il[i]] - pops[k][L.iu[i]]) for k in 1:n]
        S = [fn[k] * Bν(ν0, T[k]) for k in 1:n]
        # source at the interfaces (linear in r between zone centres; centre and surface: the zone value) — the segments
        # use S linear in τ between the entry and exit interfaces, so the flux is right in optically thick zones with a
        # gradient (constant S per zone gives the flux only through the jumps between zones)
        Si = similar(S, n + 1); Si[1] = S[1]; Si[n+1] = S[n]
        for k in 2:n
            rc0 = 0.5 * (r[k-1] + r[k]); rc1 = 0.5 * (r[k] + r[k+1]); w = (r[k] - rc0) / (rc1 - rc0)
            Si[k] = (1 - w) * S[k-1] + w * S[k]
        end
        # optically thin line (RHFC02: "If the line absorption depth is < 0.1, we use the optically thin
        # approximation"): radial line-centre depth from the centre to the surface below τthin → emission 4π κ0 S per zone
        τc = sum(κ0[k] / (sqrt(π) * vth[k] / c * ν0) * (r[k+1] - r[k]) for k in 1:n)
        if τc < τthin
            # thin: L(r_k) = Σ_{zones below} 4π κ0 S ΔV (κ0 S = n_u A hν/4π · f for κ0 ≈ n_l B…: use the exact thin emissivity)
            acc = 0.0
            Lline[i, 1] = 0.0
            for k in 1:n
                em = fn[k] * nH2[k] * pops[k][L.iu[i]] * L.A[i] * h * ν0         # erg cm⁻³ s⁻¹
                acc += em * 4π / 3 * (r[k+1]^3 - r[k]^3); Lline[i, k+1] = acc
            end
            continue
        end
        # per-line frequency grid (G6a): the zones that carry `cover` (99.9 %) of this line's emission (n_u A hν ΔV) or of its
        # line-centre opacity column (κ0 Δr/Δν_D) set the span ±(3 v_th + |u|) and the bin width (≤ ½ their smallest v_th)
        em = [fn[k] * nH2[k] * pops[k][L.iu[i]] * L.A[i] * (r[k+1]^3 - r[k]^3) for k in 1:n]
        oc = [κ0[k] / (vth[k] / c * ν0) * (r[k+1] - r[k]) for k in 1:n]
        act = falses(n)
        for w in (em, oc)
            tot = sum(abs, w); tot > 0 || continue
            p_ = sortperm(abs.(w); rev = true); acc_ = 0.0
            for k in p_; act[k] = true; acc_ += abs(w[k]); acc_ >= cover * tot && break; end
        end
        any(act) || (act .= true)
        ka = findall(act)
        vspan = 3maximum(vth[ka]) + maximum(max(abs(u[k]), abs(u[k+1])) for k in ka)
        nb = clamp(ceil(Int, 2vspan / (0.5minimum(vth[ka]))) + 1, nbin, nbinmax)
        Threads.atomic_max!(nbmax_used, nb)
        xs = range(-vspan, vspan, length = nb); dv = step(xs)
        Fν = zeros(n + 1)                          # ∫ flux over the line (erg cm⁻² s⁻¹) at each interface
        ΔI = zeros(n + 1, n)                       # I(μ) − I(−μ) at interface k on ray j (frequency bin), j ≤ k
        for ib in 1:nb
            vb = xs[ib]                            # observer-frame offset: ν = ν0 (1 + vb/c)
            fill!(ΔI, 0.0)
            # Ray j (p = r[j]) crosses interfaces k = j..n+1 at z_k = √(r_k² − p²); segment m (between interfaces
            # j+m−1 and j+m) lies in zone j+m−1. Incoming half (−ẑ) from the outer boundary (I = 0) to z = 0, then the
            # outgoing half; at interface k the outgoing and incoming intensities give ΔI(μ_k), μ_k = z_k/r_k.
            for j in 1:n
                p = r[j]
                nk = n + 1 - j + 1
                zk = [sqrt(max(r[k]^2 - p^2, 0.0)) for k in j:n+1]
                function dτ(m, dir)
                    kz = j + m - 1
                    zm = 0.5 * (zk[m] + zk[m+1]); rm = sqrt(p^2 + zm^2)
                    w = clamp((rm - r[kz]) / (r[kz+1] - r[kz]), 0.0, 1.0)
                    ur = (1 - w) * u[kz] + w * u[kz+1]
                    μ = rm > 0 ? dir * zm / rm : 0.0
                    vcom = vb - μ * ur
                    φ = exp(-(vcom / vth[kz])^2) / (sqrt(π) * vth[kz] / c * ν0)
                    κ0[kz] * φ * (zk[m+1] - zk[m])
                end
                # segment m runs between interfaces j+m−1 (inner) and j+m (outer); at the tangent point of ray j (m = 1,
                # inner end) the source is that of interface j
                Iin = zeros(nk); I = 0.0
                for m in nk-1:-1:1
                    I = sh(I, dτ(m, -1.0), Si[j+m], Si[j+m-1]); Iin[m] = I
                end
                I = Iin[1]; Iout1 = I
                for (m, k) in enumerate(j:n+1)
                    if m > 1
                        I = sh(I, dτ(m - 1, 1.0), Si[j+m-2], Si[j+m-1])
                    end
                    ΔI[k, j] = m == 1 ? 0.0 : I - Iin[m]          # at the tangent point μ = 0: ΔI = 0
                end
            end
            # F(r_k) = 2π ∫₀¹ μ ΔI(μ) dμ with ΔI linear in μ between the rays (exact in the thin limit, ΔI ∝ μ)
            for k in 2:n+1
                acc = 0.0
                for j in 1:k-1                    # interval between rays j+1 (smaller μ) and j (larger μ)
                    μa = sqrt(max(1 - (r[j+1] / r[k])^2, 0.0)); μb = sqrt(max(1 - (r[j] / r[k])^2, 0.0))
                    Ia = j + 1 <= n ? ΔI[k, j+1] : 0.0; Ib = ΔI[k, j]      # ray p = r[n+1] grazes the surface: ΔI = 0
                    d = μb - μa; d <= 0 && continue
                    s1 = (Ib - Ia) / d
                    # ∫ μ (Ia + s1 (μ − μa)) dμ from μa to μb
                    acc += Ia * (μb^2 - μa^2) / 2 + s1 * ((μb^3 - μa^3) / 3 - μa * (μb^2 - μa^2) / 2)
                end
                Fν[k] += 2π * acc * (dv / c * ν0)
            end
        end
        for k in 1:n+1; Lline[i, k] = 4π * r[k]^2 * Fν[k]; end
        Threads.atomic_add!(nthick, 1)
    end
    for k in 1:n+1; Lr[k] = sum(view(Lline, :, k)); end
    Λ = [(Lr[k+1] - Lr[k]) / (4π / 3 * (r[k+1]^3 - r[k]^3)) for k in 1:n]
    LAST[] = (nbins = nbmax_used[], nthick = nthick[])
    Λ, Lr, Lline
end

# ---------------------------------------------------------------------------------------------------------------
# G6c: non-LTE level populations (H2Levels) coupled to the line transfer
#
# Per zone and line the populations give  κ0 = (hν/4π)(n_l B_lu − n_u B_ul),  S = n_u A/(n_l B_lu − n_u B_ul)
# = (2hν³/c²)/((g_u n_l)/(g_l n_u) − 1)  (stimulated emission included). The formal solution (tangent rays, comoving-frame
# profile shifts, source linear in τ) gives, besides the flux, the profile-weighted mean intensity J̄_k = ∫ φ_k J dν per
# zone and line, and the local (diagonal) operator Λ*_k of the same weighting.
#
# Local operator: in the short-characteristics step I_exit = I_entry e^{−Δτ} + S_exit ψ(Δτ) + S_entry (…), the coefficient
# of the exit-node source is ψ(Δτ) = w0 − w1/Δτ (→ Δτ/2 for Δτ ≪ 1, → 1 for Δτ ≫ 1). Λ* at a node is the mean over
# directions of this coefficient for the segment arriving at the node: inward rays arrive through the zone outside the
# node, outward rays through the zone inside it (the tangent ray: the same segment twice). It is accumulated with the same
# μ-quadrature and the same profile weights φ(v − μ u) as J̄, so J̄ = Λ[S] and Λ* are the same linear map restricted to the
# segment adjacent to each node; the zone value is the average of its two nodes. (The exact diagonal of Λ needs one
# unit-source solution per zone, n× the cost; the local operator is the usual approximate operator of ALI. The fixed point
# does not depend on the choice of Λ*.)
# ALI step (Rybicki & Hummer 1991): with J̄ = J̄_fs[S_old] + Λ*(S_new − S_old) the net rate of line u→l is
#   n_u A (1 − Λ*) − (n_l B_lu − n_u B_ul)(J̄_fs − Λ* S_old)   — linear in the populations, so each iteration is one
# statistical-equilibrium solve per zone with β = 1 − Λ*, J̄_eff = max(J̄_fs − Λ* S_old, 0) (H2Levels `β`, `Jbar`).
# ---------------------------------------------------------------------------------------------------------------
@inline ψloc(Δτ) = Δτ < 1e-4 ? 0.5Δτ : (e = exp(-Δτ); w0 = 1 - e; w1 = w0 - Δτ * e; w0 - w1 / Δτ)

"Source at the interfaces (linear in r between zone centres, centre/surface: zone values)."
function node_source(r, S)
    n = length(S); Si = zeros(n + 1); Si[1] = S[1]; Si[n+1] = S[n]
    for k in 2:n
        rc0 = 0.5 * (r[k-1] + r[k]); rc1 = 0.5 * (r[k] + r[k+1]); w = (r[k] - rc0) / (rc1 - rc0)
        Si[k] = (1 - w) * S[k-1] + w * S[k]
    end
    Si
end

"""
    transfer_line(r, u, vth, κ0, Si, ν0, xs, lam) -> (Fν, Jz, Λz)

Formal solution of one line on the offsets `xs` (observer-frame velocity): Fν[k] = ∫ flux over the line at the interfaces
[erg cm⁻² s⁻¹]; Jz[k] profile-averaged mean intensity of zone k (average of its two interfaces, each with the zone's
comoving profile) and Λz[k] the local operator in the same weighting (zero if `lam = false`).
"""
# LINESRC=zone (S6): the source is constant across each zone's segment (the zone's own S from its populations), not linear between
# interface values interpolated from neighbouring zones. Then a segment's emergent I = I e^{−Δτ} + S_zone (1 − e^{−Δτ}) exactly: a zone thick by
# itself absorbs the incoming intensity exponentially and re-emits at its own S, and emission κS and absorption κI in a zone use that
# zone's κ and S only (photon conservation per zone; with linear interpolation a thick zone next to a thin one whose S is large, e.g. near
# inversion, radiated at the neighbour's S).
const SRCZONE = Ref(get(ENV, "LINESRC", "linear") == "zone")
# LINEQUAD=p2 (S6): flux at node k as a trapezoid in μ² (∝ p², impact parameters fixed): every ray keeps its weight from node to node, so
# 4π r²F(k+1) − 4π r²F(k) is exactly the energy the segments of zone k took from or gave to the rays (linear in μ, the default, is more
# accurate for F itself in the thin limit but reweights the rays at each node: the difference of a large flux crossing dilute zones was noise)
const QUADP2 = Ref(get(ENV, "LINEQUAD", "mu") == "p2")
@inline shz(I, Δτ, S) = Δτ < 1e-4 ? I * (1 - Δτ) + Δτ * S : (e = exp(-Δτ); I * e + S * (1 - e))
@inline ψz(Δτ) = Δτ < 1e-4 ? Δτ : -expm1(-Δτ)
@inline ψlin(Δτ) = Δτ < 1e-4 ? 0.5Δτ : (e = exp(-Δτ); w0 = 1 - e; w0 - (w0 - Δτ * e) / Δτ)
"""
Ray geometry of `transfer_line`: everything that depends on r and u only, so it is the same for every line, frequency and ALI
iteration of a call. Ray j is tangent at r[j]; its nodes m = 1..n+2−j are the interfaces j..n+1.
dz[m, j]: path length of segment m; μu[m, j]: μ·u at the segment midpoint; μku[m, j]: μ·u at node m;
μa[k, j], μb[k, j]: direction cosines at interface k of the rays tangent at r[j+1] and r[j] (the angle quadrature).
"""
struct RayGeom
    dz::Matrix{Float64}; μu::Matrix{Float64}; μku::Matrix{Float64}; μa::Matrix{Float64}; μb::Matrix{Float64}
end
function raygeom(r, u)
    n = length(r) - 1
    dz = zeros(n + 1, n); μu = zeros(n + 1, n); μku = zeros(n + 1, n); μa = zeros(n + 1, n); μb = zeros(n + 1, n)
    zk = zeros(n + 1)
    for j in 1:n
        p = r[j]; nk = n + 2 - j
        for m in 1:nk; zk[m] = sqrt(max(r[j+m-1]^2 - p^2, 0.0)); end
        for m in 1:nk-1
            kz = j + m - 1
            zm = 0.5 * (zk[m] + zk[m+1]); rm = sqrt(p^2 + zm^2)
            w = clamp((rm - r[kz]) / (r[kz+1] - r[kz]), 0.0, 1.0)
            ur = (1 - w) * u[kz] + w * u[kz+1]
            μ = rm > 0 ? zm / rm : 0.0
            dz[m, j] = zk[m+1] - zk[m]; μu[m, j] = μ * ur
        end
        for m in 1:nk
            k = j + m - 1; μk = r[k] > 0 ? zk[m] / r[k] : 0.0
            μku[m, j] = μk * u[k]
        end
    end
    for k in 2:n+1, j in 1:k-1
        μa[k, j] = sqrt(max(1 - (r[j+1] / r[k])^2, 0.0)); μb[k, j] = sqrt(max(1 - (r[j] / r[k])^2, 0.0))
    end
    RayGeom(dz, μu, μku, μa, μb)
end

# (Nov 2026: the geometry is taken from `G`, the optical-depth prefactors are formed once per line and the work arrays once per call;
# the arithmetic of every value is unchanged, so the results are the same to the last bit. Before, the three per-ray arrays were
# allocated for every ray and frequency: 16 GB per ALI call at 150 zones and 100 thick lines.)
# `bg = (Bc, Bs, kq)` (incident radiation, see `incident`): the formal solution is done for I′ = I − I_bg, I_bg the continuum the ray
# would carry without the line: Bc (the isotropic background) on the inward leg and on the outward leg of the rays that miss the star,
# Bs (the photosphere) on the outward leg from node kq on of the rays j < kq, where I′ is set to 0 (the star is opaque: I = Bs there,
# whatever came in). I_bg is constant along a segment, so I′ obeys the same equation with the source S − I_bg. The returned Jz is that
# of I′ (the caller adds the background's own mean intensity, which it knows exactly); Fν is the flux of I′: the continuum that passes
# through has r²F independent of the node (LINEQUAD=p2), so the difference of 4πr²Fν across a zone is still what the zone gave the rays.
function transfer_line(r, u, vth, κ0, Si, ν0, xs, lam; Sz = nothing, G::RayGeom = raygeom(r, u), bg = nothing)
    n = length(vth); nb = length(xs); dv = step(xs)
    Fν = zeros(n + 1); Jz = zeros(n); Λz = zeros(n)
    ΔI = zeros(n + 1, n); QJ = zeros(n + 1, n, 2); QΛ = zeros(n + 1, n, 2)
    dz = G.dz; μu = G.μu; μku = G.μku; Ma = G.μa; Mb = G.μb
    spv = [sqrt(π) * vth[k] for k in 1:n]                      # profile normalisation: φ(v, k) = exp(−(v/vth[k])²) / spv[k]
    pre = zeros(n + 1, n)
    for j in 1:n, m in 1:n+1-j
        kz = j + m - 1; pre[m, j] = κ0[kz] * dz[m, j] / (sqrt(π) * vth[kz] / c * ν0)
    end
    dti = zeros(n + 1); dto = zeros(n + 1); Iin = zeros(n + 1); ψin = zeros(n + 1); Iout = zeros(n + 1); ψout = zeros(n + 1)
    quadp2 = QUADP2[]
    @inbounds for ib in 1:nb
        vb = xs[ib]
        fill!(ΔI, 0.0)
        for j in 1:n
            nk = n + 2 - j
            for m in 1:nk-1
                kz = j + m - 1; vt = vth[kz]; q = μu[m, j]; pr = pre[m, j]
                dti[m] = pr * exp(-((vb + q) / vt)^2)          # inward (dir = −1)
                dto[m] = pr * exp(-((vb - q) / vt)^2)          # outward
            end
            Iin[nk] = 0.0; ψin[nk] = 0.0; I = 0.0
            if bg !== nothing
                mq = (bg.kq > 0 && j < bg.kq) ? bg.kq - j + 1 : 0          # this ray's node at the photosphere (0: it misses the star)
                for m in nk-1:-1:1
                    if Sz === nothing
                        I = sh(I, dti[m], Si[j+m] - bg.Bc, Si[j+m-1] - bg.Bc); ψin[m] = ψlin(dti[m])
                    else
                        I = shz(I, dti[m], Sz[j+m-1] - bg.Bc); ψin[m] = ψz(dti[m])
                    end
                    Iin[m] = I
                end
                Iout[1] = Iin[1]; ψout[1] = ψin[1]; I = Iin[1]
                for m in 2:nk
                    bq = (mq > 0 && m > mq) ? bg.Bs : bg.Bc
                    if Sz === nothing
                        I = sh(I, dto[m-1], Si[j+m-2] - bq, Si[j+m-1] - bq); ψout[m] = ψlin(dto[m-1])
                    else
                        I = shz(I, dto[m-1], Sz[j+m-2] - bq); ψout[m] = ψz(dto[m-1])
                    end
                    m == mq && (I = 0.0; ψout[m] = 0.0)
                    Iout[m] = I
                end
            elseif Sz === nothing
                for m in nk-1:-1:1
                    Δτ = dti[m]; Sin = Si[j+m]; Sout = Si[j+m-1]
                    if Δτ < 1e-4                               # `sh` and `ψloc` with their common exponential taken once
                        I = I * (1 - Δτ) + 0.5Δτ * (Sin + Sout); ψin[m] = 0.5Δτ
                    else
                        e = exp(-Δτ); w0 = 1 - e; w1 = w0 - Δτ * e
                        I = I * e + Sout * (w0 - w1 / Δτ) + Sin * (w1 / Δτ); ψin[m] = w0 - w1 / Δτ
                    end
                    Iin[m] = I
                end
                Iout[1] = Iin[1]; ψout[1] = ψin[1]; I = Iin[1]
                for m in 2:nk
                    Δτ = dto[m-1]; Sin = Si[j+m-2]; Sout = Si[j+m-1]
                    if Δτ < 1e-4
                        I = I * (1 - Δτ) + 0.5Δτ * (Sin + Sout); ψout[m] = 0.5Δτ
                    else
                        e = exp(-Δτ); w0 = 1 - e; w1 = w0 - Δτ * e
                        I = I * e + Sout * (w0 - w1 / Δτ) + Sin * (w1 / Δτ); ψout[m] = w0 - w1 / Δτ
                    end
                    Iout[m] = I
                end
            else
                for m in nk-1:-1:1
                    I = shz(I, dti[m], Sz[j+m-1]); Iin[m] = I; ψin[m] = ψz(dti[m])
                end
                Iout[1] = Iin[1]; ψout[1] = ψin[1]; I = Iin[1]
                for m in 2:nk
                    I = shz(I, dto[m-1], Sz[j+m-2]); Iout[m] = I; ψout[m] = ψz(dto[m-1])
                end
            end
            for m in 1:nk
                k = j + m - 1; ΔI[k, j] = m == 1 ? 0.0 : Iout[m] - Iin[m]
                q = μku[m, j]
                for s in 1:2
                    zs = k - 2 + s                             # the zone inside (s = 1) and outside (s = 2) interface k
                    (zs < 1 || zs > n) && continue
                    pin = exp(-((vb + q) / vth[zs])^2) / spv[zs]; pout = exp(-((vb - q) / vth[zs])^2) / spv[zs]
                    QJ[k, j, s] = pin * Iin[m] + pout * Iout[m]
                    lam && (QΛ[k, j, s] = pin * ψin[m] + pout * ψout[m])
                end
            end
        end
        for k in 2:n+1
            acc = 0.0
            for j in 1:k-1
                μa = Ma[k, j]; μb = Mb[k, j]
                Ia = j + 1 <= n ? ΔI[k, j+1] : 0.0; Ib = ΔI[k, j]
                d = μb - μa; d <= 0 && continue
                if quadp2                        # trapezoid in μ² = 1 − p²/r_k²: ray weights ∝ Δp², the same at every node
                    acc += 0.25 * (Ia + Ib) * (μb^2 - μa^2)
                else
                    s1 = (Ib - Ia) / d
                    acc += Ia * (μb^2 - μa^2) / 2 + s1 * ((μb^3 - μa^3) / 3 - μa * (μb^2 - μa^2) / 2)
                end
            end
            Fν[k] += 2π * acc * (dv / c * ν0)
        end
        # mean intensity at the nodes: J = ½ ∫₀¹ (Q_in + Q_out) dμ, trapezoid in μ between the rays (Q includes the profile)
        for k in 1:n+1, s in 1:2
            zs = k - 2 + s
            (zs < 1 || zs > n) && continue
            if k == 1
                Jz[zs] += 0.25 * QJ[1, 1, s] * dv; lam && (Λz[zs] += 0.25 * QΛ[1, 1, s] * dv)
                continue
            end
            aJ = 0.0; aΛ = 0.0
            for j in 1:k-1
                μa = Ma[k, j]; μb = Mb[k, j]
                d = μb - μa; d <= 0 && continue
                qa = j + 1 <= n ? QJ[k, j+1, s] : 0.0
                aJ += 0.5 * (qa + QJ[k, j, s]) * d
                if lam
                    la = j + 1 <= n ? QΛ[k, j+1, s] : 0.0
                    aΛ += 0.5 * (la + QΛ[k, j, s]) * d
                end
            end
            Jz[zs] += 0.25 * aJ * dv; lam && (Λz[zs] += 0.25 * aΛ * dv)
        end
    end
    Fν, Jz, Λz
end

# LINEFAST=1 (Nov 2026; default off): the formal solution of `transfer_line` for the zone-constant source (LINESRC=zone), reorganised.
# The method, the rays, the frequency bins and the quadratures are those of `transfer_line`; the results agree to round-off, not to
# the last bit. What differs:
#   - the loops run ray → segment → frequency bin (before: bin → ray → segment), so the innermost loop is over independent bins
#     and is vectorised (`@simd`), with a branch-free exponential (`expk`, within 1 ulp of `exp`);
#   - a segment acts on a bin only where the bin lies within the line profile of the segment's zone: |v − μu| ≤ v_th √(38 + ln max(τ₀, 1))
#     (τ₀ the line-centre depth of the segment), outside of which its optical depth is below 3e-17; a node's profile weight is taken
#     within |v − μu| ≤ √38 v_th (weight below 3e-17 of the peak). Segments without opacity are skipped (I passes unchanged);
#   - the sums over the frequency bins are taken first (per ray and node) and the angle quadratures once per line (before: the
#     quadratures per bin): another order of summation;
#   - the work arrays come from a pool of one set per task (`FastWS`) instead of being allocated for every line.
const FAST = Ref(get(ENV, "LINEFAST", "0") == "1")
const FASTXC2 = 38.0
const EXPC = ntuple(i -> 1.0 / factorial(i), 13)             # (e^r − 1)/r = Σ r^(i−1)/i!
"`(exp(x), expm1(x))` for x ≤ 0 without branches or tables (vectorisable): x = k ln 2 + r, |r| ≤ ½ ln 2, Taylor series to r¹³. Exactly 0 below −708."
@inline function expk(x)
    xc = max(x, -708.0)
    kf = round(xc * 1.4426950408889634)
    rr = muladd(kf, -0.6931471803691238, xc); rr = muladd(kf, -1.9082149292705877e-10, rr)
    em = rr * evalpoly(rr, EXPC)
    s = reinterpret(Float64, (Int64(unsafe_trunc(Int32, kf)) + 1023) << 52)
    e = muladd(s, em, s)
    (ifelse(x < -708.0, 0.0, e), ifelse(kf == 0.0, em, e - 1.0))
end
struct FastWS
    Iin::Matrix{Float64}; I::Vector{Float64}; ψ::Vector{Float64}; xv::Vector{Float64}
    ΔI::Matrix{Float64}; QJ::Array{Float64,3}; QΛ::Array{Float64,3}; pre::Matrix{Float64}
    Fν::Vector{Float64}; Jz::Vector{Float64}; Λz::Vector{Float64}
end
FastWS(n, nb) = FastWS(zeros(nb, n + 1), zeros(nb), zeros(nb), zeros(nb), zeros(n + 1, n), zeros(n + 1, n, 2), zeros(n + 1, n, 2),
                       zeros(n + 1, n), zeros(n + 1), zeros(n), zeros(n))
function fastpool(n, nb)
    ch = Channel{FastWS}(Threads.nthreads())
    for _ in 1:Threads.nthreads(); put!(ch, FastWS(n, nb)); end
    ch
end
# the bins lo:hi of the grid x1 + (ib − 1)/rdv that can lie within c0 ± w (one bin of margin; all bins if c0 or w is not a number)
@inline function binwin(c0, w, x1, rdv, nb)
    a = (c0 - w - x1) * rdv + 1.0; b = (c0 + w - x1) * rdv + 1.0
    (a == a && b == b) || return (1, nb)
    lo = a <= 1.0 ? 1 : a >= nb + 1.0 ? nb + 1 : floor(Int, a)
    hi = b >= nb ? nb : b < 0.0 ? 0 : ceil(Int, b)
    (lo, hi)
end
# one segment for all the bins it acts on: I ← I e^{−Δτ} + S (1 − e^{−Δτ}), ψ = 1 − e^{−Δτ}, Δτ = pr φ̂((v + qs)/vt) (`shz`, `ψz`)
@inline function fast_step!(I, ψ, xv, pr, qs, vt, S, x1, rdv, nb)
    pr == 0 && return (1, 0)
    lo, hi = binwin(-qs, vt * sqrt(FASTXC2 + log(max(abs(pr), 1.0))), x1, rdv, nb)
    @inbounds @simd for ib in lo:hi
        a = (xv[ib] + qs) / vt
        Δτ = pr * expk(-(a * a))[1]
        e, em = expk(-Δτ)
        Io = I[ib]
        I[ib] = ifelse(Δτ < 1e-4, Io * (1 - Δτ) + Δτ * S, Io * e + S * (1 - e))
        ψ[ib] = ifelse(Δτ < 1e-4, Δτ, -em)
    end
    (lo, hi)
end
# Σ over the bins of φ̂((v + qs)/vt) I and of φ̂ ψ (the profile's normalisation is applied by the caller)
@inline function fast_node(I, ψ, xv, qs, vt, x1, rdv, nb, lam)
    lo, hi = binwin(-qs, vt * sqrt(FASTXC2), x1, rdv, nb)
    sJ = 0.0; sΛ = 0.0
    if lam
        @inbounds @simd for ib in lo:hi
            a = (xv[ib] + qs) / vt; p = expk(-(a * a))[1]
            sJ += p * I[ib]; sΛ += p * ψ[ib]
        end
    else
        @inbounds @simd for ib in lo:hi
            a = (xv[ib] + qs) / vt; p = expk(-(a * a))[1]
            sJ += p * I[ib]
        end
    end
    (sJ, sΛ)
end
"`transfer_line` with `Sz` (zone-constant source), reorganised (see LINEFAST). Returns (Fν, Jz, Λz) as views of the work set `W`: copy before releasing it."
function transfer_line_fast!(W::FastWS, r, u, vth, κ0, ν0, xs, lam::Bool, Sz, G::RayGeom, bg)
    n = length(vth); nb = length(xs); dv = step(xs); rdv = 1 / dv
    Fν = W.Fν; Jz = W.Jz; Λz = W.Λz; ΔI = W.ΔI; QJ = W.QJ; QΛ = W.QΛ; pre = W.pre; Iin = W.Iin; I = W.I; ψ = W.ψ; xv = W.xv
    @assert size(ΔI, 2) == n && length(I) >= nb
    fill!(Fν, 0.0); fill!(Jz, 0.0); fill!(Λz, 0.0); fill!(ΔI, 0.0); fill!(QJ, 0.0); fill!(QΛ, 0.0)
    dz = G.dz; μu = G.μu; μku = G.μku; Ma = G.μa; Mb = G.μb
    @inbounds for ib in 1:nb; xv[ib] = xs[ib]; end
    x1 = xv[1]
    @inbounds for j in 1:n, m in 1:n+1-j
        kz = j + m - 1; pre[m, j] = κ0[kz] * dz[m, j] / (sqrt(π) * vth[kz] / c * ν0)
    end
    hasbg = bg !== nothing; Bc = hasbg ? Float64(bg.Bc) : 0.0; Bs = hasbg ? Float64(bg.Bs) : 0.0; kq = hasbg ? Int(bg.kq) : 0
    quadp2 = QUADP2[]
    @inbounds for j in 1:n
        nk = n + 2 - j
        mq = (kq > 0 && j < kq) ? kq - j + 1 : 0                  # this ray's node at the photosphere (0: it misses the star)
        for ib in 1:nb; I[ib] = 0.0; ψ[ib] = 0.0; Iin[ib, nk] = 0.0; end              # nothing comes in at the outer node
        for m in nk-1:-1:1                                         # inwards (I′ = 0 at the outer node)
            kz = j + m - 1
            lo, hi = fast_step!(I, ψ, xv, pre[m, j], μu[m, j], vth[kz], Sz[kz] - Bc, x1, rdv, nb)
            @simd for ib in 1:nb; Iin[ib, m] = I[ib]; end
            qn = μku[m, j]
            for s in 1:2
                zs = kz - 2 + s                                    # the zone inside (s = 1) and outside (s = 2) interface k = kz
                (zs < 1 || zs > n) && continue
                sJ, sΛ = fast_node(I, ψ, xv, qn, vth[zs], x1, rdv, nb, lam)
                if m == 1                                          # the turning point: I_out = I_in, ψ_out = ψ_in
                    oJ, oΛ = fast_node(I, ψ, xv, -qn, vth[zs], x1, rdv, nb, lam); sJ += oJ; sΛ += oΛ
                end
                spv = sqrt(π) * vth[zs]
                QJ[kz, j, s] = sJ / spv; QΛ[kz, j, s] = sΛ / spv
            end
            for ib in lo:hi; ψ[ib] = 0.0; end
        end
        for m in 2:nk                                              # outwards
            kz = j + m - 2; k = j + m - 1
            bq = (mq > 0 && m > mq) ? Bs : Bc
            lo, hi = fast_step!(I, ψ, xv, pre[m-1, j], -μu[m-1, j], vth[kz], Sz[kz] - bq, x1, rdv, nb)
            if m == mq                                             # the star is opaque: I′ = 0 from here on
                for ib in 1:nb; I[ib] = 0.0; end
                for ib in lo:hi; ψ[ib] = 0.0; end
            end
            d = 0.0
            @simd for ib in 1:nb; d += I[ib] - Iin[ib, m]; end
            ΔI[k, j] = d
            qn = μku[m, j]
            for s in 1:2
                zs = k - 2 + s
                (zs < 1 || zs > n) && continue
                sJ, sΛ = fast_node(I, ψ, xv, -qn, vth[zs], x1, rdv, nb, lam)
                spv = sqrt(π) * vth[zs]
                QJ[k, j, s] += sJ / spv; QΛ[k, j, s] += sΛ / spv
            end
            for ib in lo:hi; ψ[ib] = 0.0; end
        end
    end
    @inbounds for k in 2:n+1                                       # the quadratures of `transfer_line`, once, on the sums over the bins
        acc = 0.0
        for j in 1:k-1
            μa = Ma[k, j]; μb = Mb[k, j]
            Ia = j + 1 <= n ? ΔI[k, j+1] : 0.0; Ib = ΔI[k, j]
            d = μb - μa; d <= 0 && continue
            if quadp2
                acc += 0.25 * (Ia + Ib) * (μb^2 - μa^2)
            else
                s1 = (Ib - Ia) / d
                acc += Ia * (μb^2 - μa^2) / 2 + s1 * ((μb^3 - μa^3) / 3 - μa * (μb^2 - μa^2) / 2)
            end
        end
        Fν[k] = 2π * acc * (dv / c * ν0)
    end
    @inbounds for k in 1:n+1, s in 1:2
        zs = k - 2 + s
        (zs < 1 || zs > n) && continue
        if k == 1
            Jz[zs] += 0.25 * QJ[1, 1, s] * dv; lam && (Λz[zs] += 0.25 * QΛ[1, 1, s] * dv)
            continue
        end
        aJ = 0.0; aΛ = 0.0
        for j in 1:k-1
            μa = Ma[k, j]; μb = Mb[k, j]
            d = μb - μa; d <= 0 && continue
            qa = j + 1 <= n ? QJ[k, j+1, s] : 0.0
            aJ += 0.5 * (qa + QJ[k, j, s]) * d
            if lam
                la = j + 1 <= n ? QΛ[k, j+1, s] : 0.0
                aΛ += 0.5 * (la + QΛ[k, j, s]) * d
            end
        end
        Jz[zs] += 0.25 * aJ * dv; lam && (Λz[zs] += 0.25 * aΛ * dv)
    end
    Fν, Jz, Λz
end

# EPGRAD (escape-probability mode, the Sobolev part): which radial velocity gradient a zone's lines escape through.
#   zone: the difference of the zone's two interface velocities, |Δu|/Δr, and |ū|/r of the zone (the method up to Nov 2026).
#   thermal (default): v_th/ℓ, ℓ the distance from the zone centre to the nearest radius (inwards or outwards) where the velocity differs from
#     the zone's by one thermal width v_th; no radial term if it never does. Where the flow is resolved and supersonic this is |du/dr|;
#     in a subsonic core (Δu across a zone ≪ v_th, Sobolev length ≫ zone) a line photon does not see the zone's own Δu, and |Δu|/Δr
#     carries the zone-to-zone velocity noise amplified by r/Δr into β and so into the cooling: a zone compressing slightly faster
#     cooled more, lost pressure and compressed faster still (growth time a few zone sound crossings; the central noise of full1
#     from n_c = 1e10, when the lines turn thick; the ALI populations, from the real transfer, do not show it).
#     The tangential term is taken over a thermal width too (not |u|/r of the zone): in the core |u|/r is a ratio of two small numbers
#     and a 10 m/s sound wave (Mach 0.003) changed it, and with it the EP cooling of the central zones, by 40 % from one step to the next
#     (t_f7, steps 3742 → 3744, n_c = 1.7e11: Λ_EP 2.5e-11 → 3.6e-11 while the ALI cooling moved 5.65e-11 → 6.00e-11).
#   flow: no radial term at all, only the tangential one |u|/r (diagnostic).
const EPGRAD = Ref(Symbol(get(ENV, "EPGRAD", "thermal")))
"LINESCHED=zone (default): the level populations of the zones are distributed one zone at a time; `chunk`: `Threads.@threads` blocks. Same results."
const ZSCHED = Ref(get(ENV, "LINESCHED", "zone") == "zone")
"Exponential integral E₂(x) = e^{−x} − x E₁(x), x ≥ 0 (series for x ≤ 1, continued fraction above)."
function E2(x)
    x <= 0 && return 1.0
    if x <= 1
        E1 = -0.5772156649015329 - log(x); t = 1.0
        for k in 1:30; t *= -x / k; E1 -= t / k; end
        return exp(-x) - x * E1
    end
    # Lentz continued fraction for E₁
    b = x + 1.0; cc = 1e300; d = 1 / b; f = d
    for k in 1:200
        a = -k^2.0; b += 2.0; d = 1 / (a * d + b); cc = b + a / cc; Δ = cc * d; f *= Δ; abs(Δ - 1) < 1e-14 && break
    end
    exp(-x) - x * (f * exp(-x))
end
# EPTAB=1 (Nov 2026; default off): the static part of the escape probability of `epβ` from a table. That part is a function of one
# variable, F(τ) = Σ_x w_x E₂(τ e^{−x²}) over the 21-point frequency quadrature (β_static = ½[F(τ_out) + F(τ_in)]), and was summed
# with 42 E₂ evaluations (series or continued fraction) for every line and zone on every escape-probability call: after LINEFAST
# the largest single item of the line transfer between ALI refreshes. The table holds the same sum on a grid in ln τ (step 1/128,
# 1e-10 ≤ τ ≤ 1e10; 4-point Lagrange interpolation, relative error < 1e-9; outside the range the sum itself). The table is used
# only to DECIDE: β is the larger of the static and the Sobolev term, and where the Sobolev term exceeds the tabulated static
# term by more than the table's error (margin 1e-6) the sum is not needed; elsewhere the sum is evaluated as before. The result
# is therefore the same to the last bit. (In the accreting envelope the Sobolev term wins everywhere: 0 of 4.2e5 evaluations
# on four production snapshots needed the sum.)
const EPTAB = Ref(get(ENV, "EPTAB", "0") == "1")
const EPXG = 0.0:0.25:5.0
const EPWG = [(x == 0 || x == 5.0 ? 0.5 : 1.0) * 0.25 * 2 / sqrt(π) * exp(-x^2) for x in EPXG]
const EPPH = [exp(-x^2) for x in EPXG]
function epF_sum(τ)
    b = 0.0
    for j in eachindex(EPWG); b += EPWG[j] * E2(τ * EPPH[j]); end
    b
end
const EPTL0 = log(1e-10); const EPTH = 1 / 128; const EPTN = ceil(Int, (log(1e10) - EPTL0) / EPTH) + 4
const EPTABLE = Ref(Float64[])
"The table of F at ln τ = EPTL0 + (j − 2) EPTH, j = 1..EPTN (built on first use; call before a threaded loop)."
eptable() = (isempty(EPTABLE[]) && (EPTABLE[] = [epF_sum(exp(EPTL0 + (j - 2) * EPTH)) for j in 1:EPTN]); EPTABLE[])
@inline function epF(tab, τ)
    τ > 0 || return epF_sum(0.0)
    s = (log(τ) - EPTL0) / EPTH
    (s >= 0 && s < EPTN - 4) || return epF_sum(τ)
    j = floor(Int, s); t = s - j
    @inbounds f0 = tab[j+1]; @inbounds f1 = tab[j+2]; @inbounds f2 = tab[j+3]; @inbounds f3 = tab[j+4]
    -t * (t - 1) * (t - 2) / 6 * f0 + (t + 1) * (t - 1) * (t - 2) / 2 * f1 - (t + 1) * t * (t - 2) / 2 * f2 + (t + 1) * t * (t - 1) / 6 * f3
end

"Max relative change of the level populations (levels with fraction > pmin in either iterate)."
function popchange(P, Q; pmin = 1e-6)
    d = 0.0
    for idx in eachindex(P)
        m = max(P[idx], Q[idx]); m > pmin || continue
        d = max(d, abs(P[idx] - Q[idx]) / m)
    end
    d
end

function check_levels(L::LineList, LV)
    length(LV.A) == length(L.A) || error("H2Lines/H2Levels: different line counts")
    ν = [(LV.E[LV.iu[i]] - LV.E[LV.il[i]]) / h for i in eachindex(L.A)]
    (all(isapprox.(ν, L.ν; rtol = 1e-9)) && all(LV.A .== L.A) && all(isapprox.(LV.g[LV.iu], L.g[L.iu])) &&
        all(isapprox.(LV.g[LV.il], L.g[L.il]))) || error("H2Lines/H2Levels: level or line mismatch")
    true
end

"""
    incident(inc, r, ν) -> (kq, Bc, Bs, Jbg)

Incident radiation for the level populations (keyword `inc = (Tbg, Tstar, Rph)` of `cooling`, NLTE paths; driver switch LINEINC):
an isotropic blackbody background of temperature Tbg (the CMB) and the continuum of the photosphere, B_ν(Tstar) from the sphere of
radius Rph, seen diluted. Without it the lines see no radiation but their own, so a level can only be excited by collisions; with it
the levels are also pumped (and de-excited by stimulated emission), and the gas gains or loses energy through the collisions that
follow. kq: the node of the grid nearest Rph (0: no star; then the background only); Bc, Bs: B_ν of the two per line;
Jbg[i, k] = Bc + W_k (Bs − Bc), W_k the mean of the dilution factor W(r) = ½[1 − √(1 − (R/r)²)] at the zone's two interfaces (½ at
kq, 0 inside), the same weighting as the J̄ of `transfer_line`.
Not included: absorption of the continuum between the photosphere and the zone (the grey depth there is < 2/3), the line photons the
gas takes from or adds to the continuum are not debited from the grey radiation field, and the photosphere is a blackbody of T_eff.
"""
function incident(inc, r, ν)
    n = length(r) - 1; nl = length(ν)
    Tb = Float64(inc.Tbg); Ts = Float64(inc.Tstar); Rp = Float64(inc.Rph)
    Bc = [Tb > 0 ? Bν(ν[i], Tb) : 0.0 for i in 1:nl]; Bs = copy(Bc); kq = 0
    if Ts > Tb && Rp > 0
        kn = argmin([abs(r[k] - Rp) for k in 1:n+1])
        2 <= kn <= n && (kq = kn; Bs = [Bν(ν[i], Ts) for i in 1:nl])
    end
    Wn = zeros(n + 1)
    if kq > 0
        for k in kq:n+1
            x2 = (r[kq] / r[k])^2; Wn[k] = k == kq ? 0.5 : 0.5 * x2 / (1 + sqrt(1 - x2))
        end
    end
    Jbg = [Bc[i] + 0.5 * (Wn[k] + Wn[k+1]) * (Bs[i] - Bc[i]) for i in 1:nl, k in 1:n]
    (; kq, Bc, Bs, Jbg)
end

"""
    cooling_nlte(L, r, u, T, nH2, nlte; mode = :ali | :ep, ntol, maxit, init, lam)

See `cooling` (keyword `nlte = (LV, nH, nHe, nHp)`). Returns (Λ, Lr, Lline) as the LTE path (Λ = −∇·F per zone from the
flux of the formal solution of the last iterate); `LAST[]` gets `iters` (ALI iterations; the EP start is counted in `nep`),
`converged`, `hist` (max relative population change per iteration), `Λrate` (net radiative rate Σ hν[n_u A − (n_l B_lu −
n_u B_ul) J̄] per zone, ALI) and `pops` (63 × n).
"""
function cooling_nlte(L::LineList, r, u, T, nH2, nlte; nbin = 30, tol = 0.01, nbinmax = 200, cover = 0.999, mode = :ali,
                      ntol = 1e-3, maxit = 60, init = :ep, lam = :diag, inc = nothing)
    LV = nlte.LV; HV = Main.H2Levels; check_levels(L, LV)
    n = length(T); @assert length(r) == n + 1
    nl = length(L.A); nlev = length(LV.E)
    vth = [sqrt(2kB * T[k] / mH2) for k in 1:n]
    τthin = 2tol
    iu = LV.iu; il = LV.il; Bul = LV.Bul; Blu = LV.Blu; ν = L.ν
    P = zeros(nlev, n); β = ones(nl, n); Jeff = zeros(nl, n)
    # incident radiation (`incident`): in the escape-probability populations J̄ = (1 − β) S + β J_bg, in the ALI ones J̄ = J̄′ + J_bg
    bgd = inc === nothing ? nothing : incident(inc, r, ν)
    function solve!(Pn)
        if ZSCHED[]                                   # one zone at a time from a shared counter: the ortho:para iteration count varies between zones
            next = Threads.Atomic{Int}(1)
            @sync for _ in 1:Threads.nthreads()
                Threads.@spawn begin
                    ws = HV.PopWS(LV)                 # work arrays of this task
                    while true
                        kz = Threads.atomic_add!(next, 1)
                        kz > n && break
                        Pn[:, kz] .= HV.populations(LV, T[kz], nlte.nH[kz], nH2[kz], nlte.nHe[kz], nlte.nHp[kz]; β = view(β, :, kz), Jbar = view(Jeff, :, kz), ws)
                    end
                end
            end
        else
            Threads.@threads for k in 1:n
                Pn[:, k] .= HV.populations(LV, T[k], nlte.nH[k], nH2[k], nlte.nHe[k], nlte.nHp[k]; β = β[:, k], Jbar = Jeff[:, k])
            end
        end
        Pn
    end
    # opacity and source from populations
    function opac(Pc)
        κ0 = zeros(nl, n); S = zeros(nl, n); nm = zeros(Int, n)
        Threads.@threads for k in 1:n                 # (S6: threaded; zones independent)
            for i in 1:nl
                pu = Pc[iu[i], k]; pl = Pc[il[i], k]; den = pl * Blu[i] - pu * Bul[i]
                if den > 0
                    κ0[i, k] = h * ν[i] / (4π) * nH2[k] * den; S[i, k] = pu * LV.A[i] / den
                else
                    nm[k] += pu > 0 ? 1 : 0; S[i, k] = Bν(ν[i], T[k])
                end
            end
        end
        κ0, S, sum(nm)
    end
    # escape probability per line and zone from the current per-line depths (EP mode): static part = profile-averaged
    # two-direction slab escape ½[E₂(τ_out φ̂) + E₂(τ_in φ̂)], φ̂ = e^{−x²}, with τ_out the line-centre depth along the radius to
    # the surface and τ_in that along the inward radius through the centre to the surface; Sobolev part (velocity gradient
    # g = max(|du/dr|, |u/r|)): (1 − e^{−τ_S})/τ_S, τ_S = κ0 c/(ν0 g); β = the larger of the two (a photon leaves by the
    # easier route). β = 1 where κ0 ≤ 0.
    # EPGRAD=thermal: the radial gradient over the distance in which the velocity changes by one thermal width (see EPGRAD)
    gth = zeros(n); gtn = zeros(n)
    if EPGRAD[] == :thermal
        for k in 1:n
            rc = 0.5 * (r[k] + r[k+1]); uc = 0.5 * (u[k] + u[k+1]); ℓ = Inf
            for (js, rev) in ((k+1:n+1, false), (k:-1:1, true))
                rp = rc; dp = 0.0
                for j in js
                    dj = abs(u[j] - uc)
                    if dj >= vth[k]
                        ℓ = min(ℓ, abs(rp + (r[j] - rp) * (vth[k] - dp) / (dj - dp) - rc)); break
                    end
                    rp = r[j]; dp = dj
                end
            end
            gth[k] = vth[k] / ℓ
            # the same along the tangent: at path length s from the zone centre the line-of-sight velocity is u(r') s/r', r' = √(rc² + s²)
            # (zero at the tangent point); for homologous flow u = −a r this gives a, as |u|/r does, and |u|/r where |u| ≫ v_th
            sp = 0.0; vp = 0.0
            for j in k+1:n+1
                sj = sqrt(max(r[j]^2 - rc^2, 0.0)); vj = abs(u[j]) * sj / r[j]
                if vj >= vth[k]
                    gtn[k] = vth[k] / (sp + (sj - sp) * (vth[k] - vp) / (vj - vp)); break
                end
                sp = sj; vp = vj
            end
        end
    end
    function epβ(κ0)
        βn = ones(nl, n)
        xg = 0.0:0.25:5.0; wg = [(x == 0 || x == 5.0 ? 0.5 : 1.0) * 0.25 * 2 / sqrt(π) * exp(-x^2) for x in xg]
        tab = EPTAB[] ? eptable() : nothing
        Threads.@threads for i in 1:nl                # (S6: threaded; lines independent)
            τs = zeros(n); acc = 0.0
            for k in n:-1:1
                dτ = max(κ0[i, k], 0.0) / (sqrt(π) * vth[k] / c * ν[i]) * (r[k+1] - r[k])
                τs[k] = acc + 0.5dτ; acc += dτ
            end
            τtot = acc
            for k in 1:n
                κ0[i, k] > 0 || continue
                τo = τs[k]; τi = (τtot - τs[k]) + τtot
                rc = 0.5 * (r[k] + r[k+1])
                g = EPGRAD[] == :thermal ? max(gth[k], gtn[k]) :
                    max(EPGRAD[] == :flow ? 0.0 : abs(u[k+1] - u[k]) / (r[k+1] - r[k]), abs(0.5 * (u[k] + u[k+1])) / max(rc, 1e-300))
                bS = -1.0
                if g > 0
                    τS = κ0[i, k] * c / (ν[i] * g); bS = τS < 1e-6 ? 1 - τS / 2 : -expm1(-τS) / τS
                end
                b = 0.0
                if tab !== nothing && bS >= (1 + 1e-6) * 0.5 * (epF(tab, τo) + epF(tab, τi))
                    b = bS                                # (EPTAB) the static term cannot win: its sum is not needed
                else
                    for (x, w) in zip(xg, wg); φ̂ = exp(-x^2); b += w * 0.5 * (E2(τo * φ̂) + E2(τi * φ̂)); end
                    g > 0 && (b = max(b, bS))
                end
                βn[i, k] = min(b, 1.0)
            end
        end
        βn
    end
    hist = Float64[]; nep = 0; iters = 0; converged = false
    bgd !== nothing && (Jeff .= bgd.Jbg)
    solve!(P)
    if mode == :ep || init == :ep
        while nep < maxit
            κ0, S, _ = opac(P); β .= epβ(κ0); bgd !== nothing && (Jeff .= β .* bgd.Jbg); Pn = solve!(copy(P)); nep += 1
            d = popchange(P, Pn); P .= Pn; mode == :ep && push!(hist, d); d < ntol && break
        end
        mode == :ep && (converged = !isempty(hist) && hist[end] < ntol)
    end
    nbmax = 0; nthick = 0; Lline = zeros(nl, n + 1); Λrate = zeros(n); nmaser = 0
    Λem = zeros(n)                # (incident radiation, EP mode) the emission alone, n_u A β hν summed over the lines [erg cm⁻³ s⁻¹]
    if mode == :ep
        κ0, S, nmaser = opac(P)
        # local emission n_u A β hν accumulated over the volume
        if bgd === nothing
        for i in 1:nl, k in 1:n
            Lline[i, k+1] = Lline[i, k] + nH2[k] * P[iu[i], k] * LV.A[i] * β[i, k] * h * ν[i] * 4π / 3 * (r[k+1]^3 - r[k]^3)
        end
        else                      # net of the background photons absorbed: β hν [n_u A − (n_l B_lu − n_u B_ul) J_bg]
            for i in 1:nl, k in 1:n
                pu = P[iu[i], k]; pl = P[il[i], k]
                em = nH2[k] * pu * LV.A[i] * β[i, k] * h * ν[i]
                ab = nH2[k] * (pl * Blu[i] - pu * Bul[i]) * β[i, k] * bgd.Jbg[i, k] * h * ν[i]
                Λem[k] += em; Lline[i, k+1] = Lline[i, k] + (em - ab) * 4π / 3 * (r[k+1]^3 - r[k]^3)
            end
        end
        iters = nep
    else
        Jfs = zeros(nl, n); Lam = zeros(nl, n); β .= 1; Jeff .= 0
        geom = raygeom(r, u)
        wpool = FAST[] && SRCZONE[] ? fastpool(n, nbinmax) : nothing          # LINEFAST: one work set per task
        for it in 1:maxit
            κ0, S, nmaser = opac(P)
            fill!(Jfs, 0.0); fill!(Lam, 0.0); fill!(Lline, 0.0)
            nb_t = Threads.Atomic{Int}(0); nth_t = Threads.Atomic{Int}(0)
            Threads.@threads :greedy for i in 1:nl                  # one line per task: the thick lines are neighbours in the list
                ν0 = ν[i]; κl = κ0[i, :]
                τc = sum(κl[k] / (sqrt(π) * vth[k] / c * ν0) * (r[k+1] - r[k]) for k in 1:n)
                if τc < τthin
                    acc = 0.0
                    if bgd === nothing
                    for k in 1:n
                        acc += nH2[k] * P[iu[i], k] * LV.A[i] * h * ν0 * 4π / 3 * (r[k+1]^3 - r[k]^3); Lline[i, k+1] = acc
                    end
                    else          # thin line in the incident field: J̄ = J_bg, net emission
                        for k in 1:n
                            pu = P[iu[i], k]; pl = P[il[i], k]; Jfs[i, k] = bgd.Jbg[i, k]
                            acc += nH2[k] * (pu * LV.A[i] - (pl * Blu[i] - pu * Bul[i]) * bgd.Jbg[i, k]) * h * ν0 * 4π / 3 * (r[k+1]^3 - r[k]^3)
                            Lline[i, k+1] = acc
                        end
                    end
                    continue
                end
                Si = node_source(r, S[i, :])
                em = [nH2[k] * P[iu[i], k] * LV.A[i] * (r[k+1]^3 - r[k]^3) for k in 1:n]
                oc = [κl[k] / (vth[k] / c * ν0) * (r[k+1] - r[k]) for k in 1:n]
                act = falses(n)
                for w in (em, oc)
                    tot = sum(abs, w); tot > 0 || continue
                    p_ = sortperm(abs.(w); rev = true); acc_ = 0.0
                    for k in p_; act[k] = true; acc_ += abs(w[k]); acc_ >= cover * tot && break; end
                end
                any(act) || (act .= true)
                ka = findall(act)
                vspan = 3maximum(vth[ka]) + maximum(max(abs(u[k]), abs(u[k+1])) for k in ka)
                nb = clamp(ceil(Int, 2vspan / (0.5minimum(vth[ka]))) + 1, nbin, nbinmax)
                Threads.atomic_max!(nb_t, nb)
                xs = range(-vspan, vspan, length = nb)
                bgl = bgd === nothing ? nothing : (Bc = bgd.Bc[i], Bs = bgd.Bs[i], kq = bgd.kq)
                W = wpool === nothing ? nothing : take!(wpool)
                Fν, Jz, Λz = W === nothing ? transfer_line(r, u, vth, κl, Si, ν0, xs, lam == :diag; Sz = SRCZONE[] ? S[i, :] : nothing, G = geom, bg = bgl) :
                             transfer_line_fast!(W, r, u, vth, κl, ν0, xs, lam == :diag, S[i, :], geom, bgl)
                for k in 1:n+1; Lline[i, k] = 4π * r[k]^2 * Fν[k]; end
                Jfs[i, :] .= Jz; Lam[i, :] .= Λz
                W === nothing || put!(wpool, W)
                if bgd !== nothing
                    for k in 1:n; Jfs[i, k] += bgd.Jbg[i, k]; end
                end
                Threads.atomic_add!(nth_t, 1)
            end
            nbmax = max(nbmax, nb_t[]); nthick = nth_t[]
            @. β = max(1 - Lam, 1e-8); @. Jeff = max(Jfs - Lam * S, 0.0)
            Pn = solve!(copy(P)); iters = it
            d = popchange(P, Pn); push!(hist, d)
            if d < ntol; converged = true; break; end
            P .= Pn
        end
        # net radiative rate of the last iterate with the formal-solution J̄ (thin lines: J̄ = 0)
        for k in 1:n
            s = 0.0
            for i in 1:nl
                pu = P[iu[i], k]; pl = P[il[i], k]
                s += h * ν[i] * (pu * LV.A[i] - (pl * Blu[i] - pu * Bul[i]) * Jfs[i, k])
            end
            Λrate[k] = nH2[k] * s
        end
    end
    Lr = [sum(view(Lline, :, k)) for k in 1:n+1]
    Λ = [(Lr[k+1] - Lr[k]) / (4π / 3 * (r[k+1]^3 - r[k]^3)) for k in 1:n]
    LAST[] = (nbins = nbmax, nthick = nthick, iters = iters, nep = nep, converged = converged, hist = hist, Λrate = Λrate,
              pops = P, nmaser = nmaser, mode = mode, Λem = Λem, kq = bgd === nothing ? 0 : bgd.kq)
    Λ, Lr, Lline
end


"Short-characteristics step with the source linear in τ (S_in → S_out across Δτ)."
@inline function sh(I, Δτ, Sin, Sout)
    if Δτ < 1e-4
        return I * (1 - Δτ) + 0.5Δτ * (Sin + Sout)
    end
    e = exp(-Δτ); w0 = 1 - e; w1 = w0 - Δτ * e
    I * e + Sout * (w0 - w1 / Δτ) + Sin * (w1 / Δτ)       # S(t) = S_out + (S_in − S_out) t/Δτ, t from the exit
end

"Frequency-integrated flux [per unit of the S·dν weights] at each interface from a formal solution on tangent rays without\nDoppler shifts: χ[k, f], S[k, f] per zone and frequency (S interpolated to the interfaces, linear in τ along segments);\neach frequency contributes with weight dν[f] = 1 (caller multiplies)."
function formal_flux(r, χ, S)
    n = size(χ, 1); nf = size(χ, 2); Fm = zeros(n + 1, nf)
    Threads.@threads for f in 1:nf
        Ff = zeros(n + 1); ΔI = zeros(n + 1, n)
        Si = zeros(n + 1); Si[1] = S[1, f]; Si[n+1] = S[n, f]          # source at the interfaces (linear in r)
        for k in 2:n
            rc0 = 0.5 * (r[k-1] + r[k]); rc1 = 0.5 * (r[k] + r[k+1]); w = (r[k] - rc0) / (rc1 - rc0)
            Si[k] = (1 - w) * S[k-1, f] + w * S[k, f]
        end
        for j in 1:n
            p = r[j]; nk = n + 2 - j; zk = [sqrt(max(r[k]^2 - p^2, 0.0)) for k in j:n+1]
            Iin = zeros(nk); I = 0.0
            for m in nk-1:-1:1
                kz = j + m - 1; I = sh(I, χ[kz, f] * (zk[m+1] - zk[m]), Si[j+m], Si[j+m-1]); Iin[m] = I
            end
            I = Iin[1]
            for (m, k) in enumerate(j:n+1)
                if m > 1; kz = j + m - 2; I = sh(I, χ[kz, f] * (zk[m] - zk[m-1]), Si[j+m-2], Si[j+m-1]); end
                ΔI[k, j] = m == 1 ? 0.0 : I - Iin[m]
            end
        end
        for k in 2:n+1
            acc = 0.0
            for j in 1:k-1
                μa = sqrt(max(1 - (r[j+1] / r[k])^2, 0.0)); μb = sqrt(max(1 - (r[j] / r[k])^2, 0.0))
                Ia = j + 1 <= n ? ΔI[k, j+1] : 0.0; Ib = ΔI[k, j]; d = μb - μa; d <= 0 && continue
                s1 = (Ib - Ia) / d
                acc += Ia * (μb^2 - μa^2) / 2 + s1 * ((μb^3 - μa^3) / 3 - μa * (μb^2 - μa^2) / 2)
            end
            Ff[k] = 2π * acc
        end
        Fm[:, f] .= Ff                  # one column per frequency (no shared writes between threads)
    end
    vec(sum(Fm; dims = 2))
end


end # module
