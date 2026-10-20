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

function load(dir = joinpath(homedir(), "Projects", "popiii_lit", "H2lines"); vmax = 2, Jmax = 20)
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
                 nlte = nothing, nlte_mode = :ali, nlte_tol = 1e-3, nlte_maxit = 60, nlte_init = :ep, nlte_lambda = :diag)
    if nlte !== nothing
        return cooling_nlte(L, r, u, T, nH2, nlte; nbin, tol, nbinmax, cover, mode = nlte_mode, ntol = nlte_tol,
                            maxit = nlte_maxit, init = nlte_init, lam = nlte_lambda)
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
function transfer_line(r, u, vth, κ0, Si, ν0, xs, lam)
    n = length(vth); nb = length(xs); dv = step(xs)
    Fν = zeros(n + 1); Jz = zeros(n); Λz = zeros(n)
    ΔI = zeros(n + 1, n); QJ = zeros(n + 1, n, 2); QΛ = zeros(n + 1, n, 2)
    φv(v, k) = exp(-(v / vth[k])^2) / (sqrt(π) * vth[k])
    for ib in 1:nb
        vb = xs[ib]
        fill!(ΔI, 0.0)
        for j in 1:n
            p = r[j]; nk = n + 2 - j
            zk = [sqrt(max(r[k]^2 - p^2, 0.0)) for k in j:n+1]
            dti = zeros(nk - 1); dto = zeros(nk - 1)
            for m in 1:nk-1
                kz = j + m - 1
                zm = 0.5 * (zk[m] + zk[m+1]); rm = sqrt(p^2 + zm^2)
                w = clamp((rm - r[kz]) / (r[kz+1] - r[kz]), 0.0, 1.0)
                ur = (1 - w) * u[kz] + w * u[kz+1]
                μ = rm > 0 ? zm / rm : 0.0
                pre = κ0[kz] * (zk[m+1] - zk[m]) / (sqrt(π) * vth[kz] / c * ν0)
                dti[m] = pre * exp(-((vb + μ * ur) / vth[kz])^2)          # inward (dir = −1)
                dto[m] = pre * exp(-((vb - μ * ur) / vth[kz])^2)          # outward
            end
            Iin = zeros(nk); ψin = zeros(nk); Iout = zeros(nk); ψout = zeros(nk); I = 0.0
            for m in nk-1:-1:1
                I = sh(I, dti[m], Si[j+m], Si[j+m-1]); Iin[m] = I; ψin[m] = ψloc(dti[m])
            end
            Iout[1] = Iin[1]; ψout[1] = ψin[1]; I = Iin[1]
            for m in 2:nk
                I = sh(I, dto[m-1], Si[j+m-2], Si[j+m-1]); Iout[m] = I; ψout[m] = ψloc(dto[m-1])
            end
            for m in 1:nk
                k = j + m - 1; ΔI[k, j] = m == 1 ? 0.0 : Iout[m] - Iin[m]
                μk = r[k] > 0 ? zk[m] / r[k] : 0.0; uk = u[k]
                for (s, zs) in ((1, k - 1), (2, k))
                    (zs < 1 || zs > n) && continue
                    pin = φv(vb + μk * uk, zs); pout = φv(vb - μk * uk, zs)
                    QJ[k, j, s] = pin * Iin[m] + pout * Iout[m]
                    lam && (QΛ[k, j, s] = pin * ψin[m] + pout * ψout[m])
                end
            end
        end
        for k in 2:n+1
            acc = 0.0
            for j in 1:k-1
                μa = sqrt(max(1 - (r[j+1] / r[k])^2, 0.0)); μb = sqrt(max(1 - (r[j] / r[k])^2, 0.0))
                Ia = j + 1 <= n ? ΔI[k, j+1] : 0.0; Ib = ΔI[k, j]
                d = μb - μa; d <= 0 && continue
                s1 = (Ib - Ia) / d
                acc += Ia * (μb^2 - μa^2) / 2 + s1 * ((μb^3 - μa^3) / 3 - μa * (μb^2 - μa^2) / 2)
            end
            Fν[k] += 2π * acc * (dv / c * ν0)
        end
        # mean intensity at the nodes: J = ½ ∫₀¹ (Q_in + Q_out) dμ, trapezoid in μ between the rays (Q includes the profile)
        for k in 1:n+1, (s, zs) in ((1, k - 1), (2, k))
            (zs < 1 || zs > n) && continue
            if k == 1
                Jz[zs] += 0.25 * QJ[1, 1, s] * dv; lam && (Λz[zs] += 0.25 * QΛ[1, 1, s] * dv)
                continue
            end
            aJ = 0.0; aΛ = 0.0
            for j in 1:k-1
                μa = sqrt(max(1 - (r[j+1] / r[k])^2, 0.0)); μb = sqrt(max(1 - (r[j] / r[k])^2, 0.0))
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
    cooling_nlte(L, r, u, T, nH2, nlte; mode = :ali | :ep, ntol, maxit, init, lam)

See `cooling` (keyword `nlte = (LV, nH, nHe, nHp)`). Returns (Λ, Lr, Lline) as the LTE path (Λ = −∇·F per zone from the
flux of the formal solution of the last iterate); `LAST[]` gets `iters` (ALI iterations; the EP start is counted in `nep`),
`converged`, `hist` (max relative population change per iteration), `Λrate` (net radiative rate Σ hν[n_u A − (n_l B_lu −
n_u B_ul) J̄] per zone, ALI) and `pops` (63 × n).
"""
function cooling_nlte(L::LineList, r, u, T, nH2, nlte; nbin = 30, tol = 0.01, nbinmax = 200, cover = 0.999, mode = :ali,
                      ntol = 1e-3, maxit = 60, init = :ep, lam = :diag)
    LV = nlte.LV; HV = Main.H2Levels; check_levels(L, LV)
    n = length(T); @assert length(r) == n + 1
    nl = length(L.A); nlev = length(LV.E)
    vth = [sqrt(2kB * T[k] / mH2) for k in 1:n]
    τthin = 2tol
    iu = LV.iu; il = LV.il; Bul = LV.Bul; Blu = LV.Blu; ν = L.ν
    P = zeros(nlev, n); β = ones(nl, n); Jeff = zeros(nl, n)
    function solve!(Pn)
        Threads.@threads for k in 1:n
            Pn[:, k] .= HV.populations(LV, T[k], nlte.nH[k], nH2[k], nlte.nHe[k], nlte.nHp[k]; β = β[:, k], Jbar = Jeff[:, k])
        end
        Pn
    end
    # opacity and source from populations
    function opac(Pc)
        κ0 = zeros(nl, n); S = zeros(nl, n); nmaser = 0
        for k in 1:n, i in 1:nl
            pu = Pc[iu[i], k]; pl = Pc[il[i], k]; den = pl * Blu[i] - pu * Bul[i]
            if den > 0
                κ0[i, k] = h * ν[i] / (4π) * nH2[k] * den; S[i, k] = pu * LV.A[i] / den
            else
                nmaser += pu > 0 ? 1 : 0; S[i, k] = Bν(ν[i], T[k])
            end
        end
        κ0, S, nmaser
    end
    # escape probability per line and zone from the current per-line depths (EP mode): static part = profile-averaged
    # two-direction slab escape ½[E₂(τ_out φ̂) + E₂(τ_in φ̂)], φ̂ = e^{−x²}, with τ_out the line-centre depth along the radius to
    # the surface and τ_in that along the inward radius through the centre to the surface; Sobolev part (velocity gradient
    # g = max(|du/dr|, |u/r|)): (1 − e^{−τ_S})/τ_S, τ_S = κ0 c/(ν0 g); β = the larger of the two (a photon leaves by the
    # easier route). β = 1 where κ0 ≤ 0.
    function epβ(κ0)
        βn = ones(nl, n)
        xg = 0.0:0.25:5.0; wg = [(x == 0 || x == 5.0 ? 0.5 : 1.0) * 0.25 * 2 / sqrt(π) * exp(-x^2) for x in xg]
        for i in 1:nl
            τs = zeros(n); acc = 0.0
            for k in n:-1:1
                dτ = max(κ0[i, k], 0.0) / (sqrt(π) * vth[k] / c * ν[i]) * (r[k+1] - r[k])
                τs[k] = acc + 0.5dτ; acc += dτ
            end
            τtot = acc
            for k in 1:n
                κ0[i, k] > 0 || continue
                τo = τs[k]; τi = (τtot - τs[k]) + τtot
                b = 0.0
                for (x, w) in zip(xg, wg); φ̂ = exp(-x^2); b += w * 0.5 * (E2(τo * φ̂) + E2(τi * φ̂)); end
                rc = 0.5 * (r[k] + r[k+1]); g = max(abs(u[k+1] - u[k]) / (r[k+1] - r[k]), abs(0.5 * (u[k] + u[k+1])) / max(rc, 1e-300))
                if g > 0
                    τS = κ0[i, k] * c / (ν[i] * g); b = max(b, τS < 1e-6 ? 1 - τS / 2 : -expm1(-τS) / τS)
                end
                βn[i, k] = min(b, 1.0)
            end
        end
        βn
    end
    hist = Float64[]; nep = 0; iters = 0; converged = false
    solve!(P)
    if mode == :ep || init == :ep
        while nep < maxit
            κ0, S, _ = opac(P); β .= epβ(κ0); Pn = solve!(copy(P)); nep += 1
            d = popchange(P, Pn); P .= Pn; mode == :ep && push!(hist, d); d < ntol && break
        end
        mode == :ep && (converged = !isempty(hist) && hist[end] < ntol)
    end
    nbmax = 0; nthick = 0; Lline = zeros(nl, n + 1); Λrate = zeros(n); nmaser = 0
    if mode == :ep
        κ0, S, nmaser = opac(P)
        # local emission n_u A β hν accumulated over the volume
        for i in 1:nl, k in 1:n
            Lline[i, k+1] = Lline[i, k] + nH2[k] * P[iu[i], k] * LV.A[i] * β[i, k] * h * ν[i] * 4π / 3 * (r[k+1]^3 - r[k]^3)
        end
        iters = nep
    else
        Jfs = zeros(nl, n); Lam = zeros(nl, n); β .= 1; Jeff .= 0
        for it in 1:maxit
            κ0, S, nmaser = opac(P)
            fill!(Jfs, 0.0); fill!(Lam, 0.0); fill!(Lline, 0.0)
            nb_t = Threads.Atomic{Int}(0); nth_t = Threads.Atomic{Int}(0)
            Threads.@threads for i in 1:nl
                ν0 = ν[i]; κl = κ0[i, :]
                τc = sum(κl[k] / (sqrt(π) * vth[k] / c * ν0) * (r[k+1] - r[k]) for k in 1:n)
                if τc < τthin
                    acc = 0.0
                    for k in 1:n
                        acc += nH2[k] * P[iu[i], k] * LV.A[i] * h * ν0 * 4π / 3 * (r[k+1]^3 - r[k]^3); Lline[i, k+1] = acc
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
                Fν, Jz, Λz = transfer_line(r, u, vth, κl, Si, ν0, xs, lam == :diag)
                for k in 1:n+1; Lline[i, k] = 4π * r[k]^2 * Fν[k]; end
                Jfs[i, :] .= Jz; Lam[i, :] .= Λz
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
              pops = P, nmaser = nmaser, mode = mode)
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
