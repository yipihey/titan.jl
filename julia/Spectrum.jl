"""
    Spectrum

Frequency-dependent emergent spectrum of a spherical TITAN snapshot (post-processing; static, LTE).

* Opacity: the Z = 0 continuum of `z0/opacity.jl` at each zone (H I bound-free n ≤ 12, He I ground, He II n ≤ 8,
  free–free on H⁺/He⁺/He⁺⁺, H⁻ bound-free and free–free, Thomson and H Rayleigh scattering), LTE level populations
  from the zone's (ρ, T) via Saha. No lines, no gaunt factors ≠ 1, no molecules. The ground-state bound-free cross
  sections of H I, He I and He II are replaced by the Verner et al. (1996) fits (`coeffs`), which matter in X-rays.
  Zones flagged in `atm.ion` are taken as fully ionised (limiting cases).
* Source: S_ν = ε B_ν(T) + (1 − ε) J_ν, ε = κ/(κ + σ): LTE thermal emission plus coherent isotropic scattering.
  Thomson scattering is treated as coherent (no Compton energy exchange in the spectrum).
* Transfer: impact-parameter rays (one through every zone centre plus core rays), Feautrier second-order form,
  accelerated Λ-iteration with the exact diagonal of each ray's operator (Rybicki & Hummer 1991) and Ng acceleration;
  where that has not converged after `maxit` iterations (scattering-dominated hard X-rays), a direct solution with the
  full Λ matrix (`direct_J`). Emergent intensities from a formal solution with S linear in τ (`emergent`).
* Grid: zone centres, plus a pair of points either side of an interface where T or ρ jumps by more than ×1.5.
* Velocities are ignored (static frame): the infall at v/c ≤ 0.006 shifts continuum features by less than their width.
* Inner boundary: rays hitting the innermost zone see I⁺ = B_ν there (deep, thermalised).
* Beyond the TITAN grid the infall can be continued as free fall at the same Ṁ out to `rext`:
  ρ = Ṁ/(4πr²v_ff), T ∝ r^{-1/2} from the grid's outer temperature (dilute radiative equilibrium) with a floor.

`spectrum(atm, ν)` returns (L_ν [erg/s/Hz], ALI iteration counts[, frequency-integrated moments]).
`atmosphere(snapshot; ...)` builds `atm`; `lucy` is the Unsöld–Lucy temperature correction.
"""
module Spectrum

import ..Z0Physics
using LinearAlgebra
const Z0 = Z0Physics
const h = Z0.h; const kB = Z0.kB; const c = Z0.c

Bν(ν, T) = 2h * ν^3 / c^2 / expm1(h * ν / (kB * T))

"""
Optical depth of the segment between two neighbouring points. Normally the trapezoid; across an interface pair
(`jump`, the two points either side of a zone interface where T or ρ jumps) the smaller χ, so that the source
function, interpolated linearly between the two sides, does not collect optical depth from the opaque side (the
segment is 2 % of a zone wide, and a ramp from B(7000 K) to B(10⁶ K) over even that much of a dense zone's χ
would radiate more in the EUV than the hot zone itself).
"""
dtau(χa, χb, Δz, jump) = (jump ? min(χa, χb) : 0.5 * (χa + χb)) * Δz

"Zone-centred atmosphere from a snapshot (r interfaces, d, T zones). `ext = (Mdot [g/s], M [g], rmax/r_out, Tfloor)`."
function atmosphere(s; ext = nothing, nextra = 60, rmin = nothing, split = 1.5)
    rc = [0.5 * (s.r[k] + s.r[k+1]) for k in eachindex(s.d)]
    keep = rmin === nothing ? collect(eachindex(rc)) : findall(>=(rmin), rc)
    # Points: zone centres. Where T or ρ jumps by more than `split` between neighbouring zones, two more points are put
    # just either side of the interface, carrying the two zones' values: the solvers interpolate χ and S linearly
    # between points, and across a jump (e.g. 7000 K dense gas next to 1e6 K shocked gas) that would spread the hot
    # zone's source function over the cool zone's optical depth and create spurious emission (the flux inside the hot layer reached 4.2 L0, in the EUV).
    R = Float64[]; ρ = Float64[]; T = Float64[]; zone = Int[]; pair = Int[]
    for (i, k) in enumerate(keep)
        push!(R, rc[k]); push!(ρ, s.d[k]); push!(T, s.T[k]); push!(zone, k)
        if i < length(keep) && keep[i+1] == k + 1
            q = max(s.T[k+1] / s.T[k], s.T[k] / s.T[k+1], s.d[k+1] / s.d[k], s.d[k] / s.d[k+1])
            if q > split
                ri = s.r[k+1]; ε = 0.01 * min(ri - rc[k], rc[k+1] - ri)
                push!(R, ri - ε); push!(ρ, s.d[k]); push!(T, s.T[k]); push!(zone, k); push!(pair, length(R))
                push!(R, ri + ε); push!(ρ, s.d[k+1]); push!(T, s.T[k+1]); push!(zone, k + 1)
            end
        end
    end
    if ext !== nothing
        Mdot, M, fmax, Tfloor = ext
        r0 = s.r[end]; T0 = s.T[end]
        for x in exp.(range(log(r0 * 1.02), log(r0 * fmax), length = nextra))
            vff = sqrt(2 * Z0.G * M / x)
            push!(R, x); push!(ρ, Mdot / (4π * x^2 * vff)); push!(T, max(Tfloor, T0 * sqrt(r0 / x))); push!(zone, 0)
        end
    end
    jump = falses(length(R) - 1); jump[pair] .= true      # segment pair[i] → pair[i]+1 crosses an interface
    (r = R, ρ = ρ, T = T, zone = zone, jump = jump)
end

"Absorption and scattering coefficients (cm⁻¹) of every zone at frequency ν (precomputes Saha states once)."
struct Opac
    s::Vector{Any}; pH::Vector{Vector{Float64}}; pHe::Vector{Vector{Float64}}; T::Vector{Float64}
end
function Opac(atm)
    S = Any[]; PH = Vector{Float64}[]; PHe = Vector{Float64}[]
    for k in eachindex(atm.r)
        s = Z0.saha(atm.ρ[k], atm.T[k])
        if haskey(atm, :ion) && atm.ion[k]     # limiting case: H and He fully ionised (e.g. photoionised infall)
            s = merge(s, (ne = s.nH + 2s.nHe, nH0 = 0.0, nHp = s.nH, nHe0 = 0.0, nHe1 = 0.0, nHe2 = s.nHe, xH = 1.0))
        end
        pH, pHe = Z0.level_pops(atm.T[k], s)
        push!(S, s); push!(PH, pH); push!(PHe, pHe)
    end
    Opac(S, PH, PHe, copy(atm.T))
end
"""
Verner et al. (1996) ground-state photoionisation cross section (cm²) at photon energy E (eV):
σ = σ0 F(y), x = E/E0 − y0, y = √(x² + y1²), F = [(x−1)² + yw²] y^(P/2 − 5.5) (1 + √(y/ya))^(−P).
"""
function σverner(E, Eth, E0, σ0Mb, ya, P, yw, y0, y1)
    E < Eth && return 0.0
    x = E / E0 - y0; y = sqrt(x^2 + y1^2)
    1e-18 * σ0Mb * ((x - 1)^2 + yw^2) * y^(0.5P - 5.5) * (1 + sqrt(y / ya))^(-P)
end
σH1(E) = σverner(E, 13.6, 0.4298, 5.475e4, 32.88, 2.963, 0.0, 0.0, 0.0)
σHe1(E) = σverner(E, 24.59, 13.61, 949.2, 1.469, 3.188, 2.039, 0.4434, 2.136)
σHe2(E) = σH1(E / 4) / 4          # hydrogenic scaling σ_Z(E) = σ_H(E/Z²)/Z²

"""
Absorption and scattering coefficients (cm⁻¹). Same as Z0.kappa_nu (the opacities TITAN's tables are built from),
except that the ground-state bound-free cross sections of H I, He I and He II are replaced by the Verner et al. (1996)
fits. Z0.kappa_nu uses Kramers ν⁻³ (H I, He II) and 7.4e-18 (ν₀/ν)^1.7 (He I), which are adequate within a few times
threshold but too large in X-rays (He I by a factor of ~10³ at 10 keV), where the absorption by the cool infall matters.
"""
function coeffs(o::Opac, k, ν)
    ab, sc = Z0.kappa_nu(ν, o.T[k], o.s[k], o.pH[k], o.pHe[k])
    T = o.T[k]; s = o.s[k]; stim = -expm1(-Z0.h * ν / (Z0.kB * T))
    E = Z0.h * ν / Z0.eV
    old = 0.0; new = 0.0
    if ν >= Z0.χH / Z0.h; old += o.pH[k][1] * 2.815e29 / ν^3; new += o.pH[k][1] * σH1(E); end
    if ν >= Z0.χHe2 / Z0.h; old += o.pHe[k][1] * 2.815e29 * 16 / ν^3; new += o.pHe[k][1] * σHe2(E); end
    ν0 = Z0.χHe1 / Z0.h
    if ν >= ν0; old += s.nHe0 * 7.4e-18 * (ν0 / ν)^1.7; new += s.nHe0 * σHe1(E); end
    max(ab + (new - old) * stim, 1e-300), sc
end

"""
(absorption, scattering, thermal source function) at frequency ν for point k. For an LTE `Opac` the source is B_ν(T);
other opacity types (the non-LTE hydrogen of NLTEH.jl) add methods returning S_th = η_th/κ.
"""
coeffs3(o::Opac, k, ν) = (coeffs(o, k, ν)..., Bν(ν, o.T[k]))

"""
Rays: impact parameters p = 0, core rays inside r[1], and p = r[k] for every zone centre. For each ray the indices of
the zones it crosses (r ≥ p), outermost first, and the path coordinate z = √(r² − p²).
"""
function rays(r; ncore = 12)
    ps = vcat(r[1] .* sin.(range(0, π / 2, length = ncore + 1))[1:end-1], r)
    R = []
    for p in ps
        idx = findall(x -> x >= p * (1 - 1e-14), r)
        z = [sqrt(max(r[k]^2 - p^2, 0.0)) for k in idx]
        push!(R, (p = p, idx = reverse(idx), z = reverse(z), core = p < r[1]))
    end
    ps, R
end

"""
Feautrier solution along one ray for source S at the ray's points (outermost first) with extinction χ.
Returns u = (I⁺ + I⁻)/2 at the points, the diagonal of the Λ operator, and the emergent intensity I⁺ at the top.
Outer boundary: no incident radiation. Bottom: core ray → I⁺ = Ib; otherwise symmetry at z = 0.
"""
function feautrier(z, χ, S, core, Ib; jm = falses(length(z) - 1))
    n = length(z)
    if n == 1
        τ = 2χ[1] * z[1] + 1e-30
        return [S[1] * (1 - exp(-τ)) / 2], [(1 - exp(-τ)) / 2], S[1] * (1 - exp(-τ))
    end
    dτ = [dtau(χ[i], χ[i+1], abs(z[i] - z[i+1]), jm[i]) + 1e-300 for i in 1:n-1]
    a = zeros(n); b = zeros(n); cc = zeros(n); rhs = copy(S)
    # top (i = 1): u − I⁻=0 → first-order boundary with second-order correction (Auer 1967)
    b[1] = 1 + 2 / dτ[1] + 2 / dτ[1]^2; cc[1] = 2 / dτ[1]^2; rhs[1] = S[1]
    for i in 2:n-1
        d = 0.5 * (dτ[i-1] + dτ[i])
        a[i] = 1 / (dτ[i-1] * d); cc[i] = 1 / (dτ[i] * d); b[i] = 1 + a[i] + cc[i]
    end
    if core       # I⁺ = Ib at the bottom: u − (Ib − u)... diffusion-type: u = Ib − du/dτ
        a[n] = 2 / dτ[n-1]^2; b[n] = 1 + 2 / dτ[n-1] + 2 / dτ[n-1]^2; rhs[n] = S[n] + 2Ib / dτ[n-1]
    else          # symmetry at z = 0 (the last point is at z ≥ 0; for p = r[k] it is exactly 0)
        a[n] = 2 / dτ[n-1]^2; b[n] = 1 + 2 / dτ[n-1]^2
    end
    # Thomas algorithm with the Rybicki–Hummer diagonal of the inverse
    D = zeros(n); Z = zeros(n)
    D[1] = cc[1] / b[1]; Z[1] = rhs[1] / b[1]
    for i in 2:n
        den = b[i] - a[i] * D[i-1]
        D[i] = i < n ? cc[i] / den : 0.0
        Z[i] = (rhs[i] + a[i] * Z[i-1]) / den
    end
    u = zeros(n); u[n] = Z[n]
    for i in n-1:-1:1; u[i] = Z[i] + D[i] * u[i+1]; end
    E = zeros(n + 1); E[n] = a[n] / b[n]
    for i in n-1:-1:1; E[i] = a[i] / (b[i] - cc[i] * E[i+1]); end
    diag = zeros(n)
    for i in 1:n
        lower = i > 1 ? b[i] - a[i] * D[i-1] : b[1]
        diag[i] = 1 / (lower * (1 - D[i] * (i < n ? E[i+1] : 0.0)))
    end
    # emergent intensity: I⁻ = 0 at the top, so I⁺ = 2u
    return u, diag, 2u[1]
end

"""
Direct (non-iterative) solution of J = Λ[S], S = εB + (1−ε)J, for frequencies where the diagonal ALI converges too
slowly (scattering-dominated, ε ≪ 1, at large optical depth: hard X-rays inside the star and the hot layer). The
Feautrier system of every ray, T u = S + boundary term, is inverted explicitly; Λ[k′,k] = Σ_rays w_μ (T⁻¹)[k′,k] is
assembled and (I − Λ diag(1−ε)) J = Λ(εB) + J_b is solved by LU. Cost ∝ Σ_rays n², about 10× one ALI iteration
sweep per frequency here.
"""
function direct_J(atm, R, χ, ε, B, Ib)
    nr = length(atm.r); Λ = zeros(nr, nr); Jb = zeros(nr); wsum = zeros(nr)
    for ray in R
        idx = ray.idx; n = length(idx); z = ray.z; χr = χ[idx]
        if n == 1
            τ = 2χr[1] * z[1] + 1e-30; k = idx[1]; w = ray.wμ[1]
            Λ[k, k] += w * (1 - exp(-τ)) / 2; wsum[k] += w; continue
        end
        dτ = [dtau(χr[i], χr[i+1], abs(z[i] - z[i+1]), ray.jm[i]) + 1e-300 for i in 1:n-1]
        dl = zeros(n - 1); d = zeros(n); du = zeros(n - 1); bnd = zeros(n)
        d[1] = 1 + 2 / dτ[1] + 2 / dτ[1]^2; du[1] = -2 / dτ[1]^2
        for i in 2:n-1
            dd = 0.5 * (dτ[i-1] + dτ[i]); a = 1 / (dτ[i-1] * dd); c = 1 / (dτ[i] * dd)
            dl[i-1] = -a; du[i] = -c; d[i] = 1 + a + c
        end
        if ray.core
            dl[n-1] = -2 / dτ[n-1]^2; d[n] = 1 + 2 / dτ[n-1] + 2 / dτ[n-1]^2; bnd[n] = 2Ib / dτ[n-1]
        else
            dl[n-1] = -2 / dτ[n-1]^2; d[n] = 1 + 2 / dτ[n-1]^2
        end
        F = lu(Tridiagonal(dl, d, du))
        Tinv = F \ Matrix{Float64}(I, n, n)
        ub = F \ bnd
        for m2 in 1:n
            k2 = idx[m2]; w = ray.wμ[m2]; wsum[k2] += w; Jb[k2] += w * ub[m2]
            @inbounds for m in 1:n; Λ[k2, idx[m]] += w * Tinv[m2, m]; end
        end
    end
    Λ ./= wsum; Jb ./= wsum
    A = Matrix{Float64}(I, nr, nr) .- Λ .* (1 .- ε)'
    J = A \ (Λ * (ε .* B) .+ Jb)
    max.(J, 0.0)
end

"One segment of the formal solution with S linear in τ: intensity after crossing Δτ from S_a to S_b."
@inline function segment(I, Sa, Sb, dτ)
    dτ < 1e-6 && return I * (1 - dτ) + 0.5dτ * (Sa + Sb)
    e = exp(-dτ); f = (1 - e) / dτ
    I * e + Sa * (f - e) + Sb * (1 - f)
end
"""
Emergent intensity along a ray by the formal solution (S linear in τ between points): for a core ray from I = Ib at
the bottom outwards; otherwise inwards from I = 0 at the top to z = 0 and back out. Exact for any segment optical
depth, so it does not depend on resolving τ = 1 at the surface (the Feautrier u there does)."""
function emergent(z, χ, S, core, Ib; jm = falses(length(z) - 1))
    n = length(z)
    n == 1 && return S[1] * (1 - exp(-2χ[1] * z[1]))
    dτ = [dtau(χ[i], χ[i+1], abs(z[i] - z[i+1]), jm[i]) for i in 1:n-1]
    if core
        I = Ib
    else
        I = 0.0
        for i in 1:n-1; I = segment(I, S[i], S[i+1], dτ[i]); end      # inward to z = 0
    end
    for i in n-1:-1:1; I = segment(I, S[i+1], S[i], dτ[i]); end       # outward
    I
end

"""
    monochromatic(atm, op, rayset, ν; tol, maxit) -> (Lν, J, iterations)

ALI for one frequency. Angle weights: J(r_k) = Σ_rays w u with trapezoidal weights in μ = z/r at that radius.
"""
const NG_EVERY = 4
"""
Ng (1974) extrapolation of the last four iterates (oldest first), order 2, with weights 1/J² (Auer 1987). Returns
`nothing` if the 2×2 system is singular.
"""
function ng(h)
    x3, x2, x1, x0 = h[1], h[2], h[3], h[4]          # x0 newest
    w = 1 ./ (x0 .^ 2 .+ 1e-300)
    d0 = x0 .- x1; d1 = x1 .- x2; d2 = x2 .- x3
    a1 = d0 .- d1; a2 = d0 .- d2
    A11 = sum(w .* a1 .* a1); A12 = sum(w .* a1 .* a2); A22 = sum(w .* a2 .* a2)
    b1 = sum(w .* d0 .* a1); b2 = sum(w .* d0 .* a2)
    det = A11 * A22 - A12^2
    (abs(det) > 1e-14 * abs(A11 * A22) && isfinite(det)) || return nothing
    c1 = (b1 * A22 - b2 * A12) / det; c2 = (A11 * b2 - A12 * b1) / det
    (1 - c1 - c2) .* x0 .+ c1 .* x1 .+ c2 .* x2
end

function monochromatic(atm, op, R, ν; tol = 1e-5, maxit = 80, moments = false, maxit_direct = true)
    nr = length(atm.r)
    κ = zeros(nr); σ = zeros(nr); B = zeros(nr)
    for k in 1:nr
        ab, sc, Bk = coeffs3(op, k, ν); κ[k] = max(ab, 1e-300); σ[k] = sc; B[k] = Bk
    end
    χ = κ .+ σ; ε = κ ./ χ
    # μ weights per radius (trapezoid over the rays that reach it, μ from 1 (p = 0) down to 0 (p = r))
    J = copy(B); S = similar(B)
    Λd = zeros(nr)
    Ib = B[1]
    it = 0
    hist = Vector{Vector{Float64}}()
    for outer_it in 1:maxit
        it = outer_it
        S .= ε .* B .+ (1 .- ε) .* J
        Jn = zeros(nr); Λd .= 0.0; wsum = zeros(nr)
        for (j, ray) in enumerate(R)
            u, dg, _ = feautrier(ray.z, χ[ray.idx], S[ray.idx], ray.core, Ib; jm = ray.jm)
            for (m, k) in enumerate(ray.idx)
                μ = ray.z[m] / atm.r[k]
                # weight: trapezoid in μ between this ray and its neighbours (stored per ray below)
                w = ray.wμ[m]
                Jn[k] += w * u[m]; Λd[k] += w * dg[m]; wsum[k] += w
            end
        end
        Jn ./= wsum; Λd ./= wsum
        # ALI update: J_new = J_FS + Λ* (S_new − S_old) solved locally
        ΔJ = (Jn .- J) ./ (1 .- (1 .- ε) .* Λd)
        J .+= ΔJ
        # relative change, ignoring zones where J is below 1e-10 of its maximum (deep Wien tails, irrelevant to L_ν)
        maximum(abs.(ΔJ) ./ (abs.(J) .+ 1e-10 * maximum(J) .+ 1e-300)) < tol && break
        # Ng acceleration (order 2, Auer 1987) every NG_EVERY iterations once four iterates are stored
        push!(hist, copy(J)); length(hist) > 4 && popfirst!(hist)
        if length(hist) == 4 && outer_it % NG_EVERY == 0
            Jng = ng(hist)
            if Jng !== nothing && all(>(0), Jng)
                J .= Jng; empty!(hist)
            end
        end
    end
    if it == maxit && maxit_direct       # ALI did not converge: solve the scattering problem directly
        J = direct_J(atm, R, χ, ε, B, Ib); it = -1
    end
    S .= ε .* B .+ (1 .- ε) .* J
    Lν = 0.0
    ps = [ray.p for ray in R]
    Iout = [emergent(ray.z, χ[ray.idx], S[ray.idx], ray.core, Ib; jm = ray.jm) for ray in R]
    for j in 1:length(R)-1
        Lν += 0.5 * (Iout[j] * ps[j] + Iout[j+1] * ps[j+1]) * (ps[j+1] - ps[j])
    end
    Lν *= 8π^2       # L = 4π · 2π ∫ I p dp
    moments === :I && return Lν, J, it, (p = ps, I = Iout)
    if moments
        Jm, Hm = sweep_moments(atm, R, χ, S, Ib)
        return Lν, J, it, (J = Jm, H = Hm, κ = κ, χ = χ, B = B)
    end
    return Lν, J, it
end

"""
J and H (Eddington flux) at every radius from formal-solution sweeps along all rays with the converged source:
I⁻ inwards from the top (0 incident), I⁺ outwards (core rays from Ib, other rays from the reflected I⁻ at z = 0).
"""
function sweep_moments(atm, R, χ, S, Ib)
    nr = length(atm.r); J = zeros(nr); H = zeros(nr); wsum = zeros(nr)
    for ray in R
        idx = ray.idx; n = length(idx); z = ray.z
        Iin = zeros(n); Iout = zeros(n)
        dτ = [dtau(χ[idx[i]], χ[idx[i+1]], abs(z[i] - z[i+1]), ray.jm[i]) for i in 1:n-1]
        for i in 2:n; Iin[i] = segment(Iin[i-1], S[idx[i-1]], S[idx[i]], dτ[i-1]); end
        Iout[n] = ray.core ? Ib : Iin[n]
        for i in n-1:-1:1; Iout[i] = segment(Iout[i+1], S[idx[i+1]], S[idx[i]], dτ[i]); end
        for m in 1:n
            k = idx[m]; μ = z[m] / atm.r[k]; w = ray.wμ[m]
            J[k] += w * 0.5 * (Iin[m] + Iout[m]); H[k] += w * 0.5 * (Iout[m] - Iin[m]) * μ; wsum[k] += w
        end
    end
    J ./ wsum, H ./ wsum
end

"Attach trapezoidal μ weights to the rays (per radius, over the rays reaching it)."
function weights!(atm, R)
    nr = length(atm.r)
    lists = [Tuple{Int,Int,Float64}[] for _ in 1:nr]          # (ray, point, μ)
    for (j, ray) in enumerate(R), (m, k) in enumerate(ray.idx)
        push!(lists[k], (j, m, ray.z[m] / atm.r[k]))
    end
    W = [zeros(length(ray.idx)) for ray in R]
    for k in 1:nr
        L = sort(lists[k], by = x -> -x[3])                    # μ descending: 1 → 0
        μ = [x[3] for x in L]
        for (i, (j, m, _)) in enumerate(L)
            lo = i < length(L) ? μ[i] - μ[i+1] : 0.0
            hi = i > 1 ? μ[i-1] - μ[i] : 0.0
            W[j][m] = 0.5 * (lo + hi)
        end
    end
    jm(ray) = (n = length(ray.idx); haskey(atm, :jump) ? [atm.jump[ray.idx[i+1]] for i in 1:n-1] : falses(max(n - 1, 0)))
    [merge(ray, (wμ = W[j], jm = jm(ray))) for (j, ray) in enumerate(R)]
end

"Spectrum L_ν on frequency grid ν. With `moments = true` also the frequency-integrated moments for the temperature correction."
function spectrum(atm, ν; ncore = 12, tol = 1e-5, maxit = 80, verbose = false, moments = false)
    op = Opac(atm)
    _, R = rays(atm.r; ncore)
    R = weights!(atm, R)
    L = zeros(length(ν)); its = zeros(Int, length(ν))
    nr = length(atm.r)
    acc = [zeros(nr, 6) for _ in 1:Threads.nthreads()]          # ∫J, ∫H, ∫κJ, ∫κB, ∫B, ∫χH
    dν = [i == 1 ? ν[2] - ν[1] : i == length(ν) ? ν[end] - ν[end-1] : 0.5 * (ν[i+1] - ν[i-1]) for i in eachindex(ν)]
    Threads.@threads for i in eachindex(ν)
        if moments
            L[i], _, its[i], m = monochromatic(atm, op, R, ν[i]; tol, maxit, moments = true)
            a = acc[Threads.threadid()]; w = dν[i]
            @. a[:, 1] += w * m.J; @. a[:, 2] += w * m.H; @. a[:, 3] += w * m.κ * m.J
            @. a[:, 4] += w * m.κ * m.B; @. a[:, 5] += w * m.B; @. a[:, 6] += w * m.χ * m.H
        else
            L[i], _, its[i] = monochromatic(atm, op, R, ν[i]; tol, maxit)
        end
    end
    verbose && println("ALI iterations: median ", sort(its)[end÷2], ", max ", maximum(its), "; direct solutions: ", count(==(-1), its))
    moments || return L, its
    A = sum(acc)
    L, its, (J = A[:, 1], H = A[:, 2], κJ = A[:, 3] ./ A[:, 1], κB = A[:, 4] ./ A[:, 5], B = A[:, 5], χH = A[:, 6] ./ A[:, 2])
end

"""
Unsöld–Lucy temperature correction (plane-parallel form, Lucy 1964) on the zones `sel` (indices into atm, a
contiguous hydrostatic atmosphere, outermost last), towards radiative equilibrium with the target Eddington flux
H0(r) = L(r)/(16π²r²) (L may be a number or a per-zone vector, e.g. TITAN's radiative + convective luminosity):
  ΔB = (κ_J J − κ_B B)/κ_B + (κ_J/κ_B)·[3∫_r^{r_top} χ_H ΔH dr' + 2ΔH(r_top)],  ΔH = H0 − H.
Densities are held fixed. Returns the new temperatures and the largest relative flux error in `sel`.
"""
function lucy(atm, mom, sel, Ltarget; damp = 1.0, thin = Int[])
    σ = Z0.σSB; T = copy(atm.T)
    H0 = [(Ltarget isa Number ? Ltarget : Ltarget[k]) / (16π^2 * atm.r[k]^2) for k in eachindex(atm.r)]
    ΔH = H0 .- mom.H
    top = last(sel)
    integ = 0.0
    for k in reverse(sel)
        if k < top
            integ += 0.5 * (mom.χH[k] * ΔH[k] + mom.χH[k+1] * ΔH[k+1]) * (atm.r[k+1] - atm.r[k])
        end
        ΔB = (mom.κJ[k] * mom.J[k] - mom.κB[k] * mom.B[k]) / mom.κB[k] +
             mom.κJ[k] / mom.κB[k] * (3integ + 2ΔH[top])
        Bk = mom.B[k] + ΔB                          # frequency-integrated B = σT⁴/π
        Tn = (π * max(Bk, 0.25 * mom.B[k]) / σ)^0.25
        T[k] = clamp(atm.T[k] * (Tn / atm.T[k])^damp, 0.7 * atm.T[k], 1.3 * atm.T[k])
    end
    # points in `thin` (the infall): local radiative equilibrium only, ∫κ_ν J_ν = ∫κ_ν B_ν(T) (the first term above):
    # gas heated and cooled by radiation alone, without a flux condition (the infall is not a hydrostatic atmosphere)
    for k in thin
        Bk = mom.B[k] + (mom.κJ[k] * mom.J[k] - mom.κB[k] * mom.B[k]) / mom.κB[k]
        Tn = (π * max(Bk, 0.25 * mom.B[k]) / σ)^0.25
        T[k] = clamp(atm.T[k] * (Tn / atm.T[k])^damp, 0.7 * atm.T[k], 1.3 * atm.T[k])
    end
    T, isempty(sel) ? 0.0 : maximum(abs.(ΔH[sel] ./ H0[sel]))
end

"Frequency grid (Hz) from λ = λmax down to photon energy Emax [keV], with every opacity edge bracketed."
function grid(; λmax_μm = 10.0, Emax_keV = 40.0, n = 500, nedge = 6)
    ν1 = c / (λmax_μm * 1e-4); ν2 = Emax_keV * 1.602176634e-9 / h
    ν = collect(exp.(range(log(ν1), log(ν2), length = n)))
    edges = Float64[]
    for nl in 1:nedge; push!(edges, Z0.χH / (h * nl^2)); end
    for nl in 1:3; push!(edges, Z0.χHe2 / (h * nl^2)); end
    push!(edges, Z0.χHe1 / h); push!(edges, c / 1.6419e-4)
    for e in edges, f in (0.995, 1.005); ν1 < e * f < ν2 && push!(ν, e * f); end
    sort!(unique!(ν))
end

end # module
