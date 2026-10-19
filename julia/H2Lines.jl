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
characteristics):  I_out = I_in e^{−Δτ} + S_in (w0 − w1/Δτ) + S_out (w1/Δτ),  w0 = 1 − e^{−Δτ},  w1 = w0 − Δτ e^{−Δτ}
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

"""
    cooling(L, r, u, ρ, T, nH2; nbin = 30, vmax_extra = 0.0) -> (Λ, Lr)

Line cooling rate Λ_k [erg cm⁻³ s⁻¹] of each zone and the line luminosity L(r_i) [erg s⁻¹] at each interface.
r: interfaces (n+1, r[1] = 0 allowed), u: interface velocities (n+1), ρ, T, nH2: zones (n). Zone values are uniform
within the zone; the velocity along a ray inside a zone is interpolated linearly in r between the interfaces.
Frequency bins: per line, nbin bins spanning ±(3 v_th,max + max|u|) around the rest frequency (comoving shift taken
into account per segment).
"""
function cooling(L::LineList, r, u, ρ, T, nH2; nbin = 30)
    n = length(T); @assert length(r) == n + 1
    vth = [sqrt(2kB * T[k] / mH2) for k in 1:n]
    vspan = 3maximum(vth) + maximum(abs, u)
    xs = range(-vspan, vspan, length = nbin)          # velocity offsets of the bin centres [cm/s]
    dv = step(xs)
    # LTE populations per zone
    pops = [populations(L, T[k]) for k in 1:n]
    nl = length(L.A)
    Lr = zeros(n + 1)
    # rays: impact parameters p_j = r[j], j = 1..n (r[1] = 0: the central ray); intersections with interfaces k ≥ j
    Lline = zeros(nl, n + 1)
    Threads.@threads for i in 1:nl
        ν0 = L.ν[i]; gu = L.g[L.iu[i]]; gl = L.g[L.il[i]]
        # line-integrated opacity coefficient per zone: κ0_k = (hν/4π)(c²/2hν³) A (g_u/g_l n_l − n_u)   [cm⁻¹ × Hz]
        κ0 = [h * ν0 / (4π) * c^2 / (2h * ν0^3) * L.A[i] * nH2[k] * (gu / gl * pops[k][L.il[i]] - pops[k][L.iu[i]]) for k in 1:n]
        S = [Bν(ν0, T[k]) for k in 1:n]
        Fν = zeros(n + 1)                          # ∫ flux over the line (erg cm⁻² s⁻¹) at each interface
        for ib in 1:nbin
            vb = xs[ib]                            # observer-frame offset: ν = ν0 (1 + vb/c)
            # I along each ray at each interface crossing, outgoing (μ > 0) half; the incoming half from the far side
            # Ray j (p = r[j]): crosses interfaces k = j..n+1 at z_k = √(r_k² − p²); segments between consecutive k.
            for j in 1:n
                p = r[j]
                nk = n + 1 - j + 1                 # interfaces j..n+1
                zk = [sqrt(max(r[k]^2 - p^2, 0.0)) for k in j:n+1]
                # inward half (from z = +R down to z = 0): direction of propagation −ẑ; then outward half 0 → +R
                # segment m between interfaces j+m-1 and j+m lies in zone j+m-1
                function dτS(m, dir)               # optical depth and source of segment m (zone kz) along direction dir
                    kz = j + m - 1
                    zm = 0.5 * (zk[m] + zk[m+1]); rm = sqrt(p^2 + zm^2)
                    w = clamp((rm - r[kz]) / (r[kz+1] - r[kz]), 0.0, 1.0)
                    ur = (1 - w) * u[kz] + w * u[kz+1]
                    μ = rm > 0 ? dir * zm / rm : 0.0  # cosine between the ray and the radial direction
                    vcom = vb - μ * ur              # comoving-frame offset
                    φ = exp(-(vcom / vth[kz])^2) / (sqrt(π) * vth[kz] / c * ν0)    # profile per Hz
                    (κ0[kz] * φ * (zk[m+1] - zk[m]), S[kz])
                end
                # incoming: I = 0 at the outer boundary; integrate segments m = nk-1 down to 1 (towards the centre)
                Iin = zeros(nk); I = 0.0
                for m in nk-1:-1:1
                    Δτ, Sm = dτS(m, -1.0)
                    I = sh(I, Δτ, Sm, Sm)
                    Iin[m] = I
                end
                # outgoing: from the centre crossing (z = 0 at interface j) outward; at z = 0 the incoming I continues
                Iout = zeros(nk); I = Iin[1]; Iout[1] = I
                for m in 1:nk-1
                    Δτ, Sm = dτS(m, 1.0)
                    I = sh(I, Δτ, Sm, Sm)
                    Iout[m+1] = I
                end
                # contribution to the flux at each interface k ≥ j: weight from μ-quadrature (trapezoid in μ over rays)
                for (m, k) in enumerate(j:n+1)
                    k == 1 && continue
                    μk = zk[m] / r[k]
                    # trapezoid weights in μ between neighbouring rays j−1, j, j+1 at radius r[k]
                    μlo = j < k ? sqrt(max(1 - (r[j+1] / r[k])^2, 0.0)) : 0.0
                    μhi = j > 1 ? sqrt(max(1 - (r[j-1] / r[k])^2, 0.0)) : 1.0
                    wμ = 0.5 * (μhi - μlo)
                    Ip = Iout[m]; Im = m == 1 ? Iout[1] : Iin[m]
                    Fν[k] += 2π * wμ * μk * (Ip - Im) * (dv / c * ν0)
                end
            end
        end
        for k in 1:n+1; Lline[i, k] = 4π * r[k]^2 * Fν[k]; end
    end
    for k in 1:n+1; Lr[k] = sum(view(Lline, :, k)); end
    Λ = [(Lr[k+1] - Lr[k]) / (4π / 3 * (r[k+1]^3 - r[k]^3)) for k in 1:n]
    Λ, Lr, Lline
end

"Short-characteristics step with the source linear in τ (S_in → S_out across Δτ)."
@inline function sh(I, Δτ, Sin, Sout)
    if Δτ < 1e-4
        return I * (1 - Δτ) + 0.5Δτ * (Sin + Sout)
    end
    e = exp(-Δτ); w0 = 1 - e; w1 = w0 - Δτ * e
    I * e + Sin * (w0 - w1 / Δτ) + Sout * (w1 / Δτ)
end

end # module
