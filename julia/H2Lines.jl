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
"Bins used and number of lines solved by transfer in the last call."
const LAST = Ref((nbins = 0, nthick = 0))

"""
    cooling(L, r, u, ρ, T, nH2; nbin = 30, vmax_extra = 0.0) -> (Λ, Lr)

Line cooling rate Λ_k [erg cm⁻³ s⁻¹] of each zone and the line luminosity L(r_i) [erg s⁻¹] at each interface.
r: interfaces (n+1, r[1] = 0 allowed), u: interface velocities (n+1), ρ, T, nH2: zones (n). Zone values are uniform
within the zone; the velocity along a ray inside a zone is interpolated linearly in r between the interfaces.
Frequency bins: per line, nbin bins spanning ±(3 v_th,max + max|u|) around the rest frequency (comoving shift taken
into account per segment).
"""
function cooling(L::LineList, r, u, ρ, T, nH2; nbin = 30, fnlte = nothing, τthin = 0.02, nbinmax = 200)
    n = length(T); @assert length(r) == n + 1
    vth = [sqrt(2kB * T[k] / mH2) for k in 1:n]
    # span and bin width from the molecular zones (fraction of H nuclei in H₂ ≥ 0.1; μ_H = 2.2316e-24 g, metal free with
    # n_He/n_H = 1/12): the hot, mostly dissociated core (f_H₂ ≈ 0.04 at T ~ 4e4 K, v_th ~ 17 km/s) would otherwise set the
    # span and push the bins past the cap (M3d). Its lines are buried under the opaque continuum anyway.
    act = [2nH2[k] * 2.2316e-24 / ρ[k] >= 0.1 for k in 1:n]
    any(act) || (act = trues(n))
    vspan = 3maximum(vth[act]) + maximum(abs, u)
    # bins no wider than half the narrowest thermal width (cold outer zones), at least nbin, at most nbinmax
    nb = clamp(ceil(Int, 2vspan / (0.5minimum(vth[act]))) + 1, nbin, nbinmax)
    xs = range(-vspan, vspan, length = nb)            # velocity offsets of the bin centres [cm/s]
    dv = step(xs)
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
    LAST[] = (nbins = nb, nthick = nthick[])
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
