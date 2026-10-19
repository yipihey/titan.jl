# Population III project, M3c: where does each approximation hold during the collapse? Post-processing of snapshots
# (lib/popiii/snaps_*.jls, written by popiii_m2c.jl with SNAPS=1).
#   julia -t 16 --project=env popiii_regimes.jl SNAPFILE [every]
# Per snapshot and shell:
#   nLTE   n_H / n_crit(H₂), n_crit = Λ_LTE/Λ_0 (GP98/HM79): ≳ 1 → LTE level populations
#   chem   t_3b / t_ff, t_3b = 1/(k4 n_H²) (PSS83 three-body formation), and the equilibrium flag (ρ ≥ ρ_eq)
#   τline  max over the 231 H₂ lines of the radial line-centre depth from the shell to the surface (thermal profile)
#   continuum (ContinuumOpacity, 300 frequencies, 10–1e5 cm⁻¹), radial depths from the shell to the surface:
#          τP (Planck mean of absorption), τR (Rosseland mean of extinction), τRa (Rosseland of absorption);
#          fwin = fraction of ∫ ∂B/∂T dν carried at frequencies where the thermalisation depth
#                 τ_eff,ν = ∫ √(α_a (α_a + α_s)) dr < 1 (radiation there not thermalised: J_ν ≠ B_ν)
#   flux test: frequency-resolved formal solution (tangent rays, S = B_ν, scattering treated as absorption — upper
#          bracket; and absorption only — lower bracket) vs the Rosseland diffusion flux F_D = −(16σT³/3κ_R ρ) dT/dr at
#          each interface: ratio F/F_D (≈ 1 where diffusion with κ_R is justified).
using Printf, Serialization
include(joinpath(@__DIR__, "PrimordialChem.jl")); include(joinpath(@__DIR__, "ContinuumOpacity.jl")); include(joinpath(@__DIR__, "H2Lines.jl"))
const PC = PrimordialChem; const CO = ContinuumOpacity; const HL = H2Lines
const σSB = 5.670374e-5

snaps = deserialize(ARGS[1]); every = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 1
L = HL.load()
νs = [CO.c * 10^x for x in range(1.0, 5.0, length = 301)]
νm = [0.5 * (νs[i] + νs[i+1]) for i in 1:300]; dν = diff(νs)

const LiH = 4.15e-10          # ⁷Li/H, SBBN (Coc et al. 2004) as Mayer & Duschl 2005 Table 1
const USELI = get(ENV, "LI", "1") == "1"
"Saha ionisation of Li (5.3917 eV; statistical weights g(Li⁺) g(e)/g(Li) = 1·2/2): free electrons from Li at given n_e."
KLi(T) = (2π * PC.me * PC.kB * T / PC.hP^2)^1.5 * exp(-5.3917PC.eV / (PC.kB * T))
"""
Composition of zone j: number densities. Electrons: the network's (or Saha's) x_e from H, plus Li ⇌ Li⁺ + e⁻ in Saha
(LI=1, default; Li provides the electrons below ≈ 3500 K, Mayer & Duschl 2005 §2.3). H⁻: collisional steady state =
Saha with the local electrons (associative detachment and its inverse are 20–1e5× faster than photodetachment in a full
Planck field above 1e16 cm⁻³, M3c), unless the network has more.
"""
function comp(ρ, T, y)
    nHn = ρ / PC.μH; nH = PC.xH(y) * nHn; nH2 = y[3] * nHn; nHe = PC.yHe * nHn; nHp = y[1] * nHn; ne = max(PC.xe(y), 0.0) * nHn
    if USELI                                                  # n_Li⁺ = n_Li K/(K + n_e), n_e = n_e,H + n_Li⁺ (fixed point)
        nLi = LiH * nHn; K = KLi(T); ne0 = ne
        for _ in 1:50; ne = ne0 + nLi * K / (K + ne + 1e-300); end
    end
    nHm = max(y[2] * nHn, nH * ne * (PC.hP^2 / (2π * PC.me * PC.kB * T))^1.5 / 4 * exp(0.754PC.eV / (PC.kB * T)))
    (; nHn, nH, nH2, nHe, nHp, ne, nHm)
end

flux(r, χ, S) = HL.formal_flux(r, χ, S .* reshape(dν, 1, :))      # S·dν: frequency integral

@printf("%s: %d snapshots\n", ARGS[1], length(snaps))
for (is, s) in enumerate(snaps)
    (is % every == 0 || is == length(snaps)) || continue
    n = length(s.T); r = s.r
    lnc = log10(s.ρ[1] / PC.mH)
    C = [comp(s.ρ[j], s.T[j], s.y[j]) for j in 1:n]
    αa = zeros(n, 300); αs = zeros(n, 300); Bm = zeros(n, 300)
    for j in 1:n, f in 1:300
        a, sc, _ = CO.alpha(νm[f], s.T[j]; nH = C[j].nH, nH2 = C[j].nH2, nHe = C[j].nHe, nHp = C[j].nHp, ne = C[j].ne, nHm = C[j].nHm)
        αa[j, f] = a; αs[j, f] = sc; Bm[j, f] = CO.Bν(νm[f], s.T[j])
    end
    dB = [CO.dBdT(νm[f], s.T[j]) * dν[f] for j in 1:n, f in 1:300]
    κR = [sum(dB[j, :]) / sum(dB[j, :] ./ max.(αa[j, :] .+ αs[j, :], 1e-300)) / s.ρ[j] for j in 1:n]
    κRa = [sum(dB[j, :]) / sum(dB[j, :] ./ max.(αa[j, :], 1e-300)) / s.ρ[j] for j in 1:n]
    κP = [sum(αa[j, :] .* Bm[j, :] .* dν) / sum(Bm[j, :] .* dν) / s.ρ[j] for j in 1:n]
    dr = diff(r)
    # radial depths from shell j (its centre) to the surface
    cum(v) = (c = zeros(n); acc = 0.0; for j in n:-1:1; c[j] = acc + 0.5v[j] * dr[j]; acc += v[j] * dr[j]; end; c)
    τP = cum(κP .* s.ρ); τR = cum(κR .* s.ρ); τRa = cum(κRa .* s.ρ)
    τeff = zeros(n, 300)
    for f in 1:300; τeff[:, f] = cum(sqrt.(αa[:, f] .* (αa[:, f] .+ αs[:, f]))); end
    fwin = [sum(dB[j, f] for f in 1:300 if τeff[j, f] < 1; init = 0.0) / sum(dB[j, :]) for j in 1:n]
    # line-centre depth (max over lines)
    τline = zeros(n)
    pops = [HL.populations(L, s.T[j]) for j in 1:n]
    for i in eachindex(L.A)
        ν0 = L.ν[i]
        κc = [HL.h * ν0 / (4π) * HL.c^2 / (2HL.h * ν0^3) * L.A[i] * C[j].nH2 *
              (L.g[L.iu[i]] / L.g[L.il[i]] * pops[j][L.il[i]] - pops[j][L.iu[i]]) / (sqrt(π) * sqrt(2HL.kB * s.T[j] / HL.mH2) / HL.c * ν0) for j in 1:n]
        τline .= max.(τline, cum(κc))
    end
    nLTE = [C[j].nH * PC.Λ_H2_0_GP(s.T[j]) / PC.Λ_H2_LTE_HM(s.T[j]) for j in 1:n]
    tff = [sqrt(3π / (32PC.G * s.ρ[j])) for j in 1:n]; t3b = [1 / (PC.k4_pss(s.T[j]) * C[j].nH^2 + 1e-300) for j in 1:n]
    # flux test (interfaces 2..n): two brackets
    Fsa = flux(r, αa .+ αs, Bm); Fa = flux(r, αa, Bm)
    FD = zeros(n + 1)
    for k in 2:n
        κRi = 0.5 * (κR[k-1] + κR[k]); ρi = 0.5 * (s.ρ[k-1] + s.ρ[k]); Ti = 0.5 * (s.T[k-1] + s.T[k])
        dTdr = (s.T[k] - s.T[k-1]) / (0.5 * (r[k+1] - r[k-1]))
        FD[k] = -16σSB * Ti^3 / (3κRi * ρi) * dTdr
    end
    @printf("\n=== snapshot %d: step %d, t = %.6e yr, log n_c = %.2f, T_c = %.0f K\n", is, s.step, s.t / 3.156e7, lnc, s.T[1])
    @printf("%10s %9s %8s %8s %6s %9s %9s %9s %9s %9s %6s %10s %10s\n", "m [Msun]", "log n_H", "T", "n/ncrit", "eq",
            "t3b/tff", "tau_line", "tau_P", "tau_R", "tau_Rabs", "f_win", "F/F_D(s+a)", "F/F_D(a)")
    for mt in (1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0, 1000.0)
        j = findfirst(x -> x >= mt * 1.98847e33, s.m[2:end]); j === nothing && continue
        k = j + 1
        @printf("%10.3g %9.2f %8.0f %8.2g %6s %9.2g %9.2g %9.2g %9.2g %9.2g %6.3f %10.3g %10.3g\n", s.m[j+1] / 1.98847e33,
                log10(C[j].nHn), s.T[j], nLTE[j], s.eq[j] ? "yes" : "no", t3b[j] / tff[j], τline[j], τP[j], τR[j], τRa[j], fwin[j],
                FD[k] != 0 ? Fsa[k] / FD[k] : NaN, FD[k] != 0 ? Fa[k] / FD[k] : NaN)
    end
end
