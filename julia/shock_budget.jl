# Energy budget of the accretion shock and line estimates for a TITAN snapshot [estimates; see report ch. 8].
#   Z0_TABLE=ext TITAN_NOLIB=1 julia --project=env shock_budget.jl RUN TIME[s|end]
# 1. Kinetic power arriving at the shock, ½Ṁv², against the cooling of the hot layer (from the shock down to where
#    T < 1e4 K): Compton (Λ_C as in TitanResidFast), free–free (1.4e-27 √T n_e Σ Z²n_i g, g = 1.2, Rybicki & Lightman
#    1979) and TITAN's Planck-mean exchange 4πκ_Pρ(B − cE/4π); X-ray luminosity of the free–free emission above 0.5
#    and 2 keV (thermal spectrum ∝ exp(−hν/kT)).
# 2. Energy-limited bounds on hydrogen line emission powered by the accretion (case B, Osterbrock 1989).
# 3. Sobolev optical depth of Hα in the free-falling gas above the shock, with LTE n = 2 populations:
#    τ_S = (πe²/m_e c) f λ n₂ / |dv/dr|, dv/dr = v/(2r) for free fall.
include(joinpath(@__DIR__, "Z0Physics.jl"))
using Serialization, Printf
const Z0 = Z0Physics
run = ARGS[1]; tq = ARGS[2]
S = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(run)_relax_snaps.jls"))
s = tq == "end" ? S[end] : S[argmin([abs(x.t - parse(Float64, tq)) for x in S])]
R0 = 8.4352e15; L0 = S[1].Lrad[end-1]; tab = Z0.table()
kB, me, c, σT, arad, h = Z0.kB, Z0.me, Z0.c, 6.6524616e-25, Z0.arad, Z0.h
vol(k) = 4π / 3 * (s.r[k+1]^3 - s.r[k]^3)
nz = length(s.d)
du = [s.u[k+1] - s.u[k] for k in 1:length(s.u)-1]; ks = argmin(du)               # shock: interface ks+1
kc = findlast(k -> k < ks && s.T[k] < 1e4, 1:nz)                               # bottom of the hot layer
# photosphere (Rosseland τ = 2/3 from the outer boundary, as in qstar_summary)
τ = 0.0; kph = 1
for k in nz:-1:1
    global τ += Z0.kap(tab, s.d[k], s.T[k]) * s.d[k] * (s.r[k+1] - s.r[k]); τ >= 2 / 3 && (global kph = k; break)
end
j = min(ks + 3, nz); vpre = abs(s.u[j]); Mdot = -4π * s.r[j]^2 * s.d[j] * s.u[j]
Pkin = 0.5 * Mdot * vpre^2
LC = 0.0; Lff = 0.0; Lex = 0.0; LX05 = 0.0; LX2 = 0.0
for k in kc+1:ks
    ρ, T, E = s.d[k], s.T[k], s.er[k]; st = Z0.saha(ρ, T); V = vol(k)
    Tr = (E / arad)^0.25
    global LC += 4σT * st.ne * c * E * kB * (T - Tr) / (me * c^2) * V
    εff = 1.4e-27 * sqrt(T) * st.ne * (st.nHp + st.nHe1 + 4st.nHe2) * 1.2
    global Lff += εff * V
    global LX05 += εff * exp(-0.5 * 1.602e-9 / (kB * T)) * V
    global LX2 += εff * exp(-2.0 * 1.602e-9 / (kB * T)) * V
    κP = Z0.kapP(tab, ρ, T)
    global Lex += κP * ρ * (4 * Z0.σSB * T^4 - c * E) * V
end
@printf("%s t = %.4e s: shock at r/R = %.5f (zone %d), hot layer down to r/R = %.5f (zone %d), photosphere r/R = %.5f\n",
        run, s.t, s.r[ks+1] / R0, ks, s.r[kc+1] / R0, kc, s.r[kph+1] / R0)
@printf("  Ṁ = %.3f Msun/yr, v_pre = %.3e cm/s, ½Ṁv² = %.3e erg/s = %.4f L0\n", Mdot * Z0.yr / Z0.Msun, vpre, Pkin, Pkin / L0)
@printf("  hot layer: Compton %.3e (%.3f of ½Ṁv²), free–free %.3e (%.3f), TITAN Planck-mean exchange %.3e (%.3f)\n",
        LC, LC / Pkin, Lff, Lff / Pkin, Lex, Lex / Pkin)
@printf("  free–free above 0.5 keV %.3e erg/s, above 2 keV %.3e erg/s (%.2e and %.2e L0)\n", LX05, LX2, LX05 / L0, LX2 / L0)
# Line emission powered by the accretion: case B cannot be applied to the gas near the photosphere (there the stellar
# radiation field keeps hydrogen ionised and the lines are set by the source function, not by recombination counting;
# applied blindly it gives L(Hα) larger than the accretion power). Energy-limited bounds instead:
#  (a) every accreted H atom recombines once after the shock: L(Lyα) ≤ (Ṁ X/m_H)·10.2 eV, L(Hα) ≤ (Ṁ X/m_H)·0.45·1.89 eV
#      (case B: 0.45 Hα photons per recombination at 1e4 K, Osterbrock 1989);
#  (b) X-rays absorbed in the infall and reprocessed: at most L_X(> 0.5 keV) in total.
NH = Mdot * Z0.COMP.X / Z0.mH
@printf("  energy bounds: L(Lyα) ≤ %.2e, L(Hα) ≤ %.2e erg/s from one recombination per accreted H; reprocessed X-rays ≤ %.2e erg/s\n",
        NH * 10.2 * 1.602e-12, NH * 0.45 * 1.89 * 1.602e-12, LX05)
# Hα Sobolev optical depth in the infall above the shock (grid part), LTE n = 2
fHa = 0.6407; λHa = 6562.8e-8; πe2 = 0.02654      # πe²/(m_e c) [cm² s⁻¹]
println("  Hα Sobolev optical depth in the infall (LTE n = 2):")
for k in round.(Int, range(ks + 2, nz, length = 4))
    ρ, T = s.d[k], s.T[k]; st = Z0.saha(ρ, T); pH, _ = Z0.level_pops(T, st)
    r = 0.5 * (s.r[k] + s.r[k+1]); v = abs(0.5 * (s.u[k] + s.u[k+1])); dvdr = v / (2r)
    τS = πe2 * fHa * λHa * pH[2] / dvdr
    @printf("    r/R = %.4f  T = %.0f K  ρ = %.2e  n₂ = %.2e cm⁻³  v = %.0f km/s  τ_S(Hα) = %.2e\n", r / R0, T, ρ, pH[2], v / 1e5, τS)
end
println("  … continued beyond the grid as free fall (T ∝ r^-1/2 from the grid's outer T, floor 1000 K; as Spectrum.jl):")
let r0 = s.r[end], T0 = s.T[end], M = 1e6 * Z0.Msun
    for f in (1.5, 2.0, 3.0, 5.0, 10.0)
        r = r0 * f; vff = sqrt(2Z0.G * M / r); ρ = Mdot / (4π * r^2 * vff); T = max(1000.0, T0 * sqrt(r0 / r))
        st = Z0.saha(ρ, T); pH, _ = Z0.level_pops(T, st); τS = πe2 * fHa * λHa * pH[2] / (vff / (2r))
        @printf("    r/R = %6.2f  T = %.0f K  ρ = %.2e  n₂ = %.2e  v = %.0f km/s  τ_S(Hα) = %.2e\n", r / R0, T, ρ, pH[2], vff / 1e5, τS)
    end
end
# 4. The infall column outside the shock: hydrogen column, Thomson depth, LTE neutral fraction, and a photon-counting
#    test of whether the shock's X-rays could keep it ionised (decides between the LTE and fully ionised limits used for
#    the spectra). Recombinations if fully ionised: ∫ α_B n_e n_p dV, α_B = 2.59e-13 (T/1e4 K)^-0.7 cm³ s⁻¹ (Osterbrock
#    1989, power-law fit, an estimate below 1e4 K). Ionisations available: L_X(> 0.5 keV)/(36 eV) if every X-ray's
#    energy goes into ionising neutral H (36 eV per ion pair, the high-ionisation-yield limit), an upper bound.
println("  infall column outside the shock (grid, then the free-fall extension to 100 r_out):")
let r0 = s.r[end], T0 = s.T[end], M = 1e6 * Z0.Msun
    NHc = 0.0; NH0 = 0.0; τes = 0.0; Rrec = 0.0
    seg(r1, r2, ρ, T) = begin
        st = Z0.saha(ρ, T); nH = st.nH; dr = r2 - r1
        NHc += nH * dr; NH0 += st.nH0 * dr; τes += st.ne * σT * dr
        αB = 2.59e-13 * (T / 1e4)^-0.7
        Rrec += αB * (nH + 2st.nHe) * nH * 4π * (0.5 * (r1 + r2))^2 * dr     # n_e(full) · n_p
    end
    for k in ks+1:nz; seg(s.r[k], s.r[k+1], s.d[k], s.T[k]); end
    rr = exp.(range(log(r0), log(100r0), length = 400))
    for i in 1:length(rr)-1
        r = sqrt(rr[i] * rr[i+1]); vff = sqrt(2Z0.G * M / r)
        seg(rr[i], rr[i+1], Mdot / (4π * r^2 * vff), max(1000.0, T0 * sqrt(r0 / r)))
    end
    QX = LX05 / (36 * 1.602e-12)
    @printf("    N_H = %.2e cm⁻², neutral (LTE) N_H0 = %.2e cm⁻² (%.3f), Thomson depth (LTE n_e) %.3f\n", NHc, NH0, NH0 / NHc, τes)
    @printf("    recombinations if fully ionised %.2e s⁻¹; ionisations the X-rays can supply ≤ %.2e s⁻¹ (ratio %.2e)\n", Rrec, QX, QX / Rrec)
end
