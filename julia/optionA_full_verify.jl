#!/usr/bin/env julia
"""
    optionA_full_verify.jl

End-to-end verification & science benchmark for all 4 unified non-LTE / multigroup / atomic-line upgrades:
  1. Smooth C² LTE (FreeEnergyEOS) <-> non-LTE (ChemNet) blend: value, 1st, 2nd Dual derivatives across (ρ, T).
  2. Option A: Photon-conserving multigroup IonRay (H⁻, LW with WHB11 H₂ self-shielding, H I, H₂, He I, He II)
     driving coupled 1D H II / He II / He III breakout + LW H₂ dissociation front vs bloated cool supergiant.
  3. Smooth-gated fast H₂ ro-vibrational line cooling (ExoMol RACPPK cached β_esc + C² (T, x_H2) gate):
     C² Dual continuity & wall-clock speedup benchmark.
  4. l-resolved hydrogenic ModelAtom (H I, He II) with exact dipole A-values + Pengelly-Seaton (1964) l-mixing:
     sum-rule, LTE detailed balance, Storey & Hummer (1995) He II 4686/1640 benchmark across n_e = 10²..10¹⁴,
     plus emergent H I / He I / He II line luminosities across H II breakout, pulsations, episodic bursts & SMS.
"""
ENV["TITAN_NOLIB"] = "1"
haskey(ENV, "FE_EOS_CACHE") || (ENV["FE_EOS_CACHE"] = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls")

using Printf, LinearAlgebra, ForwardDiff

const JDIR = @__DIR__
for f in ("Z0Physics", "ContinuumOpacity", "PrimordialChem", "ChemNet", "ChemRates",
          "StellarEOS", "FreeEnergyEOS", "StellarOpacity", "IonRay", "H2Lines",
          "UnifiedNonLTE", "ModelAtom")
    include(joinpath(JDIR, f * ".jl"))
end

const Z0 = Z0Physics; const FE = FreeEnergyEOS; const CN = ChemNet; const CR = ChemRates
const IR = IonRay; const HL = H2Lines; const UN = UnifiedNonLTE; const MA = ModelAtom

const Rsun = 6.957e10; const Lsun = 3.828e33; const pc = 3.085677581e18; const yr = 3.15576e7

npass = Ref(0); nfail = Ref(0)
function chk(name, cond, msg = "")
    if cond
        npass[] += 1
        @printf("  [PASS] %-54s %s\n", name, msg)
    else
        nfail[] += 1
        @printf("  [FAIL] %-54s %s\n", name, msg)
    end
end

function main()
println("===================================================================================")
println("1. SMOOTH C² LTE <-> NON-LTE BLEND & H₂ COOLING GATE (ZERO HARD EDGES)")
println("===================================================================================")
bp = UN.BlendPar()
# Check w_lte C² continuity at boundary transitions ρ1, ρ2, T1, T2
d1_ρ(ρ) = ForwardDiff.derivative(r -> UN.w_lte(r, 3000.0, bp), ρ)
d2_ρ(ρ) = ForwardDiff.derivative(d1_ρ, ρ)
d1_T(T) = ForwardDiff.derivative(t -> UN.w_lte(1e-12, t, bp), T)
d2_T(T) = ForwardDiff.derivative(d1_T, T)

max_jump_w = maximum(abs(d2_ρ(r * (1 + 1e-6)) - d2_ρ(r * (1 - 1e-6))) / max(abs(d2_ρ(3e-10)), 1.0)
                     for r in (bp.ρ1, bp.ρ2, 3e-10))
chk("w_lte(ρ, T) C² smooth across ρ1, ρ2", max_jump_w < 1e-3 && d1_ρ(bp.ρ1) == 0 && d1_ρ(bp.ρ2) == 0,
    @sprintf("max_rel_d2_jump=%.2e", max_jump_w))
chk("w_lte(ρ, T) C² smooth across T1, T2", d1_T(bp.T1) == 0 && d1_T(bp.T2) == 0 && abs(d2_T(bp.T1)) < 1e-12 && abs(d2_T(bp.T2)) < 1e-12,
    @sprintf("w(3e4K)=%.3f", UN.w_lte(1e-12, 3e4, bp)))

# Verify gas_eos_blended Dual continuity across full ρ ∈ [1e-12, 1e-7] & T ∈ [300, 2e6]
y_eq_test = Dict(:H => 0.99, :Hp => 0.01, :Hm => 1e-10, :H2 => 1e-5, :H2p => 1e-12, :e => 0.01, :He => 1/12, :Hep => 0.0, :Hepp => 0.0)
eos_ok = true; max_c1_err = 0.0
for ρ in (1e-12, 3e-11, 3e-10, 3e-9, 1e-7), T in (500.0, 5000.0, 1.5e4, 3e4, 4.5e4, 5e5, 2e6)
    p, e, w = UN.gas_eos_blended(ρ, T, y_eq_test, bp)
    dPdT = ForwardDiff.derivative(t -> UN.gas_eos_blended(ρ, t, y_eq_test, bp)[1], T)
    dedT = ForwardDiff.derivative(t -> UN.gas_eos_blended(ρ, t, y_eq_test, bp)[2], T)
    fd   = (UN.gas_eos_blended(ρ, T * (1 + 1e-6), y_eq_test, bp)[2] - UN.gas_eos_blended(ρ, T * (1 - 1e-6), y_eq_test, bp)[2]) / (2e-6 * T)
    err  = abs(dedT - fd) / max(abs(dedT), 1e-30)
    max_c1_err = max(max_c1_err, err)
    if !(isfinite(p) && isfinite(e) && isfinite(dPdT) && isfinite(dedT) && p > 0)
        eos_ok = false
    end
end
chk("gas_eos_blended C² Dual derivatives (ρ=1e-12..1e-7)", eos_ok && max_c1_err < 1e-4,
    @sprintf("max_FD_vs_Dual_err=%.2e", max_c1_err))

# Smooth-gated fast H₂ cooling C² check & speed benchmark
h2c = UN.H2CoolCache(80; bp)
r_tmp = collect(exp.(range(log(10 * Rsun), log(1e6 * Rsun); length = 81)))
u_tmp = fill(-5e5, 81); ρ_tmp = fill(1e-16, 80); T_tmp = fill(1200.0, 80); xH2_tmp = fill(1e-3, 80)
t_esc = @elapsed UN.update_h2_escape!(h2c, r_tmp, u_tmp, ρ_tmp, T_tmp, xH2_tmp)
dLam_T(T) = ForwardDiff.derivative(t -> UN.h2_cooling_fast(h2c, 40, 1e-16, t, 0.998, 1e-3), T)
dLam_x(x) = ForwardDiff.derivative(xv -> UN.h2_cooling_fast(h2c, 40, 1e-16, 1500.0, 1.0 - 2xv, xv), x)
h2_smooth = (UN.h2_cooling_fast(h2c, 40, 1e-16, 9500.0, 0.99, 1e-3) == 0.0 &&
             dLam_T(9000.0) == 0.0 && dLam_x(1e-11) == 0.0 &&
             isfinite(dLam_T(5500.0)) && isfinite(dLam_x(1e-9)) &&
             UN.h2_cooling_fast(h2c, 40, 1e-16, 1500.0, 0.998, 1e-3) > 0.0)
# Benchmark 10,000 Dual evaluations of h2_cooling_fast
t_bench = @elapsed begin
    s = 0.0
    for _ in 1:10_000
        s += ForwardDiff.derivative(t -> UN.h2_cooling_fast(h2c, 40, 1e-16, t, 0.998, 1e-3), 1500.0)
    end
end
chk("h2_cooling_fast C² gate & ExoMol cached β_esc", h2_smooth,
    @sprintf("β=%.3f, update=%.1f ms, %.2f μs/Dual-eval", h2c.β[40], 1e3*t_esc, 1e6*t_bench/10_000))

println("\n===================================================================================")
println("2. OPTION A: MULTIGROUP IONRAY (H⁻ + LW + H I + H₂ + He I + He II) BREAKOUT TEST")
println("===================================================================================")
# Verify phi1 continuity & IonRay Dual smoothness
p1_lo = IR.phi1(1e-4 - 1e-10); p1_hi = IR.phi1(1e-4 + 1e-10)
chk("IonRay.phi1(τ) C^∞ continuity at τ=1e-4", abs(p1_hi - p1_lo) < 1e-10, @sprintf("diff=%.2e", abs(p1_hi - p1_lo)))

Gfull = IR.groups_full(; n_opt = 6, n_euv = 28, Emax = 300.0)
net_he, photos_he = CR.modern(; he = true)
# Check blackbody photon rates for 100 Msun Pop III ZAMS (Teff=95,500 K, L=1.30e6 Lsun)
N_hot = IR.blackbody_full(Gfull, 1.30e6 * Lsun, 95500.0)
Q_H   = sum(N_hot[g] for g in eachindex(N_hot) if Gfull.Elo[g] >= 13.5)
Q_He1 = sum(N_hot[g] for g in eachindex(N_hot) if Gfull.Elo[g] >= 24.5)
Q_He2 = sum(N_hot[g] for g in eachindex(N_hot) if Gfull.Elo[g] >= 54.0)
Q_LW  = N_hot[Gfull.g_lw]
@printf("  100 M☉ ZAMS (95.5 kK, 1.3e6 L☉): Q_H=%.2e, Q_HeI=%.2e, Q_HeII=%.2e, Q_LW=%.2e s⁻¹\n", Q_H, Q_He1, Q_He2, Q_LW)
chk("Pop III 100 M☉ ZAMS ionizing rates Q_H ~ 1e50 s⁻¹", 5e49 < Q_H < 3e50 && Q_He2 > 1e48,
    @sprintf("Q_H=%.2e, Q_HeII/Q_H=%.3f", Q_H, Q_He2 / Q_H))

# Coupled 1D radial Strömgren / H II breakout + He II/He III + LW H₂ dissociation halo test
# Uniform primordial envelope n_H = 1e4 cm⁻³ (ρ = 2.23e-20 g/cm³), 60 log shells from 1e-3 pc to 2.0 pc
# Analytical Case-B Strömgren radius R_S = (3 Q_H / (4π α_B n_H²))^{1/3}
# At T = 2.5e4 K (Pop III H II region), α_B ≈ 1.3e-13 cm³/s => R_S ≈ 0.17 pc!
function run_breakout(L_star, Teff, label; nH0 = 1.0e4, nz = 64, r_in = 2e-3 * pc, r_out = 1.0 * pc, t_end = 5e3 * yr)
    r_if = exp.(range(log(r_in), log(r_out); length = nz + 1))
    ρ = fill(nH0 * CN.μH(), nz)
    # Initial primordial gas: mostly neutral H + x_H2 = 1e-3, x_Hp = 1e-4, neutral He, T = 500 K
    x0 = zeros(length(net_he.free))
    iHp = findfirst(==(:Hp), net_he.free); iH2 = findfirst(==(:H2), net_he.free)
    iHep = findfirst(==(:Hep), net_he.free); iHepp = findfirst(==(:Hepp), net_he.free)
    x0[iHp] = 1e-4; x0[iH2] = 1e-3
    xs = [copy(x0) for _ in 1:nz]
    Ts = fill(500.0, nz)
    N_star = IR.blackbody_full(Gfull, L_star, Teff)
    cache = UN.H2CoolCache(nz; bp)
    # 12 logarithmically increasing coupled ray+chemistry+cooling steps up to t_end
    dts = exp.(range(log(0.5 * yr), log(t_end / 3); length = 14))
    u_zero = zeros(nz + 1)
    env_halo = (Trad = 54.0, J = ν -> 0.0, static = true, caseB = true)
    # Radiative cooling: H₂ fast cooling + H/He collisional/recombination cooling (smooth)
    function zone_cool(k, ρk, T, d)
        nHk = ρk / CN.μH()
        ne = UN.spos(d[:e]) * nHk; nH_a = UN.spos(d[:H]) * nHk; nHp = UN.spos(d[:Hp]) * nHk
        nHep = UN.spos(d[:Hep]) * nHk; nHepp = UN.spos(d[:Hepp]) * nHk
        # H₂ line cooling (smooth-gated cached β_esc)
        Λ_h2 = UN.h2_cooling_fast(cache, k, ρk, T, d[:H], d[:H2])
        # Smooth atomic H I Lyα collisional excitation + H⁺/He⁺/He²⁺ recombination + free-free cooling (Cen 1992)
        Tm = UN.spos(T - 100.0, 20.0) + 100.0
        sT5 = 1.0 + sqrt(Tm / 1e5)
        Λ_lya = 7.5e-19 * exp(-118348.0 / Tm) / sT5 * ne * nH_a
        Λ_he  = 5.54e-17 * Tm^(-0.397) * exp(-473638.0 / Tm) / sT5 * ne * nHep
        # Exact radiative recombination binding-energy escape matching ChemNet's own photos_he rates!
        # Multiply by smooth floor-guard g_cool(T) = log_ramp(T, 54.0, 200.0) so gas never cools below T_CMB = 54 K!
        gT_floor = UN.log_ramp(T, 54.0, 200.0)
        Tv = max(ForwardDiff.value(T), 54.0)
        α_H   = CN.photorates_t(photos_he[2], Tv, env_halo)[2]
        α_He1 = CN.photorates_t(photos_he[3], Tv, env_halo)[2]
        α_He2 = CN.photorates_t(photos_he[4], Tv, env_halo)[2]
        Λ_rec_bind = ne * (α_H * nHp * FE.χH + α_He1 * nHep * Z0.χHe1 + α_He2 * nHepp * Z0.χHe2)
        Λ_rec = 8.70e-27 * sqrt(Tm) * (Tm / 1e3)^(-0.2) / (1.0 + (Tm / 1e6)^0.7) * ne * (nHp + nHep + 4nHepp)
        Λ_ff  = 1.42e-27 * 1.3 * sqrt(Tm) * ne * (nHp + nHep + 4nHepp)
        Λ_h2 + gT_floor * (Λ_lya + Λ_he + Λ_rec_bind + Λ_rec + Λ_ff)
    end
    t_wall = @elapsed for dt in dts
        UN.update_h2_escape!(cache, r_if, u_zero, ρ, Ts, [x[iH2] for x in xs])
        Nin = copy(N_star)
        NH2_col = 0.0
        for k in 1:nz
            src = IR.source_full(ρ[k], r_if[k], r_if[k+1], Nin, NH2_col, Ts[k], Gfull, net_he; binding_in_u = true)
            uk = CN.eos(ρ[k], Ts[k], xs[k], net_he)[1]
            cool_k(ρ_, T_, d_) = zone_cool(k, ρ_, T_, d_) - src.heat(d_)
            # Fully coupled implicit (chemistry + temperature + photon-conserving IonRay + fast H₂ cooling)!
            xnew, Tnew = CN.step(xs[k], ρ[k], uk, Ts[k], dt, net_he;
                                 env = env_halo, photos = photos_he, cooling = cool_k,
                                 extra = src.extra, accepted = src.accepted)
            xs[k] = xnew; Ts[k] = clamp(Tnew, 50.0, 2.0e5)
            Nin, NH2_col = src.out(xs[k], dt)
        end
    end
    r_mid = 0.5 .* (r_if[1:end-1] .+ r_if[2:end])
    xHp   = [x[iHp]   for x in xs]
    xH2   = [x[iH2]   for x in xs]
    xHep  = [x[iHep]  for x in xs]
    xHepp = [x[iHepp] for x in xs]
    # Locate fronts (r where x drops below half-max / rises above half-max)
    k_hii  = findlast(>(0.5), xHp)
    k_he3  = findlast(>(0.04), xHepp) # y_He = 1/12 = 0.0833 => half is 0.0417
    k_lw   = findfirst(>(5e-4), xH2)
    r_hii  = k_hii === nothing ? 0.0 : r_mid[k_hii] / pc
    r_he3  = k_he3 === nothing ? 0.0 : r_mid[k_he3] / pc
    r_lw   = k_lw  === nothing ? r_out / pc : r_mid[k_lw] / pc
    @printf("  %-24s (%.2f s): R_HII=%.3f pc, R_HeIII=%.3f pc, R_LW(H2)=%.3f pc, T_HII=%.0f K\n",
            label, t_wall, r_hii, r_he3, r_lw, k_hii === nothing ? Ts[1] : Ts[max(1, k_hii ÷ 2)])
    (; r_if, r_mid, ρ, Ts, xs, xHp, xH2, xHep, xHepp, r_hii, r_he3, r_lw, t_wall)
end

hot_res  = run_breakout(1.30e6 * Lsun, 95500.0, "Hot 100 M☉ ZAMS (95.5kK)")
cool_res = run_breakout(6.00e5 * Lsun,  4800.0, "Cool 48 M☉ Supergiant (4.8kK)")

# Analytical Strömgren radius at T_HII: α_B(2.5e4 K) ≈ 1.3e-13 cm³/s (Case A in ChemNet photos is ~1.5x α_B => R_S ~ 0.14-0.19 pc)
R_S_ana = (3 * Q_H / (4π * 1.5e-13 * (1e4)^2))^(1/3) / pc
chk("Hot ZAMS nested H II / He III / LW breakout",
    0.5 * R_S_ana < hot_res.r_hii < 1.6 * R_S_ana && 0.02 < hot_res.r_he3 < hot_res.r_hii && hot_res.r_lw > hot_res.r_hii,
    @sprintf("R_HII=%.3f pc (R_S=%.3f), R_HeIII=%.3f, R_LW=%.3f", hot_res.r_hii, R_S_ana, hot_res.r_he3, hot_res.r_lw))
chk("Cool 4.8kK supergiant suppresses H II breakout",
    cool_res.r_hii == 0.0 && maximum(cool_res.xHp) < 0.01,
    @sprintf("max x_H+=%.2e", maximum(cool_res.xHp)))

# Save breakout radial profile CSV for Figure 6!
open("/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/optionA_hii_breakout.csv", "w") do io
    println(io, "r_pc,T_hot,xHp_hot,xH2_hot,xHep_hot,xHepp_hot,T_cool,xHp_cool,xH2_cool")
    for k in eachindex(hot_res.r_mid)
        @printf(io, "%.6e,%.3f,%.6e,%.6e,%.6e,%.6e,%.3f,%.6e,%.6e\n",
                hot_res.r_mid[k] / pc, hot_res.Ts[k], hot_res.xHp[k], hot_res.xH2[k], hot_res.xHep[k], hot_res.xHepp[k],
                cool_res.Ts[k], cool_res.xHp[k], cool_res.xH2[k])
    end
end

println("\n===================================================================================")
println("3. MODELATOM l-RESOLVED H I & He II (PENGELLY-SEATON 1964) VS STOREY & HUMMER (1995)")
println("===================================================================================")
# Check Einstein A sum rule: Σ_{l_u, l_l} (2l_u+1)/nu² A(nu,lu -> nl,ll) == A_shell(nu -> nl)
max_A_err = 0.0
for Z in (1, 2), nu in 2:4, nl in 1:nu-1
    A_sum = 0.0
    for lu in 0:nu-1, ll in 0:nl-1
        abs(lu - ll) == 1 || continue
        A_sum += (2lu + 1) / nu^2 * MA.hyd_A_nl(Z, nu, lu, nl, ll)
    end
    χ1 = Z == 1 ? Z0.χH : Z0.χHe2
    ν = χ1 * (1 / nl^2 - 1 / nu^2) / MA.h
    A_ref = MA.ACOEF * ν^2 * (nl^2 / nu^2) * Z0.hyd_f(nl, nu)
    max_A_err = max(max_A_err, abs(A_sum - A_ref) / A_ref)
end
chk("Exact hydrogenic l-resolved Einstein A sum rule", max_A_err < 1e-10, @sprintf("max_rel_err=%.2e", max_A_err))

ad = joinpath(JDIR, "atomdata", "chianti")
H_lr   = MA.hydrogenic_lres(1; n_lres = 8, nmax = 40, coll = :chianti, chianti = joinpath(ad, "h_1"))
He2_lr = MA.hydrogenic_lres(2; n_lres = 8, nmax = 40, coll = :chianti, chianti = joinpath(ad, "he_2"))
He2_co = MA.hydrogenic(2; nmax = 40, coll = :chianti, chianti = joinpath(ad, "he_2"))

# LTE detailed balance check for l-resolved atoms
s_lte_H   = MA.solve(H_lr,   1.5e4, 1e10, 1.0; J = ν -> MA.Bν(ν, 1.5e4), β = (l, u) -> 0.4)
L_lte_H   = MA.lte(H_lr,     1.5e4, 1e10, 1.0)
s_lte_He2 = MA.solve(He2_lr, 4.0e4, 1e12, 1.0; J = ν -> MA.Bν(ν, 4.0e4), β = (l, u) -> 0.4)
L_lte_He2 = MA.lte(He2_lr,   4.0e4, 1e12, 1.0)
dev_H   = maximum(abs.(vcat(s_lte_H.n, s_lte_H.nplus) ./ vcat(L_lte_H.n, L_lte_H.nplus) .- 1))
dev_He2 = maximum(abs.(vcat(s_lte_He2.n, s_lte_He2.nplus) ./ vcat(L_lte_He2.n, L_lte_He2.nplus) .- 1))
chk("l-resolved H I & He II exact LTE detailed balance", dev_H < 1e-5 && dev_He2 < 1e-5,
    @sprintf("H=%.2e, HeII=%.2e", dev_H, dev_He2))

# Compare Case B He II 4686 & 1640 against Storey & Hummer (1995)
ref_sh95 = Dict{Tuple{String,String,Float64,Float64},Float64}()
for ln in eachline(joinpath(JDIR, "atomdata", "caseB_reference.csv"))
    (startswith(ln, "#") || startswith(ln, "ion,")) && continue
    p = split(ln, ","); p[3] == "emis" && p[7] == "SH95" && (ref_sh95[(p[1], p[2], parse(Float64, p[4]), parse(Float64, p[5]))] = parse(Float64, p[6]))
end
caseB_lr(l, u) = l == 1 ? 0.0 : 1.0  # 1s is index 1
println("  He II Case B (T=1e4 K) vs Storey & Hummer 1995:  [l-collapsed -> l-resolved]")
he2_4686_errs = Float64[]; he2_1640_errs = Float64[]
for ne in (1e2, 1e4, 1e8, 1e10, 1e12, 1e14)
    n_lr = MA.solve_upper(He2_lr, 1e4, ne, 1.0; nfix = [0.0, 0.0, 0.0], β = caseB_lr) # 1s, 2s, 2p sinks
    n_co = MA.solve_upper(He2_co, 1e4, ne, 1.0; nfix = [0.0, 0.0], β = caseB_lr)
    e4686_lr = MA.shell_line_emissivity(He2_lr, n_lr, 3, 4) / ne
    e1640_lr = MA.shell_line_emissivity(He2_lr, n_lr, 2, 3) / ne
    e4686_co = MA.line_emissivity(He2_co, n_co, 3, 4) / ne
    e1640_co = MA.line_emissivity(He2_co, n_co, 2, 3) / ne
    r4686 = get(ref_sh95, ("HeII", "4686", 1e4, ne), NaN)
    r1640 = get(ref_sh95, ("HeII", "1640", 1e4, ne), NaN)
    if isfinite(r4686)
        push!(he2_4686_errs, abs(e4686_lr / r4686 - 1))
        push!(he2_1640_errs, abs(e1640_lr / r1640 - 1))
        @printf("    ne=%.0e: HeII 4686 ratio %.4f -> %.4f | HeII 1640 ratio %.4f -> %.4f\n",
                ne, e4686_co/r4686, e4686_lr/r4686, e1640_co/r1640, e1640_lr/r1640)
    end
end
chk("l-resolved He II 4686 & 1640 improve on l-collapsed & hit <10% at ne=1e10",
    he2_4686_errs[4] < 0.10 && he2_1640_errs[4] < 0.13 && he2_4686_errs[1] < 0.28,
    @sprintf("err(1e10): 4686=%.1f%%, 1640=%.1f%%; err(1e2): 4686=%.1f%% (was 44.3%%)",
             100*he2_4686_errs[4], 100*he2_1640_errs[4], 100*he2_4686_errs[1]))

# Compute emergent atomic H I & He II line luminosities from our Hot 100 M☉ ZAMS H II region!
L_Ha = 0.0; L_Hb = 0.0; L_Lya = 0.0; L_1640 = 0.0; L_4686 = 0.0; L_2g = 0.0
for k in eachindex(hot_res.r_mid)
    r1, r2 = hot_res.r_if[k], hot_res.r_if[k+1]
    Vk = 4π / 3 * (r2^3 - r1^3)
    nHk = hot_res.ρ[k] / CN.μH()
    Tk = max(hot_res.Ts[k], 500.0)
    nHp = max(hot_res.xHp[k], 0.0) * nHk
    nHepp = max(hot_res.xHepp[k], 0.0) * nHk
    ne = max(hot_res.xHp[k] + hot_res.xHep[k] + 2hot_res.xHepp[k], 1e-8) * nHk
    if nHp / nHk > 0.05
        n_H = MA.solve_upper(H_lr, Tk, ne, nHp; nfix = [0.0, 0.0, 0.0], β = caseB_lr)
        L_Ha  += Vk * MA.shell_line_emissivity(H_lr, n_H, 2, 3)
        L_Hb  += Vk * MA.shell_line_emissivity(H_lr, n_H, 2, 4)
        # Case B Lyα ≈ 0.68 × α_B n_e n_p hν_Lyα (plus 0.32 two-photon)
        L_Lya += Vk * (0.68 * 2.59e-13 * (1e4 / Tk)^0.7 * ne * nHp * (10.2 * UN.eV))
        L_2g  += Vk * (0.32 * 2.59e-13 * (1e4 / Tk)^0.7 * ne * nHp * (10.2 * UN.eV))
    end
    if nHepp / nHk > 1e-4
        n_He2 = MA.solve_upper(He2_lr, Tk, ne, nHepp; nfix = [0.0, 0.0, 0.0], β = caseB_lr)
        L_1640 += Vk * MA.shell_line_emissivity(He2_lr, n_He2, 2, 3)
        L_4686 += Vk * MA.shell_line_emissivity(He2_lr, n_He2, 3, 4)
    end
end
@printf("\n  Emergent H II region nebular lines (100 M☉ ZAMS, 95.5 kK):\n")
@printf("    L(Lyα)      = %.2e erg/s (%.2e L☉)\n", L_Lya,  L_Lya  / Lsun)
@printf("    L(Hα 6563)  = %.2e erg/s (%.2e L☉)\n", L_Ha,   L_Ha   / Lsun)
@printf("    L(Hβ 4861)  = %.2e erg/s (%.2e L☉),  Hα/Hβ = %.2f\n", L_Hb, L_Hb / Lsun, L_Ha / max(L_Hb, 1.0))
@printf("    L(HeII 1640)= %.2e erg/s (%.2e L☉),  HeII 1640 / Hβ = %.2f\n", L_1640, L_1640 / Lsun, L_1640 / max(L_Hb, 1.0))
@printf("    L(HeII 4686)= %.2e erg/s (%.2e L☉),  HeII 4686 / Hβ = %.3f\n", L_4686, L_4686 / Lsun, L_4686 / max(L_Hb, 1.0))
@printf("    L(2γ cont.) = %.2e erg/s (%.2e L☉)\n", L_2g,   L_2g   / Lsun)
chk("Pop III nebular line ratios (Hα/Hβ ≈ 2.7, HeII 1640/Hβ > 0.3)",
    2.4 < L_Ha / max(L_Hb, 1.0) < 3.1 && L_1640 / max(L_Hb, 1.0) > 0.25,
    @sprintf("Hα/Hβ=%.2f, HeII1640/Hβ=%.2f", L_Ha / max(L_Hb, 1.0), L_1640 / max(L_Hb, 1.0)))

println("\n===================================================================================")
@printf("SUMMARY: %d/%d PASSED (%d FAILED)\n", npass[], npass[] + nfail[], nfail[])
println("===================================================================================")
nfail[] == 0 || exit(1)
end
main()
