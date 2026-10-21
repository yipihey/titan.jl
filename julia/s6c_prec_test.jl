# S6b Part 2 tests of the thick precursor / ionisation-thermostat surface condition (s6c_precursor.jl).
#   TITAN_NOLIB=1 julia --project=env s6c_prec_test.jl
include(joinpath(@__DIR__, "s4_common.jl"))
include(joinpath(@__DIR__, "s6_star.jl"))
include(joinpath(@__DIR__, "s6_setup.jl"))
include(joinpath(@__DIR__, "s6b_setup.jl"))
include(joinpath(@__DIR__, "s6c_precursor.jl"))
ok = true
check(name, cond) = (global ok &= cond; @printf("  %-86s %s\n", name, cond ? "PASS" : "FAIL"))
# 1. thin limit with no energy sink: E_s equals S6a's Lucy value
let M = 5Msun, R = 30Rsun, Md = 1e-8Msun / YR
    p = precursor_thick(R, M, Md, 1e3Lsun); l = precursor(R, M, Md, 1e3Lsun + GRAV * M * Md / R)
    check("thin limit (Mdot 1e-8): E_s = Lucy to 1e-3", abs(p.Es / l.Es - 1) < 1e-3)
end
# 2. the energy of the shock just ionises the accreted gas at Mdot_fid, M = 0.22, R ≈ 22–28 Rsun (G M/R = 1.5–1.9e13 erg/g ≈ ionisation energy of H + He):
#    S6a's Lucy law gave T_s = 3–5e4 K, tau = 400–1400; with the enthalpy sink the precursor is thin-ish (tau < 100) and T_s ≈ 1e4 K
let M = 0.2174Msun, Md = 4.41e-3Msun / YR
    for R_ in (22.0, 25.0, 28.0)
        R = R_ * Rsun; p = precursor_thick(R, M, Md, 30Lsun)
        @printf("  Mdot_fid, R = %.0f Rsun: T_s %.0f K (Lucy %.0f), tau %.1f, h_R/(GM/R) = %.2f, L_inf/L_in = %.2f\n", R_, p.Ts, (p.Es_lucy / Z0Physics.arad)^0.25, p.tau, p.hR / (GRAV * M / R), p.Linf / (30Lsun + GRAV * M * Md / R))
        check("  energy budget: h_R <= G M/R + L*/Mdot and L_inf >= 0", p.hR <= GRAV * M / R + 30Lsun / Md && p.Linf >= 0)
        check("  T_s below Lucy's (the energy to ionise the gas is not radiated)", p.Ts < (p.Es_lucy / Z0Physics.arad)^0.25)
    end
end
# 3. no cliff: T_s(R) smooth and monotonically decreasing across the ionisation limit R = G M/h_ion ≈ 28 Rsun (steps of 0.5 Rsun), at fixed L*
let M = 0.2174Msun, Md = 4.41e-3Msun / YR
    Rs = 18.0:0.5:34.0; Ts = [precursor_thick(R_ * Rsun, M, Md, 30Lsun).Ts for R_ in Rs]
    check("T_s(R) non-increasing with R to 3 % (no cliff, 0.5 Rsun steps)", all(diff(log.(Ts)) .< 0.03))
    check("largest step of ln T_s per 0.5 Rsun < 0.25", maximum(abs.(diff(log.(Ts)))) < 0.25)
end
# 4. monotone and weak in L* (the lagged L* of the star must not make the boundary hunt): T_s increases with L*, d ln T_s / d ln L* < 0.5
let M = 0.2174Msun, Md = 4.41e-3Msun / YR, R = 21.45Rsun
    Ls = [1.0, 3.0, 10.0, 30.0, 60.0, 100.0, 200.0]; Ts = [precursor_thick(R, M, Md, l * Lsun).Ts for l in Ls]
    check("T_s non-decreasing in L*", all(diff(Ts) .>= -1e-3 .* Ts[1:end-1]))
    check("L* from 10 to 30 Lsun changes T_s by < 5 %", abs(Ts[4] / Ts[3] - 1) < 0.05)
end
# 5. low accretion rate / large R (Mdot 1.1e-3, M = 9, R = 120 Rsun): the enthalpy sink is negligible, E_s = Lucy to a few %
let M = 9Msun, Md = 1.1e-3Msun / YR, R = 120Rsun
    p = precursor_thick(R, M, Md, 1500Lsun); check("Mdot 1.1e-3, 9 Msun, 120 Rsun: T_s within 2 % of Lucy", abs(p.Ts / (p.Es_lucy / Z0Physics.arad)^0.25 - 1) < 0.02)
end
# 6. cost per surface condition
let M = 0.2174Msun, Md = 4.41e-3Msun / YR, R = 25Rsun
    t0 = time(); for _ in 1:20; precursor_thick(R, M, Md, 30Lsun); end; dt = (time() - t0) / 20
    @printf("  time per precursor_thick: %.1f ms\n", 1e3dt); check("cost < 0.5 s", dt < 0.5)
end
println(ok ? "ALL PASS" : "FAILURES")
