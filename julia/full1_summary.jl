# FULL-1: tables of the coupled stage (star track, envelope records, inline 21 cm) of the chained run.  julia --project=env full1_summary.jl [tag]
using DelimitedFiles, Printf
tag = isempty(ARGS) ? "abn02z19_fb" : ARGS[1]
L = joinpath(@__DIR__, "..", "lib", "popiii")
function rd(f)
    A, h = readdlm(joinpath(L, f), ','; header = true); Dict(String(h[i]) => (eltype(A[:, i]) <: Number ? Float64.(A[:, i]) : A[:, i]) for i in eachindex(h))
end
st = rd("s6b_e2e_$(tag)_star.csv"); en = rd("s7b_$(tag)_env.csv")
M = st["M"]; t = st["t_yr"] .- st["t_yr"][1]
println("# protostar track: M[Msun] t[yr] Mdot R*[Rsun] L* Lacc Teff_ph Tc fD")
for m in (0.5, 1, 3, 10, 20, 30, 40, 60, 80, 100, 110, 115)
    i = argmin(abs.(M .- m)); abs(M[i] - m) / m < 0.03 && @printf("%7.2f %9.1f %.2e %7.1f %.2e %.2e %7.0f %.3e %.3f\n", M[i], t[i], st["Mdot"][i], st["R_Rsun"][i], st["Lstar_Lsun"][i], st["Lacc_Lsun"][i], st["Teff"][i], st["Tc"][i], st["fD"][i])
end
i = argmax(st["R_Rsun"]); @printf("R* peak %.0f Rsun at M=%.1f, t=%.0f yr\n", st["R_Rsun"][i], M[i], t[i])
for thr in (1e2, 1e3, 1e4); j = findfirst(>(thr), st["Lstar_Lsun"]); @printf("L*>%g at M=%.1f\n", thr, M[j]); end
for thr in (0.99, 0.5, 0.1); j = findfirst(<(thr), st["fD"]); j === nothing || @printf("fD<%g at M=%.1f\n", thr, M[j]); end
lc = st["Lcno"] ./ max.(st["Lstar_Lsun"], 1e-30)
for thr in (0.01, 0.5, 0.9); j = findfirst(>(thr), lc); println(j === nothing ? "CNO never > $(100thr)% of L*" : @sprintf("CNO>%g%% of L* at M=%.1f", 100thr, M[j])); end
@printf("last row: t=%.1f M=%.3f R=%.2f L*=%.3e Lacc=%.3e Teff=%.0f Tc=%.4e rho_c=%.2f fD=%.4f Sion=%.2e Mdot=%.2e\n", t[end], M[end], st["R_Rsun"][end], st["Lstar_Lsun"][end], st["Lacc_Lsun"][end], st["Teff"][end], st["Tc"][end], st["rho_c"][end], st["fD"][end], st["Sion"][end], st["Mdot"][end])
println("\n# envelope records: j t_c M_pt Mdot_win Mstar R* L* Teff* Q* RII/R* Q/Qcrit kLW_max kLW_z1 arad/g F/(Mdot v) Tz1 fH2z1 redge[AU]")
for k in 1:20:length(en["j"])
    @printf("%4d %9.1f %7.2f %.2e %7.2f %7.1f %.2e %.2e %.2e %.2f %.1e %.1e %.1e %.2e %.1e %.0f %.1e %.0f\n", en["j"][k], en["tc_yr"][k], en["Mpt"][k], en["Mdot_win"][k], en["Mstar"][k], en["Rstar_Rsun"][k], en["Lstar"][k], en["Teffstar"][k], en["Qstar"][k], en["RII_over_R"][k], en["Q_over_Qcrit"][k], en["kLW_max"][k], en["kLW_zone1"][k], en["arad_over_g_max"][k], en["Fsum_over_Mdotv"][k], en["Tzone1"][k], en["fH2_zone1"][k], en["r_edge_AU"][k])
end
k = length(en["j"]); @printf("final barrier j=%d: t_c=%.1f M_pt=%.2f Mstar=%.2f Mdot_win=%.2e wall=%.0f s\n", en["j"][k], en["tc_yr"][k], en["Mpt"][k], en["Mstar"][k], en["Mdot_win"][k], en["wall_s"][k])
kk = findlast(>(0), diff(en["Mpt"])); @printf("M_pt last growth at t_c = %.0f yr, M_pt = %.3f\n", en["tc_yr"][kk+1], en["Mpt"][kk+1])
for m in (10, 30, 50, 75, 100, 105)
    j = findfirst(>=(m), en["Mpt"]); j === nothing || @printf("M_pt>=%d: t_c=%.0f Mdot_win=%.2e redge=%.0f AU Tz1=%.0f fH2z1=%.1e F/Mv=%.2e Q*=%.2e RII/R*=%.2f\n", m, en["tc_yr"][j], en["Mdot_win"][j], en["r_edge_AU"][j], en["Tzone1"][j], en["fH2_zone1"][j], en["Fsum_over_Mdotv"][j], en["Qstar"][j], en["RII_over_R"][j])
end
@printf("max Q/Qcrit %.3g, max RII/R* %.3g, max esc fraction %.3g, max arad/g %.3g, max F/(Mdot v) before halt %.3g\n", maximum(en["Q_over_Qcrit"]), maximum(en["RII_over_R"]), maximum(en["Qesc_in"]), maximum(en["arad_over_g_max"]), maximum(en["Fsum_over_Mdotv"][en["Mpt"] .< 100]))
if isfile(joinpath(L, "hi21_e2e_$(tag)_coupled.csv"))
    h = rd("hi21_e2e_$(tag)_coupled.csv")
    println("\n# HI21 inline: step t_yr ncen Tcen Rout MHI Tk_mw Ts_mw Ts_min tau_max dTb_disc_min dTb_disc_max S_core_nJy")
    for k in vcat(1:8:length(h["step"]), length(h["step"]))
        @printf("%6d %.5e %.2e %.0f %.2e %.0f %.0f %.0f %.0f %.2e %.1f %.1f %.2e\n", h["step"][k], h["t_yr"][k], h["ncen"][k], h["Tcen"][k], h["Rout_AU"][k], h["MHI"][k], h["Tk_mw"][k], h["Ts_mw"][k], h["Ts_min"][k], h["tau_max"][k], h["dTb_disc_min"][k], h["dTb_disc_max"][k], h["S_core_nJy"][k])
    end
    k = length(h["step"]); println("last: ", Dict(c => h[c][k] for c in ("step", "t_yr", "Tk_mw", "Ts_mw", "tau_max", "dTb_disc_min", "dTb_disc_max", "dTb_cen_min", "dTb_cen_max", "S_core_nJy", "EW", "teq_max_over_tff")))
end
