# S6b Part 1: infall state just outside the sink (zone above the interface that crossed r_s) from lib/popiii/s6c_<TAG>_events.csv.
#   julia s6c_infall.jl <TAG> [M list]
using Printf
const G = 6.6743e-8; const Msun = 1.98847e33; const AU = 1.495978707e13; const kB = 1.380649e-16; const mH = 1.6726e-24
tag = ARGS[1]
Ms = length(ARGS) > 1 ? parse.(Float64, ARGS[2:end]) : [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100]
L = readlines(joinpath(@__DIR__, "..", "lib", "popiii", "s6c_$(tag)_events.csv"))[2:end]
ev = sort([parse.(Float64, split(l, ",")) for l in L if !isempty(strip(l))], by = e -> e[3])
# columns: step, t_cross_yr, m_interface, r_s_AU, u_kms, rho, T, fH2, xe, e, p, cs_kms, rho_up, T_up, u_up_kms
println(tag, ": state of the zone just outside the sink at the crossing of the shell M (events within a factor 1.1 in M averaged in ln)")
@printf("%9s %8s %10s %9s %9s %8s %9s %9s %8s %8s %9s %10s\n", "M [Msun]", "r_s [AU]", "rho", "T [K]", "u [km/s]", "v_ff", "u/v_ff", "c_s", "Mach", "f_H2", "x_e", "n_H [cm-3]")
for M in Ms
    sel = [e for e in ev if abs(log(e[3] / M)) < log(1.12)]
    isempty(sel) && continue
    avg(i) = exp(sum(log(abs(e[i]) + 1e-300) for e in sel) / length(sel))
    am(i) = sum(e[i] for e in sel) / length(sel)
    rs = avg(4); u = avg(5); vff = sqrt(2G * M * Msun / (rs * AU)) / 1e5; cs = avg(12)
    @printf("%9.3g %8.3g %10.3e %9.0f %9.2f %8.2f %9.2f %9.2f %8.2f %8.3f %9.2e %10.2e\n", M, rs, avg(6), avg(7), u, vff, u / vff, cs, u / cs, am(8), am(9), avg(6) * 0.75 / mH)
end
