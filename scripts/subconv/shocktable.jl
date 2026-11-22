# Table of the convection diagnostics (cd) of one snapshot around the shock: julia shocktable.jl <snaps.jls> [step] [nbelow=10] [nabove=2]
# One row per interface j (between zones j-1 and j); T = mean of the two zones; branch = zone j's; R* = first interface from the centre with u < -2 km/s.
using Serialization, Printf
function main()
    v = Any[]; open(ARGS[1]) do io; while !eof(io); push!(v, deserialize(io)); end; end
    st = length(ARGS) > 1 ? parse(Int, ARGS[2]) : v[end].step
    s = v[findlast(x -> x.step == st, v)]; cd = s.cd; nb = length(ARGS) > 2 ? parse(Int, ARGS[3]) : 10; na = length(ARGS) > 3 ? parse(Int, ARGS[4]) : 2
    n = length(s.ρ); js = findfirst(<(-2e5), s.u); js === nothing && error("no shock (no u < -2 km/s)")
    Rst = s.r[js]
    @printf("step %d  t = %.6e yr  zones %d  R* = %.4g Rsun (interface %d)  layer [r_in/R*, r_out/R*, m_layer/Msun, max(-N2)tdyn2, max Fconv/Ftot, max Fpot/Ftot, n, R*/Rsun] = %s\n",
            s.step, s.t / 3.156e7, n, Rst / 6.957e10, js, join((@sprintf("%.4g", x) for x in cd.layer), ", "))
    @printf("%4s %7s %9s %6s %4s %10s %10s %9s %9s %8s %8s %6s %6s %5s %10s %10s\n", "j", "r/R*", "T[K]", "Mach", "br", "N2_eq", "N2_L", "nab", "nad_used", "nad_eq", "nad_fr", "w_tau", "w_inf", "w_th", "Fconv/F", "Fpot/F")
    for j in max(2, js - nb):min(n, js + na)
        T = 0.5 * (s.T[j-1] + s.T[j])
        @printf("%4d %7.4f %9.1f %6.3f %4.0f %10.3e %10.3e %9.3g %9.3g %8.4f %8.4f %6.3f %6.3f %5.0f %10.3e %10.3e%s\n", j, s.r[j] / Rst, T, cd.mach_i[j], cd.branch[j], cd.N2_eq[j], cd.N2_L[j],
                cd.nab[j], cd.nad_used[j], cd.nad_eq[j], cd.nad_fr[j], cd.w_tau[j], cd.w_inf[j], cd.w_th[j], cd.Fconv_over_Ftot[j], cd.Fpot_over_Ftot[j], j == js ? "   <- shock" : "")
    end
end
main()
