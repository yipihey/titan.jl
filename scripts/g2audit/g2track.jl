# g2track.jl <tree> <rundir> [every]: per snapshot of a RUNREC run: step, dt, M*, R*, L_ph, T_eff, the sharpest T jump (the cliff:
# largest |Δ log T| between neighbouring zones with T > 2500 K) and whether it sits on a table|network interface; summary at the end.
include(joinpath(ARGS[1], "julia", "runrec.jl"))
function main(dir, every)
L = lineage(dir)
@printf("%6s %10s %9s %8s %7s %9s %7s %7s | %6s %5s %6s %6s %6s %s\n", "step", "t_yr", "dt_s", "M*", "R*", "L_ph", "Teff_L", "Teff_sl", "cliffk", "dlgT", "T_in", "T_out", "tau_in", "on table|network interface? (nearest interface)")
ons = 0; nn = 0; dts = Float64[]
for (i, x) in enumerate(L)
    F, S = frame(x)
    T = 10 .^ F["logT"]; eq = F["eq"]; τ = 10 .^ F["logtau"]; n = length(T)
    best = 0.0; kc = 0
    for k in 1:n-1
        (T[k] > 2500 && T[k+1] > 2500) || continue
        d = abs(log10(T[k] / T[k+1])); d > best && (best = d; kc = k)
    end
    ifs = [k for k in 1:n-1 if eq[k] != eq[k+1]]
    on = kc > 0 && eq[kc] != eq[kc+1]
    dist = isempty(ifs) ? -1 : minimum(abs.(ifs .- kc))
    push!(dts, x.dt); nn += 1; ons += on
    (i % every == 0 || i == length(L)) && @printf("%6d %10.4f %9.3e %8.5f %7.3f %9.3e %7.0f %7.0f | %6d %5.2f %6.0f %6.0f %6.1e %s (%s)\n",
        S[1], S[2] - 238933.0, S[3], S[6], S[7], S[11], S[9], S[10], kc, best, kc > 0 ? T[kc] : NaN, kc > 0 ? T[kc+1] : NaN, kc > 0 ? τ[kc] : NaN,
        on ? "YES" : "no", isempty(ifs) ? "no interface" : "interfaces at $(ifs), $(dist) zones away")
end
@printf("snapshots %d: cliff on a table|network interface in %d; dt min %.3e median %.3e last %.3e s\n", nn, ons, minimum(dts), sort(dts)[max(1, length(dts) ÷ 2)], dts[end])
end
main(ARGS[2], length(ARGS) > 2 ? parse(Int, ARGS[3]) : 1)
