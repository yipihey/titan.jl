# Along run profiles: network (HM88 + H2PF, LTE) vs table (SCVHXE) EOS in the network domain (s < 1) and in the blend (1 ≤ s < 3).
include(joinpath(@__DIR__, "profiles.jl"))    # (prints the profiles too; defines readstate, sval, wval, nonideal, net_pu_hm)
println("\n==== verification")
for (lab, path) in ((ARGS[i], ARGS[i+1]) for i in 1:2:length(ARGS))
    st = readstate(path); ρ = vec(st["dn"]); T = vec(st["tn"])
    res = Dict("net" => Float64[], "blend" => Float64[]); worst = Dict{String,Any}()
    for k in eachindex(ρ)
        (ρ[k] > 2e-20 && 1e3 <= T[k] <= 3e5) || continue
        nd = nonideal(ρ[k], T[k]); s = sval(nd); s >= 3 && continue
        key = s < 1 ? "net" : "blend"
        tt = thermo(tab_pu, ρ[k], T[k]); tn = thermo(net_pu_hm, ρ[k], T[k])
        e = (abs(log(tn.p / tt.p)), abs(tn.u - tt.u) * ρ[k] / tt.p, abs(log(tn.Γ1 / tt.Γ1)), abs(log(tn.cP / tt.cP)))
        push!(res[key], maximum(e[1:2]))
        w = get(worst, key, nothing)
        (w === nothing || maximum(e[1:2]) > maximum(w[2][1:2])) && (worst[key] = (k, e, log10(ρ[k]), log10(T[k])))
    end
    for key in ("net", "blend")
        isempty(res[key]) && continue
        w = worst[key]
        @printf("%-14s %-6s zones %3d: max(|Δln p|, |Δu| kT) median %.3f max %.3f; worst zone %d (log ρ %.2f log T %.2f): Δln p %.3f Δu %.3f kT Δln Γ1 %.3f Δln cP %.3f\n",
                lab, key, length(res[key]), sort(res[key])[cld(end, 2)], maximum(res[key]), w[1], w[3], w[4], w[2]...)
    end
end
