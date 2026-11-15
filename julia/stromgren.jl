# Strömgren budget of a RUNREC snapshot: cumulative case-B recombinations from r0 outward, with the snapshot's own n_e n_p and fully
# ionized (n_e = 1.08 n_H, n_p = n_H), against a photon rate Q.   julia --project=<env> stromgren.jl <snaps.jls> <r0 R☉> <Q s⁻¹> [last|step]
using Serialization, Printf
include(joinpath(@__DIR__, "runrec.jl"))
v = readstream(ARGS[1]); x = (length(ARGS) < 4 || ARGS[4] == "last") ? v[end] : v[findfirst(y -> y.step == parse(Int, ARGS[4]), v)]
r0 = parse(Float64, ARGS[2]) * Rs; Q = parse(Float64, ARGS[3]); μH = 1.6735e-24 * (1 + 4 / 12)
αB(T) = (t = T / 1e4; 4.309e-13 * t^-0.6166 / (1 + 0.6703 * t^0.53))          # RECFAST fit, as PopIIIMaterial.alphaB_recfast
ie = findfirst(==(:e), x.spn); ip = findfirst(==(:Hp), x.spn)
s1 = 0.0; s2 = 0.0; done1 = false; done2 = false
println("step $(x.step): Q = $Q from r0 = $(ARGS[2]) R☉")
for k in eachindex(x.ρ)
    x.r[k] >= r0 || continue
    nH = x.ρ[k] / μH; V = 4π / 3 * (x.r[k+1]^3 - x.r[k]^3)
    global s1 += αB(x.T[k]) * nH^2 * x.Y[ie, k] * x.Y[ip, k] * V; global s2 += αB(max(x.T[k], 1e4)) * nH^2 * 1.083 * V
    !done1 && s1 >= Q && (@printf("  snapshot ionization: recombinations reach Q at r = %.4g R☉ (%.3g pc), zone %d\n", x.r[k+1] / Rs, x.r[k+1] / 3.0857e18, k); global done1 = true)
    !done2 && s2 >= Q && (@printf("  fully ionized (T ≥ 1e4 K): recombinations reach Q at r = %.4g R☉ (%.3g pc), zone %d\n", x.r[k+1] / Rs, x.r[k+1] / 3.0857e18, k); global done2 = true)
end
@printf("  totals to the edge (%.3g pc): snapshot %.3e, fully ionized %.3e s⁻¹\n", x.r[end] / 3.0857e18, s1, s2)
