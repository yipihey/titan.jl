# Network (HM88, LTE) vs table (SCVHXE) EOS on (log ρ, log T) with the non-ideality parameters; OUT/gridnd.bin + gridnd.txt.
ENV["HM88"] = "1"
include(joinpath(@__DIR__, "load.jl")); include(joinpath(@__DIR__, "mats.jl")); include(joinpath(@__DIR__, "mats_hm.jl"))
const OUT = get(ENV, "OUT", joinpath(@__DIR__, "out_nd")); mkpath(OUT)
const LR = collect(-14.0:0.1:1.0); const LTg = collect(3.0:0.02:5.5)
const FIELDS = [:p_t, :u_t, :G1_t, :cP_t, :p_n, :u_n, :G1_n, :cP_n, :xe_t, :xe_n, :Gam, :Deg, :PHM, :w_scvh, :Pneu]
EQ.table(); net_lte_hm(1e-6, 1e4)
const A = fill(NaN, length(FIELDS), length(LR), length(LTg))
t0 = time()
Threads.@threads for i in eachindex(LR)
    for j in eachindex(LTg)
        ρ = 10^LR[i]; T = 10^LTg[j]
        try
            tt = thermo(tab_pu, ρ, T); tn = thermo(net_pu_hm, ρ, T); s = net_lte_hm(ρ, T); nd = nonideal(ρ, T)
            y = EQ.composition(ρ, T); z = EQ.composition_he(ρ, T)
            tb = EQ.table(); x = EQ.gx(tb, LR[i]); yy = EQ.gy(tb, LTg[j])
            ii = clamp(floor(Int, x), 1, tb.nρ - 1); jj = clamp(floor(Int, yy), 1, tb.nT - 1); a = x - ii; b = yy - jj
            w = (1 - a) * (1 - b) * tb.w[ii, jj] + a * (1 - b) * tb.w[ii+1, jj] + (1 - a) * b * tb.w[ii, jj+1] + a * b * tb.w[ii+1, jj+1]
            v = (tt.p, tt.u, tt.Γ1, tt.cP, tn.p, tn.u, tn.Γ1, tn.cP, y[1] + z[1] + 2z[2], s.ne / s.NH, nd.Γ, nd.D, nd.PHM, w, nd.Pneu)
            for (k, val) in enumerate(v); A[k, i, j] = val; end
        catch e
        end
    end
end
@printf("gridnd %d × %d in %.0f s\n", length(LR), length(LTg), time() - t0)
open(joinpath(OUT, "gridnd.bin"), "w") do io; write(io, A); end
open(joinpath(OUT, "gridnd.txt"), "w") do io
    println(io, join(string.(FIELDS), " ")); println(io, join(string.(LR), " ")); println(io, join(string.(LTg), " "))
end
