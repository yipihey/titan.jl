# FULL-1: max x_alpha / x_c per snapshot of a hi21_snap profile file.  julia full1_xa.jl hi21_<run>_prof.csv
using DelimitedFiles, Printf
A, h = readdlm(ARGS[1], ','; header = true); h = vec(String.(h))
c(n) = A[:, findfirst(==(n), h)]
lab = c("label"); xa = Float64.(c("x_a")); xc = Float64.(c("x_c")); ts = Float64.(c("T_s")); tk = Float64.(c("T_k")); tau = Float64.(c("tau21")); m = Float64.(c("m_Msun")); r = Float64.(c("r_AU"))
for l in unique(lab)
    i = findall(==(l), lab); q = xa[i] ./ xc[i]; k = argmax(q)
    @printf("%s: max x_a/x_c = %.2e at r = %.3g AU (x_c = %.2e, T_s = %.0f, T_k = %.0f); max tau21 = %.2e\n", l, q[k], r[i][k], xc[i][k], ts[i][k], tk[i][k], maximum(tau[i]))
end
