# S6c: Mdot(M*) of sink runs (smoothed tables of s6c_mdot.jl) side by side, with OHY10 (Z = 0, Fig. B16 read off).
#   julia s6c_atomic_cmp.jl TAG1 TAG2 ...   (files lib/popiii/s6c_<TAG>_mdot_smooth.csv)
using Printf
const OHY = [(0.003, 0.5), (0.01, 0.2), (0.03, 0.12), (0.1, 0.085), (0.3, 0.068), (1.0, 0.055), (3.0, 0.045), (10.0, 0.035), (30.0, 0.03), (100.0, 0.02)]
function interp(tab, M)
    (M < tab[1][1] || M > tab[end][1]) && return NaN
    i = findlast(p -> p[1] <= M, tab); i == length(tab) && return tab[end][2]
    w = log(M / tab[i][1]) / log(tab[i+1][1] / tab[i][1]); exp((1 - w) * log(tab[i][2]) + w * log(tab[i+1][2]))
end
load(tag) = [Tuple(parse.(Float64, split(l, ","))) for l in readlines(joinpath(@__DIR__, "..", "lib", "popiii", "s6c_$(tag)_mdot_smooth.csv"))[2:end] if !isempty(strip(l))]
tags = ARGS; tabs = [load(t) for t in tags]
@printf("%8s", "M [Msun]"); for t in tags; @printf(" %12s", t); end; @printf(" %10s\n", "OHY10 Z=0")
for M in (0.03, 0.1, 0.3, 1, 3, 10, 20, 40, 100)
    @printf("%8.2g", M); for tb in tabs; @printf(" %12.3e", interp(tb, M)); end; @printf(" %10.3e\n", interp(OHY, M))
end
