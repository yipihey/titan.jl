# Bound on the Balmer break for lines above n = NMAX: the region between the Balmer edge and the n = NMAX line of a
# balmer_series.jl result set to at most f × the continuum (f = 0.18, the depth reached just redward; f = 0: black).
#   julia --project=env balmer_bound.jl [NAME_TAG=W30_late_ext_n40_b12] [nmax=40]
using Serialization, Printf
nmax = parse(Int, get(ARGS, 2, "40"))
b = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "balmer_$(get(ARGS, 1, "W30_late_ext_n40_b12")).jls")); λ = b.λ
win(y, a, c) = (s = (λ .>= a) .& (λ .<= c); sum(y[s]) / count(s))
brk(y) = win(y, 4000, 4100) / win(y, 3620, 3720)
λe = 1e8 / (109677.58 / 4); λ40 = λe / (1 - 4 / nmax^2)
@printf("n=%d last line at %.1f Å; break full %.2f\n", nmax, λ40, brk(b.Lnu))
for f in (0.18, 0.0)
    y = copy(b.Lnu); s = (λ .> λe) .& (λ .< λ40); y[s] .= min.(y[s], f .* b.Lc[s])
    @printf("upper bound, %.2f of the continuum between the edge and the n=%d line: break %.2f\n", f, nmax, brk(y))
end
