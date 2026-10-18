# The outer infall of a wide-domain run (W30d): temperature at fixed radii, mean and spread over the snapshots after
# TMIN, against dilute radiative equilibrium with the photosphere, T⁴ = T_eff⁴ W(r), W = ½[1 − √(1 − (R/r)²)], and
# against the T ∝ r^-1/2 law from T(1.2R) assumed for the free-fall extension in Spectrum.atmosphere.
#   julia --project=env outer_infall.jl RUN [TMIN=2e8]
using Serialization, Printf
S = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(ARGS[1])_relax_snaps.jls"))
tmin = length(ARGS) > 1 ? parse(Float64, ARGS[2]) : 2e8
S = filter(s -> s.t >= tmin, S); R = 8.435236237964494e15; Teff = 7310.0
W(x) = 0.5 * (1 - sqrt(1 - 1 / x^2))
xs = [1.1, 1.2, 1.5, 2, 3, 4, 6, 8, 12, 16, 20, 25, 29]
@printf("%d snapshots, t = %.3e – %.3e s\n", length(S), S[1].t, S[end].t)
Ts = [[begin rc = [0.5(s.r[k] + s.r[k+1]) / R for k in eachindex(s.d)]
        i = clamp(searchsortedlast(rc, x), 1, length(rc) - 1); t = (log(x) - log(rc[i])) / (log(rc[i+1]) - log(rc[i]))
        exp(log(s.T[i]) * (1 - t) + log(s.T[i+1]) * t) end for s in S] for x in xs]
T12 = sum(Ts[2]) / length(S)
@printf("%6s %8s %6s %8s %8s %8s %8s\n", "r/R", "T mean", "sd", "min", "max", "diluteRE", "r^-1/2")
for (x, T) in zip(xs, Ts)
    m = sum(T) / length(T); sd = sqrt(sum((T .- m) .^ 2) / length(T))
    @printf("%6.1f %8.0f %6.0f %8.0f %8.0f %8.0f %8.0f\n", x, m, sd, minimum(T), maximum(T), Teff * W(x)^0.25, T12 * sqrt(1.2 / x))
end
