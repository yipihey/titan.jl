# Check of mc_escatter.jl: with electron scattering off (es=0, Sobolev line resonances on), the Monte Carlo Hα profile
# should match the ray-integrated profile of line_profile.jl (LP_NLTE) for the same non-LTE solution, F/F_c − 1
# normalized to unit area. Prints FWHM, centroid, the fractions beyond ±500 km/s and the L1 distance of the two shapes.
#   julia --project=env mc_validate.jl NAME
using Serialization, Printf
name = ARGS[1]
l = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "lines_$(name)_nlte.jls")); ha = l[findfirst(z -> z.name == "Halpha", l)]
m = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "mc_$(name)_Halpha_es0_res1.jls"))
v = ha.v; e = max.(ha.nlte .- 1, 0.0); dv = v[2] - v[1]; e ./= sum(e) * dv          # ray profile (emission part), unit area
pm = [m.V[1] <= x <= m.V[end] ? (i = clamp(searchsortedlast(m.V, x), 1, length(m.V) - 1); t = (x - m.V[i]) / (m.V[i+1] - m.V[i]); m.prof[i] * (1 - t) + m.prof[i+1] * t) : 0.0 for x in v]
pm ./= sum(pm) * dv
stats(p) = (k = argmax(p); h = findall(>(p[k] / 2), p); (v[h[end]] - v[h[1]], sum(v .* p) * dv, sum(p[v .< -500]) * dv, sum(p[v .> 500]) * dv))
a, b = stats(e), stats(pm)
@printf("ray-integrated: FWHM %.0f km/s, centroid %+.0f, beyond −500 %.3f, beyond +500 %.3f\n", a...)
@printf("Monte Carlo   : FWHM %.0f km/s, centroid %+.0f, beyond −500 %.3f, beyond +500 %.3f\n", b...)
@printf("L1 distance of the normalized shapes: %.3f (0 = identical, 2 = disjoint)\n", sum(abs.(e .- pm)) * dv)
