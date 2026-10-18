# The fed quasi-star in the metrics used for Little Red Dots: Balmer-break flux ratio of the continuum, ⟨f_ν⟩ over
# 4000–4100 Å / ⟨f_ν⟩ over 3620–3720 Å (the windows of de Graaff et al. 2025 and Naidu et al. 2025), L_bol/λL_λ(5100 Å), L_bol/L(Hα), and the electron-scattering depth outside the radius where Hα forms
# (which sets how much electron scattering could broaden the line).
#   Z0_TABLE=ext TITAN_NOLIB=1 julia --project=env lrd_metrics.jl NAME   (spectra NAME and NAME_nlte, nlte_NAME, linespec_NAME)
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf
const Z0 = Z0Physics
name = ARGS[1]; sp(n) = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", n))
Lnu(d, λ) = (ν = Z0.c / (λ * 1e-8); i = clamp(searchsortedlast(d.ν, ν), 1, length(d.ν) - 1); t = (log(ν) - log(d.ν[i])) / (log(d.ν[i+1]) - log(d.ν[i]));
             exp(log(d.L[i]) * (1 - t) + log(d.L[i+1]) * t))
# the jump at the edge itself: L_ν at the grid point just redward of 3646 Å over the one just blueward (the continuum
# grid has a point pair at the edge; interpolating across it, or averaging over a window that straddles it, mixes the two)
edge(d) = (λs = 2.99792458e18 ./ d.ν; i = findlast(l -> l > 3646.0, λs); d.L[i] / d.L[i+1])
win(d, a, b) = sum(Lnu(d, λ) for λ in range(a, b, length = 41)) / 41
Lbol(d) = sum(0.5 * (d.L[i] + d.L[i+1]) * (d.ν[i+1] - d.ν[i]) for i in 1:length(d.ν)-1)
ls = sp("linespec_$name.jls"); LHa = ls.lines[findfirst(x -> x.name == "Hα", ls.lines)].Lline
Lines = sum(x.Lline for x in ls.lines)
for (lab, n) in (("star alone (star_t0)", "star_t0"), ("LTE infall", name), ("non-LTE infall", "$(name)_nlte"))
    d = sp("$n.jls"); Lb = Lbol(d) + (endswith(n, "_nlte") ? Lines : 0.0)
    L5100 = Z0.c / 5100e-8 * Lnu(d, 5100.0)
    @printf("%-22s jump at the edge %5.2f, break ⟨f_ν⟩(4000–4100)/⟨f_ν⟩(3620–3720) = %5.2f   L_bol = %.3e (continuum%s)   L_bol/λL_5100 = %.2f%s\n", lab, edge(d), win(d, 4000, 4100) / win(d, 3620, 3720), Lb,
            endswith(n, "_nlte") ? " + 23 lines" : "", Lb / L5100, endswith(n, "_nlte") ? @sprintf("   L_bol/L(Hα) = %.1f", Lb / LHa) : "")
end
dn = sp("nlte_$name.jls"); R0 = 8.4352e15
ne = [st.ne for st in dn.nl]; τ = zeros(length(ne))
for k in length(ne)-1:-1:1; τ[k] = τ[k+1] + 0.5 * (ne[k] + ne[k+1]) * 6.6524e-25 * (dn.r[k+1] - dn.r[k]); end
# Hα-power-weighted mean of the electron-scattering depth outside the emitting radius, and of the gas temperature
ν32 = (NLTEH.En[3] - NLTEH.En[2]) / Z0.h
w = [4π * dn.r[q]^2 * max(dn.nl[q].P, 0.0) * dn.dV[q] / (4π * dn.r[q]^2) for q in eachindex(dn.nl)]
@printf("line-power-weighted: ⟨τ_es outside⟩ = %.2f, ⟨T⟩ = %.0f K; τ_es outside 5R = %.2f, outside 12R = %.2f\n",
        sum(w .* τ) / sum(w), sum(w .* dn.T) / sum(w), τ[argmin(abs.(dn.r ./ R0 .- 5))], τ[argmin(abs.(dn.r ./ R0 .- 12))])
