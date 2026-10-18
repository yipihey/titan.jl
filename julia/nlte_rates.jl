# What ionizes the infall in a non-LTE solution (nlte_run.jl): one continuum transfer with the final populations,
# then at selected radii the photoionization rate from each level, n_i ∫ 4π σ_i J_ν/(hν) dν, as a fraction of the
# total, and the mean intensity relative to the geometrically diluted star, J_ν / [W(r) I_ν], at a few wavelengths
# (I_ν: blackbody at T_eff, the reference of the simple model, analytic_infall.jl).
#   Z0_TABLE=ext TITAN_NOLIB=1 julia -t 2 --project=env nlte_rates.jl NAME
# Writes lib/spectra/nlterates_NAME.jls with J_ν on the frequency grid at every non-LTE point.
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf
const SP = Spectrum; const Z0 = Z0Physics; const NH = NLTEH
name = ARGS[1]
d = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "$name.jls"))
dn = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$name.jls"))
atm = merge(d.atm, (T = dn.Tstar,)); sel = dn.sel; R0 = 8.4352e15
nl = Any[nothing for _ in eachindex(atm.r)]; for (q, j) in enumerate(sel); nl[j] = dn.nl[q]; end
op = NH.NOpac(SP.Opac(atm), nl)
_, R = SP.rays(atm.r; ncore = 12); R = SP.weights!(atm, R)
ν = dn.ν
Jm = zeros(length(atm.r), length(ν))
Threads.@threads for i in eachindex(ν)
    _, J, _ = SP.monochromatic(atm, op, R, ν[i])
    Jm[:, i] .= J
end
wν = [i == 1 ? 0.5(ν[2] - ν[1]) : i == length(ν) ? 0.5(ν[end] - ν[end-1]) : 0.5(ν[i+1] - ν[i-1]) for i in eachindex(ν)]
Teff = 7310.0; L = 1.4204e44; Rst = sqrt(L / (4π * Z0.σSB * Teff^4))
W(r) = 0.5 * (1 - sqrt(max(1 - (Rst / r)^2, 0.0)))
Jat(j, λ) = (x = Z0.c / (λ * 1e-8); i = clamp(searchsortedlast(ν, x), 1, length(ν) - 1);
             Jm[j, i] + (Jm[j, i+1] - Jm[j, i]) * (x - ν[i]) / (ν[i+1] - ν[i]))
@printf("%7s %7s | %-34s | %s\n", "r/R0", "x_H", "photoionizations from n = 1, 2, 3, 4, ≥5", "J/(W B(T_eff)) at 1000, 3000, 3600, 4000, 6563, 8000 Å")
NET = []
for x in (0.995, 1.05, 1.2, 1.6, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 50.0, 80.0, 115.0)
    q = argmin(abs.(dn.r ./ R0 .- x)); j = sel[q]; st = dn.nl[q]; r = atm.r[j]
    rate = zeros(NH.NL)
    for i in 1:NH.NL, k in eachindex(ν)
        Z0.h * ν[k] > 100Z0.eV && continue
        σ = NH.σbf(i, ν[k]); σ == 0 && continue
        rate[i] += st.n[i] * 4π * wν[k] * σ * Jm[j, k] / (Z0.h * ν[k])
    end
    f = rate ./ sum(rate)
    # net radiative ionization per level: photoionization − radiative (spontaneous + induced) recombination, per volume
    T = atm.T[j]; ne = st.ne
    recb = [st.np * ne * NH.Φ(i, T) * sum((Z0.h * ν[k] > 100Z0.eV ? 0.0 : 4π * wν[k] * NH.σbf(i, ν[k]) * (2Z0.h * ν[k]^3 / Z0.c^2 + Jm[j, k]) *
            exp(-Z0.h * ν[k] / (Z0.kB * T)) / (Z0.h * ν[k])) for k in eachindex(ν)) for i in 1:NH.NL]
    net = rate .- recb
    push!(NET, (r = r / R0, net = net, tot = sum(recb)))
    @printf("%7.3f %7.4f | %6.3f %6.3f %6.3f %6.3f %6.3f | %s\n", r / R0, st.np / (st.np + sum(st.n)), f[1], f[2], f[3], f[4], sum(f[5:end]),
            join([@sprintf("%7.3f", Jat(j, λ) / (W(r) * SP.Bν(Z0.c / (λ * 1e-8), Teff))) for λ in (1000.0, 3000.0, 3600.0, 4000.0, 6563.0, 8000.0)], " "))
end
println("net radiative ionizations per level (photoionization − radiative recombination) / total radiative recombinations:")
for q in NET
    @printf("  r/R0 %7.3f  n=1 %+8.4f  n=2 %+8.4f  n=3 %+8.4f  n=4 %+8.4f  n≥5 %+8.4f\n", q.r, (q.net[1:4] ./ q.tot)..., sum(q.net[5:end]) / q.tot)
end
serialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlterates_$name.jls"), (ν = ν, r = atm.r[sel], J = Jm[sel, :]))
