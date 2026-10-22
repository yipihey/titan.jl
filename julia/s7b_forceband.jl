# S7b: which frequencies carry the radiation force on the neutral envelope? Acceleration per energy band for a star record (M, R*, L*, L_inf, T_ph) in a cold neutral envelope zone.
#   TITAN_NOLIB=1 julia --project=env s7b_forceband.jl <M> <R_Rsun> <Lstar_Lsun> <Linf_Lsun> <Tph> [r_AU] [nH]
using Printf
const J = @__DIR__
include(joinpath(J, "PrimordialChem.jl")); include(joinpath(J, "ContinuumOpacity.jl")); include(joinpath(J, "ChemNet.jl")); include(joinpath(J, "ChemRates.jl"))
include(joinpath(J, "AtomicHHe.jl")); include(joinpath(J, "RadTransfer.jl")); include(joinpath(J, "s7b_exchange.jl")); include(joinpath(J, "s7b_feedback.jl"))
const S = S7b; const RT = RadTransfer; const CN = ChemNet
M, R, Ls, Linf, Tph = parse.(Float64, ARGS[1:5]); rAU = length(ARGS) > 5 ? parse(Float64, ARGS[6]) : 3000.0; nH = length(ARGS) > 6 ? parse(Float64, ARGS[7]) : 1e7
net, photos = AtomicHHe.atomic(); g = RT.default_grid(dIR = 0.3, dopt = 0.25, dLW = 0.02, dUV = 0.05); ia = RT.ion_atoms(photos, g; dr = AtomicHHe.k_di_HeII)
rec = Dict("Lstar" => Ls * S.Lsun, "R" => R * S.Rsun, "Linf" => Linf * S.Lsun, "Tph" => Tph)
src = S.star_source(g, rec); inn = S.inner_column(ia, src, M * S.Msun, R * S.Rsun, 1000 * S.AU, 3e-17, 1e6; reemit = (g, Tph))
ab = RT.absorbers(photos, g)
n = 1; r = [rAU * S.AU, 1.05 * rAU * S.AU]; X = Dict{Symbol,Vector{Float64}}(:H => [1.0], :H2 => [1e-3], :He => [CN.yHe], :e => [1e-5], :Hp => [1e-5])
for T in (1800.0,)
    med = RT.Medium(r, [T], [nH], X)
    # the direct beam enters with the spectrum Lout diluted geometrically to r (the transfer treats a source at the centre: use the unattenuated inner-column output)
    f = RT.transfer(g, ab, RT.Source(inn.Lout, 0.0, sum(inn.Lout)), med, RT.Opts(diffuse = false); rho = [nH * CN.μH()])
    op = RT.new_opacity(1, RT.nfreq(g), length(ab)); RT.fill_opacity!(op, g, ab, med, RT.Opts(diffuse = false))
    emis = sum(4π * op.kother[1, k] * RT.Bν(g.ν[k], T) * g.dν[k] for k in 1:RT.nfreq(g))
    @printf("heating by the non-network continua: absorption of the star %.3e, thermal emission %.3e erg cm^-3 s^-1 (per H nucleus: %.3e / %.3e erg/s)\n", f.heat_other[1], emis, f.heat_other[1] / nH, emis / nH)
    bandsH = [(0.1, 1.0), (1.0, 4.0), (4.0, 10.2), (10.2, 13.6), (13.6, 100.0)]
    for (lo, hi) in bandsH
        h = sum(4π * op.kother[1, k] * f.J[1, k] * g.dν[k] for k in 1:RT.nfreq(g) if g.E[k] >= lo && g.E[k] < hi; init = 0.0)
        @printf("   absorbed power in %5.1f-%5.1f eV: %.3e erg cm^-3 s^-1\n", lo, hi, h)
    end
    tot = f.arad[1]; gl = S.G * M * S.Msun / (rAU * S.AU)^2
    @printf("M = %.1f, r = %.0f AU, n_H = %.1e: a_rad = %.3e, g = %.3e (a/g = %.3f); thin-limit flux F = L/(4π r²) with total extinction\n", M, rAU, nH, tot, gl, tot / gl)
    # band decomposition of (L_in - L_out)/c per mass in the zone
    bands = [(0.1, 4.0), (4.0, 10.2), (10.2, 13.6), (13.6, 24.6), (24.6, 100.0)]
    for (lo, hi) in bands
        s = 0.0; Lb = 0.0; kb = 0.0
        for k in 1:RT.nfreq(g)
            E = g.E[k]; (E >= lo && E < hi) || continue
            s += (f.Ldir[1, k] - f.Ldir[2, k]) / (CN.c * nH * CN.μH() * RT.volume(r, 1)); Lb += inn.Lout[k]; kb += (op.kabs[1, k] + op.ssc[1, k]) / (nH * CN.μH()) * inn.Lout[k]
        end
        @printf("  %5.1f - %5.1f eV: L_band = %.3e Lsun (%.1f%%), a = %.3e (%.1f%% of total), mean extinction %.2e cm²/g\n", lo, hi, Lb / S.Lsun, 100 * Lb / sum(inn.Lout), s, 100 * s / tot, Lb > 0 ? kb / Lb : 0.0)
    end
end
