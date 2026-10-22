# S7b offline evaluation of the ionising radiation that leaves the star+precursor through the free-fall column inside the sink radius, along a stored S6b protostar track and the
# infall state at the sink of the matching collapse run (no feedback, so these are the uncoupled values the coupled run starts from; the coupled run then closes the loop).
#   TITAN_NOLIB=1 julia --project=env s7b_offline.jl <s6b csv> <events csv> <out csv> [T_inner]
# For each track row: Q_* (blackbody of the stellar surface), Q_ph (S6b photosphere blackbody), the fraction escaping r_s = 5 AU M/Msun in both modes, R_II/R*, Q/Q_crit (OI02), Eddington ratio,
# the Lyα bound on the radiation force relative to the infall momentum flux, and the full radiation momentum rate L/c relative to Ṁ v_ff(R*).
using Printf
const J = @__DIR__
include(joinpath(J, "PrimordialChem.jl")); include(joinpath(J, "ContinuumOpacity.jl")); include(joinpath(J, "ChemNet.jl")); include(joinpath(J, "ChemRates.jl"))
include(joinpath(J, "AtomicHHe.jl")); include(joinpath(J, "RadTransfer.jl")); include(joinpath(J, "s7b_feedback.jl"))
const S = S7b; const RT = RadTransfer; const CN = ChemNet
trackf, evf, outf = ARGS[1], ARGS[2], ARGS[3]; Tin = length(ARGS) > 3 ? parse(Float64, ARGS[4]) : 1.5e4; rhomode = length(ARGS) > 4 ? ARGS[5] : "grid"      # grid: ρ_s of the zone outside the sink in the collapse run; mdot: ρ_s = Ṁ_dyn/(4π r_s² u_s)
rd(f) = [parse.(Float64, split(l, ",")[1:end]) for l in readlines(f)[2:end] if !isempty(strip(l))]
track = rd(trackf); ev = rd(evf)
# columns of the track: 1 t_yr, 3 M, 4 R_Rsun, 5 Lstar, 6 Lacc, 7 Teff(ph), 23 Linf, 25 Mdot
evM = [e[3] for e in ev]; evrho = [e[6] for e in ev]; evu = [abs(e[5]) * 1e5 for e in ev]
function interp_ev(M, v)
    M <= evM[1] && return v[1]; M >= evM[end] && return v[end]
    j = searchsortedlast(evM, M); w = (log(M) - log(evM[j])) / (log(evM[j+1]) - log(evM[j])); exp((1 - w) * log(v[j]) + w * log(v[j+1]))
end
net, photos = AtomicHHe.atomic()
g = RT.default_grid(dUV = 0.05); ia = RT.ion_atoms(photos, g; dr = AtomicHHe.k_di_HeII)
io = open(outf, "w")
println(io, "t_yr,M,Mdot,R_Rsun,Lstar,Linf,Teff_ph,Teff_star,Q_star,Q_esc_star,RII_Rstar_star,Qcrit,Q_over_Qcrit,Q_ph,Q_esc_ph,RII_Rstar_ph,Gamma_Edd,Lc_over_Mdotvff,Lya_bound_over_ram,rs_AU,rho_s,u_s_kms")
function main()
last_M = 0.0
for (i, row) in enumerate(track)
    M = row[3]; M < 0.9 && continue; M < last_M * 1.03 && continue; last_M = M
    Ms = M * S.Msun; R = row[4] * S.Rsun; Ls = max(row[5], 1e-3) * S.Lsun; Lacc = row[6] * S.Lsun; Tph = row[7]; Linf = max(row[23], row[5]) * S.Lsun; Mdot = row[25] * S.Msun / S.YR
    rs = max(5.0 * M, 0.05) * S.AU; u_s = interp_ev(M, evu); rho_s = rhomode == "mdot" ? Mdot / (4π * rs^2 * u_s) : interp_ev(M, evrho)
    rec = Dict("Lstar" => Ls, "R" => R, "Linf" => Linf, "Tph" => Tph)
    srcs = S.star_source(g, rec; mode = :star); srcp = S.star_source(g, rec; mode = :photosphere)
    cs = S.inner_column(ia, srcs, Ms, R, rs, rho_s, u_s; T = Tin, reemit = (g, Tph)); cp = S.inner_column(ia, srcp, Ms, R, rs, rho_s, u_s; T = Tin, reemit = (g, Tph))
    vff = sqrt(2 * S.G * Ms / R); Mdotenv = cs.Mdot
    # Lyα bound at the H II edge: proton density of the free-fall gas at R_II
    nRII = rho_s * (cs.RII / rs)^(-1.5) / CN.μH()
    LLy, Nsc, FLy = S.lya_force_bound(cs.Qstar - cs.Qesc, nRII; T = Tin)
    ram = Mdotenv * sqrt(2 * S.G * Ms / cs.RII)
    println(io, join([row[1], M, Mdotenv * S.YR / S.Msun, row[4], row[5], Linf / S.Lsun, Tph, (Ls / (4π * R^2 * S.σSB))^0.25, cs.Qstar, cs.Qesc, cs.RII / R, cs.Qcrit, cs.Qstar / cs.Qcrit, cp.Qstar, cp.Qesc, cp.RII / R,
                      S.Γedd(Ls + Lacc, Ms), Linf / S.c / (Mdotenv * vff), FLy / ram, rs / S.AU, rho_s, u_s / 1e5], ","))
    @printf("M = %7.2f  R = %7.1f Rsun  Mdot(env) = %.2e  Q* = %.2e (esc %.1e)  R_II/R* = %.3f  Q/Qcrit = %.2e | photosphere: Q = %.2e (esc %.1e)  Γ = %.2f  L/c/(Ṁ v_ff) = %.3f  Lyα/ram ≤ %.2e\n",
            M, row[4], Mdotenv * S.YR / S.Msun, cs.Qstar, cs.Qesc / max(cs.Qstar, 1e-300), cs.RII / R, cs.Qstar / cs.Qcrit, cp.Qstar, cp.Qesc / max(cp.Qstar, 1e-300), S.Γedd(Ls + Lacc, Ms), Linf / S.c / (Mdotenv * vff), FLy / ram)
    flush(io)
end
end
main()
close(io)
