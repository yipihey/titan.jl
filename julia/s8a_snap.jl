# S8a: Monte Carlo Lyα field of stored envelope snapshots (env_snap_<j>.jls + star_<j>.txt of a coupled run) -> s8a_<tag>_J_<j>.csv for the 21 cm post-processor (hi21_snap.jl, LYA=mc:<dir>:<tag>).
#   julia --project=env -t 4 s8a_snap.jl <dir with env_snap_*.jls and star_*.txt> <tag> [out dir = lib/popiii]
# Same transfer as the coupled run (inner column in the scaled medium, grid Monte Carlo on the snapshot's zones); the snapshot's own zone-1 density is used for the inner column.
using Printf, Serialization
const J_ = @__DIR__
for f in ("PrimordialChem", "ContinuumOpacity", "ChemNet", "ChemRates", "AtomicHHe", "RadTransfer", "s7b_exchange", "s7b_feedback", "LyaMC", "s8a_couple"); include(joinpath(J_, f * ".jl")); end
const RT = RadTransfer; const CN = ChemNet; const S = S7b; const S8 = S8a; const Lm = LyaMC
dir = ARGS[1]; tag = ARGS[2]; outdir = length(ARGS) > 2 ? ARGS[3] : joinpath(J_, "..", "lib", "popiii")
files = sort(filter(f -> occursin(r"^env_snap_\d+\.jls$", f), readdir(dir)))
g = RT.default_grid(dIR = 0.3, dopt = 0.25, dLW = 0.02, dUV = 0.05)
net, photos = AtomicHHe.atomic(); ia = RT.ion_atoms(photos, g; dr = AtomicHHe.k_di_HeII)
function read_rec(path)
    d = Dict{String,Float64}()
    for l in eachline(path); m = match(r"^(\w+) = (\S+)", l); m === nothing || (d[m.captures[1]] = parse(Float64, m.captures[2])); end
    d
end
TIN = 1.5e4
log = open(joinpath(outdir, "s8a_$(tag)_snaps.csv"), "w")
println(log, "j,t_yr,Mpt,Mstar,rate_line,rate_win,inner_Ncol,F_mc_over_Lc,f_out,f_in,f_2g,f_cut,nsc_per_photon,relerr_F,mc_s")
for f in files
    j = match(r"env_snap_(\d+)\.jls", f).captures[1]
    s = deserialize(joinpath(dir, f))[1]; rec = read_rec(joinpath(dir, "star_$(j).txt"))
    n = length(s.T); nH = s.ρ ./ CN.μH(); X = Dict{Symbol,Vector{Float64}}()
    for (nm, i) in ((:Hp, 1), (:Hm, 2), (:H2, 3), (:H2p, 4)); X[nm] = [s.y[k][i] for k in 1:n]; end
    X[:H] = [max(1 - s.y[k][1] - s.y[k][2] - 2s.y[k][3] - 2s.y[k][4], 0.0) for k in 1:n]
    X[:e] = [max(s.y[k][1] + s.y[k][4] - s.y[k][2], 0.0) for k in 1:n]
    src = S.star_source(g, rec; mode = :star)
    inner = S.inner_column(ia, src, rec["M"], rec["R"], s.r[1], s.ρ[1], max(abs(s.u[1]), 1e3); nz = 300, T = TIN, reemit = (g, max(rec["Tph"], 3000.0)))
    Rst = rec["R"]; RII = max(inner.RII, 1.001 * Rst)
    rate_line = Lm.f_lya_caseB(TIN) * max(inner.Qstar - inner.Qesc, 0.0); rate_win = S8.window_rates(g, inner.Lout, RT.LYA_MC_WIN[])
    nHe = s.ρ[1] / CN.μH()
    inn = S8.inner_response!(Rst, RII, s.r[1], nHe, max(abs(s.u[1]), 1e3), TIN)
    Ncol = S8.column(nHe * (1 - S8.XE_IN), RII, s.r[1]); tau0 = Ncol * Lm.sig0(TIN)
    emerge = S8.emergent(inn, rate_line, rate_win, tau0)
    mz = diff(s.m)
    t0 = time(); acc, res = S8.grid_mc!(collect(float.(s.r)), collect(float.(s.T)), collect(float.(s.ρ)), collect(float.(s.u)), X, emerge, Ncol, TIN, mz; seed = parse(Int, j) + 1); tmc = time() - t0
    Lt = (rate_line + sum(rate_win)) * Lm.E0
    big = res === nothing ? Int[] : findall(>(0.02 * maximum(abs.(res.force))), abs.(res.force))
    relerr = isempty(big) ? 0.0 : sqrt(sum(res.force_err[big] .^ 2)) / max(abs(sum(res.force[big])), 1e-300)
    @printf(log, "%s,%.6e,%.4f,%.4f,%.4e,%.4e,%.4e,%.4f,%.4f,%.4f,%.4f,%.4f,%.4e,%.4f,%.1f\n", j, s.t / 3.15576e7, s.m[1] / 1.98847e33, rec["M"] / 1.98847e33, rate_line, sum(rate_win), Ncol, res === nothing ? 0.0 : sum(res.force) * Lm.c / max(Lt, 1e-300),
            res.f_out, res.f_in, res.f_2g, res.f_lost, res.nsc_per_photon, relerr, tmc); flush(log)
    open(joinpath(outdir, "s8a_$(tag)_J_$(parse(Int, j)).csv"), "w") do io
        println(io, "# t = ", s.t); println(io, "zone,r,T,nHI,J_erg,Tc,Tc_valid,Rsc,force,accel")
        for k in 1:n
            println(io, join((k, 0.5 * (s.r[k] + s.r[k+1]), s.T[k], res.nHI[k], Lm.E0 * res.Jalpha[k], res.Tcolor[k], res.Tcolor_valid[k] ? 1 : 0, res.Rsc[k], res.force[k], res.force[k] / mz[k]), ","))
        end
    end
    @printf("snapshot %s: M_pt %.1f, star %.1f Msun, N_H(inner) %.2e, rate_line %.2e, window %.2e s^-1, MC %.1f s, F c/L %.3f\n", j, s.m[1] / 1.98847e33, rec["M"] / 1.98847e33, Ncol, rate_line, sum(rate_win), tmc, sum(res.force) * Lm.c / max(Lt, 1e-300))
end
close(log)
