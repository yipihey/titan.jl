# S8b: 21 cm diagnostic of stored snapshots (post-processor).
#   julia --project=env hi21_snap.jl SNAPFILE[,SNAPFILE2,...] [every]
# Environment:
#   RUN=<name>        output lib/popiii/hi21_<RUN>.csv (summary, one row per snapshot), hi21_<RUN>_prof.csv (per zone), hi21_<RUN>_spec_<k>.csv (velocity-resolved spectra); default: file name
#   TCMB (54.5), Z (19)    background temperature at the source and its redshift (z = 19 ↔ T_CMB = 2.7255 (1+z) = 54.5 K)
#   LYA=none|scan|rt|file   none: no Lyα (collisions + CMB); scan:<f>: uniform J_α = f × J_αC (the critical flux, x_α = 1 for S_α = 1); rt: RadTransfer's zone mean intensity at the Lyα bin
#                           from the star record next to the snapshot (STARDIR, env_snap_<j>.jls ↔ star_<j>.txt), continuum photons (Hirata S̃_α, T_c); file:<csv> per-zone columns J (erg…) or P, Tc
#   NEQ=1             time-dependent level population carried from snapshot to snapshot (Lagrangian zones)
#   SPEC=1            write spectra files;  RAYS=0 skips the ray tracing (zone quantities only);  KH2: H–H₂ coefficient as a multiple of κ_HH (default 0)
#   RADIO_AR, RADIO_BETA   radio excess at 78.1 MHz (default none)
using Printf, Serialization
const J_ = @__DIR__
include(joinpath(J_, "HI21.jl")); include(joinpath(J_, "Perf.jl"))
const H = HI21
const ROOT = joinpath(J_, "..")
files = split(ARGS[1], ","); file = files[1]; every = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 1
run = get(ENV, "RUN", replace(basename(file), r"\.jls$" => "", r"^snaps_" => ""))
Tcmb = parse(Float64, get(ENV, "TCMB", "54.5")); z = parse(Float64, get(ENV, "Z", "19"))
lyamode = get(ENV, "LYA", "none"); neq_on = get(ENV, "NEQ", "0") == "1"; spec_on = get(ENV, "SPEC", "0") == "1"; rays_on = get(ENV, "RAYS", "1") == "1"
KH2 = parse(Float64, get(ENV, "KH2", "0")); radio = (parse(Float64, get(ENV, "RADIO_AR", "0")), parse(Float64, get(ENV, "RADIO_BETA", "-2.6")))
Perf.reset!("s8b_" * run; stage = "hi21_snap", file = basename(file), lya = lyamode, neq = neq_on)

snaps = Any[]; jlabel = String[]       # jlabel: exchange index of an env_snap_<j>.jls file (for the star record of LYA=rt)
for f in files
    ss = deserialize(f); m = match(r"env_snap_(\w+)\.jls", basename(f))
    append!(snaps, ss); append!(jlabel, fill(m === nothing ? "" : m.captures[1], length(ss)))
end
@printf("%s: %d snapshots; T_CMB = %.2f K, z = %g, Lyα = %s\n", basename(file), length(snaps), Tcmb, z, lyamode)

# ---- Lyα providers
const JC = 27 * H.A10 * H.Tgamma(Tcmb, z; Ar = radio[1], β = radio[2]) / (16π * H.SIGINT_A * H.TSTAR)      # critical J_α [photons cm⁻² s⁻¹ Hz⁻¹ sr⁻¹]: x_α = 1 at S = 1
function lya_for(s, k)
    n = length(s.T)
    if lyamode == "none"; return nothing
    elseif startswith(lyamode, "scan:")
        f = parse(Float64, split(lyamode, ":")[2]); return H.lya_from_Jnu(fill(f * JC * H.hP * H.NU_A, n))
    elseif lyamode == "rt"
        return get(ENV, "STARDIR", "") == "" ? error("LYA=rt needs STARDIR and snapshot names env_snap_<j>.jls") : rt_lya(s, jlabel[k])
    else
        error("unknown LYA = $lyamode")
    end
end

# RadTransfer fallback for the Lyα input: the star's spectrum (S7b star_source, inner free-fall column) transferred through the snapshot's envelope (no re-ionisation sweep: the snapshot's
# state is the one the run's sweep produced); J at the bin containing Lyα: the continuum near the line, not the trapped line intensity (that belongs to the Monte Carlo transfer, LyaMC).
if lyamode == "rt"
    for f in ("PrimordialChem", "ContinuumOpacity", "ChemNet", "ChemRates", "AtomicHHe", "RadTransfer", "s7b_exchange", "s7b_feedback"); include(joinpath(J_, f * ".jl")); end
end
function read_rec(path)
    d = Dict{String,Float64}()
    for l in eachline(path); m = match(r"^(\w+) = (\S+)", l); m === nothing || (d[m.captures[1]] = parse(Float64, m.captures[2])); end
    d
end
function rt_lya(s, jj)
    RT = Main.RadTransfer; CN = Main.ChemNet; S = Main.S7b; AH = Main.AtomicHHe
    rec = read_rec(joinpath(ENV["STARDIR"], "star_$(jj).txt"))
    g = RT.default_grid(); net, photos = AH.atomic(); ab = RT.absorbers(photos, g); ia = RT.ion_atoms(photos, g; dr = AH.k_di_HeII)
    src = S.star_source(g, rec; mode = :star)
    inner = S.inner_column(ia, src, rec["M"], rec["R"], s.r[1], s.ρ[1], max(abs(s.u[1]), 1e3); nz = 300, T = 1.5e4, reemit = (g, max(rec["Tph"], 3000.0)))
    n = length(s.T); nH = s.ρ ./ CN.μH(); X = Dict{Symbol,Vector{Float64}}()
    for (nm, i) in ((:Hp, 1), (:Hm, 2), (:H2, 3), (:H2p, 4)); X[nm] = [s.y[j][i] for j in 1:n]; end
    X[:H] = [max(1 - s.y[j][1] - s.y[j][2] - 2s.y[j][3] - 2s.y[j][4], 0.0) for j in 1:n]
    X[:e] = [max(s.y[j][1] + s.y[j][4] - s.y[j][2], 0.0) for j in 1:n]
    yHe = 1 / 12; X[:He] = fill(yHe, n); X[:Hep] = zeros(n); X[:Hepp] = zeros(n)
    med = RT.Medium(collect(float.(s.r)), collect(float.(s.T)), nH, X)
    rf = RT.transfer(g, ab, RT.Source(inner.Lout, 0.0, sum(inner.Lout)), med, RT.Opts(diffuse = false); rho = s.ρ)
    fa = argmin(abs.(g.ν .- H.NU_A)); @printf("   RT Lyα bin %d: ν = %.4e Hz (Lyα %.4e); star M = %.1f, L* = %.2e L☉, Q_ion* = %.2e s⁻¹, escaping the inner column %.2e\n",
                                           fa, g.ν[fa], H.NU_A, rec["M"] / 1.98847e33, rec["Lstar"] / 3.828e33, inner.Qstar, inner.Qesc)
    H.lya_from_Jnu(rf.J[:, fa])
end

mkpath(joinpath(ROOT, "lib", "popiii"))
fsum = open(joinpath(ROOT, "lib", "popiii", "hi21_$(run).csv"), "w"); println(fsum, H.summary_header())
fprof = open(joinpath(ROOT, "lib", "popiii", "hi21_$(run)_prof.csv"), "w"); println(fprof, H.PROFILE_HEADER)
yprev = Ref{Any}(nothing); tprev = Ref(0.0)
for (k, s) in enumerate(snaps)
    (k % every == 0 || k == length(snaps)) || continue
    lya = Perf.@section :lya lya_for(s, k)
    dt = yprev[] === nothing ? 0.0 : s.t - tprev[]
    o_, sp, sm = Perf.@section :hi21 H.analyze(s; Tcmb = Tcmb, z = z, lya = lya, neq = neq_on && yprev[] !== nothing ? yprev[][end-length(s.T)+1:end] : nothing, dt = dt, rays = rays_on, radio = radio, KH2 = KH2)
    yprev[] = o_.y; tprev[] = s.t
    label = "$(run)_$(k)"
    println(fsum, H.summary_row(label, s.step, s.t, sm)); foreach(l -> println(fprof, l), H.profile_rows(label, s, o_, sp)); flush(fsum); flush(fprof)
    spec_on && sp !== nothing && open(joinpath(ROOT, "lib", "popiii", "hi21_$(run)_spec_$(k).csv"), "w") do io
        hdr, rows = H.spectrum_rows(sp); println(io, hdr); foreach(l -> println(io, l), rows)
    end
    Perf.count!(:snapshots)
    @printf("%3d step %6d  log n_c %6.2f T_c %8.0f  M_HI %9.2f Msun  <Ts>_m %8.1f K (<Tk>_m %8.1f)  max tau21 %9.3e  dTb(disc) %9.3f .. %8.3f mK  centre %9.3f .. %8.3f  EW %9.3e  t_eq/t_ff max %8.1e  theta %.2e\"\n",
            k, s.step, sm.ncen, sm.Tcen, sm.MHI, sm.Ts_mw, sm.Tk_mw, sm.tau_max, sm.dTb_disc_min, sm.dTb_disc_max, sm.dTb_cen_min, sm.dTb_cen_max, sm.EW, sm.teq_max_over_tff, sm.theta_arcsec)
end
close(fsum); close(fprof)
Perf.write()
