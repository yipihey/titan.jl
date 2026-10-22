# S7b envelope side of the coupled protostellar feedback: hooks for popiii_m2c.jl (SINK=1 collapse driver), activated with S7B=1. Everything is opt-in.
#
# Scheme (see docs/POPIII_NOTES.md, S7b): two processes with a common clock t_c = 0 at the synchronisation point. This process (the envelope) is the master of the clock: it chooses the next
# barrier t_{j+1} = t_j + Δτ_j (Δτ = FEX × accretion time, FEX default 0.01), runs its hydrodynamics to it (the time step is cut to land on it), then
#   (1) waits for the star's record at t_j (star_<j>.txt, written by s7b_star.jl in the protostar driver),
#   (2) computes the radiation field of that star on the envelope: inner free-fall column between R* and the grid edge → spectrum entering the grid → equilibrium sweep + transfer on the grid
#       (S7b.envfields): photo-rates J_ν per zone for the ChemNet step (ionisation, H⁻ detachment, LW dissociation, photo-heating), radiation acceleration for the momentum row (TitanResidFast.AEXT),
#   (3) writes env_<j>.txt: the accretion rate through the sink over the last window, the accreted mass, the next barrier, the infall state at the grid edge.
# The star advances in parallel over the same interval with a constant Ṁ predicted from the envelope's last window, corrected for the mass difference (S7b star hook).
# One-way mode (S7B_TRACK = S6b csv): the star is read from a stored track as a function of the accreted mass instead of a partner process (no exchange).
#
# Environment: S7B_DIR (exchange directory), S7B_TAG, S7B_MSYNC [Msun, accretion since the start of the run at which the coupling starts], S7B_FEX [0.01], S7B_MODE [star|photosphere],
# S7B_RHOMODE [mdot|grid], S7B_TINNER [1.5e4], S7B_NZIN [300], S7B_DIFFUSE [0], S7B_GRID [coarse|default], S7B_LW [1], S7B_FORCE [1], S7B_FEEDBACK [1] (0: transfer evaluated, nothing applied),
# S7B_MCAP [Msun], S7B_MDOTSTOP [Msun/yr], S7B_MDOTWIN [0.3] (window in ln M of the accretion-rate estimate), S7B_TRACK, S7B_SAVEAT [list of Msun: snapshot of the envelope], S7B_RESUME.
module S7bEnv
using Printf, Serialization
const S = Main.S7b; const RT = Main.RadTransfer; const CN = Main.ChemNet
envf(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
envi(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
const Msun = S.Msun; const YR = S.YR; const AU = S.AU; const Rsun = S.Rsun; const G = S.G
const TAG = get(ENV, "S7B_TAG", "run"); const DIR = get(ENV, "S7B_DIR", joinpath(Main.TF.ROOT, "lib", "popiii", "s7b_x_$(TAG)"))
const FEX = envf("S7B_FEX", 0.01); const MSYNC = envf("S7B_MSYNC", 10.0) * Msun; const MODE = Symbol(get(ENV, "S7B_MODE", "star")); const RHOMODE = get(ENV, "S7B_RHOMODE", "mdot")
const TINNER = envf("S7B_TINNER", 1.5e4); const NZIN = envi("S7B_NZIN", 300); const DIFFUSE = get(ENV, "S7B_DIFFUSE", "0") == "1"; const USELW = get(ENV, "S7B_LW", "1") == "1"
const HEATFAC = envf("S7B_HEAT", 1.0); const HEATMODE = get(ENV, "S7B_HEATMODE", "abs"); const FORCE = get(ENV, "S7B_FORCE", "1") == "1"; const FEEDBACK = get(ENV, "S7B_FEEDBACK", "1") == "1"; const FORCEFAC = envf("S7B_FORCEFAC", 1.0)
const MCAP = envf("S7B_MCAP", 1e9) * Msun; const MDOTSTOP = envf("S7B_MDOTSTOP", 0.0) * Msun / YR; const MDOTWIN = envf("S7B_MDOTWIN", 0.3)
const NPROF = envi("S7B_NPROF", 0)
const TRACK = get(ENV, "S7B_TRACK", ""); const EXCH = TRACK == "" && get(ENV, "S7B_EXCH", "1") == "1"
const SAVEAT = haskey(ENV, "S7B_SAVEAT") ? sort(parse.(Float64, split(ENV["S7B_SAVEAT"], ","))) .* Msun : Float64[]
const PERFD = Dict{Symbol,Tuple{Int,Float64}}()
function perf(name, f)
    t0 = time_ns(); v = f(); cs = get(PERFD, name, (0, 0.0)); PERFD[name] = (cs[1] + 1, cs[2] + (time_ns() - t0) * 1e-9); v
end
const GRID = Ref{Any}(nothing); const IA = Ref{Any}(nothing); const AB = Ref{Any}(nothing); const ZE = Ref{Any}(nothing); const NET = Ref{Any}(nothing); const PHOTOS = Ref{Any}(nothing)
const JCMB = Ref(Float64[])

Base.@kwdef mutable struct State
    phase::Symbol = :pre                 # :pre | :coupled | :done
    events::Vector{Tuple{Float64,Float64}} = Tuple{Float64,Float64}[]      # (crossing time [s], interface mass [g]) of the interfaces removed by the sink
    last_a::Int = 0; prev_r::Vector{Float64} = Float64[]; prev_t::Float64 = 0.0
    j::Int = 0; t0::Float64 = 0.0; Macc0::Float64 = 0.0; tb::Float64 = 0.0; dtau::Float64 = 0.0; Mstar0::Float64 = 0.0
    Mdot::Float64 = 0.0; Macc::Float64 = 0.0
    a0::Int = 0; J::Matrix{Float64} = zeros(0, 0); kLW::Vector{Float64} = Float64[]; heat::Vector{Float64} = Float64[]; active::Bool = false; heatnow::Float64 = 0.0   # fields in force
    log::Any = nothing; prof::Any = nothing; track::Any = nothing; texch::Float64 = 0.0
end
const ST = State()

# ------------------------------------------------------------------ set-up (called from popiii_m2c.jl while the network is being built)
function setup_grid()
    GRID[] === nothing || return
    Main.AtomicHHe.FJ_FLOOR[] = envf("S7B_FJFLOOR", 0.0)
    GRID[] = get(ENV, "S7B_GRID", "coarse") == "default" ? RT.default_grid() : RT.default_grid(dIR = 0.3, dopt = 0.25, dLW = 0.02, dUV = 0.05)
    mkpath(DIR)
end
"Network with the Lyman–Werner reaction (RadTransfer.with_lw) when S7B_LW = 1."
function prepare_net(netw, phot)
    setup_grid(); NET[] = USELW ? RT.with_lw(netw) : netw; PHOTOS[] = phot
    DIFFUSE && (Main.AtomicHHe.ONSPOT[] = false)
    (NET[], phot)
end
"The per-zone radiation environment of ChemNet (RadTransfer.ZoneEnv) with the CMB as the initial field."
function make_env(phot, tcmb)
    setup_grid(); g = GRID[]
    ze = RT.ZoneEnv(g, phot; Trad = tcmb)
    JCMB[] = tcmb > 0 ? [RT.Bν(g.ν[f], tcmb) for f in 1:RT.nfreq(g)] : zeros(RT.nfreq(g))
    RT.set_zone!(ze, JCMB[], 0.0)
    ZE[] = ze; AB[] = RT.absorbers(phot, g); IA[] = RT.ion_atoms(phot, g; dr = Main.AtomicHHe.k_di_HeII)
    ze
end
"Wrap the zone hook (escape probabilities) so that the zone's radiation field is loaded before its chemistry step, and the other photo-heating enters the cooling."
function install!(mat)
    old = mat.zonehook; cool0 = mat.cooling
    mat.zonehook = function (k, d)
        old === nothing || old(k, d)
        j = k - ST.a0 + 1
        if ST.active && 1 <= j <= size(ST.J, 1)
            RT.set_zone_light!(ZE[], view(ST.J, j, :), ST.kLW[j]); ST.heatnow = ST.heat[j]
        else
            RT.set_zone_light!(ZE[], JCMB[], 0.0); ST.heatnow = 0.0
        end
        nothing
    end
    mat.cooling = (ρ, T, y) -> cool0(ρ, T, y) - ST.heatnow          # photo-heating of the continua that the network does not carry (H⁻ ff, …): negative cooling
    ST.log = open(joinpath(Main.TF.ROOT, "lib", "popiii", "s7b_$(TAG)_env.csv"), "w")
    println(ST.log, "j,t_yr,tc_yr,dtau_yr,Macc,Mpt,Mdot_win,Mstar,Rstar_Rsun,Lstar,Linf,Tph,Teffstar,Qstar,Qesc_in,RII_over_R,Qcrit,Q_over_Qcrit,Labs_in_frac,rho_s,u_s_kms,r_edge_AU,Q_in_grid_abs,n_lit,rfront_AU,xe_max,kLW_max,kLW_zone1,arad_max,arad_over_g_max,Fsum_dyn,Fsum_over_Mdotv,Tzone1,fH2_zone1,wall_s")
    if NPROF > 0
        ST.prof = open(joinpath(Main.TF.ROOT, "lib", "popiii", "s7b_$(TAG)_prof.csv"), "w")
        println(ST.prof, "j,t_yr,zone,r_AU,rho,T,u_kms,p,fH2,xHp,arad,kLW,g,x_sweep,m_int")
    end
    if TRACK != ""
        rows = [parse.(Float64, split(l, ",")) for l in readlines(TRACK)[2:end] if !isempty(strip(l))]
        ST.track = (M = [r[3] for r in rows], R = [r[4] for r in rows], Ls = [r[5] for r in rows], Lacc = [r[6] for r in rows], Tph = [r[7] for r in rows], Linf = [r[23] for r in rows])
    end
    nothing
end

# ------------------------------------------------------------------ accretion bookkeeping
rsink(M) = Main.rsink(M)
"Record the interfaces that the sink removed in the last step (crossing time interpolated inside the step, as the sink driver's events)."
function track_events!(t)
    TF = Main.TF; a = Int(TF.getscalar(:ngrs)); rn = TF.var(:rn); xm = Main.TX.XMINT[]
    if ST.last_a == 0
        ST.last_a = a; ST.prev_r = copy(rn); ST.prev_t = t; return
    end
    for i in ST.last_a:a-1
        m = 4π * xm[i]; rs = rsink(m); rp = ST.prev_r[i]
        fr = rp > rn[i] ? clamp((rp - rs) / (rp - rn[i]), 0.0, 1.0) : 1.0
        push!(ST.events, (ST.prev_t + fr * (t - ST.prev_t), m))
    end
    ST.last_a = a; ST.prev_r = copy(rn); ST.prev_t = t
    length(ST.events) > 2000 && deleteat!(ST.events, 1:1000)
end
"""
Accretion rate through the sink and the accreted mass at time t from the sink events (crossing time, interface mass): least-squares line m = m̄ + Ṁ (t' − t̄) through the events of the last window
(ln M span ≥ MDOTWIN, at least 4 events), evaluated at t. The events are a staircase (one zone mass each, 5 % of M, with jitter in the crossing times: single pairs scatter by ×4); the line gives a
smooth rate and a smooth mass whose level, unlike a window-averaged rate integrated in time, has no lag offset. Returns (Ṁ, M(t)) or (NaN, NaN) with fewer than 4 events.
"""
function mdot_window(t)
    ev = ST.events; n = length(ev)
    n < 4 && return (NaN, NaN)
    m1 = ev[n][2]; k = n - 3
    while k > 1 && log(m1 / ev[k][2]) < MDOTWIN; k -= 1; end
    ts = [e[1] for e in ev[k:n]]; ms = [e[2] for e in ev[k:n]]
    tb = sum(ts) / length(ts); mb = sum(ms) / length(ms)
    den = sum((ts .- tb) .^ 2)
    den > 0 || return (NaN, NaN)
    Md = sum((ts .- tb) .* (ms .- mb)) / den
    Md > 0 || return (NaN, NaN)
    # the sink may starve (the infall halts): the line through the staircase would keep growing. Bound the rate by the average since the first event of the window, and the mass by the point mass
    t1, mlast = ev[n]; t0w, m0w = ev[k]
    t > t1 && (Md = min(Md, (mlast - m0w) / max(t - t0w, 1e-30)))
    (Md, min(mb + Md * (t - tb), 4π * Main.TX.XMINT[][Int(Main.TF.getscalar(:ngrs))]))
end

# ------------------------------------------------------------------ the field on the grid
"Star record → (spectrum entering the grid, diagnostics); the star record `rec` has the keys of s7b_star.jl."
function field!(rec, Mdot_dyn)
    TF = Main.TF; TX = Main.TX; a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); n = b - a + 1
    rn = TF.var(:rn)[a:b+1]; dn = TF.var(:dn)[a:b]; tn = TF.var(:tn)[a:b]; un = TF.var(:un)[a:b+1]; dv = TF.var(:dvoln)[a:b]
    g = GRID[]; ia = IA[]
    src = S.star_source(g, rec; mode = MODE)
    Mstar = rec["M"]
    redge = rn[1]; u_edge = abs(un[1])
    ρs = RHOMODE == "mdot" && isfinite(Mdot_dyn) && Mdot_dyn > 0 ? Mdot_dyn / (4π * redge^2 * max(u_edge, 1e3)) : dn[1]
    inner = perf(:rt_inner, () -> S.inner_column(ia, src, Mstar, rec["R"], redge, ρs, max(u_edge, 1e3); nz = NZIN, T = TINNER, reemit = (g, max(rec["Tph"], 3000.0))))
    mz = 4π .* dn .* dv
    rmid = 0.5 .* (rn[1:end-1] .+ rn[2:end]); ubar = 0.5 .* (abs.(un[1:end-1]) .+ abs.(un[2:end]))
    mdotH = 4π .* rmid .^ 2 .* dn .* ubar ./ CN.μH()
    xnet = [Main.mat.x[a+j-1] for j in 1:n]
    opts = RT.Opts(diffuse = DIFFUSE)
    f = S.envfields(g, ia, AB[], NET[], PHOTOS[], rn, tn, dn, mdotH, xnet, inner.Lout; opts = opts, Jcmb = JCMB[], perf = perf)
    # apply
    if FEEDBACK
        ST.a0 = a; ST.J = f.J; ST.kLW = copy(f.kLW); ST.heat = (HEATMODE == "net" ? f.heat_net : f.heat_other) .* HEATFAC; ST.active = true
        if FORCE
            acc = RT.interface_accel(f.arad .* FORCEFAC, mz); full = zeros(length(un_full())); full[a:b+1] .= acc; TX.AEXT[] = full
        end
    end
    # diagnostics
    x_sw = f.sw === nothing ? zeros(n) : f.sw.x
    kf = findlast(>=(0.5), x_sw)
    vff = sqrt.(2 * G * (4π .* TX.XMINT[][a:b]) ./ rmid)
    gloc = G .* (4π .* TX.XMINT[][a:b]) ./ rmid .^ 2
    Fsum = sum(f.arad .* mz)
    Nabs = f.sw === nothing ? 0.0 : f.sw.ledger.Nion
    (; inner, src, f, ρs, redge, u_edge, n, kf, rn, x_sw, Nabs, Fsum, Fmax = maximum(f.arad ./ gloc), amax = maximum(abs.(f.arad)), mz)
end
un_full() = Main.TF.var(:un)

# ------------------------------------------------------------------ the barrier
function star_record_track(Macc)
    tr = ST.track; M = (ST.Mstar0 + Macc) / Msun
    M = clamp(M, tr.M[1], tr.M[end]); j = clamp(searchsortedlast(tr.M, M), 1, length(tr.M) - 1); w = (M - tr.M[j]) / (tr.M[j+1] - tr.M[j])
    lerp(v) = (1 - w) * v[j] + w * v[j+1]
    Dict("M" => M * Msun, "R" => lerp(tr.R) * Rsun, "Lstar" => lerp(tr.Ls) * S.Lsun, "Lacc" => lerp(tr.Lacc) * S.Lsun, "Tph" => lerp(tr.Tph), "Linf" => max(lerp(tr.Linf), lerp(tr.Ls)) * S.Lsun, "t" => 0.0)
end
fname(kind, j) = joinpath(DIR, @sprintf("%s_%05d.txt", kind, j))

function barrier!(t)
    TF = Main.TF; t_w0 = time()
    j = ST.j; tc = t - ST.t0
    Md, Macc_t = mdot_window(t)
    if !isfinite(Md); Md = ST.Mdot; Macc_t = ST.Macc0 + ST.Mdot * tc; end
    Macc = Macc_t - ST.Macc0; ST.Macc = Macc; ST.Mdot = Md
    Mpt = 4π * Main.TX.XMINT[][Int(TF.getscalar(:ngrs))]
    # star record at t_j
    rec = nothing
    if TRACK != ""
        rec = star_record_track(Macc)
    else
        r = S.wait_rec(fname("star", j)); rec = Dict(r)
    end
    j == 0 && get(ENV, "S7B_SAVESYNC", "0") == "1" && save_snapshot!(0)       # restart point of the coupled run (with the star's checkpoint at its sync)
    fd = perf(:rt_total, () -> field!(rec, Md))
    t_rt = time() - t_w0
    # next barrier
    tacc = max(rec["M"], Mpt) / max(Md, 1e-30)
    dtau = FEX * tacc
    ST.dtau = dtau; ST.tb = t + dtau
    # records
    Rr = Main.TF.var(:rn)
    row = (j, t / YR, tc / YR, dtau / YR, Macc / Msun, Mpt / Msun, Md * YR / Msun, rec["M"] / Msun, rec["R"] / Rsun, rec["Lstar"] / S.Lsun, rec["Linf"] / S.Lsun, rec["Tph"],
           (max(rec["Lstar"], 1e-3 * S.Lsun) / (4π * rec["R"]^2 * S.σSB))^0.25, fd.inner.Qstar, fd.inner.Qesc, fd.inner.RII / rec["R"], fd.inner.Qcrit, fd.inner.Qstar / fd.inner.Qcrit,
           fd.inner.Labs / max(sum(fd.src.L), 1e-300), fd.ρs, fd.u_edge / 1e5, fd.redge / AU, fd.Nabs, count(>=(0.5), fd.x_sw), fd.kf === nothing ? 0.0 : fd.rn[fd.kf+1] / AU,
           maximum(get(fd.f.med.X, :e, [0.0])), maximum(fd.f.kLW), fd.f.kLW[1], fd.amax, fd.Fmax, fd.Fsum, abs(fd.Fsum) / max(Md * sqrt(2 * G * rec["M"] / rec["R"]), 1e-300),
           Main.TF.var(:tn)[Int(TF.getscalar(:ngrs))], 2 * Main.ycomp(Int(TF.getscalar(:ngrs)))[3], time() - Main.T0WALL)
    println(ST.log, join(row, ",")); flush(ST.log)
    if NPROF > 0
        io = ST.prof; a_ = Int(TF.getscalar(:ngrs)); rn_ = TF.var(:rn); dn_ = TF.var(:dn); tn_ = TF.var(:tn); un_ = TF.var(:un); pg_ = TF.var(:pgn); xm_ = Main.TX.XMINT[]
        for jj in 1:min(NPROF, length(fd.f.arad))
            kk = a_ + jj - 1; yc = Main.ycomp(kk)
            println(io, join((j, t / YR, jj, rn_[kk] / AU, dn_[kk], tn_[kk], un_[kk] / 1e5, pg_[kk], 2yc[3], yc[1], fd.f.arad[jj], fd.f.kLW[jj], G * 4π * xm_[kk] / rn_[kk]^2, fd.x_sw[jj], 4π * xm_[kk] / Msun), ","))
        end
        flush(io)
    end
    @printf("S7B j=%d t_c=%.4f yr  M_acc=%.3f M_pt=%.3f Mdot=%.3e | star M=%.3f R=%.1f L*=%.3e Teff*=%.0f | Q*=%.2e esc_in=%.2e R_II/R*=%.2f Q/Qcrit=%.2e | kLW_max=%.2e a_max=%.2e F/g=%.2e | RT %.1f s\n",
            j, tc / YR, Macc / Msun, Mpt / Msun, Md * YR / Msun, rec["M"] / Msun, rec["R"] / Rsun, rec["Lstar"] / S.Lsun, row[13], fd.inner.Qstar, fd.inner.Qesc, fd.inner.RII / rec["R"], fd.inner.Qstar / fd.inner.Qcrit,
            maximum(fd.f.kLW), fd.amax, fd.Fmax, t_rt)
    flush(stdout)
    # envelope record for the star
    done = Mpt >= MCAP || (MDOTSTOP > 0 && Md < MDOTSTOP && j > 20)
    if EXCH
        save = !isempty(SAVEAT) && Macc / Msun + ST.Mstar0 / Msun >= SAVEAT[1] / Msun
        save && popfirst!(SAVEAT)
        save && save_snapshot!(j)
        S.write_rec(fname("env", j), Dict("j" => j, "t" => t, "tc" => tc, "tnext" => ST.tb, "dtau" => dtau, "Mdot" => Md, "Macc" => Macc, "Macc_pred" => Macc + Md * dtau, "Mpt" => Mpt,
                                          "done" => done ? 1.0 : 0.0, "save" => save ? 1.0 : 0.0, "rho_s" => fd.ρs, "u_s" => fd.u_edge, "r_edge" => fd.redge))
    end
    ST.j += 1; ST.texch += time() - t_w0
    done && (ST.phase = :done)
    nothing
end

"Snapshot of the envelope in the format of popiii_m2c.jl SNAPS (restartable with RESTART=…, RESIDX=0)."
function save_snapshot!(j)
    TF = Main.TF; a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); nsh = b - a + 1; mat = Main.mat
    snap = (; step = 0, t = TF.getscalar(:timen), r = TF.var(:rn)[a:b+1], u = TF.var(:un)[a:b+1], ρ = TF.var(:dn)[a:b], T = TF.var(:tn)[a:b], e = TF.var(:egn)[a:b],
            y = [Main.ycomp(a + i - 1) for i in 1:nsh], m = 4π .* Main.TX.XMINT[][a:b+1], Λext = copy(mat.Λext[a:b]), eq = [Main.eqflag(a + i - 1) for i in 1:nsh])
    serialize(joinpath(DIR, @sprintf("env_snap_%05d.jls", j)), [snap])
end

# ------------------------------------------------------------------ per-step hook
function after_step(k)
    TF = Main.TF; t = TF.getscalar(:timen)
    perf(:s7b_hook, () -> begin
        track_events!(t)
        if ST.phase === :pre
            Md, Macc_t = mdot_window(t)
            if isfinite(Md) && Macc_t >= MSYNC
                # start of the coupling: the clock zero is here; the accreted mass is counted from here
                ST.phase = :coupled; ST.t0 = t; ST.Macc0 = Macc_t; ST.Mdot = Md; ST.Mstar0 = Macc_t
                if EXCH
                    s0 = S.wait_rec(fname("star", 0)); ST.Mstar0 = Macc_t
                end
                @printf("S7B: coupling starts at t = %.6e s, M_acc = %.4f Msun, Mdot = %.3e Msun/yr\n", t, Macc_t / Msun, Md * YR / Msun)
                barrier!(t)
            end
        elseif ST.phase === :coupled
            if t >= ST.tb * (1 - 1e-9)
                barrier!(t)
            end
        end
        if ST.phase === :coupled
            rem = ST.tb - t; dt = TF.getscalar(:dtime)
            if rem < dt * 1.5 && rem > 0
                TF.setscalar!(:dtime, rem < 1.25 * dt && rem > dt ? rem : (rem <= dt ? rem : rem / 2))
            end
        elseif ST.phase === :done
            TF.setscalar!(:tmax, t)
        end
        nothing
    end)
end

function finish()
    ST.log === nothing || close(ST.log)
    if EXCH && ST.phase !== :pre
        S.write_rec(joinpath(DIR, "env_final.txt"), Dict("j" => ST.j, "done" => 1.0))
    end
    # performance record (own file, independent of the driver's Perf instrumentation)
    mkpath(joinpath(Main.TF.ROOT, "lib", "popiii", "perf"))
    open(joinpath(Main.TF.ROOT, "lib", "popiii", "perf", "s7b_$(TAG)_env_sections.csv"), "w") do io
        println(io, "section,calls,total_s,mean_ms")
        for (kk, (c, s)) in sort(collect(PERFD); by = x -> -x[2][2]); @printf(io, "%s,%d,%.4f,%.4f\n", kk, c, s, 1e3 * s / max(c, 1)); end
        @printf(io, "exchange_wait_and_rt_total,%d,%.4f,%.4f\n", ST.j, ST.texch, 1e3 * ST.texch / max(ST.j, 1))
        @printf(io, "wall_total,1,%.4f,%.4f\n", time() - Main.T0WALL, 1e3 * (time() - Main.T0WALL))
    end
    nothing
end
end # module
