# Full emergent spectra (continuum + H I / He I / He II / H₂ lines, observer frame) of recorded snapshots: module SpecFull.
# Post-processing; nothing here feeds back into a run.
#
#   julia -t 8 --project=<env> snapfull.jl <run dir> [stepmod]        # cache <run dir>/specfull/<step>_<t>.jls, export/specfull.{json,f32}
#   julia -t 8 --project=<env> snapfull.jl --one <run dir> [step|last]  # one snapshot, diagnostics printed
#   julia -t 8 --project=<env> snapfull.jl --test                       # LTE slab: no line features
#
# SPECMAX=<n>: stop after n new spectra. The selection is coarse-to-fine as snapspec.jl's (stepmod default 640: ~4 × fewer than
# snapspec, the line pass costs ~1 min per snapshot).
using Serialization, Printf
include(joinpath(@__DIR__, "snapspec.jl"))                 # runrec, Z0Physics, Spectrum, ContinuumOpacity, H2Lines, H2Levels, SnapSpec
isdefined(Main, :ModelAtom) || include(joinpath(@__DIR__, "ModelAtom.jl"))
isdefined(Main, :SpecFull) || include(joinpath(@__DIR__, "SpecFull.jl"))
const SF = SpecFull; const SS = SnapSpec

"Full spectrum of snapshot x, with the run's ALI H₂ populations (snapspec h2lines) for the H₂ lines."
function full_spectrum(x)
    τ = SS.greytau(x); kph = findlast(>(2 / 3), τ); d = SS.densities(x)
    h2 = try SS.h2lines(x, τ, kph) catch e; @warn "H₂ line transfer failed" e; nothing end
    HL = SS.HLV[]                                              # H2Levels: the 63 levels of the ALI populations
    P = h2 === nothing ? nothing : H2Lines.LAST[].pops
    ok = P !== nothing && HL !== nothing && size(P, 1) == length(HL.E)
    jc = h2 === nothing ? 0 : h2.jc; nc = jc >= 2 ? 1 : 0
    h2pop(z) = (ok && z > jc && z - jc + nc <= size(P, 2)) ? view(P, :, z - jc + nc) : nothing
    s = SF.spectrum(x, d, τ, kph; h2pop, HL)
    merge(s, (; h2ok = ok, Lrun_lines = x.Llines))
end

const LSUN = 3.83e33; const RSUN = 6.957e10
function report(s)
    q(Lν) = (ν = s.ν; sum(0.5 * (Lν[i] + Lν[i+1]) * (ν[i+1] - ν[i]) for i in 1:length(ν)-1 if isfinite(Lν[i]) && isfinite(Lν[i+1]); init = 0.0))
    @printf("step %d: %d points (cut zone %d, %d refined LTE sub-points), H₂ populations %s\n", s.step, s.np, s.kc, s.nrefined, s.h2ok ? "ALI" : "none")
    isfinite(s.Ltop) && @printf("  refined photosphere: run L at the zone top (r %.1f R☉) %.4e L☉; T⁴ factor fT %.4g (iterates %s)\n", s.rtop / RSUN, s.Ltop / LSUN, s.fT, join(round.(s.fhist, sigdigits = 4), ", "))
    @printf("  times [s]: %s\n", join(["$k $(round(v, digits = 1))" for (k, v) in pairs(s.times)], ", "))
    @printf("  ∫L_ν: edge %.4e L☉ (run %.4e), photosphere %.4e L☉ (run %.4e); T_eff %.0f K, r_ph %.1f R☉\n",
            q(s.cont.Le) / LSUN, s.Lrun_edge / LSUN, q(s.cont.Lp) / LSUN, s.Lrun_ph / LSUN, s.Teff, s.rph / RSUN)
    @printf("  %-14s %12s %12s %12s %10s %10s\n", "line", "λ0 [Å]", "L_line [L☉]", "EW [Å]", "τ0(r_in)", "τ0(mid)")
    for (j, ln) in enumerate(s.lines)
        @printf("  %-14s %12.2f %12.4e %12.4g %10.3e %10.3e\n", ln.name, 1e8 * SF.c / ln.ν0, s.Lline[j] / LSUN, s.EW[j], s.τr[j, 1], s.τr[j, SF.NRT ÷ 2])
    end
end

"""
LTE test: an isothermal slab-like sphere (T, ρ uniform, LTE composition, J = B for the populations, the static case and with
infall): the line windows must show no features (L/L_c = 1 to the solver error), since S = B everywhere.
"""
function lte_test()
    T = 8000.0; ρ = 1e-10; n = 40; r = collect(exp.(range(log(1e12), log(3e12), length = n + 1)))
    s = Z0Physics.opacity_state(ρ, T)
    pts = Any[]
    for v0 in (0.0, -1e7)
        empty!(pts)
        for k in 1:n
            rc = 0.5 * (r[k] + r[k+1])
            push!(pts, (; r = rc, ρ, T, v = v0 * sqrt(r[1] / rc), dvdr = 0.0, zone = k, Δr = r[k+1] - r[k], τg = 100.0, lte = true,
                        nH0 = s.nH0, nHp = s.nHp, nHm = s.pH[1] * s.ne * exp(Z0Physics.χHm / (Z0Physics.kB * T)) / (4 * Z0Physics.λe3(T)),
                        nH2 = 0.0, ne = s.ne, nHe0 = s.nHe0, nHe1 = s.nHe1, nHe2 = s.nHe2, nHD = 0.0))
        end
        atm = (r = [p.r for p in pts], ρ = fill(ρ, n), T = fill(T, n), zone = collect(1:n), jump = falses(n - 1))
        _, R = SF.SP.rays(atm.r; ncore = 8); R = SF.SP.weights!(atm, R)
        pop = SF.populations(pts, k -> (ν -> SF.MA.Bν(ν, T)))
        op = SF.FOpac(pts, pop); LL = SF.linelist(nothing)
        for j in (1, 3, 4, 14, 21, 23)
            j > length(LL) && continue
            ln = LL[j]; ν0 = ln.ν0; νo = ν0 .* (1 .- collect(range(-3e7, 3e7, length = 41)) ./ SF.c)
            χc = zeros(n); ηc = zeros(n); K = zeros(n, 1); E = zeros(n, 1); pr = [Union{Nothing,Z0Physics.LineProf}[nothing] for _ in 1:n]
            for k in 1:n
                ab, sc, St = SF.SP.coeffs3(op, k, ν0); χc[k] = ab + sc; ηc[k] = ab * St + sc * SF.MA.Bν(ν0, T)
                K[k, 1], E[k, 1] = SF.strength(ln, k, pts, pop, z -> nothing); pr[k][1] = SF.profile(ln, pts[k])
            end
            L = SF.window_L(atm, R, pts, νo, χc, ηc, K, E, pr, [j], LL); Lc = SF.window_L(atm, R, pts, [ν0], χc, ηc, K, E, pr, [j], LL; lines = false)[1]
            rmid = atm.r[n ÷ 2] * 1.0001                                                # the photosphere-sphere output at mid-slab
            _, LP = SF.window_L(atm, R, pts, νo, χc, ηc, K, E, pr, [j], LL; rph = rmid); _, LcP = SF.window_L(atm, R, pts, [ν0], χc, ηc, K, E, pr, [j], LL; lines = false, rph = rmid)
            Sl = E[n ÷ 2, 1] / K[n ÷ 2, 1]
            @printf("v0 %6.0f km/s %-12s S_line/B %.6f  edge max|L/Lc − 1| %.2e  at r_mid max|L − Lc|/(4π²r²B) %.2e  (τ0 centre %.2e)\n", v0 / 1e5, ln.name, Sl / SF.MA.Bν(ν0, T),
                    maximum(abs.(L ./ Lc .- 1)), maximum(abs.(LP .- LcP[1])) / (4π^2 * rmid^2 * SF.MA.Bν(ν0, T)), K[n ÷ 2, 1] * Z0Physics.prof_eval(pr[n ÷ 2][1], 0.0) * (r[end] - r[1]))
        end
    end
end

const SPSCALF = ["step", "t_yr", "r_ph", "T_eff", "Lrun_edge", "Lrun_ph", "Lc_edge", "Lc_ph", "np", "nrefined"]
"""
Export: specfull.json (grids, line names, per spectrum the scalars SPSCALF) and specfull.f32, per spectrum: νL_ν/L☉ edge (nν),
νL_ν/L☉ photosphere (nν), r(τ_ν = ⅔)/R☉ (nν); per line (nl): L_ν/L_c on the velocity grid (nv), then per line L_line/L☉, EW [Å],
L_c(ν0) ν0/L☉ (3 nl); then r grid/R☉ (NRT) and log10 τ0 per line (nl × NRT); ver 2 ("ver":2 in the json) appends the same for the
light through the photosphere sphere r_ph: per line L_ν/L_c on the velocity grid (nv), then per line L_line/L☉, EW [Å], L_c(ν0) ν0/L☉
(3 nl), all NaN for spectra cached before ver 2.
"""
function export_full(L, cache, out)
    mkpath(out); rows = Any[]; meta = nothing
    open(joinpath(out, "specfull.f32.tmp"), "w") do io
        for x in L
            f = joinpath(cache, SS.key(x) * ".jls"); isfile(f) || continue
            s = try deserialize(f) catch; continue end
            meta === nothing && (meta = (ν = s.ν, names = [l.name for l in s.lines], λ0 = [1e8 * SF.c / l.ν0 for l in s.lines], v = s.vgrid))
            length(s.lines) == length(meta.names) || continue
            q(Lν) = sum(0.5 * (Lν[i] + Lν[i+1]) * (s.ν[i+1] - s.ν[i]) for i in 1:length(s.ν)-1 if isfinite(Lν[i]) && isfinite(Lν[i+1]); init = 0.0)
            write(io, Float32.(s.cont.Le .* s.ν ./ LSUN)); write(io, Float32.(s.cont.Lp .* s.ν ./ LSUN)); write(io, Float32.(s.cont.rτ ./ RSUN))
            for j in eachindex(s.lines); write(io, Float32.(s.Lprof[j, :] ./ max(s.Lcont[j], 1e-300))); end
            for j in eachindex(s.lines); write(io, Float32[s.Lline[j] / LSUN, s.EW[j], s.Lcont[j] * s.lines[j].ν0 / LSUN]); end
            write(io, Float32.(s.rgrid ./ RSUN)); for j in eachindex(s.lines); write(io, Float32.(log10.(max.(s.τr[j, :], 1e-30)))); end
            # ver 2: the lines through the photosphere sphere r_ph (NaN for spectra cached before they were computed)
            hasP = hasproperty(s, :LprofP); nv_ = length(s.vgrid)
            for j in eachindex(s.lines); write(io, hasP ? Float32.(s.LprofP[j, :] ./ (abs(s.LcontP[j]) > 0 ? s.LcontP[j] : NaN)) : fill(Float32(NaN), nv_)); end
            for j in eachindex(s.lines); write(io, hasP ? Float32[s.LlineP[j] / LSUN, s.EWP[j], s.LcontP[j] * s.lines[j].ν0 / LSUN] : fill(Float32(NaN), 3)); end
            push!(rows, [s.step, s.t / 3.15576e7, s.rph / RSUN, s.Teff, s.Lrun_edge / LSUN, s.Lrun_ph / LSUN, q(s.cont.Le) / LSUN, q(s.cont.Lp) / LSUN, s.np, s.nrefined])
        end
    end
    meta === nothing && return println("no full spectra in $cache")
    js(v) = isfinite(v) ? @sprintf("%.6g", v) : "null"
    open(joinpath(out, "specfull.json.tmp"), "w") do io
        print(io, "{\"ver\":2,\"nnu\":", length(meta.ν), ",\"nl\":", length(meta.names), ",\"nv\":", length(meta.v), ",\"nrt\":", SF.NRT,
              ",\"lam\":[", join((@sprintf("%.6g", 1e4 * SF.c / ν) for ν in meta.ν), ","), "],\"lines\":[", join(("\"$n\"" for n in meta.names), ","),
              "],\"lam0\":[", join((@sprintf("%.6g", v) for v in meta.λ0), ","), "],\"v\":[", join((@sprintf("%.4g", v / 1e5) for v in meta.v), ","),
              "],\"scalars\":[", join(("\"$s\"" for s in SPSCALF), ","), "],\"spectra\":[", join(("[" * join(js.(r), ",") * "]" for r in rows), ","), "]}")
    end
    for f in ("specfull.f32", "specfull.json"); mv(joinpath(out, f * ".tmp"), joinpath(out, f); force = true); end
    println("exported $(length(rows)) full spectra to $out")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if ARGS[1] == "--test"
        lte_test()
    elseif ARGS[1] == "--one"
        L = lineage(ARGS[2]); w = length(ARGS) >= 3 ? ARGS[3] : "last"
        x = w == "last" ? L[end] : L[findlast(y -> y.step <= parse(Int, w), L)]
        s = full_spectrum(x); report(s)
        length(ARGS) >= 4 && serialize(ARGS[4], s)
    else
        dir = ARGS[1]; L = lineage(dir); smod = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 640
        cache = joinpath(dir, "specfull"); mkpath(cache)
        sel = [i for i in eachindex(L) if L[i].step % smod == 0]; (isempty(sel) || sel[end] != length(L)) && !isempty(L) && push!(sel, length(L))
        order = isempty(sel) ? Int[] : [sel[end]]                     # the newest first (monitoring), then coarse to fine
        for st in (64, 32, 16, 8, 4, 2, 1), j in 1:st:length(sel); sel[j] in order || push!(order, sel[j]); end
        nmax = parse(Int, get(ENV, "SPECMAX", "100000")); ndone = 0
        for i in order
            ndone >= nmax && break
            f = joinpath(cache, SS.key(L[i]) * ".jls"); isfile(f) && continue
            s = try full_spectrum(L[i]) catch e; @warn "full spectrum failed" step = L[i].step e; continue end
            tmp = f * ".tmp"; serialize(tmp, s); mv(tmp, f; force = true)
            @printf("step %d: %s\n", s.step, join(["$k $(round(v, digits = 1))" for (k, v) in pairs(s.times)], ", ")); flush(stdout)
            global ndone += 1
            ndone % 5 == 0 && export_full(L[sort(sel)], cache, joinpath(dir, "export"))
        end
        export_full(L[sel], cache, joinpath(dir, "export"))
    end
end
