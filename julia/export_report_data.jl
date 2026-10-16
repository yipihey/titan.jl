# Data files for the report's interactive figures (report/app/data/*.json).
#   julia --project=env export_report_data.jl bench      benchmark profiles: Fortran 1995-code final dumps vs pure-Julia runs
#   julia --project=env export_report_data.jl qstar NAME[+CONT] [PROFILE.csv] [OUT]   a quasi-star run (lib/qs_runs/NAME_relax_snaps.jls,
#                                                     with the snapshots of its continuation CONT appended)
using Serialization, Printf
include(joinpath(@__DIR__, "Metrics.jl")); using .Metrics
const OUT = joinpath(@__DIR__, "..", "report", "app", "data"); mkpath(OUT)

js(x::Real) = isfinite(x) ? @sprintf("%.6g", x) : "null"
js(x::Integer) = string(x)
js(x::AbstractString) = "\"" * x * "\""
js(v::AbstractVector) = "[" * join(js.(v), ",") * "]"
js(d::NamedTuple) = "{" * join(["\"$k\":" * js(v) for (k, v) in pairs(d)], ",") * "}"
writejs(name, x) = (open(io -> print(io, js(x)), joinpath(OUT, name), "w"); println("wrote ", joinpath(OUT, name)))

"Fortran final dump CSV → named columns (zone values on rows 1:n-1, interface values on 1:n)."
function readdump(path)
    lines = filter(l -> !startswith(l, "#"), readlines(path))
    hdr = split(lines[1], ","); rows = [parse.(Float64, split(l, ",")) for l in lines[2:end]]
    col(n) = [r[findfirst(==(n), hdr)] for r in rows]
    t = parse(Float64, match(r"timen=([-0-9.eE+]+)", readline(path))[1])
    (; t, r = col("r"), u = col("u"), rho = col("rho")[1:end-1], T = col("T")[1:end-1], Trad = col("Trad")[1:end-1])
end

mode = ARGS[1]
if mode == "bench"
    for (name, nsub) in (("SOD", 30), ("SUBC", 40), ("SUPC", 40))
        F = readdump(joinpath(@__DIR__, "..", "reference_runs", name, "$(name)_final.csv"))
        local snaps = deserialize(joinpath(@__DIR__, "..", "lib", "p2_runs", "$(name)_pure", "snaps.jls"))
        J = snaps[end]
        idx = unique(round.(Int, range(1, length(snaps), length = nsub)))
        series = [(t = s.time, r = s.r, rho = s.rho, T = s.T, Trad = Metrics.Trad(s), u = s.u) for s in snaps[idx]]
        writejs("bench_$(name).json", (name = name, fortran = F,
            julia = (t = J.time, r = J.r, u = J.u, rho = J.rho, T = J.T, Trad = Metrics.Trad(J)),
            nsteps_julia = length(snaps), series = series))
    end
elseif mode == "qstar"
    name = ARGS[2]                     # "A+B": run B continued A (QS_FROM); their snapshots are joined
    parts = split(name, "+"); name = String(parts[1])
    snaps = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(parts[1])_relax_snaps.jls"))
    for p in parts[2:end]
        append!(snaps, deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(p)_relax_snaps.jls"))[2:end])
    end
    snaps = filter(s -> all(isfinite, s.r) && all(isfinite, s.T), snaps)
    prof = length(ARGS) >= 3 ? ARGS[3] : joinpath(@__DIR__, "..", "reference_runs", "qs_profile_f0.1_md0.csv")
    lines = readlines(prof); hdr = Dict(Symbol(split(kv, "=")[1]) => parse(Float64, split(kv, "=")[2]) for kv in split(lstrip(lines[1], ['#', ' '])))
    C = [parse.(Float64, split(l, ",")) for l in lines[3:end]]
    names = split(lines[2], ","); ci(n) = findfirst(==(n), names)
    steady = (r = [c[ci("r")] for c in C], T = [c[ci("T")] for c in C], rho = [c[ci("rho")] for c in C],
              Lconv = [max(0.0, 1 - c[ci("nabla")] / c[ci("nabla_rad")]) for c in C])
    out = length(ARGS) >= 4 ? ARGS[4] : "qstar_$(name).json"
    # time history of every snapshot (shock, photosphere, accretion rate; qstar_summary.jl), profiles of at most 150
    hist = nothing
    if length(snaps) > 2 && length(snaps[end].u) > 0 && any(s -> minimum(s.u) < -1e7, snaps)      # infall runs only
        include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "qstar_summary.jl"))
        M0 = Base.invokelatest(QSum.mass, snaps[1])
        H = [Base.invokelatest(QSum.summary, s, hdr[:R], hdr[:L], M0) for s in snaps]
        hist = (t = [h.t for h in H], rsh = [h.rsh for h in H], Tpk = [h.Tpk for h in H], rph = [h.rph for h in H],
                Teff = [h.Teff for h in H], Lph = [h.Lph for h in H], Mdot_out = [h.Mdot for h in H], Mdot_sh = [h.Mdsh for h in H])
    end
    idx = unique(round.(Int, range(1, length(snaps), length = min(length(snaps), 150))))
    prof = (name = name, R = hdr[:R], L = hdr[:L], M = hdr[:M], steady = steady,
        snaps = [(t = s.t, r = s.r, u = s.u, rho = s.d, T = s.T, pg = s.pg, Er = s.er, Lrad = s.Lrad ./ hdr[:L], Lconv = s.Lconv ./ hdr[:L])
                 for s in snaps[idx]])
    writejs(out, hist === nothing ? prof : merge(prof, (hist = hist,)))
end
