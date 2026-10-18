# Data files for the report's interactive figures (report/app/data/*.json).
#   julia --project=env export_report_data.jl bench      benchmark profiles: Fortran 1995-code final dumps vs pure-Julia runs
#   julia --project=env export_report_data.jl qstar NAME[+CONT] [PROFILE.csv] [OUT]   a quasi-star run (lib/qs_runs/NAME_relax_snaps.jls,
#                                                     with the snapshots of its continuation CONT appended)
#   Z0_TABLE=ext julia --project=env export_report_data.jl spectra NAME...   emergent spectra lib/spectra/NAME.jls (spectra_run.jl)
#                                                     → spectra.json, with band luminosities
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
elseif mode == "spectra"
    include(joinpath(@__DIR__, "Z0Physics.jl"))
    names = ARGS[2:end]
    LABEL = Dict("star_t0" => "steady star before infall (t = 0)",
                 "A220_late_ext" => "infall, 4.2e8 s, LTE infall + free-fall extension",
                 "A220_late_ion" => "infall, 4.2e8 s, infall fully ionized (limit)",
                 "A220_late_noext" => "infall, 4.2e8 s, grid only (to 1.20 R)",
                 "A220_late_ext_nlte" => "infall, 4.2e8 s, non-LTE hydrogen in the infall",
                 "A220_late_bare" => "infall, 4.2e8 s, star only (infall above the shock removed)",
                 "A220_late_hotion" => "infall, 4.2e8 s, cooling layer fully ionized (test)",
                 "R2_late_ext" => "Eulerian run R2, 5.8e7 s, LTE infall + extension",
                 "A220_corona_ext" => "early ramp, 7.9e6 s (Ṁ = 0.003 M⊙/yr at the shock, 0.14 at the grid edge), LTE + extension")
    hP, kB, cc = 6.62607015e-27, 1.380649e-16, 2.99792458e10; keV = 1.602176634e-9
    band(ν, L, ν1, ν2) = (s = 0.0; for i in 1:length(ν)-1
        a, b = max(ν[i], ν1), min(ν[i+1], ν2); b > a || continue
        f(x) = (L[i] > 0 && L[i+1] > 0) ? L[i] * (L[i+1] / L[i])^((log(x) - log(ν[i])) / (log(ν[i+1]) - log(ν[i]))) : L[i] + (L[i+1] - L[i]) * (x - ν[i]) / (ν[i+1] - ν[i]); s += 0.5 * (f(a) + f(b)) * (b - a); end; s)   # power-law interpolation (edges)
    Lat(ν, L, λÅ) = (x = cc / (λÅ * 1e-8); i = argmin(abs.(ν .- x)); L[i])
    runs = []
    for n in names
        d = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "$n.jls"))
        ν, L = d.ν, d.L
        Ltot = band(ν, L, 0.0, Inf)
        bands = (X_05_2 = band(ν, L, 0.5keV / hP, 2keV / hP), X_2_10 = band(ν, L, 2keV / hP, 10keV / hP),
                 X_10_40 = band(ν, L, 10keV / hP, 40keV / hP), EUV = band(ν, L, 13.6 * 1.602176634e-12 / hP, 0.1keV / hP),
                 FUV = band(ν, L, cc / 2000e-8, cc / 912e-8), opt = band(ν, L, cc / 7000e-8, cc / 3646e-8), total = Ltot)
        Dbalmer = log10(Lat(ν, L, 3664.0) / Lat(ν, L, 3628.0))       # Balmer jump: L_ν just red / just blue of 3646 Å
        R0 = 8.4352e15
        prof = (r = d.atm.r ./ R0, T = d.atm.T, T_titan = d.T_titan, rho = d.atm.ρ,
                xH = haskey(d, :xH) ? d.xH : [Z0Physics.saha(d.atm.ρ[k], d.atm.T[k]).xH for k in eachindex(d.atm.r)])
        tab = Z0Physics.table(); τ = 0.0; kph = length(d.atm.r)          # TITAN-style photosphere: Rosseland τ = 2/3
        for k in length(d.atm.r):-1:2
            τ += Z0Physics.kap(tab, d.atm.ρ[k], d.T_titan[k]) * d.atm.ρ[k] * (d.atm.r[k] - d.atm.r[k-1])
            τ >= 2 / 3 && (kph = k; break)
        end
        rph = d.atm.r[kph]; Teff = (d.Ltitan / (4π * 5.670374419e-5 * rph^2))^0.25
        lf = joinpath(@__DIR__, "..", "lib", "spectra", "lines_$n.jls")
        lines = isfile(lf) ? deserialize(lf) : []
        push!(runs, (name = n, lines = lines, Lline = haskey(d, :Lline) ? d.Lline : 0.0, label = get(LABEL, n, n), run = d.run, t = d.t, Mdot = d.Mdot * Z0Physics.yr / Z0Physics.Msun,
            Ltitan = d.Ltitan, rph = rph / R0, Teff = Teff, lam = cc ./ ν .* 1e8, E = hP .* ν ./ keV, nuLnu = ν .* L, bands = bands, Dbalmer = Dbalmer,
            nlucy = d.nlucy, ion = haskey(d, :ion) ? d.ion : false, prof = prof))
        @printf("%-18s ∫L = %.4e (TITAN %.4e, ratio %.4f)  X(0.5–2) %.2e  X(2–10) %.2e  X(10–40) %.2e  EUV %.2e  FUV %.3e  opt %.3e  D_B = %.3f  r_ph/R = %.5f  T_eff = %.0f K\n",
                n, Ltot, d.Ltitan, Ltot / d.Ltitan, bands.X_05_2, bands.X_2_10, bands.X_10_40, bands.EUV, bands.FUV, bands.opt, Dbalmer, rph / R0, Teff)
    end
    writejs("spectra.json", (R = 8.4352e15, runs = runs))
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
