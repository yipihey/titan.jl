# S7c-1: Perf.jl wrapper of the protostar driver (s6b_run.jl): the whole run is section :run, per-step samples every PERFEVERY [100] steps from the driver's `series`, counters from the driver.
# Records lib/popiii/perf/<PERF>_{meta.txt,sections.csv,counters.csv,samples.csv}. Usage: PERF=s7c1_x julia -t 3 --project=env s7c1_perf.jl  (environment of s6b_run.jl)
include(joinpath(@__DIR__, "Perf.jl"))
const PERFNAME = get(ENV, "PERF", "s7c1_run")
Perf.reset!(PERFNAME; stage = "protostar", driver = "s6b_run.jl (s7c1_perf.jl)", out = get(ENV, "OUT", "run"), mdot = get(ENV, "MDOT", get(ENV, "MDOTFILE", "")), prec = get(ENV, "PRECTHICK", "0"),
            extrap = get(ENV, "EXTRAP", "default"), estau = get(ENV, "ESTAU", "default"), lstarq = get(ENV, "LSTARQ", "default"))
Perf.@section :run include(joinpath(@__DIR__, "s6b_run.jl"))
let n = length(Main.series), every = parse(Int, get(ENV, "PERFEVERY", "100"))
    for (i, row) in enumerate(Main.series)
        i % every == 0 && Perf.sample!(step = i, t_yr = row[1], dt_yr = row[2], M_Msun = row[3], R_Rsun = row[4], Lstar_Lsun = row[5], Tc = row[8])
    end
    Perf.count!(:steps, Main.nstep[]); Perf.count!(:step_fail, Main.FAILS[])
    if isdefined(Main, :S7bStar)
        for (k, (c, s)) in Main.S7bStar.PERFD; Perf.add!(Symbol("s7b_", k), round(UInt64, s * 1e9), 0); Perf.count!(Symbol("s7b_", k, "_calls"), c); end
        Perf.count!(:stall_recoveries, Main.S7bStar.NRECOVER[])
    end
end
Perf.write(); Perf.summary()
