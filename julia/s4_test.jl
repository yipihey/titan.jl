# S4 test driver: runs the S4 tests (each in its own process, the TITAN state is global) and prints PASS/FAIL.
#   TITAN_NOLIB=1 julia --project=env s4_test.jl [a|b|c|d|all]     (a, b, c ≈ 5 min total; d is the long star run, see s4_d.jl)
# a  no-burning regression: StarGas (Z0 opacity, no network) vs StellarGas, bit for bit, 25 steps of a rad-hydro star
# b  closed box: TITAN step sequence vs NucNet.burn_coupled (T(t), Y_D(t)), energy conservation 1e-10
# c  mixing operator: conservation 1e-14, uniform convective region, radiative zones untouched
# d  100 Msun Eddington-model star with burning and convective core (long; reads lib/popiii/s4_d_*.log)
using Printf
const J = @__DIR__
which = isempty(ARGS) ? "all" : ARGS[1]
julia(args...; env = String[]) = addenv(`julia -t 4 --project=$(joinpath(J, "env")) $args`, "TITAN_NOLIB" => "1", env...)
function run1(name, script; env = String[])
    out = try read(pipeline(julia(joinpath(J, script); env); stderr = stdout), String) catch e; "ERROR: " * sprint(showerror, e) end
    write(joinpath(J, "..", "lib", "popiii", "tests", "s4_$name.log"), out)
    println("---- $name ($script)"); print(out)
    out
end
mkpath(joinpath(J, "..", "lib", "popiii", "tests"))
res = Pair{String,Bool}[]
(which in ("a", "all")) && (o = run1("a", "s4_a.jl"); push!(res, "a: no-burning regression, StarGas vs StellarGas bit-for-bit" => occursin(r"\nPASS", "\n" * o)))
if which in ("b", "all")
    for ef in ("0.05", "0.01")
        o = run1("b_efrac$ef", "s4_b.jl"; env = ["EFRAC" => ef])
        push!(res, "b: closed-box D burning, efrac $ef (energy 1e-10, sequence vs burn_coupled, time-step error vs fine reference)" =>
              occursin("PASS energy", o) && occursin("PASS sequence", o) && occursin("PASS time-step error", o))
    end
end
(which in ("c", "all")) && (o = run1("c", "s4_c.jl"); push!(res, "c: mixing operator" => occursin(r"\nPASS", "\n" * o)))
println("\n==== S4 summary")
for (n, ok) in res; println(ok ? "PASS  " : "FAIL  ", n); end
println(all(last, res) ? "ALL PASS" : "FAILED")
