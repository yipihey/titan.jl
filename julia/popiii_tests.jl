# Population III project, G0: one entry point for every acceptance test so far. Each test runs as a separate process
# (the TITAN state is global), its output is parsed for the number that decides PASS/FAIL.
#   julia --project=env popiii_tests.jl [quick|full]       quick ≈ 5 min on many cores; full adds the M1 collapse tests
# Logs: lib/popiii/tests/<name>.log; summary printed and written to lib/popiii/tests/summary.txt.
using Printf
const J = @__DIR__
const LOG = joinpath(J, "..", "lib", "popiii", "tests"); mkpath(LOG)
tier = isempty(ARGS) ? "quick" : ARGS[1]

struct Test
    name::String; cmd::Cmd; tier::String
    check::Function            # (output::String) -> (pass::Bool, detail::String)
end
num(rx, s; k = 1) = (m = match(rx, s); m === nothing ? NaN : parse(Float64, m.captures[k]))
allnums(rx, s) = [parse(Float64, m.captures[1]) for m in eachmatch(rx, s)]
"The first match of rx in s as a string (\"?\" if none)."
grab(rx, s) = (m = match(rx, s); m === nothing ? "?" : strip(m.match))
env(c::Cmd, kv...) = addenv(c, kv...)
julia(args...; t = 1) = `julia -t $t --project=$(joinpath(J, "env")) $args`

tests = Test[
    # --- regression of the verified TITAN paths (PORTING_NOTES P6, M0.9) --------------------------------------
    Test("regress_SOD", env(julia(joinpath(J, "phase6_run.jl"), "pure", "SOD"), "Z0_TABLE" => "ext", "TITAN_NOLIB" => "1"), "quick",
         o -> (n = num(r"steps = (\d+)", o); (n == 464, "steps $n (reference 464 on this machine)"))),
    Test("regress_SUBC", env(julia(joinpath(J, "phase6_run.jl"), "pure", "SUBC"), "Z0_TABLE" => "ext", "TITAN_NOLIB" => "1"), "quick",
         o -> (n = num(r"steps = (\d+)", o); (n == 587, "steps $n (587)"))),
    Test("regress_SUPC", env(julia(joinpath(J, "phase6_run.jl"), "pure", "SUPC"), "Z0_TABLE" => "ext", "TITAN_NOLIB" => "1"), "quick",
         o -> (n = num(r"steps = (\d+)", o); (n == 713, "steps $n (713)"))),
    # --- M1 ---------------------------------------------------------------------------------------------------
    Test("m1_static", env(julia(joinpath(J, "popiii_m1.jl"), "static"), "TITAN_NOLIB" => "1", "TAG" => "_test"), "quick",
         o -> (occursin("PASS", o), grab(r"max \|Δr/r\|[^\n]*", o))),
    Test("m1_dust", env(julia(joinpath(J, "popiii_m1.jl"), "dust", "40"), "TITAN_NOLIB" => "1", "TAG" => "_test",
                        "T0" => "1e-6", "THETA" => "0.5", "CDYN" => "0.003", "STOL" => "0.1"), "full",
         o -> (t = num(r"t/t_ff = ([\d.]+)", o); e = num(r"max \|r/r_exact − 1\| = ([\d.e+-]+)", o);
               (abs(t / 0.999574 - 1) < 1e-3, @sprintf("arrival t/t_ff = %.6f (exact 0.99958 at R/R0 = 0.01); max |Δr/r| to 0.01: %.2e", t, e)))),
    Test("m1_be_disc", env(julia(joinpath(J, "popiii_m1.jl"), "be", "60", "discrete"), "TITAN_NOLIB" => "1", "TAG" => "_test"), "quick",
         o -> (occursin("PASS", o), grab(r"max \|Δρ/ρ\| = [^\n]*", o))),
    Test("m1_be_le120", env(julia(joinpath(J, "popiii_m1.jl"), "be", "120", "le"), "TITAN_NOLIB" => "1", "TAG" => "_test"), "full",
         o -> (occursin("PASS", o), grab(r"max \|Δρ/ρ\| = [^\n]*", o))),
    Test("m1_lp", env(julia(joinpath(J, "popiii_m1.jl"), "lp", "200"), "TITAN_NOLIB" => "1", "TAG" => "_test"), "full",
         o -> (d = allnums(r"\n +[\d.e+-]+ +[\d.e+-]+ +([\d.]+) +-?[\d.]+\n", o); m = isempty(d) ? NaN : maximum(d);
               (7.5 < m < 7.9, @sprintf("max 4πGρr²/c_s² after 8 decades %.3f (M1b: 7.68; LP 8.84)", m)))),
    # --- M2/M3 physics modules ----------------------------------------------------------------------------------
    Test("h2lines", julia(joinpath(J, "h2lines_test.jl"); t = 8), "quick",
         o -> (r2 = allnums(r"ratio ([\d.]+) +\(Λ", o); r3 = allnums(r"total L ratio ([\d.]+)", o);
               r4 = allnums(r"F/F_D = ([\d.]+)", o) ; fd = filter(x -> x > 0.98, allnums(r"(0\.9[89]\d+|1\.00\d+)", o));
               ok = all(x -> abs(x - 1) < 1e-4, r2) && all(x -> abs(x - 1) < 0.01, r3);
               (ok, @sprintf("thin %s; sphere total L ratio ∈ [%.4f, %.4f]", isempty(r2) ? "?" : string(round(minimum(r2), digits = 6)),
                             isempty(r3) ? NaN : minimum(r3), isempty(r3) ? NaN : maximum(r3))))),
    Test("h2levels", julia(joinpath(J, "h2levels_test.jl")), "full",
         o -> (occursin("ALL PASS", o), grab(r"(ALL PASS|FAILED)[^\n]*", o))),
    Test("stellareos", env(julia(joinpath(J, "stellareos_test.jl")), "TITAN_NOLIB" => "1"), "full",
         o -> (occursin("ALL PASS", o), grab(r"(ALL PASS|FAILED)[^\n]*", o))),
    Test("chem_saha_db", julia("-e", """include("$(joinpath(J, "PrimordialChem.jl"))"); const PC = PrimordialChem
        ρ = 1e16 * PC.μH; T = 2200.0; y = (1e-11, 1e-18, 0.46, 0.0); u = PC.eos_rhfc02(ρ, T, y)[1]
        a = PC.step_coupled(y, ρ, u, T, 1e6, PC.rates_rhfc02(; k5 = T -> PC.k5_db(T))); b = PC.step_equilibrium(ρ, u, T, 1e6)
        println("Tdb=", a[2], " Tsaha=", b[2], " fdb=", 2a[1][3], " fsaha=", 2b[1][3])"""), "quick",
         o -> (a = num(r"Tdb=([\d.e+-]+)", o); b = num(r"Tsaha=([\d.e+-]+)", o); (abs(a / b - 1) < 1e-6, @sprintf("detailed-balance network T = %.4f K vs Saha %.4f K", a, b)))),
    Test("onezone_ra04", env(julia(joinpath(J, "popiii_m2_onezone.jl"), "ra04", "pss"), "NMAX" => "1e14"), "quick",
         o -> (f = num(r"\n +12\.\d+ +12\.\d+ +[\d.]+ +([\d.e+-]+)", o); (0.95 < f < 1.0, @sprintf("f_H2 at log ρ/m_H ≈ 12.4: %.4f (M2a: 0.971)", f)))),
]

sel = filter(t -> tier == "full" || t.tier == "quick", tests)
println("running $(length(sel)) tests ($tier)")
results = Vector{Any}(undef, length(sel))
@sync for (i, t) in enumerate(sel)
    @async begin
        out = try read(pipeline(t.cmd; stderr = stdout), String) catch e; "ERROR: " * sprint(showerror, e) end
        write(joinpath(LOG, t.name * ".log"), out)
        results[i] = try t.check(out) catch e; (false, "check failed: " * sprint(showerror, e)[1:min(end, 120)]) end
    end
end
io = open(joinpath(LOG, "summary.txt"), "w")
for (t, (ok, det)) in zip(sel, results)
    line = @sprintf("%-14s %s  %s", t.name, ok ? "PASS" : "FAIL", det)
    println(line); println(io, line)
end
close(io)
exit(all(first, results) ? 0 : 1)
