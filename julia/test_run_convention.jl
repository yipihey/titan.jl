# Run-convention tests for qstar_run.jl:  julia --project=env test_run_convention.jl
# (a) schema lists every QS_* env; (b) env sources are recorded; (c) a tiny relax run and a resume pass run-check.
# (c) needs the Fortran library only if you run without TITAN_NOLIB; qstar_run.jl is the pure-Julia backend (TITAN_NOLIB=1).
using Test, RunConvention
const HERE = @__DIR__
const ENVP = joinpath(HERE, "env")
const RUNCHECK = get(ENV, "RUN_CHECK", "/Users/tabel/Projects/koine/target/debug/run-check")
julia_cmd(args...; env = Dict{String,String}()) =
    setenv(`$(Base.julia_cmd()) --project=$ENVP $(joinpath(HERE, "qstar_run.jl")) $args`,
           merge(Dict(k => v for (k, v) in ENV if !startswith(k, "QS_")), Dict("TITAN_NOLIB" => "1"), env))

@testset "run convention: qstar_run.jl" begin
    @testset "schema" begin
        doc = RunConvention.parse_json(read(julia_cmd("--schema"), String))
        envs = Set(String(get(p["x-run"], "env", "")) for (_, p) in doc["properties"])
        src = read(joinpath(HERE, "qstar_run.jl"), String)
        used = Set(m.match for m in eachmatch(r"QS_[A-Z0-9_]*[A-Z0-9]", src))
        delete!(used, "QS_SET")                      # the open group `common`, not a single parameter
        @test !isempty(used) && issubset(used, envs)
        @test haskey(doc["x-run"]["groups"], "common") && doc["x-run"]["groups"]["common"]["open"] == true
        @test haskey(doc["properties"], "common.ql1")
        @test doc["properties"]["sink.tau"]["x-run"]["unit"] == "s"
        @test doc["properties"]["setup.mdot"]["x-run"]["unit"] == "Msun/yr"
        @test doc["properties"]["io.every"]["x-run"]["mutable"] == "always"
    end
    @testset "env sources" begin
        out = read(julia_cmd("relax", "--show-params"; env = Dict("QS_SINK_TAU" => "3e6", "QS_GASROW" => "1", "QS_SET" => "ql1=1e-4")), String)
        rows = Dict(split(l, '\t')[1] => split(l, '\t')[2:3] for l in split(out, '\n'; keepempty = false))
        @test parse(Float64, rows["sink.tau"][1]) == 3e6 && rows["sink.tau"][2] == "env:QS_SINK_TAU"
        @test rows["common.ql1"][2] == "env:QS_SET" || rows["common.ql1"][2] == "cli"
        @test rows["io.every"][2] == "default"
        # sink.tau without physics.gasrow is refused at parse time (x-run.requires), not at the sink.
        @test_throws ProcessFailedException read(julia_cmd("relax", "--show-params"; env = Dict("QS_SINK_TAU" => "3e6")), String)
    end
    @testset "relax run + resume conform" begin
        if !isfile(RUNCHECK)
            @info "SKIPPED: run-check binary not found at $RUNCHECK (set RUN_CHECK)"
        elseif !isfile(joinpath(HERE, "..", "reference_runs", "init", "QSTAR.tstate.gz"))
            @info "SKIPPED: reference_runs/init/QSTAR.tstate.gz missing"
        else
            name = "R3TEST$(getpid())"
            qs = joinpath(HERE, "..", "lib", "qs_runs")
            e1 = Dict("QS_NAME" => name, "QS_INIT" => "QSTAR", "QS_EVERY" => "1")
            @test success(pipeline(julia_cmd("relax", "1e30", "5"; env = e1); stdout = devnull))
            d1 = joinpath(qs, name)
            @test success(`$RUNCHECK $d1`)
            ck = sort(filter(f -> startswith(f, "checkpoint-"), readdir(d1)))
            @test !isempty(ck)
            @test isfile(joinpath(qs, "$(name)_relaxed.tstate.gz"))          # the old flat file is still written
            e2 = Dict("QS_NAME" => name, "QS_EVERY" => "1", "QS_FROM" => joinpath(d1, last(ck)))
            @test success(pipeline(julia_cmd("relax", "1e30", "2"; env = e2); stdout = devnull))
            d2 = joinpath(qs, "$name.2")
            @test success(`$RUNCHECK $d2`)
            rj = RunConvention.parse_json(read(joinpath(d2, "run.json"), String))
            p1 = RunConvention.parse_json(read(joinpath(d1, "run.json"), String))
            @test rj["lineage"]["parent"]["run"] == p1["id"]
            @test rj["lineage"]["parent"]["relation"] == "resume"
            for d in (d1, d2); rm(d; recursive = true); end
            for f in ("$(name)_relaxed.tstate.gz", "$(name)_relax_snaps.jls"); rm(joinpath(qs, f); force = true); end
        end
    end
end
