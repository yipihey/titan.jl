# usage: julia verify_phase1.jl PROBLEM_NAME NN [ntrace]
# Runs one benchmark through libtitan.so via TitanFortran (one model per process),
# then checks BITWISE equality of the final COMMON state with the reference dump
# and, if ntrace > 0, per-step equality of COMMON /new/ with trace.bin written by
# the standalone titan_trace.x (same sources, same flags).
include(joinpath(@__DIR__, "TitanFortran.jl"))
using .TitanFortran
const TF = TitanFortran

name, nn = ARGS[1], parse(Int, ARGS[2])
ntrace = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0
root = TF.ROOT
rundir = joinpath(root, "lib", "jl_runs", name)
ref = TF.read_dump(joinpath(root, "reference_runs", name, "$name.dmp"))

TF.init(nn; rundir)
recs = Tuple{Int32,Int32,Float64,Float64,Vector{UInt8}}[]
nsteps = 0
t0 = time()
while true
    js = TF.getscalar(:jstep)
    s = TF.step!(); global nsteps += 1
    if nsteps <= ntrace
        push!(recs, (js, TF.getscalar(:iter), TF.getscalar(:timen), TF.getscalar(:dtime),
                     copy(TF.blockbytes(:new))))
    end
    s != 0 && break
end
wall = time() - t0
cmpres = TF.compare_with_dump(ref)
st = TF.get_state()
TF.finalize()

nbad = count(!r[2] for r in cmpres)
println("$name: steps driven = $nsteps, final jstep = $(st.jstep) (reference $(ref[:jstep])), ",
        "t = $(st.time), wall = $(round(wall, digits=2)) s")
println("$name: final-state bitwise comparison over $(length(cmpres)) dumped quantities ",
        "(full mgr-length arrays): ", nbad == 0 ? "ALL IDENTICAL" : "$nbad DIFFER")
for r in cmpres
    r[2] && continue
    a = r[1] in TF.DUMP_ARRS || r[1] in TF.DUMP_ARR0 ? findall(reinterpret(UInt64, ref[r[1]]) .!=
            reinterpret(UInt64, TF.var(r[1]) |> collect)) : Int[]
    println("   differs: ", r, " at Fortran index ", a, " (ngrs=$(st.ngrs), ngre=$(st.ngre))")
end

if ntrace > 0 let
    trace = joinpath(rundir, "..", "$(name)_trace", "trace.bin")
    nb = TF.block_size(:new); reclen = 4 + 4 + 8 + 8 + nb
    data = read(trace); nrec = length(data) ÷ reclen
    @assert length(data) == nrec * reclen
    nok = 0; firstbad = 0
    for i in 1:min(nrec, length(recs))
        o = (i - 1) * reclen
        js = reinterpret(Int32, data[o+1:o+4])[1]; it = reinterpret(Int32, data[o+5:o+8])[1]
        tn = reinterpret(Float64, data[o+9:o+16])[1]; dt = reinterpret(Float64, data[o+17:o+24])[1]
        blk = data[o+25:o+24+nb]
        r = recs[i]
        ok = js == r[1] && it == r[2] && TF.bitsequal(tn, r[3]) && TF.bitsequal(dt, r[4]) && blk == r[5]
        ok ? (nok += 1) : (firstbad == 0 && (firstbad = i))
    end
    println("$name: per-step check of jstep, iter, timen, dtime and all of COMMON /new/ ($nb bytes) ",
            "for $(min(nrec, length(recs))) steps: $nok identical",
            firstbad == 0 ? "" : ", first mismatch at step $firstbad")
end end
