"""
    Z0Physics

Zero-metallicity (X = 0.76, Y = 0.24) LTE equation of state and Rosseland/Planck opacities,
AD-ready (ForwardDiff Duals pass through every function). Copied from QuasiStarInfall.jl v2;
see docs/Z0_EOS_OPACITY.md for the validation record. The opacity table is stored as gzipped
raw Float64 (`z0/kaptable_Z0.bin.gz`), independent of the Julia version; `build_table` regenerates it.
"""
module Z0Physics
include(joinpath(@__DIR__, "z0", "constants.jl"))
include(joinpath(@__DIR__, "z0", "eos.jl"))
include(joinpath(@__DIR__, "z0", "opacity.jl"))

# Default table: log T 3.3–8.0, log R −12–1 (the steady model's). `Z0_TABLE=ext` selects the extended table
# (log T 3.0–8.0, log R −20–1, same nodes inside the default range) needed by cold, rarefied inflow.
const TABLEFILE = joinpath(@__DIR__, "z0", get(ENV, "Z0_TABLE", "") == "ext" ? "kaptable_Z0_ext.bin.gz" : "kaptable_Z0.bin.gz")
"Write a table: Int64 nT, nR, then lT, lR, lk (Rosseland), lkP (Planck), all Float64, gzipped."
function save_table(tab::KapTable, path = TABLEFILE)
    buf = IOBuffer()
    write(buf, Int64(length(tab.lT)), Int64(length(tab.lR)), tab.lT, tab.lR, vec(tab.lk), vec(tab.lkP))
    open(pipeline(`gzip -9 -c`, stdout = path), "w") do io; write(io, take!(buf)); end
    path
end
function load_table(path = TABLEFILE)
    io = IOBuffer(read(`gzip -dc $path`))
    nT, nR = read(io, Int64), read(io, Int64)
    rd(n) = collect(reinterpret(Float64, read(io, 8n)))
    lT = rd(nT); lR = rd(nR); lk = reshape(rd(nT * nR), nT, nR); lkP = reshape(rd(nT * nR), nT, nR)
    KapTable(lT, lR, lk, lkP)
end
const TAB = Ref{Any}(nothing)
table() = (TAB[] === nothing && (TAB[] = load_table()); TAB[]::KapTable)

"Gas-only pressure and specific internal energy (ionization included; radiation excluded — TITAN evolves E_r separately)."
gas_p(ρ, T) = Pgas(ρ, T)
gas_e(ρ, T) = _egas(saha(ρ, T), ρ, T)
end
