# G2 material audit: load the material modules of a test tree (git archive of the branch), standalone (no TITAN).
# usage: TREE=<dir with julia/> include("load.jl")
ENV["HECHEM"] = get(ENV, "HECHEM", "1")          # as zdiagB (run.env HECHEM=1): table with He Saha, network with He chemistry
const JD = joinpath(ENV["TREE"], "julia")
using ForwardDiff, Printf
include(joinpath(JD, "PrimordialChem.jl")); include(joinpath(JD, "SCVH.jl")); include(joinpath(JD, "PopIIIEqEOS.jl"))
include(joinpath(JD, "ChemNet.jl")); include(joinpath(JD, "ContinuumOpacity.jl")); include(joinpath(JD, "Z0Physics.jl"))
include(joinpath(JD, "StellarOpacity.jl")); include(joinpath(JD, "StellarOpacityLow.jl"))
const PC = PrimordialChem; const EQ = PopIIIEqEOS; const CN = ChemNet; const COp = ContinuumOpacity; const SO = StellarOpacity
const Z0 = Z0Physics; const SL = StellarOpacityLow
const XH = 1 / (1 + 4PC.yHe)
const EVERG = 1.602176634e-12; const E21 = 10.19881 * EVERG
