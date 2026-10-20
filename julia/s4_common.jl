# Shared set-up for the S4 tests (included by the s4_*.jl scripts). Run with TITAN_NOLIB=1.
using Printf, ForwardDiff
for f in ("TitanDriver", "JuliaPhysics", "TitanResidFast", "BlockSolver", "TitanCore", "Z0Physics", "TitanMaterials",
          "ProblemSetup", "StellarEOS", "StellarOpacity", "NucNet", "StellarMaterial", "StarMaterial")
    include(joinpath(@__DIR__, f * ".jl"))
end
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TC = TitanCore
const TM = TitanMaterials; const Z0 = Z0Physics; const SE = StellarEOS; const NN = NucNet; const SM = StarMaterial
const Msun = 1.98847e33; const Lsun = 3.828e33; const Rsun = 6.957e10
