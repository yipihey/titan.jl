# Population III protostar formation and accretion from the ABN02 cloud: 1D radiation hydrodynamics (TITAN, pure Julia) with the
# coupled implicit primordial chemistry, H₂ line transfer, star light, the G2 network material and the GRIDEQ moving mesh — the one
# supported configuration is julia/popiii_stack.jl (docs/POPIII_STACK.md); launch with scripts/popiii_run/iceulrun.sh.
#   julia --project=<env> popiii_m2c.jl lines <IC shells>      (start from the IC; RESTART=<checkpoint> continues a run)
using Printf
include(joinpath(@__DIR__, "popiii_stack.jl"))
# the one supported configuration: defaults for every unset stack key
# TITAN_PKG (default off): load the module files as the precompiled package julia/TitanStack (no compile of the modules) instead of
# including them; "1": fall back to the includes if the image was built with a different load-time environment, "strict": stop then.
# The package needs the environment to have TitanStack developed (docs/POPIII_PRECOMPILE.md). Everything below is the same either way.
const PKGMODE = get(ENV, "TITAN_PKG", "0")
PKGMODE in ("0", "1", "strict") || error("TITAN_PKG=$PKGMODE: 0 (include the module files), 1 (the package, include path if its image does not match), strict (the package or stop)")
function pkg_load(mode)       # true: the modules are Main's constants, from the package image
    mode == "0" && return false
    t0 = time()
    try
        Core.eval(Main, :(using TitanStack))
        for nm in Core.eval(Main, :(TitanStack.MODULES)); Core.eval(Main, :(const $nm = getfield(TitanStack, $(QuoteNode(nm))))); end     # Main.<module>, as the includes define them (the H2Levels lookups of H2Lines/PopIIIMaterial go through Main)
        println("TITAN_PKG: TitanStack image loaded in ", round(time() - t0, digits = 1), " s (", length(Core.eval(Main, :(TitanStack.BAKED))), " guarded environment keys, ", Core.eval(Main, :(TitanStack.PC_APPLIED))[][1], " precompile statements)")
        return true
    catch e
        while !(e isa ErrorException) && (hasproperty(e, :error) || hasproperty(e, :ex)); e = hasproperty(e, :error) ? e.error : e.ex; end     # InitError / LoadError wrappers
        msg = e isa ErrorException ? e.msg : sprint(showerror, e)
        mode == "strict" && error("TITAN_PKG=strict: " * msg)
        println("TITAN_PKG=1: not using the TitanStack image, falling back to the include path. Reason: ", msg)
        return false
    end
end
const PKGON = pkg_load(PKGMODE)
if !PKGON
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl")); include(joinpath(@__DIR__, "PrimordialChem.jl"))
include(joinpath(@__DIR__, "SCVH.jl")); include(joinpath(@__DIR__, "PopIIIEqEOS.jl"))
include(joinpath(@__DIR__, "ContinuumOpacity.jl")); include(joinpath(@__DIR__, "G2Opacity.jl")); include(joinpath(@__DIR__, "ChemNet.jl")); include(joinpath(@__DIR__, "ChemRates.jl"))
include(joinpath(@__DIR__, "PopIIIMaterial.jl")); include(joinpath(@__DIR__, "H2Lines.jl"))
# LINEGPU=cuda|metal|cpu: formal solution of the thick H2 lines through KernelAbstractions (H2LinesKA.jl); needs an
# environment with KernelAbstractions and CUDA or Metal. Off unless set.
haskey(ENV, "LINEGPU") && include(joinpath(@__DIR__, "h2lines_gpu.jl"))
include(joinpath(@__DIR__, "StellarOpacity.jl")); include(joinpath(@__DIR__, "IonRay.jl"))
# MGROUPS (docs/POPIII_MULTIGROUP.md; default off): the pivoted banded solver for the wider block system
haskey(ENV, "MGROUPS") && ENV["MGROUPS"] != "0" && include(joinpath(@__DIR__, "JuliaSolver.jl"))
end
# S6-C1: HOT=1: stellar EOS, nuclear network and convective mixing in the hot core (T > HOTT1), see PopIIIMaterial.HOT
const HOTON = get(ENV, "HOT", "0") == "1"
if HOTON && !PKGON
    for f in ("StellarEOS", "NucNet", "StellarMaterial", "StarMaterial"); include(joinpath(@__DIR__, f * ".jl")); end
end
# INTEGRATOR=trbdf2 (default) | theta: second-order L-stable TR-BDF2 (Trbdf2.jl, G2) or TITAN's theta scheme
PKGON || (include(joinpath(@__DIR__, "Trbdf2.jl")); include(joinpath(@__DIR__, "NewtonKrylov.jl")); include(joinpath(@__DIR__, "StageCons.jl")))
# The driver: PopIIIRun.jl holds its definitions (part of the TitanStack image), popiii_run_body.jl the run itself, which PopIIIRun.main
# evaluates now (ENV and ARGS are read then, not when the module is loaded).
PKGON || (include(joinpath(@__DIR__, "H2Levels.jl")); include(joinpath(@__DIR__, "GridEq.jl")); include(joinpath(@__DIR__, "PopIIIRun.jl")))      # (with TITAN_PKG: Main.H2Levels, Main.GridEq and Main.PopIIIRun are TitanStack's, bound above)
PopIIIRun.main(ARGS)
