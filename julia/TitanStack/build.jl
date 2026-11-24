# Build (or rebuild) the TitanStack image with the environment of THIS process. The image freezes the load-time environment (see
# src/envguard.jl); build it from the configuration it will be used for, e.g.
#   env $(cat run.env) MGROUPS=4 julia --project=<env with TitanStack developed> julia/TitanStack/build.jl
# (popiii_stack.jl's defaults are applied first, as popiii_m2c.jl does). Forces a recompile even if an image exists.
include(joinpath(@__DIR__, "..", "popiii_stack.jl"))
using Pkg
id = Base.identify_package("TitanStack")
id === nothing && error("TitanStack is not in the active project: Pkg.develop(path = \"julia/TitanStack\") in your environment first")
t = @elapsed Base.compilecache(id)
println("TitanStack image built in ", round(t, digits = 1), " s")
