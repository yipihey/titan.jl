# LINEGPU=cuda|metal|cpu : the formal solution of the optically thick H₂ lines of an ALI call (H2Lines.cooling_nlte) as one batch
# on a KernelAbstractions backend (H2LinesKA.jl) instead of the threaded CPU loop over `transfer_line_fast!`.
#
#   include after H2Lines.jl (popiii_m2c.jl does, when LINEGPU is set), or call `h2lines_gpu!("cuda")` by hand.
#
# Needs KernelAbstractions and the backend package (CUDA.jl on NVIDIA, Metal.jl on Apple silicon) in the active project; the
# production environment has neither, so nothing here is loaded unless LINEGPU is set.
# LINEGPUFLOAT=32|64 : float type on the device (default 32; Metal has no Float64, and consumer NVIDIA cards are ~30x slower in
# Float64). In Float32 the net line cooling differs from the CPU path by ~1e-8 of its peak value (H2LinesKA.jl header);
# "cpu" (KernelAbstractions' CPU backend) in Float64 reproduces the CPU path to round-off but is slow: for validation only.
# Acts only on the LINEFAST=1 LINESRC=zone path; everything else (opacities, rate equations, thin lines, the escape-probability
# call of every step) is unchanged and stays on the host.

"""
    h2lines_gpu!(name; float = nothing)

Switch the formal solution of H2Lines to the KernelAbstractions backend `name` ("cuda", "metal", "cpu"; "off"/"" switches back to
the CPU path). The backend package must be loadable. Returns the backend object (or `nothing`).
"""
function h2lines_gpu!(name::AbstractString; float = nothing)
    name = lowercase(name)
    if name in ("", "0", "off", "none")
        H2Lines.KA_BACKEND[] = nothing
        return nothing
    end
    name in ("cuda", "metal", "cpu") || error("LINEGPU=$name: expected cuda, metal or cpu")
    isdefined(Main, :H2LinesKA) || Base.include(Main, joinpath(@__DIR__, "H2LinesKA.jl"))
    name == "cuda" && (isdefined(Main, :CUDA) || Core.eval(Main, :(using CUDA)))
    name == "metal" && (isdefined(Main, :Metal) || Core.eval(Main, :(using Metal)))
    Base.invokelatest(_h2lines_gpu_set!, name, float)
end

function _h2lines_gpu_set!(name, float)
    KA = Main.H2LinesKA.KA
    T = float === nothing ? (get(ENV, "LINEGPUFLOAT", "32") == "64" ? Float64 : Float32) : float
    name == "metal" && T !== Float32 && error("LINEGPU=metal: Metal has no Float64; use LINEGPUFLOAT=32")
    be = name == "cuda" ? Main.CUDA.CUDABackend() : name == "metal" ? Main.Metal.MetalBackend() : KA.CPU()
    name == "cuda" && !Main.CUDA.functional() && error("LINEGPU=cuda: CUDA.jl is not functional on this machine")
    H2Lines.KA_MOD[] = Main.H2LinesKA; H2Lines.KA_FLOAT[] = T; H2Lines.KA_BACKEND[] = be
    (H2Lines.FAST[] && H2Lines.SRCZONE[]) || @warn "LINEGPU acts only with LINEFAST=1 LINESRC=zone; the CPU path will be used"
    println("LINEGPU: formal solution of the thick H2 lines on ", typeof(be), " in ", T)
    be
end

haskey(ENV, "LINEGPU") && h2lines_gpu!(ENV["LINEGPU"])
