"""
    TitanDriver

Phase-2 scaffold: TITAN's time-step control flow (step.f + the titan.f loop body)
re-implemented in Julia.  Every physics/numerics routine is still the Fortran
one, called through ccall, with all state in the mapped COMMON blocks.  Routines
can be swapped at two levels:

* top level (`ROUTINES[:matgen] = f` …): anything called directly by step;
* deep level (`install_hooks!`): `eos`/`opac`, which Fortran routines call
  internally, via the C shim in lib/src/hooks.c (eos_ -> Julia or eos_f_).
"""
module TitanDriver

include(joinpath(@__DIR__, "TitanFortran.jl"))
using .TitanFortran
const TF = TitanFortran
export TF

# E2E-A: performance recording (Perf.jl, no dependencies). `TitanDriver.Perf` is the single instance the whole stack and the drivers share.
include(joinpath(@__DIR__, "Perf.jl"))

# ------------------------------------------------------------------ ccalls
const PTR = Dict{Symbol,Ptr{Cvoid}}()
fptr(s::Symbol) = get!(() -> TF.Libdl.dlsym(TF.load(), s), PTR, s)
fcall(s::Symbol) = ccall(fptr(s), Cvoid, ())
fmsg(code, xs = 0.0) = ccall(fptr(:port_msg_), Cvoid, (Ref{Int32}, Ref{Float64}), Int32(code), xs)
"Where Julia-backend messages go (TITAN writes them to the terminal and the .out file)."
const MSGIO = Ref{IO}(stdout)
"step.f messages; Fortran port_msg on the COMMON backend, Julia text (and errors instead of STOP) on the Julia backend."
function msg(code, xs = 0.0)
    TF.BACKEND[] === :common && return fmsg(code, xs)
    io = MSGIO[]
    code == 1 && println(io, "  divergence. jtry = ", g(:jtry))
    code == 2 && error("TITAN step01: failure to converge after $(g(:ntry)) attempts (jstep $(g(:jstep)), t = $(g(:timen)))")
    code == 3 && println(io, "  changes over timestep ", xs, " times larger than limit")
    code == 4 && error("TITAN step02: total changes still too large after $(g(:nback)) attempts")
    code == 5 && println(io, "  repeat integration step with dtime = ", g(:dtime))
    nothing
end
"EOS / opacity entry points used by the Julia routines (Fortran through the hook shim by default)."
const EOS = Ref{Function}(() -> fcall(:eos_))
const OPAC = Ref{Function}(() -> fcall(:opac_))

"Routines called by step; replace an entry to swap in a Julia version."
const ROUTINES = Dict{Symbol,Function}(
    s => (() -> fcall(Symbol(s, "_"))) for s in
    (:peq, :matgen, :matslv, :maxdel, :update, :timstp, :reset, :totnrg))
"Run routine `s` of the step. `matgen` (residual + Jacobian assembly) and `matslv` (block solve) are timed by Perf (`:matgen`, `:linsolve`)."
r(s) = s === :matgen ? Perf.@section(:matgen, ROUTINES[s]()) : s === :matslv ? Perf.@section(:linsolve, ROUTINES[s]()) : ROUTINES[s]()

# scalar accessors (zero-copy 0-dim views, cached)
const SV = Dict{Symbol,Array}()
sv(s::Symbol) = get!(() -> TF.var(s), SV, s)
g(s) = sv(s)[]
s!(s, x) = (a = sv(s); a[] = convert(eltype(a), x))

"""
Bookkeeping fixes (Phase 4), ON by default in the port:

* `FIX_CLOCK[]`: after a Newton-divergence retry (and in the rarely used
  "changes too large" branch) step.f calls reset, which sets timen = timeo, so the
  retried step is integrated with dtime but the clock does not advance.  Fixed by
  setting timen = timeo + dtime after reset.
* `FIX_ENERGY[]`: step.f calls timstp BEFORE totnrg, so the energy bookkeeping
  (boundary work, luminosity, surface and dissipation terms) is multiplied by the
  dt planned for the NEXT step (20x too large before a divergence retry).  Fixed
  by calling totnrg before timstp, i.e. with the dt actually taken.

`bitwise!(true)` (or env TITAN_BITWISE=1) switches both off for a bitwise
reproduction of TITAN 1.2.  Neither fix changes the solution of any step; they
change the time labels (hence tmax termination) and the energy diagnostic.
"""
const FIX_CLOCK = Ref(true)
const FIX_ENERGY = Ref(true)
"""
`EDD_MATERIAL_SOURCE[]`: source function of the Eddington-factor formal solution (eddfac).
false (TITAN): thermal absorption = χ_F ρ − σ_T n_e, i.e. the flux mean minus electron scattering. Exact for
TITAN's grey built-in opacities. For a tabulated material (χ_F a Rosseland mean, κ_P a Planck mean) that
difference is round-off/table noise wherever scattering dominates, and in hot thin gas it multiplies a Planck
function 1e9 times J, so the formal solution returns garbage (f clamped to 0 or 1).
true (set by `TitanMaterials.use!` for non-built-in materials): S = (κ_P ρ B + σ_T n_e J)/(κ_P ρ + σ_T n_e).
"""
const EDD_MATERIAL_SOURCE = Ref(false)
bitwise!(on::Bool = true) = (FIX_CLOCK[] = !on; FIX_ENERGY[] = !on; nothing)
get(ENV, "TITAN_BITWISE", "0") == "1" && bitwise!(true)

"""
`STAGE_HOOK[]` (S5): `f(kind::Symbol, dt)` called by the step functions right after every converged implicit solve, while
the new arrays hold the stage state and the old arrays its starting state: `:single` (TITAN's θ-scheme step, dt = dtime),
and from `Trbdf2` `:tr` (stage 1), `:bdf` (stage 2) and `:be` (backward-Euler fallback). Used by StarMaterial to book the
energy injected by the nuclear source and the energy lost through the boundary with the weights of the scheme.
"""
const STAGE_HOOK = Ref{Any}(nothing)
"The full step size Δ of the step being solved (TR-BDF2 solves stages of γΔ, (1−γ)/(2−γ)Δ; `dtime` is the stage size). Set by `step!` and `Trbdf2.step!`; read by StarMaterial's source."
const FULL_DT = Ref(0.0)
"Factor on dtime after a Newton failure in `step!` (TITAN: 0.05; S5 star runs: 0.5, so that dt hovers near the largest step Newton can take instead of dropping 20×)."
const FAIL_TFAC = Ref(5.0e-2)
"S6a: accretion rate [g/s] of the prescribed self-similar mass grid (see TitanResidFast.ALE_MDOT; 0 = off). Lives here so that TitanCore and TitanResidFast share it."
const ALE_MDOT = Ref(0.0)
"S6a: per-interface weights ψ_k ∈ [0,1] of the prescribed mass grid (index k of the interface; empty = all 1): m_k(t+Δt) = m_k (1 + ψ_k Ṁ Δt/M). ψ_k = 0: Lagrangian interface."
const ALE_PSI = Ref(Float64[])
"S6a-fix (tests): false = the prescribed mass grid moves through the matter without inflow at the surface (relabelling only; the surface interface is Lagrangian)."
const ALE_INFLOW = Ref(true)
"""
S6a-fix: `CONV_STALL[]` (default 0 = off, TITAN): Newton is also accepted when the largest correction is ≤ CONV_STALL[] and repeats itself
(|dmax − dmax two iterations earlier| < 1e-3 dmax: a 2-cycle of a non-differentiable switch, e.g. the MLT flux switching on and off at the edge of a
convective zone, with corrections of 1e-7…1e-6 in ln E_r that alternate in sign forever).
"""
const CONV_STALL = Ref(0.0)
"""
S7c-1: `SAMEDT_RETRY[]` (default false = TITAN): after the FIRST Newton failure of a step, retry once at the same dt from the old state (`reset` restores it: no extrapolated radii, no predicted T, E) before the time step
is reduced. The extrapolated initial guess of the interface radii (EXTRAP_FRAC, `predict!`) fails where the surface layers of a luminous, extended star oscillate (first Newton correction ×14 in the H-ionisation zone) and the
dt reduction then repeated at every step (dt → seconds: "the S6b stall"); from the old state the same dt converges.
"""
const SAMEDT_RETRY = Ref(false)
"""
S6b: `INNER_GHOST_GEOM[]` (default false = TITAN): with the Lagrangian inner boundary (llibc = 1) the ghost interface below the first zone is the
geometric reflection r(a−1) = r(a)²/r(a+1) instead of the radius of an equal-volume zone r(a−1)³ = r(a)³ − 3Δm/ρ (which has no real root
when the zone is wider than a factor 2^(1/3) in r). Used by the mass sink of popiii_m2c.jl.
"""
const INNER_GHOST_GEOM = Ref(false)
"Julia version of subroutine step (step.f), same control flow and arithmetic. (Timed by Perf as `:integrator`.)"
step!() = Perf.@section(:integrator, step_inner!())
function step_inner!()
    r(:peq)
    s!(:jtry, 1); s!(:jback, 0)
    @label two
    FULL_DT[] = g(:dtime)
    s!(:iter, 1)
    dm1 = Inf; dm2 = Inf
    @label three
    r(:matgen); r(:matslv); r(:maxdel); r(:update); Perf.count!(:newton_iter)
    g(:dmax) <= g(:conv) && @goto seven
    (g(:dmax) <= CONV_STALL[] && g(:iter) >= 4 && abs(g(:dmax) - dm2) < 1e-3 * g(:dmax)) && @goto seven
    dm2 = dm1; dm1 = g(:dmax)
    s!(:iter, g(:iter) + 1)
    g(:iter) <= g(:niter) && @goto three
    # failure to converge: reduce timestep
    if SAMEDT_RETRY[] && g(:jtry) == 1
        s!(:jtry, 2); r(:reset); FIX_CLOCK[] && s!(:timen, g(:timeo) + g(:dtime))
        @goto two
    end
    msg(1); Perf.count!(:step_fail)
    s!(:tfac, FAIL_TFAC[])
    s!(:dtime, g(:tfac) * g(:dtime))
    s!(:timen, g(:timeo) + g(:dtime))
    s!(:jtry, g(:jtry) + 1)
    g(:jtry) > g(:ntry) && msg(2)          # Fortran writes, writout, clzfil, STOP
    r(:reset)
    FIX_CLOCK[] && s!(:timen, g(:timeo) + g(:dtime))
    @goto two
    @label seven
    s!(:lboos, 1)
    r(:matgen); r(:matslv); r(:maxdel); r(:update)
    s!(:lboos, 0)
    STAGE_HOOK[] === nothing || STAGE_HOOK[](:single, g(:dtime))
    FIX_ENERGY[] && r(:totnrg)              # book with the dt actually taken
    r(:timstp)
    xs = g(:smax) / g(:stol)
    xchek = 1.0e6
    if xs < xchek
        FIX_ENERGY[] || r(:totnrg)
        return nothing
    end
    msg(3, xs)
    s!(:jback, g(:jback) + 1)
    g(:jback) > g(:nback) && msg(4)
    msg(5)
    r(:reset)                               # (restores tes/tew/tel/teq, so an early totnrg is undone)
    FIX_CLOCK[] && s!(:timen, g(:timeo) + g(:dtime))
    @goto two
end

const LDONE = Ref{Ptr{Int32}}(C_NULL)
ldone() = (LDONE[] == C_NULL && (LDONE[] = Ptr{Int32}(fptr(:portctl_))); LDONE[])

"The step implementation used by titan_step! (TITAN's by default; e.g. DampedNewton.step!)."
const STEPFN = Ref{Function}(step!)

"Julia version of titan_step (one iteration of the titan.f loop). 0 = continue."
function titan_step!()
    st = unsafe_load(ldone()); st != 0 && return Int(st)
    STEPFN[]()
    fcall(:archive_); fcall(:writout_); fcall(:writsty_)
    if g(:timen) >= g(:tmax)
        fcall(:clzfil_); unsafe_store!(ldone(), 1); return 1
    end
    s!(:jstep, g(:jstep) + 1)
    if g(:jstep) > g(:jstepe)
        fcall(:clzfil_); unsafe_store!(ldone(), 2); return 2
    end
    0
end

"""
    run_julia!(; output = nothing, maxsteps = typemax(Int))

Pure-Julia main loop (titan.f without archive/writout/writsty/clzfil): one time step, then
`output(n)`, then TITAN's termination tests (timen ≥ tmax, or jstep > jstepe). No library calls
if every routine in ROUTINES, EOS/OPAC and TitanCore.EDDFAC is a Julia function.
"""
function run_julia!(; output = nothing, maxsteps = typemax(Int))
    n = 0
    while n < maxsteps
        STEPFN[](); n += 1
        output === nothing || output(n)
        g(:timen) >= g(:tmax) && return (n, :tmax)
        s!(:jstep, g(:jstep) + 1)
        g(:jstep) > g(:jstepe) && return (n, :jstepe)
    end
    (n, :maxsteps)
end

"Run to completion with the Julia driver; `callback(n)` after every step."
function run!(; callback = nothing)
    n = 0
    while true
        st = titan_step!(); n += 1
        callback === nothing || callback(n)
        st != 0 && return n
    end
end

# --------------------------------------------- COMMON snapshot / restore
const ALLBLOCKS = sort(collect(keys(TF.BLOCKS)))
snapshot() = TF.BACKEND[] === :julia ? Dict{Symbol,Array}(n => copy(a) for (n, a) in TF.JSTATE) :
                                       Dict(b => copy(TF.blockbytes(b)) for b in ALLBLOCKS)
function restore!(snap)
    if TF.BACKEND[] === :julia
        for (n, a) in snap; copyto!(TF.JSTATE[n], a); end
    else
        for b in ALLBLOCKS; copyto!(TF.blockbytes(b), snap[b]); end
    end
    nothing
end
push!(TF.ON_SWITCH, () -> empty!(SV))

# ------------------------------------------------------- deep-level hooks
const HOOKFUN = Dict{Int,Function}()
const HOOKPTR = Dict{Int,Ptr{Cvoid}}()
_hook_eos()  = (HOOKFUN[1](); nothing)
_hook_opac() = (HOOKFUN[2](); nothing)

"""
    install_hooks!(; eos = f_or_nothing, opac = f_or_nothing)

`f()` is called in place of the Fortran routine whenever any Fortran code calls
eos/opac.  Pass `nothing` to restore the Fortran routine.
"""
function install_hooks!(; eos = nothing, opac = nothing)
    for (which, f, jf) in ((1, eos, _hook_eos), (2, opac, _hook_opac))
        if f === nothing
            ccall(fptr(:titan_set_hook), Cvoid, (Cint, Ptr{Cvoid}), which, C_NULL)
        else
            HOOKFUN[which] = f
            p = which == 1 ? @cfunction(_hook_eos, Cvoid, ()) : @cfunction(_hook_opac, Cvoid, ())
            HOOKPTR[which] = p
            ccall(fptr(:titan_set_hook), Cvoid, (Cint, Ptr{Cvoid}), which, p)
        end
    end
end
fortran_eos()  = fcall(:eos_f_)
fortran_opac() = fcall(:opac_f_)

end # module
