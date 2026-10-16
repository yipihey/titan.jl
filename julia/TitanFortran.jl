"""
    TitanFortran

Phase-1 Julia driver for TITAN 1.2 (Gehmeyr & Mihalas; UIUC licence, notice
retained in ../fortran/README).  Loads `libtitan.so`, drives it through the
entry points `titan_init_/titan_step_/titan_finalize_` and maps every COMMON
block of `titan.com` onto Julia arrays (no copies) using offsets computed by
parsing `titan.com` + `titan.par` (implicit typing from `titan.imp`:
i-n -> Int32, otherwise Float64; `header` is character*80).

RE-ENTRANCY: all solver state lives in global COMMON blocks and Fortran I/O
units, so there is exactly ONE model per process.  A Fortran STOP (e.g.
non-convergence 'step01') terminates the Julia process too.
"""
module TitanFortran

using Libdl, Serialization

const ROOT    = normpath(joinpath(@__DIR__, ".."))
# run-output directories are git-ignored; create them so a fresh checkout works
for d in ("jl_runs", "p2_runs", "p3_runs", "p4_runs", "p4_legs"); mkpath(joinpath(ROOT, "lib", d)); end
const FSRC    = joinpath(ROOT, "fortran")
const LIBPATH = Ref(joinpath(ROOT, "lib", "libtitan.so"))
const HANDLE  = Ref{Ptr{Cvoid}}(C_NULL)
const LIBTITAN = joinpath(ROOT, "lib", "libtitan.so")   # for literal ccall tuples

# --------------------------------------------------------------------------
# COMMON-block layout from titan.par / titan.com
# --------------------------------------------------------------------------
struct Var
    name::Symbol
    T::DataType          # Float64, Int32 or NTuple{80,UInt8} (character*80)
    dims::Vector{Int}
    offset::Int          # bytes from the start of the block
end
nbytes(v::Var) = sizeof(v.T) * prod(v.dims; init = 1)

"Read a fixed-form include file, drop comments, join continuation lines."
function fixed_form_statements(path)
    stmts = String[]
    for raw in eachline(path)
        isempty(raw) && continue
        c = raw[1]
        (c == 'c' || c == 'C' || c == '*' || c == '!') && continue
        line = replace(raw, r"!.*$" => "")
        if length(line) >= 6 && !isspace(line[6]) && line[6] != '0' && !isempty(stmts)
            stmts[end] *= " " * strip(line[7:end])
        else
            push!(stmts, strip(line))
        end
    end
    filter(!isempty, stmts)
end

"Split on commas that are not inside parentheses."
function split_top(s)
    out = String[]; depth = 0; buf = IOBuffer()
    for ch in s
        if ch == '(' ; depth += 1 elseif ch == ')' ; depth -= 1 end
        if ch == ',' && depth == 0
            push!(out, strip(String(take!(buf))))
        else
            write(buf, ch)
        end
    end
    push!(out, strip(String(take!(buf))))
    filter(!isempty, out)
end

function int_eval(ex, P::Dict{Symbol,Any})
    ex isa Integer && return Int(ex)
    ex isa Symbol && return Int(P[ex])
    if ex isa Expr && ex.head == :call
        a = [int_eval(x, P) for x in ex.args[2:end]]
        f = ex.args[1]
        f == :+ && return sum(a)
        f == :- && return length(a) == 1 ? -a[1] : a[1] - a[2]
        f == :* && return prod(a)
        f == :÷ && return a[1] ÷ a[2]
    end
    error("cannot evaluate dimension expression $ex")
end

"Integer PARAMETERs of titan.par."
function parse_params(path = joinpath(FSRC, "titan.par"))
    P = Dict{Symbol,Any}()
    for st in fixed_form_statements(path)
        m = match(r"^parameter\s*\((.*)\)\s*$"i, st)
        m === nothing && continue
        for item in split_top(m.captures[1])
            k, v = strip.(split(item, '='; limit = 2))
            vs = lowercase(v)
            if occursin('.', vs)            # real*8 parameter (D exponent)
                P[Symbol(lowercase(k))] = parse(Float64, replace(vs, 'd' => 'e'))
            else
                P[Symbol(lowercase(k))] = int_eval(Meta.parse(vs), P)
            end
        end
    end
    P
end

"Block name => Vector{Var} (declaration order, packed offsets)."
function parse_commons(path = joinpath(FSRC, "titan.com"), P = parse_params())
    chartypes = Dict{Symbol,DataType}()
    blocks = Dict{Symbol,Vector{Var}}()
    for st in fixed_form_statements(path)
        s = lowercase(st)
        m = match(r"^character\s*\*\s*(\d+)\s+(.*)$", s)
        if m !== nothing
            n = parse(Int, m.captures[1]); @assert n == 80
            for nm in split_top(m.captures[2]); chartypes[Symbol(nm)] = NTuple{80,UInt8}; end
            continue
        end
        m = match(r"^common\s*/\s*(\w+)\s*/\s*(.*)$", s)
        m === nothing && continue
        blk = Symbol(m.captures[1]); vars = Var[]; off = 0
        for item in split_top(m.captures[2])
            mm = match(r"^(\w+)\s*(?:\((.*)\))?$", item)
            name = Symbol(mm.captures[1])
            dims = mm.captures[2] === nothing ? Int[] :
                   [int_eval(Meta.parse(d), P) for d in split_top(mm.captures[2])]
            T = get(chartypes, name, first(String(name)) in 'i':'n' ? Int32 : Float64)
            v = Var(name, T, dims, off); push!(vars, v); off += nbytes(v)
        end
        blocks[blk] = vars
    end
    blocks
end

const PARAMS  = parse_params()
const BLOCKS  = parse_commons()
const VARMAP  = Dict{Symbol,Tuple{Symbol,Var}}(v.name => (b, v) for (b, vs) in BLOCKS for v in vs)
const MGR     = PARAMS[:mgr]
block_size(b::Symbol) = (v = BLOCKS[b][end]; v.offset + nbytes(v))

# --------------------------------------------------------------------------
# library handle and raw access
# --------------------------------------------------------------------------
function load(path = LIBPATH[])
    HANDLE[] == C_NULL || return HANDLE[]
    LIBPATH[] = path
    # RTLD_DEEPBIND: resolve libtitan's own undefined symbols (pow, exp, log …)
    # against ITS dependencies (glibc libm) first.  Without it, the binding
    # depended on load order: if Julia's libopenlibm was already in the global
    # scope, libtitan's pow/exp/log bound to OpenLibm instead of glibc and the
    # trajectory changed in the last bits (found in phase 3, see PORTING_NOTES).
    HANDLE[] = Libdl.dlopen(path, Libdl.RTLD_GLOBAL | Libdl.RTLD_NOW | Libdl.RTLD_DEEPBIND)
end
# load (and bind) before any ccall can trigger a plain dlopen; TITAN_NOLIB=1 runs without libtitan (julia backend only)
__init__() = get(ENV, "TITAN_NOLIB", "0") == "1" || load()
"Address of COMMON block `blk` (gfortran symbol `blk_`)."
blockptr(blk::Symbol) = Ptr{UInt8}(Libdl.dlsym(load(), Symbol(blk, "_")))

"Zero-copy Julia view of a COMMON variable (array or 0-dim)."
# --------------------------------------------------------------------------
# state backends: :common (Fortran COMMON blocks in libtitan) or :julia (Julia-owned arrays,
# same names and shapes). Everything in the port reaches state through `var`, so the routines
# run unchanged on either backend.
# --------------------------------------------------------------------------
const BACKEND = Ref(:common)
const JSTATE = Dict{Symbol,Array}()
"Functions run after a backend switch (modules register their view-cache clearers here)."
const ON_SWITCH = Function[]
_switched() = foreach(f -> f(), ON_SWITCH)

"Copy every COMMON variable into Julia-owned arrays and switch to the :julia backend."
function to_julia!()
    BACKEND[] === :julia && return nothing
    empty!(JSTATE)
    for (_, vs) in BLOCKS, v in vs
        JSTATE[v.name] = copy(var(v.name))
    end
    BACKEND[] = :julia; _switched(); nothing
end
"Copy the Julia state back into the COMMON blocks and switch to :common (for Fortran output or comparisons)."
function to_common!()
    BACKEND[] === :common && return nothing
    BACKEND[] = :common
    for (n, a) in JSTATE; copyto!(var(n), a); end
    _switched(); nothing
end
"""
    save_state(path) / load_state!(path)

Complete model state (every COMMON variable) in a simple, Julia-version-independent binary format
("TITANST1": per variable name, element type, dims, raw little-endian data). A path ending in
`.gz` is written/read through the system `gzip` (initial states are mostly zeros: ~5.6 MB → ~0.1 MB).
`load_state!` switches to the :julia backend (no libtitan needed).
"""
function save_state(path)
    st = BACKEND[] === :julia ? JSTATE : Dict{Symbol,Array}(v.name => copy(var(v.name)) for (_, vs) in BLOCKS for v in vs)
    buf = IOBuffer()
    write(buf, b"TITANST1"); write(buf, Int64(length(st)))
    for n in sort(collect(keys(st)))
        a = st[n]; nb = Vector{UInt8}(String(n))
        code = eltype(a) === Float64 ? 0x01 : eltype(a) === Int32 ? 0x02 : eltype(a) === NTuple{80,UInt8} ? 0x03 : error("type $(eltype(a))")
        write(buf, Int64(length(nb)), nb, code, Int64(ndims(a))); for d in size(a); write(buf, Int64(d)); end
        write(buf, reinterpret(UInt8, vec(collect(a))))
    end
    data = take!(buf)
    if endswith(path, ".gz")
        open(pipeline(`gzip -9 -c`, stdout = path), "w") do io; write(io, data); end
    else
        write(path, data)
    end
    path
end
function load_state!(path)
    data = endswith(path, ".gz") ? read(`gzip -dc $path`) : read(path)
    io = IOBuffer(data)
    String(read(io, 8)) == "TITANST1" || error("load_state!: $path is not a TITANST1 file")
    n = read(io, Int64); st = Dict{Symbol,Array}()
    for _ in 1:n
        name = Symbol(String(read(io, read(io, Int64)))); code = read(io, UInt8)
        dims = Tuple(read(io, Int64) for _ in 1:read(io, Int64))
        T = code == 0x01 ? Float64 : code == 0x02 ? Int32 : NTuple{80,UInt8}
        raw = read(io, prod(dims; init = 1) * sizeof(T))
        st[name] = reshape(collect(reinterpret(T, raw)), dims)
    end
    for (_, vs) in BLOCKS, v in vs
        haskey(st, v.name) || error("load_state!: $(v.name) missing in $path")
        size(st[v.name]) == (isempty(v.dims) ? () : Tuple(v.dims)) || error("load_state!: $(v.name) has wrong shape")
    end
    empty!(JSTATE); merge!(JSTATE, st)
    BACKEND[] = :julia; _switched(); nothing
end

function var(name::Symbol)
    BACKEND[] === :julia && return JSTATE[name]
    blk, v = VARMAP[name]
    p = blockptr(blk) + v.offset
    isempty(v.dims) ? unsafe_wrap(Array, Ptr{v.T}(p), ()) : unsafe_wrap(Array, Ptr{v.T}(p), Tuple(v.dims))
end
getscalar(name::Symbol) = var(name)[]
setscalar!(name::Symbol, x) = (var(name)[] = x; x)
"Raw bytes of a whole COMMON block (zero-copy)."
blockbytes(blk::Symbol) = unsafe_wrap(Array, blockptr(blk), block_size(blk))

"Compare parsed block sizes with the ELF symbol sizes (nm -S)."
function check_layout(; nm = `nm -S --defined-only $(LIBPATH[])`)
    sizes = Dict{Symbol,Int}()
    for l in eachline(nm)
        f = split(l); length(f) == 4 || continue
        endswith(f[4], "_") && (sizes[Symbol(f[4][1:end-1])] = parse(Int, f[2]; base = 16))
    end
    res = Tuple{Symbol,Int,Int,Bool}[]
    for b in sort(collect(keys(BLOCKS)))
        p = block_size(b); e = sizes[b]
        push!(res, (b, p, e, p == e || (e - p) in 1:7 && e % 8 == 0))   # trailing pad to 8
    end
    res
end

# --------------------------------------------------------------------------
# driving the model
# --------------------------------------------------------------------------
const RUNDIR = Ref{String}("")
const OLDDIR = Ref{String}("")

"""
    init(input; rundir, predump=nothing)

`input` is an input deck name in ../fortran (e.g. "titan.in08") or a problem
number 3..10 (-> "titan.inNN").  The deck is copied to `rundir/titan.in`; an
optional dump file (`predump`) is copied in for restarts (COOL, RLST leg 2).
The working directory stays in `rundir` until `finalize()`.
"""
function init(input; rundir::AbstractString, predump = nothing)
    load()
    deck = input isa Integer ? "titan.in" * lpad(input, 2, '0') : String(input)
    mkpath(rundir)
    cp(joinpath(FSRC, deck), joinpath(rundir, "titan.in"); force = true)
    predump === nothing || cp(predump, joinpath(rundir, basename(predump)); force = true)
    OLDDIR[] = pwd(); RUNDIR[] = abspath(rundir); cd(RUNDIR[])
    st = Ref{Int32}(0)
    ccall((:titan_init_, LIBTITAN), Cvoid, (Ref{Int32},), st)
    st[]
end

"One TITAN time step (step, archive, writout, writsty). Returns 0 = continue, 1/2 = finished."
function step!()
    st = Ref{Int32}(0)
    ccall((:titan_step_, LIBTITAN), Cvoid, (Ref{Int32},), st)
    st[]
end

function finalize()
    ccall((:titan_finalize_, LIBTITAN), Cvoid, ())
    isempty(OLDDIR[]) || cd(OLDDIR[])
    nothing
end

"Run to completion; `callback(nstep)` is called after each step. Returns the number of steps."
function run!(; callback = nothing, maxsteps = typemax(Int))
    n = 0
    while n < maxsteps
        s = step!(); n += 1
        callback === nothing || callback(n)
        s != 0 && break
    end
    n
end

"""
    get_state()

Copy of the current (new-time-level) solution on the active grid.
Interfaces k = ngrs..ngre+1: r, xme (exterior mass/4π), u, F.
Zones      k = ngrs..ngre  : rho, T, eg, pg, Er, fedd, chiF, kapE, kapP, cs.
"""
function get_state()
    ngrs = Int(getscalar(:ngrs)); ngre = Int(getscalar(:ngre))
    I = ngrs:ngre+1; Z = ngrs:ngre
    (; jstep = Int(getscalar(:jstep)), iter = Int(getscalar(:iter)),
       time = getscalar(:timen), dtime = getscalar(:dtime), ngrs, ngre,
       r = var(:rn)[I], xme = var(:xmen)[I], u = var(:un)[I], F = var(:frn)[I],
       rho = var(:dn)[Z], T = var(:tn)[Z], eg = var(:egn)[Z], pg = var(:pgn)[Z],
       Er = var(:ern)[Z], fedd = var(:fedd)[Z], chiF = var(:chifn)[Z],
       kapE = var(:xken)[Z], kapP = var(:xkpn)[Z], cs = var(:asn)[Z])
end

# --------------------------------------------------------------------------
# reader for the direct-access dump files (layout from archive.f)
# --------------------------------------------------------------------------
const DUMP_SCAL1 = (:cmax, :dmax, :smax, :timen, :dtime, :tfac, :tetot, :tee, :tes, :tew, :tel, :teq)
const DUMP_ARR0  = (:r0, :xm0, :d0, :u0, :egrt0, :dvol0, :tescr0)
const DUMP_SCAL2 = (:geddl, :geddr, :xmass, :xmext, :xmtot)
const DUMP_ARRS  = (:rn, :rmun, :rmup1n, :rmum1n, :avchin, :urel, :dvoln, :xmen, :dn, :dsn,
                    :chifn, :drdt, :un, :usn, :tn, :asn, :xken, :dmdt, :egn, :egsn, :pgn,
                    :xnen, :xkpn, :xnt, :ern, :ersn, :egrn, :egrsn, :plfn, :etn, :etsn,
                    :egtn, :egtsn, :fedd, :frn, :frsn, :egrtn, :egrtsn)
const LDUMP = 8 * (45 * MGR + 19)

function read_dump(file, jdump = nothing)
    open(file) do io
        if jdump === nothing
            seek(io, LDUMP); jdump = read(io, Int32)
        end
        seek(io, (jdump + 1) * LDUMP)
        d = Dict{Symbol,Any}(:jdump => read(io, Int32), :jstep => read(io, Int32))
        for n in DUMP_SCAL1; d[n] = read(io, Float64); end
        for n in DUMP_ARR0;  d[n] = read!(io, Vector{Float64}(undef, MGR)); end
        for n in DUMP_SCAL2; d[n] = read(io, Float64); end
        for n in DUMP_ARRS;  d[n] = read!(io, Vector{Float64}(undef, MGR)); end
        d
    end
end

"Bitwise equality (distinguishes ±0 and NaN payloads)."
bitsequal(a::AbstractArray{Float64}, b::AbstractArray{Float64}) =
    size(a) == size(b) && reinterpret(UInt64, vec(collect(a))) == reinterpret(UInt64, vec(collect(b)))
bitsequal(a::Float64, b::Float64) = reinterpret(UInt64, a) == reinterpret(UInt64, b)

"""
Compare the live COMMON state with a dump record bitwise; returns a vector of
(name, equal::Bool, nmismatch).  Every dumped variable has the same name in
COMMON (e.g. rn in /new/, fedd in /mid/, r0 in /energy/, xnt in /agrid/).
"""
function compare_with_dump(d)
    out = Tuple{Symbol,Bool,Int}[]
    for n in (DUMP_SCAL1..., DUMP_SCAL2...)
        push!(out, (n, bitsequal(getscalar(n), d[n]), Int(!bitsequal(getscalar(n), d[n]))))
    end
    for n in (DUMP_ARR0..., DUMP_ARRS...)
        a = reinterpret(UInt64, collect(var(n))); b = reinterpret(UInt64, d[n])
        nm = count(a .!= b); push!(out, (n, nm == 0, nm))
    end
    push!(out, (:jstep, Int(getscalar(:jstep)) == d[:jstep], 0))
    out
end

end # module
