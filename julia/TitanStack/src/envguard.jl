# Load-time environment guard of TitanStack.
#
# A precompiled image freezes everything that was computed when its modules were loaded: every `const X = get(ENV, ...)`, every
# conditional `include`, every default that is read once. This file finds those environment dependences mechanically from the
# sources, at precompile time, and records the values the image was built with (BAKED). `check()` compares them with the
# environment of the process that loads the image.
#
# Method (see `scan`): every source file the package includes (the module files, and the files they include in turn) is parsed.
#  1. The load-time code is every expression that is not inside a function, macro or lambda body (module bodies, `const`, `if`,
#     `for`, `let`, macro arguments, string interpolations, ...).
#  2. In it, every use of `ENV` yields a key: `get(ENV, "K", ..)`, `haskey(ENV, "K")`, `ENV["K"]`, `get!(ENV, "K", ..)`, ... An `ENV`
#     use without a literal key (a variable key, the whole ENV) makes the scan fail: the keys cannot be known.
#  3. Every identifier in the load-time code that names a function defined in the scanned files adds the keys read in that function's
#     body, and (transitively) in the bodies of the functions it mentions. Matching is by name only (over-strict: all methods and
#     all modules of that name count). `__init__` bodies are not load-time code of the image (they run at every load).
#  4. Other things that fix a value in the image: `homedir()` (key HOME), `Threads.nthreads()` (pseudo-key #nthreads), `pwd()`
#     (#pwd). Calls of time/getpid/gethostname/rand/tempname/ccall in load-time code are reported by `scan` and fail it unless they
#     are listed in ALLOWED_IMPURE (the reviewed ones).
# Being over-strict is acceptable (a run is sent to the include path); missing a key is not.

const ENVREADERS = (:get, :haskey, :get!, :getindex, :pop!, :delete!, :setindex!, :in, :getkey, :haskey)

# `f(args) = body`, `function f ... end`, `macro m ... end`
function _isfuncdef(ex)
    ex isa Expr || return false
    ex.head in (:function, :macro) && return true
    if ex.head === :(=) && length(ex.args) == 2
        l = ex.args[1]
        while l isa Expr && l.head in (:where, :(::)); l = l.args[1]; end
        return l isa Expr && l.head === :call
    end
    false
end
function _funcname(ex)::Union{Symbol,Nothing}
    sig = ex.args[1]
    while sig isa Expr && sig.head in (:where, :(::)); sig = sig.args[1]; end
    sig isa Expr && sig.head === :call || return nothing
    f = sig.args[1]
    while f isa Expr && f.head in (:curly, :(::)); f = f.args[1]; end
    f isa Expr && f.head === :. && (f = f.args[end]; f isa QuoteNode && (f = f.value))
    f isa Symbol ? f : (f isa Expr && f.head === :. ? nothing : nothing)
end

mutable struct ScanState
    loadkeys::Dict{String,Set{String}}     # key -> where (file:line) it is read in load-time code
    fkeys::Dict{Tuple{String,Symbol},Set{String}}   # (file, function name) -> keys read directly in its bodies
    fsyms::Dict{Tuple{String,Symbol},Set{Symbol}}   # (file, function name) -> identifiers in its bodies
    loadsyms::Set{Tuple{String,Symbol}}             # (file, identifier) in the load-time code
    loadsyms_where::Dict{Tuple{String,Symbol},String}
    initkeys::Dict{String,Set{String}}     # file -> keys read in its __init__ bodies (re-read at every load)
    initsyms::Dict{String,Set{Symbol}}
    includes::Vector{String}               # nested include() files found (absolute)
    problems::Vector{String}
    runtime::Set{String}                   # files whose function and lambda bodies are run-time only (see `scan`)
end
ScanState() = ScanState(Dict{String,Set{String}}(), Dict{Tuple{String,Symbol},Set{String}}(), Dict{Tuple{String,Symbol},Set{Symbol}}(), Set{Tuple{String,Symbol}}(), Dict{Tuple{String,Symbol},String}(), Dict{String,Set{String}}(), Dict{String,Set{Symbol}}(), String[], String[], Set{String}())

# collect identifiers and ENV keys in `ex` (a body: nothing is skipped). `keys` receives (key, line) pairs.
function _collect!(ex, syms::Set{Symbol}, keys::Vector{Tuple{String,String}}, probs::Vector{String}, where_::String)
    if ex isa Symbol
        push!(syms, ex)
    elseif ex isa QuoteNode
        _collect!(ex.value, syms, keys, probs, where_)
    elseif ex isa Expr
        if ex.head === :ref && length(ex.args) >= 2 && ex.args[1] === :ENV
            k = ex.args[2]
            k isa String ? push!(keys, (k, where_)) : push!(probs, "ENV[...] with a non-literal key at $where_: $(ex)")
        elseif ex.head === :call && any(a -> a === :ENV, ex.args[2:end])
            i = findfirst(a -> a === :ENV, ex.args[2:end]) + 1
            if i < length(ex.args)
                k = ex.args[i+1]
                k isa String ? push!(keys, (k, where_)) : push!(probs, "$(ex.args[1])(ENV, <non-literal key>) at $where_: $(ex)")
            else
                push!(probs, "ENV passed whole at $where_: $(ex)")
            end
        elseif ex.head === :call && ex.args[1] === :ENV
            push!(probs, "ENV called at $where_")
        elseif ex.head === :. && length(ex.args) >= 1 && ex.args[1] === :ENV
            push!(probs, "ENV property at $where_: $(ex)")
        elseif ex.head === :for || ex.head === :generator || ex.head === :comprehension
            # `for x in ENV` or iteration over ENV
            any(a -> a === :ENV, ex.args) && push!(probs, "iteration over ENV at $where_")
        end
        for a in ex.args
            if a === :ENV && ex.head !== :ref && !(ex.head === :call)
                ex.head in (:for, :generator, :comprehension, :.) || push!(probs, "bare ENV (head $(ex.head)) at $where_: $(ex)")
            end
            _collect!(a, syms, keys, probs, where_)
        end
    end
    nothing
end

const LINE_RE = Ref{Int}(0)

# walk the load-time code of one file; function bodies are recorded separately
function _walk!(S::ScanState, ex, file::String, line::Int; allow_dynamic_include::Bool = false)
    if ex isa LineNumberNode
        return ex.line
    end
    if ex isa QuoteNode
        ex.value isa Symbol && (push!(S.loadsyms, (file, ex.value)); get!(S.loadsyms_where, (file, ex.value), "$file:$line"))
        return line
    end
    ex isa Expr || (ex isa Symbol && (push!(S.loadsyms, (file, ex)); get!(S.loadsyms_where, (file, ex), "$file:$line")); return line)
    where_ = "$(basename(file)):$line"
    if ex.head === :function && length(ex.args) == 1
        return line                                        # `function f end`: a declaration without a body
    elseif _isfuncdef(ex) && basename(file) in S.runtime
        return line                                        # a run-time-only file (`scan`'s `runtime`): its function bodies are not scanned
    elseif _isfuncdef(ex)
        nm = _funcname(ex)
        syms = Set{Symbol}(); ks = Tuple{String,String}[]
        # the signature's default arguments are evaluated per call, not at load: the whole definition is a body
        _collect!(ex, syms, ks, S.problems, where_)
        if nm === nothing
            # a method on a callable object or a qualified name: attribute to every identifier of the signature (over-strict)
            nm = :__anonymous_method__
        end
        if nm === :__init__
            union!(get!(S.initkeys, basename(file), Set{String}()), (k for (k, _) in ks)); union!(get!(S.initsyms, basename(file), Set{Symbol}()), syms)
        else
            union!(get!(S.fkeys, (file, nm), Set{String}()), (k for (k, _) in ks))
            union!(get!(S.fsyms, (file, nm), Set{Symbol}()), syms)
        end
        return line
    elseif (ex.head === :-> || ex.head === :quote) && basename(file) in S.runtime
        return line
    elseif ex.head === :-> || ex.head === :quote
        # lambda or quoted code: not executed at load; its keys count as the keys of an anonymous function that a load-time call may reach
        syms = Set{Symbol}(); ks = Tuple{String,String}[]
        _collect!(ex, syms, ks, S.problems, where_)
        union!(get!(S.fkeys, (file, :__lambda__), Set{String}()), (k for (k, _) in ks))
        union!(get!(S.fsyms, (file, :__lambda__), Set{Symbol}()), syms)
        return line
    elseif ex.head === :macro
        return line
    elseif ex.head === :call && ex.args[1] === :include
        a = length(ex.args) >= 2 ? ex.args[end] : nothing
        p = _static_path(a, file)
        if p === nothing
            allow_dynamic_include || push!(S.problems, "include with a non-static path at $where_: $ex")
        else
            push!(S.includes, p)
        end
        return line
    end
    # a load-time expression: its own ENV uses and identifiers (the sub-expressions that are function bodies are handled by recursion below)
    if ex.head === :module
        for a in ex.args[3:end]; line = _walk!(S, a, file, line; allow_dynamic_include); end
        return line
    elseif ex.head === :struct
        # the type's own name (and its constructors' names) are definitions, not load-time uses; field types are identifiers of types only
        for blk in ex.args[3:end]
            if blk isa Expr && blk.head === :block
                for st in blk.args
                    if st isa Symbol; continue                                            # an untyped field name
                    elseif st isa Expr && st.head === :(::) && length(st.args) == 2 && st.args[1] isa Symbol
                        line = _walk!(S, st.args[2], file, line; allow_dynamic_include)    # `name::Type`: the field name is not a use
                    else
                        line = _walk!(S, st, file, line; allow_dynamic_include)
                    end
                end
            else
                line = _walk!(S, blk, file, line; allow_dynamic_include)
            end
        end
        return line
    end
    # ENV uses directly in this node
    ks = Tuple{String,String}[]; syms = Set{Symbol}()
    if (ex.head === :ref && length(ex.args) >= 2 && ex.args[1] === :ENV) ||
       (ex.head === :call && any(a -> a === :ENV, ex.args[2:end])) || ex.head === :for
        _collect_node!(ex, ks, S.problems, where_)
    end
    for (k, w) in ks; push!(get!(S.loadkeys, k, Set{String}()), w); end
    for a in ex.args; line = _walk!(S, a, file, line; allow_dynamic_include); end
    line
end
# the ENV keys of this node only (not of its children)
function _collect_node!(ex::Expr, keys, probs, where_)
    if ex.head === :ref
        k = ex.args[2]; k isa String ? push!(keys, (k, where_)) : push!(probs, "ENV[...] with a non-literal key at $where_: $(ex)")
    elseif ex.head === :call
        i = findfirst(a -> a === :ENV, ex.args[2:end]) + 1
        if i < length(ex.args)
            k = ex.args[i+1]; k isa String ? push!(keys, (k, where_)) : push!(probs, "$(ex.args[1])(ENV, <non-literal key>) at $where_: $(ex)")
        else
            push!(probs, "ENV passed whole at $where_: $(ex)")
        end
    elseif ex.head === :for
        any(a -> a === :ENV, ex.args[1].args) && push!(probs, "iteration over ENV at $where_")
    end
end
# `joinpath(@__DIR__, "a", "b.jl")` or a string literal -> absolute path; anything else -> nothing
function _static_path(a, file)
    dir = dirname(file)
    a isa String && return normpath(dir, a)
    if a isa Expr && a.head === :call && a.args[1] === :joinpath
        parts = String[]
        for (i, p) in enumerate(a.args[2:end])
            if p isa String; push!(parts, p)
            elseif p isa Expr && p.head === :macrocall && p.args[1] === Symbol("@__DIR__") && i == 1; push!(parts, dir)
            else return nothing end
        end
        return normpath(parts...)
    end
    nothing
end

# the scan: parse `roots` and everything they include. Returns (keys, report)
function scan(roots::Vector{String}; dynamic_include_ok::Set{String} = Set{String}(), skip::Set{String} = Set{String}(), reread = Tuple{String,String,Union{String,Nothing}}[], runtime::Set{String} = Set{String}())
    S = ScanState(); S.runtime = runtime; seen = Set{String}(); queue = copy(roots); files = String[]
    while !isempty(queue)
        f = popfirst!(queue); f in seen && continue; basename(f) in skip && continue; push!(seen, f); push!(files, f)
        isfile(f) || (push!(S.problems, "missing file $f"); continue)
        ast = Meta.parseall(read(f, String); filename = f)
        nbefore = length(S.includes)
        _walk!(S, ast, f, 0; allow_dynamic_include = basename(f) in dynamic_include_ok)
        append!(queue, S.includes[nbefore+1:end])
    end
    keys = Set{String}(Base.keys(S.loadkeys))
    # transitive closure over functions mentioned in load-time code
    # Functions are matched by name; a name that the referencing file defines itself resolves to that file's definitions only (names
    # like `table` or `load` exist in many modules), otherwise to all files that define it (imports across modules).
    defs = Dict{Symbol,Vector{String}}()
    for (f, n) in Base.keys(S.fkeys); f in get!(defs, n, String[]) || push!(defs[n], f); end
    for (f, n) in Base.keys(S.fsyms); f in get!(defs, n, String[]) || push!(defs[n], f); end
    resolve(f, n) = haskey(defs, n) ? (f in defs[n] ? [(f, n)] : [(g, n) for g in defs[n]]) : Tuple{String,Symbol}[]
    reached = Set{Tuple{String,Symbol}}(); work = Tuple{String,Symbol}[]; via = Dict{Tuple{String,Symbol},String}()
    for (f, s) in S.loadsyms, d in resolve(f, s)
        push!(work, d); haskey(via, d) || (via[d] = "load-time code at " * S.loadsyms_where[(f, s)])
    end
    # lambdas (and anonymous methods) are always counted (over-strict)
    for f in files, n in (:__lambda__, :__anonymous_method__)
        haskey(S.fsyms, (f, n)) && (push!(work, (f, n)); via[(f, n)] = "always")
    end
    why = Dict{String,Set{String}}(k => copy(w) for (k, w) in S.loadkeys)
    while !isempty(work)
        d = pop!(work); d in reached && continue; push!(reached, d)
        for k in get(S.fkeys, d, ()); push!(keys, k); push!(get!(why, k, Set{String}()), "fn $(d[2]) in $(basename(d[1])) (via $(via[d]))"); end
        for t in get(S.fsyms, d, ()), e in resolve(d[1], t)
            e in reached || (push!(work, e); haskey(via, e) || (via[e] = "fn $(d[2]) <- " * via[d]))
        end
    end
    # Value-fixing calls in load-time code itself. `homedir()` reads HOME (a guarded key). Thread counts and ids fix the size of
    # per-thread tables to the BUILD process's single thread: a use is accepted only at a reviewed site (`reread`: file, symbol, ENV key
    # that feeds it) whose module recomputes it in an `__init__` in the same file, which is checked here. Anything else fails the scan.
    impure = Dict{Symbol,String}()
    IMPURE = (:homedir, :nthreads, :maxthreadid, :threadid, :nthreadpools, :threadpoolsize, :pwd, :gethostname, :getpid, :time, :time_ns, :rand, :randn,
              :tempname, :mktemp, :mktempdir, :ccall, :cglobal, :dlopen, :dlsym, :now, :CPU_THREADS, :nprocs, :myid)
    for (f, sym) in S.loadsyms
        sym in IMPURE || continue
        w = S.loadsyms_where[(f, sym)]; impure[sym] = w
        if sym === :homedir
            push!(keys, "HOME")
        elseif (i = findfirst(r -> r[1] == basename(f) && Symbol(r[2]) === sym, reread)) !== nothing
            k = reread[i][3]
            if sym in get(S.initsyms, basename(f), Set{Symbol}()) && (k === nothing || k in get(S.initkeys, basename(f), Set{String}()))
                k === nothing || delete!(keys, k)
            else
                push!(S.problems, "$sym at $w is listed as recomputed in an __init__ of $(basename(f)), but no __init__ there reads it")
            end
        else
            push!(S.problems, "load-time use of $sym at $w: its value would be frozen in the image (guard it, or recompute it in an __init__ and list the site in `reread`)")
        end
    end
    # also report the ones that only occur in functions the name closure reached (information only: not load-time uses)
    for r in reached, sym in get(S.fsyms, r, ())
        sym in IMPURE && !haskey(impure, sym) && (impure[sym] = "(in function $(r[2]) of $(basename(r[1])), reached by name)")
    end
    (; keys = sort!(collect(keys)), why, files, impure, problems = S.problems, reached = sort!(collect(string(basename(f), ":", n) for (f, n) in reached)))
end


# ---- the guard -----------------------------------------------------------------------------------------------------------------
"The value the process has for a guarded key now (nothing: absent)."
current(k::AbstractString) = k == "#nthreads" ? string(Threads.nthreads()) : get(ENV, k, nothing)

"The keys whose current value differs from the one the image was built with: (key, baked, current) triples."
function mismatches(baked)
    bad = Tuple{String,Union{String,Nothing},Union{String,Nothing}}[]
    for (k, v) in baked
        c = current(k)
        c == v || push!(bad, (k, v, c))
    end
    bad
end
show_val(v) = v === nothing ? "(unset)" : repr(v)
function mismatch_message(bad, nkeys)
    io = IOBuffer()
    println(io, "TitanStack: the precompiled image was built with a different load-time environment than this process has ($(length(bad)) of $nkeys guarded keys differ):")
    for (k, b, c) in bad; println(io, "  ", k, ": image ", show_val(b), ", environment ", show_val(c)); end
    print(io, "The image froze these values when it was built; running on it would silently use the image's. Rebuild it with the wanted settings (julia/TitanStack/build.jl) or use the include path (unset TITAN_PKG).")
    String(take!(io))
end
