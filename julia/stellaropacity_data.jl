# S3 data helpers shared by stellaropacity_build.jl and stellaropacity_test.jl: GN93hz and Cassisi parsers, nonuniform Hermite.
const OPALFILE = expanduser("~/Projects/popiii_lit/opal/ex/kap_input_data/opal/GN93hz")
const CONDFILE = expanduser("~/Projects/popiii_lit/opal/raw/condtabl.data")
smoothstep(x) = (x = clamp(x, 0.0, 1.0); x * x * (3 - 2x))

# ------------------------------------------------------------------ parsers
"Z = 0 tables of GN93hz: Dict X => (nT × 19 matrix of log10 κ_R, NaN where the file has 9.999 / blank), the log T and log R nodes."
function read_gn93(path)
    L = readlines(path); tabs = Dict{Float64,Matrix{Float64}}(); lT = Float64[]; lR = collect(-8.0:0.5:1.0)
    i = 242
    while i <= length(L)
        m = match(r"^TABLE #\s*(\d+)\s+.*X=([\d.]+) Y=([\d.]+) Z=([\d.]+)", L[i])
        if m === nothing; i += 1; continue; end
        X = parse(Float64, m[2]); Z = parse(Float64, m[4]); rows = Vector{Float64}[]; lt = Float64[]; j = i + 1
        while j <= length(L) && !startswith(L[j], "TABLE")
            l = L[j]
            if occursin(r"^\d\.\d\d", l)
                push!(lt, parse(Float64, l[1:4])); v = fill(NaN, 19)
                for k in 1:19
                    a = 5 + (k - 1) * 7 + 1; a > length(l) && break
                    s = strip(l[a:min(a + 6, length(l))]); isempty(s) && continue
                    x = parse(Float64, s); x >= 9.99 || (v[k] = x)
                end
                push!(rows, v)
            end
            j += 1
        end
        Z == 0.0 && (tabs[X] = Matrix(reduce(hcat, rows)'); lT = lt)
        i = j
    end
    tabs, lT, lR
end
"Cassisi et al. 2007 conduction: log10 K [erg cm⁻¹ s⁻¹ K⁻¹] for Z = 1 and 2 on (log ρ −6…11.5, log T 3…10), step 0.25."
function read_cond(path)
    L = readlines(path); out = Dict{Float64,Matrix{Float64}}()
    for z in (1.0, 2.0)
        i0 = findfirst(l -> (v = split(l); length(v) == 30 && parse(Float64, v[1]) == z && parse(Float64, v[2]) == 3.0), L)
        rows = Vector{Float64}[]; i = i0 + 1
        while i <= length(L) && length(split(L[i])) == 30 && parse(Float64, split(L[i])[2]) != 3.0
            push!(rows, parse.(Float64, split(L[i]))[2:end]); i += 1
        end
        out[z] = Matrix(reduce(hcat, rows)')            # [iρ, iT]
    end
    out
end
# ------------------------------------------------------------------ nonuniform Hermite interpolation
@inline function slope(x, y, k)
    n = length(x)
    k == 1 && return (y[2] - y[1]) / (x[2] - x[1])
    k == n && return (y[n] - y[n-1]) / (x[n] - x[n-1])
    hl = x[k] - x[k-1]; hr = x[k+1] - x[k]
    (hr * (y[k] - y[k-1]) / hl + hl * (y[k+1] - y[k]) / hr) / (hl + hr)
end
function herm(x, y, xq)
    n = length(x); xq = clamp(xq, x[1], x[n])
    i = clamp(searchsortedlast(x, xq), 1, n - 1); h = x[i+1] - x[i]; t = (xq - x[i]) / h
    d0 = slope(x, y, i); d1 = slope(x, y, i + 1)
    t2 = t * t; t3 = t2 * t
    (2t3 - 3t2 + 1) * y[i] + (t3 - 2t2 + t) * h * d0 + (-2t3 + 3t2) * y[i+1] + (t3 - t2) * h * d1
end

