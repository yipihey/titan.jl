# Phase 3: quantify the Jacobian terms TITAN's hand-coded matgen omits.
# A colored central finite-difference Jacobian of the ACTUAL Fortran residual
# (perturb unknowns in COMMON -> update_ with zero correction -> [eddfac_] -> matgen_)
# is compared with TITAN's analytic blocks em2..ep2 on a mid-Newton state.
# usage: julia fd_jacobian.jl NAME NN STEP
using Printf, Serialization
include(joinpath(@__DIR__, "TitanDriver.jl")); using .TitanDriver
const TD = TitanDriver; const TF = TD.TF

name, nn, nstep = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3])
TF.init(nn; rundir = joinpath(TF.ROOT, "lib", "p3_runs", "fd_$(name)"))
for i in 1:nstep-1; TD.titan_step!(); end
TD.r(:peq); TD.r(:matgen); TD.r(:matslv); TD.r(:maxdel); TD.r(:update)   # mid-Newton state

ngrs = Int(TF.getscalar(:ngrs)); ngre = Int(TF.getscalar(:ngre)); neqn = Int(TF.getscalar(:neqn))
idx = Dict(s => Int(TF.getscalar(s)) for s in (:ir, :im, :id, :iu, :it, :ie, :if))
lrad = TF.getscalar(:lrad) == 1
# unknown j -> (array, kind, zone range)
VARS = [(idx[:ir], :rn, :log, ngrs:ngre+1), (idx[:im], :xmen, :log, ngrs:ngre+1),
        (idx[:id], :dn, :log, ngrs:ngre), (idx[:iu], :un, :lin, ngrs:ngre+1),
        (idx[:it], :tn, :log, ngrs:ngre)]
lrad && append!(VARS, [(idx[:ie], :ern, :log, ngrs:ngre), (idx[:if], :frn, :lin, ngrs:ngre+1)])
eqname = Dict(idx[:ir] => "grid", idx[:im] => "mass", idx[:id] => "contin", idx[:iu] => "gasmom",
              idx[:it] => "gasnrg")
lrad && (eqname[idx[:ie]] = "radnrg"; eqname[idx[:if]] = "radmom")
varname = Dict(v[1] => String(v[2]) for v in VARS)

rhs = TF.var(:rhs)
base = TD.snapshot()
function residual(; edd = false)
    rhs .= 0.0
    TD.fcall(:update_)
    edd && TD.fcall(:eddfac_)
    TD.fcall(:matgen_)
    R = copy(rhs[1:neqn, :])
    R
end
# analytic Jacobian at the base point (update with zero correction is idempotent)
R0 = residual()
Jan = Dict(o => copy(TF.var(s)) for (o, s) in ((-2, :em2), (-1, :em1), (0, :e00), (1, :ep1), (2, :ep2)))
TD.restore!(base)
nom = Dict(:un => copy(TF.var(:unom)), :frn => copy(TF.var(:frnom)))

"colored central FD; returns Dict((i,j,o) => Vector over row zone k)"
const HBASE = Dict(:rn => 1e-10, :xmen => 1e-11, :dn => 1e-6, :un => 1e-6, :tn => 1e-6, :ern => 1e-6, :frn => 1e-6)
function fdjac(hs; edd = false)
    J = Dict{Tuple{Int,Int,Int},Vector{Float64}}()
    for (j, arr, kind, kr) in VARS, c in 0:4
        h = HBASE[arr] * hs
        ks = [k for k in kr if mod(k, 5) == c]
        Rp = Rm = nothing
        for sgn in (+1, -1)
            TD.restore!(base)
            a = TF.var(arr)
            for k in ks
                a[k] = kind == :log ? a[k] * exp(sgn * h) : a[k] + sgn * h * nom[arr][k]
            end
            R = residual(; edd)
            sgn == 1 ? (Rp = R) : (Rm = R)
        end
        D = (Rp .- Rm) ./ (2h)
        for i in 1:neqn, kc in ks, o in -2:2
            k = kc - o                                  # row zone whose column offset o hits kc
            (ngrs <= k <= ngre + 1) || continue
            get!(() -> zeros(TF.MGR), J, (i, j, o))[k] = D[i, k]
        end
    end
    TD.restore!(base)
    J
end

J1 = fdjac(1.0); J2 = fdjac(4.0)
Je = lrad ? fdjac(1.0; edd = true) : nothing

rows = (ngrs + 2):(ngre - 2)           # rows whose ±2 stencil avoids phantom zones
rowmax = Dict((i, k) => maximum(abs(Jan[o][i, j, k]) for o in -2:2, j in 1:neqn) for i in 1:neqn, k in rows)
function compare(J; label)
    println("\n== $label  (rows k = $(first(rows))..$(last(rows)); |ΔJ| / max|J_an(row)|)")
    out = Dict{Tuple{Int,Int},Float64}()
    for i in 1:neqn, (j, _...) in VARS, o in -2:2, k in rows
        a = Jan[o][i, j, k]; f = haskey(J, (i, j, o)) ? J[(i, j, o)][k] : 0.0
        d = abs(f - a) / max(rowmax[(i, k)], floatmin())
        out[(i, j)] = max(get(out, (i, j), 0.0), d)
    end
    out
end
c1 = compare(J1; label = "FD(h) vs analytic")
c2 = compare(J2; label = "FD(4h) vs analytic")
ce = Je === nothing ? nothing : compare(Je; label = "FD incl. Eddington-factor response vs analytic")
@printf("%s step %d: neqn=%d, t=%.5e, dt=%.4e\n", name, nstep, neqn, TF.getscalar(:timen), TF.getscalar(:dtime))
println("equation \\ variable : max |ΔJ|/rowmax   [FD h | FD 4h | FD h incl. Eddington response]; h: r 1e-10, m 1e-11, others 1e-6")
for i in sort(collect(keys(eqname)))
    line = @sprintf("  %-7s", eqname[i])
    for (j, arr, _...) in VARS
        line *= @sprintf("  %s: %.1e|%.1e%s", varname[j], c1[(i, j)], c2[(i, j)],
                         ce === nothing ? "" : @sprintf("|%.1e", ce[(i, j)]))
    end
    println(line)
end
serialize(joinpath(TF.ROOT, "lib", "p3_runs", "fdjac_$(name)_$nstep.jls"), (; c1, c2, ce, eqname, varname))
TF.finalize()
