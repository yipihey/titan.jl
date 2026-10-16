# Effect of the dropped Eddington-factor Jacobian terms on Newton convergence,
# measured inside TITAN's own full block system (phantom rows included).
# usage: julia eddington_newton2.jl NAME NN STEP
#  ltran = 0  -> update_ recomputes f, g every iteration (TITAN's "live f" option)
#  (A) J = TITAN analytic (f frozen in J)          [what TITAN does for ltran = 0]
#  (B) J = analytic + ΔJ_f, ΔJ_f = FD[R(ltran=0)] − FD[R(ltran=1)] per non-phantom column
#  (C) ltran = 1 (f frozen in R and J)             [what the benchmarks do]
using Printf, LinearAlgebra
include(joinpath(@__DIR__, "TitanDriver.jl")); using .TitanDriver
const TD = TitanDriver; const TF = TD.TF

name, nn, nstep = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3])
TF.init(nn; rundir = joinpath(TF.ROOT, "lib", "p3_runs", "edd2_$(name)"))
for i in 1:nstep-1; TD.titan_step!(); end
ngrs = Int(TF.getscalar(:ngrs)); ngre = Int(TF.getscalar(:ngre)); neqn = Int(TF.getscalar(:neqn))
k0 = ngrs - 2; nb = ngre - ngrs + 5; N = neqn * nb
gidx(i, k) = (k - k0) * neqn + i
ix(s) = Int(TF.getscalar(s))
VARS = [(ix(:ir), :rn, :log, ngrs:ngre+1, 1e-10), (ix(:im), :xmen, :log, ngrs:ngre+1, 1e-11),
        (ix(:id), :dn, :log, ngrs:ngre, 1e-6), (ix(:iu), :un, :lin, ngrs:ngre+1, 1e-6),
        (ix(:it), :tn, :log, ngrs:ngre, 1e-6), (ix(:ie), :ern, :log, ngrs:ngre, 1e-6),
        (ix(:if), :frn, :lin, ngrs:ngre+1, 1e-6)]
rhs = TF.var(:rhs)
E = (TF.var(:em2), TF.var(:em1), TF.var(:e00), TF.var(:ep1), TF.var(:ep2))
function dense_system()
    A = zeros(N, N); b = zeros(N)
    for kb in 1:nb, o in -2:2
        cb = kb + o; (1 <= cb <= nb) || continue
        k = k0 + kb - 1
        for j in 1:neqn, i in 1:neqn
            A[(kb-1)*neqn+i, (cb-1)*neqn+j] = E[o+3][i, j, k]
        end
    end
    for kb in 1:nb, i in 1:neqn; b[(kb-1)*neqn+i] = rhs[i, k0+kb-1]; end
    A, b
end
function resid_full(ltran)
    TF.setscalar!(:ltran, Int32(ltran)); rhs .= 0.0
    TD.fcall(:update_); TD.fcall(:matgen_)
    [rhs[i, k] for k in k0:k0+nb-1 for i in 1:neqn]
end
"ΔJ_f at the current COMMON state (state restored afterwards)."
function dJ_edd()
    S = TD.snapshot(); D = zeros(N, N)
    nom = Dict(:un => copy(TF.var(:unom)), :frn => copy(TF.var(:frnom)))
    for (j, arr, kind, kr, h) in VARS, k in kr
        col = gidx(j, k); R = Dict{Tuple{Int,Int},Vector{Float64}}()
        for sgn in (1, -1), lt in (0, 1)
            TD.restore!(S); a = TF.var(arr)
            a[k] = kind == :log ? a[k] * exp(sgn * h) : a[k] + sgn * h * nom[arr][k]
            R[(sgn, lt)] = resid_full(lt)
        end
        D[:, col] = ((R[(1, 0)] .- R[(-1, 0)]) .- (R[(1, 1)] .- R[(-1, 1)])) ./ (2h)
    end
    TD.restore!(S); D
end
function newton(variant; maxit = 15)
    TD.restore!(S0)
    TF.setscalar!(:ltran, Int32(variant == :C ? 1 : 0))
    variant != :C && TD.fcall(:eddfac_)               # start from f of the predictor state
    hist = Float64[]
    for it in 1:maxit
        TF.setscalar!(:ltran, Int32(variant == :C ? 1 : 0))
        TD.fcall(:matgen_)
        A, b = dense_system()
        if variant == :B
            S = TD.snapshot(); D = dJ_edd(); TD.restore!(S)
            TF.setscalar!(:ltran, Int32(0)); A .+= D
        end
        sc = [let m = maximum(abs, view(A, r, :)); m > 0 ? 1 / m : 1.0 end for r in 1:N]   # row equilibration
        Δ = -(lu(sc .* A) \ (sc .* b))
        for kb in 1:nb, i in 1:neqn; rhs[i, k0+kb-1] = Δ[(kb-1)*neqn+i]; end
        TD.fcall(:maxdel_); TD.fcall(:update_)
        push!(hist, TF.getscalar(:dmax))
        TF.getscalar(:dmax) <= TF.getscalar(:conv) && return it, hist
    end
    maxit + 1, hist
end

TD.r(:peq)
S0 = TD.snapshot()
@printf("%s step %d (t=%.5e, dt=%.4e), N=%d (full system incl. phantom rows)\n", name, nstep,
        TF.getscalar(:timen), TF.getscalar(:dtime), N)
for (v, lbl) in ((:C, "C: ltran=1 (benchmarks): f frozen in R and J"),
                 (:A, "A: ltran=0, TITAN Jacobian (∂f/∂x dropped)"),
                 (:B, "B: ltran=0, Jacobian + ∂f/∂x (FD)"))
    it, h = newton(v)
    println(@sprintf("  %-48s iterations: %s  dmax: ", lbl, it > 15 ? ">15" : string(it)),
            join([@sprintf("%.1e", x) for x in h], " "))
end
TF.finalize()
