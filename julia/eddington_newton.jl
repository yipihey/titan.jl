# Phase 3: size and effect of the Jacobian terms TITAN drops by freezing the
# Eddington factors (f, g) during the Newton iteration.
# usage: julia eddington_newton.jl NAME NN STEP
#  * unknowns: TITAN's Newton variables on the non-phantom points
#    (ln r, ln m, ln ρ, ln T, ln E, u/u_nom, F/F_nom); phantom values are slaved by
#    update_'s boundary formulas.
#  * R_live(x): update_ (zero correction) + eddfac_ + matgen_  (f follows the iterate)
#    R_froz(x): update_ + matgen_                               (f frozen, TITAN ltran=1)
#  * Jacobians by uncolored central FD (no stencil assumption -> captures the
#    nonlocal f-response), J_froz ≈ TITAN's analytic matrix, J_live = full.
using Printf, LinearAlgebra
include(joinpath(@__DIR__, "TitanDriver.jl")); using .TitanDriver
const TD = TitanDriver; const TF = TD.TF

name, nn, nstep = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3])
TF.init(nn; rundir = joinpath(TF.ROOT, "lib", "p3_runs", "edd_$(name)"))
for i in 1:nstep-1; TD.titan_step!(); end
TD.r(:peq)                                     # predictor for step `nstep` = Newton start
ngrs = Int(TF.getscalar(:ngrs)); ngre = Int(TF.getscalar(:ngre)); conv = TF.getscalar(:conv)
ix(s) = Int(TF.getscalar(s))
VARS = [(ix(:ir), :rn, :log, ngrs:ngre+1, 1e-10), (ix(:im), :xmen, :log, ngrs:ngre+1, 1e-11),
        (ix(:id), :dn, :log, ngrs:ngre, 1e-6), (ix(:iu), :un, :lin, ngrs:ngre+1, 1e-6),
        (ix(:it), :tn, :log, ngrs:ngre, 1e-6), (ix(:ie), :ern, :log, ngrs:ngre, 1e-6),
        (ix(:if), :frn, :lin, ngrs:ngre+1, 1e-6)]
UNK = [(v, k) for v in VARS for k in v[4]]                 # unknown list
ROW = [(v[1], k) for v in VARS for k in v[4]]              # residual rows (eq i ↔ var i)
N = length(UNK)
base = TD.snapshot()
nom = Dict(:un => copy(TF.var(:unom)), :frn => copy(TF.var(:frnom)))
x0 = zeros(N)                                    # x = displacement in Newton variables
rhs = TF.var(:rhs)
function setx!(x)
    TD.restore!(base)
    for (n, ((j, arr, kind, kr, h), k)) in enumerate(UNK)
        a = TF.var(arr)
        a[k] = kind == :log ? a[k] * exp(x[n]) : a[k] + x[n] * nom[arr][k]
    end
end
function R(x; live)
    setx!(x); rhs .= 0.0
    TD.fcall(:update_); live && TD.fcall(:eddfac_); TD.fcall(:matgen_)
    [rhs[i, k] for (i, k) in ROW]
end
function J(x; live)
    M = zeros(N, N)
    for n in 1:N
        h = UNK[n][1][5]
        xp = copy(x); xp[n] += h; xm = copy(x); xm[n] -= h
        M[:, n] = (R(xp; live) .- R(xm; live)) ./ (2h)
    end
    # rows whose residual does not depend on x (prescribed-value BC rows: piston u,
    # reflecting-boundary F, outer wall u): TITAN uses rhs = 0, e00 = 1, i.e. Δ = 0
    for n in 1:N
        all(iszero, view(M, n, :)) && (M[n, n] = 1.0)
    end
    M
end
logvar = [UNK[n][1][3] == :log for n in 1:N]
function newton(; Rlive, Jlive, maxit = 15)
    x = copy(x0); hist = Float64[]
    for it in 1:maxit
        r = R(x; live = Rlive)
        Jm = J(x; live = Jlive)
        if !all(isfinite, Jm) || !all(isfinite, r)
            bad = findall(n -> !isfinite(r[n]) || any(!isfinite, view(Jm, n, :)), 1:N)
            println("    non-finite residual/Jacobian rows at iteration $it: ",
                    [(eqn[ROW[n][1]], ROW[n][2]) for n in bad[1:min(end, 8)]], " … total ", length(bad))
            return (-1, hist)
        end
        Δ = -(lu(Jm) \ r)
        nb = argmax(abs.(Δ))
        @printf("    it %d: |r|max=%.2e  max|Δ|=%.2e at %s k=%d\n", it, maximum(abs, r), abs(Δ[nb]), String(UNK[nb][1][2]), UNK[nb][2])
        x .+= Δ
        dmax = maximum(abs, Δ[logvar])
        push!(hist, dmax)
        dmax <= conv && return (it, hist)
    end
    (maxit + 1, hist)
end

# --- size of the dropped terms at the Newton start point
Jf = J(x0; live = false); Jl = J(x0; live = true)
D = Jl .- Jf
band(n, m) = abs(UNK[n][2] - UNK[m][2]) <= 2
@printf("%s, Newton start of step %d (t=%.5e, dt=%.4e), N = %d unknowns\n", name, nstep,
        TF.getscalar(:timen), TF.getscalar(:dtime), N)
eqn = Dict(ix(:ir) => "grid", ix(:im) => "mass", ix(:id) => "contin", ix(:iu) => "gasmom",
           ix(:it) => "gasnrg", ix(:ie) => "radnrg", ix(:if) => "radmom")
println("Eddington-factor response (J_live − J_frozen), per equation, relative to max|J_frozen| of the row:")
for i in sort(collect(keys(eqn)))
    rows = [n for n in 1:N if ROW[n][1] == i && count(!iszero, view(Jf, n, :)) > 1]
    inb = maximum(abs(D[n, m]) / maximum(abs, Jf[n, :]) for n in rows for m in 1:N if band(n, m); init = 0.0)
    outb = maximum(abs(D[n, m]) / maximum(abs, Jf[n, :]) for n in rows for m in 1:N if !band(n, m); init = 0.0)
    frac = sum(abs(D[n, m]) for n in rows for m in 1:N if !band(n, m); init = 0.0) /
           max(sum(abs(D[n, m]) for n in rows for m in 1:N; init = 0.0), floatmin())
    @printf("  %-7s in-band max %.2e   out-of-band max %.2e   out-of-band share of |ΔJ| %.2f\n", eqn[i], inb, outb, frac)
end
offband_froz = maximum(abs(Jf[n, m]) / maximum(abs, Jf[n, :]) for n in 1:N for m in 1:N if !band(n, m) && count(!iszero, view(Jf, n, :)) > 1)
@printf("  (check: frozen-f FD Jacobian outside the ±2-zone band: max %.2e of row max)\n", offband_froz)

# --- Newton iteration counts
for (lbl, Rl, Jl_) in (("TITAN ltran=1: f frozen in R and J", false, false),
                       ("f live in R, frozen-f Jacobian (ltran=0 as in TITAN)", true, false),
                       ("f live in R, full Jacobian incl. ∂f/∂x", true, true))
    it, hist = newton(; Rlive = Rl, Jlive = Jl_)
    println(@sprintf("  %-55s iterations to dmax ≤ %.0e: %s   dmax history: ", lbl, conv, it > 15 ? ">15" : string(it)),
            join([@sprintf("%.1e", h) for h in hist], " "))
end
TF.finalize()
