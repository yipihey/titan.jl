"""
    DampedNewtonS (S7c-1: DampedNewton with the pieces of TitanDriver.step! that the star driver needs: STAGE_HOOK, FULL_DT, FAIL_TFAC, Julia messages)

Alternative to TitanDriver.step!: TITAN's Newton iteration with (1) step-length clipping
of the multiplicative variables and (2) a backtracking line search on the row-equilibrated
residual norm, so that a hard step (wall reflection, shock formation) is first attacked by
damping instead of TITAN's "15 failed iterations → dt × 0.05 → restart".
The dt-cut retry stays as a last resort (factor `TFAC[]`, default TITAN's 0.05).

Damping is only engaged when the full Newton correction is large (dmax > `DTRIG[]`), so
well-behaved iterations are exactly TITAN's (same updates, same convergence test).
Requires the Julia residual (TitanResidFast) for the merit function and a matslv that does
not overwrite em2..ep2 (BlockSolver.Penta).
"""
module DampedNewtonS

using LinearAlgebra
import ..TitanDriver, ..TitanResidFast
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast
g(s) = TD.g(s); s!(s, x) = TD.s!(s, x); r(s) = TD.r(s)

const ALPHA = Ref(0.5)       # max |relative change| of r, m, ρ, T, E per iterate
const DTRIG = Ref(0.05)      # engage damping only when dmax of the full step exceeds this
const LMIN = Ref(1 / 32)     # smallest line-search step
const TFAC = Ref(0.05)       # last-resort dt factor (TITAN: 0.05)
const KFREEZE = Ref(0)       # >0: freeze branch decisions from Newton iteration KFREEZE on (0 = never)
const TRACE = Ref(false); const FAILS = Vector{Vector{Float64}}(); const HIST = Float64[]
const STATS = Dict{Symbol,Int}(:damped => 0, :backtracks => 0, :dtcuts => 0, :iters => 0)
resetstats!() = (for k in keys(STATS); STATS[k] = 0; end)

const PRIM = (:rn, :xmen, :un, :dn, :tn, :ern, :frn)
const SAVE = Dict{Symbol,Vector{Float64}}()
save!() = (for s in PRIM; copyto!(get!(() -> similar(TF.var(s)), SAVE, s), TF.var(s)); end)
restore!() = (for s in PRIM; copyto!(TF.var(s), SAVE[s]); end)

"row scales 1/max|J_row| from em2..ep2 (current Jacobian) -> S[i,k]"
function rowscales!(S, neqn, k0, nb)
    E = (TF.var(:em2), TF.var(:em1), TF.var(:e00), TF.var(:ep1), TF.var(:ep2))
    @inbounds for kb in 1:nb, i in 1:neqn
        k = k0 + kb - 1; m = 0.0
        for M in E, j in 1:neqn; m = max(m, abs(M[i, j, k])); end
        S[i, kb] = m > 0 ? inv(m) : 1.0
    end
    S
end
meritof(R, S) = sqrt(sum(abs2, R .* S))     # R, S: neqn × nb
function merit_now(S)
    c = TX.context()
    R = zeros(c.neqn, c.nb); TX.residual!(R, TX.getws(Float64, c.mgr), c, zeros(c.neqn * c.nb))
    meritof(R, S)
end

"Newton iteration for one step: returns true if converged."
function newton!()
    neqn = Int(g(:neqn)); ngrs = Int(g(:ngrs)); ngre = Int(g(:ngre)); k0 = ngrs - 2; nb = ngre - ngrs + 5
    S = zeros(neqn, nb); R0 = zeros(neqn, nb); δ = zeros(neqn, nb)
    rhs = TF.var(:rhs)
    mult = [Int(g(s)) for s in (:ir, :im, :id, :it)]
    Int(g(:lrad)) == 1 && push!(mult, Int(g(:ie)))
    s!(:iter, 1); empty!(HIST); dm1 = Inf; dm2 = Inf
    while true
        STATS[:iters] += 1
        kf = KFREEZE[]
        TX.FREEZE[] = kf == 0 || g(:iter) < kf ? 0 : (g(:iter) == kf ? 1 : 2)
        r(:matgen)
        TX.FREEZE[] == 1 && (TX.FREEZE[] = 2)        # merit evaluations replay the same branches
        rowscales!(S, neqn, k0, nb)
        @inbounds for kb in 1:nb, i in 1:neqn; R0[i, kb] = rhs[i, k0+kb-1]; end
        r(:matslv)
        r(:maxdel)
        dfull = g(:dmax); TRACE[] && push!(HIST, dfull)
        if dfull <= DTRIG[]                       # TITAN's undamped iteration
            r(:update)
            g(:dmax) <= g(:conv) && (TX.FREEZE[] = 0; return true)
            (g(:dmax) <= TD.CONV_STALL[] && g(:iter) >= 4 && abs(g(:dmax) - dm2) < 1e-3 * g(:dmax)) && (TX.FREEZE[] = 0; return true)
            dm2 = dm1; dm1 = g(:dmax)
        else
            STATS[:damped] += 1
            @inbounds for kb in 1:nb, i in 1:neqn; δ[i, kb] = rhs[i, k0+kb-1]; end
            dm = maximum(abs(δ[i, kb]) for i in mult, kb in 1:nb)
            λ = min(1.0, ALPHA[] / dm)
            f0 = meritof(R0, S)
            save!()
            while true
                @inbounds for kb in 1:nb, i in 1:neqn; rhs[i, k0+kb-1] = λ * δ[i, kb]; end
                r(:update)
                f = merit_now(S)
                (isfinite(f) && f <= (1 - 1e-4λ) * f0) && break
                λ <= LMIN[] && break                  # accept the smallest step anyway
                restore!(); λ /= 2; STATS[:backtracks] += 1
            end
            # convergence only on an undamped step (TITAN's criterion)
            λ == 1.0 && dfull <= g(:conv) && (TX.FREEZE[] = 0; return true)
        end
        s!(:iter, g(:iter) + 1)
        g(:iter) > g(:niter) && (TRACE[] && push!(FAILS, copy(HIST)); TX.FREEZE[] = 0; return false)
    end
end

"Drop-in for TitanDriver.step! (same control flow apart from the Newton iteration; STAGE_HOOK, FULL_DT, FAIL_TFAC and the Julia `msg` as in TitanDriver.step!)."
function step!()
    r(:peq)
    s!(:jtry, 1); s!(:jback, 0)
    while true
        TD.FULL_DT[] = g(:dtime)
        if newton!()
            s!(:lboos, 1)
            r(:matgen); r(:matslv); r(:maxdel); r(:update)
            s!(:lboos, 0)
            TD.STAGE_HOOK[] === nothing || TD.STAGE_HOOK[](:single, g(:dtime))
            TD.FIX_ENERGY[] && r(:totnrg)
            r(:timstp)
            xs = g(:smax) / g(:stol)
            if xs < 1.0e6
                TD.FIX_ENERGY[] || r(:totnrg)
                return nothing
            end
            TD.msg(3, xs); s!(:jback, g(:jback) + 1); g(:jback) > g(:nback) && TD.msg(4); TD.msg(5)
            r(:reset); TD.FIX_CLOCK[] && s!(:timen, g(:timeo) + g(:dtime))
        else
            STATS[:dtcuts] += 1
            TD.msg(1)
            s!(:tfac, TFAC[]); s!(:dtime, g(:tfac) * g(:dtime)); s!(:timen, g(:timeo) + g(:dtime))
            s!(:jtry, g(:jtry) + 1); g(:jtry) > g(:ntry) && TD.msg(2)
            r(:reset); TD.FIX_CLOCK[] && s!(:timen, g(:timeo) + g(:dtime))
        end
    end
end

end # module
