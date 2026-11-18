"""
    NewtonKrylov

NKRYLOV (popiii_m2c `NKRYLOV=1[,η[,mmax[,ε[,central]]]]`, provisional, default off = bit-identical; docs/POPIII_NEWTON_KRYLOV.md).

Newton's linear system J δ = −R(x) is solved on the TRUE Jacobian J of the residual after a full `update!`, not on the assembled band Jacobian J_b
(which lacks ∂/∂x of everything `update!` refreshes outside `residual!`: the Eddington factor f, h, the flux-tie weight, ...). The band Jacobian is the
preconditioner: GMRES on A δ = δ0, A v = P J_b⁻¹ (J v), δ0 = P(−J_b⁻¹ R) the correction the code's own matslv has just written into `rhs`.

* J v = [R(x ⊕ εv) − R(x ⊖ εv)]/2ε (central; one-sided [R(x ⊕ εv) − R(x)]/ε with the 5th field of NKRYLOV = 0). x ⊕ εv: εv is written into `rhs` and `update!` is called
  (the Newton update: EOS, opacity, Eddington factors, tie weights, species), then a fresh `context()` and the Float64 `residual!` at zero correction. Afterwards every COMMON
  array and the closure caches are restored from a snapshot (exact copies) and a residual evaluation at the restored state is compared with R(x) (the drift counter).
* ε = EPS/max|Wv| with the GMRES weights W (below), EPS = 1e-5: the largest weighted entry of εv is 1e-5. The residual carries round-off noise of the formal solution
  (1 − (1 + τ)e^{−τ} in thin zones, up to 1e-6 relative in f and h); a one-sided difference with max|εv| = 1e-7 divides that by 1e-7 and GMRES stagnates, the central difference
  has O(ε²) truncation and tolerates 1e-5; measuring ε in the weighted norm keeps the entries the norm ignores (F in the opaque interior, relative changes of 1e2) from forcing ε down.
* The branch decisions of the residual (`decide`: limiter, sign, upwind, viscosity) are recorded at R(x) and replayed (TitanResidFast.FREEZE = 1, 2) in the
  perturbed evaluations, as the assembled Jacobian differentiates the base state's branches.
* Only the entries `update!` applies are unknowns (P projects onto them): interface unknowns r, m, u, F of zones ngrs..ngre+1, zone unknowns ρ, T, E_r of ngrs..ngre,
  species in coupled zones. The ghost entries are slaved to them by the boundary conditions inside `update!`.
* GMRES without restart, modified Gram-Schmidt, Givens rotations; stops when the preconditioned residual is ≤ η ‖Pδ0‖ in the weighted 2-norm ‖W·‖₂, W the weights of the DTERR error norm (`weights`).
  If it does not get there in `mmax` iterations the last iterate (the best) is used and the event counted. With an exact band Jacobian Pδ0 is an eigenvector with eigenvalue 1 and the solve ends after one product.
Then the unchanged `maxdel!`, `update!`.
"""
module NewtonKrylov

using LinearAlgebra, Printf, Statistics
import ..TitanDriver, ..TitanCore, ..TitanResidFast, ..Trbdf2
const TD = TitanDriver; const TF = TD.TF; const TC = TitanCore; const TX = TitanResidFast
vf(s) = TF.var(s)::Vector{Float64}
gi(s) = Int(TD.g(s))
gf(s) = Float64(TD.g(s))

mutable struct Opts
    eta::Float64; mmax::Int; eps::Float64
    central::Bool                    # central (default) or one-sided differences for J v
    check::Int                       # NKCHECK=<Newton iteration>: matrix-vector check at that iteration of the first attempt, then exit (diagnostic)
end
const OPT = Opts(1e-3, 30, 1e-5, true, 0)
mutable struct Count
    solves::Int; mv::Int; nonconv::Int; fallback::Int; ndrift::Int; maxdrift::Float64; secs::Float64
end
const CNT = Count(0, 0, 0, 0, 0, 0.0, 0.0)
"NKLS=1[,trigger[,λmin]]: backtracking line search on T(λ) = ‖P M⁻¹ R(x ⊕ λ x_gmres)‖_W (true residual, branches re-evaluated); default off = bit-identical."
mutable struct LSOpts
    on::Bool; trigger::Float64; lmin::Float64
end
const LS = LSOpts(false, 1e-2, 1 / 64)
mutable struct LSCount
    engaged::Int; backtracks::Int; fail::Int; evals::Int; maxrel0::Float64
end
const LSC = LSCount(0, 0, 0, 0, 0.0)
const MARK = Ref((0, 0, 0, 0, 0))        # counters and step count at the last trace line

# ---------------------------------------------------------------- snapshot of everything update! writes
const SKIP = (:em2, :em1, :e00, :ep1, :ep2, :rhs)       # the Jacobian and the correction are not state
mutable struct Snap
    arrs::Dict{Symbol,Array}; extra::Vector{Vector{Float64}}; gx::Float64
end
Snap() = Snap(Dict{Symbol,Array}(), [Float64[] for _ in 1:6], 0.0)
extras() = (TD.HEDD, TD.HXEDD, TD.FXEDD, TD.TIEW, TX.XSN[][1], TX.XSN[][2])      # closure h, h_x, f_x, tie weight; species iterate
function take!(S::Snap)
    for (n, a) in TF.JSTATE
        n in SKIP && continue
        haskey(S.arrs, n) ? copyto!(S.arrs[n], a) : (S.arrs[n] = copy(a))
    end
    for (i, a) in enumerate(extras()); S.extra[i] = copy(a); end
    S.gx = TD.GXEDD[]
    nothing
end
function restore!(S::Snap)
    for (n, a) in S.arrs; copyto!(TF.JSTATE[n], a); end
    for (b, a) in zip(extras(), S.extra)
        length(b) == length(a) || resize!(b, length(a))
        copyto!(b, a)
    end
    TD.GXEDD[] = S.gx
    nothing
end

# ---------------------------------------------------------------- the unknowns, residual, preconditioner
const LOGCODES = (0, 3, 5, 6)       # equation codes of the unknowns applied as x e^δ: species, ρ, T, E_r
"Lower bound of the trust-region weights of the unknowns that are not logarithmic (u, F_r, m, r): an entry the GMRES norm all but ignores (weight ~1e-6) is still limited."
const TRFLOOR = 1e-2
"A finite-difference step changes a logarithmic unknown by at most this (δ ln x): an entry of weight ~1e-12 (a trace species) otherwise gets ε v ~ 1e5, x → 0 and NaN in the EOS."
const FDCLIP = 0.1
struct Sys
    ne::Int; nb::Int; k0::Int; n::Int
    code::Vector{Int}; act::BitVector            # act[(kb-1)ne+i]: the entry update! applies
    wt::Vector{Float64}                         # weight of each entry in the GMRES norm (Trbdf2.wnorm's)
    wtr::Vector{Float64}                        # weight in the trust-region measure: 1 for the log variables (ρ, T, E_r, species), max(wt, TRFLOOR) for the others
    R0::Vector{Float64}                         # residual at the base state (Float64 path)
    snap::Snap
    AB::Matrix{Float64}; ipiv::Vector{LinearAlgebra.BlasInt}; sc::Vector{Float64}; kl::Int; ku::Int
end

function resid(ne, nb)
    c = TX.context(); ws = TX.getws(Float64, c.mgr); R = zeros(ne, nb)
    TX.residual!(R, ws, c, zeros(ne * nb)); vec(R)
end

function activeset(c)
    ne, nb, k0 = c.neqn, c.nb, c.k0; a, b = c.ngrs, c.ngre
    act = falses(ne * nb)
    for kb in 1:nb, i in 1:ne
        k = k0 + kb - 1; cd = c.code[i]
        act[(kb-1)*ne+i] = cd in (1, 2, 4, 7) ? a <= k <= b + 1 :
                           cd in (3, 5, 6) ? a <= k <= b :
                           (cd == 0 && i in (8, 9) && TC.CPLH[] !== nothing) ? (a <= k <= b && TX.cplz(k)) : false
    end
    TX.MG[] && error("NKRYLOV: MGROUPS is not supported")
    act
end

"""
Weights of the correction entries in the GMRES norm: those of the DTERR error norm `Trbdf2.wnorm` (ρ, T relative; E_r against E + efl aT⁴; species against x + xfl;
u against max(|u|, c_s); m against the zone mass; r against the smaller zone width on a moving mesh; F against max(|F|, ffl c E)), so that corrections of a flux far below
c E (relative changes of 1e2 in the opaque interior) are not measured in units of that flux. The linear tolerance is a fraction of the weighted 2-norm of δ0.
"""
function weights(c, act)
    ne, nb, k0 = c.neqn, c.nb, c.k0; a, b = c.ngrs, c.ngre; o = Trbdf2.EOPT
    lag = gi(:lgrid) == 3 || TC.GRIDMOVE[]
    rn, xmen, dn, un, tn, unom, asn, dvn = (vf(s) for s in (:rn, :xmen, :dn, :un, :tn, :unom, :asn, :dvoln))
    ern, frn, frnom = vf(:ern), vf(:frn), vf(:frnom); wt = ones(ne * nb)
    for kb in 1:nb, i in 1:ne
        r = (kb - 1) * ne + i; act[r] || continue; k = k0 + kb - 1; cd = c.code[i]
        wt[r] = cd == 3 || cd == 5 ? 1.0 :
                cd == 6 ? ern[k] / (ern[k] + o.efl * 7.565723e-15 * tn[k]^4) :
                cd == 0 ? (q = TX.XSN[][i-7][k]; q / (q + o.xfl)) :
                cd == 4 ? unom[k] / max(abs(un[k]), max(asn[max(k-1, a)], asn[min(k, b)])) :
                cd == 2 ? xmen[k] / max(abs(xmen[k]), dn[max(k-1, a)] * dvn[max(k-1, a)], dn[min(k, b)] * dvn[min(k, b)]) :
                cd == 1 ? (lag ? rn[k] / min(rn[k] - rn[max(k-1, a)], rn[min(k+1, b+1)] - rn[k]) : 1.0) :
                frnom[k] / max(abs(frn[k]), o.ffl * 2.99792458e10 * 0.5 * (ern[max(k-1, a)] + ern[min(k, b)]))
    end
    [isfinite(x) ? max(abs(x), 1e-12) : 1.0 for x in wt]       # (a 0/0 at a boundary interface: weight 1)
end

"Pivoted banded LU of the (row-equilibrated) band Jacobian in COMMON (em2..ep2), as JuliaSolver.matslv!: the preconditioner."
function bandlu(ne, nb, k0)
    E = (vf2(:em2), vf2(:em1), vf2(:e00), vf2(:ep1), vf2(:ep2)); n = ne * nb; kl = ku = 3ne - 1
    AB = zeros(2kl + ku + 1, n); sc = zeros(n)
    @inbounds for kb in 1:nb, o in -2:2
        cb = kb + o; (1 <= cb <= nb) || continue
        M = E[o+3]; k = k0 + kb - 1
        for j in 1:ne, i in 1:ne
            r = (kb - 1) * ne + i; sc[r] = max(sc[r], abs(M[i, j, k]))
        end
    end
    @inbounds for r in 1:n; sc[r] = sc[r] > 0 ? inv(sc[r]) : 1.0; end
    @inbounds for kb in 1:nb, o in -2:2
        cb = kb + o; (1 <= cb <= nb) || continue
        M = E[o+3]; k = k0 + kb - 1
        for j in 1:ne, i in 1:ne
            r = (kb - 1) * ne + i; cc = (cb - 1) * ne + j
            AB[kl+ku+1+r-cc, cc] = M[i, j, k] * sc[r]
        end
    end
    AB, ipiv = LAPACK.gbtrf!(kl, ku, n, AB)
    AB, ipiv, sc, kl, ku
end
vf2(s) = TF.var(s)::Array{Float64,3}

"w = J_b⁻¹ y"
function precond!(w, S::Sys, y)
    @inbounds for r in 1:S.n; w[r] = y[r] * S.sc[r]; end
    LAPACK.gbtrs!('N', S.kl, S.ku, S.n, S.AB, S.ipiv, reshape(w, S.n, 1))
    w
end

"J_b w from the band matrices in COMMON (diagnostic)."
function bandmul(S::Sys, w)
    E = (vf2(:em2), vf2(:em1), vf2(:e00), vf2(:ep1), vf2(:ep2)); ne, nb, k0 = S.ne, S.nb, S.k0
    z = zeros(S.n); sa = zeros(S.n)
    for kb in 1:nb, o in -2:2
        cb = kb + o; (1 <= cb <= nb) || continue
        M = E[o+3]; k = k0 + kb - 1
        for j in 1:ne, i in 1:ne
            t = M[i, j, k] * w[(cb-1)*ne+j]; z[(kb-1)*ne+i] += t; sa[(kb-1)*ne+i] += abs(t)
        end
    end
    z, sa
end

"""
Residual at x ⊕ s·v: s·v is written into `rhs`, `update!` applied (branch decisions replayed), a fresh context and the Float64 residual evaluated; the state is restored.
`realized` (optional) receives the change the update actually made to every unknown including the ghost entries, per unit s.
"""
function perturbed(S::Sys, v, s; realized = nothing, freeze = 2, inspect = nothing, clip = Inf)
    ne, nb, k0 = S.ne, S.nb, S.k0; rhs = TF.var(:rhs)
    base = realized === nothing ? nothing : (c = TX.context(); (copy(c.base), copy(vf(:unom)), copy(vf(:frnom))))
    @inbounds for kb in 1:nb, i in 1:ne; rhs[i, k0+kb-1] = S.act[(kb-1)*ne+i] ? (S.code[i] in LOGCODES ? clamp(s * v[(kb-1)*ne+i], -clip, clip) : s * v[(kb-1)*ne+i]) : 0.0; end
    TX.FREEZE[] = freeze; updated = false
    try
        TD.r(:update); updated = true
        inspect === nothing || inspect()          # (NKSCAN: the state after the trial update, before the restore)
        R1 = resid(ne, nb)
        if realized !== nothing
            cnew = TX.context()
            for kb in 1:nb, i in 1:ne
                cd = S.code[i]; k = k0 + kb - 1; r = (kb - 1) * ne + i
                cd == 0 && (realized[r] = S.act[r] ? v[r] : 0.0; continue)       # species: applied as given
                nw = cnew.base[k, cd]; b0 = base[1][k, cd]
                realized[r] = nw == b0 ? 0.0 : cd == 4 ? (nw - b0) / base[2][k] / s : cd == 7 ? (nw - b0) / base[3][k] / s : log(nw / b0) / s
            end
        end
        R1
    catch
        (inspect === nothing || updated) || inspect()
        rethrow()
    finally
        restore!(S.snap); TX.FREEZE[] = 0
    end
end
"""
Backtracking line search along the GMRES step x. T(λ) = ‖W · P M⁻¹ R(x_k ⊕ λx)‖₂ with the TRUE residual at the trial state (`freeze = 0`: branch decisions re-evaluated, not replayed).
T(0) is computed from S.R0. Accept the first λ = 1, 1/2, … ≥ λmin with T ≤ (1 − λ/4) T(0); otherwise the finite trial of smallest T if < T(0), else λmin (counted as a failure).
Returns (λ, T0, T(λ), trials, ‖sc R‖ at 0, ‖sc R‖ at λ, relative difference of T(0) from ‖W δ0‖, failed).
"""
function linesearch(S::Sys, x, b)
    w = zeros(S.n)
    function level(R)
        precond!(w, S, R); w[.!S.act] .= 0.0; w .*= S.wt; norm(w)
    end
    scn(R) = norm(S.sc .* R)
    T0 = level(S.R0); sc0 = scn(S.R0); nb0 = norm(b .* S.wt)
    rel0 = abs(T0 - nb0) / max(nb0, 1e-300); LSC.maxrel0 = max(LSC.maxrel0, rel0)
    λ = 1.0; ntr = 0; best = (Inf, 0.0, NaN); acc = nothing; Tlast = NaN; sclast = NaN
    while λ >= LS.lmin * (1 - 1e-12)
        ntr += 1; LSC.evals += 1
        T, scl = try
            R = perturbed(S, x, λ; freeze = 0)
            level(R), scn(R)
        catch ex
            ex isa InterruptException && rethrow()
            Inf, NaN
        end
        Tlast = T; sclast = scl
        isfinite(T) && T < best[1] && (best = (T, λ, scl))
        if isfinite(T) && T <= (1 - λ / 4) * T0
            acc = (λ, T, scl); break
        end
        λ /= 2
    end
    acc === nothing || (acc[1] < 1 && (LSC.backtracks += 1); return acc[1], T0, acc[2], ntr, sc0, acc[3], rel0, false)
    LSC.fail += 1; LSC.backtracks += 1
    if best[1] < T0
        return best[2], T0, best[1], ntr, sc0, best[3], rel0, true
    end
    LS.lmin, T0, Tlast, ntr, sc0, sclast, rel0, true          # (no finite trial below T(0): lmin; the printed T is that of the last (smallest-lambda) trial)
end
"""
Directional derivative of the full residual, y = J v: central difference [R(x ⊕ εv) − R(x ⊖ εv)]/2ε (one-sided [R(x ⊕ εv) − R(x)]/ε if `OPT.central` is false). v is a vector over the active entries.
"""
function fdjv!(y, S::Sys, v, ε; realized = nothing, freeze = 2)
    R1 = perturbed(S, v, ε; realized, freeze, clip = FDCLIP)
    if OPT.central
        R2 = perturbed(S, v, -ε; freeze, clip = FDCLIP)
        @inbounds for r in 1:S.n; y[r] = (R1[r] - R2[r]) / (2ε); end
    else
        @inbounds for r in 1:S.n; y[r] = (R1[r] - S.R0[r]) / ε; end
    end
    y
end

function setup()
    TF.BACKEND[] === :julia || error("NKRYLOV needs the Julia state backend")
    c = TX.context(); ne, nb, k0 = c.neqn, c.nb, c.k0
    act = activeset(c); snap = Snap(); take!(snap)
    TX.FREEZE[] = 1
    R0 = try resid(ne, nb) finally TX.FREEZE[] = 0 end
    AB, ipiv, sc, kl, ku = bandlu(ne, nb, k0)
    wt = weights(c, act)
    wtr = [c.code[(r-1)%ne+1] in LOGCODES ? 1.0 : max(wt[r], TRFLOOR) for r in 1:ne*nb]
    Sys(ne, nb, k0, ne * nb, copy(c.code), act, wt, wtr, R0, snap, AB, ipiv, sc, kl, ku)
end

# ---------------------------------------------------------------- GMRES
"Solve A x = b (A v = P J_b⁻¹ J v) from x = 0 by GMRES in the norm ‖W·‖₂ (weights S.wt); returns (x, iterations, converged, relative residual)."
function gmres(S::Sys, b)
    n = S.n; mmax = OPT.mmax; η = OPT.eta; wt = S.wt; bs = b .* wt           # scaled space: ṽ = W v
    β = norm(bs); x = zeros(n)
    β == 0 && return x, 0, true, 0.0
    V = zeros(n, mmax + 1); H = zeros(mmax + 1, mmax); cs = zeros(mmax); sn = zeros(mmax); gv = zeros(mmax + 1)
    V[:, 1] = bs ./ β; gv[1] = β
    y = zeros(n); w = zeros(n); v = zeros(n); m = 0; conv = false; res = β
    for j in 1:mmax
        v .= (@view V[:, j]) ./ wt; v[.!S.act] .= 0.0
        fdjv!(y, S, v, OPT.eps / maximum(abs, @view V[:, j])); CNT.mv += 1       # ε v: the largest weighted entry changes by OPT.eps
        all(isfinite, y) || (rb = findfirst(!isfinite, y); error(@sprintf("NKRYLOV: non-finite J v (%d entries; first row %s; ε %.3g)", count(!isfinite, y), zname(S, rb), OPT.eps / maximum(abs, @view V[:, j]))))
        precond!(w, S, y); w[.!S.act] .= 0.0; w .*= wt
        for i in 1:j                                       # modified Gram-Schmidt
            h = dot(w, @view V[:, i]); H[i, j] = h; w .-= h .* @view V[:, i]
        end
        H[j+1, j] = norm(w)
        for i in 1:j-1                                     # earlier rotations
            t = cs[i] * H[i, j] + sn[i] * H[i+1, j]; H[i+1, j] = -sn[i] * H[i, j] + cs[i] * H[i+1, j]; H[i, j] = t
        end
        ρ = hypot(H[j, j], H[j+1, j]); cs[j] = ρ == 0 ? 1.0 : H[j, j] / ρ; sn[j] = ρ == 0 ? 0.0 : H[j+1, j] / ρ
        H[j, j] = ρ; hn = H[j+1, j]; H[j+1, j] = 0.0
        gv[j+1] = -sn[j] * gv[j]; gv[j] = cs[j] * gv[j]
        m = j; res = abs(gv[j+1])
        if res <= η * β || hn == 0.0 || !isfinite(res)
            conv = res <= η * β || hn == 0.0; break
        end
        V[:, j+1] = w ./ hn
    end
    yk = UpperTriangular(H[1:m, 1:m]) \ gv[1:m]
    for j in 1:m; x .+= yk[j] .* @view V[:, j]; end
    x ./= wt; x[.!S.act] .= 0.0
    x, m, conv, res / β
end

function wtop(b, S::Sys, ne, k0)
    r = argmax(abs.(b .* S.wt)); kb = (r - 1) ÷ ne + 1; i = (r - 1) % ne + 1
    @sprintf("%d@%d:%.3g(w %.2g)", i, k0 + kb - 1 - gi(:ngrs) + 1, b[r], S.wt[r])
end

# ---------------------------------------------------------------- NKSCAN (diagnostic, default off): what the true residual does along the GMRES step
"NKSCAN=<att>,<it>[,<it>...]: at Newton iteration(s) <it> of attempt <att> (needs NWTRACE for the attempt counter) print the line scan, column probes and state; no effect on the solve."
const SCAN = Ref((0, Int[]))
const NMS = ("r", "m", "rho", "u", "T", "E_r", "F_r", "x_H+", "x_H2")
zname(S::Sys, r) = @sprintf("%d@%d", (r - 1) % S.ne + 1, (r - 1) ÷ S.ne + 1 + S.k0 - 1 - gi(:ngrs) + 1)
toplist(S::Sys, v, n) = join((@sprintf("%s:%.4g", zname(S, r), v[r]) for r in sortperm(abs.(v); rev = true)[1:min(n, length(v))]), " ")
function levelT(S::Sys, R)
    w = zeros(S.n); precond!(w, S, R); w[.!S.act] .= 0.0; w .*= S.wt
    norm(w), w
end
"Unknowns of the current (trial) state that are non-finite or non-positive where they must be positive."
function badstate(S::Sys)
    a, b = gi(:ngrs), gi(:ngre); out = String[]
    rn, xmen, un, frn, dn, tn, ern = (vf(s) for s in (:rn, :xmen, :un, :frn, :dn, :tn, :ern)); x1, x2 = TX.XSN[][1], TX.XSN[][2]
    for k in a:b+1, (nm, v) in (("r", rn), ("m", xmen), ("u", un), ("F_r", frn))       # (r, m: only finiteness; m is an offset that is negative in the base state)
        q = v[k]; isfinite(q) || push!(out, @sprintf("%s@%d=%.5g", nm, k - a + 1, q))
    end
    for k in a:b, (nm, v) in (("rho", dn), ("T", tn), ("E_r", ern), ("x_H+", x1), ("x_H2", x2))
        q = v[k]; (!isfinite(q) || q <= 0) && push!(out, @sprintf("%s@%d=%.5g", nm, k - a + 1, q))
    end
    out
end
function scantrial(S::Sys, x, λ, fr)
    bad = String[]; exc = nothing; R = nothing
    try R = perturbed(S, x, λ; freeze = fr, inspect = () -> (bad = badstate(S)))
    catch ex
        ex isa InterruptException && rethrow()
        exc = ex
    end
    R, exc, bad
end
"Number of branch decisions of the trial state (recorded with FREEZE = 1, the values are those of FREEZE = 0) that differ from the base record; the base record is restored."
function decdiff(S::Sys, x, λ)
    base = Dict(k => copy(v) for (k, v) in TX.DEC); a = gi(:ngrs)
    try perturbed(S, x, λ; freeze = 1)
    catch ex; ex isa InterruptException && rethrow(); end
    nd = 0; parts = String[]
    for (k, v) in TX.DEC
        haskey(base, k) || (push!(parts, "newkey $(k[1])/$(k[2])"); continue)
        d = findall(v .!= base[k]); nd += length(d)
        isempty(d) || push!(parts, @sprintf("%s/%d:%d[k %s]", k[1], k[2], length(d), join(first(d, 4) .- a .+ 1, ",")))
    end
    for k in collect(keys(TX.DEC)); haskey(base, k) ? copyto!(TX.DEC[k], base[k]) : delete!(TX.DEC, k); end
    nd, join(parts, " ")
end
function scan!(S::Sys, b, x)
    ne, nb, k0 = S.ne, S.nb, S.k0; a = gi(:ngrs); wb = b .* S.wt; wx = x .* S.wt; pre = "NKSCAN "
    @printf("%sA att %d it %d jstep %d dtime %.6g theta %g |b|_W %.6g |x|_W %.6g\n", pre, Trbdf2.NWATT[], gi(:iter), gi(:jstep), TD.g(:dtime), TD.g(:thet), norm(wb), norm(wx))
    println(pre, "A largest weighted x (eq@zone:value): ", toplist(S, wx, 8)); println(pre, "A largest weighted b: ", toplist(S, wb, 8))
    rn, dn, tn, ern, frn, un = (vf(s) for s in (:rn, :dn, :tn, :ern, :frn, :un)); x1, x2 = TX.XSN[][1], TX.XSN[][2]
    println(pre, "B state (zone: r rho T E_r F_r x_H+ x_H2 u)")
    for z in 26:33
        k = a + z - 1; k > gi(:ngre) + 1 && continue
        @printf("%sB zone %d: r %.6e rho %.6e T %.6e E_r %.6e F_r %.6e x_H+ %.6e x_H2 %.6e u %.6e\n", pre, z, rn[k], k <= gi(:ngre) ? dn[k] : NaN, k <= gi(:ngre) ? tn[k] : NaN, k <= gi(:ngre) ? ern[k] : NaN, frn[k], k <= gi(:ngre) ? x1[k] : NaN, k <= gi(:ngre) ? x2[k] : NaN, un[k])
    end
    for z in (29, 30)
        kb = a + z - 1 - k0 + 1
        println(pre, "B weights zone ", z, ": ", join((@sprintf("%s %.4g", NMS[i], S.wt[(kb-1)*ne+i]) for i in 1:min(ne, 9)), " "))
    end
    println(pre, "B raw x (r m rho u T E_r F_r x_H+ x_H2)")
    for z in 27:32
        kb = a + z - 1 - k0 + 1
        println(pre, "B raw x zone ", z, ": ", join((@sprintf("%.4e", x[(kb-1)*ne+i]) for i in 1:min(ne, 9)), " "))
    end
    # C: line scan
    T0, _ = levelT(S, S.R0); @printf("%sC T0 %.6g (norm(b W) %.6g), |sc R0| %.6g\n", pre, T0, norm(wb), norm(S.sc .* S.R0))
    lams = [1e-6, 1e-5, 1e-4, 1e-3, 1/256, 1/128, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1.0]
    for sg in (-1, 1), l in lams
        λ = sg * l; line = ""
        for fr in (0, 2)
            R, exc, bad = scantrial(S, x, λ, fr)
            if exc !== nothing
                println(pre, "C lambda ", λ, " freeze ", fr, " EXCEPTION ", typeof(exc), ": ", first(split(sprint(showerror, exc), '\n')))
            else
                T, w = levelT(S, R)
                if isfinite(T)
                    @printf("%sC lambda %+.6g freeze %d T %.6g T/T0 %.6g linear |1-l| %.6g |sc R| %.6g worst row %s\n", pre, λ, fr, T, T / T0, abs(1 - λ), norm(S.sc .* R), zname(S, argmax(abs.(w))))
                else
                    nf = findall(!isfinite, R)
                    @printf("%sC lambda %+.6g freeze %d T non-finite (%d non-finite residual rows): %s\n", pre, λ, fr, length(nf), join((zname(S, r) * "=" * string(R[r]) for r in first(nf, 10)), " "))
                end
            end
            isempty(bad) || println(pre, "C   state after update (lambda ", λ, " freeze ", fr, "): ", length(bad), " flagged: ", join(first(bad, 12), " "))
        end
        nd, parts = decdiff(S, x, λ)
        println(pre, "C   lambda ", λ, " branch decisions differing from the base record: ", nd, " ", parts)
    end
    # D: column probes
    r1 = argmax(abs.(wx)); e = zeros(S.n); e[r1] = 1.0; xn = norm(wx); v = x ./ xn
    y = zeros(S.n); out = zeros(S.n)
    function Aop(vv, eps, fr)
        fdjv!(y, S, vv, eps / maximum(abs, vv .* S.wt); freeze = fr)
        precond!(out, S, y); out[.!S.act] .= 0.0; copy(out)
    end
    println(pre, "D unit entry e = ", zname(S, r1), " (weight ", S.wt[r1], "); v = x/|x|_W")
    for (fr, epss) in ((2, (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7)), (0, (1e-5,)))
        for eps in epss
            for (nmv, vv) in (("e", e), ("v", v))
                try
                    o = Aop(vv, eps, fr); d = nmv == "e" ? o .- e : o
                    @printf("%sD freeze %d eps %.0e %s: |Av|_W %.6g |Av - %s|_W %.6g  top W*(Av%s): %s\n", pre, fr, eps, nmv, norm(o .* S.wt), nmv == "e" ? "e" : "0", norm(d .* S.wt), nmv == "e" ? "-e" : "", toplist(S, d .* S.wt, 10))
                catch ex
                    ex isa InterruptException && rethrow()
                    println(pre, "D freeze ", fr, " eps ", eps, " ", nmv, " EXCEPTION ", typeof(ex), ": ", first(split(sprint(showerror, ex), '\n')))
                end
            end
        end
    end
    # E: raw column against the band column
    fdjv!(y, S, e, 1e-5 / maximum(abs, e .* S.wt)); z, _ = bandmul(S, e)
    d = abs.(y .- z) .* S.sc
    println(pre, "E rows with the largest row-scaled |J e - J_b e| (e = ", zname(S, r1), "; values times S.sc): row, J e, J_b e, diff")
    for r in sortperm(d; rev = true)[1:12]
        @printf("%sE %s: Je %+.5e Jb e %+.5e diff %+.5e (sc %.3e)\n", pre, zname(S, r), y[r] * S.sc[r], z[r] * S.sc[r], (y[r] - z[r]) * S.sc[r], S.sc[r])
    end
    nothing
end

# ---------------------------------------------------------------- globalisation: trust region + Armijo line search on the weighted true residual
"""
Newton globalisation (NKGLOB, default on with NKRYLOV; docs/POPIII_SWITCH_AUDIT2.md). One method for every unknown:

* Positivity by construction: `update!` applies ρ, T, E_r (and the species) as x ← x e^δ with δ the relative Newton correction (the step in ln x), r and m as x(1 + δ);
  a trial state that is non-finite, has x ≤ 0, or has interfaces that cross (r not increasing, m not decreasing) has merit Inf (`trial_ok`).
* Merit function T(λ) = ‖W J_b⁻¹ R(x ⊕ λ s)‖₂ (J_b the band Jacobian of the iteration's base state, W the GMRES weights): the residual of all rows with the base state's branch
  decisions (limiter, sign, upwind, viscosity) replayed, as in J v and J_b: the model the step solves. The residual with the decisions re-evaluated is discontinuous where one flips
  (the regime switch of zone 1), and a step that is exact for the branch it differentiates then shows no decrease however small it is (tested: stall at |step| 2e-4 for 25 iterations);
  a flip is resolved by the next iteration, whose base state records the new decisions. It is the residual in the units of the corrections, i.e. the size of the Newton step the base Jacobian would take from the trial state. s is the GMRES step.
* Trust region: the largest weighted entry of the applied step, ‖W λ s‖∞ (δ ln x for ρ, T, E_r; species weighted x/(x + xfl); u against max(|u|, c_s); r against the smaller
  adjacent zone width; ...), is limited to the radius Δ (start `delta0` = 0.5 of every Newton solve, at most `deltamax`): λ_cap = min(1, Δ/‖Ws‖∞).
* Armijo backtracking from λ_cap: halve λ until T(λ) ≤ (1 − c λ) T(0) (c = 1e-4), at most down to smin·λ_cap (smin = 1e-3); no such λ is a Newton failure (the caller retries
  with a smaller dt). A step whose weighted size is below the Newton tolerance `conv` is applied without the test (the iteration ends with it).
* Radius update from the ratio ρ_TR of the actual to the predicted reduction (T(0) − T(λ))/(λ (1 − η_lin) T(0)), η_lin the GMRES relative residual: ρ_TR < 1/4 or a backtracked step
  sets Δ to (1/2 of / the size of) the accepted step; ρ_TR > 3/4 with the step on the boundary doubles Δ.
* `maxdel!` then measures convergence on the full step s (TitanCore.STEP_FRAC = λ): a cut step is never converged.
"""
mutable struct GlobOpts
    on::Bool; delta0::Float64; deltamax::Float64; armijo::Float64; smin::Float64
end
const GLOB = GlobOpts(true, 0.5, 4.0, 1e-4, 1e-3)
mutable struct GlobCount
    solves::Int; capped::Int; backtracked::Int; fail::Int; evals::Int; rejected::Int
end
const GC = GlobCount(0, 0, 0, 0, 0, 0)
const TRR = Ref(0.5)           # current trust radius (weighted ‖·‖∞)

"The trial state in the arrays: finite, ρ, T, E_r and the species positive, r increasing and m decreasing across the interfaces."
function trial_ok()
    a, b = gi(:ngrs), gi(:ngre)
    rn, xmen, un, frn, dn, tn, ern = (vf(s) for s in (:rn, :xmen, :un, :frn, :dn, :tn, :ern)); x1, x2 = TX.XSN[][1], TX.XSN[][2]
    @inbounds for k in a:b+1
        (isfinite(rn[k]) && isfinite(xmen[k]) && isfinite(un[k]) && isfinite(frn[k])) || return false
    end
    @inbounds for k in a:b
        rn[k+1] > rn[k] && xmen[k+1] <= xmen[k] || return false
        (isfinite(dn[k]) && dn[k] > 0 && isfinite(tn[k]) && tn[k] > 0 && isfinite(ern[k]) && ern[k] > 0) || return false
        TX.cplz(k) && !(x1[k] > 0 && x2[k] > 0 && isfinite(x1[k]) && isfinite(x2[k])) && return false
    end
    true
end

"Merit T(λ) of the trial state x ⊕ λ s (Inf if the update fails, or the state is not acceptable)."
function merit(S::Sys, s, λ)
    GC.evals += 1; ok = true
    try
        R = perturbed(S, s, λ; freeze = 2, inspect = () -> (ok = trial_ok()))
        ok || (GC.rejected += 1; return Inf)
        T = levelT(S, R)[1]; isfinite(T) ? T : Inf
    catch ex
        ex isa InterruptException && rethrow()
        GC.rejected += 1; Inf
    end
end

"Fraction λ of the GMRES step s to apply (see `GlobOpts`); throws if no acceptable λ exists. `rel`: relative GMRES residual."
function globalize(S::Sys, s, rel)
    GC.solves += 1
    wx = 0.0
    @inbounds for r in 1:S.n; S.act[r] && (wx = max(wx, abs(s[r] * S.wtr[r]))); end
    isfinite(wx) || error("NKRYLOV: non-finite Newton step")
    wx <= gf(:conv) && return 1.0
    T0 = levelT(S, S.R0)[1]
    λcap = min(1.0, TRR[] / wx); λ = λcap; λcap < 1 && (GC.capped += 1)
    while λ >= GLOB.smin * λcap * (1 - 1e-12)
        T = merit(S, s, λ)
        if T <= (1 - GLOB.armijo * λ) * T0
            sz = λ * wx; ratio = (T0 - T) / (λ * (1 - min(rel, 0.99)) * T0)
            if λ < λcap
                GC.backtracked += 1; TRR[] = sz
            elseif ratio < 0.25
                TRR[] = 0.5 * sz
            elseif ratio > 0.75 && sz >= 0.9 * TRR[]
                TRR[] = min(2TRR[], GLOB.deltamax)
            end
            return λ
        end
        λ /= 2
    end
    GC.fail += 1; TRR[] = 0.5 * λ * wx
    error(@sprintf("NKRYLOV: no decrease of the merit function along the Newton step (‖Ws‖∞ %.3g, λ from %.3g to %.3g)", wx, λcap, 2λ))
end

# ---------------------------------------------------------------- the hook (TD.KRYLOV[]), called after matslv
function krylov!()
    t0 = time(); CNT.solves += 1
    ne = gi(:neqn); k0 = gi(:ngrs) - 2; nb = gi(:ngre) - gi(:ngrs) + 5
    rhs = TF.var(:rhs); blk = view(rhs, 1:ne, k0:k0+nb-1)
    δ0 = copy(blk); S = nothing
    try
        S = setup()
        OPT.check > 0 && gi(:iter) == OPT.check && (check!(S); flush(stdout); exit(0))
        b = zeros(S.n)
        @inbounds for r in 1:S.n; S.act[r] && (b[r] = δ0[(r-1)%ne+1, (r-1)÷ne+1]); end
        x, m, conv, rel = gmres(S, b)
        conv || (CNT.nonconv += 1)
        if SCAN[][1] == Trbdf2.NWATT[] && gi(:iter) in SCAN[][2]
            try scan!(S, b, x)
            catch ex
                ex isa InterruptException && rethrow()
                println("NKSCAN: diagnostic failed: ", first(sprint(showerror, ex), 300))
            end
            flush(stdout)
        end
        lsinfo = " ls off"; λ = 1.0
        if GLOB.on
            gi(:iter) == 1 && (TRR[] = GLOB.delta0)
            λ = globalize(S, x, rel)
            lsinfo = @sprintf(" glob lambda %.6g radius %.4g", λ, TRR[])
        elseif LS.on && maximum(abs, x .* S.wt) > LS.trigger
            LSC.engaged += 1
            λ, T0, Tl, ntr, sc0, scl, rel0, failed = linesearch(S, x, b)
            lsinfo = @sprintf(" ls lambda %.6g T0 %.6g T %.6g trials %d scnorm0 %.6g scnorm %.6g T0rel %.2e%s", λ, T0, Tl, ntr, sc0, scl, rel0, failed ? " lsfail" : "")
        end
        TitanCore.STEP_FRAC[] = λ
        blk .= δ0           # (the perturbed updates used `rhs` as scratch); entries update! ignores keep the band solution
        @inbounds for r in 1:S.n; S.act[r] && (blk[(r-1)%ne+1, (r-1)÷ne+1] = λ * x[r]); end
        Rf = resid(ne, nb); d = maximum(abs, Rf .- S.R0) / (maximum(abs, S.R0) + 1e-300)       # drift of the restored state
        CNT.maxdrift = max(CNT.maxdrift, d); d > 1e-10 && (CNT.ndrift += 1)
        Trbdf2.NWTRACE[] === nothing || println(Trbdf2.NWTRACE[], "NKR att ", Trbdf2.NWATT[], " it ", gi(:iter), " krylov its ", m, " converged ", conv, " rel residual ", rel, " |delta0|_W ", norm(b .* S.wt), " |delta|_W ", norm(x .* S.wt), " drift ", d, " largest weighted delta0 ", wtop(b, S, ne, k0), " largest weighted delta ", wtop(x, S, ne, k0), LS.on ? lsinfo : "")
    catch ex
        ex isa InterruptException && rethrow()
        S === nothing || restore!(S.snap)
        TX.FREEZE[] = 0; blk .= δ0; TitanCore.STEP_FRAC[] = 1.0
        GLOB.on && rethrow()          # a failed linear solve or globalisation is a Newton failure (the caller retries with a smaller dt)
        CNT.fallback += 1
        CNT.fallback <= 5 && println("NKRYLOV: fallback to the band-Jacobian correction: ", first(sprint(showerror, ex), 300))
    end
    CNT.secs += time() - t0
    nothing
end

# ---------------------------------------------------------------- matrix-vector check (NKCHECK, diagnostic)
"""
At Newton iteration `OPT.check` (with NWPROBE): for a random v over the active entries, the finite-difference J v against J_b applied to the change `update!` actually
made (the ghost entries follow the boundary conditions), per row relative to Σ_j |J_b,ij w_j|; the sensitivity to ε ×10 and ÷10; the residual of the AD pass
against the Float64 residual.
"""
function check!(S::Sys)
    ne, nb, k0 = S.ne, S.nb, S.k0; a = gi(:ngrs)
    rnd = let st = UInt64(12345); () -> (st = st * 6364136223846793005 + 1442695040888963407; (Float64(st >> 11) / 2^53) * 2 - 1); end
    v = [S.act[r] ? rnd() : 0.0 for r in 1:S.n]
    y = zeros(S.n); wr = zeros(S.n)
    ε = OPT.eps / maximum(abs, v)
    fdjv!(y, S, v, ε; realized = wr)
    z, sa = bandmul(S, wr)
    ok = isfinite.(z) .& (sa .> 0) .& .!((y .== 0) .& (z .!= 0))       # (rows that are the correction itself, R = δ: the Dirichlet rows of the boundary interfaces, are not functions of the state)
    nd = count(isfinite.(z) .& (sa .> 0)) - count(ok)
    relof(yy) = [ok[r] ? abs(yy[r] - z[r]) / sa[r] : 0.0 for r in 1:S.n]
    rel = relof(y)
    @printf("NKCHECK ne %d nb %d active %d eps %.3e; rows compared %d (%d Dirichlet rows R = delta left out); |Jv - J_b w| / sum_j |J_b,ij w_j| per row: max %.3e, median %.3e, rows > 1e-4: %d, > 1e-2: %d\n",
            ne, nb, count(S.act), ε, count(ok), nd, maximum(rel), median(rel[ok]), count(>(1e-4), rel), count(>(1e-2), rel))
    nm = ("r", "m", "rho", "u", "T", "E_r", "F_r", "x_H+", "x_H2")
    for r in sortperm(rel; rev = true)[1:30]
        kb = (r - 1) ÷ ne + 1; i = (r - 1) % ne + 1
        @printf("   row %-5s zone %3d: Jv %+.5e  J_b w %+.5e  rowscale %.3e  rel %.3e\n", i <= length(nm) ? nm[i] : string(i), k0 + kb - 1 - a + 1, y[r], z[r], sa[r], rel[r])
    end
    for f in (1e-2, 1e-1, 10.0, 100.0)
        y2 = zeros(S.n); fdjv!(y2, S, v, ε * f); r2 = relof(y2); d2 = [ok[r] ? abs(y2[r] - y[r]) / sa[r] : 0.0 for r in 1:S.n]
        @printf("NKCHECK eps x%g (max|eps v| %.1e): vs J_b w: max %.3e median %.3e rows>1e-4 %d >1e-2 %d | vs the eps run: |Jv - Jv(eps)|/|Jv(eps)| (2-norm) %.3e, worst row (rel to rowscale) %.3e\n",
                f, OPT.eps * f, maximum(r2), median(r2[ok]), count(>(1e-4), r2), count(>(1e-2), r2), norm(y2 .- y) / norm(y), maximum(d2))
    end
    rows = [argmax(rel)]
    if haskey(ENV, "NKROW")           # NKROW=eq,zone[;eq,zone...] (eq: name or number; zone counted from the inner boundary as in the listing): probe these rows instead of the worst one
        rows = Int[]
        for spec in split(ENV["NKROW"], ";")
            f = split(spec, ","); ie = something(findfirst(==(f[1]), nm), tryparse(Int, f[1]), 0); zs = parse(Int, f[2]); kbs = zs + a - k0
            (1 <= ie <= ne && 1 <= kbs <= nb) || error("NKROW=$(ENV["NKROW"]): no such row $spec")
            push!(rows, (kbs - 1) * ne + ie)
        end
    end
    @printf("NKCHECK row probe: finite-difference rows against the band rows, one unit vector per active column (%d products)\n", count(S.act))
    E = (vf2(:em2), vf2(:em1), vf2(:e00), vf2(:ep1), vf2(:ep2)); fdcols = Dict(r => Tuple{Float64,Int,Float64,Float64}[] for r in rows)
    for c in findall(S.act)
        e = zeros(S.n); e[c] = 1.0; yy = zeros(S.n); fdjv!(yy, S, e, OPT.eps)
        kc = (c - 1) ÷ ne + 1; jc = (c - 1) % ne + 1
        for rstar in rows
            kbs = (rstar - 1) ÷ ne + 1; is = (rstar - 1) % ne + 1; o = kc - kbs
            jb = abs(o) <= 2 ? E[o+3][is, jc, k0 + kbs - 1] : 0.0         # (the band row: offset kb - kc, column jc)
            push!(fdcols[rstar], (abs(yy[rstar] - jb), c, yy[rstar], jb))
        end
    end
    for rstar in rows
        kbs = (rstar - 1) ÷ ne + 1; is = (rstar - 1) % ne + 1; cols = sort!(fdcols[rstar]; rev = true)
        @printf("NKCHECK row probe: row %s of zone %d\n", is <= length(nm) ? nm[is] : string(is), k0 + kbs - 1 - a + 1)
        for (d, c, fd, jb) in cols[1:min(12, length(cols))]
            kc = (c - 1) ÷ ne + 1; jc = (c - 1) % ne + 1
            @printf("   column %-5s zone %3d (offset %+d): FD %+.4e  band %+.4e  |diff| %.3e\n", jc <= length(nm) ? nm[jc] : string(jc), k0 + kc - 1 - a + 1, kc - kbs, fd, jb, d)
        end
    end
    TD.r(:matgen)                        # AD pass at the same state: its residual against the Float64 one
    rad = vec(TF.var(:rhs)[1:ne, k0:k0+nb-1])
    @printf("NKCHECK AD residual vs Float64 residual: max |diff| / max |R| = %.3e\n", maximum(abs, rad .- S.R0) / maximum(abs, S.R0))
    nothing
end

# ---------------------------------------------------------------- switch, reporting
"NKRYLOV=1[,η[,mmax[,ε[,central]]]] (ε: target max|εv| of the scaled vector, central: 1 central / 0 one-sided differences; test parameters)."
function configure!(spec::AbstractString; activate = true)
    f = split(spec, ",")
    length(f) >= 2 && (OPT.eta = parse(Float64, f[2]))
    length(f) >= 3 && (OPT.mmax = parse(Int, f[3]))
    length(f) >= 4 && (OPT.eps = parse(Float64, f[4]))
    length(f) >= 5 && (OPT.central = f[5] != "0")
    GLOB.on = get(ENV, "NKGLOB", "1") != "0"; TitanCore.GLOBAL_NEWTON[] = GLOB.on
    activate && (TD.KRYLOV[] = krylov!)
    nothing
end
function configure_ls!(spec::AbstractString)
    f = split(spec, ","); LS.on = f[1] != "0"
    length(f) >= 2 && (LS.trigger = parse(Float64, f[2]))
    length(f) >= 3 && (LS.lmin = parse(Float64, f[3]))
    nothing
end
"Counters since the last call (`steps`: accepted steps so far), for the DTERR trace line."
function trace_text(steps)
    c = CNT; m = MARK[]; ds = c.solves - m[1]; dm = c.mv - m[2]; dn = c.nonconv - m[3]; df = c.fallback - m[4]; nsteps = steps - m[5]
    MARK[] = (c.solves, c.mv, c.nonconv, c.fallback, steps)
    @sprintf(" | NK %.2f Krylov its/Newton it, %.1f mv/step, nonconv %d, fallback %d", ds == 0 ? 0.0 : dm / ds, nsteps == 0 ? 0.0 : dm / nsteps, dn, df)
end
ls_text() = LS.on ? @sprintf(" | NKLS (trigger %g, lmin %g): %d line searches engaged, %d backtracked, %d failed, %d extra residual evaluations, max T(0) rel. diff %.2e",
                          LS.trigger, LS.lmin, LSC.engaged, LSC.backtracks, LSC.fail, LSC.evals, LSC.maxrel0) : ""
glob_text() = GLOB.on ? @sprintf(" | NKGLOB: %d globalised steps, %d trust-region capped, %d backtracked, %d failed, %d merit evaluations (%d trial states rejected)", GC.solves, GC.capped, GC.backtracked, GC.fail, GC.evals, GC.rejected) : ""
summary_text() = @sprintf("NKRYLOV (eta %g, mmax %d, eps %g): %d Newton linear solves, %d matrix-vector products (%.2f per solve), %d GMRES non-convergences, %d fallbacks, %d state-drift events (max %.2e), %.1f s",
                          OPT.eta, OPT.mmax, OPT.eps, CNT.solves, CNT.mv, CNT.solves == 0 ? 0.0 : CNT.mv / CNT.solves, CNT.nonconv, CNT.fallback, CNT.ndrift, CNT.maxdrift, CNT.secs) * ls_text() * glob_text()

end # module
