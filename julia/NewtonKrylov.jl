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

using LinearAlgebra, Printf, Statistics, ForwardDiff
import ..TitanDriver, ..TitanCore, ..TitanResidFast, ..Trbdf2, ..BlockSolver
const TD = TitanDriver; const TF = TD.TF; const TC = TitanCore; const TX = TitanResidFast
TX.EDDFAC_DUAL[] = TC.eddfac_dual
vf(s) = TF.var(s)::Vector{Float64}
gi(s) = Int(TD.g(s))
gf(s) = Float64(TD.g(s))

mutable struct Opts
    eta::Float64; mmax::Int; eps::Float64
    central::Bool                    # central (default) or one-sided differences for J v
    check::Int                       # NKCHECK=<Newton iteration>: matrix-vector check at that iteration of the first attempt, then exit (diagnostic)
    checkn::Int                      # NKCHECKN=<n>: the check is made at the n first such iterations (default 1), then exit
    ad::Bool                         # NKAD=1 (default): exact 1-partial forward-mode Dual J v instead of central finite differences
    ew::Bool                         # NKEW=1 (default): Eisenstat-Walker adaptive GMRES tolerance on non-boost iterations
end
const OPT = Opts(1e-3, 30, 1e-5, true, 0, parse(Int, get(ENV, "NKCHECKN", "1")), get(ENV, "NKAD", "1") != "0", get(ENV, "NKEW", "1") != "0")
mutable struct Count
    solves::Int; mv::Int; nonconv::Int; ndrift::Int; maxdrift::Float64; secs::Float64; jb_reused::Int
end
const CNT = Count(0, 0, 0, 0, 0.0, 0.0, 0)
const MARK = Ref((0, 0, 0, 0))        # counters and step count at the last trace line
const LUCACHE = Ref{Any}(nothing)
const LAST_SOLVE = Ref((false, 0, 0.0))

# ---------------------------------------------------------------- snapshot of everything update! writes
const SKIP = (:em2, :em1, :e00, :ep1, :ep2, :rhs)       # the Jacobian and the correction are not state
mutable struct Snap
    arrs::Dict{Symbol,Array}
    pairs1::Vector{Tuple{Vector{Float64},Vector{Float64}}}
    pairs_other::Vector{Tuple{Array,Array}}
    extra::Vector{Vector{Float64}}
    gx::Float64
    taken::Bool
    mg::Vector{Any}                  # MGROUPS: copies of the group arrays `update!` writes (mgarrs), the flag MGVEFOK
    mgflag::Bool
end
Snap() = Snap(Dict{Symbol,Array}(), Tuple{Vector{Float64},Vector{Float64}}[], Tuple{Array,Array}[], [Float64[] for _ in 1:6], 0.0, false, Any[], false)
"MGROUPS: the group arrays that `update!` writes and `restore!` must give back: E_g, F_g, the group Eddington factors, tie h and weights, h⁺, g_g, the opacity records."
mgarrs() = TX.MG[] ? vcat(Any[a for R in (TX.MGEN, TX.MGFN, TX.MGFEDD, TX.MGHEDD, TX.MGTIEW, TX.MGHPLUS) for a in R[]], Any[TX.MGGEDDR[], TX.MGVN[]]) : Any[]
const SNAP = Snap()
push!(TF.ON_SWITCH, () -> (empty!(SNAP.pairs1); empty!(SNAP.pairs_other); empty!(SNAP.arrs); SNAP.taken = false))
extras() = (TD.HEDD, TD.HXEDD, TD.FXEDD, TD.TIEW, TX.XSN[][1], TX.XSN[][2])      # closure h, h_x, f_x, tie weight; species iterate
function take!(S::Snap)
    if length(S.pairs1) + length(S.pairs_other) + length(SKIP) != length(TF.JSTATE)
        empty!(S.pairs1); empty!(S.pairs_other)
        for (n, a) in TF.JSTATE
            n in SKIP && continue
            b = haskey(S.arrs, n) && size(S.arrs[n]) == size(a) ? S.arrs[n] : (S.arrs[n] = copy(a))
            if a isa Vector{Float64} && b isa Vector{Float64}
                push!(S.pairs1, (a, b))
            else
                push!(S.pairs_other, (a, b))
            end
        end
    end
    @inbounds for (dst, src) in S.pairs1
        length(src) == length(dst) || resize!(src, length(dst))
        copyto!(src, dst)
    end
    @inbounds for (dst, src) in S.pairs_other
        copyto!(src, dst)
    end
    for (i, a) in enumerate(extras())
        b = S.extra[i]
        length(b) == length(a) || resize!(b, length(a))
        copyto!(b, a)
    end
    S.gx = TD.GXEDD[]
    if TX.MG[]
        A = mgarrs(); length(S.mg) == length(A) && all(i -> size(S.mg[i]) == size(A[i]), eachindex(A)) || (S.mg = Any[similar(a) for a in A])
        for (b, a) in zip(S.mg, A); copyto!(b, a); end
        S.mgflag = TX.MGVEFOK[]
    end
    S.taken = true
    nothing
end
function restore!(S::Snap)
    @inbounds for (dst, src) in S.pairs1
        copyto!(dst, src)
    end
    @inbounds for (dst, src) in S.pairs_other
        copyto!(dst, src)
    end
    for (b, a) in zip(extras(), S.extra)
        length(b) == length(a) || resize!(b, length(a))
        copyto!(b, a)
    end
    TD.GXEDD[] = S.gx
    if TX.MG[]
        for (a, b) in zip(mgarrs(), S.mg); copyto!(a, b); end
        TX.MGVEFOK[] = S.mgflag
    end
    nothing
end

# ---------------------------------------------------------------- the unknowns, residual, preconditioner
const LOGCODES = (0, 3, 5, 6, 8)    # equation codes of the strictly positive unknowns: species, ρ, T, E_r, E_g (MGROUPS)
const POSCODES = (3, 5, 6)          # (ρ, T, E_r: x ψ(δ); the species: x e^δ)
"A finite-difference step changes a logarithmic unknown by at most this (δ ln x): an entry of weight ~1e-12 (a trace species) otherwise gets ε v ~ 1e5, x → 0 and NaN in the EOS."
const FDCLIP = 0.1
"""
Equation codes of the unknowns by row: those of the context (1–7 the state variables, 0 the species), and with MGROUPS 8 for the groups' y_g = asinh(E_g/s_g) and 9 for their
fluxes F_g (linear, scaled by the nominal flux): rows 10 … 9 + 2G (TitanResidFast.mgiE, mgiF). `update!` applies them through MGRad.mg_update!.
"""
function nkcode(c)
    code = copy(c.code)
    if TX.MG[]
        for g in 1:TX.MGN[]; code[TX.mgiE(g)] = 8; code[TX.mgiF(g)] = 9; end
    end
    code
end
struct Sys
    ne::Int; nb::Int; k0::Int; n::Int
    c::TX.Ctx
    code::Vector{Int}; act::BitVector            # act[(kb-1)ne+i]: the entry update! applies
    wt::Vector{Float64}                         # weight of each entry in the GMRES norm (Trbdf2.wnorm's)
    wtr::Vector{Float64}                        # weight in the trust-region measure: 1 for the log variables (ρ, T, E_r, species), the GMRES weight for the others
    R0::Vector{Float64}                         # residual at the base state (Float64 path)
    snap::Snap
    AB::Matrix{Float64}; ipiv::Vector{LinearAlgebra.BlasInt}; sc::Vector{Float64}; kl::Int; ku::Int
end

const XZ_BUF = Float64[]
const R1_BUF = Ref(zeros(0, 0))
const R2_BUF = Ref(zeros(0, 0))
@inline function _get_rbuf(r::Base.RefValue{Matrix{Float64}}, ne::Int, nb::Int)
    M = r[]
    size(M) == (ne, nb) || (M = zeros(ne, nb); r[] = M)
    M
end
function resid(ne, nb; rbuf::Union{Nothing,Matrix{Float64}} = nothing)
    c = TX.context(); ws = TX.getws(Float64, c.mgr)
    R = (rbuf !== nothing && size(rbuf) == (ne, nb)) ? rbuf : zeros(ne, nb)
    fill!(R, 0.0)
    n = ne * nb
    length(XZ_BUF) == n || (resize!(XZ_BUF, n); fill!(XZ_BUF, 0.0))
    TX.residual!(R, ws, c, XZ_BUF); vec(R)
end

function activeset(c)
    ne, nb, k0 = c.neqn, c.nb, c.k0; a, b = c.ngrs, c.ngre
    act = falses(ne * nb); code = nkcode(c)
    for kb in 1:nb, i in 1:ne
        k = k0 + kb - 1; cd = code[i]
        act[(kb-1)*ne+i] = cd in (1, 2, 4, 7, 9) ? a <= k <= b + 1 :
                           cd in (3, 5, 6, 8) ? a <= k <= b :
                           (cd == 0 && i in (8, 9) && TC.CPLH[] !== nothing) ? (a <= k <= b && TX.cplz(k)) : false
    end
    act
end

"""
Weights of the correction entries in the GMRES norm: those of the DTERR error norm `Trbdf2.wnorm` (ρ, T relative; E_r against E + efl aT⁴; species against x + xfl;
u against max(|u|, c_s); m against the zone mass; r against the smaller zone width on a moving mesh; F against max(|F|, ffl c E)), so that corrections of a flux far below
c E (relative changes of 1e2 in the opaque interior) are not measured in units of that flux. The linear tolerance is a fraction of the weighted 2-norm of δ0.
Floor 1e-60: the F_r weight of the initial state (F ≈ 0, frnom ~ 1e-30) is ~1e-40 and the band step has δF ~ 1e28 there; a floor of 1e-12 made that entry the whole norm, J v = 0 for the normalised v and a singular Hessenberg matrix (SingularException at step 0).
"""
function weights(c, act)
    ne, nb, k0 = c.neqn, c.nb, c.k0; a, b = c.ngrs, c.ngre; o = Trbdf2.EOPT
    lag = gi(:lgrid) == 3 || TC.GRIDMOVE[]
    rn, xmen, dn, un, tn, unom, asn, dvn = (vf(s) for s in (:rn, :xmen, :dn, :un, :tn, :unom, :asn, :dvoln))
    ern, frn, frnom = vf(:ern), vf(:frn), vf(:frnom); wt = ones(ne * nb); code = nkcode(c)
    for kb in 1:nb, i in 1:ne
        r = (kb - 1) * ne + i; act[r] || continue; k = k0 + kb - 1; cd = code[i]
        wt[r] = cd == 3 || cd == 5 ? 1.0 :
                cd == 8 ? 1.0 :
                cd == 9 ? (g = (i - TX.MGI0[]) ÷ 2; Fg = TX.MGFN[][g][k]; TX.MGFNOM[][g][k] / max(abs(Fg), o.ffl * 2.99792458e10 * 0.5 * (ern[max(k-1, a)] + ern[min(k, b)]), TX.MGEFL * abs(frn[k]))) :
                cd == 6 ? ern[k] / (ern[k] + o.efl * 7.565723e-15 * tn[k]^4) :
                cd == 0 ? (q = TX.XSN[][i-7][k]; q / (q + o.xfl)) :
                cd == 4 ? unom[k] / max(abs(un[k]), max(asn[max(k-1, a)], asn[min(k, b)])) :
                cd == 2 ? xmen[k] / max(abs(xmen[k]), dn[max(k-1, a)] * dvn[max(k-1, a)], dn[min(k, b)] * dvn[min(k, b)]) :
                cd == 1 ? (lag ? rn[k] / min(rn[k] - rn[max(k-1, a)], rn[min(k+1, b+1)] - rn[k]) : 1.0) :
                frnom[k] / max(abs(frn[k]), o.ffl * 2.99792458e10 * 0.5 * (ern[max(k-1, a)] + ern[min(k, b)]))
    end
    [isfinite(x) ? max(abs(x), 1e-60) : 1.0 for x in wt]       # (a 0/0 at a boundary interface: weight 1)
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
    # An unknown no row depends on (an all-zero column: the ghost entries of the first step, whose boundary rows do not reference them) makes the band matrix exactly singular
    # (SingularException(1) at step 0 from the IC); P projects such entries out of every vector, so a unit on the diagonal changes nothing in the preconditioner's action on the active ones.
    @inbounds for cc in 1:n
        z = true
        for r in max(1, cc - ku):min(n, cc + kl); AB[kl+ku+1+r-cc, cc] != 0 && (z = false; break); end
        z && (AB[kl+ku+1, cc] = 1.0)
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
function perturbed(S::Sys, v, s; realized = nothing, freeze = 2, inspect = nothing, clip = Inf, rbuf = nothing)
    S.snap.taken || take!(S.snap)
    ne, nb, k0 = S.ne, S.nb, S.k0; rhs = TF.var(:rhs)
    base = realized === nothing ? nothing : (c = TX.context(); (copy(c.base), copy(vf(:unom)), copy(vf(:frnom))))
    @inbounds for kb in 1:nb, i in 1:ne; rhs[i, k0+kb-1] = S.act[(kb-1)*ne+i] ? (S.code[i] in LOGCODES ? clamp(s * v[(kb-1)*ne+i], -clip, clip) : s * v[(kb-1)*ne+i]) : 0.0; end
    TX.FREEZE[] = freeze; updated = false
    try
        TD.r(:update); updated = true
        inspect === nothing || inspect()          # (NKSCAN: the state after the trial update, before the restore)
        R1 = resid(ne, nb; rbuf)
        if realized !== nothing
            cnew = TX.context()
            for kb in 1:nb, i in 1:ne
                cd = S.code[i]; k = k0 + kb - 1; r = (kb - 1) * ne + i
                (cd == 0 || cd >= 8) && (realized[r] = S.act[r] ? v[r] : 0.0; continue)       # species, group fields: applied as given
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
Directional derivative of the full residual, y = J v: central difference [R(x ⊕ εv) − R(x ⊖ εv)]/2ε (one-sided [R(x ⊕ εv) − R(x)]/ε if `OPT.central` is false). v is a vector over the active entries.
"""
function fdjv!(y, S::Sys, v, ε; realized = nothing, freeze = 2)
    R1 = perturbed(S, v, ε; realized, freeze, clip = FDCLIP, rbuf = _get_rbuf(R1_BUF, S.ne, S.nb))
    if OPT.central
        R2 = perturbed(S, v, -ε; freeze, clip = FDCLIP, rbuf = _get_rbuf(R2_BUF, S.ne, S.nb))
        @inbounds for r in 1:S.n; y[r] = (R1[r] - R2[r]) / (2ε); end
    else
        @inbounds for r in 1:S.n; y[r] = (R1[r] - S.R0[r]) / ε; end
    end
    y
end

struct NKTag end
ForwardDiff.:≺(::Type{NKTag}, ::Type{<:ForwardDiff.Tag}) = true
ForwardDiff.:≺(::Type{<:ForwardDiff.Tag}, ::Type{NKTag}) = false
"Exact 1-partial forward-mode directional derivative y = J v at the base state."
function adjv!(y, S::Sys, v)
    D = ForwardDiff.Dual{NKTag,Float64,1}
    xd = get!(() -> Vector{D}(undef, S.n), TX.DBUF, (:nk_x, D, S.n))::Vector{D}
    Rd = get!(() -> Matrix{D}(undef, S.ne, S.nb), TX.DBUF, (:nk_R, D, S.ne, S.nb))::Matrix{D}
    @inbounds for r in 1:S.n
        xd[r] = D(0.0, ForwardDiff.Partials((S.act[r] ? v[r] : 0.0,)))
    end
    ws = TX.getws(D, S.c.mgr)
    TX.AD_JV[] = true; TX.FREEZE[] = 2
    try
        TX.residual!(Rd, ws, S.c, xd)
    finally
        TX.AD_JV[] = false; TX.FREEZE[] = 0
    end
    @inbounds for r in 1:S.n
        y[r] = ForwardDiff.partials(Rd[r], 1)
    end
    y
end

function setup()
    TF.BACKEND[] === :julia || error("NKRYLOV needs the Julia state backend")
    c = TX.context(); ne, nb, k0 = c.neqn, c.nb, c.k0
    act = activeset(c); snap = SNAP; snap.taken = false; OPT.ad || take!(snap)
    TX.FREEZE[] = 1
    R0 = try resid(ne, nb) finally TX.FREEZE[] = 0 end
    AB, ipiv, sc, kl, ku = if TX.RESID_ONLY[] && LUCACHE[] !== nothing
        CNT.jb_reused += 1
        LUCACHE[]::Tuple{Matrix{Float64},Vector{LinearAlgebra.BlasInt},Vector{Float64},Int,Int}
    else
        lu_ = bandlu(ne, nb, k0)
        LUCACHE[] = lu_
        lu_
    end
    wt = weights(c, act)
    code = nkcode(c); wtr = [code[(r-1)%ne+1] in LOGCODES ? 1.0 : wt[r] for r in 1:ne*nb]
    Sys(ne, nb, k0, ne * nb, c, code, act, wt, wtr, R0, snap, AB, ipiv, sc, kl, ku)
end

# ---------------------------------------------------------------- GMRES
"Solve A x = b (A v = P J_b⁻¹ J v) from x = 0 by GMRES in the norm ‖W·‖₂ (weights S.wt); returns (x, iterations, converged, relative residual, truncated). With a finite `radius` the iteration stops (truncated Steihaug-style) at the first iterate whose step is outside the trust region (`stepsize` > radius) and keeps the one before."
function gmres(S::Sys, b, radius)
    n = S.n; mmax = OPT.mmax; wt = S.wt; bs = b .* wt           # scaled space: ṽ = W v
    η = (OPT.ew && gi(:lboos) == 0) ? (gi(:iter) == 1 ? max(OPT.eta, 1e-2) : clamp(0.1 * gf(:dmax), OPT.eta, 1e-2)) : OPT.eta
    β = norm(bs); x = zeros(n)
    β == 0 && return x, 0, true, 0.0, false
    V = zeros(n, mmax + 1); H = zeros(mmax + 1, mmax); cs = zeros(mmax); sn = zeros(mmax); gv = zeros(mmax + 1)
    V[:, 1] = bs ./ β; gv[1] = β
    y = zeros(n); w = zeros(n); v = zeros(n); m = 0; conv = false; res = β; trunc = false; mgood = 1; ygood = Float64[0.0]; xprev = zeros(n); xcut = zeros(n)
    for j in 1:mmax
        v .= (@view V[:, j]) ./ wt; v[.!S.act] .= 0.0
        OPT.ad ? adjv!(y, S, v) : fdjv!(y, S, v, OPT.eps / maximum(abs, @view V[:, j])); CNT.mv += 1
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
        # trust-region truncation (Steihaug): the iterate j is the combination V y with the GMRES coefficients y; stop at the first one outside the radius and keep the one before it
        # (the first iterate is kept and scaled by the line search): the later iterates amplify the near-null modes of A = P J_b⁻¹ J that the band step b does not contain
        if radius < Inf
            yj = UpperTriangular(H[1:j, 1:j]) \ gv[1:j]; xt = (@view V[:, 1:j]) * yj; xt ./= wt; xt[.!S.act] .= 0.0
            if j >= 2 && (!all(isfinite, xt) || stepsize(S, boxclip!(S, copy(xt)), 1.0) > radius)
                trunc = true
                # Steihaug–Toint (Steihaug 1983; Nocedal & Wright, Alg. 7.2): the first iterate that leaves the trust region is replaced by the point where the segment
                # from the last iterate inside it crosses the boundary (bisection on the step measure), not by the last iterate itself
                if all(isfinite, xt)
                    lo_, hi_ = 0.0, 1.0
                    for _ in 1:40
                        mid_ = 0.5 * (lo_ + hi_)
                        stepsize(S, boxclip!(S, xprev .+ mid_ .* (xt .- xprev)), 1.0) > radius ? (hi_ = mid_) : (lo_ = mid_)
                    end
                    xcut = xprev .+ lo_ .* (xt .- xprev)
                else
                    xcut = copy(xprev)
                end
                break
            end
            mgood = j; ygood = yj; xprev = copy(xt)
        end
        V[:, j+1] = w ./ hn
    end
    if trunc
        m = mgood; x .= xcut; x[.!S.act] .= 0.0
    else
        yk = UpperTriangular(H[1:m, 1:m]) \ gv[1:m]
        for j in 1:m; x .+= yk[j] .* @view V[:, j]; end
        x ./= wt; x[.!S.act] .= 0.0
    end
    x, m, trunc ? false : conv, (trunc ? abs(gv[m+1]) : res) / β, trunc
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
const SCANPREV = Ref{Any}(nothing)
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
    lams = haskey(ENV, "NKSCANX") ? collect(0.5:0.05:1.0) : [1e-6, 1e-5, 1e-4, 1e-3, 1/256, 1/128, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1.0]
    r937 = (a + 37 - 1 - k0) * ne + 9
    for sg in (-1, 1), l in lams
        λ = sg * l; line = ""; (haskey(ENV, "NKSCANX") && sg < 0) && continue
        for fr in (0, 2)
            R, exc, bad = scantrial(S, x, λ, fr)
            if exc !== nothing
                println(pre, "C lambda ", λ, " freeze ", fr, " EXCEPTION ", typeof(exc), ": ", first(split(sprint(showerror, exc), '\n')))
            else
                T, w = levelT(S, R)
                if fr == 0 && haskey(ENV, "NKSCANX")
                    SCANPREV[] !== nothing && println(pre, "C   lambda ", λ, " rows jumping vs previous lambda (sc R diff): ", toplist(S, (S.sc .* R) .- SCANPREV[], 8), " | R_sc now: ", join((@sprintf("%s:%.3e", zname(S, r), S.sc[r] * R[r]) for r in sortperm(abs.(S.sc .* R); rev = true)[1:6]), " "))
                    SCANPREV[] = copy(S.sc .* R)
                end
                if isfinite(T)
                    @printf("%sC lambda %+.6g freeze %d T %.6g T/T0 %.6g linear |1-l| %.6g |sc R| %.6g worst row %s | R_sc[9@37] %+.4e P R[9@37] %+.4e P R[5@37] %+.4e P R[8@37] %+.4e\n", pre, λ, fr, T, T / T0, abs(1 - λ), norm(S.sc .* R), zname(S, argmax(abs.(w))), S.sc[r937] * R[r937], w[r937] / S.wt[r937], w[r937-4] / S.wt[r937-4], w[r937-1] / S.wt[r937-1])
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
    # E: raw column against the band column, plus J v vs J_b v on v = x/|x|_W and top 5 columns of v
    for (lbl, vec_) in vcat([("v=x/|x|_W", v)], [(zname(S, rk), (ek = zeros(S.n); ek[rk] = 1.0; ek)) for rk in sortperm(abs.(wx); rev = true)[1:5]])
        fdjv!(y, S, vec_, 1e-5 / maximum(abs, vec_ .* S.wt)); z, _ = bandmul(S, vec_)
        d = abs.(y .- z) .* S.sc
        println(pre, "E rows with largest row-scaled |J u - J_b u| (u = ", lbl, "; values times S.sc): row, J u, J_b u, diff")
        for r in sortperm(d; rev = true)[1:8]
            @printf("%sE [%s] %s: Ju %+.5e Jb u %+.5e diff %+.5e (sc %.3e)\n", pre, lbl, zname(S, r), y[r] * S.sc[r], z[r] * S.sc[r], (y[r] - z[r]) * S.sc[r], S.sc[r])
        end
    end
    nothing
end

# ---------------------------------------------------------------- globalisation: trust region + Armijo line search on the weighted true residual
"""
Newton globalisation (always on with NKRYLOV; docs/POPIII_SWITCH_AUDIT2.md). One method for every unknown:

* Positivity by construction: `update!` applies ρ, T, E_r (and the species) as x ← x e^δ with δ the relative Newton correction (the step in ln x), r and m as x(1 + δ);
  a trial state that is non-finite, has x ≤ 0, or has interfaces that cross (r not increasing, m not decreasing) has merit Inf (`trial_ok`).
* Merit function T(λ) = ‖W J_b⁻¹ R(x ⊕ λ s)‖₂ (J_b the band Jacobian of the iteration's base state, W the GMRES weights): the residual of all rows with the base state's branch
  decisions (limiter, sign, upwind, viscosity) replayed, as in J v and J_b: the model the step solves. The residual with the decisions re-evaluated is discontinuous where one flips
  (the regime switch of zone 1), and a step that is exact for the branch it differentiates then shows no decrease however small it is (tested: stall at |step| 2e-4 for 25 iterations);
  a flip is resolved by the next iteration, whose base state records the new decisions. It is the residual in the units of the corrections, i.e. the size of the Newton step the base Jacobian would take from the trial state. s is the GMRES step.
* Trust region: the largest entry of the applied step λ s, measured as ln ψ(λ s) for ρ, T, E_r and |λ s| for the species (both δ ln x), and in the GMRES weights for the others (u against
  max(|u|, c_s), F_r against cE, m against the zone mass, r against the zone width), is limited to the radius Δ (start `delta0` = 0.5 in every Newton solve, at most `deltamax`):
  λ_cap = min(1, largest λ with size ≤ Δ) (`capfraction`).
* The step is first truncated componentwise to the box |δ ln x| ≤ `deltamax` for the logarithmic unknowns (`boxclip!`).
* Armijo backtracking from λ_cap: halve λ until T(λ) ≤ (1 − c λ + ν) T(0) (c = 1e-4; ν = 1e-6 the round-off of the residual: a step that changes only entries the weights all but ignore leaves T(λ) = T(0) to 9 digits and must not be a failure), at most down to smin·λ_cap (smin = 1e-3); no such λ is a Newton failure (the caller retries
  with a smaller dt). A step of size ≤ `tol` = 1e-4 is applied without the test: the merit cannot resolve a decrease there (round-off of the formal solution is 1e-6 in R and is amplified by J_b⁻¹; with the test every step of 1e-5 backtracked and shrank the radius), and it is the quadratic regime.
* Radius update from the ratio ρ_TR of the actual to the predicted reduction (T(0) − T(λ))/(λ (1 − η_lin) T(0)), η_lin the GMRES relative residual: ρ_TR < 1/4 or a backtracked step
  sets Δ to (1/2 of / the size of) the accepted step; ρ_TR > 1/2 with the step on the boundary doubles Δ.
* `maxdel!` then measures convergence on the full step s (TitanCore.STEP_FRAC = λ): a cut step is never converged.
"""
mutable struct GlobOpts
    delta0::Float64; deltamax::Float64; armijo::Float64; smin::Float64; tol::Float64; noise::Float64; rmin::Float64
end
const GLOB = GlobOpts(0.5, 4.0, 1e-4, 1e-3, 1e-4, 1e-6, 0.05)
mutable struct GlobCount
    solves::Int; capped::Int; backtracked::Int; fail::Int; evals::Int; rejected::Int; truncated::Int
end
const GC = GlobCount(0, 0, 0, 0, 0, 0, 0)
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


"Size of the step λ s in the trust-region measure: ln ψ(λ s_i) (`TitanCore.posfac`) for ρ, T, E_r, |λ s_i| for the species (their update is e^δ), weighted |λ s_i| for the others."
function stepsize(S::Sys, s, λ)
    m = 0.0
    @inbounds for r in 1:S.n
        S.act[r] || continue
        c = S.code[(r-1)%S.ne+1]; d = λ * s[r]
        m = max(m, c in POSCODES ? (d >= 0 ? log1p(d) : -log(TitanCore.posfac(d))) : (c == 0 || c == 8) ? abs(d) : abs(d) * S.wtr[r])
    end
    m
end
"Largest λ with `stepsize(S, s, λ) ≤ Δ`."
function capfraction(S::Sys, s, Δ)
    f = Inf
    @inbounds for r in 1:S.n
        S.act[r] || continue
        c = S.code[(r-1)%S.ne+1]; d = s[r]; d == 0 && continue
        f = min(f, c in POSCODES ? (d > 0 ? expm1(Δ) / d : TitanCore.psiinv(exp(-Δ)) / d) : (c == 0 || c == 8) ? Δ / abs(d) : Δ / (abs(d) * S.wtr[r]))
    end
    f
end

"""
Box on the logarithmic unknowns: the step of ρ, T, E_r and the species is truncated to ln ψ(δ) ≤ `deltamax` (growth) and δ ≥ −`deltamax` (decay) componentwise, before the radius scales the whole vector.
A uniform scaling cannot handle a step whose size is set by entries the norm all but ignores: a trace species (x ~ 1e-30) asks for δ ln x ~ 1e12 and λ_cap = Δ/1e12 freezes every other variable;
the merit is blind to those entries (weight ~1e-12), so truncating them leaves the descent property of the rest. (This replaces CPLCLAMP, which clipped the species alone.)
"""
function boxclip!(S::Sys, s)
    hi = expm1(GLOB.deltamax); lo = TitanCore.psiinv(exp(-GLOB.deltamax))
    @inbounds for r in 1:S.n
        S.act[r] || continue
        c = S.code[(r-1)%S.ne+1]
        c in POSCODES ? (s[r] = clamp(s[r], lo, hi)) : (c == 0 || c == 8) ? (s[r] = clamp(s[r], -GLOB.deltamax, GLOB.deltamax)) : nothing
    end
    s
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

"Fraction λ of the GMRES step s to apply (trust-region cap without Armijo backtracking)."
function globalize(S::Sys, s, rel)
    GC.solves += 1
    wx = stepsize(S, s, 1.0)
    isfinite(wx) || error("NKRYLOV: non-finite Newton step")
    λ = min(1.0, capfraction(S, s, TRR[]))
    λ < 1 && (GC.capped += 1)
    return λ
end

# ---------------------------------------------------------------- the hook (TD.KRYLOV[]), called after matslv
function krylov!()
    t0 = time(); CNT.solves += 1
    ne = gi(:neqn); k0 = gi(:ngrs) - 2; nb = gi(:ngre) - gi(:ngrs) + 5
    rhs = TF.var(:rhs); blk = view(rhs, 1:ne, k0:k0+nb-1)
    δ0 = copy(blk); S = nothing
    try
        S = setup()
        if OPT.check > 0 && gi(:iter) == OPT.check
            check!(S); flush(stdout); OPT.checkn -= 1
            OPT.checkn <= 0 && exit(0)
        end
        b = zeros(S.n)
        @inbounds for r in 1:S.n; S.act[r] && (b[r] = δ0[(r-1)%ne+1, (r-1)÷ne+1]); end
        gi(:iter) == 1 && (TRR[] = GLOB.delta0)
        x, m, conv, rel, trunc = gmres(S, b, TRR[])
        if TX.RESID_ONLY[] && (!conv || m > 6)
            # the cached band factorisation no longer preconditions J at this state (the estimate's solves run at a perturbed old state):
            # refactorise at the current state (matgen + matslv with a fresh LU, which also replaces the cache) and solve again
            CNT.solves -= 1; CNT.jb_reused -= 1
            S.snap.taken && restore!(S.snap)
            TX.FREEZE[] = 0; TitanCore.STEP_FRAC[] = 1.0
            TX.RESID_ONLY[] = false; BlockSolver.Penta.REUSE_LU[] = false
            TitanDriver.r(:matgen); TitanDriver.r(:matslv)
            krylov!()
            TX.RESID_ONLY[] = true; BlockSolver.Penta.REUSE_LU[] = true       # (the caller's later solves reuse the fresh factorisation)
            return nothing
        end
        conv || trunc || (CNT.nonconv += 1); trunc && (GC.truncated += 1)
        if SCAN[][1] == Trbdf2.NWATT[] && gi(:iter) in SCAN[][2]
            try scan!(S, b, x)
            catch ex
                ex isa InterruptException && rethrow()
                println("NKSCAN: diagnostic failed: ", first(sprint(showerror, ex), 300))
            end
            flush(stdout)
        end
        boxclip!(S, x)
        λ = globalize(S, x, rel)
        LAST_SOLVE[] = (conv, m, λ); Trbdf2.LAST_KRYLOV_SOLVE[] = (conv, m, λ)
        lsinfo = @sprintf(" glob lambda %.6g radius %.4g", λ, TRR[])
        TitanCore.STEP_FRAC[] = λ
        blk .= δ0           # (the perturbed updates used `rhs` as scratch); entries update! ignores keep the band solution
        @inbounds for r in 1:S.n; S.act[r] && (blk[(r-1)%ne+1, (r-1)÷ne+1] = λ * x[r]); end
        d = if S.snap.taken
            Rf = resid(ne, nb)
            d_ = maximum(abs, Rf .- S.R0) / (maximum(abs, S.R0) + 1e-300)       # drift of the restored state
            CNT.maxdrift = max(CNT.maxdrift, d_); d_ > 1e-10 && (CNT.ndrift += 1)
            d_
        else
            0.0
        end
        OPT_DIAG && diag!(S, b, x, λ, rel, m, conv, trunc)
        Trbdf2.NWTRACE[] === nothing || println(Trbdf2.NWTRACE[], "NKR att ", Trbdf2.NWATT[], " it ", gi(:iter), " krylov its ", m, " converged ", conv, " rel residual ", rel, " |delta0|_W ", norm(b .* S.wt), " |delta|_W ", norm(x .* S.wt), " drift ", d, " largest weighted delta0 ", wtop(b, S, ne, k0), " largest weighted delta ", wtop(x, S, ne, k0), lsinfo)
    catch ex
        ex isa InterruptException && rethrow()
        (S === nothing || !S.snap.taken) || restore!(S.snap)
        TX.FREEZE[] = 0; blk .= δ0; TitanCore.STEP_FRAC[] = 1.0
        rethrow()          # a failed linear solve or globalisation is a Newton failure (the caller retries with a smaller dt), not a silent band solve
    end
    CNT.secs += time() - t0
    nothing
end


# ---------------------------------------------------------------- NKDIAG (diagnostic, default off; no behaviour change)
const OPT_DIAG = get(ENV, "NKDIAG", "") != ""
const DIAGZ = [z for z0 in split(get(ENV, "NKDIAGZ", "35:39,99:101"), ",") for z in (occursin(":", z0) ? collect(parse(Int, split(z0, ":")[1]):parse(Int, split(z0, ":")[2])) : [parse(Int, z0)])]
function diag!(S::Sys, b, x, λ, rel, m, conv, trunc)
    io = Trbdf2.NWTRACE[] === nothing ? stdout : Trbdf2.NWTRACE[]; ne = S.ne; a = gi(:ngrs)
    zn(kb) = S.k0 + kb - 1 - a + 1
    fcap = Inf; rcap = 0
    @inbounds for r in 1:S.n
        S.act[r] || continue
        c = S.code[(r-1)%ne+1]; d = x[r]; d == 0 && continue
        f = c in POSCODES ? (d > 0 ? expm1(TRR[]) / d : TitanCore.psiinv(exp(-TRR[])) / d) : (c == 0 || c == 8) ? TRR[] / abs(d) : TRR[] / (abs(d) * S.wtr[r])
        f < fcap && (fcap = f; rcap = r)
    end
    Rs = S.R0 .* S.sc
    @printf(io, "NKD att %d it %d jstep %d dt %.5g | lam %.4g radius %.4g cap-entry %s (x %.3g wtr %.2g) | gmres its %d conv %s trunc %s rel %.2e | |R|_sc inf %.3e %s | |b|_W %.3e |x|_W %.3e\n",
            Trbdf2.NWATT[], gi(:iter), gi(:jstep), TD.g(:dtime), λ, TRR[], rcap == 0 ? "-" : zname(S, rcap), rcap == 0 ? 0.0 : x[rcap], rcap == 0 ? 0.0 : S.wtr[rcap], m, conv, trunc, rel,
            maximum(abs, Rs), zname(S, argmax(abs.(Rs))), norm(b .* S.wt), norm(x .* S.wt))
    for i in 1:ne
        rr = [(kb-1)*ne+i for kb in 1:S.nb if S.act[(kb-1)*ne+i]]; isempty(rr) && continue
        rb = rr[argmax(abs.(b[rr]))]; rx = rr[argmax(abs.(x[rr] .* S.wt[rr]))]; rR = rr[argmax(abs.(Rs[rr]))]
        @printf(io, "NKD   eq %-4s max|R_sc| %+.3e @%d | max|band corr| %+.3e @%d | max W|gmres step| %+.3e @%d (raw %+.3e)\n", NMS[min(i, 9)], Rs[rR], zn((rR-1)÷ne+1), b[rb], zn((rb-1)÷ne+1), x[rx] * S.wt[rx], zn((rx-1)÷ne+1), x[rx])
    end
    for z in DIAGZ
        kb = z + a - S.k0; (1 <= kb <= S.nb) || continue
        @printf(io, "NKD   zone %3d R_sc:", z); for i in 1:ne; @printf(io, " %+.2e", Rs[(kb-1)*ne+i]); end
        @printf(io, " | step:"); for i in 1:ne; @printf(io, " %+.2e", x[(kb-1)*ne+i]); end; println(io)
    end
    flush(io)
end

# ---------------------------------------------------------------- matrix-vector check (NKCHECK, diagnostic)
"""
At Newton iteration `OPT.check` (with NWPROBE): for a random v over the active entries, the finite-difference J v against J_b applied to the change `update!` actually
made (the ghost entries follow the boundary conditions), per row relative to Σ_j |J_b,ij w_j|; the sensitivity to ε ×10 and ÷10; the residual of the AD pass
against the Float64 residual.
"""
function check!(S::Sys)
    ne, nb, k0 = S.ne, S.nb, S.k0; a = gi(:ngrs)
    # The finite differences leave the lagged Float data of the residual (what the last update! wrote and the restore does not give back: opacities and star-light
    # state of the last trial state) at x − εv: the exact product is evaluated after a zero-perturbation update!, i.e. at the state the solver's own products see.
    adj!(y, S_, v_) = (perturbed(S_, zeros(S_.n), 0.0); adjv!(y, S_, v_))
    rnd = let st = UInt64(12345); () -> (st = st * 6364136223846793005 + 1442695040888963407; (Float64(st >> 11) / 2^53) * 2 - 1); end
    v = [S.act[r] ? rnd() : 0.0 for r in 1:S.n]
    y = zeros(S.n); wr = zeros(S.n)
    ε = OPT.eps / maximum(abs, v)
    # NKXS (default 1e-2): the species perturbations of the finite difference stay within NKXS × the neutral fraction n_H/n_H,tot = 1 − x_H⁺ − 2x_H₂: the species are x e^δ
    # and the neutral fraction of the ionised gas (1 − x_H⁺ ~ 1e-6 … 1e-9) is a difference, so a relative step of 1e-5 in x_H⁺ leaves the linear range of the derivative that the
    # check compares (the difference converges to the AD product as the step shrinks: NKROWPROBE=2). NKXS=0: no restriction.
    xsc = parse(Float64, get(ENV, "NKXS", "1e-2"))
    function xcap(r, vr, ε)
        xsc > 0 && (r - 1) % ne + 1 in (8, 9) && S.code[(r-1)%ne+1] == 0 || return vr
        k = k0 + (r - 1) ÷ ne; x1 = TX.XSN[][1][k]; x2 = TX.XSN[][2][k]; xi = (r - 1) % ne + 1 == 8 ? x1 : 2x2
        xi > 0 || return vr
        sign(vr) * min(abs(vr), xsc * max(1 - x1 - 2x2, 1e-30) / (xi * ε))
    end
    for r in 1:S.n; S.act[r] && (v[r] = xcap(r, v[r], ε)); end
    fdjv!(y, S, v, ε; realized = wr)
    z, sa = bandmul(S, wr)
    ok = isfinite.(z) .& (sa .> 0) .& .!((y .== 0) .& (z .!= 0))       # (rows that are the correction itself, R = δ: the Dirichlet rows of the boundary interfaces, are not functions of the state)
    nd = count(isfinite.(z) .& (sa .> 0)) - count(ok)
    relof(yy) = [ok[r] ? abs(yy[r] - z[r]) / sa[r] : 0.0 for r in 1:S.n]
    rel = relof(y)
    @printf("NKCHECK ne %d nb %d active %d eps %.3e; rows compared %d (%d Dirichlet rows R = delta left out); |Jv - J_b w| / sum_j |J_b,ij w_j| per row: max %.3e, median %.3e, rows > 1e-4: %d, > 1e-2: %d\n",
            ne, nb, count(S.act), ε, count(ok), nd, maximum(rel), median(rel[ok]), count(>(1e-4), rel), count(>(1e-2), rel))
    nm = ("r", "m", "rho", "u", "T", "E_r", "F_r", "x_H+", "x_H2")
    rr1 = zeros(S.n)
    if OPT.ad                          # the exact forward-mode product against the central finite difference through update! (same v): per row relative to the band row scale
        ya = zeros(S.n); adj!(ya, S, v)
        let ya2 = zeros(S.n); adj!(ya2, S, v); ya3 = zeros(S.n); adj!(ya3, S, v)         # (the product must not depend on how often it has been evaluated: repeat it)
            d12 = [ok[r] ? abs(ya[r] - ya2[r]) / sa[r] : 0.0 for r in 1:S.n]; d23 = [ok[r] ? abs(ya2[r] - ya3[r]) / sa[r] : 0.0 for r in 1:S.n]
            wr_ = sortperm(d12; rev = true)[1:5]
            @printf("NKCHECK AD repeat: 1st vs 2nd call of J v: max %.3e (rows > 1e-4: %d); 2nd vs 3rd: max %.3e; worst rows of the first pair: %s\n", maximum(d12), count(>(1e-4), d12), maximum(d23),
                    join((@sprintf("%s z%d %.2e", (i = (r - 1) % ne + 1; i <= length(nm) ? nm[i] : string(i)), k0 + (r - 1) ÷ ne - a + 1, d12[r]) for r in wr_), "; "))
        end
        worst1 = Int[]
        for ef in (1.0, 0.1, 10.0, 100.0, 1e-2, 1e-3, 1e-4)      # (the FD at ε, ε/10, 10ε, … ε/1e4: the truncation and round-off of the difference itself; a kink of the residual at the iterate does not go away)
            yf = ef == 1.0 ? y : (yy = zeros(S.n); fdjv!(yy, S, v, ε * ef); yy)
            rr = [ok[r] ? abs(ya[r] - yf[r]) / sa[r] : 0.0 for r in 1:S.n]
            @printf("NKCHECK AD J v vs FD (eps x%g): per row |J v(AD) - J v(FD)| / row scale: max %.3e, median %.3e, rows > 1e-4: %d, > 1e-2: %d | 2-norm rel %.3e\n", ef, maximum(rr), median(rr[ok]),
                    count(>(1e-4), rr), count(>(1e-2), rr), norm(ya .- yf) / max(norm(yf), 1e-300))
            ef == 1.0 && (worst1 = sortperm(rr; rev = true)[1:6])
            ef == 1.0 || @printf("NKCHECK AD vs FD (eps x%g): the worst rows of the eps x1 list: %s\n", ef, join((@sprintf("%s z%d %.2e", (i = (r - 1) % ne + 1; i <= length(nm) ? nm[i] : string(i)), k0 + (r - 1) ÷ ne - a + 1, rr[r]) for r in worst1), "; "))
            if ef == 1.0                                 # the same outside the front zones NKZX (default 36:46: photosphere, H₂ and ionisation fronts, where the state is not smooth)
                zx1, zx2 = parse.(Int, split(get(ENV, "NKZX", "36:46"), ":"))
                out = [ok[r] && !(zx1 <= k0 + (r - 1) ÷ ne - a + 1 <= zx2) for r in 1:S.n]
                if any(out)
                    ro = [out[r] ? rr[r] : 0.0 for r in 1:S.n]
                    @printf("NKCHECK AD vs FD outside zones %d:%d: %d rows, max %.3e, median %.3e, rows > 1e-4: %d, > 1e-2: %d\n", zx1, zx2, count(out), maximum(ro), median(rr[out]), count(>(1e-4), ro), count(>(1e-2), ro))
                    for i in 1:ne
                        ri = [(kb - 1) * ne + i for kb in 1:nb]; kall = argmax([ro[r] for r in ri])
                        nmi = i <= length(nm) ? nm[i] : (i - 9) % 2 == 1 ? "E_g$((i - 9 + 1) ÷ 2)" : "F_g$((i - 9) ÷ 2)"
                        @printf("NKCHECK AD vs FD outside eq %-6s: worst row %.3e (zone %d; J v AD %+.5e FD %+.5e)\n", nmi, ro[ri[kall]], k0 + kall - 1 - a + 1, ya[ri[kall]], yf[ri[kall]])
                    end
                end
            end
            ef == 1.0 || continue
            rr1 = rr
            for i in 1:ne
                ri = [(kb - 1) * ne + i for kb in 1:nb]; kall = argmax([rr[r] for r in ri]); zn(kb) = k0 + kb - 1 - a + 1
                nmi = i <= length(nm) ? nm[i] : (i - 9) % 2 == 1 ? "E_g$((i - 9 + 1) ÷ 2)" : "F_g$((i - 9) ÷ 2)"
                @printf("NKCHECK AD vs FD eq %-6s: worst row %.3e (zone %d; J v AD %+.5e FD %+.5e; row scale %.3e)\n", nmi, rr[ri[kall]], zn(kall), ya[ri[kall]], yf[ri[kall]], sa[ri[kall]])
            end
        end
    end
    if get(ENV, "NKROWPROBE", "1") == "2" && OPT.ad       # AD against FD column by column (the exact product's own columns): which entries of the worst AD-vs-FD rows disagree
        zp1, zp2 = get(ENV, "NKPZ", "34:48") == "all" ? (1, 0) : parse.(Int, split(get(ENV, "NKPZ", "34:48"), ":")); zq(kb) = k0 + kb - 1 - a + 1
        rows = Int[]; for r in sortperm(rr1; rev = true); length(rows) >= parse(Int, get(ENV, "NKPROWS", "6")) && break; rr1[r] > 1e-4 && push!(rows, r); end
        if haskey(ENV, "NKROW")                  # NKROW=eq,zone[;eq,zone...] (as the band probe): these rows too
            for spec in split(ENV["NKROW"], ";")
                f = split(spec, ","); ie = something(findfirst(==(f[1]), nm), tryparse(Int, f[1]), 0); zs = parse(Int, f[2]); kbs = zs + a - k0
                (1 <= ie <= ne && 1 <= kbs <= nb) || error("NKROW=$(ENV["NKROW"]): no such row $spec")
                push!(rows, (kbs - 1) * ne + ie)
            end
        end
        let ya = zeros(S.n); adj!(ya, S, v)           # one-sided differences along v: a kink of the residual at the iterate shows as forward ≠ backward
            R1 = copy(perturbed(S, v, ε)); R2 = copy(perturbed(S, v, -ε))
            for r in rows
                i = (r - 1) % ne + 1; @printf("NKCHECK AD/FD one-sided: row %s zone %d: AD %+.5e  forward %+.5e  backward %+.5e  central %+.5e (row scale %.3e)\n", i <= length(nm) ? nm[i] : string(i), zq((r - 1) ÷ ne + 1),
                        ya[r], (R1[r] - S.R0[r]) / ε, (S.R0[r] - R2[r]) / ε, (R1[r] - R2[r]) / 2ε, sa[r])
            end
        end
        if haskey(ENV, "NKPW")                    # NKPW=n: windows of n zones, the check vector restricted to each (all variables), then per variable in the worst window of each row: where are the columns that disagree?
            nw = parse(Int, ENV["NKPW"]); zall = [zq(kb) for kb in 1:nb]; wins = [(z, min(z + nw - 1, maximum(zall))) for z in minimum(zall):nw:maximum(zall)]
            wres = Dict(r => Tuple{Float64,Int}[] for r in rows)
            for (iw, (w1, w2)) in enumerate(wins)
                vw = [S.act[r] && w1 <= zq((r - 1) ÷ ne + 1) <= w2 ? v[r] : 0.0 for r in 1:S.n]; any(!iszero, vw) || continue
                yA = zeros(S.n); adj!(yA, S, vw); yF = zeros(S.n); fdjv!(yF, S, vw, ε)
                for r in rows; push!(wres[r], (abs(yA[r] - yF[r]) / sa[r], iw)); end
            end
            for r in rows
                kb = (r - 1) ÷ ne + 1; i = (r - 1) % ne + 1; ws = sort!(wres[r]; rev = true)[1:min(4, length(wres[r]))]
                @printf("NKCHECK AD/FD windows: row %s zone %d (rel %.3e): worst windows %s\n", i <= length(nm) ? nm[i] : string(i), zq(kb), rr1[r], join((@sprintf("%d:%d %.2e", wins[iw][1], wins[iw][2], d) for (d, iw) in ws), ", "))
                w1, w2 = wins[ws[1][2]]
                for j in 1:ne
                    vj = [S.act[rr_] && w1 <= zq((rr_ - 1) ÷ ne + 1) <= w2 && (rr_ - 1) % ne + 1 == j ? v[rr_] : 0.0 for rr_ in 1:S.n]; any(!iszero, vj) || continue
                    yA = zeros(S.n); adj!(yA, S, vj); yF = zeros(S.n); fdjv!(yF, S, vj, ε)
                    d = abs(yA[r] - yF[r]) / sa[r]; d > 1e-6 && @printf("      variable %-5s in zones %d:%d: |AD-FD|/rowscale %.3e (AD %+.4e FD %+.4e)\n", j <= length(nm) ? nm[j] : string(j), w1, w2, d, yA[r], yF[r])
                end
            end
            flush(stdout); return nothing
        end
        cols = [c for c in findall(S.act) if get(ENV, "NKPZ", "34:48") == "all" || zp1 <= zq((c - 1) ÷ ne + 1) <= zp2]
        @printf("NKCHECK AD/FD column probe: %d rows, %d columns (zones %d:%d), FD at eps and eps/10\n", length(rows), length(cols), zp1, zp2); flush(stdout)
        tab = Dict(r => Tuple{Float64,Int,Float64,Float64,Float64}[] for r in rows)
        for c in cols
            e = zeros(S.n); e[c] = xcap(c, 1.0, OPT.eps); yA = zeros(S.n); adj!(yA, S, e); yF = zeros(S.n); fdjv!(yF, S, e, OPT.eps); yG = zeros(S.n); fdjv!(yG, S, e, OPT.eps / 10)
            for r in rows; push!(tab[r], (abs(yA[r] - yF[r]), c, yA[r], yF[r], yG[r])); end
        end
        for r in rows
            kb = (r - 1) ÷ ne + 1; i = (r - 1) % ne + 1
            @printf("NKCHECK AD/FD column probe: row %s zone %d (rowscale %.3e, rel %.3e)\n", i <= length(nm) ? nm[i] : string(i), zq(kb), sa[r], rr1[r])
            for (d, c, yA_, yF_, yG_) in sort!(tab[r]; rev = true)[1:min(10, length(tab[r]))]
                kc = (c - 1) ÷ ne + 1; jc = (c - 1) % ne + 1
                @printf("   column %-5s zone %3d: AD %+.5e  FD(eps) %+.5e  FD(eps/10) %+.5e  |AD-FD| %.3e\n", jc <= length(nm) ? nm[jc] : string(jc), zq(kc), yA_, yF_, yG_, d)
            end
        end
        flush(stdout)
    end
    for r in sortperm(rel; rev = true)[1:30]
        kb = (r - 1) ÷ ne + 1; i = (r - 1) % ne + 1
        @printf("   row %-5s zone %3d: Jv %+.5e  J_b w %+.5e  rowscale %.3e  rel %.3e\n", i <= length(nm) ? nm[i] : string(i), k0 + kb - 1 - a + 1, y[r], z[r], sa[r], rel[r])
    end
    # per equation: the largest mismatch over all zones and over zones NKZ1:NKZ2 (default 28:34, the shock / H2 front / photosphere / ionization front of the 100-zone Pop III run)
    z1, z2 = parse.(Int, split(get(ENV, "NKZ", "28:34"), ":"))
    for i in 1:ne
        ri = [(kb - 1) * ne + i for kb in 1:nb]; zn(kb) = k0 + kb - 1 - a + 1
        kall = argmax([rel[r] for r in ri]); ks = [kb for kb in 1:nb if z1 <= zn(kb) <= z2]
        kz = isempty(ks) ? 0 : ks[argmax([rel[ri[kb]] for kb in ks])]
        @printf("NKCHECK eq %-5s: max rel %.3e (zone %d) | zones %d..%d: max rel %.3e (zone %d)\n", i <= length(nm) ? nm[i] : string(i), rel[ri[kall]], zn(kall), z1, z2, kz == 0 ? NaN : rel[ri[kz]], kz == 0 ? 0 : zn(kz))
    end
    for f in (1e-2, 1e-1, 10.0, 100.0)
        y2 = zeros(S.n); fdjv!(y2, S, v, ε * f); r2 = relof(y2); d2 = [ok[r] ? abs(y2[r] - y[r]) / sa[r] : 0.0 for r in 1:S.n]
        @printf("NKCHECK eps x%g (max|eps v| %.1e): vs J_b w: max %.3e median %.3e rows>1e-4 %d >1e-2 %d | vs the eps run: |Jv - Jv(eps)|/|Jv(eps)| (2-norm) %.3e, worst row (rel to rowscale) %.3e\n",
                f, OPT.eps * f, maximum(r2), median(r2[ok]), count(>(1e-4), r2), count(>(1e-2), r2), norm(y2 .- y) / norm(y), maximum(d2))
    end
    get(ENV, "NKROWPROBE", "1") == "0" && return nothing        # NKROWPROBE=0: skip the column-by-column probe (one finite-difference product per active column)
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
    if haskey(ENV, "NKLINE")             # NKLINE=smax,n,z1:z2: continuity of the residual along ln ρ → ln ρ + s on the zones z1:z2 (zone counted from the inner boundary), s ∈ [−smax, smax]
        f = split(ENV["NKLINE"], ","); smax = parse(Float64, f[1]); nl = parse(Int, f[2]); zz = parse.(Int, split(f[3], ":"))
        vl = zeros(S.n); for kb in clamp(zz[1] + a - k0, 1, nb):clamp(zz[2] + a - k0, 1, nb); S.act[(kb-1)*ne+3] && (vl[(kb-1)*ne+3] = 1.0); end
        scan(n) = [perturbed(S, vl, smax * (2i / n - 1)) for i in 0:n]
        Rc = scan(nl); Rf = scan(2nl); sc = [maximum(abs(R[r]) for R in Rf) + 1e-300 for r in 1:S.n]
        jump(Rs) = (m = 0.0; at = (0, 0); for i in 1:length(Rs)-1, r in 1:S.n; d = abs(Rs[i+1][r] - Rs[i][r]) / sc[r]; d > m && (m = d; at = (r, i)); end; (m, at))
        (jc, ac), (jf, af) = jump(Rc), jump(Rf)
        row(r) = (i = (r - 1) % ne + 1; kb = (r - 1) ÷ ne + 1; "$(i <= length(nm) ? nm[i] : i) zone $(k0 + kb - 1 - a + 1)")
        @printf("NKCHECK line scan s in [-%g, %g] (ln rho on zones %d:%d): largest row jump between neighbours, n %d: %.3e (%s, step %d); n %d: %.3e (%s, step %d); ratio %.3f (2 = Lipschitz, 1 = jump)\n",
                smax, smax, zz[1], zz[2], nl, jc, row(ac[1]), ac[2], 2nl, jf, row(af[1]), af[2], jc / jf)
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
    TitanCore.GLOBAL_NEWTON[] = activate
    activate && (TD.KRYLOV[] = krylov!)
    nothing
end
"Counters since the last call (`steps`: accepted steps so far), for the DTERR trace line."
function trace_text(steps)
    c = CNT; m = MARK[]; ds = c.solves - m[1]; dm = c.mv - m[2]; dn = c.nonconv - m[3]; nsteps = steps - m[4]
    MARK[] = (c.solves, c.mv, c.nonconv, steps)
    @sprintf(" | NK %.2f Krylov its/Newton it, %.1f mv/step, nonconv %d", ds == 0 ? 0.0 : dm / ds, nsteps == 0 ? 0.0 : dm / nsteps, dn)
end
glob_text() = @sprintf(" | NKGLOB: %d globalised steps, %d Krylov-truncated, %d trust-region capped, %d backtracked, %d failed, %d merit evaluations (%d trial states rejected)", GC.solves, GC.truncated, GC.capped, GC.backtracked, GC.fail, GC.evals, GC.rejected)
summary_text() = @sprintf("NKRYLOV (eta %g, mmax %d, eps %g): %d Newton linear solves (%d J_b reused), %d matrix-vector products (%.2f per solve), %d GMRES non-convergences, %d state-drift events (max %.2e), %.1f s",
                          OPT.eta, OPT.mmax, OPT.eps, CNT.solves, CNT.jb_reused, CNT.mv, CNT.solves == 0 ? 0.0 : CNT.mv / CNT.solves, CNT.nonconv, CNT.ndrift, CNT.maxdrift, CNT.secs) * glob_text()

end # module
