"""
    NewtonKrylov

NKRYLOV (popiii_m2c `NKRYLOV=1[,η[,mmax[,ε]]]`, provisional, default off = bit-identical; docs/POPIII_NEWTON_KRYLOV.md).

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
* GMRES without restart, modified Gram-Schmidt, Givens rotations; stops when the preconditioned residual is ≤ η ‖δ0‖ in the weighted 2-norm ‖W·‖₂, W the weights of the DTERR error norm (`weights`). If it does not get there in `mmax` iterations the last iterate (the best) is used and the event counted. Exact band Jacobian: δ0 is
  an eigenvector with eigenvalue 1 and the solve ends after one product.
Then the unchanged `maxdel!`, `update!`.
"""
module NewtonKrylov

using LinearAlgebra, Printf, Statistics
import ..TitanDriver, ..TitanCore, ..TitanResidFast, ..Trbdf2
const TD = TitanDriver; const TF = TD.TF; const TC = TitanCore; const TX = TitanResidFast
vf(s) = TF.var(s)::Vector{Float64}
gi(s) = Int(TD.g(s))

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
struct Sys
    ne::Int; nb::Int; k0::Int; n::Int
    code::Vector{Int}; act::BitVector            # act[(kb-1)ne+i]: the entry update! applies
    wt::Vector{Float64}                         # weight of each entry in the GMRES norm (Trbdf2.wnorm's)
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
function perturbed(S::Sys, v, s; realized = nothing)
    ne, nb, k0 = S.ne, S.nb, S.k0; rhs = TF.var(:rhs)
    base = realized === nothing ? nothing : (c = TX.context(); (copy(c.base), copy(vf(:unom)), copy(vf(:frnom))))
    @inbounds for kb in 1:nb, i in 1:ne; rhs[i, k0+kb-1] = S.act[(kb-1)*ne+i] ? s * v[(kb-1)*ne+i] : 0.0; end
    TX.FREEZE[] = 2
    try
        TD.r(:update)
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
    finally
        restore!(S.snap); TX.FREEZE[] = 0
    end
end
"""
Directional derivative of the full residual, y = J v: central difference [R(x ⊕ εv) − R(x ⊖ εv)]/2ε (one-sided [R(x ⊕ εv) − R(x)]/ε if `OPT.central` is false). v is a vector over the active entries.
"""
function fdjv!(y, S::Sys, v, ε; realized = nothing)
    R1 = perturbed(S, v, ε; realized)
    if OPT.central
        R2 = perturbed(S, v, -ε)
        @inbounds for r in 1:S.n; y[r] = (R1[r] - R2[r]) / (2ε); end
    else
        @inbounds for r in 1:S.n; y[r] = (R1[r] - S.R0[r]) / ε; end
    end
    y
end

function setup(; need_lu = true)
    TF.BACKEND[] === :julia || error("NKRYLOV needs the Julia state backend")
    c = TX.context(); ne, nb, k0 = c.neqn, c.nb, c.k0
    act = activeset(c); snap = Snap(); take!(snap)
    TX.FREEZE[] = 1
    R0 = try resid(ne, nb) finally TX.FREEZE[] = 0 end
    AB, ipiv, sc, kl, ku = bandlu(ne, nb, k0)
    Sys(ne, nb, k0, ne * nb, copy(c.code), act, weights(c, act), R0, snap, AB, ipiv, sc, kl, ku)
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
        all(isfinite, y) || error("NKRYLOV: non-finite J v")
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
        blk .= δ0           # (the perturbed updates used `rhs` as scratch); entries update! ignores keep the band solution
        @inbounds for r in 1:S.n; S.act[r] && (blk[(r-1)%ne+1, (r-1)÷ne+1] = x[r]); end
        Rf = resid(ne, nb); d = maximum(abs, Rf .- S.R0) / (maximum(abs, S.R0) + 1e-300)       # drift of the restored state
        CNT.maxdrift = max(CNT.maxdrift, d); d > 1e-10 && (CNT.ndrift += 1)
        Trbdf2.NWTRACE[] === nothing || println(Trbdf2.NWTRACE[], "NKR att ", Trbdf2.NWATT[], " it ", gi(:iter), " krylov its ", m, " converged ", conv, " rel residual ", rel, " |delta0|_W ", norm(b .* S.wt), " |delta|_W ", norm(x .* S.wt), " drift ", d, " largest weighted delta0 ", wtop(b, S, ne, k0), " largest weighted delta ", wtop(x, S, ne, k0))
    catch ex
        ex isa InterruptException && rethrow()
        S === nothing || restore!(S.snap)
        TX.FREEZE[] = 0; blk .= δ0; CNT.fallback += 1
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
    rstar = argmax(rel); kbs = (rstar - 1) ÷ ne + 1; is = (rstar - 1) % ne + 1
    @printf("NKCHECK row probe: finite-difference row %s of zone %d against the band row, one unit vector per active column (%d products)\n", is <= length(nm) ? nm[is] : string(is), k0 + kbs - 1 - a + 1, count(S.act))
    E = (vf2(:em2), vf2(:em1), vf2(:e00), vf2(:ep1), vf2(:ep2)); cols = Tuple{Float64,Int,Float64,Float64}[]
    for c in findall(S.act)
        e = zeros(S.n); e[c] = 1.0; yy = zeros(S.n); fdjv!(yy, S, e, OPT.eps)
        kc = (c - 1) ÷ ne + 1; jc = (c - 1) % ne + 1; o = kc - kbs
        jb = abs(o) <= 2 ? E[o+3][is, jc, k0 + kbs - 1] : 0.0         # (the band row: offset kb - kc, column jc)
        push!(cols, (abs(yy[rstar] - jb), c, yy[rstar], jb))
    end
    sort!(cols; rev = true)
    for (d, c, fd, jb) in cols[1:min(12, length(cols))]
        kc = (c - 1) ÷ ne + 1; jc = (c - 1) % ne + 1
        @printf("   column %-5s zone %3d (offset %+d): FD %+.4e  band %+.4e  |diff| %.3e\n", jc <= length(nm) ? nm[jc] : string(jc), k0 + kc - 1 - a + 1, kc - kbs, fd, jb, d)
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
    activate && (TD.KRYLOV[] = krylov!)
    nothing
end
"Counters since the last call (`steps`: accepted steps so far), for the DTERR trace line."
function trace_text(steps)
    c = CNT; m = MARK[]; ds = c.solves - m[1]; dm = c.mv - m[2]; dn = c.nonconv - m[3]; df = c.fallback - m[4]; nsteps = steps - m[5]
    MARK[] = (c.solves, c.mv, c.nonconv, c.fallback, steps)
    @sprintf(" | NK %.2f Krylov its/Newton it, %.1f mv/step, nonconv %d, fallback %d", ds == 0 ? 0.0 : dm / ds, nsteps == 0 ? 0.0 : dm / nsteps, dn, df)
end
summary_text() = @sprintf("NKRYLOV (eta %g, mmax %d, eps %g): %d Newton linear solves, %d matrix-vector products (%.2f per solve), %d GMRES non-convergences, %d fallbacks, %d state-drift events (max %.2e), %.1f s",
                          OPT.eta, OPT.mmax, OPT.eps, CNT.solves, CNT.mv, CNT.solves == 0 ? 0.0 : CNT.mv / CNT.solves, CNT.nonconv, CNT.fallback, CNT.ndrift, CNT.maxdrift, CNT.secs)

end # module
