"""
    BlockSolver

Replacement for TITAN's matslv: solves the block-pentadiagonal Newton system
(em2..ep2, rhs in COMMON; solution overwrites rhs, as matslv does) by pairing zones so
that the system becomes block-TRIdiagonal with (2·neqn)×(2·neqn) blocks, then block
Thomas elimination with partial pivoting inside each diagonal (Schur-complement) block.
Rows are equilibrated (scaled by 1/max|row|) first, as in JuliaSolver.matslv!.
No pivoting across blocks (TITAN's own matslv does none either); all buffers are
preallocated and reused; zero allocations in steady state.
"""
module BlockSolver

using StaticArrays
import ..TitanDriver
const TD = TitanDriver
const TF = TD.TF

mutable struct Buf
    m::Int; np::Int
    A::Array{Float64,3}; B::Array{Float64,3}; C::Array{Float64,3}   # diag, super, sub blocks
    b::Matrix{Float64}; sc::Vector{Float64}; piv::Matrix{Int}; w::Matrix{Float64}
end
const BUF = Ref{Any}(nothing)
function getbuf(m, np)
    B = BUF[]
    if B === nothing || B.m != m || B.np != np
        B = Buf(m, np, zeros(m, m, np), zeros(m, m, np), zeros(m, m, np), zeros(m, np), zeros(m), zeros(Int, m, np), zeros(m, m))
        BUF[] = B
    end
    B::Buf
end

"in-place LU with partial pivoting of A[:,:,p] (m×m); pivots in piv[:,p]. Returns min |pivot|."
@inline function lu!(A, piv, p, m)
    mn = Inf
    @inbounds for j in 1:m
        # pivot search
        ip = j; amax = abs(A[j, j, p])
        for i in j+1:m
            a = abs(A[i, j, p]); a > amax && (amax = a; ip = i)
        end
        piv[j, p] = ip
        if ip != j
            for c in 1:m; A[j, c, p], A[ip, c, p] = A[ip, c, p], A[j, c, p]; end
        end
        d = A[j, j, p]; mn = min(mn, abs(d))
        id = inv(d)
        for i in j+1:m; A[i, j, p] *= id; end
        for c in j+1:m
            ajc = A[j, c, p]
            ajc == 0 && continue
            @simd for i in j+1:m; A[i, c, p] -= A[i, j, p] * ajc; end
        end
    end
    mn
end
"solve with the LU in A[:,:,p] for the m×ncol right-hand side X[:, 1:ncol, q] (in place)."
@inline function lus!(X, q, ncol, A, piv, p, m)
    @inbounds for c in 1:ncol
        for j in 1:m
            ip = piv[j, p]; ip != j && ((X[j, c, q], X[ip, c, q]) = (X[ip, c, q], X[j, c, q]))
        end
        for j in 1:m
            xj = X[j, c, q]; xj == 0 && continue
            @simd for i in j+1:m; X[i, c, q] -= A[i, j, p] * xj; end
        end
        for j in m:-1:1
            X[j, c, q] /= A[j, j, p]; xj = X[j, c, q]
            @simd for i in 1:j-1; X[i, c, q] -= A[i, j, p] * xj; end
        end
    end
    nothing
end
@inline function lusv!(y, q, A, piv, p, m)
    @inbounds begin
        for j in 1:m
            ip = piv[j, p]; ip != j && ((y[j, q], y[ip, q]) = (y[ip, q], y[j, q]))
        end
        for j in 1:m
            yj = y[j, q]; @simd for i in j+1:m; y[i, q] -= A[i, j, p] * yj; end
        end
        for j in m:-1:1
            y[j, q] /= A[j, j, p]; yj = y[j, q]
            @simd for i in 1:j-1; y[i, q] -= A[i, j, p] * yj; end
        end
    end
    nothing
end

const MINPIV = Ref(Inf)

function matslv!()
    neqn = Int(TF.getscalar(:neqn)); ngrs = Int(TF.getscalar(:ngrs)); ngre = Int(TF.getscalar(:ngre))
    k0 = ngrs - 2; nb = ngre - ngrs + 5; m = 2neqn; np = (nb + 1) ÷ 2
    E = (TF.var(:em2), TF.var(:em1), TF.var(:e00), TF.var(:ep1), TF.var(:ep2)); rhs = TF.var(:rhs)
    S = getbuf(m, np)
    _solve!(S, E, rhs, neqn, k0, nb, m, np)
end

function _solve!(S::Buf, E, rhs, neqn, k0, nb, m, np)
    A, B, C, b, sc, piv = S.A, S.B, S.C, S.b, S.sc, S.piv
    fill!(A, 0.0); fill!(B, 0.0); fill!(C, 0.0)
    @inbounds for kb in 1:2np
        p = (kb + 1) >> 1; lr = (kb - 1) & 1        # pair, local zone 0/1
        if kb > nb                                    # padding zone: identity
            for i in 1:neqn; A[lr*neqn+i, lr*neqn+i, p] = 1.0; b[lr*neqn+i, p] = 0.0; end
            continue
        end
        k = k0 + kb - 1
        # row scale
        for i in 1:neqn
            s = 0.0
            for o in -2:2
                cb = kb + o; (1 <= cb <= nb) || continue
                M = E[o+3]
                for j in 1:neqn; s = max(s, abs(M[i, j, k])); end
            end
            sc[i] = s > 0 ? inv(s) : 1.0
            b[lr*neqn+i, p] = -rhs[i, k] * sc[i]
        end
        for o in -2:2
            cb = kb + o; (1 <= cb <= nb) || continue
            q = (cb + 1) >> 1; lc = (cb - 1) & 1
            T = q == p ? A : q == p + 1 ? B : C
            M = E[o+3]
            for j in 1:neqn, i in 1:neqn
                T[lr*neqn+i, lc*neqn+j, p] = M[i, j, k] * sc[i]
            end
        end
    end
    # forward elimination: D_p = A_p - C_p X_{p-1}, X_p = D_p^{-1} B_p (stored in B), y_p in b
    mn = Inf
    @inbounds for p in 1:np
        if p > 1
            # A[:,:,p] -= C[:,:,p] * B[:,:,p-1]; b[:,p] -= C[:,:,p] * b[:,p-1]
            for c in 1:m
                for l in 1:m
                    x = B[l, c, p-1]; x == 0 && continue
                    @simd for i in 1:m; A[i, c, p] -= C[i, l, p] * x; end
                end
            end
            for l in 1:m
                x = b[l, p-1]; x == 0 && continue
                @simd for i in 1:m; b[i, p] -= C[i, l, p] * x; end
            end
        end
        mn = min(mn, lu!(A, piv, p, m))
        p < np && lus!(B, p, m, A, piv, p, m)
        lusv!(b, p, A, piv, p, m)
    end
    MINPIV[] = mn
    # back substitution: x_p = y_p - X_p x_{p+1}
    @inbounds for p in np-1:-1:1
        for l in 1:m
            x = b[l, p+1]; x == 0 && continue
            @simd for i in 1:m; b[i, p] -= B[i, l, p] * x; end
        end
    end
    @inbounds for kb in 1:nb
        p = (kb + 1) >> 1; lr = (kb - 1) & 1
        for i in 1:neqn; rhs[i, k0+kb-1] = b[lr*neqn+i, p]; end
    end
    nothing
end



# ======================================================================
# Direct block-PENTAdiagonal elimination with neqn×neqn static blocks (the default).
module Penta
using StaticArrays, LinearAlgebra
import ...TitanDriver
const TF = TitanDriver.TF
const CACHE = Dict{Any,Any}()

function matslv!()
    neqn = Int(TF.getscalar(:neqn)); ngrs = Int(TF.getscalar(:ngrs)); ngre = Int(TF.getscalar(:ngre))
    k0 = ngrs - 2; nb = ngre - ngrs + 5
    E = (TF.var(:em2), TF.var(:em1), TF.var(:e00), TF.var(:ep1), TF.var(:ep2)); rhs = TF.var(:rhs)
    neqn == 7 ? _solve!(Val(7), E, rhs, k0, nb) :
    neqn == 5 ? _solve!(Val(5), E, rhs, k0, nb) : _solve!(Val(neqn), E, rhs, k0, nb)
end

@inline function loadblk(::Val{N}, M, k, sc) where {N}
    SMatrix{N,N,Float64}(ntuple(l -> (j = (l - 1) ÷ N + 1; i = l - (j - 1) * N; @inbounds M[i, j, k] * sc[i]), Val(N * N)))
end

function _solve!(::Val{N}, E, rhs, k0, nb) where {N}
    SM = SMatrix{N,N,Float64,N*N}; SV = SVector{N,Float64}
    em2, em1, e00, ep1, ep2 = E
    bufs = get!(CACHE, (N, nb)) do
        (Vector{SM}(undef, nb), Vector{SM}(undef, nb), Vector{SM}(undef, nb), Vector{SM}(undef, nb),
         Vector{SM}(undef, nb), Vector{SV}(undef, nb))
    end::Tuple{Vector{SM},Vector{SM},Vector{SM},Vector{SM},Vector{SM},Vector{SV}}
    L2, L1, D, U1, U2, y = bufs
    Z = zero(SM)
    @inbounds for kb in 1:nb
        k = k0 + kb - 1
        sc = SV(ntuple(i -> begin
            s = 0.0
            for o in -2:2
                (1 <= kb + o <= nb) || continue
                M = E[o+3]
                for j in 1:N; s = max(s, abs(M[i, j, k])); end
            end
            s > 0 ? inv(s) : 1.0
        end, Val(N)))
        L2[kb] = kb > 2 ? loadblk(Val(N), em2, k, sc) : Z
        L1[kb] = kb > 1 ? loadblk(Val(N), em1, k, sc) : Z
        D[kb] = loadblk(Val(N), e00, k, sc)
        U1[kb] = kb < nb ? loadblk(Val(N), ep1, k, sc) : Z
        U2[kb] = kb < nb - 1 ? loadblk(Val(N), ep2, k, sc) : Z
        y[kb] = SV(ntuple(i -> -rhs[i, k] * sc[i], Val(N)))
    end
    @inbounds for kb in 1:nb
        F = lu(D[kb])
        u1 = F \ U1[kb]; u2 = F \ U2[kb]; yy = F \ y[kb]
        U1[kb] = u1; U2[kb] = u2; y[kb] = yy
        if kb + 1 <= nb
            L = L1[kb+1]
            D[kb+1] -= L * u1; U1[kb+1] -= L * u2; y[kb+1] -= L * yy
        end
        if kb + 2 <= nb
            L = L2[kb+2]
            L1[kb+2] -= L * u1; D[kb+2] -= L * u2; y[kb+2] -= L * yy
        end
    end
    @inbounds for kb in nb:-1:1
        x = y[kb]
        kb + 1 <= nb && (x -= U1[kb] * y[kb+1])
        kb + 2 <= nb && (x -= U2[kb] * y[kb+2])
        y[kb] = x
        k = k0 + kb - 1
        for i in 1:N; rhs[i, k] = x[i]; end
    end
    nothing
end
end # Penta
end # module
