"""
    H2LinesKA

KernelAbstractions port of the formal solution `H2Lines.transfer_line_fast!` (titan.jl, julia/H2Lines.jl, LINEFAST=1,
LINESRC=zone): tangent rays, zone-constant source, Gaussian thermal profile, comoving Doppler shifts. A batch of lines is solved in
three kernel launches:

  1. `sweep_kernel!`, one work item per (line, ray, chunk of the frequency bins): the inward and outward sweeps along the ray with
     the loop over the chunk's bins inside (`nchunk` = 1: one item per (line, ray)). Per item it writes, per ray node, the sums over the bins that the angle quadratures need
     (profile-weighted I and ψ for the two zones next to the node, and the net exchange of the ray with each segment).
  2. `csum_kernel!`, one work item per (line, ray, node): the chunks' partial sums added up.
  3. `quad_kernel!`, one work item per (line, zone): the angle quadrature of J̄ and Λ.

The flux quadrature (and everything that involves ν₀ dv/c, the source scale and r) is done on the host in Float64.

Float type: everything that runs on the device is generic in `T = eltype` of the device arrays; with `T = Float32` no Float64
value is formed on the device (literals are written as `T(…)`/`Float32`, see `psi`).

What is passed to the device is normalised on the host in Float64:
  - sources: `(S − B_bg)/scale`, `scale` a power of two near the largest |S − B_bg| of the line's zones that have opacity
    (so Planck-like values far below 1e-38 do not underflow as long as they are within 1e-38 of the line's largest);
  - lengths in units of a power of two near r[end]; the line-centre opacity as κ₀/(√π v_th ν₀/c) per unit length;
  - ν₀ dv/c and 4π-type factors never reach the device.

`form` (the flux, i.e. the bin-summed I_out − I_in at every ray node):
  - `form = 1` (default): the device accumulates, per ray segment, δ = Σ_bins (S − I) ψ for the inward and the outward passage
    (what the segment gave the ray), and the host forms ΔI(node) by a running sum in Float64. No difference of nearly equal
    numbers is taken on the device: (S − I) is exact in the float type when I is within a factor 2 of S.
  - `form = 0`: the device stores Σ_bins I_in per node and returns Σ_bins I_out − Σ_bins I_in (a difference of two sums; in
    Float32 this cancels where the ray is in equilibrium with the source). Kept to measure the cancellation.
"""
module H2LinesKA

using KernelAbstractions
const KA = KernelAbstractions

export Workspace, formal!, formal_batch, free!, setgeom!, cached_workspace, raygeom

const C_LIGHT = 2.99792458e10          # as H2Lines.c

# ------------------------------------------------------------------------------------------------------------------------
# device code (generic in T)
# ------------------------------------------------------------------------------------------------------------------------

"ψ = 1 − e^{−x} as the source computes it: x below 1e-4, −expm1(−x) above."
@inline psi(x::Float64) = ifelse(x < 1e-4, x, -expm1(-x))
# Float32: no expm1 (not available on every backend): x < 1e-4 → x (as the source); x < 0.5 → Taylor series to x¹⁰ (truncation
# 1e-11 relative); above 1 − e^{−x} (rounding ≤ 2.6 ulp).
@inline function psi(x::Float32)
    p = x * (1f0 + x * (-0.5f0 + x * (0.16666667f0 + x * (-0.041666668f0 + x * (0.008333334f0 + x * (-0.0013888889f0 +
        x * (0.0001984127f0 + x * (-2.4801588f-5 + x * (2.7557319f-6 + x * (-2.755732f-7))))))))))
    ifelse(x < 1f-4, x, ifelse(x < 0.5f0, p, 1f0 - exp(-x)))
end
"the weight of S in I ← I e + S (1 − e): Float64 as the source (1 − e); Float32 the accurate ψ (1 − e loses 6e-4 at Δτ = 1e-4)"
@inline om_exp(e::Float64, p::Float64) = 1.0 - e
@inline om_exp(e::Float32, p::Float32) = p

# scratch arrays I, ψ: `Val(true)` = [item, bin] (item is the fast index), `Val(false)` = [bin, item]
@inline gS(A, it, ib, ::Val{true}) = @inbounds A[it, ib]
@inline gS(A, it, ib, ::Val{false}) = @inbounds A[ib, it]
@inline sS!(A, v, it, ib, ::Val{true}) = @inbounds A[it, ib] = v
@inline sS!(A, v, it, ib, ::Val{false}) = @inbounds A[ib, it] = v

# sqrt without the DomainError branch of Base.sqrt (arguments here are ≥ 38): the plain LLVM sqrt on every backend
@inline usqrt(x) = Base.sqrt_llvm(x)

"`binwin` of the source: the bins lo:hi of the grid x1 + (ib − 1)/rdv that can lie within c0 ± w (no exceptions, no Int conversion of NaN)"
@inline function binwin(c0::T, w::T, x1::T, rdv::T, nb::Int32) where {T}
    a = (c0 - w - x1) * rdv + one(T); b = (c0 + w - x1) * rdv + one(T)
    nbf = T(nb)
    ok = (a == a) & (b == b)
    a = ifelse(ok, a, one(T)); b = ifelse(ok, b, nbf)
    lo = unsafe_trunc(Int32, clamp(floor(a), one(T), nbf + one(T)))
    hi = unsafe_trunc(Int32, clamp(ceil(b), zero(T), nbf))
    (lo, hi)
end

"the first bin ≥ lo of chunk c0 (0-based) of C (a power of two): the chunk's bins are those with (ib − 1) mod C == c0"
@inline first_bin(lo::Int32, c0::Int32, C::Int32) = lo + ((c0 - (lo - Int32(1))) & (C - Int32(1)))

"`fast_step!` of the source for the bins of item `it` that belong to the chunk; also returns δ = Σ_bins (S − I_before) ψ."
@inline function step!(I, ψ, xv, it, l, pr::T, qs::T, vt::T, S::T, x1::T, rdv::T, nb::Int32, lay, c0::Int32, C::Int32) where {T}
    pr == zero(T) && return (Int32(1), Int32(0), zero(T))
    lo, hi = binwin(-qs, vt * usqrt(T(38) + log(max(abs(pr), one(T)))), x1, rdv, nb)
    δ = zero(T); thin = T(1e-4)
    ib = first_bin(lo, c0, C)
    @inbounds while ib <= hi
        a = (xv[ib, l] + qs) / vt
        Δτ = pr * exp(-(a * a))
        e = exp(-Δτ); p = psi(Δτ)
        Io = gS(I, it, ib, lay)
        sS!(I, ifelse(Δτ < thin, Io * (one(T) - Δτ) + Δτ * S, Io * e + S * om_exp(e, p)), it, ib, lay)
        sS!(ψ, p, it, ib, lay)
        δ += (S - Io) * p
        ib += C
    end
    (lo, hi, δ)
end

"`fast_node` of the source: Σ_bins φ̂ I and Σ_bins φ̂ ψ (not normalised)"
@inline function node(I, ψ, xv, it, l, qs::T, vt::T, x1::T, rdv::T, nb::Int32, lay, c0::Int32, C::Int32) where {T}
    lo, hi = binwin(-qs, vt * usqrt(T(38)), x1, rdv, nb)
    sJ = zero(T); sΛ = zero(T)
    ib = first_bin(lo, c0, C)
    @inbounds while ib <= hi
        a = (xv[ib, l] + qs) / vt; p = exp(-(a * a))
        sJ += p * gS(I, it, ib, lay); sΛ += p * gS(ψ, it, ib, lay)
        ib += C
    end
    (sJ, sΛ)
end

@inline function binsum(I, it, nb::Int32, lay, c0::Int32, C::Int32)
    T = eltype(I); s = zero(T); ib = c0 + Int32(1)
    @inbounds while ib <= nb; s += gS(I, it, ib, lay); ib += C; end
    s
end
@inline function binfill!(I, v, it, lo::Int32, hi::Int32, lay, c0::Int32, C::Int32)
    ib = first_bin(lo, c0, C)
    @inbounds while ib <= hi; sS!(I, v, it, ib, lay); ib += C; end
    nothing
end

"""
One (line l, ray j, chunk c of C): the two sweeps of `transfer_line_fast!` for the bins ib with (ib − 1) mod C == c − 1 (the bins
are independent along a ray; only the sums over the bins couple them, and these are added up afterwards). Outputs for item
it = (j − 1) nl + l, node m of the ray (interface k = j + m − 1), partial sums over the chunk's bins:
  Eo[c, m, it]    form 1: δ_in + δ_out of segment m (between nodes m and m + 1);  form 0: Σ_bins (I_out − I_in) at node m
  QJ[c, s, m, it] Σ_bins φ̂ (I_in + I_out) with the profile of zone k − 2 + s;  QL the same for ψ
  Jmp[c, it]      form 1: −Σ_bins I_out at the photosphere node (the star is opaque: I′ is reset there), else 0
"""
@inline function sweep_item!(Eo, QJ, QL, Jmp, I, ψ, xv, rdv, dz, μu, μku, vth, kap, sA, sB, ipar, n::Int32, nl::Int32, form::Int32, lay, C::Int32, c, l, j)
    T = eltype(I)
    @inbounds begin
        it = (j - 1) * nl + l; c0 = c - Int32(1)
        nb = ipar[1, l]; lam = ipar[2, l] != Int32(0); kq = ipar[3, l]
        nk = n + Int32(2) - j
        mq = (kq > 0 && j < kq) ? kq - j + Int32(1) : Int32(0)
        x1 = xv[1, l]; rd = rdv[l]
        binfill!(I, zero(T), it, Int32(1), nb, lay, c0, C); binfill!(ψ, zero(T), it, Int32(1), nb, lay, c0, C)
        Jmp[c, it] = zero(T)
        Eo[c, nk, it] = zero(T); QJ[c, 1, nk, it] = zero(T); QJ[c, 2, nk, it] = zero(T); QL[c, 1, nk, it] = zero(T); QL[c, 2, nk, it] = zero(T)
        m = nk - Int32(1)
        while m >= 1                                              # inwards
            kz = j + m - Int32(1)
            lo, hi, δ = step!(I, ψ, xv, it, l, kap[kz, l] * dz[m, j], μu[m, j], vth[kz], sA[kz, l], x1, rd, nb, lay, c0, C)
            Eo[c, m, it] = form == Int32(0) ? binsum(I, it, nb, lay, c0, C) : δ
            qn = μku[m, j]
            for s in Int32(1):Int32(2)
                zs = kz - Int32(2) + s
                sJ = zero(T); sΛ = zero(T)
                if zs >= 1 && zs <= n
                    sJ, sΛ = node(I, ψ, xv, it, l, qn, vth[zs], x1, rd, nb, lay, c0, C)
                    if m == 1                                     # the turning point: I_out = I_in, ψ_out = ψ_in
                        oJ, oΛ = node(I, ψ, xv, it, l, -qn, vth[zs], x1, rd, nb, lay, c0, C); sJ += oJ; sΛ += oΛ
                    end
                end
                QJ[c, s, m, it] = sJ; QL[c, s, m, it] = ifelse(lam, sΛ, zero(T))
            end
            binfill!(ψ, zero(T), it, lo, hi, lay, c0, C)
            m -= Int32(1)
        end
        form == Int32(0) && (Eo[c, 1, it] = zero(T))              # ΔI = 0 at the turning point
        for m in Int32(2):nk                                      # outwards
            kz = j + m - Int32(2); k = j + m - Int32(1)
            S = (mq > 0 && m > mq) ? sB[kz, l] : sA[kz, l]
            lo, hi, δ = step!(I, ψ, xv, it, l, kap[kz, l] * dz[m-1, j], -μu[m-1, j], vth[kz], S, x1, rd, nb, lay, c0, C)
            form == Int32(0) || (Eo[c, m-1, it] += δ)
            if m == mq                                            # the star is opaque: I′ = 0 from here on
                form == Int32(0) || (Jmp[c, it] = -binsum(I, it, nb, lay, c0, C))
                binfill!(I, zero(T), it, Int32(1), nb, lay, c0, C)
                binfill!(ψ, zero(T), it, lo, hi, lay, c0, C)
            end
            form == Int32(0) && (Eo[c, m, it] = binsum(I, it, nb, lay, c0, C) - Eo[c, m, it])
            qn = μku[m, j]
            for s in Int32(1):Int32(2)
                zs = k - Int32(2) + s
                if zs >= 1 && zs <= n
                    sJ, sΛ = node(I, ψ, xv, it, l, -qn, vth[zs], x1, rd, nb, lay, c0, C)
                    QJ[c, s, m, it] += sJ; QL[c, s, m, it] += ifelse(lam, sΛ, zero(T))
                end
            end
            binfill!(ψ, zero(T), it, lo, hi, lay, c0, C)
        end
    end
    nothing
end

"One (node m, line l, ray j): the chunks' partial sums added up (E1, J1 go to the host for the flux; Q1J, Q1L to `quad_kernel!`)."
@inline function csum_item!(E1, J1, Q1J, Q1L, Eo, Jmp, QJ, QL, n::Int32, nl::Int32, C::Int32, m, l, j)
    T = eltype(E1)
    @inbounds if m <= n + Int32(2) - j
        it = (j - 1) * nl + l
        e = zero(T); a1 = zero(T); a2 = zero(T); b1 = zero(T); b2 = zero(T)
        for c in Int32(1):C
            e += Eo[c, m, it]; a1 += QJ[c, 1, m, it]; a2 += QJ[c, 2, m, it]; b1 += QL[c, 1, m, it]; b2 += QL[c, 2, m, it]
        end
        E1[m, it] = e; Q1J[1, m, it] = a1; Q1J[2, m, it] = a2; Q1L[1, m, it] = b1; Q1L[2, m, it] = b2
        if m == 1
            g = zero(T)
            for c in Int32(1):C; g += Jmp[c, it]; end
            J1[it] = g
        end
    end
    nothing
end

"One (line l, zone z): J̄ and Λ of the zone, the μ quadratures of `transfer_line_fast!` at its two interfaces. wJ = ¼ dv/(√π v_th)."
@inline function quad_item!(Jz, Lz, QJ, QL, dμ, wJ, n::Int32, nl::Int32, l, z)
    T = eltype(Jz); half = T(0.5)
    @inbounds begin
        aJ = zero(T); aL = zero(T)
        for s in Int32(1):Int32(2)
            k = z + Int32(2) - s                                   # zone z is outside interface z (s = 2) and inside z + 1 (s = 1)
            if k == 1
                aJ += QJ[s, 1, l]; aL += QL[s, 1, l]
            else
                cJ = zero(T); cL = zero(T)
                for j in Int32(1):(k-Int32(1))
                    d = dμ[k, j]
                    if d > zero(T)
                        itb = (j - 1) * nl + l; mb = k - j + 1
                        qa = zero(T); la = zero(T)
                        if j + 1 <= n
                            qa = QJ[s, mb-1, itb+nl]; la = QL[s, mb-1, itb+nl]
                        end
                        cJ += half * (qa + QJ[s, mb, itb]) * d; cL += half * (la + QL[s, mb, itb]) * d
                    end
                end
                aJ += cJ; aL += cL
            end
        end
        Jz[z, l] = wJ[z, l] * aJ; Lz[z, l] = wJ[z, l] * aL
    end
    nothing
end

@kernel function sweep_kernel!(Eo, QJ, QL, Jmp, I, ψ, @Const(xv), @Const(rdv), @Const(dz), @Const(μu), @Const(μku), @Const(vth),
                               @Const(kap), @Const(sA), @Const(sB), @Const(ipar), n, nl, form, lay, C)
    c, l, j = @index(Global, NTuple)
    sweep_item!(Eo, QJ, QL, Jmp, I, ψ, xv, rdv, dz, μu, μku, vth, kap, sA, sB, ipar, n, nl, form, lay, C, c % Int32, l % Int32, j % Int32)
end

@kernel function csum_kernel!(E1, J1, Q1J, Q1L, @Const(Eo), @Const(Jmp), @Const(QJ), @Const(QL), n, nl, C)
    m, l, j = @index(Global, NTuple)
    csum_item!(E1, J1, Q1J, Q1L, Eo, Jmp, QJ, QL, n, nl, C, m % Int32, l % Int32, j % Int32)
end

@kernel function quad_kernel!(Jz, Lz, @Const(QJ), @Const(QL), @Const(dμ), @Const(wJ), n, nl)
    l, z = @index(Global, NTuple)
    quad_item!(Jz, Lz, QJ, QL, dμ, wJ, n, nl, l % Int32, z % Int32)
end

# ------------------------------------------------------------------------------------------------------------------------
# host side
# ------------------------------------------------------------------------------------------------------------------------

"Ray geometry, a copy of `H2Lines.raygeom` (any object with the fields dz, μu, μku, μa, μb of size (n+1, n) can be passed instead)."
function raygeom(r, u)
    n = length(r) - 1
    dz = zeros(n + 1, n); μu = zeros(n + 1, n); μku = zeros(n + 1, n); μa = zeros(n + 1, n); μb = zeros(n + 1, n)
    zk = zeros(n + 1)
    for j in 1:n
        p = r[j]; nk = n + 2 - j
        for m in 1:nk; zk[m] = sqrt(max(r[j+m-1]^2 - p^2, 0.0)); end
        for m in 1:nk-1
            kz = j + m - 1
            zm = 0.5 * (zk[m] + zk[m+1]); rm = sqrt(p^2 + zm^2)
            w = clamp((rm - r[kz]) / (r[kz+1] - r[kz]), 0.0, 1.0)
            ur = (1 - w) * u[kz] + w * u[kz+1]
            μ = rm > 0 ? zm / rm : 0.0
            dz[m, j] = zk[m+1] - zk[m]; μu[m, j] = μ * ur
        end
        for m in 1:nk
            k = j + m - 1; μk = r[k] > 0 ? zk[m] / r[k] : 0.0
            μku[m, j] = μk * u[k]
        end
    end
    for k in 2:n+1, j in 1:k-1
        μa[k, j] = sqrt(max(1 - (r[j+1] / r[k])^2, 0.0)); μb[k, j] = sqrt(max(1 - (r[j] / r[k])^2, 0.0))
    end
    (; dz, μu, μku, μa, μb)
end

"""
    Workspace(backend, T, r, u, vth, G; nlmax, nbmax = 200, nchunk = default, lay = default)

Device arrays for batches of up to `nlmax` lines with up to `nbmax` bins on the grid `r` (n + 1 interfaces), `u`, `vth` (n zones),
`G` the ray geometry (`H2Lines.raygeom(r, u)` or `raygeom(r, u)` of this module). The geometry is uploaded here, once.
`nchunk` (a power of two; default 1 on the CPU, 8 on a GPU): each (line, ray) is solved by `nchunk` work items that share its
frequency bins (bin ib belongs to chunk (ib − 1) mod nchunk); a GPU thread is slow in serial work and the sweep is a long serial
loop, so splitting the bins shortens the kernel until the device is saturated. The per-node sums are stored per chunk.
`lay`: layout of the per-item bin scratch, `true` = [item, bin], `false` = [bin, item] (default; on the A6000 the faster one for
nchunk ≥ 4, the slower one for nchunk = 1).
"""
mutable struct Workspace{T,B,D,H}
    backend::B
    n::Int; nlmax::Int; nbmax::Int; nchunk::Int; lay::Bool
    Lu::Float64                      # length unit
    vth::Vector{Float64}; μa::Matrix{Float64}; μb::Matrix{Float64}; w2::Matrix{Float64}
    d::D                             # device arrays
    h::H                             # host staging arrays (eltype T)
    scale::Vector{Float64}; ffac::Vector{Float64}; mq::Vector{Int}
    timings::Dict{Symbol,Float64}
end

function Workspace(backend, ::Type{T}, r, u, vth, G; nlmax::Integer, nbmax::Integer = 200, nchunk::Integer = backend isa KA.CPU ? 1 : 8,
                   lay::Bool = false) where {T}
    n = length(vth); @assert length(r) == n + 1 && size(G.dz) == (n + 1, n)
    ispow2(nchunk) || error("H2LinesKA: nchunk must be a power of two"); C = Int(nchunk)
    nit = nlmax * n
    al(dims...) = KA.zeros(backend, T, dims...)
    d = (dz = al(n + 1, n), μu = al(n + 1, n), μku = al(n + 1, n), dμ = al(n + 1, n), vth = al(n),
         xv = al(nbmax, nlmax), rdv = al(nlmax), kap = al(n, nlmax), sA = al(n, nlmax), sB = al(n, nlmax), wJ = al(n, nlmax),
         ipar = KA.zeros(backend, Int32, 3, nlmax),
         I = lay ? al(nit, nbmax) : al(nbmax, nit), ψ = lay ? al(nit, nbmax) : al(nbmax, nit),
         Eo = al(C, n + 1, nit), QJ = al(C, 2, n + 1, nit), QL = al(C, 2, n + 1, nit), Jmp = al(C, nit), E1 = al(n + 1, nit), J1 = al(nit),
         Q1J = al(2, n + 1, nit), Q1L = al(2, n + 1, nit),
         Jz = al(n, nlmax), Lz = al(n, nlmax))
    h = (xv = zeros(T, nbmax, nlmax), rdv = zeros(T, nlmax), kap = zeros(T, n, nlmax), sA = zeros(T, n, nlmax), sB = zeros(T, n, nlmax),
         wJ = zeros(T, n, nlmax), ipar = zeros(Int32, 3, nlmax))
    ws = Workspace{T,typeof(backend),typeof(d),typeof(h)}(backend, n, nlmax, nbmax, C, lay, 1.0, zeros(n), zeros(n + 1, n), zeros(n + 1, n), zeros(n + 1, n),
                                                          d, h, zeros(nlmax), zeros(nlmax), zeros(Int, nlmax), Dict{Symbol,Float64}())
    setgeom!(ws, r, u, vth, G)
end

"""
    setgeom!(ws, r, u, vth, G)

A new grid state (same number of zones) into an existing workspace: uploads the ray geometry `G = raygeom(r, u)` and `vth`.
"""
function setgeom!(ws::Workspace{T}, r, u, vth, G) where {T}
    n = ws.n; d = ws.d
    (length(vth) == n && length(r) == n + 1 && size(G.dz) == (n + 1, n)) || error("H2LinesKA: grid size differs from the workspace's")
    Lu = 2.0^exponent(Float64(r[end])); ws.Lu = Lu
    ws.vth .= vth; ws.μa .= G.μa; ws.μb .= G.μb
    dμ = zeros(n + 1, n); fill!(ws.w2, 0.0)
    for k in 2:n+1, j in 1:k-1
        dμ[k, j] = G.μb[k, j] - G.μa[k, j]
        dμ[k, j] > 0 && (ws.w2[k, j] = 0.25 * (G.μb[k, j]^2 - G.μa[k, j]^2))   # LINEQUAD=p2: ¼ (μb² − μa²), 0 where the source skips the term
    end
    copyto!(d.dz, T.(G.dz ./ Lu)); copyto!(d.μu, T.(G.μu)); copyto!(d.μku, T.(G.μku)); copyto!(d.dμ, T.(dμ)); copyto!(d.vth, T.(vth))
    ws
end

const WSCACHE = Ref{Any}(nothing)
"free the workspace kept by `cached_workspace`"
clear_cache!() = (WSCACHE[] === nothing || free!(WSCACHE[]); WSCACHE[] = nothing)
"""
    cached_workspace(backend, T, r, u, vth, G; nlmax, nbmax = 200, kw...)

The workspace of the previous call if it fits (same backend, float type, zones, `nlmax`, `nbmax`), with the new geometry
uploaded; otherwise the old one is freed and a new one built. For a caller that solves one grid state after another
(not thread-safe: one caller at a time).
"""
function cached_workspace(backend, ::Type{T}, r, u, vth, G; nlmax::Integer, nbmax::Integer = 200, kw...) where {T}
    ws = WSCACHE[]
    if ws isa Workspace{T} && ws.backend === backend && ws.n == length(vth) && ws.nlmax == nlmax && ws.nbmax == nbmax
        return setgeom!(ws, r, u, vth, G)
    end
    ws === nothing || free!(ws)
    WSCACHE[] = Workspace(backend, T, r, u, vth, G; nlmax, nbmax, kw...)
end

"give the workspace's device arrays back right away (instead of waiting for the garbage collector); the workspace is dead afterwards"
free!(ws::Workspace) = (foreach(a -> KA.unsafe_free!(a), values(ws.d)); nothing)

"bytes held on the device by a workspace"
devbytes(ws::Workspace) = sum(sizeof(eltype(a)) * length(a) for a in values(ws.d))

# one line into the staging arrays (Float64 → T): returns nothing
function pack_line!(ws::Workspace{T}, l, κ0, ν0, xs, lam::Bool, Sz, bg) where {T}
    n = ws.n; h = ws.h; vth = ws.vth
    nb = length(xs); nb <= ws.nbmax || error("H2LinesKA: line with $nb bins, workspace has nbmax = $(ws.nbmax)")
    dv = Float64(step(xs))
    hasbg = bg !== nothing; Bc = hasbg ? Float64(bg.Bc) : 0.0; Bs = hasbg ? Float64(bg.Bs) : 0.0; kq = hasbg ? Int(bg.kq) : 0
    smax = 0.0
    @inbounds for k in 1:n
        κ0[k] == 0 && continue
        smax = max(smax, abs(Sz[k] - Bc), abs(Sz[k] - Bs))
    end
    scale = (smax > 0 && isfinite(smax)) ? 2.0^exponent(smax) : 1.0
    @inbounds for k in 1:n
        if κ0[k] == 0
            h.kap[k, l] = zero(T); h.sA[k, l] = zero(T); h.sB[k, l] = zero(T)
        else
            h.kap[k, l] = T(κ0[k] / (sqrt(π) * vth[k] / C_LIGHT * ν0) * ws.Lu)
            h.sA[k, l] = T((Sz[k] - Bc) / scale); h.sB[k, l] = T((Sz[k] - Bs) / scale)
        end
        h.wJ[k, l] = T(0.25 * dv / (sqrt(π) * vth[k]))
    end
    @inbounds for ib in 1:nb; h.xv[ib, l] = T(xs[ib]); end
    h.rdv[l] = T(1 / dv)
    h.ipar[1, l] = nb; h.ipar[2, l] = lam ? 1 : 0; h.ipar[3, l] = kq
    ws.scale[l] = scale; ws.ffac[l] = (dv / C_LIGHT * ν0) * scale; ws.mq[l] = kq
    nothing
end

"the flux quadrature of `transfer_line_fast!` for one line from the device's per-node sums (Float64)"
function flux_line!(Fν, ΔI, Eo, Jmp, l, nl, n, kq, form, quadp2, Ma, Mb, w2, ffac)
    fill!(ΔI, 0.0)
    @inbounds for j in 1:n
        it = (j - 1) * nl + l; nk = n + 2 - j
        if form == 0
            for m in 2:nk; ΔI[j, j+m-1] = Eo[m, it]; end
        else
            mq = (kq > 0 && j < kq) ? kq - j + 1 : 0
            acc = 0.0
            for m in 2:nk
                acc += Float64(Eo[m-1, it]); m == mq && (acc += Float64(Jmp[it]))
                ΔI[j, j+m-1] = acc
            end
        end
    end
    Fν[1] = 0.0                                               # ΔI[j, k]: ray j at interface k (row n + 1 stays 0: no ray n + 1)
    @inbounds for k in 2:n+1
        acc = 0.0
        if quadp2
            for j in 1:k-1; acc += (ΔI[j+1, k] + ΔI[j, k]) * w2[k, j]; end
        else
            for j in 1:k-1
                μa = Ma[k, j]; μb = Mb[k, j]
                Ia = ΔI[j+1, k]; Ib = ΔI[j, k]
                d = μb - μa; d <= 0 && continue
                s1 = (Ib - Ia) / d
                acc += Ia * (μb^2 - μa^2) / 2 + s1 * ((μb^3 - μa^3) / 3 - μa * (μb^2 - μa^2) / 2)
            end
        end
        Fν[k] = 2π * acc * ffac
    end
    nothing
end

"""
    formal!(ws, lines; quadp2 = true, form = 1) -> (Fν, Jz, Λz)

The formal solution of all `lines` (an indexable collection; each element has the fields `κ0` (n), `ν0`, `xs` (range of bin
velocities), `lam::Bool`, `Sz` (n), `bg` (`nothing` or `(Bc, Bs, kq)`): the per-line arguments of `transfer_line_fast!`).
Returns Float64 matrices: Fν (n+1, nl), Jz (n, nl), Λz (n, nl); column l = the tuple `transfer_line_fast!` returns for line l.
`quadp2`: the flux quadrature (`H2Lines.QUADP2[]`). `ws.timings` gets the wall time of the parts [s].
"""
function formal!(ws::Workspace{T}, lines; quadp2::Bool = true, form::Integer = 1) where {T}
    nl = length(lines); n = ws.n; d = ws.d; h = ws.h; be = ws.backend
    0 < nl <= ws.nlmax || error("H2LinesKA: $nl lines, workspace has nlmax = $(ws.nlmax)")
    t0 = time_ns()
    for (l, L) in enumerate(lines); pack_line!(ws, l, L.κ0, L.ν0, L.xs, L.lam, L.Sz, L.bg); end
    t1 = time_ns()
    copyto!(d.xv, h.xv); copyto!(d.rdv, h.rdv); copyto!(d.kap, h.kap); copyto!(d.sA, h.sA); copyto!(d.sB, h.sB)
    copyto!(d.wJ, h.wJ); copyto!(d.ipar, h.ipar)
    KA.synchronize(be); t2 = time_ns()
    C = Int32(ws.nchunk)
    sweep_kernel!(be)(d.Eo, d.QJ, d.QL, d.Jmp, d.I, d.ψ, d.xv, d.rdv, d.dz, d.μu, d.μku, d.vth, d.kap, d.sA, d.sB, d.ipar,
                      Int32(n), Int32(nl), Int32(form), Val(ws.lay), C; ndrange = (ws.nchunk, nl, n))
    KA.synchronize(be); t3 = time_ns()
    csum_kernel!(be)(d.E1, d.J1, d.Q1J, d.Q1L, d.Eo, d.Jmp, d.QJ, d.QL, Int32(n), Int32(nl), C; ndrange = (n + 1, nl, n))
    quad_kernel!(be)(d.Jz, d.Lz, d.Q1J, d.Q1L, d.dμ, d.wJ, Int32(n), Int32(nl); ndrange = (nl, n))
    KA.synchronize(be); t4 = time_ns()
    nit = nl * n
    Eo = Array(view(d.E1, :, 1:nit)); Jmp = Array(view(d.J1, 1:nit)); Jd = Array(view(d.Jz, :, 1:nl)); Ld = Array(view(d.Lz, :, 1:nl))
    t5 = time_ns()
    Fν = zeros(n + 1, nl); Jz = zeros(n, nl); Λz = zeros(n, nl)
    nt = min(Threads.nthreads(), nl)
    @sync for t in 1:nt
        Threads.@spawn begin
            ΔI = zeros(n + 1, n + 1)
            for l in t:nt:nl
                flux_line!(view(Fν, :, l), ΔI, Eo, Jmp, l, nl, n, ws.mq[l], Int(form), quadp2, ws.μa, ws.μb, ws.w2, ws.ffac[l])
                sc = ws.scale[l]
                @inbounds for k in 1:n; Jz[k, l] = Float64(Jd[k, l]) * sc; Λz[k, l] = Float64(Ld[k, l]); end
            end
        end
    end
    t6 = time_ns()
    tm = ws.timings
    tm[:pack] = (t1 - t0) / 1e9; tm[:upload] = (t2 - t1) / 1e9; tm[:sweep] = (t3 - t2) / 1e9; tm[:quad] = (t4 - t3) / 1e9
    tm[:download] = (t5 - t4) / 1e9; tm[:post] = (t6 - t5) / 1e9; tm[:total] = (t6 - t0) / 1e9
    (Fν, Jz, Λz)
end

"""
    formal_batch(backend, T, r, u, vth, G, lines; quadp2 = true, form = 1, nbmax = 200) -> (Fν, Jz, Λz)

One-shot entry point: builds a `Workspace` for `lines` and solves them (see `formal!`). In a loop over ALI iterations build the
workspace once and call `formal!`.
"""
function formal_batch(backend, ::Type{T}, r, u, vth, G, lines; quadp2::Bool = true, form::Integer = 1, nbmax::Integer = 200, kw...) where {T}
    ws = Workspace(backend, T, r, u, vth, G; nlmax = length(lines), nbmax = max(nbmax, maximum(length(L.xs) for L in lines)), kw...)
    res = formal!(ws, lines; quadp2, form)
    free!(ws)
    res
end

end # module
