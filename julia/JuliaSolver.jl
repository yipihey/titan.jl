"""
    JuliaSolver

Julia replacements for parts of the Newton step:

* `matslv!()`: solves TITAN's block-pentadiagonal Newton system with a
  *pivoted* banded LU (LAPACK gbtrf/gbtrs via LinearAlgebra), in place of the
  Fortran unpivoted block elimination.  Same interface: reads em2..ep2/rhs from
  COMMON, writes the correction into rhs(1:neqn, ngrs-2:ngre+2).
* `contin_rows!()`: overwrites the continuity-equation rows (id, k = ngrs..ngre)
  of the Fortran-assembled system with the Julia residual and its ForwardDiff
  Jacobian (columns r, ρ, u), i.e. a partial Julia assembly.
"""
module JuliaSolver

using LinearAlgebra, ForwardDiff
import ..TitanDriver
const TD = TitanDriver
const TF = TD.TF

const V = Dict{Symbol,Array}()
v(s::Symbol) = get!(() -> TF.var(s), V, s)
gi(s) = Int(TD.g(s))

# ------------------------------------------------------------- banded solver
const BUF = Dict{Any,Any}()
const A3 = Ref{NTuple{5,Array{Float64,3}}}()
const RHS = Ref{Matrix{Float64}}()
function bindarrays!()
    isassigned(A3) && return
    A3[] = (TF.var(:em2), TF.var(:em1), TF.var(:e00), TF.var(:ep1), TF.var(:ep2))
    RHS[] = TF.var(:rhs)
end

"Pivoted banded LU of TITAN's block-pentadiagonal system with row equilibration."
function matslv!()
    bindarrays!()
    neqn = gi(:neqn); ngrs = gi(:ngrs); ngre = gi(:ngre)
    k0 = ngrs - 2; nb = ngre - ngrs + 5; n = neqn * nb
    kl = ku = 3neqn - 1
    AB = get!(() -> zeros(2kl + ku + 1, n), BUF, (:AB, n))::Matrix{Float64}
    b  = get!(() -> zeros(n), BUF, (:b, n))::Vector{Float64}
    sc = get!(() -> zeros(n), BUF, (:sc, n))::Vector{Float64}
    _assemble!(AB, b, sc, A3[], RHS[], neqn, k0, nb, kl, ku)
    AB2, ipiv = LAPACK.gbtrf!(kl, ku, n, AB)
    LAPACK.gbtrs!('N', kl, ku, n, AB2, ipiv, reshape(b, n, 1))
    rhs = RHS[]
    @inbounds for kb in 1:nb, i in 1:neqn
        rhs[i, k0 + kb - 1] = b[(kb - 1) * neqn + i]
    end
    nothing
end

function _assemble!(AB, b, sc, E, rhs, neqn, k0, nb, kl, ku)
    fill!(AB, 0.0); fill!(sc, 0.0)
    @inbounds for kb in 1:nb, o in -2:2                 # row maxima for equilibration
        cb = kb + o; (1 <= cb <= nb) || continue
        k = k0 + kb - 1; M = E[o+3]
        for j in 1:neqn, i in 1:neqn
            r = (kb - 1) * neqn + i
            sc[r] = max(sc[r], abs(M[i, j, k]))
        end
    end
    @inbounds for r in eachindex(sc); sc[r] = sc[r] > 0 ? inv(sc[r]) : 1.0; end
    @inbounds for kb in 1:nb, o in -2:2
        cb = kb + o; (1 <= cb <= nb) || continue
        k = k0 + kb - 1; M = E[o+3]
        for j in 1:neqn, i in 1:neqn
            r = (kb - 1) * neqn + i; c = (cb - 1) * neqn + j
            AB[kl + ku + 1 + r - c, c] = M[i, j, k] * sc[r]
        end
    end
    @inbounds for kb in 1:nb, i in 1:neqn
        r = (kb - 1) * neqn + i
        b[r] = -rhs[i, k0 + kb - 1] * sc[r]
    end
    nothing
end

# ------------------------------------------- continuity rows in Julia (AD)
"Continuity residual (contin.f + advectc.f + diffuse.f, zet = 0) for zones ngrs..ngre."
function contin_res(x::AbstractVector{T}, c) where {T}
    (; ngrs, ngre, K, nK, rn0, dn0, un0, ro, dO, uo, dso, dvolo, unom, θ, dt, cadv, eps,
       sig, mu, phil0, phir0, leobc) = c
    mgr = length(rn0)
    rn = Vector{T}(rn0); dn = Vector{T}(dn0); un = Vector{T}(un0)
    @inbounds for (n, k) in enumerate(K)
        rn[k] = exp(x[n]); dn[k] = exp(x[nK+n]); un[k] = x[2nK+n] * unom[k]
    end
    r = θ .* rn .+ (1 - θ) .* ro
    u = θ .* un .+ (1 - θ) .* uo
    urel = u .- (rn .- ro) ./ dt
    rmu = r .^ mu
    R = zeros(T, ngre - ngrs + 1)
    dq = zeros(T, mgr); qsn = zeros(T, mgr); db = zeros(T, mgr); df = zeros(T, mgr)
    @inbounds for k in ngrs-1:ngre+1; dq[k] = dn[k+1] - dn[k]; end
    @inbounds for k in ngrs:ngre+1
        sdq = dq[k-1] + dq[k]
        den = sdq + (ForwardDiff.value(sdq) >= 0 ? eps : -eps)
        qsn[k] = ForwardDiff.value(dq[k-1] * dq[k]) >= 0 ? cadv * dq[k-1] * dq[k] / den : zero(T)
    end
    @inbounds for k in ngrs+1:ngre+1
        s = ForwardDiff.value(urel[k]) >= 0 ? 0.5 : -0.5
        db[k] = (0.5 + s) * (θ * (dn[k-1] + 0.5qsn[k-1]) + (1 - θ) * (dO[k-1] + 0.5dso[k-1])) +
                (0.5 - s) * (θ * (dn[k]   - 0.5qsn[k])   + (1 - θ) * (dO[k]   - 0.5dso[k]))
    end
    if sig != 0
        d = θ .* dn .+ (1 - θ) .* dO
        @inbounds for k in ngrs:ngre+1
            df[k] = 2sig * rmu[k] * (d[k] - d[k-1]) / (r[k+1] - r[k-1])
        end
    end
    @inbounds for k in ngrs:ngre
        dv = (rn[k+1]^(mu+1) - rn[k]^(mu+1)) / (mu + 1)
        fin  = k == ngrs ? rmu[k] * phil0 : rmu[k] * urel[k] * db[k]
        fout = k == ngre ? (leobc == 3 ? rmu[k+1] * phir0 + rmu[k+1] * urel[k+1] * db[k+1] :
                            rmu[k+1] * phir0) : rmu[k+1] * urel[k+1] * db[k+1]
        R[k-ngrs+1] = (dn[k] * dv - dO[k] * dvolo[k]) / dt + fout - fin - df[k+1] + df[k]
    end
    R
end

function contin_ctx()
    ngrs = gi(:ngrs); ngre = gi(:ngre); K = (ngrs - 2):(ngre + 2)
    (; ngrs, ngre, K, nK = length(K), rn0 = copy(v(:rn)), dn0 = copy(v(:dn)), un0 = copy(v(:un)),
       ro = v(:ro), dO = v(:do), uo = v(:uo), dso = v(:dso), dvolo = v(:dvolo), unom = v(:unom),
       θ = TD.g(:thet), dt = TD.g(:dtime), cadv = TD.g(:cadv), eps = TD.g(:epsadv),
       sig = TD.g(:sigd), mu = gi(:mu), phil0 = TD.g(:phil0), phir0 = TD.g(:phir0), leobc = gi(:leobc))
end

struct ContinTag end
"Overwrite rows (id, k = ngrs..ngre) of the assembled system with Julia R and ∂R/∂x.
The Jacobian uses stencil coloring: unknowns (var v, zone index n) are seeded in
partial slot 5v + mod(n-1, 5) + 1, so one evaluation with 15 partials recovers the
whole ±2-zone band exactly (each row touches each color at most once per variable)."
function contin_rows!()
    c = contin_ctx()
    nK = c.nK
    x0 = vcat(log.(c.rn0[c.K]), log.(c.dn0[c.K]), c.un0[c.K] ./ c.unom[c.K])
    D = ForwardDiff.Dual{ContinTag,Float64,15}
    xd = Vector{D}(undef, 3nK)
    for v in 0:2, n in 1:nK
        slot = 5v + mod(n - 1, 5) + 1
        xd[v*nK + n] = D(x0[v*nK + n], ForwardDiff.Partials(ntuple(i -> i == slot ? 1.0 : 0.0, 15)))
    end
    Rd = contin_res(xd, c)
    id = gi(:id); cols = (gi(:jr), gi(:jd), gi(:ju))
    E = (v(:em2), v(:em1), v(:e00), v(:ep1), v(:ep2)); rhs = v(:rhs); neqn = gi(:neqn)
    @inbounds for (ir, k) in enumerate(c.ngrs:c.ngre)
        rhs[id, k] = ForwardDiff.value(Rd[ir])
        P = ForwardDiff.partials(Rd[ir])
        for o in -2:2, j in 1:neqn; E[o+3][id, j, k] = 0.0; end
        for (vv, j) in enumerate(cols), o in -2:2
            n = k + o - first(c.K) + 1
            (1 <= n <= nK) || continue
            E[o+3][id, j, k] = P[5(vv - 1) + mod(n - 1, 5) + 1]
        end
    end
    nothing
end

"Dense ForwardDiff reference for the colored version (testing only)."
function contin_jacobian_dense()
    c = contin_ctx()
    x0 = vcat(log.(c.rn0[c.K]), log.(c.dn0[c.K]), c.un0[c.K] ./ c.unom[c.K])
    ForwardDiff.jacobian(x -> contin_res(x, c), x0), c
end

# ------------------------------------------------ mass-equation rows in Julia (AD)
"Cumulative-mass residual (mass.f) for interfaces k = ngrs..ngre+1 (unknowns δln r, δln m_ext, δln ρ about the current state)."
function mass_res(x::AbstractVector{T}, c) where {T}
    (; ngrs, ngre, K, nK, rn0, xme0, dn0, ro, xmeo, θ, dt, mu, phil0) = c
    rn = Vector{T}(rn0); xm = Vector{T}(xme0); dn = Vector{T}(dn0)
    @inbounds for (n, k) in enumerate(K)
        # multiplicative perturbation about the current values (TITAN's δln variables;
        # works for the negative exterior mass that can occur in the outer phantom zone)
        rn[k] = rn0[k] * exp(x[n]); xm[k] = xme0[k] * exp(x[nK+n]); dn[k] = dn0[k] * exp(x[2nK+n])
    end
    R = zeros(T, ngre - ngrs + 2)
    k = ngrs                                          # inner boundary (bc14 / bc65)
    rk = θ * rn[k] + (1 - θ) * ro[k]
    R[1] = xm[k] - xmeo[k] + rk^mu * phil0 * dt
    @inbounds for k in ngrs+1:ngre+1                 # m14
        dv = (rn[k]^(mu+1) - rn[k-1]^(mu+1)) / (mu + 1)
        R[k-ngrs+1] = xm[k-1] - xm[k] - dn[k-1] * dv
    end
    R
end
struct MassTag end
function mass_rows!()
    ngrs = gi(:ngrs); ngre = gi(:ngre); K = (ngrs - 2):(ngre + 2)
    c = (; ngrs, ngre, K, nK = length(K), rn0 = copy(v(:rn)), xme0 = copy(v(:xmen)), dn0 = copy(v(:dn)),
           ro = v(:ro), xmeo = v(:xmeo), θ = TD.g(:thet), dt = TD.g(:dtime), mu = gi(:mu), phil0 = TD.g(:phil0))
    nK = c.nK
    x0 = zeros(3nK)
    D = ForwardDiff.Dual{MassTag,Float64,15}
    xd = [D(x0[i], ForwardDiff.Partials(ntuple(j -> j == 5((i - 1) ÷ nK) + mod((i - 1) % nK, 5) + 1 ? 1.0 : 0.0, 15)))
          for i in 1:3nK]
    Rd = mass_res(xd, c)
    im = gi(:im); cols = (gi(:jr), gi(:jm), gi(:jd))
    E = (v(:em2), v(:em1), v(:e00), v(:ep1), v(:ep2)); rhs = v(:rhs); neqn = gi(:neqn)
    @inbounds for (ir, k) in enumerate(ngrs:ngre+1)
        rhs[im, k] = ForwardDiff.value(Rd[ir]); P = ForwardDiff.partials(Rd[ir])
        for o in -2:2, j in 1:neqn; E[o+3][im, j, k] = 0.0; end
        for (vv, j) in enumerate(cols), o in -2:2
            n = k + o - first(K) + 1
            (1 <= n <= nK) || continue
            E[o+3][im, j, k] = P[5(vv - 1) + mod(n - 1, 5) + 1]
        end
    end
    nothing
end

end # module
