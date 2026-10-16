# Step 4 of phase 2: Julia port of the continuity residual (contin.f + advectc.f
# + the mass-diffusion branch of diffuse.f) with a ForwardDiff Jacobian, compared
# with TITAN's hand-coded Jacobian blocks on a mid-Newton benchmark state.
# usage: julia --project=env contin_ad.jl NAME NN STEP
using ForwardDiff, Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); using .TitanDriver
const TD = TitanDriver; const TF = TD.TF

name, nn, nstep = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3])
TF.init(nn; rundir = joinpath(TF.ROOT, "lib", "p2_runs", "$(name)_contin"))
for i in 1:nstep-1; TD.titan_step!(); end
# one Newton iteration of the next step -> a consistent, unconverged state
TD.r(:peq); TD.r(:matgen); TD.r(:matslv); TD.r(:maxdel); TD.r(:update)
iter_state_t = TF.getscalar(:timen)
TD.fcall(:matgen_)                                   # Fortran residual + Jacobian here
ix = (id = Int(TF.getscalar(:id)), jr = Int(TF.getscalar(:jr)), jd = Int(TF.getscalar(:jd)),
      ju = Int(TF.getscalar(:ju)))
Fblk = Dict(o => copy(TF.var(s)) for (o, s) in ((-2, :em2), (-1, :em1), (0, :e00), (1, :ep1), (2, :ep2)))
Frhs = copy(TF.var(:rhs))

# ------------------------------------------------------------ inputs (COMMON)
ngrs = Int(TF.getscalar(:ngrs)); ngre = Int(TF.getscalar(:ngre))
c = (; thet = TF.getscalar(:thet), dt = TF.getscalar(:dtime), cadv = TF.getscalar(:cadv),
      eps = TF.getscalar(:epsadv), sig = TF.getscalar(:sigd), mu = Int(TF.getscalar(:mu)),
      phil0 = TF.getscalar(:phil0), phir0 = TF.getscalar(:phir0), leobc = Int(TF.getscalar(:leobc)))
rn0, dn0, un0 = copy(TF.var(:rn)), copy(TF.var(:dn)), copy(TF.var(:un))
ro, dO, uo, dso, dvolo, unom = (copy(TF.var(s)) for s in (:ro, :do, :uo, :dso, :dvolo, :unom))
K = (ngrs - 2):(ngre + 2); nK = length(K)
pack() = vcat(log.(rn0[K]), log.(dn0[K]), un0[K] ./ unom[K])

"Continuity residual rows k = ngrs..ngre as a function of (ln r, ln ρ, u/unom) on K."
function contin_res(x::AbstractVector{T}) where {T}
    mgr = length(rn0)
    rn = Vector{T}(undef, mgr); dn = similar(rn); un = similar(rn)
    rn .= rn0; dn .= dn0; un .= un0
    rn[K] = exp.(x[1:nK]); dn[K] = exp.(x[nK+1:2nK]); un[K] = x[2nK+1:3nK] .* unom[K]
    θ = c.thet; mu = c.mu
    r  = θ .* rn .+ (1 - θ) .* ro
    u  = θ .* un .+ (1 - θ) .* uo
    d  = θ .* dn .+ (1 - θ) .* dO
    urel = u .- (rn .- ro) ./ c.dt
    rmu = r .^ mu
    dvoln = similar(rn); for k in ngrs-1:ngre+1; dvoln[k] = (rn[k+1]^(mu+1) - rn[k]^(mu+1)) / (mu + 1); end
    # advectc: van Leer interface densities db(k), k = ngrs+1..ngre+1
    dq = zeros(T, mgr); for k in ngrs-1:ngre+1; dq[k] = dn[k+1] - dn[k]; end
    qsn = zeros(T, mgr)
    for k in ngrs:ngre+1
        sdq = dq[k-1] + dq[k]
        den = sdq + (sdq >= 0 ? c.eps : -c.eps)
        qr = c.cadv * dq[k-1] * dq[k] / den
        qsn[k] = dq[k-1] * dq[k] >= 0 ? qr : zero(T)
    end
    db = zeros(T, mgr)
    for k in ngrs+1:ngre+1
        s = urel[k] >= 0 ? 0.5 : -0.5
        db[k] = (0.5 + s) * (θ * (dn[k-1] + 0.5qsn[k-1]) + (1 - θ) * (dO[k-1] + 0.5dso[k-1])) +
                (0.5 - s) * (θ * (dn[k]   - 0.5qsn[k])   + (1 - θ) * (dO[k]   - 0.5dso[k]))
    end
    # diffuse (mass branch, zet = 0)
    df = zeros(T, mgr)
    if c.sig != 0
        for k in ngrs:ngre+1
            df[k] = 2c.sig * rmu[k] * (d[k] - d[k-1]) / (r[k+1] - r[k-1])
        end
    end
    R = zeros(T, ngre - ngrs + 1)
    for k in ngrs:ngre
        flux_in  = k == ngrs ? rmu[k] * c.phil0 : rmu[k] * urel[k] * db[k]
        flux_out = k == ngre ? (c.leobc == 3 ? rmu[k+1] * c.phir0 + rmu[k+1] * urel[k+1] * db[k+1] :
                                rmu[k+1] * c.phir0) : rmu[k+1] * urel[k+1] * db[k+1]
        R[k-ngrs+1] = (dn[k] * dvoln[k] - dO[k] * dvolo[k]) / c.dt + flux_out - flux_in - df[k+1] + df[k]
    end
    R
end

x0 = pack()
Rj = contin_res(x0)
Rf = [Frhs[ix.id, k] for k in ngrs:ngre]
scale = maximum(abs, Rf)
@printf("%s step %d (state after 1st Newton iteration, t = %.6e, dt = %.3e, sigd = %g)\n",
        name, nstep, iter_state_t, c.dt, c.sig)
termscale = maximum(abs(dn0[k] * (rn0[k+1]^(c.mu+1) - rn0[k]^(c.mu+1)) / (c.mu+1) / c.dt) for k in ngrs:ngre)
@printf("residual (tier a): max |R_jl - R_f| / max|R_f| = %.3e   (max|R_f| = %.3e); relative to the largest term (ρΔV/Δt = %.3e): %.3e\n",
        maximum(abs.(Rj .- Rf)) / scale, scale, termscale, maximum(abs.(Rj .- Rf)) / termscale)

J = ForwardDiff.jacobian(contin_res, x0)
# compare every entry of rows k = ngrs..ngre, columns (var, zone) on K
vars = ((:r, ix.jr, 0), (:d, ix.jd, nK), (:u, ix.ju, 2nK))
worst = 0.0; nent = 0; dropped = Tuple[]; extra = Tuple[]; big = Tuple[]
for (ir, k) in enumerate(ngrs:ngre)
    rowscale = maximum(abs, J[ir, :])
    for (vn, j, off) in vars, (ic, kk) in enumerate(K)
        a = J[ir, off+ic]
        o = kk - k
        f = abs(o) <= 2 ? Fblk[o][ix.id, j, k] : 0.0
        (a == 0 && f == 0) && continue
        global nent += 1
        d = abs(a - f) / max(abs(a), abs(f))
        if f == 0 && abs(a) > 1e-12rowscale
            push!(dropped, (k, vn, o, a, abs(a) / rowscale))
        elseif a == 0 && abs(f) > 1e-12rowscale
            push!(extra, (k, vn, o, f))
        elseif abs(a - f) > 1e-10 * rowscale
            push!(big, (k, vn, o, a, f, d, abs(a - f) / rowscale))
        end
        abs(a - f) <= 1e-10rowscale && (global worst = max(worst, abs(a - f) / rowscale))
    end
end
@printf("Jacobian: %d nonzero entries compared (rows id, cols r/ρ/u, offsets -2..+2)\n", nent)
@printf("  agreeing entries: max |J_ad - J_f| / rowmax = %.3e\n", worst)
println("  entries AD≠0 but Fortran=0 (dropped by TITAN): ", length(dropped))
for e in first(dropped, 12); @printf("     k=%d d/d%s offset %+d : AD %.4e  (|.|/rowmax = %.2e)\n", e...); end
println("  entries Fortran≠0 but AD=0: ", length(extra))
for e in first(extra, 12); println("     ", e); end
println("  entries differing by > 1e-10 rowmax: ", length(big))
for e in first(sort(big, by = e -> -e[end]), 12)
    @printf("     k=%d d/d%s offset %+d : AD %.6e  F %.6e  rel %.2e  |Δ|/rowmax %.2e\n", e...)
end
TF.finalize()
