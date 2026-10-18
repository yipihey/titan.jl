# Line profiles of an axisymmetric envelope: equatorial infall plus bipolar outflow, seen at several inclinations.
# A kinematic model on top of the spherical non-LTE solution (nlte_run.jl), not a two-dimensional hydrodynamics run:
# inside the two polar cones |cos θ| > cos θ_o the gas moves outwards at the free-fall speed (u = +v_ff r̂), elsewhere
# it falls in (u = −v_ff r̂). Densities, level populations, line opacities and source functions are those of the
# spherical solution at the same radius. In the Sobolev approximation these depend on the velocity field only through
# |dv/dr| and v/r, which reversing the flow leaves unchanged; that the radiation field is the same is an assumption.
# The continuum is the spherical one. Rays: the impact parameters of the spherical transfer × NPHI azimuths.
#   Z0_TABLE=ext TITAN_NOLIB=1 julia -t 2 --project=env axisym_lines.jl NAME [lines=Halpha,Hbeta] [theta=0,30,60,90] [incl=0,45,90] [nphi=16]
# θ_o = 0: pure infall (the spherical model); θ_o = 90: pure outflow. Inclination i: angle between the line of sight
# and the outflow axis. Writes lib/spectra/axisym_NAME.jls with every profile (F/F_c vs v, positive = redshift).
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf
const SP = Spectrum; const Z0 = Z0Physics; const NH = NLTEH
name = ARGS[1]
opt = Dict(Symbol(split(a, "=")[1]) => split(a, "=")[2] for a in ARGS[2:end])
lnames = split(get(opt, :lines, "Halpha,Hbeta"), ","); thetas = parse.(Float64, split(get(opt, :theta, "0,30,60,90"), ","))
incls = parse.(Float64, split(get(opt, :incl, "0,45,90"), ",")); nphi = parse(Int, get(opt, :nphi, "16"))
d = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "$name.jls"))
dn = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$name.jls"))
S = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(d.run)_relax_snaps.jls"))
s = S[argmin([abs(x.t - d.t) for x in S])]
atm = haskey(dn, :Tstar) ? merge(d.atm, (T = dn.Tstar,)) : d.atm
c = Z0.c; kB = Z0.kB; ξturb = 10e5
rsh = s.r[d.ksh+1]; inf = dn.sel; rI = atm.r[inf]; TI = atm.T[inf]; rmax = rI[end]; uabs = dn.vabs
function interp(x, xs, ys)
    x <= xs[1] && return ys[1]; x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x); t = (log(x) - log(xs[i])) / (log(xs[i+1]) - log(xs[i])); ys[i] + t * (ys[i+1] - ys[i])
end
nlfull = Any[nothing for _ in eachindex(atm.r)]; for (q, j) in enumerate(inf); nlfull[j] = dn.nl[q]; end
op = NH.NOpac(SP.Opac(atm), nlfull)
_, R = SP.rays(atm.r; ncore = 12); R = SP.weights!(atm, R)
LU = Dict("Halpha" => (2, 3), "Hbeta" => (2, 4), "Paalpha" => (3, 4))
vgrid = collect(range(-2500e5, 2500e5, length = 201))
n = 1500; gq = [(exp(5 * (i - 1) / (n - 1)) - 1) / (exp(5) - 1) for i in 1:n]
results = []
for ln in lnames
    l, u = LU[ln]; ν0 = (NH.En[u] - NH.En[l]) / Z0.h; λ0 = c / ν0; f = NH.fosc(l, u); gl, gu = NH.gn[l], NH.gn[u]
    _, _, _, ci = SP.monochromatic(atm, op, R, ν0; moments = :I); ps, Ic = ci.p, ci.I
    wp = [ps[j] * ((j > 1 ? ps[j] - ps[j-1] : 0.0) + (ps[j+1] - ps[j])) / 2 for j in 1:length(ps)-1]   # the same ray weights for line and continuum
    Fc = sum(Ic[j] * wp[j] for j in 1:length(ps)-1)
    Kr = [max(NH.πe2mc * f * λ0 * (st.n[l] - gl / gu * st.n[u]), 0.0) for st in dn.nl]
    Sr = [K > 0 ? 2Z0.h * ν0^3 / c^2 / (st.n[l] * gu / (st.n[u] * gl) - 1) : 0.0 for (K, st) in zip(Kr, dn.nl)]
    for θo in thetas, inc in incls
        (θo == 0 || θo == 90) && inc != incls[1] && continue          # spherical cases: one inclination suffices
        ci_ = cosd(θo); ax = (sind(inc), 0.0, cosd(inc))                 # outflow axis in observer coordinates (z → observer)
        F = zeros(length(vgrid))
        Threads.@threads for iv in eachindex(vgrid)
            v = vgrid[iv]; acc = 0.0
            for j in 1:length(ps)-1
                p = ps[j]
                znear = sqrt(max(rmax^2 - p^2, 0.0)); zfar = p < rsh ? sqrt(rsh^2 - p^2) : -znear
                Ij = 0.0
                for k in 1:nphi
                    φ = 2π * (k - 0.5) / nphi; x, y = p * cos(φ), p * sin(φ)
                    I = Ic[j]
                    if znear > zfar && p < rmax
                        zs = p < rsh ? zfar .+ (znear - zfar) .* gq : (zp = znear .* gq; vcat(-reverse(zp), zp[2:end]))
                        prev = nothing
                        for zz in zs
                            r = sqrt(p^2 + zz^2)
                            if r < rsh || r < rI[1]; prev = nothing; continue; end
                            sgn = abs((x * ax[1] + y * ax[2] + zz * ax[3]) / r) > ci_ ? -1.0 : 1.0     # −1: outflow cone
                            w = sgn * interp(r, rI, uabs) * zz / r                    # line-of-sight velocity, + = redshift
                            Δv = sqrt(2kB * interp(r, rI, TI) / Z0.mH + ξturb^2)
                            K = interp(r, rI, Kr); Sl = interp(r, rI, Sr)
                            if prev !== nothing
                                zp_, wp, Δp, Kp, Sp = prev
                                lo, hi = minmax(wp, w); D = max(Δp, Δv)
                                if !(v < lo - 5D || v > hi + 5D) && (Kp > 0 || K > 0)
                                    m = max(1, ceil(Int, (hi - lo) / (0.25 * min(Δp, Δv)))); dz = (zz - zp_) / m
                                    for q in 1:m
                                        t = (q - 0.5) / m; wt = wp + t * (w - wp); Dt = Δp + t * (Δv - Δp)
                                        xx = (v - wt) / Dt; abs(xx) > 6 && continue
                                        τ = (Kp + t * (K - Kp)) * exp(-xx^2) / (sqrt(π) * Dt) * dz
                                        I = I * exp(-τ) + (Sp + t * (Sl - Sp)) * (-expm1(-τ))
                                    end
                                end
                            end
                            prev = (zz, w, Δv, K, Sl)
                        end
                    end
                    Ij += I / nphi
                end
                # trapezoid in p: accumulate this ray's I with half-weights to its two intervals
                acc += Ij * wp[j]
            end
            F[iv] = acc / Fc
        end
        EW = sum(0.5 * ((1 - F[i]) + (1 - F[i+1])) * (vgrid[i+1] - vgrid[i]) for i in 1:length(vgrid)-1) / c * λ0 * 1e8
        k = argmax(F); km = argmin(F)
        @printf("%-6s θ_o = %2.0f° i = %2.0f°: peak %7.2f at %+5.0f km/s, min %.3f at %+5.0f km/s, EW %+8.1f Å\n",
                ln, θo, inc, F[k], vgrid[k] / 1e5, F[km], vgrid[km] / 1e5, EW)
        flush(stdout)
        push!(results, (line = ln, θo = θo, incl = inc, v = vgrid ./ 1e5, F = F, EW = EW))
    end
end
serialize(joinpath(@__DIR__, "..", "lib", "spectra", "axisym_$name.jls"), results)
println("saved lib/spectra/axisym_$name.jls")
