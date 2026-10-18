# Hα and Hβ profiles of the infall, over the continuum computed by spectra_run.jl.
#   Z0_TABLE=ext TITAN_NOLIB=1 julia -t 2 --project=env line_profile.jl NAME     (NAME: lib/spectra/NAME.jls, ext=1 runs)
# [estimate] Assumptions:
#  * Only the free-falling gas outside the shock (TITAN zones above the shock zone, then the free-fall extension) forms
#    the lines; the photosphere's own line (absorption in a static atmosphere) is not computed, so the profile shows the
#    infall's contribution only.
#  * Lower-level population: LTE (Saha–Boltzmann) at the zone's (ρ, T), times a departure coefficient b₂ = 0.1, 1, 10
#    to bracket the unknown NLTE populations (Lyα trapping in this column would push b₂ up).
#  * Line source function, two limits: S_L = B_ν(T) (thermal, LTE) or pure scattering of the photospheric light,
#    S_L = β_c I_*/β (Sobolev escape probabilities, Castor 1970), with I_* the mean continuum intensity of the disc.
#  * Line optical depth integrated along each ray with a Gaussian profile (Doppler width √(2kT/m_H + ξ²), assumed
#    microturbulence ξ = 10 km/s, LP_XI), from the far side to the near side, I ← I e^{−Δτ} + S_L (1 − e^{−Δτ}),
#    starting from the continuum intensity I_c(p) of the full spectrum calculation (the continuum is treated as
#    background). A Sobolev root-finding version was not converged in the ray sampling at the velocity turning points.
#  * Rays with p below the shock radius see only the near side (the star occults the far side).
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl"))
using Serialization, Printf
const SP = Spectrum; const Z0 = Z0Physics
name = ARGS[1]
d = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "$name.jls"))
S = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(d.run)_relax_snaps.jls"))
s = S[argmin([abs(x.t - d.t) for x in S])]
atm = d.atm; M = 1e6 * Z0.Msun; R0 = 8.4352e15
rsh = s.r[d.ksh+1]
# infall points: TITAN zones above the shock, and the extension; |u| at each
# line-forming gas: the infall points with T < 2e4 K (the shock-heated zone just above the shock is ionised and is left
# out, so that S_L is never interpolated between the shocked gas and the cool infall)
inf = findall(j -> (atm.zone[j] == 0 || atm.zone[j] > d.ksh) && atm.T[j] < 2e4, eachindex(atm.r))
uabs = [atm.zone[j] == 0 ? sqrt(2Z0.G * M / atm.r[j]) : abs(0.5 * (s.u[atm.zone[j]] + s.u[atm.zone[j]+1])) for j in inf]
rI = atm.r[inf]; TI = atm.T[inf]; ρI = atm.ρ[inf]
rmax = rI[end]
"interpolate in log r over the infall points (constant beyond the ends)"
function interp(x, xs, ys)
    x <= xs[1] && return ys[1]; x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x); t = (log(x) - log(xs[i])) / (log(xs[i+1]) - log(xs[i]))
    ys[i] + t * (ys[i+1] - ys[i])
end
pops = [begin st = Z0.saha(ρI[k], TI[k]); pH, _ = Z0.level_pops(TI[k], st); pH end for k in eachindex(rI)]
const πe2mc = 0.02654
# line profile: Gaussian, Doppler width Δv = √(2kT/m_H + ξ²) with an assumed microturbulence ξ (default 10 km/s)
const ξturb = parse(Float64, get(ENV, "LP_XI", "10")) * 1e5
LINES = [(name = "Halpha", λ = 6562.8, l = 2, u = 3, f = 0.6407), (name = "Hbeta", λ = 4861.3, l = 2, u = 4, f = 0.1193)]
vgrid = collect(range(-2500e5, 2500e5, length = 251))
bvals = (0.1, 1.0, 10.0)
out = []
op = SP.Opac(atm); _, R = SP.rays(atm.r; ncore = 12); R = SP.weights!(atm, R)
for L in LINES
    ν0 = Z0.c / (L.λ * 1e-8)
    _, _, _, ci = SP.monochromatic(atm, op, R, ν0; moments = :I)
    ps, Ic = ci.p, ci.I
    gl, gu = 2L.l^2, 2L.u^2
    # line opacity factor per unit b: K(r) = πe²/mc f λ n_l (1 − (g_l n_u)/(g_u n_l)) [cm² s⁻¹ cm⁻³ · cm = s⁻¹]
    Kr = [πe2mc * L.f * L.λ * 1e-8 * pops[k][L.l] * (1 - gl * pops[k][L.u] / (gu * pops[k][L.l])) for k in eachindex(rI)]
    SLlte = [SP.Bν(ν0, TI[k]) for k in eachindex(rI)]
    Fc = sum(0.5 * (Ic[j] * ps[j] + Ic[j+1] * ps[j+1]) * (ps[j+1] - ps[j]) for j in 1:length(ps)-1)
    # mean intensity of the stellar disc (rays inside the shock radius), the light the infall scatters
    jin = findall(<(rsh), ps)
    Istar = sum(0.5 * (Ic[j] * ps[j] + Ic[j+1] * ps[j+1]) * (ps[j+1] - ps[j]) for j in jin[1:end-1]) / (0.5 * ps[jin[end]]^2)
    dlnu = [k == 1 ? 0.0 : (log(uabs[k]) - log(uabs[k-1])) / (log(rI[k]) - log(rI[k-1])) for k in eachindex(rI)]
    function SLscat(b)
        μs = range(-1, 1, length = 401)
        map(eachindex(rI)) do k
            r = rI[k]; μc = sqrt(max(1 - (rsh / r)^2, 0.0))
            esc(μ) = (τ = b * Kr[k] / max(abs(uabs[k] / r * ((1 - μ^2) + μ^2 * dlnu[k])), 1e-30); τ < 1e-8 ? 1.0 : -expm1(-τ) / τ)
            β = sum(esc, μs) / length(μs)
            βc = sum(μ -> μ >= μc ? esc(μ) : 0.0, μs) / length(μs)
            βc * Istar / β
        end
    end
    prof = Dict{Tuple{Float64,Symbol},Vector{Float64}}()
    # rays: samples along z (far side to near side; the far side only if the ray misses the star), crowded towards
    # the star, with the redshift velocity w = |u| z/r, the line opacity factor and the LTE source function there
    n = parse(Int, get(ENV, "LP_NZ", "3000")); gq = [(exp(5 * (i - 1) / (n - 1)) - 1) / (exp(5) - 1) for i in 1:n]
    RAYS = map(ps) do p
        znear = sqrt(max(rmax^2 - p^2, 0.0))
        zfar = p < rsh ? sqrt(rsh^2 - p^2) : -znear
        (znear > zfar && p < rmax) || return nothing
        zs = p < rsh ? zfar .+ (znear - zfar) .* gq : (zp = znear .* gq; vcat(-reverse(zp), zp[2:end]))
        rz = sqrt.(p^2 .+ zs .^ 2)
        (z = zs, r = rz, w = [r < rsh ? NaN : interp(r, rI, uabs) * z / r for (z, r) in zip(zs, rz)],
         K = [r < rI[1] ? 0.0 : interp(r, rI, Kr) for r in rz], Slte = [interp(r, rI, SLlte) for r in rz],
         Δv = [sqrt(2Z0.kB * interp(r, rI, TI) / Z0.mH + ξturb^2) for r in rz])
    end
    for b in bvals, src in (:lte, :scat)
        SL = src == :lte ? SLlte : SLscat(b)
        SZ = [ry === nothing ? nothing : [interp(r, rI, SL) for r in ry.r] for ry in RAYS]
        F = zeros(length(vgrid))
        Threads.@threads for iv in eachindex(vgrid)
            v = vgrid[iv]
            Iv = similar(Ic)
            for j in eachindex(ps)
                I = Ic[j]; ry = RAYS[j]
                if ry !== nothing
                    zs, wz, Kz = ry.z, ry.w, ry.K
                    Sz = src == :lte ? ry.Slte : SZ[j]
                    Δvz = ry.Δv
                    for i in 1:length(zs)-1
                        (isfinite(wz[i]) && isfinite(wz[i+1])) || continue
                        # skip segments far from resonance (|v − w| > 5 Doppler widths at both ends and not straddling)
                        lo, hi = minmax(wz[i], wz[i+1]); dv = max(Δvz[i], Δvz[i+1])
                        (v < lo - 5dv || v > hi + 5dv) && continue
                        m = max(1, ceil(Int, (hi - lo) / (0.25 * min(Δvz[i], Δvz[i+1]))))
                        dz = (zs[i+1] - zs[i]) / m
                        for q in 1:m
                            t = (q - 0.5) / m
                            w = wz[i] + t * (wz[i+1] - wz[i]); D = Δvz[i] + t * (Δvz[i+1] - Δvz[i])
                            x = (v - w) / D
                            abs(x) > 6 && continue
                            τ = b * (Kz[i] + t * (Kz[i+1] - Kz[i])) * exp(-x^2) / (sqrt(π) * D) * dz
                            S_ = Sz[i] + t * (Sz[i+1] - Sz[i])
                            I = I * exp(-τ) + S_ * (-expm1(-τ))
                        end
                    end
                end
                Iv[j] = I
            end
            F[iv] = sum(0.5 * (Iv[j] * ps[j] + Iv[j+1] * ps[j+1]) * (ps[j+1] - ps[j]) for j in 1:length(ps)-1) / Fc
        end
        prof[(b, src)] = F
        k = argmin(F); kx = argmax(F)
        EW = sum(0.5 * ((1 - F[i]) + (1 - F[i+1])) * (vgrid[i+1] - vgrid[i]) for i in 1:length(vgrid)-1) / Z0.c * L.λ
        @printf("%-6s %-4s b₂ = %4.1f: minimum %.3f at %+.0f km/s, maximum %.3f at %+.0f km/s, equivalent width %.2f Å (positive = net absorption)\n",
                L.name, src, b, F[k], vgrid[k] / 1e5, F[kx], vgrid[kx] / 1e5, EW)
    end
    push!(out, (name = L.name, λ = L.λ, v = vgrid ./ 1e5,
                lte_b01 = prof[(0.1, :lte)], lte_b1 = prof[(1.0, :lte)], lte_b10 = prof[(10.0, :lte)],
                scat_b01 = prof[(0.1, :scat)], scat_b1 = prof[(1.0, :scat)], scat_b10 = prof[(10.0, :scat)]))
end
serialize(joinpath(@__DIR__, "..", "lib", "spectra", "lines_$(name)$(get(ENV, "LP_TAG", "")).jls"), out)
println("saved lib/spectra/lines_$name.jls")
