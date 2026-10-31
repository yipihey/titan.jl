"""
    Rezone

Population III project, S6-C2: conservative rezoning of a spherical Lagrangian mass grid (regular centre). Works on a
zone profile — interface radii `rI` (n+1), zone densities `ρ`, specific gas energies `e`, temperatures `T`, radiation energy
densities `Er` (n), interface velocities `u` and radiative fluxes `Fr` (n+1) — and returns a new profile plus, for every new
zone, the old zones it came from with mass weights (for the material's per-zone state). The TITAN state is then rebuilt
from the new profile by the caller (as for the RAD rebuild), so no COMMON array is edited here.

**Split** of zone k into two zones of equal mass with a slope-limited (minmod) linear reconstruction: ρ_in/ρ_out from the
limited ln ρ slope (the interface where the inner half holds Δm/2), e_in,out = e(1 ∓ δ), E_r volume-weighted with the limited
ln E_r slope; T of each half from the material EOS (`Tofe`). A piecewise-constant split leaves a staircase whose second
differences re-trigger the indicator at the edge of the refined band. Mass, internal and radiation energy are conserved exactly; u and F_r at the new
interface are interpolated linearly in mass. The kinetic energy of the staggered grid changes by the new interface's
½ m u² minus the change of the neighbours' interface masses: the difference is taken from (or given to) the internal
energy of the two halves, so the total energy is conserved to round-off.
**Merge** of zones k, k+1: the interface between them is removed; ρ = (Δm_k + Δm_{k+1})/ΔV, e mass-weighted, E_r
volume-weighted (radiation energy conserved); the kinetic energy of the removed interface (and the change of the
neighbours' interface masses) goes into e; T from the material's EOS at the new (ρ, e, composition) via `Tofe`.
Momentum Σ m_i u_i (interface masses) changes at second order in the velocity variation across the zones involved:
merges are allowed only where that variation is small (see `plan`).

**Criteria** (`plan`): a zone is split when its indicator J_k (`jumps`: second differences of ln ρ, ln T; compression
/c_s) exceeds `jsplit` and its mass exceeds both `dmmin` and `fmin` × the mass it encloses (at most `nhard` zones in all, largest J first); then, until
none is left, any zone heavier than `qmax` × a neighbour is split (a smooth mass grid). Pairs are merged, smoothest first,
only while the zone count exceeds `nmax`, when J of both and of their neighbours is below `jmerge`, the merged mass stays
within `qmax` of both neighbours, and the pair is outside the innermost `keepcentre` zones.
"""
module Rezone

export Profile, plan, apply, totals

minmod(a, b) = a * b <= 0 ? 0.0 : (abs(a) < abs(b) ? a : b)

struct Profile
    rI::Vector{Float64}; u::Vector{Float64}; Fr::Vector{Float64}     # n + 1 interfaces (rI[1] = 0: regular centre)
    ρ::Vector{Float64}; e::Vector{Float64}; T::Vector{Float64}; Er::Vector{Float64}   # n zones
    cs::Vector{Float64}                                              # sound speed per zone (criteria only)
    P::Vector{Float64}                                               # gas pressure per zone (pressure-consistent split)
end
Profile(rI, u, Fr, ρ, e, T, Er, cs) = Profile(rI, u, Fr, ρ, e, T, Er, cs, fill(NaN, length(ρ)))
nz(p::Profile) = length(p.ρ)
vol(rI, k) = 4π / 3 * (rI[k+1]^3 - rI[k]^3)
dmass(p::Profile) = [p.ρ[k] * vol(p.rI, k) for k in 1:nz(p)]
"Interface masses of the staggered grid: ½(Δm_{k−1} + Δm_k); the centre and the surface get half a zone."
imass(dm) = (n = length(dm); [0.5 * ((i > 1 ? dm[i-1] : 0.0) + (i <= n ? dm[i] : 0.0)) for i in 1:n+1])

"Totals: mass, internal energy, radiation energy, kinetic energy (interface masses), momentum Σ m_i |u_i| r̂ (radial, signed)."
function totals(p::Profile)
    dm = dmass(p); mi = imass(dm)
    (; M = sum(dm), Eint = sum(dm .* p.e), Erad = sum(p.Er[k] * vol(p.rI, k) for k in 1:nz(p)),
       Ekin = sum(0.5 .* mi .* p.u .^ 2), P = sum(mi .* p.u))
end

"Totals over zones lo:hi (interfaces lo:hi+1 for the kinetic energy and momentum)."
function totals(p::Profile, lo, hi)
    dm = dmass(p); mi = imass(dm)
    (; M = sum(dm[lo:hi]), Eint = sum(dm[lo:hi] .* p.e[lo:hi]), Erad = sum(p.Er[k] * vol(p.rI, k) for k in lo:hi),
       Ekin = sum(0.5 .* mi[lo:hi+1] .* p.u[lo:hi+1] .^ 2), P = sum(mi[lo:hi+1] .* p.u[lo:hi+1]))
end

"""
Refinement indicator per zone: the largest of |second difference| of ln ρ and of ln T over the zone and its neighbours
(≈ 0 on a smooth power law, whatever the zone-to-zone ratio of the geometric grid) and the compression across the zone
max(0, u_k − u_{k+1})/c_s (shocks); the indicator of a zone is the max over itself and its two neighbours.
"""
function jumps(p::Profile)
    n = nz(p); J0 = zeros(n); lρ = log.(p.ρ); lT = log.(p.T)
    for k in 2:n-1
        J0[k] = max(abs(lρ[k+1] - 2lρ[k] + lρ[k-1]), abs(lT[k+1] - 2lT[k] + lT[k-1]))
    end
    for k in 1:n; J0[k] = max(J0[k], max(0.0, p.u[k] - p.u[k+1]) / max(p.cs[k], 1e-30)); end
    [maximum(J0[max(1, k - 1):min(n, k + 1)]) for k in 1:n]
end

"""
    plan(p; jsplit = 0.3, jmerge = 0.05, dmmin, qmax = 1.5, nmax, nhard, keepcentre = 5) -> (splits::Vector{Int}, merges::Vector{Int})

Zones to split (indices in `p`) and the left zones of the pairs to merge (indices in `p`, disjoint from each other and from
the splits). See the module docstring.
"""
function plan(p::Profile; jsplit = 0.3, jmerge = 0.05, dmmin = 0.0, qmax = 1.5, nmax = typemax(Int), nhard = typemax(Int), keepcentre = 5, evenn = true, fmin = 0.0)
    n = nz(p); dm = dmass(p); J = jumps(p)
    split = falses(n)
    menc = cumsum(dm)                                # mass enclosed by the outer interface of each zone
    for k in 1:n; (J[k] > jsplit && dm[k] > max(2dmmin, fmin * menc[k])) && (split[k] = true); end
    # smooth mass grid: a zone heavier than qmax × a neighbour (after the neighbour's own split) is split too
    for _ in 1:50
        eff = [split[k] ? dm[k] / 2 : dm[k] for k in 1:n]; changed = false
        for k in 1:n
            split[k] && continue
            for l in (k - 1, k + 1)
                (1 <= l <= n) || continue
                if eff[k] > qmax * eff[l] && dm[k] > 2dmmin
                    split[k] = true; changed = true; break
                end
            end
        end
        changed || break
    end
    nnew = n + count(split)
    if nnew > nhard                                  # hard limit (array size): keep the splits of the largest indicators
        ks = sort(findall(split); by = k -> -J[k])
        for k in ks[nhard - n + 1:end]; split[k] = false; end
        nnew = n + count(split)
    end
    merges = Int[]
    if nnew > nmax
        eff = [split[k] ? dm[k] / 2 : dm[k] for k in 1:n]
        cand = [k for k in max(keepcentre + 1, 2):n-2 if !split[k] && !split[k+1] && max(J[k-1], J[k], J[k+1], J[k+2]) < jmerge]
        sort!(cand; by = k -> max(J[k], J[k+1]))
        used = falses(n)
        for k in cand
            nnew - length(merges) <= nmax && break
            (used[k] || used[k+1] || used[k-1] || used[k+2]) && continue
            mm = eff[k] + eff[k+1]
            (mm <= qmax * eff[k-1] && mm <= qmax * eff[k+2]) || continue
            push!(merges, k); used[k] = used[k+1] = true
        end
    end
    sp = findall(split)
    if isodd(n + length(sp) - length(merges)) && evenn         # TITAN's builder needs an even zone count
        if !isempty(sp)
            deleteat!(sp, argmin(J[sp]))
        elseif !isempty(merges)
            pop!(merges)
        end
    end
    sp, sort!(merges)
end

"""
    apply(p, splits, merges; Tofe) -> (q::Profile, src)

New profile `q` and, for each new zone j, `src[j]` = vector of (old zone, mass weight) pairs (weights sum to 1).
`Tofe(ρ, e, T0, src_j)` returns T for a merged or split zone (material EOS with the mass-weighted composition of `src_j`);
`eosT(ρ, T, src_j) -> (e, P)` (optional) switches on the pressure-consistent split.
Total mass, internal + kinetic + radiation energy are conserved to round-off (the kinetic-energy change of the staggered
grid is moved into the internal energy of the zones involved).
"""
function apply(p::Profile, splits, merges; Tofe, eosT = nothing)
    n = nz(p); dm = dmass(p); isS = falses(n); isS[splits] .= true; isM = falses(n); isM[merges] .= true
    rI = [p.rI[1]]; u = [p.u[1]]; Fr = [p.Fr[1]]
    ρ = Float64[]; e = Float64[]; T = Float64[]; Er = Float64[]; cs = Float64[]; P = Float64[]; src = Vector{Vector{Tuple{Int,Float64}}}()
    kept = Int[]                                    # for every new interface: the old interface it is (0 = new)
    push!(kept, 1)
    k = 1
    while k <= n
        if isM[k]
            r0, r2 = p.rI[k], p.rI[k+2]; V = 4π / 3 * (r2^3 - r0^3); m = dm[k] + dm[k+1]
            push!(ρ, m / V); push!(e, (dm[k] * p.e[k] + dm[k+1] * p.e[k+1]) / m)
            push!(Er, (p.Er[k] * vol(p.rI, k) + p.Er[k+1] * vol(p.rI, k + 1)) / V)
            push!(cs, max(p.cs[k], p.cs[k+1])); push!(P, NaN); push!(src, [(k, dm[k] / m), (k + 1, dm[k+1] / m)])
            push!(T, NaN)                           # from Tofe below, after the energy correction
            push!(rI, r2); push!(u, p.u[k+2]); push!(Fr, p.Fr[k+2]); push!(kept, k + 2)
            k += 2
        elseif isS[k]
            # slope-limited (minmod) linear reconstruction, exactly conservative: halves of equal mass Δm/2 with
            # ρ_in/ρ_out = exp(−g) (g: limited ln ρ difference per half zone), e_in,out = e (1 ± δ), E_r volume-weighted
            V = vol(p.rI, k); m = dm[k]
            g = k == 1 || k == n ? 0.0 : 0.5 * minmod(log(p.ρ[k+1] / p.ρ[k]), log(p.ρ[k] / p.ρ[k-1]))
            x = exp(-g); ρout = m * (1 + 1 / x) / (2V); ρin = x * ρout          # Δm/(2ρin) + Δm/(2ρout) = V
            rh = cbrt(p.rI[k]^3 + 3 / (4π) * m / (2ρin))
            δ = k == 1 || k == n ? 0.0 : 0.25 * minmod((p.e[k+1] - p.e[k]) / p.e[k], (p.e[k] - p.e[k-1]) / p.e[k])
            Tin = Tout = NaN
            if eosT !== nothing && 1 < k < n
                # pressure-consistent split: P_in,out = P̄ ∓ ½ΔP with ΔP from the limited pressure slope in mass across the
                # neighbours (the new interface sees the gradient of its surroundings: no hydrostatic kick), P̄ from the
                # energy constraint e_in + e_out = 2e (equal masses); T of each half from P(ρ_half, T) = P_half
                mc(i) = sum(dm[1:i-1]) + 0.5 * dm[i]
                sP = minmod((p.P[k+1] - p.P[k]) / (mc(k + 1) - mc(k)), (p.P[k] - p.P[k-1]) / (mc(k) - mc(k - 1)))
                ΔP = -sP * m / 2                              # P_in − P_out (half-zone centres Δm/2 apart)
                srck = [(k, 1.0)]
                Tof(ρh, Ph, T0) = (lo = log(T0) - 2; hi = log(T0) + 2;
                    for _ in 1:80; mid = 0.5 * (lo + hi); eosT(ρh, exp(mid), srck)[2] > Ph ? (hi = mid) : (lo = mid); end; exp(0.5 * (lo + hi)))
                function f(Pb)
                    ti = Tof(ρin, Pb + ΔP / 2, p.T[k]); to = Tof(ρout, Pb - ΔP / 2, p.T[k])
                    (eosT(ρin, ti, srck)[1] + eosT(ρout, to, srck)[1] - 2p.e[k], ti, to)
                end
                a_ = 0.98p.P[k]; b_ = 1.02p.P[k]; fa = f(a_)[1]; fb = f(b_)[1]
                if fa * fb < 0 && abs(ΔP) < 0.5p.P[k]
                    for _ in 1:60
                        c_ = 0.5 * (a_ + b_); fc = f(c_)[1]
                        fa * fc <= 0 ? (b_ = c_; fb = fc) : (a_ = c_; fa = fc)
                    end
                    _, Tin, Tout = f(0.5 * (a_ + b_))
                end
            end
            Vin = 4π / 3 * (rh^3 - p.rI[k]^3); Vout = V - Vin
            ε = k == 1 || k == n ? 0.0 : 0.25 * minmod(log(p.Er[k+1] / p.Er[k]), log(p.Er[k] / p.Er[k-1]))
            Ein = p.Er[k] * V / (Vin + Vout * exp(2ε)); Eout = Ein * exp(2ε)  # E_in V_in + E_out V_out = E_r V
            uh = 0.5 * (p.u[k] + p.u[k+1]); Fh = 0.5 * (p.Fr[k] + p.Fr[k+1])   # equal masses: linear in mass = mean
            if isnan(Tin)
                halves = ((ρin, p.e[k] * (1 - δ), Ein, NaN), (ρout, p.e[k] * (1 + δ), Eout, NaN))
            else                                         # e of the halves from the EOS (their sum is 2e to the bisection tolerance)
                ei = eosT(ρin, Tin, [(k, 1.0)])[1]; eo = eosT(ρout, Tout, [(k, 1.0)])[1]; c_ = 2p.e[k] / (ei + eo)
                halves = ((ρin, ei * c_, Ein, NaN), (ρout, eo * c_, Eout, NaN))   # exact energy; T re-inverted below
            end
            for (ρh, eh, Eh, Th) in halves
                push!(ρ, ρh); push!(e, eh); push!(T, Th); push!(Er, Eh); push!(cs, p.cs[k]); push!(P, NaN); push!(src, [(k, 1.0)])
            end
            push!(rI, rh); push!(u, uh); push!(Fr, Fh); push!(kept, 0)
            push!(rI, p.rI[k+1]); push!(u, p.u[k+1]); push!(Fr, p.Fr[k+1]); push!(kept, k + 1)
            k += 1
        else
            push!(ρ, p.ρ[k]); push!(e, p.e[k]); push!(T, p.T[k]); push!(Er, p.Er[k]); push!(cs, p.cs[k]); push!(P, p.P[k]); push!(src, [(k, 1.0)])
            push!(rI, p.rI[k+1]); push!(u, p.u[k+1]); push!(Fr, p.Fr[k+1]); push!(kept, k + 1)
            k += 1
        end
    end
    q = Profile(rI, u, Fr, ρ, e, T, Er, cs, P)
    # kinetic energy: total before − after, moved into the internal energy of the zones next to every changed interface
    # (in proportion to their masses), so Eint + Ekin is conserved
    dmq = dmass(q); miq = imass(dmq); mip = imass(dm)
    ΔK = sum(0.5 .* mip .* p.u .^ 2) - sum(0.5 .* miq .* q.u .^ 2)
    touched = Int[]
    for j in 1:length(q.ρ)
        (length(src[j]) > 1 || isS[src[j][1][1]]) && push!(touched, j)
    end
    if !isempty(touched) && ΔK != 0
        Mt = sum(dmq[touched])
        for j in touched; q.e[j] += ΔK / Mt; end
    end
    for j in 1:length(q.ρ)
        if isnan(q.T[j])                            # merged and split zones: T from the material EOS
            q.T[j] = Tofe(q.ρ[j], q.e[j], sum(w * p.T[i] for (i, w) in src[j]), src[j])
        end
    end
    q, src
end

end # module
