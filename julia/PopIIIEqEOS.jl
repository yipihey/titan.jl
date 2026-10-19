"""
    PopIIIEqEOS

Equilibrium-chemistry equation of state for ρ ≥ ρ_eq (RHFC02 App. A1.2–A1.3): Saha (H₂ ⇌ 2H, H ⇌ H⁺ + e⁻; ideal gas,
PrimordialChem.eos_rhfc02) for p ≤ 1e4 dyn cm⁻², SCVH95 (SCVH.jl) for p ≥ 1e6, interpolated in between with a smooth
weight in log p (RHFC02: "appropriately interpolate ... to ensure a smooth transition"; the form of the weight is ours).

Tabulated once on a regular grid in (log ρ, log T) and interpolated bicubically (Catmull–Rom), so that p(ρ, T), u(ρ, T)
are smooth and ForwardDiff Duals pass through (TITAN's Newton sees the softening by dissociation and ionisation).
Stored: log p, log(u + u_off) with u_off = (χ_H + χ_H₂/2)/μ_H (u + u_off = energy relative to cold H₂ and neutral He,
SCVH's zero; positive), f_H₂ (fraction of H nuclei in H₂), x_H⁺.

SCVH composition columns read as number fractions of H₂, H among the hydrogen particles (H₂, H, H⁺) [interp, checked
against Saha where both apply: X_H₂ = 0.604 vs Saha's 0.579 at 3000 K, log P = 5]: f_H₂ = 2X_H₂/(1 + X_H₂),
x_H⁺ = (1 − X_H₂ − X_H)/(1 + X_H₂).
"""
module PopIIIEqEOS

import ..PrimordialChem, ..SCVH
const PC = PrimordialChem
using Serialization, ForwardDiff

const UOFF = (PC.χH + PC.χH2 / 2) / PC.μH

struct EqTable
    lρ0::Float64; dlρ::Float64; nρ::Int
    lT0::Float64; dlT::Float64; nT::Int
    lp::Matrix{Float64}; lus::Matrix{Float64}; fH2::Matrix{Float64}; xHp::Matrix{Float64}; w::Matrix{Float64}
end

smoothstep(x) = (x = clamp(x, 0.0, 1.0); x * x * (3 - 2x))

"One node: Saha-ideal, SCVH where the ideal pressure is ≥ 1e4 and the point is inside SCVH's T range, blended."
function node(lρ, lT)
    ρ = 10^lρ; T = 10^lT
    y = PC.y_saha(ρ, T); u1, p1 = PC.eos_rhfc02(ρ, T, y); f1 = 2y[3]; xp1 = y[1]
    w = (lT >= 2.10 && lT <= 7.06) ? smoothstep((log10(p1) - 4) / 2) : 0.0
    if w > 0
        s = try SCVH.state(ρ, T) catch; nothing end
        s === nothing && (w = 0.0)
        if w > 0
            X2, X1 = s.xH2, s.xH; f2 = 2X2 / (1 + X2); xp2 = max(1 - X2 - X1, 0.0) / (1 + X2)
            p = (1 - w) * p1 + w * s.p; us = (1 - w) * (u1 + UOFF) + w * s.u
            return (log10(p), log10(us), (1 - w) * f1 + w * f2, (1 - w) * xp1 + w * xp2, w)
        end
    end
    (log10(p1), log10(u1 + UOFF), f1, xp1, 0.0)
end

function build(; lρ0 = -12.0, lρ1 = 1.5, dlρ = 0.05, lT0 = 1.7, lT1 = 6.0, dlT = 0.02)
    nρ = round(Int, (lρ1 - lρ0) / dlρ) + 1; nT = round(Int, (lT1 - lT0) / dlT) + 1
    lp = zeros(nρ, nT); lus = zeros(nρ, nT); f = zeros(nρ, nT); xp = zeros(nρ, nT); w = zeros(nρ, nT)
    Threads.@threads for i in 1:nρ
        for j in 1:nT
            lp[i, j], lus[i, j], f[i, j], xp[i, j], w[i, j] = node(lρ0 + (i - 1) * dlρ, lT0 + (j - 1) * dlT)
        end
    end
    EqTable(lρ0, dlρ, nρ, lT0, dlT, nT, lp, lus, f, xp, w)
end

const CACHE = joinpath(homedir(), "Projects", "popiii_lit", "SCVH95", "eqeos_table.jls")
const TAB = Ref{Any}(nothing)
"The table (built once and cached outside the repository; SCVH95 tables are not redistributable from here)."
function table()
    TAB[] === nothing || return TAB[]::EqTable
    TAB[] = isfile(CACHE) ? deserialize(CACHE) : (t = build(); serialize(CACHE, t); t)
    TAB[]::EqTable
end

# Catmull–Rom bicubic interpolation (clamped at the edges); x, y in grid units, generic in the number type
@inline cr(p0, p1, p2, p3, t) = p1 + 0.5t * (p2 - p0 + t * (2p0 - 5p1 + 4p2 - p3 + t * (3(p1 - p2) + p3 - p0)))
function bicubic(A, x, y)
    n, m = size(A)
    xv = ForwardDiff.value(x); yv = ForwardDiff.value(y)
    i = clamp(floor(Int, xv), 1, n - 2); j = clamp(floor(Int, yv), 1, m - 2)
    tx = x - i; ty = y - j
    col(jj) = (jc = clamp(jj, 1, m); cr(A[clamp(i - 1, 1, n), jc], A[i, jc], A[i+1, jc], A[min(i + 2, n), jc], tx))
    cr(col(j - 1), col(j), col(j + 1), col(j + 2), ty)
end
gx(t, lρ) = (lρ - t.lρ0) / t.dlρ + 1
gy(t, lT) = (lT - t.lT0) / t.dlT + 1

"(u [erg/g, our zero], p) at (ρ, T) in equilibrium."
function eos(ρ, T)
    t = table(); x = gx(t, log10(ρ)); y = gy(t, log10(T))
    (10^bicubic(t.lus, x, y) - UOFF, 10^bicubic(t.lp, x, y))
end
"Equilibrium composition y = (x_H⁺, 0, x_H₂, 0) at (ρ, T) (bilinear, Float64)."
function composition(ρ, T)
    t = table(); x = gx(t, log10(ρ)); y = gy(t, log10(T))
    i = clamp(floor(Int, x), 1, t.nρ - 1); j = clamp(floor(Int, y), 1, t.nT - 1); a = clamp(x - i, 0, 1); b = clamp(y - j, 0, 1)
    bl(A) = (1 - a) * (1 - b) * A[i, j] + a * (1 - b) * A[i+1, j] + (1 - a) * b * A[i, j+1] + a * b * A[i+1, j+1]
    (bl(t.xHp), 0.0, bl(t.fH2) / 2, 0.0)
end
"T at fixed (ρ, u − dt Λ/ρ) in equilibrium (bisection in ln T); returns (y, T)."
function step(ρ, u, T0, dt; cooling = nothing)
    F(T) = (y = composition(ρ, T); eos(ρ, T)[1] - u + (cooling === nothing ? 0.0 : dt * cooling(ρ, T, y) / ρ))
    lo, hi = T0 / 1.5, T0 * 1.5
    while F(lo) > 0; lo /= 1.5; lo < 30.0 && error("PopIIIEqEOS.step: no bracket below"); end
    while F(hi) < 0; hi *= 1.5; hi > 1e6 && error("PopIIIEqEOS.step: no bracket above"); end
    for _ in 1:200; m = sqrt(lo * hi); F(m) > 0 ? (hi = m) : (lo = m); hi / lo - 1 < 1e-13 && break; end
    T = sqrt(lo * hi); (composition(ρ, T), T)
end

end # module
