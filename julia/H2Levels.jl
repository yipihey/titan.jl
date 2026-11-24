isdefined(@__MODULE__, :H2Lines) || include(joinpath(@__DIR__, "H2Lines.jl"))
"""
    H2Levels

Population III project, G6b: non-LTE level populations of H₂ (the 63 levels v = 0–2, J = 0–20 of `H2Lines`' RHFC02
subset) from statistical equilibrium with collisional and radiative rates.

Levels, degeneracies, A-values: ExoMol RACPPK as in `H2Lines` (the line list is `H2Lines.load()` remapped onto the 63
levels, which are sorted by energy). Collisional de-excitation rate coefficients k(u→l) [cm³ s⁻¹] (Cloudy 25.00
`data/h2`, see `~/Projects/popiii_lit/H2coll/README.md`); excitation by detailed balance
k(l→u) = k(u→l) (g_u/g_l) exp(−(E_u−E_l)/kT) with the ExoMol energies. Every row's direction is checked against
those energies (upper level = the higher one; rows written the other way would be swapped and counted in
`CollTable.nswap`; none occur in these files).

| collider | file | temperature grid | coverage (63-level box) |
| H        | coll_rates_H_07.dat (Wrathmall, Gusdorf & Flower 2007; first 3 rows Sun & Dalgarno 1994 reactive o↔p) | 100–6000 K, 60 pts | 844 transitions (v ≤ 2, J ≤ 20) |
| H₂ (p/o) | coll_rates_H2{para,ortho}_LeBourlot.dat (Le Bourlot, Pineau des Forêts & Flower 1999) | 9 pts 100–6000 K | v ≤ 2, J ≤ 16 (para) / 15 (ortho) |
| He       | coll_rates_He_ORNL.dat (Lee et al. 2007, in prep.) | 41 pts 1–10⁴ K | 529 transitions |
| H⁺       | coll_rates_Hp.dat (Gerlich 1990) | 3, 10, 30, 100, 300, 6000 K | v = 0, J ≤ 9 (incl. o↔p) |
| e⁻       | none: electron collisions are NOT included |

(`coverage(L)` prints the measured numbers.) Temperature interpolation: linear in log k versus log T; outside the table
the end values are held (k(T < T_min) = k(T_min), k(T > T_max) = k(T_max); detailed balance still uses the actual T).
Transitions absent from a collider's file have k = 0 (in both directions); levels without any listed partner are
populated only radiatively/by other colliders. Ortho↔para changes are only present in the H and H⁺ tables (reactive
collisions); the H₂ and He tables conserve the nuclear spin. The 3:1 ratio is NOT imposed.

Statistical equilibrium: dn_i/dt = Σ_j n_j R_ji − n_i Σ_j R_ij = 0, Σ n_i = 1, solved with the GTH algorithm (state
reduction: no subtractions, so the nearly decoupled ortho and para chains of the low-density limit keep full relative
accuracy; equivalent to the linear solve with one row replaced by the normalisation). Radiative rates: u→l
A_ul β + B_ul J̄, l→u B_lu J̄ (β = escape probability per line, J̄ = line-averaged mean intensity [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹]
per line, B_ul = c²A/2hν³, B_lu = g_u B_ul/g_l); defaults β = 1, J̄ = 0. If the ortho and para chains are not coupled at
all (no H or H⁺ density), each chain is normalised to the fraction `opr/(1+opr)` (default 3:1); `opr = 3.0` forces a
fixed ortho:para ratio at any density (reactive transitions dropped; used for the comparison with the 3:1 fits of GA08).
For the H₂ collider the density is split into ortho and para using the computed populations (fixed-point iteration,
converged to 1e-9, starts at 3:1).
"""
module H2Levels

using Printf
import ..H2Lines as HL

const h = HL.h; const c = HL.c; const kB = HL.kB
const DIR = joinpath(get(ENV, "POPIII_LIT", joinpath(homedir(), "Projects", "popiii_lit")), "H2coll", "cloudy", "c25.00_extract", "data", "h2")

"One collider's de-excitation table on the 63 levels (pairs u > l in energy)."
struct CollTable
    name::String
    Tg::Vector{Float64}          # temperature grid [K]
    lT::Vector{Float64}          # log T grid
    pu::Vector{Int}; pl::Vector{Int}
    lk::Matrix{Float64}          # (npairs, nT): log k(u→l)
    nrows::Int; nswap::Int; ndup::Int
end

struct Levels
    E::Vector{Float64}; g::Vector{Float64}; v::Vector{Int}; J::Vector{Int}
    ortho::BitVector
    iu::Vector{Int}; il::Vector{Int}; A::Vector{Float64}; ν::Vector{Float64}     # radiative lines
    Bul::Vector{Float64}; Blu::Vector{Float64}
    colls::Dict{Symbol,CollTable}                                                  # :H, :H2o, :H2p, :He, :Hp
end

function read_coll(path, name, idx, E)
    ls = readlines(path)
    Tg = Float64[]; rows = Tuple{Int,Int,Vector{Float64}}[]; nrows = 0; nswap = 0
    pairs = Dict{Tuple{Int,Int},Vector{Float64}}(); ndup = 0
    for l in Iterators.drop(ls, 1)                       # line 1: magic number
        s = strip(l); (isempty(s) || startswith(s, "#")) && continue
        tok = split(s)
        if isempty(Tg); Tg = [parse(Float64, replace(t, r"[dD]" => "e")) for t in tok]; continue; end
        length(tok) == 4 + length(Tg) || error("$name: row with $(length(tok)) fields (expected $(4 + length(Tg))): $s")
        nrows += 1
        vu, Ju, vl, Jl = parse.(Int, tok[1:4])
        haskey(idx, (vu, Ju)) && haskey(idx, (vl, Jl)) || continue
        u = idx[(vu, Ju)]; lo = idx[(vl, Jl)]
        if E[u] < E[lo]; u, lo = lo, u; nswap += 1; end
        u == lo && continue
        k = [parse(Float64, replace(t, r"[dD]" => "e")) for t in tok[5:end]]
        haskey(pairs, (u, lo)) && (ndup += 1)
        pairs[(u, lo)] = k
    end
    ks = sort!(collect(keys(pairs)))
    lk = zeros(length(ks), length(Tg))
    for (p, key) in enumerate(ks); lk[p, :] .= log.(max.(pairs[key], 1e-300)); end
    CollTable(name, Tg, log.(Tg), first.(ks), last.(ks), lk, nrows, nswap, ndup)
end

"""
    load(; he = :ORNL, dir_lines, dir_coll) -> Levels

The 63 levels v ≤ 2, J ≤ 20 (sorted by energy), the 231 radiative lines of `H2Lines.load()` and the collision tables.
`he = :LeBourlot` selects the Le Bourlot et al. 1999 He–H₂ table (9 T points, v ≤ 3, J ≤ 16) instead of the ORNL one.
"""
function load(; he = :ORNL, dir_lines = joinpath(get(ENV, "POPIII_LIT", joinpath(homedir(), "Projects", "popiii_lit")), "H2lines"), dir_coll = DIR)
    HLL = HL.load(dir_lines; vmax = 2, Jmax = 20)
    # v, J of the ExoMol states (H2Lines.LineList does not keep them)
    vJ = Dict{Int,Tuple{Int,Int}}()
    for l in readlines(joinpath(dir_lines, "1H2__RACPPK.states"))
        s = split(l); isempty(s) && continue
        vJ[parse(Int, s[1])] = (parse(Int, s[5]), parse(Int, s[4]))
    end
    keep = sort!([k for (k, (v, J)) in vJ if v <= 2 && J <= 20]; by = k -> HLL.E[k])
    new = Dict(k => i for (i, k) in enumerate(keep)); n = length(keep)
    E = [HLL.E[k] for k in keep]; g = [HLL.g[k] for k in keep]
    v = [vJ[k][1] for k in keep]; J = [vJ[k][2] for k in keep]
    idx = Dict((v[i], J[i]) => i for i in 1:n)
    iu = [new[k] for k in HLL.iu]; il = [new[k] for k in HLL.il]
    Bul = [c^2 * HLL.A[i] / (2h * HLL.ν[i]^3) for i in eachindex(HLL.A)]
    Blu = [g[iu[i]] / g[il[i]] * Bul[i] for i in eachindex(HLL.A)]
    hef = he === :ORNL ? "coll_rates_He_ORNL.dat" : "coll_rates_He_LeBourlot.dat"
    colls = Dict{Symbol,CollTable}(
        :H   => read_coll(joinpath(dir_coll, "coll_rates_H_07.dat"), "H", idx, E),
        :H2p => read_coll(joinpath(dir_coll, "coll_rates_H2para_LeBourlot.dat"), "H2(para)", idx, E),
        :H2o => read_coll(joinpath(dir_coll, "coll_rates_H2ortho_LeBourlot.dat"), "H2(ortho)", idx, E),
        :He  => read_coll(joinpath(dir_coll, hef), "He", idx, E),
        :Hp  => read_coll(joinpath(dir_coll, "coll_rates_Hp.dat"), "H+", idx, E))
    Levels(E, g, v, J, BitVector(J .% 2 .== 1), iu, il, HLL.A, HLL.ν, Bul, Blu, colls)
end

"Print the coverage of every collision table on the level set: transitions found, of the n(n−1)/2 pairs, and J/v reach."
function coverage(L::Levels; io = stdout)
    n = length(L.E); npair = n * (n - 1) ÷ 2
    for key in (:H, :H2p, :H2o, :He, :Hp)
        t = L.colls[key]; np = length(t.pu)
        conn = falses(n); conn[t.pu] .= true; conn[t.pl] .= true
        reactive = count(p -> L.ortho[t.pu[p]] != L.ortho[t.pl[p]], 1:np)
        nrad = count(i -> any(p -> t.pu[p] == L.iu[i] && t.pl[p] == L.il[i], 1:np), eachindex(L.A))
        @printf(io, "  %-10s %4d of %d pairs (rows read %d, swapped %d, duplicates %d); levels linked %d/%d; o<->p pairs %d; radiative lines with a collision partner %d/%d; T %g-%g K (%d pts); max J_u %d, max v_u %d\n",
                t.name, np, npair, t.nrows, t.nswap, t.ndup, count(conn), n, reactive, nrad, length(L.A), t.Tg[1], t.Tg[end], length(t.Tg),
                np == 0 ? 0 : maximum(L.J[t.pu]), np == 0 ? 0 : maximum(L.v[t.pu]))
    end
end

"Interpolated de-excitation coefficients of a table at T (linear in log k vs log T, end values held outside the grid)."
function kde(t::CollTable, T)
    x = clamp(log(T), t.lT[1], t.lT[end]); nT = length(t.lT)
    i = clamp(searchsortedlast(t.lT, x), 1, nT - 1); w = (x - t.lT[i]) / (t.lT[i+1] - t.lT[i])
    @views exp.((1 - w) .* t.lk[:, i] .+ w .* t.lk[:, i+1])
end

"Collisional rate matrix Q[i, j] (rate i→j [s⁻¹]) for the colliders with densities `dens` (name => n [cm⁻³], H2 handled by the caller)."
function collrates!(Q, L::Levels, T, dens::AbstractVector{<:Pair{Symbol,<:Real}}; reactive = true)
    for (key, n) in dens
        n > 0 || continue
        t = L.colls[key]; k = kde(t, T)
        for p in eachindex(k)
            u = t.pu[p]; lo = t.pl[p]
            (!reactive && L.ortho[u] != L.ortho[lo]) && continue
            Q[u, lo] += n * k[p]
            Q[lo, u] += n * k[p] * L.g[u] / L.g[lo] * exp(-(L.E[u] - L.E[lo]) / (kB * T))
        end
    end
    Q
end

const COLLKEYS = (:H, :He, :Hp, :H2o, :H2p)          # the collider order of `populations`
"""
Work arrays of `populations` for one task (`ws = PopWS(L)`): the rate matrix, the copy that `gth` reduces, and per collider
the de-excitation coefficients and Boltzmann factors at the zone's temperature, which the ortho:para iteration reuses
(they do not depend on the ortho fraction). The arithmetic and its order are those of the path without `ws`: the
populations are bit-identical.
"""
struct PopWS{T}
    Q::Matrix{T}; P::Matrix{T}; s::Vector{T}; π::Vector{T}
    Qo::Matrix{T}; Po::Matrix{T}; so::Vector{T}; πo::Vector{T}
    Qp::Matrix{T}; Pp::Matrix{T}; sp::Vector{T}; πp::Vector{T}
    k::Vector{Vector{T}}; ex::Vector{Vector{T}}; have::Vector{Bool}
    o::Vector{Int}; p::Vector{Int}
    κ0::Vector{T}; S::Vector{T}; β::Vector{T}; J::Vector{T}
end
function PopWS{T}(L::Levels) where {T}
    n = length(L.E); o = findall(L.ortho); p = findall(.!L.ortho); no = length(o); np = length(p); nl = length(L.A)
    PopWS{T}(zeros(T, n, n), zeros(T, n, n), zeros(T, n), zeros(T, n),
             zeros(T, no, no), zeros(T, no, no), zeros(T, no), zeros(T, no),
             zeros(T, np, np), zeros(T, np, np), zeros(T, np), zeros(T, np),
             [zeros(T, length(L.colls[c].pu)) for c in COLLKEYS],
             [zeros(T, length(L.colls[c].pu)) for c in COLLKEYS], fill(false, length(COLLKEYS)), o, p,
             zeros(T, nl), zeros(T, nl), zeros(T, nl), zeros(T, nl))
end
PopWS(L::Levels) = PopWS{Float64}(L)
const POPWS_POOL = [IdDict{Levels,Dict{DataType,Any}}() for _ in 1:Threads.maxthreadid()]
__init__() = (resize!(POPWS_POOL, Threads.maxthreadid()); for i in eachindex(POPWS_POOL); isassigned(POPWS_POOL, i) || (POPWS_POOL[i] = IdDict{Levels,Dict{DataType,Any}}()); end)   # (when a precompiled TitanStack image is loaded: the build process has one thread)
function popws(L::Levels, ::Type{T} = Float64) where {T}
    d0 = POPWS_POOL[Threads.threadid()]
    d = get!(d0, L) do; Dict{DataType,Any}(); end
    w = get(d, T, nothing)
    w === nothing ? (w2 = PopWS{T}(L); d[T] = w2; w2) : (w::PopWS{T})
end

"`kde` into k (the same interpolation, element by element)."
function kde!(k::Vector{T}, t::CollTable, Tval) where {T}
    x = clamp(log(Tval), t.lT[1], t.lT[end]); nT = length(t.lT)
    i = clamp(searchsortedlast(t.lT, x), 1, nT - 1); w = (x - t.lT[i]) / (t.lT[i+1] - t.lT[i])
    @inbounds for p in eachindex(k); k[p] = exp((1 - w) * t.lk[p, i] + w * t.lk[p, i+1]); end
    k
end

"`collrates!` with the densities in the order of COLLKEYS and the coefficients of `ws` (filled on first use at this T)."
function collrates!(Q::Matrix{T}, L::Levels, Tval, ns::NTuple{5,T}, ws::PopWS{T}; reactive = true) where {T}
    for (c, key) in enumerate(COLLKEYS)
        n = ns[c]; n > 0 || continue
        t = L.colls[key]; k = ws.k[c]; ex = ws.ex[c]
        if !ws.have[c]
            kde!(k, t, Tval)
            @inbounds for p in eachindex(k); ex[p] = exp(-(L.E[t.pu[p]] - L.E[t.pl[p]]) / (kB * Tval)); end
            ws.have[c] = true
        end
        @inbounds for p in eachindex(k)
            u = t.pu[p]; lo = t.pl[p]
            (!reactive && L.ortho[u] != L.ortho[lo]) && continue
            Q[u, lo] += n * k[p]
            Q[lo, u] += n * k[p] * L.g[u] / L.g[lo] * ex[p]
        end
    end
    Q
end

"Stationary distribution of the Markov generator with off-diagonal rates Q[i,j] (i→j), by GTH state reduction."
gth(Q0::AbstractMatrix) = gth!(Matrix{eltype(Q0)}(Q0), zeros(eltype(Q0), size(Q0, 1)))
"`gth` on work arrays: P is overwritten with Q0 and reduced, s is scratch (same arithmetic as `gth(Q0)`)."
gth(Q0::AbstractMatrix, P::Matrix{T}, s::Vector{T}) where {T<:Real} = gth!(copyto!(P, Q0), fill!(s, zero(T)))
gth(Q0::AbstractMatrix, P::Matrix{T}, s::Vector{T}, π::Vector{T}) where {T<:Real} = gth!(copyto!(P, Q0), fill!(s, zero(T)), π)
function gth!(P::Matrix{T}, s::Vector{T}, π::Vector{T} = zeros(T, size(P, 1))) where {T<:Real}
    n = size(P, 1)
    for i in 1:n; P[i, i] = 0.0; end
    for k in n:-1:2
        sk = 0.0; for j in 1:k-1; sk += P[k, j]; end
        sk > 0 || error("H2Levels.gth: state $k has no route to lower states (reducible chain)")
        s[k] = sk
        for i in 1:k-1
            pik = P[i, k]; pik == 0 && continue
            for j in 1:k-1; P[i, j] += pik * P[k, j] / sk; end
        end
    end
    fill!(π, zero(T)); π[1] = 1.0
    for k in 2:n
        a = zero(T); for i in 1:k-1; a += π[i] * P[i, k]; end
        π[k] = a / s[k]
    end
    sp = sum(π)
    @inbounds for i in 1:n; π[i] = π[i] / sp; end
    π
end

"Radiative rates added to Q: u→l A β + B_ul J̄, l→u B_lu J̄."
function radrates!(Q, L::Levels, β, Jbar)
    for i in eachindex(L.A)
        u = L.iu[i]; lo = L.il[i]
        b = β === nothing ? 1.0 : (β isa Number ? Float64(β) : β[i])
        J = Jbar === nothing ? 0.0 : (Jbar isa Number ? Float64(Jbar) : Jbar[i])
        Q[u, lo] += L.A[i] * b + L.Bul[i] * J
        Q[lo, u] += L.Blu[i] * J
    end
    Q
end

function solve_Q(L::Levels, Q::Matrix{T}, opr, ws::Union{Nothing,PopWS{T}} = nothing) where {T}
    n = length(L.E); o = ws === nothing ? findall(L.ortho) : ws.o; p = ws === nothing ? findall(.!L.ortho) : ws.p
    cross = false
    @inbounds for j in p, i in o
        if Q[i, j] > 0 || Q[j, i] > 0; cross = true; break; end
    end
    if cross && opr === nothing
        return ws === nothing ? gth(Q) : gth(Q, ws.P, ws.s, ws.π)
    end
    xo = opr === nothing ? 0.75 : opr / (1 + opr)
    if ws === nothing
        π = zeros(T, n)
        π[o] .= xo .* gth(Q[o, o]); π[p] .= (1 - xo) .* gth(Q[p, p])
        return π
    end
    @inbounds for (jj, j) in enumerate(o), (ii, i) in enumerate(o); ws.Qo[ii, jj] = Q[i, j]; end
    @inbounds for (jj, j) in enumerate(p), (ii, i) in enumerate(p); ws.Qp[ii, jj] = Q[i, j]; end
    go = gth(ws.Qo, ws.Po, ws.so, ws.πo); gp = gth(ws.Qp, ws.Pp, ws.sp, ws.πp)
    π = fill!(ws.π, zero(T))
    @inbounds for (ii, i) in enumerate(o); π[i] = xo * go[ii]; end
    @inbounds for (ii, i) in enumerate(p); π[i] = (1 - xo) * gp[ii]; end
    π
end

"""
    populations(L, T, n_H, n_H2, n_He, n_Hp; β = nothing, Jbar = nothing, opr = nothing) -> Vector (63, Σ = 1)

Statistical-equilibrium level populations (fractions of the H₂ molecules) at gas temperature T [K] with collider
densities in cm⁻³ (no electron collisions). β: escape probability (scalar or per line, default 1), Jbar: mean intensity
per line [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹] (scalar or per line, default 0). `opr`: fixed ortho:para ratio (conversion rates dropped), see
the module docs.
"""
@inline function _solve_for!(w::PopWS{ET}, L::Levels, T, n_H, n_H2, n_He, n_Hp, xo, β, Jbar, opr) where {ET}
    Q = fill!(w.Q, zero(ET))
    collrates!(Q, L, ET(T), (ET(n_H), ET(n_He), ET(n_Hp), ET(n_H2 * xo), ET(n_H2 * (1 - xo))), w; reactive = opr === nothing)
    radrates!(Q, L, β, Jbar)
    solve_Q(L, Q, opr, w)
end
function populations(L::Levels, T, n_H, n_H2, n_He, n_Hp; β = nothing, Jbar = nothing, opr = nothing, ws::Union{Nothing,PopWS} = nothing)
    ET = promote_type(typeof(T), typeof(n_H), typeof(n_H2), typeof(n_He), typeof(n_Hp), typeof(opr === nothing ? 0.75 : opr / (1 + opr)),
                      β === nothing ? Float64 : eltype(β), Jbar === nothing ? Float64 : eltype(Jbar))
    w = (ws isa PopWS{ET}) ? ws : popws(L, ET)
    fill!(w.have, false)
    if opr !== nothing
        π = _solve_for!(w, L, T, n_H, n_H2, n_He, n_Hp, opr / (1 + opr), β, Jbar, opr)
        return ws === nothing ? copy(π) : π
    end
    if !(n_H2 > 0)
        π = _solve_for!(w, L, T, n_H, n_H2, n_He, n_Hp, 0.75, β, Jbar, opr)
        return ws === nothing ? copy(π) : π
    end
    xo = 0.75 + zero(ET); π = _solve_for!(w, L, T, n_H, n_H2, n_He, n_Hp, xo, β, Jbar, opr)
    for it in 1:50
        @inbounds for (ii, i) in enumerate(w.o); w.πo[ii] = π[i]; end
        xn = sum(w.πo)
        abs(xn - xo) < 1e-9 && break
        xo = xn; π = _solve_for!(w, L, T, n_H, n_H2, n_He, n_Hp, xo, β, Jbar, opr)
    end
    ws === nothing ? copy(π) : π
end

"Line cooling Σ n_u A β hν per volume [erg cm⁻³ s⁻¹] for H₂ number density n_H2 (pops: fractions)."
function cooling(L::Levels, pops, n_H2; β = nothing)
    s = 0.0
    for i in eachindex(L.A)
        b = β === nothing ? 1.0 : (β isa Number ? Float64(β) : β[i])
        s += pops[L.iu[i]] * L.A[i] * b * h * L.ν[i]
    end
    n_H2 * s
end

"""
    collcool(L, T, pops, n_H, n_H2, n_He, n_Hp) -> erg cm⁻³ s⁻¹

Net energy the gas loses to the H₂ levels by collisions, n_H2 Σ_{i≠j} n_i C_ij (E_j − E_i), at the populations `pops`
(fractions; the ortho share taken from them for the H₂ collider, as `populations`). This is the gas's line cooling
(negative: heating). In statistical equilibrium it equals the net radiative loss Σ hν [n_u A − (n_l B_lu − n_u B_ul) J̄]
of the same populations, however the levels are pumped: radiative pumping and fluorescence move no gas energy.
"""
function collcool(L::Levels, T, pops, n_H, n_H2, n_He, n_Hp)
    n_H2 > 0 || return 0.0
    n = length(L.E); Q = zeros(n, n); xo = sum(pops[L.ortho])
    collrates!(Q, L, T, [:H => Float64(n_H), :He => Float64(n_He), :Hp => Float64(n_Hp), :H2o => n_H2 * xo, :H2p => n_H2 * (1 - xo)])
    s = 0.0
    @inbounds for j in 1:n, i in 1:n
        i == j && continue
        q = Q[i, j]; q == 0 && continue
        s += pops[i] * q * (L.E[j] - L.E[i])
    end
    n_H2 * s
end

"""
    netrad(L, pops, n_H2, β, J) -> erg cm⁻³ s⁻¹

Net radiative loss n_H2 Σ hν [n_u (A β + B_ul J̄) − n_l B_lu J̄] of populations `pops` under the radiative rates of `radrates!`
(β, J̄ per line). For populations in statistical equilibrium under those rates (`populations(...; β, Jbar = J)`) it equals
`collcool`, the gas's net collisional loss (energy conservation of the level system), at O(lines) cost.
"""
function netrad(L::Levels, pops, n_H2, β, J)
    s = 0.0
    @inbounds for i in eachindex(L.A)
        u = L.iu[i]; lo = L.il[i]
        s += h * L.ν[i] * (pops[u] * (L.A[i] * β[i] + L.Bul[i] * J[i]) - pops[lo] * L.Blu[i] * J[i])
    end
    n_H2 * s
end

"""
    linekappa(L, pops, n_H2) -> κ0 per line

Line-centre opacity coefficient hν/(4π) n_H2 (n_l B_lu − n_u B_ul) [cm⁻¹ Hz, the profile factor excluded] of populations `pops` (as H2Lines.cooling_nlte's opac).
"""
linekappa(L::Levels, P, n_H2) = [h * L.ν[i] / (4π) * n_H2 * (P[L.il[i]] * L.Blu[i] - P[L.iu[i]] * L.Bul[i]) for i in eachindex(L.A)]
function linekappa!(κ0::Vector{T}, L::Levels, P::Vector{T}, n_H2) where {T}
    @inbounds for i in eachindex(L.A)
        κ0[i] = h * L.ν[i] / (4π) * n_H2 * (P[L.il[i]] * L.Blu[i] - P[L.iu[i]] * L.Bul[i])
    end
    κ0
end

"Line source function n_u A/(n_l B_lu − n_u B_ul) of populations `pops` (B_ν(T) where the line is inverted), per line."
linesource(L::Levels, P, T) = [(den = P[L.il[i]] * L.Blu[i] - P[L.iu[i]] * L.Bul[i]; den > 0 ? P[L.iu[i]] * L.A[i] / den : HL.Bν(L.ν[i], T)) for i in eachindex(L.A)]
function linesource!(S::Vector{T}, L::Levels, P::Vector{T}, Tval) where {T}
    @inbounds for i in eachindex(L.A)
        den = P[L.il[i]] * L.Blu[i] - P[L.iu[i]] * L.Bul[i]
        S[i] = den > 0 ? P[L.iu[i]] * L.A[i] / den : HL.Bν(L.ν[i], Tval)
    end
    S
end

"""
Allocation-free LINEIMP cooling kernel (`PopIIIMaterial.li_lam`) evaluated inside `H2Levels` so that `populations`, `linekappa!`,
`linesource!`, and `netrad` use concrete per-thread work buffers (`PopWS{ET}`) for both `Float64` and `ForwardDiff.Dual`.
The populations are re-solved at the iterate with the radiative rates (β, J̄) of the last (ALI) transfer call.
"""
function li_lam_eval(LV::Levels, β0::Vector{Float64}, J0::Vector{Float64}, T, nH, nH2, nHe, nHp)
    Tm = max(T, 1.0)
    ET = promote_type(typeof(Tm), typeof(nH), typeof(nH2), typeof(nHe), typeof(nHp))
    ws = popws(LV, ET)
    P = populations(LV, Tm, nH, nH2, nHe, nHp; β = β0, Jbar = J0, ws)
    netrad(LV, P, nH2, β0, J0)
end

"LTE (Boltzmann) fractions over the 63 levels."
function lte(L::Levels, T)
    w = L.g .* exp.(-(L.E .- L.E[1]) ./ (kB * T)); w ./ sum(w)
end

end # module
