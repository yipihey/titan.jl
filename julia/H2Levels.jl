isdefined(Main, :H2Lines) || include(joinpath(@__DIR__, "H2Lines.jl"))
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
import Main.H2Lines as HL

const h = HL.h; const c = HL.c; const kB = HL.kB
const DIR = joinpath(homedir(), "Projects", "popiii_lit", "H2coll", "cloudy", "c25.00_extract", "data", "h2")

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
function load(; he = :ORNL, dir_lines = joinpath(homedir(), "Projects", "popiii_lit", "H2lines"), dir_coll = DIR)
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
function collrates!(Q, L::Levels, T, dens::Vector{Pair{Symbol,Float64}}; reactive = true)
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

"Stationary distribution of the Markov generator with off-diagonal rates Q[i,j] (i→j), by GTH state reduction."
function gth(Q0::AbstractMatrix)
    n = size(Q0, 1); P = Matrix{Float64}(Q0); s = zeros(n)
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
    π = zeros(n); π[1] = 1.0
    for k in 2:n
        a = 0.0; for i in 1:k-1; a += π[i] * P[i, k]; end
        π[k] = a / s[k]
    end
    π ./ sum(π)
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

function solve_Q(L::Levels, Q, opr)
    n = length(L.E); o = findall(L.ortho); p = findall(.!L.ortho)
    cross = (sum(@view Q[o, p]) + sum(@view Q[p, o])) > 0
    if cross && opr === nothing
        return gth(Q)
    end
    xo = opr === nothing ? 0.75 : opr / (1 + opr)
    π = zeros(n)
    π[o] .= xo .* gth(Q[o, o]); π[p] .= (1 - xo) .* gth(Q[p, p])
    π
end

"""
    populations(L, T, n_H, n_H2, n_He, n_Hp; β = nothing, Jbar = nothing, opr = nothing) -> Vector (63, Σ = 1)

Statistical-equilibrium level populations (fractions of the H₂ molecules) at gas temperature T [K] with collider
densities in cm⁻³ (no electron collisions). β: escape probability (scalar or per line, default 1), Jbar: mean intensity
per line [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹] (scalar or per line, default 0). `opr`: fixed ortho:para ratio (conversion rates dropped), see
the module docs.
"""
function populations(L::Levels, T, n_H, n_H2, n_He, n_Hp; β = nothing, Jbar = nothing, opr = nothing)
    n = length(L.E)
    function solve_for(xo)
        Q = zeros(n, n)
        collrates!(Q, L, T, [:H => Float64(n_H), :He => Float64(n_He), :Hp => Float64(n_Hp),
                             :H2o => n_H2 * xo, :H2p => n_H2 * (1 - xo)]; reactive = opr === nothing)
        radrates!(Q, L, β, Jbar)
        solve_Q(L, Q, opr)
    end
    opr !== nothing && return solve_for(opr / (1 + opr))
    n_H2 > 0 || return solve_for(0.75)
    xo = 0.75; π = solve_for(xo)
    for it in 1:50
        xn = sum(π[L.ortho])
        abs(xn - xo) < 1e-9 && break
        xo = xn; π = solve_for(xo)
    end
    π
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

"LTE (Boltzmann) fractions over the 63 levels."
function lte(L::Levels, T)
    w = L.g .* exp.(-(L.E .- L.E[1]) ./ (kB * T)); w ./ sum(w)
end

end # module
