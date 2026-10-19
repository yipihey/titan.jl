"""
    SCVH

Saumon, Chabrier & Van Horn (1995, ApJS 99, 713) equation-of-state tables for H and He (the interpolated "i" tables
h_tab_i.dat, he_tab_i.dat; copies from YREC, see POPIII_NOTES M2e). Each table: 63 isotherms, log T = 2.10 … 7.06 in
steps of 0.08; per isotherm a header `logT npts`, then rows
  log P, X_H2, X_H, log ρ, log S, log U, (∂log ρ/∂log T)_P, (∂log ρ/∂log P)_T, (∂log S/∂log T)_P, (∂log S/∂log P)_T, ∇_ad
(cgs; X = number concentrations of H₂ and H for the H table, and of He and He⁺ for the He table). The column order
is confirmed by an independent labelled copy: github.com/chreinhardt/scvh, eos-tables-pt/hydrogen_scvh_extended.data,
header "log T, logP, x_H2, x_H, logRho, logE, logS". At log T = 2.10 and log P = 4.0 its logE = 10.1201 and
logS = 8.8910 are this table's columns 6 and 5. That copy extends to log T = 1.06 and low P with Kovetz et al.
(2009); not used here.

Mixture (SCVH95 additive-volume rule): at fixed (P, T), 1/ρ = X/ρ_H + Y/ρ_He, U = X U_H + Y U_He (the entropy of
mixing is ignored here: it does not enter p or U).

`state(ρ, T; X)` → (; p, u, xH2, xH) at given ρ, T: Newton on log P (bilinear interpolation of log ρ and log U in
(log T, log P) on each table).
"""
module SCVH

struct Table
    logT::Vector{Float64}
    logP::Vector{Vector{Float64}}
    cols::Vector{Matrix{Float64}}        # per isotherm: npts × 10 (columns 2–11)
end

function read_table(path)
    lines = readlines(path); i = 1
    logT = Float64[]; logP = Vector{Float64}[]; cols = Matrix{Float64}[]
    while i <= length(lines)
        h = split(lines[i]); i += 1
        isempty(h) && continue
        t = parse(Float64, h[1]); n = parse(Int, h[2])
        P = zeros(n); C = zeros(n, 10)
        for j in 1:n
            v = parse.(Float64, split(lines[i])); i += 1
            P[j] = v[1]; C[j, :] .= v[2:11]
        end
        push!(logT, t); push!(logP, P); push!(cols, C)
    end
    Table(logT, logP, cols)
end

"Column c (2 = X_H2, 3 = X_H, 4 = log ρ, 5 = log S, 6 = log U, …) on isotherm j at log P (linear; clamped)."
function oniso(tb::Table, j, lP, c)
    P = tb.logP[j]; C = tb.cols[j]
    lP <= P[1] && return C[1, c-1]
    lP >= P[end] && return C[end, c-1]
    k = searchsortedfirst(P, lP)
    w = (lP - P[k-1]) / (P[k] - P[k-1])
    (1 - w) * C[k-1, c-1] + w * C[k, c-1]
end
"Column c at (log T, log P): linear in log T between the bracketing isotherms."
function at(tb::Table, lT, lP, c)
    j = clamp(searchsortedfirst(tb.logT, lT), 2, length(tb.logT))
    w = clamp((lT - tb.logT[j-1]) / (tb.logT[j] - tb.logT[j-1]), 0.0, 1.0)
    (1 - w) * oniso(tb, j - 1, lP, c) + w * oniso(tb, j, lP, c)
end

const TABS = Ref{Any}(nothing)
const DIR = Ref(joinpath(homedir(), "Projects", "popiii_lit", "SCVH95"))
tables() = (TABS[] === nothing && (TABS[] = (read_table(joinpath(DIR[], "h_tab_i.dat")),
                                             read_table(joinpath(DIR[], "he_tab_i.dat")))); TABS[])

"Mixture at (log T, log P): (log ρ, U [erg/g], X_H2, X_H)."
function mix(lT, lP; X = 0.75)
    H, He = tables(); Y = 1 - X
    ρH = 10^at(H, lT, lP, 4); ρHe = 10^at(He, lT, lP, 4)
    ρ = 1 / (X / ρH + Y / ρHe)
    U = X * 10^at(H, lT, lP, 6) + Y * 10^at(He, lT, lP, 6)
    (log10(ρ), U, at(H, lT, lP, 2), at(H, lT, lP, 3))
end

"(ρ, T) → (; p, u, xH2, xH): bisection on log P (log ρ increases monotonically with log P on every isotherm) so that
the mixture density equals ρ; the bracket is the tables' log P range."
function state(ρ, T; X = 0.75)
    lT = log10(T); lρ = log10(ρ); H, _ = tables()
    lo, hi = minimum(first.(H.logP)), maximum(last.(H.logP))
    f(x) = mix(lT, x; X)[1] - lρ
    (f(lo) > 0 || f(hi) < 0) && error("SCVH.state: (ρ, T) = ($ρ, $T) outside the table")
    for _ in 1:80
        m = (lo + hi) / 2; f(m) > 0 ? (hi = m) : (lo = m)
        hi - lo < 1e-12 && break
    end
    x = (lo + hi) / 2
    _, U, xH2, xH = mix(lT, x; X)
    (; p = 10^x, u = U, xH2, xH)
end

end # module
