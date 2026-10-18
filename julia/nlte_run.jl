# Non-LTE hydrogen in the infall of a TITAN snapshot (NLTEH.jl), iterated with the continuum transfer (Spectrum.jl).
#   Z0_TABLE=ext TITAN_NOLIB=1 julia -t 2 --project=env nlte_run.jl NAME [maxit=40] [cscale=1] [caseB=0] [xray=1] [out=…]
# NAME: a spectrum from spectra_run.jl (lib/spectra/NAME.jls, with the free-fall extension); its atmosphere, including
# the Unsöld–Lucy temperatures below the shock, is used unchanged. Non-LTE points: the infall outside the shock zone
# (TITAN zones above it and the extension) with T < 2e4 K. Everything else stays in LTE.
# Each iteration: continuum J_ν at every point and frequency with the current populations, then the statistical
# equilibrium of every non-LTE point with those J_ν (Λ-iteration on the continua; the lines are treated by Sobolev
# escape probabilities, which are already "preconditioned"). Writes lib/spectra/nlte_OUT.jls.
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf
const SP = Spectrum; const Z0 = Z0Physics; const NH = NLTEH
name = ARGS[1]
opt = Dict(Symbol(split(a, "=")[1]) => split(a, "=")[2] for a in ARGS[2:end])
maxit = parse(Int, get(opt, :maxit, "40")); cscale = parse(Float64, get(opt, :cscale, "1"))
caseB = get(opt, :caseB, "0") == "1"; xray = get(opt, :xray, "1") == "1"
out = get(opt, :out, name * (cscale == 1 ? "" : "_c$(opt[:cscale])") * (caseB ? "_caseB" : "") * (xray ? "" : "_noX"))
d = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "$name.jls"))
S = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(d.run)_relax_snaps.jls"))
s = S[argmin([abs(x.t - d.t) for x in S])]
atm = d.atm; R0 = 8.4352e15; M = 1e6 * Z0.Msun
sel = findall(j -> (atm.zone[j] == 0 || atm.zone[j] > d.ksh) && atm.T[j] < 2e4, eachindex(atm.r))
vabs = [atm.zone[j] == 0 ? sqrt(2Z0.G * M / atm.r[j]) : abs(0.5 * (s.u[atm.zone[j]] + s.u[atm.zone[j]+1])) for j in sel]
dlnv = similar(vabs)
for (q, j) in enumerate(sel)
    a, b = max(q - 1, 1), min(q + 1, length(sel))
    dlnv[q] = a == b ? -0.5 : (log(vabs[b]) - log(vabs[a])) / (log(atm.r[sel[b]]) - log(atm.r[sel[a]]))
end
# shell volume of each non-LTE point (between the midpoints to its neighbours)
rs = atm.r[sel]
dV = [4π / 3 * ((q < length(rs) ? 0.5(rs[q] + rs[q+1]) : rs[q])^3 - (q > 1 ? 0.5(rs[q-1] + rs[q]) : rs[q])^3) for q in eachindex(rs)]
ν = SP.grid(; λmax_μm = 30.0, n = parse(Int, get(opt, :n, "260")), nedge = 12)
wν = [i == 1 ? 0.5(ν[2] - ν[1]) : i == length(ν) ? 0.5(ν[end] - ν[end-1]) : 0.5(ν[i+1] - ν[i-1]) for i in eachindex(ν)]
base = SP.Opac(atm)
# Unsöld–Lucy target and zones (as spectra_run.jl: from 0.95 R to zone kc − 2, TITAN's L_rad + L_conv)
dolucy = get(opt, :lucy, "1") == "1"
Lz = [0.5 * (s.Lrad[k] + s.Lconv[k] + s.Lrad[k+1] + s.Lconv[k+1]) for k in eachindex(s.d)]
Ltarget = [z > 0 ? Lz[z] : s.Lrad[end] for z in atm.zone]
lsel = findall(j -> 0 < atm.zone[j] <= d.kc - 2 && atm.r[j] >= 0.95R0, eachindex(atm.r))
_, R = SP.rays(atm.r; ncore = 12); R = SP.weights!(atm, R)
# start: LTE (Saha–Boltzmann) populations
nl = Any[nothing for _ in eachindex(atm.r)]
for j in sel
    st = base.s[j]; nl[j] = (n = copy(base.pH[j]), np = st.nHp, ne = st.ne)
end
lte = [(xH = base.s[j].xH, n2 = base.pH[j][2], ne = base.s[j].ne) for j in sel]
# init=OUT: start from an earlier solution lib/spectra/nlte_OUT.jls (same atmosphere)
if haskey(opt, :init)
    d0 = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$(opt[:init]).jls"))
    d0.sel == sel || error("init: different non-LTE points")
    for (q, j) in enumerate(sel); nl[j] = d0.nl[q]; end
    if haskey(d0, :Tstar)          # and from its corrected stellar temperatures
        global atm = merge(atm, (T = d0.Tstar,)); global base = SP.Opac(atm)
    end
    println("started from nlte_$(opt[:init])")
end
# Ng acceleration of the Λ-iteration: every NGEV iterations, extrapolate log(populations) from the last four iterates
NGEV = parse(Int, get(opt, :ng, "5")); minit = parse(Int, get(opt, :minit, "1"))
logpop() = reduce(vcat, [log.(vcat(nl[j].n, nl[j].np)) for j in sel])
function setpop!(v)
    for (q, j) in enumerate(sel)
        x = exp.(v[(q-1)*(NH.NL+1)+1:q*(NH.NL+1)])
        ne = x[end] + (nl[j].ne - nl[j].np)
        nl[j] = merge(nl[j], (n = x[1:NH.NL], np = x[end], ne = ne))
    end
end
nghist = Vector{Vector{Float64}}()
"Ng (order 2) extrapolation with equal weights (the iterates are logarithms)."
function ngl(h)
    x3, x2, x1, x0 = h
    d0 = x0 .- x1; a1 = d0 .- (x1 .- x2); a2 = d0 .- (x2 .- x3)
    A11 = sum(a1 .* a1); A12 = sum(a1 .* a2); A22 = sum(a2 .* a2); b1 = sum(d0 .* a1); b2 = sum(d0 .* a2)
    det = A11 * A22 - A12^2
    (abs(det) > 1e-14 * abs(A11 * A22) && isfinite(det)) || return nothing
    c1 = (b1 * A22 - b2 * A12) / det; c2 = (A11 * b2 - A12 * b1) / det
    (1 - c1 - c2) .* x0 .+ c1 .* x1 .+ c2 .* x2
end
println("$(length(sel)) non-LTE points, r/R = $(round(atm.r[sel[1]]/R0, digits = 5))–$(round(atm.r[sel[end]]/R0, digits = 1)); ",
        "$(length(ν)) frequencies; cscale = $cscale, caseB = $caseB, xray = $xray")
Jm = zeros(length(atm.r), length(ν)); L = zeros(length(ν)); Lline_last = Ref(0.0)
hist = []
t0 = time()
for it in 1:maxit
    op = NH.NOpac(base, nl)
    acc = zeros(length(atm.r), 6); lk = ReentrantLock()
    Threads.@threads for i in eachindex(ν)
        L[i], J, _, m = SP.monochromatic(atm, op, R, ν[i]; moments = true)
        Jm[:, i] .= J
        w = wν[i]
        lock(lk) do
            @. acc[:, 1] += w * m.J; @. acc[:, 2] += w * m.H; @. acc[:, 3] += w * m.κ * m.J
            @. acc[:, 4] += w * m.κ * m.B; @. acc[:, 5] += w * m.B; @. acc[:, 6] += w * m.χ * m.H
        end
    end
    # Unsöld–Lucy correction of the (LTE) star below the shock, as in spectra_run.jl, so that light the infall sends
    # back is re-emitted instead of being absorbed at fixed temperature (backwarming)
    if dolucy
        mom = (J = acc[:, 1], H = acc[:, 2], κJ = acc[:, 3] ./ acc[:, 1], κB = acc[:, 4] ./ acc[:, 5], B = acc[:, 5], χH = acc[:, 6] ./ acc[:, 2])
        Tn, lerr = SP.lucy(atm, mom, lsel, Ltarget; damp = 0.5)
        global atm = merge(atm, (T = Tn,)); global base = SP.Opac(atm)
        @printf("   Lucy: max |ΔH/H| %.4f in the star\n", lerr)
    end
    maxch = 0.0
    for (q, j) in enumerate(sel)
        T = atm.T[j]; st = base.s[j]; cur = nl[j]   # (non-LTE points keep their temperatures)
        # X-ray energy absorbed above 100 eV by H and He bound–free → secondary H ionisations (36 eV per ion pair).
        # Known artifact: credited to H even when n₁ → 0, which inverts Hα in the few points right above the shock.
        Γx = 0.0
        if xray
            dep = 0.0
            for i in eachindex(ν)
                E = Z0.h * ν[i] / Z0.eV; E < 100 && continue
                κx = cur.n[1] * NH.σbf(1, ν[i]) + st.nHe0 * SP.σHe1(E) + base.pHe[j][1] * SP.σHe2(E)
                dep += 4π * wν[i] * κx * Jm[j, i]
            end
            Γx = dep / (36Z0.eV) / max(cur.n[1], 1e-300)
        end
        n, np, ne, β, P = NH.solve_point(T, st.nH, st.ne - st.nHp, vabs[q], atm.r[j], dlnv[q], ν, wν, Jm[j, :], Γx,
                                      cur.n, cur.np; cscale, caseB)
        # electrons from helium: LTE value (neutral He at these temperatures)
        maxch = max(maxch, abs(np / cur.np - 1), abs(n[2] / cur.n[2] - 1))
        nl[j] = (n = n, np = np, ne = ne, Γx = Γx, β = β, P = P)
    end
    Ltot = sum(0.5 * (L[i] + L[i+1]) * (ν[i+1] - ν[i]) for i in 1:length(ν)-1)
    Lline = sum(nl[sel[q]].P * dV[q] for q in eachindex(sel)); Lline_last[] = Lline
    push!(nghist, logpop()); length(nghist) > 4 && popfirst!(nghist)
    if NGEV > 0 && length(nghist) == 4 && it % NGEV == 0
        v = ngl(nghist)
        if v !== nothing && all(isfinite, v)
            # limit the extrapolation to a factor 10 per population
            v = nghist[end] .+ clamp.(v .- nghist[end], -log(10), log(10)); setpop!(v); empty!(nghist)
            println("  Ng extrapolation")
        end
    end
    push!(hist, maxch)
    @printf("it %2d: max |Δn_p/n_p|, |Δn₂/n₂| = %.2e   ∫L_ν = %.4e erg/s, + lines %.3e   (%.0f s)\n", it, maxch, Ltot, Lline, time() - t0)
    flush(stdout)
    maxch < 1e-3 && it >= minit && break
end
# summary along the infall
ν32 = (NH.En[3] - NH.En[2]) / Z0.h
@printf("%8s %7s %9s %9s %9s %9s %9s %9s %9s %9s\n", "r/R", "T", "x_H LTE", "x_H NLTE", "b1", "b2", "b3", "S/B(Hα)", "S/J(Hα)", "Γx")
for (q, j) in enumerate(sel)
    (q % max(1, length(sel) ÷ 25) == 0 || q == 1) || continue
    T = atm.T[j]; st = nl[j]
    b = [st.n[i] / (st.np * st.ne * NH.Φ(i, T)) for i in 1:3]
    SL = 2Z0.h * ν32^3 / Z0.c^2 / (st.n[2] * NH.gn[3] / (st.n[3] * NH.gn[2]) - 1)
    i = searchsortedlast(ν, ν32); Jc = Jm[j, i] + (Jm[j, i+1] - Jm[j, i]) * (ν32 - ν[i]) / (ν[i+1] - ν[i])
    @printf("%8.4f %7.0f %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.2e\n", atm.r[j] / R0, T, lte[q].xH,
            st.np / base.s[j].nH, b..., SL / SP.Bν(ν32, T), SL / Jc, st.Γx)
end
serialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$out.jls"),
          (name = name, out = out, ν = ν, L = L, sel = sel, Tstar = atm.T, lucy = dolucy, nl = nl[sel], lte = lte, r = atm.r[sel], T = atm.T[sel],
           vabs = vabs, dlnv = dlnv, Jhα = [Jm[j, searchsortedlast(ν, ν32)] for j in sel], hist = hist,
           dV = dV, cscale = cscale, caseB = caseB, xray = xray))
println("saved lib/spectra/nlte_$out.jls")
# the same result in the format of spectra_run.jl (for export_report_data.jl): non-LTE ionisation in the infall
xH = [base.s[j].xH for j in eachindex(atm.r)]
for (q, j) in enumerate(sel); xH[j] = nl[j].np / base.s[j].nH; end
serialize(joinpath(@__DIR__, "..", "lib", "spectra", "$(name)_nlte$(out == name ? "" : "_" * out).jls"),
          (run = d.run, t = d.t, ν = ν, L = L, its = Int[], atm = atm, T_titan = d.T_titan, nlucy = -1, Mdot = d.Mdot,
           Ltitan = d.Ltitan, ksh = d.ksh, kc = d.kc, keep = d.keep, ion = false, xH = xH, Lline = Lline_last[]))
