# E2E-A: hand-off of a star between the stages of the chained run (e2e_chain.jl).
#   export_star(path; mat, a, b)     (s6b_run.jl, at the end of the accretion): the state of the star on its own grid, zone by zone
#   star_from_handoff(path; ...)     (s5_run.jl, HANDOFF=path): a TITAN state of an isolated star (no accretion, free surface at P = 0)
#                                    with the same profiles, built exactly as `build_star` of s4_star.jl builds the polytrope
# The file holds: interface interior mass m [g], r, u (cm, cm/s), zone ρ, T, E_r, interface F_r, the full composition Y[19, zone] (mol/g, NucNet order) and scalars.
# Needs s4_common.jl (TD, TF, TX, PS, SM, ...) to be included first.
using Serialization

"Write the star in the TITAN state to `path` (zones a:b, interfaces a:b+1)."
function export_star(path; mat, a = Int(TF.getscalar(:ngrs)), b = Int(TF.getscalar(:ngre)), extra = (;))
    rn, dn, tn, un, ern, frn, xmen = (TF.var(s) for s in (:rn, :dn, :tn, :un, :ern, :frn, :xmen))
    xmtot = TF.getscalar(:xmtot)
    m = [(xmtot - xmen[k]) * 4π for k in a:b+1]
    h = (; format = "e2e-handoff-v1", t = TF.getscalar(:timen), m, r = copy(rn[a:b+1]), u = copy(un[a:b+1]), rho = copy(dn[a:b]), T = copy(tn[a:b]), Er = copy(ern[a:b]),
         Fr = copy(frn[a:b+1]), Y = copy(mat.Y[:, a:b]), alpha = TX.ALPHA_MLT[], extra...)
    serialize(path, h)
    h
end

"""
    star_from_handoff(path; alpha, efrac, dXmax, burnT, kcap, fcoupling, opac) -> (; mat, M, R, β, Tc, ρc, n, h)

Builds the TITAN state (lgrid = 3, regular centre, free outer boundary at pext = 0, StarGas with burning, MLT) from an exported star. The grid, the profiles and the
composition are those of the file; nothing is remeshed. The returned tuple has the fields `s5_run.jl` takes from `build_star` (β is NaN: no polytrope).
"""
function star_from_handoff(path; alpha = 1.5, efrac = 0.5, dXmax = 2e-3, burnT = :path, kcap = Inf, fcoupling = true, opac = :stellar, xmext_fac = 1e-3, screen = :chugunov)
    h = deserialize(path); n = length(h.rho)
    M = h.m[end]
    mat = SM.StarGas(TF.MGR; X = 0.75, Y = 0.25, Z = 0.0, nuc = true, opac, efrac, dXmax, burnT, screen, mixing = true, kcap, fcoupling)
    PS.build!(; rI = copy(h.r), rho = copy(h.rho), T = copy(h.T), u = copy(h.u), geom = 2, rad = true, Er = copy(h.Er), Fr = copy(h.Fr), lgrid = 3, inner = :wall, outer = :free,
              pext = 0.0, selfgrav = true, material = mat, convection = true, gam = 5 / 3, gmu = 1.0, xmext = xmext_fac * M / (4π))
    TX.REGULAR_CENTRE[] = true
    TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    TX.ALPHA_MLT[] = alpha
    PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    b - a + 1 == n || error("hand-off: $n zones in the file, $(b - a + 1) in the state")
    for j in 1:n; SM.setY!(mat, a + j - 1, view(h.Y, :, j)); end
    SM.ghosts!(mat)
    mat.Yprev .= mat.Y
    TD.EOS[](); TD.OPAC[]()
    TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    (; mat, M, R = h.r[end], β = NaN, Tc = h.T[1], ρc = h.rho[1], n, h)
end
