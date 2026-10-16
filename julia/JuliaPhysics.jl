"""
    JuliaPhysics

Julia-native replacements for TITAN's eos.f and opac.f, written for clarity and
speed (precomputed constants, algebraic simplifications, Julia's own sqrt/^),
NOT for bit parity.  They operate in place on the mapped COMMON arrays over
k = ngrs-1 .. ngre+2, exactly the range the Fortran routines cover.

Implemented: leos = 2 (ideal gas), lopac = 3 (constant χρ), 4 (Thomson),
6 (Kramers + Thomson).  Anything else falls back to the Fortran routine.
"""
module JuliaPhysics

import ..TitanDriver
const TD = TitanDriver
const TF = TD.TF

# cached zero-copy views of the COMMON arrays we touch
const V = Dict{Symbol,Vector{Float64}}()   # typed: keeps the kernels type-stable
v(s::Symbol) = get!(() -> TF.var(s), V, s)
gf(s::Symbol) = Float64(TD.g(s))::Float64
range_k() = (Int(TD.g(:ngrs))::Int - 1):(Int(TD.g(:ngre))::Int + 2)

"Ideal-gas EOS (leos = 2) in Julia; falls back to Fortran otherwise."
function eos!()
    Int(TD.g(:leos))::Int == 2 || return TD.fortran_eos()
    gam, gmu = gf(:gam), gf(:gmu)
    Rg  = gf(:ck) / (gf(:cm0) * gmu)          # gas constant per unit mass
    cv  = Rg / (gam - 1)
    dels = (gam - 1) / gam
    cp  = Rg / dels
    nef = 0.5 / (gf(:cm0) * gmu)
    gR  = gam * Rg
    dn, tn = v(:dn), v(:tn)
    pgn, egn, xnen, delsn, cpn, betan, gamman, asn =
        v(:pgn), v(:egn), v(:xnen), v(:delsn), v(:cpn), v(:betan), v(:gamman), v(:asn)
    d = (v(:dlpgdltn), v(:dlpgdldn), v(:dlegdltn), v(:dlegdldn), v(:dlpedltn), v(:dlpedldn),
         v(:dldsdltn), v(:dldsdldn), v(:dlcpdltn), v(:dlcpdldn), v(:dlbedltn), v(:dlbedldn))
    dval = (1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0)
    @inbounds for k in range_k()
        ρ, T = dn[k], tn[k]
        pgn[k] = Rg * ρ * T
        egn[k] = cv * T
        xnen[k] = nef * ρ
        delsn[k] = dels; cpn[k] = cp; gamman[k] = gam
        betan[k] = inv(ρ * T)
        asn[k] = sqrt(gR * T)                  # = sqrt(γ p/ρ), simplified
        for j in 1:12; d[j][k] = dval[j]; end
    end
    nothing
end

"Opacity (lopac = 3, 4, 6) in Julia; falls back to Fortran otherwise."
function opac!()
    lop = Int(TD.g(:lopac))::Int
    lop in (3, 4, 6) || return TD.fortran_opac()
    dn, tn = v(:dn), v(:tn)
    chifn, xken, xkpn, xnen = v(:chifn), v(:xken), v(:xkpn), v(:xnen)
    dcfd, dcft, dked, dket, dkpd, dkpt, dped, dpet =
        v(:dlcfdldn), v(:dlcfdltn), v(:dlkedldn), v(:dlkedltn),
        v(:dlkpdldn), v(:dlkpdltn), v(:dlpedldn), v(:dlpedltn)
    ratio, csige, cm0, gmu = gf(:ratio), gf(:csige), gf(:cm0), gf(:gmu)
    if lop == 3                                  # constant χρ = chif0
        chif0 = gf(:chif0)
        xnefac = (1 - ratio) * chif0 / csige     # (1-ratio)·χ·ρ/σ_T, ρ cancels
        @inbounds for k in range_k()
            χ = chif0 / dn[k]
            chifn[k] = χ; xken[k] = ratio * χ; xkpn[k] = ratio * χ; xnen[k] = xnefac
            dcfd[k] = -1.0; dcft[k] = 0.0; dked[k] = -1.0; dket[k] = 0.0
            dkpd[k] = -1.0; dkpt[k] = 0.0; dped[k] = 0.0; dpet[k] = 1.0
        end
    elseif lop == 4                              # Thomson
        χ = csige / cm0
        @inbounds for k in range_k()
            chifn[k] = χ; xken[k] = ratio * χ; xkpn[k] = ratio * χ
            xnen[k] = (1 - ratio) * χ * dn[k] / csige
            dcfd[k] = 0.0; dcft[k] = 0.0; dked[k] = 0.0; dket[k] = 0.0
            dkpd[k] = 0.0; dkpt[k] = 0.0; dped[k] = 1.0; dpet[k] = 1.0
        end
    else                                         # lop == 6: Kramers ff + Thomson
        χT = csige / cm0
        nef = 0.5 / (gmu * cm0)
        @inbounds for k in range_k()
            T = tn[k]
            xkff = 8.6e23 * dn[k] / (T^3 * sqrt(T))      # T^-3.5 without pow()
            χ = xkff + χT
            chifn[k] = χ; xken[k] = xkff; xkpn[k] = xkff; xnen[k] = nef * dn[k]
            dcft[k] = -3.5 * xkff / χ; dcfd[k] = xkff / χ
            dket[k] = -3.5; dked[k] = 1.0; dkpt[k] = -3.5; dkpd[k] = 1.0
            dpet[k] = 1.0; dped[k] = 1.0
        end
    end
    nothing
end

end # module
