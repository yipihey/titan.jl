# Reference solution + accuracy metric for the SOD shock tube regression (popiii_regress.jl).
# Planar 1D Euler, gamma 1.4, domain [1,2], reflecting walls, Sod IC at 1.5 (rho,p = 1,1 | 0.125,0.1, u = 0);
# MUSCL-Hancock (minmod) + HLLC with N = 4000 cells, forward-Euler update, CFL 0.4: a line-by-line port of the
# sodref.py used for the 9.37e-3 -> 8.20e-3 comparison (Nov 2026). Metric = sodcmp.py's: L1 error of the TITAN cell
# values against the reference cell-averaged over each TITAN cell, relative: sum|v-ref|dx / sum|ref|dx.
module SodReference
export sod_reference, sod_l1
const G = 1.4

@inline function hllc(rL, uL, pL, rR, uR, pR)
    EL = pL/(G-1) + 0.5rL*uL^2; ER = pR/(G-1) + 0.5rR*uR^2
    cL = sqrt(G*pL/rL); cR = sqrt(G*pR/rR)
    SL = min(uL-cL, uR-cR); SR = max(uL+cL, uR+cR)
    Ss = (pR-pL + rL*uL*(SL-uL) - rR*uR*(SR-uR)) / (rL*(SL-uL) - rR*(SR-uR))
    FL = (rL*uL, rL*uL^2+pL, (EL+pL)*uL); FR = (rR*uR, rR*uR^2+pR, (ER+pR)*uR)
    UL = (rL, rL*uL, EL); UR = (rR, rR*uR, ER)
    SL >= 0 && return FL
    SR <= 0 && return FR
    if Ss >= 0
        f = rL*(SL-uL)/(SL-Ss); st = (f, f*Ss, f*(EL/rL + (Ss-uL)*(Ss + pL/(rL*(SL-uL)))))
        return ntuple(k -> FL[k] + SL*(st[k]-UL[k]), 3)
    else
        f = rR*(SR-uR)/(SR-Ss); st = (f, f*Ss, f*(ER/rR + (Ss-uR)*(Ss + pR/(rR*(SR-uR)))))
        return ntuple(k -> FR[k] + SR*(st[k]-UR[k]), 3)
    end
end
@inline minmod(a, b) = a*b > 0 ? sign(a)*min(abs(a), abs(b)) : zero(a)
@inline eulerflux(r, u, p) = (r*u, r*u^2+p, (p/(G-1)+0.5r*u^2+p)*u)
@inline cons(r, u, p) = (r, r*u, p/(G-1)+0.5r*u^2)
@inline prim(U) = (r = U[1]; u = U[2]/r; (r, u, (G-1)*(U[3]-0.5r*u^2)))

"Cell-centre primitives (xc, rho, u, p) at time T on N cells of [1,2] (cell averages)."
function sod_reference(N::Int = 4000, T::Float64 = 1.0)
    dx = 1.0/N; xc = 1.0 .+ ((1:N) .- 0.5) .* dx
    U = [cons(x < 1.5 ? 1.0 : 0.125, 0.0, x < 1.5 ? 1.0 : 0.1) for x in xc]
    Ng = N + 4; W = Vector{NTuple{3,Float64}}(undef, Ng); S = similar(W)
    WL = similar(W); WR = similar(W); Wp = similar(W); Wm = similar(W); fl = Vector{NTuple{3,Float64}}(undef, N+1)
    t = 0.0
    while t < T
        dt = Inf
        for i in 1:N
            r, u, p = prim(U[i]); dt = min(dt, 0.4dx/(abs(u) + sqrt(G*p/r)))
        end
        dt = min(dt, T - t)
        for i in 1:N; W[i+2] = prim(U[i]); end
        mir(w) = (w[1], -w[2], w[3])
        W[2] = mir(W[3]); W[1] = mir(W[4]); W[N+3] = mir(W[N+2]); W[N+4] = mir(W[N+1])
        S[1] = S[Ng] = (0.0, 0.0, 0.0)
        for i in 2:Ng-1; S[i] = ntuple(k -> minmod(W[i][k]-W[i-1][k], W[i+1][k]-W[i][k]), 3); end
        for i in 1:Ng
            wr = ntuple(k -> W[i][k] + 0.5S[i][k], 3); wl = ntuple(k -> W[i][k] - 0.5S[i][k], 3)
            Fr = eulerflux(wr...); Fl = eulerflux(wl...); h = 0.5dt/dx
            Ur = cons(wr...); Ul = cons(wl...)
            WL[i] = prim(ntuple(k -> Ur[k] - h*(Fr[k]-Fl[k]), 3))   # right-face value, evolved
            WR[i] = prim(ntuple(k -> Ul[k] - h*(Fr[k]-Fl[k]), 3))   # left-face value, evolved
        end
        for j in 1:N+1                                              # interface between Wg cells j+1, j+2 (1-based j+1..)
            a = WL[j+1]; b = WR[j+2]; fl[j] = hllc(a[1], a[2], a[3], b[1], b[2], b[3])
        end
        for i in 1:N; U[i] = ntuple(k -> U[i][k] - dt*(fl[i+1][k]-fl[i][k])/dx, 3); end
        t += dt
    end
    P = prim.(U)
    (xc = xc, rho = [q[1] for q in P], u = [q[2] for q in P], p = [q[3] for q in P], dx = dx)
end

"Relative L1 errors (rho, u, p) of TITAN cell values on interfaces `r` against reference `ref` averaged over each cell."
function sod_l1(r, rho, u_iface, p, ref)
    h = ref.dx; N = length(ref.xc); edges = vcat(1.0, ref.xc .+ h/2)
    cum(v) = vcat(0.0, cumsum(v)) .* h
    function interp(C, x)                                    # linear interpolation of the cumulative integral
        k = clamp(searchsortedlast(edges, x), 1, N); f = (x - edges[k]) / (edges[k+1] - edges[k])
        C[k] + f*(C[k+1] - C[k])
    end
    dx = diff(r); uc = 0.5 .* (u_iface[1:end-1] .+ u_iface[2:end])
    out = map(((v, W),) -> begin
        C = cum(W); rf = [(interp(C, r[i+1]) - interp(C, r[i])) / dx[i] for i in eachindex(dx)]
        sum(abs.(v .- rf) .* dx) / sum(abs.(rf) .* dx)
    end, ((rho, ref.rho), (uc, ref.u), (p, ref.p)))
    (rho = out[1], u = out[2], p = out[3])
end
end
