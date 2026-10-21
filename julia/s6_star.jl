# S6a: initial models of the accreting protostar. `isentropic_star`: a fully convective (isentropic) gas+radiation star of mass M and central
# temperature Tc, integrated outward in mass with ∇ = ∇_ad (StellarEOS, composition X, Y, Z = 0) up to the point where its pressure falls to the ram
# pressure of the free-fall inflow, P_ram = Ṁ v_ff/(4π r²) (the contact is then in pressure balance); ρ_c is found for the required M.
using Printf
const GRAV = Z0Physics.G
"RHS of dy/dm for y = (r, ln ρ, ln T) on the isentrope."
function isentrope_rhs(m, y, comp)
    r, lρ, lT = y
    st = SE.eos(exp(lρ), exp(lT); comp, radiation = true)
    dlP = -GRAV * m / (4π * r^4) / st.P
    dlT = st.∇ad * dlP
    dlρ = (dlP - st.χT * dlT) / st.χρ
    (1 / (4π * r^2 * exp(lρ)), dlρ, dlT), st
end
"""
    isentropic_star(M, Tc, ρc; X = 0.75, Mdot = 0.0, dlm = 0.01) -> (; m, r, ρ, T, P, R, Mstop, stopped)

Integrates (RK4 in ln m) from m0 = 1e-9 M to m = M or until P ≤ P_ram(r, m) (Mdot > 0).
"""
function isentropic_star(M, Tc, ρc; X = 0.75, Mdot = 0.0, dlm = 0.01, comp = SE.composition(X, 1 - X, 0.0))
    m = 1e-9 * M; r = (3m / (4π * ρc))^(1 / 3); y = [r, log(ρc), log(Tc)]
    ms = [m]; ys = [copy(y)]; stopped = false
    while m < M * (1 - 1e-12)
        f0, st0 = isentrope_rhs(m, y, comp)
        h = min(m * (exp(dlm) - 1), M - m, 0.03 * st0.P / (GRAV * m / (4π * y[1]^4)))       # also |Δ ln P| ≤ 3 % per step
        f(mm, yy) = collect(isentrope_rhs(mm, yy, comp)[1])
        k1 = f(m, y); k2 = f(m + h / 2, y .+ h / 2 .* k1); k3 = f(m + h / 2, y .+ h / 2 .* k2); k4 = f(m + h, y .+ h .* k3)
        yn = y .+ h / 6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        m += h; y = yn; push!(ms, m); push!(ys, copy(y))
        if !(isfinite(y[1]) && y[1] > 0 && all(isfinite, y)) || y[3] < log(3500.0)
            return (; m = ms, r = [yy[1] for yy in ys], ρ = exp.([yy[2] for yy in ys]), T = exp.([yy[3] for yy in ys]), R = NaN, Mstop = m, stopped = true, bad = true)
        end
        if Mdot > 0
            st = SE.eos(exp(y[2]), exp(y[3]); comp, radiation = true)
            Pram = Mdot * sqrt(2GRAV * m / y[1]) / (4π * y[1]^2)
            if st.P <= Pram; stopped = true; break; end
        end
    end
    P = [SE.eos(exp(yy[2]), exp(yy[3]); comp, radiation = true).P for yy in ys]
    (; m = ms, r = [yy[1] for yy in ys], ρ = exp.([yy[2] for yy in ys]), T = exp.([yy[3] for yy in ys]), P, R = ys[end][1], Mstop = m, stopped, bad = false)
end
"ρc such that the isentrope of central temperature Tc contains mass M at the contact (bisection in ln ρc)."
function find_isentropic(M, Tc; X = 0.75, Mdot = 0.0, lo = 1e-6, hi = 1e3)
    f(lρ) = (s = isentropic_star(M, Tc, exp(lρ); X, Mdot); s.bad ? NaN : log(s.Mstop / M))
    a, b = log(lo), log(hi); fa = f(a)
    for _ in 1:60
        mid = 0.5(a + b); fm = f(mid)
        isnan(fm) && (fm = -1.0)
        (fm < 0) == (fa < 0) ? (a = mid; fa = fm) : (b = mid)
        abs(b - a) < 1e-12 && break
    end
    ρc = exp(0.5(a + b)); isentropic_star(M, Tc, ρc; X, Mdot)
end
