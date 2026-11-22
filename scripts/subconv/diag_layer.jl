# M1 b-d diagnosis tables of the sub-shock layer from RUNREC snapshots with cd.  julia diag_layer.jl <snaps.jls> <step> [nbuf=2]
#   prints (1) interface table (gates, N2, nabla, opacity), (2) zone energy-row terms as specific-entropy rates, (3) time scales.
#   julia diag_layer.jl <snaps.jls> series   prints the layer scalars at every snapshot.
using Serialization, Printf
const MH = 1.6726e-24; const KB = 1.380649e-16; const YR = 3.156e7; const RS = 6.957e10; const MS = 1.989e33
function load(f); v = Any[]; open(f) do io; while !eof(io); push!(v, deserialize(io)); end; end; v; end
shockif(s) = findfirst(<(-2e5), s.u)
# unstable run just inside the shock: interfaces with Mach<0.5 and (N2_L<0 or N2_eq<0), contiguous inwards from the outermost found within 4 interfaces of the shock
function layer(s)
    cd = s.cd; js = shockif(s); js === nothing && return (js, 0, 0)
    bad(j) = cd.mach_i[j] < 0.5 && (cd.N2_L[j] < 0 || cd.N2_eq[j] < 0)
    jo = 0; for j in js-1:-1:max(2, js-4); bad(j) && (jo = j; break); end
    jo == 0 && return (js, 0, 0)
    ji = jo; while ji - 1 >= 2 && bad(ji - 1); ji -= 1; end
    (js, ji, jo)
end
function itable(s, ji, jo, nb)
    cd = s.cd; js = shockif(s); R = s.r[js]; n = length(s.ρ)
    @printf("# step %d t=%.6e yr  zones %d  R*=%.4g Rsun (shock interface %d)  layer interfaces %d..%d  layer[]=%s\n", s.step, s.t/YR, n, R/RS, js, ji, jo, join((@sprintf("%.4g", x) for x in cd.layer), ","))
    @printf("%4s %6s %8s %9s %6s %3s %10s %10s %10s %7s %7s %7s %7s %7s %5s %5s %5s %9s %9s %9s %9s %8s %11s\n", "j", "r/R*", "T", "rho", "Mach", "br", "N2_eq", "N2_S", "N2_L", "nab", "nab_rad", "nab_L", "nad_fr", "nad_eq", "w_tau", "w_inf", "w_th", "kap", "Fpot/F", "Fconv/F", "dr/HP", "dlnT/dr", "L_rad[Lsun]")
    for j in max(2, ji - nb):min(n, jo + nb + 1)
        T = 0.5 * (s.T[j-1] + s.T[j]); ρ = 0.5 * (s.ρ[j-1] + s.ρ[j]); dr = 0.5 * (s.r[min(j + 1, n + 1)] - s.r[j-1])
        @printf("%4d %6.4f %8.1f %9.3e %6.3f %3.0f %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.4f %7.4f %5.2f %5.2f %5.0f %9.3e %9.3e %9.3e %9.3e %8.3f %10.3e%s\n", j, s.r[j]/R, T, ρ, cd.mach_i[j], cd.branch[j], cd.N2_eq[j], cd.N2_S[j], cd.N2_L[j],
                cd.nab[j], cd.nab_rad[j], cd.nab_L[j], cd.nad_fr[j], cd.nad_eq[j], cd.w_tau[j], cd.w_inf[j], cd.w_th[j], cd.kap[j], cd.Fpot_over_Ftot[j], cd.Fconv_over_Ftot[j], dr/cd.HP[j],
                (log(s.T[j]) - log(s.T[j-1])) / (0.5 * (s.r[j+1] - s.r[j-1])) * R, 4pi * s.r[j]^2 * cd.frad[j] / 3.828e33, j == js ? " <-shock" : (ji <= j <= jo ? " *" : ""))
    end
end
# zone k lies between interfaces k and k+1; contributions c_i = -e_i/(m T) in k_B per baryon per yr (heating positive); net = e_stor+e_pdv = m T ds/dt
function etable(s, ji, jo, nb)
    cd = s.cd; n = length(s.ρ); R = s.r[shockif(s)]
    names = (:e_av, :e_kin, :e_adv, :e_dfe, :e_rad, :e_conv, :e_chem, :e_src)
    @printf("# energy-row terms as entropy rates [k_B/baryon/yr] = -e_i mu mH/(kB T mz) yr  (heating>0); net = (e_stor+e_pdv) same units = ds/dt; mz = rho dV/sr\n")
    @printf("%4s %6s %8s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %8s\n", "k", "r/R*", "T", "rho", "mu", "stor+pdv", "sum_terms", "av", "kin", "adv", "dfe", "rad", "conv", "chem", "src", "res/sum", "|rad|/|tot|", "t_rad[yr]", "t_settle", "dr/HP")
    for k in max(1, ji - 1 - nb):min(n, jo + nb)
        mz = s.ρ[k] * (s.r[k+1]^3 - s.r[k]^3) / 3; f = cd.mu[k] * MH / (KB * s.T[k] * mz) * YR
        c = [-getfield(cd, nm)[k] * f for nm in names]; net = (cd.e_stor[k] + cd.e_pdv[k]) * f
        tot = sum(abs, c)
        @printf("%4d %6.4f %8.1f %9.3e %9.4f %9.2e %9.2e %s %9.2e %9.2e %9.2e %9.2e %8.3f%s\n", k, 0.5 * (s.r[k] + s.r[k+1]) / R, s.T[k], s.ρ[k], cd.mu[k], net, sum(c), join((@sprintf("%9.2e", x) for x in c), " "), cd.e_res[k] * f / (abs(sum(c)) + 1e-300), abs(c[5]) / tot, s.e[k] * mz / abs(cd.e_rad[k] + 1e-300) / YR, 0.5 * (cd.t_settle[k] + cd.t_settle[k+1]) / YR,
                (s.r[k+1] - s.r[k]) / cd.HP[k], (ji - 1 <= k <= jo) ? " *" : "")
    end
end
function ttable(s, ji, jo, nb)
    cd = s.cd; n = length(s.ρ); arad = 7.5657e-15; cc = 2.9979e10
    @printf("# time scales: t_settle=H_P/|u|, t_N=1/sqrt|N2_L| [yr], t_chem [SECONDS] (coupled zone j); potential MLT, all gates open, same formula as fpot_th:\n")
    @printf("#   fr = frozen (c_P,fr, delta_fr, nad_fr; U with c_P,fr): sup_fr = nab-nad_fr, zeta, v_pot, t_turn=lmix/v, Fpot/Ftot;\n")
    @printf("#   eq = equilibrium chemistry (delta_eq, nad_eq; U->0 limit zeta=sqrt(sup_eq); c_P,fr used for the flux = LOWER BOUND since c_P,eq >= c_P,fr): sup_eq, v_eq, t_turn_eq, Feq_lb/Ftot\n")
    @printf("%4s %9s %9s %9s %8s %8s %9s %9s %9s %9s %9s %9s %9s %9s\n", "j", "t_settle", "t_N", "t_chem[s]", "U", "sup_fr", "v_pot_fr", "t_turn_fr", "Fpot/Ft", "sup_eq", "v_eq", "t_turn_eq", "Feq_lb/Ft", "tau")
    for j in max(2, ji - nb):min(n, jo + nb)
        ρ = 0.5 * (s.ρ[j-1] + s.ρ[j]); T = 0.5 * (s.T[j-1] + s.T[j]); open_ = cd.w_th[j] == 1.0
        cp = open_ ? cd.cp_used[j] : cd.cp_fr[j]; dl = open_ ? cd.del_used[j] : cd.del_fr[j]; na = open_ ? cd.nad_used[j] : cd.nad_fr[j]
        sup = cd.nab[j] - na; U = 3arad * cc * T^3 / (cp * ρ^2 * cd.kap[j] * cd.lmix[j]^2) * sqrt(8cd.HP[j] / (cd.g[j] * dl))
        ζ = sup > 0 ? sqrt(U^2 + sup) - U : 0.0; v = cd.lmix[j] * sqrt(cd.g[j] * dl / (8cd.HP[j])) * ζ
        se = cd.nab[j] - cd.nad_eq[j]; ze = se > 0 ? sqrt(se) : 0.0; ve = cd.lmix[j] * sqrt(cd.g[j] * cd.del_eq[j] / (8cd.HP[j])) * ze
        Fe = se > 0 ? ρ * cd.cp_fr[j] * T * sqrt(cd.g[j] * cd.del_eq[j]) * cd.lmix[j]^2 * ze^3 / (4sqrt(2.0) * cd.HP[j] * sqrt(cd.HP[j])) : 0.0
        @printf("%4d %9.2e %9.2e %9.2e %8.1e %8.3f %9.2e %9.2e %9.2e %9.3f %9.2e %9.2e %9.2e %9.2e\n", j, cd.t_settle[j]/YR, cd.t_N[j]/YR, cd.t_chem[j], U, sup, v, v > 0 ? cd.lmix[j]/v/YR : NaN, cd.Fpot_over_Ftot[j],
                se, ve, ve > 0 ? cd.lmix[j]/ve/YR : NaN, Fe / cd.Ftot[j], cd.tau[j])
    end
end
function series(v)
    @printf("%6s %11s %3s %3s %3s %8s %8s %8s %9s %9s %9s %9s\n", "step", "t[yr]", "ji", "jo", "n", "r_in/R*", "r_out/R*", "m_lay", "max-N2td2", "Fconv/F", "Fpot/F", "R*[Rsun]")
    for s in v
        js, ji, jo = layer(s); cd = s.cd
        ji == 0 ? @printf("%6d %11.5e  none (layer[]=%s)\n", s.step, s.t/YR, join((@sprintf("%.3g", x) for x in cd.layer), ",")) :
        @printf("%6d %11.5e %3d %3d %3d %8.4f %8.4f %8.4f %9.3e %9.3e %9.3e %9.2f\n", s.step, s.t/YR, ji, jo, jo-ji+1, s.r[ji]/s.r[js], s.r[jo]/s.r[js], cd.layer[3], cd.layer[4], cd.layer[5], cd.layer[6], s.r[js]/RS)
    end
end
# per unstable interface j (zones j-1 inside, j outside): Delta s = s_j - s_{j-1} [k_B/baryon] from N2_L (ds/dr = c_P N2/(g delta)), and the rate of change of Delta s from each energy-row term
# (rate_i = c_i(zone j) - c_i(zone j-1)); tau_s = |Delta s| / |net rate|. Negative Delta s = unstable (entropy decreasing outwards).
function summary(v)
    names = (:e_av, :e_kin, :e_adv, :e_dfe, :e_rad, :e_conv, :e_chem)
    hd = ("step", "j", "r/R*", "N2_L", "N2_eq", "Ds_L", "Ds_eq", "dDs/dt", "av", "kin", "adv", "dfe", "rad", "conv", "chem", "tauS_L/tN", "nabrad/nab", "kap_j", "kap_max", "kapmax_j")
    @printf("# Ds = s_j - s_(j-1) [k_B/baryon] from N2 (ds/dr = c_P N2/(g_eff delta)); Ds_eq uses c_P,fr (LOWER BOUND |Ds_eq|); dDs/dt and its terms = rate_i(zone j) - rate_i(zone j-1) [k_B/baryon/yr]; tauS_L = |Ds_L|/|dDs/dt|\n")
    println(join(hd, " "))
    for s in v
        js, ji, jo = layer(s); ji == 0 && continue; cd = s.cd; R = s.r[js]; n = length(s.ρ)
        rate(k) = (mz = s.ρ[k] * (s.r[k+1]^3 - s.r[k]^3) / 3; f = cd.mu[k] * MH / (KB * s.T[k] * mz) * YR; [-getfield(cd, nm)[k] * f for nm in names])
        for j in ji:jo
            dr = 0.5 * (s.r[j+1] - s.r[j-1]); fm = 0.5 * (cd.mu[j-1] + cd.mu[j]) * MH / KB
            dsL = cd.cp_fr[j] * cd.N2_L[j] * dr / (cd.g_eff[j] * cd.del_fr[j]) * fm; dsE = cd.cp_fr[j] * cd.N2_eq[j] * dr / (cd.g_eff[j] * cd.del_eq[j]) * fm
            dc = rate(j) - rate(j-1); net = sum(dc); lo = max(2, ji - 3); km = argmax(cd.kap[lo:min(n, jo + 3)]) + lo - 1
            @printf("%d %d %.4f %.2e %.2e %.2e %.2e %.2e | %s | %.2e %.2e %.2e %.2e %d\n", s.step, j, s.r[j]/R, cd.N2_L[j], cd.N2_eq[j], dsL, dsE, net, join((@sprintf("%.2e", x) for x in dc), " "),
                    abs(dsL) / abs(net) / (cd.t_N[j] / YR), cd.nab_rad[j] / cd.nab[j], cd.kap[j], cd.kap[km], km)
        end
    end
end
# gate statistics over all subsonic (Mach<0.5) interfaces inside the shock with N2_L<0 or N2_eq<0 (not only the contiguous layer): which gate blocks MLT
function gates(v)
    @printf("%5s %4s %5s %5s %6s %6s %6s %6s %6s %6s %9s %9s\n", "step", "n", "N2L<0", "N2eq<0", "w_th=0", "w_tau<1", "w_inf<1", "br=1", "Fconv>0", "Fpot>0", "maxFconv/F", "maxFpot/F")
    tot = zeros(Int, 9)
    for s in v
        cd = s.cd; js = shockif(s); js === nothing && continue
        c = zeros(Int, 9); mc = 0.0; mp = 0.0
        for j in 2:js-1
            cd.mach_i[j] < 0.5 && (cd.N2_L[j] < 0 || cd.N2_eq[j] < 0) || continue
            c .+= [1, cd.N2_L[j] < 0, cd.N2_eq[j] < 0, cd.w_th[j] == 0, cd.w_tau[j] < 0.99, cd.w_inf[j] < 0.99, cd.branch[j] == 1, cd.fconv[j] > 0, cd.fpot_th[j] > 0]
            mc = max(mc, cd.Fconv_over_Ftot[j]); mp = max(mp, cd.Fpot_over_Ftot[j])
        end
        tot .+= c
        @printf("%5d %4d %5d %5d %6d %6d %6d %6d %6d %6d %9.2e %9.2e\n", s.step, c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], c[9], mc, mp)
    end
    @printf("TOTAL over %d snapshots: unstable subsonic interfaces %d; N2_L<0 %d; N2_eq<0 %d; w_th=0 %d; w_tau<1 %d; w_inf<1 %d; branch 1 %d; Fconv>0 %d; fpot_th>0 %d\n", length(v), tot...)
end
function main()
    v = load(ARGS[1])
    if ARGS[2] == "series"; series(v); return; end
    if ARGS[2] == "summary"; summary(v); return; end
    if ARGS[2] == "gates"; gates(v); return; end
    st = parse(Int, ARGS[2]); nb = length(ARGS) > 2 ? parse(Int, ARGS[3]) : 2
    s = v[findlast(x -> x.step == st, v)]; js, ji, jo = layer(s)
    if ji == 0; println("step $st: no unstable subsonic layer (js=$js); table around the shock instead"); ji = jo = js - 1; nb = 4; end
    itable(s, ji, jo, nb); println(); etable(s, ji, jo, nb); println(); ttable(s, ji, jo, nb)
    println("# species: ", s.spn)
end
main()
