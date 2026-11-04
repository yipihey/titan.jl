# M3 test: Iliev et al. 2006 test 1. n_H = 1e-3, T = 1e4 K fixed, 5e48 photons/s at 13.6 eV, static gas, no hydro. Analytic: r_I = r_S (1 - exp(-t/t_rec))^(1/3),
# with alpha_B = 2.59e-13: r_S = 5.4 kpc, t_rec = 122.4 Myr; the comparison's ~5 % scatter is the bar.  julia --project=env ionfront_test.jl [nzones] [dt_Myr]
using Printf, ForwardDiff
for f in ("PrimordialChem", "ContinuumOpacity", "ChemNet", "ChemRates", "IonRay"); include(joinpath(@__DIR__, f * ".jl")); end
const CN = ChemNet; const IR = IonRay
nz = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 128
dtMyr = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1.0
kpc = 3.0857e21; Myr = 3.15576e13
net, photos = ChemRates.modern()
env = (Trad = 0.0, J = ν -> 0.0, static = true, caseB = true)
nH = 1e-3; ρ = nH * CN.μH(); T = 1e4
# recombination coefficient the network really uses
ph = [(p, CN.photorates_t(p, T, env)...) for p in photos]
iH = findfirst(p -> p.R === :H, photos); αB = ph[iH][3]
envA = (Trad = 0.0, J = ν -> 0.0, static = true, caseB = false); αA = CN.photorates_t(photos[iH], T, envA)[2]
@printf("network recombination at 1e4 K: case B %.3e (Iliev 2.59e-13), case A %.3e (4.2e-13 standard); levels to n = %d\n", αB, αA, CN.NH_LEVELS)
Ndot = 5e48; L = 5.5e-3 * kpc
# Test 1 is monochromatic at 13.6 eV (SPEC=bb: the 1e5 K blackbody of Test 2, still with T fixed here)
if get(ENV, "SPEC", "mono") == "bb"
    G = IR.groups(); N0 = IR.blackbody(G, 1.0, 1e5); N0 .*= Ndot / sum(N0)
else
    G = IR.groups(; Elo = 13.6, Emax = 13.601, n = 1); N0 = [Ndot]
end
@printf("%d groups, photons above 13.6 eV normalised to %.2e; mean energy %.1f eV; photons in the last group %.1e of total\n", length(N0), sum(N0), sum(N0 .* G.Ec) / sum(N0) / CN.eV, N0[end] / sum(N0))
L = 6.6 * kpc; rI = collect(range(0.0, L, length = nz + 1))
x0 = zeros(length(net.free)); x0[findfirst(==(:Hp), net.free)] = 1.2e-3
xs = [copy(x0) for _ in 1:nz]
function front(xs, rI)                                    # 50 % neutral, linear between zone centres
    xp = [x[findfirst(==(:Hp), net.free)] for x in xs]; rc = 0.5 .* (rI[1:end-1] .+ rI[2:end])
    j = findfirst(<(0.5), xp); j === nothing && return rI[end]; j == 1 && return rc[1]
    rc[j-1] + (0.5 - xp[j-1]) / (xp[j] - xp[j-1]) * (rc[j] - rc[j-1])
end
rS = (3Ndot / (4π * αB * nH^2))^(1 / 3); trec = 1 / (αB * nH)
rS26 = (3Ndot / (4π * 2.59e-13 * nH^2))^(1 / 3); trec26 = 1 / (2.59e-13 * nH)
t = 0.0; dt = dtMyr * Myr; nbad = 0; tw = time()
outs = (10.0, 30.0, 100.0, 200.0, 500.0); io = 1
@printf("  t [Myr]  r_I [kpc]  analytic (own alpha)  ratio  analytic(2.59e-13)  ratio  photons: in vs ionized+recombined\n")
global t, io
let t = 0.0, io = 1, ion = 0.0, esc = 0.0
    while t < 500Myr * (1 - 1e-9)
        Nin = copy(N0); rec = 0.0
        for k in 1:nz
            src = IR.source(ρ, rI[k], rI[k+1], Nin, G, net)
            x1, _ = CN.step(xs[k], ρ, 0.0, T, dt, net; env, photos, extra = src.extra, accepted = src.accepted, Tfix = T)
            xp0 = xs[k][1]; xs[k] = x1; Nin = src.out(x1, dt)
            V = 4π / 3 * (rI[k+1]^3 - rI[k]^3)
            ion += (x1[findfirst(==(:Hp), net.free)] - xp0) * nH * V
            xp = x1[findfirst(==(:Hp), net.free)]
            rec += αB * (nH * xp)^2 * V * dt
        end
        ion += rec; t += dt; esc += sum(Nin) * dt
        if io <= length(outs) && t >= outs[io] * Myr * (1 - 1e-9)
            r = front(xs, rI); a = rS * (1 - exp(-t / trec))^(1 / 3); b = rS26 * (1 - exp(-t / trec26))^(1 / 3)
            @printf("  %6.0f   %7.3f    %7.3f              %.3f     %7.3f           %.3f  %.3f  (+ escaped %.3f)\n", t / Myr, r / kpc, a / kpc, r / a, b / kpc, r / b, ion / (Ndot * t), esc / (Ndot * t))
            io += 1
        end
    end
end
xp = [x[findfirst(==(:Hp), net.free)] for x in xs]
@printf("neutral fraction at 0.2, 0.5, 0.7, 0.9 r_S: %s (wall %.1f s)\n", join([@sprintf("%.2e", 1 - xp[clamp(searchsortedfirst(rI, f * rS) - 1, 1, nz)]) for f in (0.2, 0.5, 0.7, 0.9)], " "), time() - tw)
