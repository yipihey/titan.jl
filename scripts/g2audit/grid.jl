# G2 gap maps: table vs network-in-LTE on a (log ρ, log T) grid. Output: OUT/grid.bin (Float64, [field, iρ, iT]) + grid.txt (field names).
isdefined(Main, :net_lte) || (include(joinpath(@__DIR__, "load.jl")); include(joinpath(@__DIR__, "mats.jl")))
isdefined(Main, :PROTO) || (const PROTO = nothing)       # prototype opacity terms (proto.jl), or nothing
const OUT = get(ENV, "OUT", joinpath(@__DIR__, "out")); mkpath(OUT)
const LR = collect(-14.0:0.1:0.0); const LTg = collect(3.0:0.02:6.0)
const FIELDS = [:p_t, :u_t, :G1_t, :nab_t, :cP_t, :p_n, :u_n, :G1_n, :nab_n, :cP_n, :kR_t, :kP_t, :kR_n, :kP_n, :xe_t, :xe_n, :fH2_t, :fH2_n,
                :w_scvh, :kR_z, :kP_z, :kR_p, :kP_p, :xe_tfix, :uS, :pS, :xe_w]
# warm the lazy caches before threading
EQ.table(); SO.table(); COp.tables(); nlf(); ztab = Z0.load_table(joinpath(JD, "z0", isfile(joinpath(JD, "z0", "kaptable_Z0_ext.bin.gz")) ? "kaptable_Z0_ext.bin.gz" : "kaptable_Z0.bin.gz"))
net_means(1e-8, 5000.0, net_lte(1e-8, 5000.0))
const A = fill(NaN, length(FIELDS), length(LR), length(LTg))
fi(s) = findfirst(==(s), FIELDS)
t0 = time()
Threads.@threads for i in eachindex(LR)
    for j in eachindex(LTg)
        ρ = 10^LR[i]; T = 10^LTg[j]
        try
            tt = thermo(tab_pu, ρ, T); s = net_lte(ρ, T); tn = thermo(net_pu, ρ, T)
            κt = tab_kap(ρ, T); m = net_means(ρ, T, s)
            y = EQ.composition(ρ, T); z = EQ.composition_he(ρ, T)
            tb = EQ.table(); x = EQ.gx(tb, LR[i]); yy = EQ.gy(tb, LTg[j])
            ii = clamp(floor(Int, x), 1, tb.nρ - 1); jj = clamp(floor(Int, yy), 1, tb.nT - 1); a = x - ii; b = yy - jj
            w = (1 - a) * (1 - b) * tb.w[ii, jj] + a * (1 - b) * tb.w[ii+1, jj] + (1 - a) * b * tb.w[ii, jj+1] + a * b * tb.w[ii+1, jj+1]
            # Saha-ideal part of the table (PopIIIEqEOS.node without SCVH) and the SCVH x_e under the X_e reading (see doc)
            ys, z1, z2, _ = EQ.saha_he(ρ, T); uS, pS = EQ.eos_he(ρ, T, ys, z1, z2)
            xfix = NaN
            if w > 0
                sc = try SCVH.state(ρ, T) catch; nothing end
                if sc !== nothing
                    X2, X1 = sc.xH2, sc.xH; Xp = max(1 - X2 - X1, 0.0) / 2; xp2 = Xp / (2X2 + X1 + Xp)
                    xfix = (1 - w) * ys[1] + w * xp2 + z[1] + 2z[2]
                end
            else
                xfix = ys[1] + z[1] + 2z[2]
            end
            kz = (Z0.kap(ztab, ρ, T), Z0.kapP(ztab, ρ, T))
            kp = PROTO === nothing ? (NaN, NaN) : (mp = PROTO(ρ, T, s); (mp.κR, mp.κP))
            v = (tt.p, tt.u, tt.Γ1, tt.∇ad, tt.cP, tn.p, tn.u, tn.Γ1, tn.∇ad, tn.cP, κt[1], κt[2], m.κR, m.κP,
                 y[1] + z[1] + 2z[2], s.ne / s.NH, 2y[3], 2s.nH2 / s.NH, w, kz[1], kz[2], kp[1], kp[2], xfix, uS, pS, (T < 2500 ? NaN : try os = Z0.opacity_state(ρ, T; comp = Z0.Comp(XH, 1 - XH)); os.ne / os.nH catch; NaN end))
            for (k, val) in enumerate(v); A[k, i, j] = val; end
        catch e
            i == 1 && j == 1 && @warn "point failed" LR[i] LTg[j] e
        end
    end
end
@printf("grid %d × %d in %.0f s on %d threads\n", length(LR), length(LTg), time() - t0, Threads.nthreads())
open(joinpath(OUT, "grid.bin"), "w") do io; write(io, A); end
open(joinpath(OUT, "grid.txt"), "w") do io
    println(io, join(string.(FIELDS), " ")); println(io, join(string.(LR), " ")); println(io, join(string.(LTg), " "))
end
