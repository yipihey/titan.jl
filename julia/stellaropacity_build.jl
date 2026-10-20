# Population III project, S3: build julia/stellar/kaptable_stellar.bin.gz (the StellarOpacity run-time table).
#   julia -t 4 --project=env stellaropacity_build.jl            (≈ 10–15 min)
# Inputs (data only, outside the repository; provenance in ~/Projects/popiii_lit/opal/README.md):
#   ~/Projects/popiii_lit/opal/ex/kap_input_data/opal/GN93hz      OPAL GN93hz (Iglesias & Rogers 1996), Z = 0 tables
#   ~/Projects/popiii_lit/opal/raw/condtabl.data                  Cassisi et al. 2007 electron conduction (MESA copy)
# Output table: log10 κ_R and log10 κ_P at (log T, log R, X nodes of OPAL) = 321 × 81 × 10, log T 1.5–9.5, log R −15…+5.
#
# Construction of log κ_R at a fine node (all weights are smoothsteps, hence the result is C¹ in the data):
#   OPf  = (1−s(T)) F_lin + s(T) F_anch,  s = smoothstep((log T − 5.5)/1)
#     F_lin  : OPAL (nonuniform cubic Hermite in log T, log R); holes at high log R filled by linear extrapolation with the
#              slope of the last interval (clamped to [0, 1.2]); below log R = −8 saturating-exponential extension;
#     F_anch : kappa_high + Δ, Δ = OPAL − kappa_high on OPAL's nodes, held constant beyond OPAL's last valid node of
#              each row (in log R) and beyond log T = 8.7 (so the high-T model carries the T and R shape, OPAL the level);
#   base = w_OP OPf + (1 − w_OP) LOW   (w_OP < 1 only near the cold-corner holes of OPAL, log T < 4.0, low R)
#   lkR  = w_lo base' ...  = w_lo LOW + (1 − w_lo) base,  w_lo = 1 − smoothstep((log T − 3.8)/0.3)
#   the table holds the radiative κ_R; conduction (Cassisi tables, embedded) is added at run time: 1/κ = 1/κ_rad + w_c/κ_cond.
# κ_P: w_lo LOW_P + (1 − w_lo) Z0-continuum Planck mean (general X) × degeneracy factor.
using Printf, Dates
const J = @__DIR__
include(joinpath(J, "PrimordialChem.jl")); include(joinpath(J, "ContinuumOpacity.jl")); include(joinpath(J, "Z0Physics.jl"))
include(joinpath(J, "StellarOpacity.jl")); include(joinpath(J, "StellarOpacityLow.jl"))
const SO = StellarOpacity; const SL = StellarOpacityLow
include(joinpath(J, "stellaropacity_data.jl"))
# ------------------------------------------------------------------ main
function build()
    t0 = time()
    tabs, lTo, lRo = read_gn93(OPALFILE); cond = read_cond(CONDFILE)
    nTo = length(lTo)
    lTf = [SO.LT0 + SO.DLT * (i - 1) for i in 1:SO.NT]; lRf = [SO.LR0 + SO.DLR * (j - 1) for j in 1:SO.NR]
    nX = length(SO.XNODES)
    ilowmax = findlast(x -> x <= 4.15 + 1e-9, lTf); izmin = findfirst(x -> x >= 3.8 - 1e-9, lTf)
    @printf("grid %d×%d×%d; LOW to log T %.3f (iT ≤ %d), Z0-continuum from log T %.3f (iT ≥ %d); threads %d\n",
            SO.NT, SO.NR, nX, lTf[ilowmax], ilowmax, lTf[izmin], izmin, Threads.nthreads())
    # --- own low-T means (LOW) and Z0-continuum Planck means ---------------------------------------------------------
    LOWR = fill(NaN, SO.NT, SO.NR, nX); LOWP = fill(NaN, SO.NT, SO.NR, nX)
    Threads.@threads for it in 1:ilowmax
        T = 10^lTf[it]; p = SL.prep(T)
        for ix in 1:nX, ir in 1:SO.NR
            ρ = 10^(lRf[ir] + 3 * (lTf[it] - 6)); X = SO.XNODES[ix]
            m = SL.means(p, ρ, SL.equilibrium(ρ, T, X))
            LOWR[it, ir, ix] = log10(max(m.κR, 1e-40)); LOWP[it, ir, ix] = log10(max(m.κP, 1e-40))
        end
    end
    @printf("LOW done %.0f s\n", time() - t0)
    PZ = fill(NaN, SO.NT, SO.NR, nX)
    jobs = [(it, ix) for it in izmin:SO.NT for ix in 1:nX]
    Threads.@threads for jb in jobs
        it, ix = jb; T = 10^lTf[it]; X = SO.XNODES[ix]
        for ir in 1:SO.NR
            ρ = 10^(lRf[ir] + 3 * (lTf[it] - 6))
            PZ[it, ir, ix] = log10(max(SL.z0_means(ρ, T, X).κP, 1e-40))
        end
    end
    @printf("Z0 continuum Planck done %.0f s\n", time() - t0)
    # --- OPAL, anchoring, blends -------------------------------------------------------------------------------------
    lkR = zeros(SO.NT, SO.NR, nX); lkP = zeros(SO.NT, SO.NR, nX); flag = zeros(UInt8, SO.NT, SO.NR, nX)
    rep = IOBuffer()
    stats = Dict{String,Any}()
    for (ix, X) in enumerate(SO.XNODES)
        A = copy(tabs[X])
        # OPAL's logT = 8.7, log R ≤ −7 entries are the same placeholders in every table (-0.668, -0.931, -1.004 vs the
        # pure-scattering value ≈ −1.017): treat as missing
        i87 = findfirst(x -> abs(x - 8.7) < 1e-9, lTo); for j in 1:3; A[i87, j] = NaN; end
        valid = .!isnan.(A)
        # F_lin: fill holes of every row by extrapolation
        Flin = copy(A); Hn = zeros(nTo, 19)
        for i in 1:nTo, j in 1:19; Hn[i, j] = log10(SO.kappa_high(10^(lRo[j] + 3 * (lTo[i] - 6)), 10^lTo[i], X)); end
        Dl = A .- Hn
        for i in 1:nTo
            v = findall(valid[i, :]); isempty(v) && continue
            j1, j2 = first(v), last(v)
            for j in 1:j1-1; Flin[i, j] = A[i, j1]; Dl[i, j] = Dl[i, j1]; end
            s = j2 > j1 ? clamp((A[i, j2] - A[i, j2-1]) / 0.5, 0.0, 1.2) : 0.0
            for j in j2+1:19; Flin[i, j] = A[i, j2] + s * 0.5 * (j - j2); Dl[i, j] = Dl[i, j2]; end
        end
        # the 8.7 row: the three masked entries take the first valid
        holes = [(i, j) for i in 1:nTo, j in 1:19 if !valid[i, j] && lTo[i] < 4.3]
        # evaluation, R outermost
        for (ir, lR) in enumerate(lRf)
            rows_lin = zeros(nTo); rows_del = zeros(nTo)
            lRc = clamp(lR, -8.0, 1.0); d = lR - lRc
            for i in 1:nTo
                rows_lin[i] = herm(lRo, Flin[i, :], lRc); rows_del[i] = herm(lRo, Dl[i, :], lRc)
                if d != 0
                    δ = d > 0 ? -0.01 : 0.01                                             # step inward from the edge
                    dlin = (herm(lRo, Flin[i, :], lRc + δ) - rows_lin[i]) / δ            # d/dlogR at the edge
                    ddel = (herm(lRo, Dl[i, :], lRc + δ) - rows_del[i]) / δ
                    φ = sign(d) * 0.5 * (1 - exp(-abs(d) / 0.5))                         # saturating (τ = 0.5)
                    rows_del[i] += ddel * φ
                    rows_lin[i] += d < 0 ? dlin * φ : clamp(dlin, 0.0, 1.2) * d          # linear in log R beyond +1
                end
            end
            for (it, lT) in enumerate(lTf)
                lin = herm(lTo, rows_lin, lT); del = herm(lTo, rows_del, lT)
                ρ = 10^(lR + 3 * (lT - 6)); Tt = 10^lT
                hi = log10(SO.kappa_high(ρ, Tt, X))
                s = smoothstep((lT - 5.5) / 1.0)
                OPf = (1 - s) * lin + s * (hi + del)
                wlo = 1 - smoothstep((lT - 3.8) / 0.3)
                wop = 1.0
                if lT <= 4.3 && !isempty(holes)
                    lRh = lRc
                    dmin = minimum(sqrt(((lT - lTo[i]) / 0.05)^2 + ((lRh - lRo[j]) / 0.5)^2) for (i, j) in holes)
                    wop = smoothstep(dmin / 2)
                end
                low = isnan(LOWR[it, ir, ix]) ? 0.0 : LOWR[it, ir, ix]
                base = wop * OPf + (1 - wop) * low
                lkR[it, ir, ix] = wlo * low + (1 - wlo) * base
                inbox = lT >= 3.75 - 1e-9 && lT <= 8.7 + 1e-9 && lR >= -8 - 1e-9 && lR <= 1 + 1e-9
                flag[it, ir, ix] = wlo >= 1 ? 4 : (wlo > 0 || wop < 1 ? 3 : (inbox ? 1 : 2))
                # Planck mean
                pz = PZ[it, ir, ix]
                lT >= 3.8 - 1e-9 && (pz -= SO.FF.p * log10(1 + smoothstep((lT - 6.0) / 1.0) * SO.FF.b * 2.7e11 * ρ / Tt^2))   # Pauli blocking of free–free, as in kappa_ff, only above 1e6 K (the 2.7e11ρ/T² form is a high-T fit)
                lowp = isnan(LOWP[it, ir, ix]) ? 0.0 : LOWP[it, ir, ix]
                lkP[it, ir, ix] = wlo >= 1 ? lowp : (wlo > 0 ? wlo * lowp + (1 - wlo) * pz : pz)
            end
        end
        # diagnostics at this X: OPAL vs LOW in the overlap, anchoring offset Δ on OPAL's last row
        @printf(rep, "X=%.2f  OPAL − LOW (log10 κ_R) on the fine grid where OPAL has data:\n", X)
        for lT in (3.80, 3.85, 3.90, 3.95, 4.00, 4.05, 4.10)
            it = round(Int, (lT - SO.LT0) / SO.DLT) + 1; ds = Float64[]
            io_ = findfirst(x -> abs(x - lT) < 1e-6, lTo)
            for (ir, lR) in enumerate(lRf)
                (lR < -8 || lR > 1 || abs(lR * 2 - round(lR * 2)) > 1e-9) && continue
                jo = round(Int, (lR + 8) / 0.5) + 1
                valid[io_, jo] || continue
                push!(ds, A[io_, jo] - LOWR[it, ir, ix])
            end
            isempty(ds) || @printf(rep, "   log T=%.2f: n=%2d  mean %+.3f  min %+.3f  max %+.3f\n", lT, length(ds), sum(ds) / length(ds), minimum(ds), maximum(ds))
        end
        i87 = findfirst(x -> abs(x - 8.7) < 1e-9, lTo)
        @printf(rep, "   high-T model vs OPAL: Δ = OPAL − kappa_high at log T = 8.7: %s\n",
                join([@sprintf("R=%.1f:%+.3f", lRo[j], Dl[i87, j]) for j in 4:2:19 if valid[i87, j]], " "))
        for lT in (7.0, 8.0, 8.5)
            i = findfirst(x -> abs(x - lT) < 1e-9, lTo)
            @printf(rep, "   Δ at log T = %.1f: %s\n", lT, join([@sprintf("R=%.1f:%+.3f", lRo[j], Dl[i, j]) for j in 1:2:19 if valid[i, j]], " "))
        end
    end
    @printf("OPAL blends done %.0f s\n", time() - t0)
    (; lkR, lkP, flag, cond, report = String(take!(rep)))
end

function save_table(t, path)
    buf = IOBuffer()
    write(buf, Int64(SO.NT), Int64(SO.NR), Int64(length(SO.XNODES)))
    write(buf, vec(t.lkR), vec(t.lkP)); write(buf, vec(t.flag)); write(buf, vec(t.cond[1.0]), vec(t.cond[2.0]))
    open(pipeline(`gzip -9 -c`, stdout = path), "w") do io; write(io, take!(buf)); end
    path
end

if abspath(PROGRAM_FILE) == @__FILE__
    t = build()
    println(t.report)
    @assert all(isfinite, t.lkR) && all(isfinite, t.lkP) "non-finite table entries"
    mkpath(joinpath(J, "stellar"))
    p = save_table(t, SO.TABLEFILE)
    @printf("wrote %s (%.2f MB)\n", p, filesize(p) / 1e6)
end
