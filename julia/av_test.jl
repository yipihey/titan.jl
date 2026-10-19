# G1 test: artificial-viscosity length ∝ local zone width (TitanResidFast.AV_ZONE) vs the old l = ql1·r̄.
# Spherical Noh problem on the Lagrangian grid with a regular centre (as popiii_m1.jl): cold uniform gas, ρ0, u = −v0,
# γ = 5/3, wall (reflecting) centre at r = 0, free outer boundary (p_ext = 0, so the outer shell keeps falling in at v0).
# Exact: shock speed v0/3, post-shock density 64 ρ0, pre-shock ρ = ρ0 (1 + v0 t/r)². Measured at the time the shock is
# at r_s = 0.1 R0: shock-position error (half-jump crossing of (ρ − ρ_pre)/(64ρ0 − ρ_pre)) and width = number of zones
# between the 10 % and 90 % crossings of u/u_pre ("wu"; the ρ-based 10–90 % width is ill-defined because the post-shock plateau is not flat). Grids: uniform in r, and Δr_{k+1} = s Δr_k (s = 1.02 … 1.2) beyond r = 0.3 R0.
#   TITAN_NOLIB=1 julia -t 4 --project=env av_test.jl                (runs the whole table, one process per case)
#   TITAN_NOLIB=1 julia --project=env av_test.jl one GEOM GRID S MODE VAL  (GEOM 0 planar | 2 spherical; GRID uniform|geom; MODE ql1|zone)
# Env: QL1S/QL1P (old-viscosity ql1, spherical/planar), AVC (default 2.0), THET (default 0.55), STOL (0.05).
using Printf, Statistics
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl"))
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup

const R0 = 1e12; const V0 = 1e7; const RHO0 = 1e-10; const GAM = 5 / 3; const GMU = 0.6
const XS_TARGET = 0.1R0                 # shock radius at the measurement time

function grid(kind, s; dx = 0.004R0, r_stretch = 0.3R0)
    # uniform zones dx up to r_stretch, then Δr_{k+1} = s Δr_k (stretch ratio s; kind "uniform" = s = 1) up to R0
    # (the shock reaches the stretched part when it has crossed half of the run; the last zone is trimmed; n even)
    kind == "uniform" && (s = 1.0)
    r = [0.0]; Δ = dx
    while r[end] < R0 - 1e-9R0
        r[end] >= r_stretch - 1e-9R0 && (Δ *= s)
        push!(r, min(r[end] + Δ, R0))
    end
    isodd(length(r) - 1) && deleteat!(r, length(r) - 1)
    r[end] = R0; r
end

"Run one case; returns a NamedTuple of measurements."
function run_case(geom, kind, s, mode, val; thet = parse(Float64, get(ENV, "THET", "0.55")))
    rmin = geom == 0 ? R0 : 0.0                  # planar: TITAN's unknowns are ln r, so the wall sits at r = R0 (x = r − R0)
    rI = rmin .+ grid(kind, s); n = length(rI) - 1
    Rg = 1.380658e-16 / (1.6605402e-24 * GMU); T0 = (V0 / 200)^2 / (GAM * Rg)
    PS.build!(; rI, rho = fill(RHO0, n), T = fill(T0, n), u = [k == 1 ? 0.0 : -V0 for k in 1:n+1], geom, rad = false,
              inner = :wall, outer = :free, pext = 0.0, gam = GAM, gmu = GMU, lgrid = 3, xscale = R0,
              av = (ql0 = 0.0, ql1 = mode == "ql1" ? val : 0.0, cq1 = 0.0, cq2 = 1.0))
    TX.AV_ZONE[] = mode == "zone" ? val : 0.0
    if geom == 2; TX.REGULAR_CENTRE[] = true; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64}); end
    PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
    tend = XS_TARGET / (V0 / 3)
    TF.setscalar!(:dtime, 1e-3 * (rI[2] - rI[1]) / V0); TF.setscalar!(:tmax, tend); TF.setscalar!(:thet, thet)
    TF.setscalar!(:stol, parse(Float64, get(ENV, "STOL", "0.05")))
    nst, why = try TD.run_julia!(maxsteps = 400000) catch e; (0, :failed) end
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); t = TF.getscalar(:timen)
    r = TF.var(:rn)[a:b+1] .- rmin; ρ = TF.var(:dn)[a:b]
    rc = 0.5 .* (r[1:end-1] .+ r[2:end])
    ρpre = geom == 0 ? fill(RHO0, length(rc)) : @. RHO0 * (1 + V0 * t / rc)^2
    jump = geom == 0 ? 4RHO0 : 64RHO0
    # outermost upward crossing of level L of f = (ρ − ρ_pre)/(ρ_post − ρ_pre) scanning inward: (fractional zone index, radius)
    function cross(L, ρpost)
        f = (ρ .- ρpre) ./ (ρpost .- ρpre)
        i = findlast(j -> f[j] >= L, eachindex(f)); (i === nothing || i == length(f)) && return (NaN, NaN)
        w = (f[i] - L) / (f[i] - f[i+1]); (i + w, rc[i] + w * (rc[i+1] - rc[i]))
    end
    # post-shock level: the finite-resolution plateau (median ρ in 0.6–0.9 of the shock radius, found first with the
    # analytic 64ρ0), because the numerical plateau sits ~10 % below 64ρ0 (wall heating) and would make the 90 % level ill-defined
    r50a = cross(0.5, jump)[2]
    sel = findall(x -> 0.6r50a <= x <= 0.9r50a, rc); ρpost = isempty(sel) ? jump : median(ρ[sel])
    # shock position and width from the velocity profile g = u/u_pre (1 pre-shock, 0 behind the shock; robust against the
    # density plateau being unresolved on coarse grids): outermost downward crossing of level L scanning inward
    xs = V0 / 3 * t
    uu = TF.var(:un)[a:b+1]; upre = uu[something(findfirst(>=(1.6xs), r), length(r))]   # pre-shock speed just ahead of the shock
    g = uu ./ upre
    function ucross(L)
        i = findlast(j -> g[j] <= L, eachindex(g)); (i === nothing || i == length(g)) && return (NaN, NaN)
        w = (L - g[i]) / (g[i+1] - g[i]); (i + w, r[i] + w * (r[i+1] - r[i]))
    end
    (iu50, ru50) = ucross(0.5); (iu10, _) = ucross(0.1); (iu90, _) = ucross(0.9)
    (i50, r50) = cross(0.5, ρpost); (i10, _) = cross(0.1, ρpost); (i90, _) = cross(0.9, ρpost)
    j = clamp(floor(Int, i50), 1, length(rc) - 1)
    (; n, steps = nst, why, poserr = (r50 - xs) / xs, width = i10 - i90, wu = iu90 - iu10, plateau = ρpost / jump, tfrac = t / tend,
       dr_over_r = (rc[j+1] - rc[j]) / rc[j])
end

if !isempty(ARGS) && ARGS[1] == "one"
    geom = parse(Int, ARGS[2])
    x = run_case(geom, ARGS[3], parse(Float64, ARGS[4]), ARGS[5], parse(Float64, ARGS[6]))
    if get(ENV, "PROFILE", "0") == "1"             # zone-by-zone ρ and u at the end
        a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); rmin = geom == 0 ? R0 : 0.0
        r = TF.var(:rn)[a:b+1] .- rmin; ρ = TF.var(:dn)[a:b]; u = TF.var(:un)[a:b+1]
        for k in eachindex(ρ); @printf("%3d r=%.4f rho/jump=%.3f u/v0=%.3f\n", k, r[k] / R0, ρ[k] / (geom == 0 ? 4RHO0 : 64RHO0), u[k] / V0); end
    end
    @printf("RESULT %d %s s=%s %s=%s n=%d steps=%d poserr=%.4f wu=%.3f plateau=%.3f why=%s t/tend=%.4f dr/r_shock=%.4f\n",
            geom, ARGS[3], ARGS[4], ARGS[5], ARGS[6], x.n, x.steps, x.poserr, x.wu, x.plateau, x.why, x.tfrac, x.dr_over_r)
else
    # old viscosity: ql1 chosen so that l = ql1 r̄ equals 2 × (current zone width in the shocked gas) on the uniform grid
    # (spherical: Δr/r ≈ 0.023 there → 0.05; planar, with the wall at r = R0: shocked zones are dx/4 wide → 2·dx/4/(1.1 R0))
    ql1 = Dict(2 => parse(Float64, get(ENV, "QL1S", "0.05")), 0 => parse(Float64, get(ENV, "QL1P", "0.0018")))
    avc = parse(Float64, get(ENV, "AVC", "2.0"))
    cases = Tuple[]
    for geom in (0, 2), (kind, s) in (("uniform", 1.0), ("geom", 1.02), ("geom", 1.05), ("geom", 1.1), ("geom", 1.2)),
        (m, v) in (("ql1", ql1[geom]), ("zone", avc))
        push!(cases, (geom, kind, s, m, v))
    end
    me = @__FILE__
    outs = Vector{String}(undef, length(cases))
    Threads.@threads for i in eachindex(cases)
        g, k, s, m, v = cases[i]
        cmd = addenv(`julia -t 1 --project=$(joinpath(@__DIR__, "env")) $me one $g $k $s $m $v`, "TITAN_NOLIB" => "1", "TITAN_MGR" => "700")
        outs[i] = try read(pipeline(cmd; stderr = devnull), String) catch e; "FAILED $e" end
    end
    for o in outs; foreach(l -> startswith(l, "RESULT") && println(l), split(o, "\n")); end
end
