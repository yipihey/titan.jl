# Tests of ModelAtom (R2, docs/POPIII_RADIATION_PLAN.md):  julia --project=<env> modelatom_test.jl
# 1. LTE: J = B_ν(T), any β, any collisions -> Saha–Boltzmann populations to round-off.
# 2. Case B (Lyman lines β = 0, others 1, no radiation, n_+ given): Hβ, Hα, He II 4686/1640 emissivities against
#    Storey & Hummer 1995 (atomdata/caseB_reference.csv). The l-collapsed atom is expected to deviate at n_e = 1e2–1e4
#    (no l-mixing) and to converge to SH95 at n_e ≳ 1e8–1e10.
using Printf
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "ModelAtom.jl"))
const MA = ModelAtom
ref = Dict{Tuple{String,String,Float64,Float64},Float64}()
for ln in eachline(joinpath(@__DIR__, "atomdata", "caseB_reference.csv"))
    (startswith(ln, "#") || startswith(ln, "ion,")) && continue
    p = split(ln, ","); p[3] == "emis" && p[7] == "SH95" && (ref[(p[1], p[2], parse(Float64, p[4]), parse(Float64, p[5]))] = parse(Float64, p[6]))
end
function main()
    H = MA.hydrogenic(1; nmax = 20); He2 = MA.hydrogenic(2; nmax = 20)
    println("1. LTE test (J = B, β = 0.3, collisions on): max |n/n_LTE − 1|")
    for (a, T, ne) in ((H, 6000.0, 1e12), (H, 1.2e4, 1e10), (H, 3e4, 1e14), (He2, 3e4, 1e12), (He2, 8e4, 1e10))
        s = MA.solve(a, T, ne, 1.0; J = ν -> MA.Bν(ν, T), β = (l, u) -> 0.3); L = MA.lte(a, T, ne, 1.0)
        dev = maximum(abs.(vcat(s.n, s.nplus) ./ vcat(L.n, L.nplus) .- 1))
        @printf("   %-5s T %.1e ne %.0e: %.2e %s\n", a.name, T, ne, dev, dev < 1e-7 ? "PASS" : "FAIL")
    end
    println("2. Case B emissivities 4πj/(n_e n_+) against SH95")
    caseB(l, u) = l == 1 ? 0.0 : 1.0
    println("   (n = 1 and 2 held as empty sinks, as in SH95's case B; nmax = 20 levels, l-collapsed)")
    for (a, ion, lines) in ((H, "HI", (("Hbeta", 2, 4), ("Halpha", 2, 3))), (He2, "HeII", (("4686", 3, 4), ("1640", 2, 3))))
        for (T, ne) in ((1e4, 1e2), (1e4, 1e4), (2e4, 1e2), (1e4, 1e6), (1e4, 1e8), (1e4, 1e10), (1e4, 1e11), (1e4, 1e12), (1e4, 1e13), (1e4, 1e14))
            n = MA.solve_upper(a, T, ne, 1.0; nfix = [0.0, 0.0], β = caseB)
            for (nm, l, u) in lines
                e = MA.line_emissivity(a, n, l, u) / ne
                r = get(ref, (ion, nm, T, ne), NaN)
                @printf("   %-4s %-6s T %.0e ne %.0e: %.4e  SH95 %.4e  ratio %.4f\n", ion, nm, T, ne, e, r, e / r)
            end
        end
    end
end
main()
function test_collisions()
    println("3. Collision rates: CHIANTI (H I Anderson et al. 2002; He II Ballance) n-collapsed vs Johnson 1972 and Vriens & Smeets (memory)")
    for (Z, stem) in ((1, joinpath(@__DIR__, "atomdata", "chianti", "h_1")), (2, joinpath(@__DIR__, "atomdata", "chianti", "he_2")))
        cu, nmx = MA.chianti_collapsed(stem); a = MA.hydrogenic(Z; nmax = 6)
        for (p, q) in ((1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5)), T in (1e4 * Z^2, 3e4 * Z^2)
            Υ = cu(p, q, T); qj = MA.johnson_exc(p, q, T; Z); qv = MA.vs_exc(p, q, T; Z)
            if Υ === nothing
                @printf("   Z %d %d-%d T %.1e: not in CHIANTI;  Johnson %.3e  V&S(memory) %.3e\n", Z, p, q, T, qj, qv)
            else
                qc = MA.q_from_ups(Υ, a.levels[p].g, a.levels[q].E - a.levels[p].E, T)
                @printf("   Z %d %d-%d T %.1e: CHIANTI Υ %.3f q %.3e  Johnson/CHIANTI %.2f  V&S(memory)/CHIANTI %.2f\n", Z, p, q, T, Υ, qc, qj / qc, qv / qc)
            end
        end
    end
end
function test_caseB2()
    println("4. Case B with CHIANTI (n = 1 → n′) + Johnson 1972 collisions, nmax = 40 (n = 1, 2 sinks)")
    caseB(l, u) = l == 1 ? 0.0 : 1.0
    for (Z, ion, stem, lines) in ((1, "HI", "h_1", (("Hbeta", 2, 4), ("Halpha", 2, 3))), (2, "HeII", "he_2", (("4686", 3, 4), ("1640", 2, 3))))
        a = MA.hydrogenic(Z; nmax = 40, coll = :chianti, chianti = joinpath(@__DIR__, "atomdata", "chianti", stem))
        for ne in (1e2, 1e8, 1e10, 1e11, 1e12, 1e13, 1e14)
            n = MA.solve_upper(a, 1e4, ne, 1.0; nfix = [0.0, 0.0], β = caseB)
            for (nm, l, u) in lines
                e = MA.line_emissivity(a, n, l, u) / ne; r = get(ref, (ion, nm, 1e4, ne), NaN)
                @printf("   %-4s %-6s T 1e4 ne %.0e: %.4e  SH95 %.4e  ratio %.4f\n", ion, nm, ne, e, r, e / r)
            end
        end
    end
end
test_collisions(); test_caseB2()
function test_heI()
    println("5. He I case B (1¹S held empty; singlet resonance lines n¹P → 1¹S thick) against Porter et al. 2012/13, nmax = 40")
    a = MA.heI(; nmax = 40)
    lab = [l.label for l in a.levels]; t(s) = findfirst(==(s), lab)
    lines = (("4471", t("1s.2p 3P"), t("1s.4d 3D")), ("5876", t("1s.2p 3P"), t("1s.3d 3D")), ("6678", t("1s.2p 1P"), t("1s.3d 1D")), ("10830", t("1s.2s 3S"), t("1s.2p 3P")))
    @printf("   %d levels; LTE check: ", MA.nlev(a))
    s = MA.solve(a, 2e4, 1e12, 1.0; J = ν -> MA.Bν(ν, 2e4), β = (l, u) -> 0.5); L = MA.lte(a, 2e4, 1e12, 1.0)
    @printf("max |n/n_LTE − 1| = %.2e\n", maximum(abs.(vcat(s.n, s.nplus) ./ vcat(L.n, L.nplus) .- 1)))
    caseB(l, u) = (l == 1 && a.levels[u].g == 3 && occursin("1P", a.levels[u].label)) || (l == 1 && occursin("singlet", a.levels[u].label)) ? 0.0 : 1.0
    pal = Dict{Tuple{String,Float64,Float64},Float64}()
    for ln in eachline(joinpath(@__DIR__, "atomdata", "caseB_reference.csv"))
        p = split(ln, ","); length(p) >= 7 && p[1] == "HeI" && p[7] == "Pal13" && (pal[(p[2], parse(Float64, p[4]), parse(Float64, p[5]))] = parse(Float64, p[6]))
    end
    for (T, ne) in ((1e4, 1e2), (1e4, 1e4), (1e4, 1e6), (1e4, 1e8), (1e4, 1e10), (1e4, 1e12), (1e4, 1e14), (2e4, 1e10), (2e4, 1e12), (2e4, 1e14))
        n = MA.solve_upper(a, T, ne, 1.0; nfix = [0.0], β = caseB)
        @printf("   T %.0e ne %.0e:", T, ne)
        for (nm, l, u) in lines
            e = MA.line_emissivity(a, n, l, u) / ne; r = get(pal, (nm, T, ne), NaN)
            @printf("  %s %.3f", nm, e / r)
        end
        println()
    end
end
test_heI()
function test_energy()
    println("6. Radiative energy balance: J = B_ν(T) and LTE populations -> lines + bf ≈ 0 (relative to the gross bf emission)")
    for (a, T, ne) in ((MA.hydrogenic(1; nmax = 20), 8000.0, 1e12), (MA.hydrogenic(2; nmax = 20), 4e4, 1e13), (MA.heI(; nmax = 15), 2e4, 1e12))
        L = MA.lte(a, T, ne, 1.0); J = ν -> MA.Bν(ν, T)
        rc = MA.radiative_cooling(a, L.n, L.nplus, ne, T; J, β = (l, u) -> 0.4)
        gross = sum(ne * L.nplus * MA.bf_energy(a, i, T; J)[2] for i in 1:MA.nlev(a))
        @printf("   %-5s T %.0e ne %.0e: lines %.3e  bf %.3e  (gross bf %.3e)  2γ %.3e  -> |lines+bf|/gross %.2e\n", a.name, T, ne, rc.lines, rc.bf, gross, rc.twophoton, abs(rc.lines + rc.bf) / gross)
    end
end
test_energy()
