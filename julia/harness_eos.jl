# Per-routine harness (phase 1 demonstration).
#   1. Run SUBC to completion through libtitan.so to get a realistic state.
#   2. ccall the single Fortran subroutine `eos` (ideal-gas path, leos = 2) on
#      that snapshot and on a perturbed snapshot; compare bitwise with a
#      hand-written Julia translation (`eos_ideal!`).
#   3. libm demonstration: same procedure for `opac` with lopac = 6 (Kramers +
#      Thomson), whose only libm call is pow(T, -3.5).  `parity = true` routes
#      pow/sqrt through glibc via ccall; `parity = false` uses Julia's own.
#   4. Micro-checks of the other libm/libgcc entry points the Fortran uses:
#      exp, log, pow, and __powidf2 (gfortran's x**n for INTEGER-variable n).
include(joinpath(@__DIR__, "TitanFortran.jl"))
using .TitanFortran
const TF = TitanFortran

# ---------------------------------------------------------------- math switch
const LIBM = "libm.so.6"
const LIBGCC = "libgcc_s.so.1"
c_sqrt(x::Float64) = ccall((:sqrt, LIBM), Float64, (Float64,), x)
c_pow(x::Float64, y::Float64) = ccall((:pow, LIBM), Float64, (Float64, Float64), x, y)
c_exp(x::Float64) = ccall((:exp, LIBM), Float64, (Float64,), x)
c_log(x::Float64) = ccall((:log, LIBM), Float64, (Float64,), x)
c_powi(x::Float64, n::Integer) = ccall((:__powidf2, LIBGCC), Float64, (Float64, Cint), x, n)
msqrt(x, parity) = parity ? c_sqrt(x) : sqrt(x)
mpow(x, y, parity) = parity ? c_pow(x, y) : x^y

# ------------------------------------------------- Julia translation of eos.f
"eos.f, leos = 2 branch, k = ngrs-1 .. ngre+2. Operation order copied from the Fortran."
function eos_ideal!(o, dn, tn, k1, k2, gam, gmu, ck, cm0; parity = false)
    for k in k1:k2
        o[:pgn][k]    = ck / (cm0 * gmu) * dn[k] * tn[k]
        o[:egn][k]    = ck / (cm0 * gmu) * tn[k] / (gam - 1.0)
        o[:xnen][k]   = 0.5 / (cm0 * gmu) * dn[k]
        o[:delsn][k]  = (gam - 1.0) / gam
        o[:cpn][k]    = ck / (cm0 * gmu) / o[:delsn][k]
        o[:betan][k]  = 1.0 / (tn[k] * dn[k])
        o[:gamman][k] = gam
        o[:asn][k]    = msqrt(gam * o[:pgn][k] / dn[k], parity)
        o[:dlpgdltn][k] = 1.0; o[:dlpgdldn][k] = 1.0
        o[:dlegdltn][k] = 1.0; o[:dlegdldn][k] = 0.0
        o[:dlpedltn][k] = 1.0; o[:dlpedldn][k] = 1.0
        o[:dldsdltn][k] = 0.0; o[:dldsdldn][k] = 0.0
        o[:dlcpdltn][k] = 0.0; o[:dlcpdldn][k] = 0.0
        o[:dlbedltn][k] = -1.0; o[:dlbedldn][k] = -1.0
    end
    o
end
const EOS_OUT = (:pgn, :egn, :xnen, :delsn, :cpn, :betan, :gamman, :asn,
                 :dlpgdltn, :dlpgdldn, :dlegdltn, :dlegdldn, :dlpedltn, :dlpedldn,
                 :dldsdltn, :dldsdldn, :dlcpdltn, :dlcpdldn, :dlbedltn, :dlbedldn)

"opac.f, lopac = 6 branch."
function opac_kramers!(o, dn, tn, k1, k2, gmu, cm0, csige; parity = false)
    for k in k1:k2
        xkff = 8.6e23 * dn[k] * mpow(tn[k], -3.5, parity)
        o[:chifn][k] = xkff + csige / cm0
        o[:xken][k]  = xkff
        o[:xkpn][k]  = xkff
        o[:xnen][k]  = 0.5 / (gmu * cm0) * dn[k]
        o[:dlcfdltn][k] = -3.5 * xkff / o[:chifn][k]
        o[:dlcfdldn][k] = 1.0 * xkff / o[:chifn][k]
        o[:dlkedltn][k] = -3.5; o[:dlkedldn][k] = 1.0
        o[:dlkpdltn][k] = -3.5; o[:dlkpdldn][k] = 1.0
        o[:dlpedltn][k] = 1.0;  o[:dlpedldn][k] = 1.0
    end
    o
end
const OPAC_OUT = (:chifn, :xken, :xkpn, :xnen, :dlcfdltn, :dlcfdldn, :dlkedltn, :dlkedldn,
                  :dlkpdltn, :dlkpdldn, :dlpedltn, :dlpedldn)

# --------------------------------------------------------------- comparison
"Zero the outputs, call the Fortran routine, then the Julia one; bitwise compare on k1:k2."
function compare_routine(fsym::Symbol, outs, jl!; parity)
    k1 = Int(TF.getscalar(:ngrs)) - 1; k2 = Int(TF.getscalar(:ngre)) + 2
    for n in outs; TF.var(n)[k1:k2] .= 0.0; end
    ccall(Libdl_dlsym(fsym), Cvoid, ())
    F = Dict(n => copy(TF.var(n)) for n in outs)
    J = Dict(n => zeros(TF.MGR) for n in outs)
    jl!(J, k1, k2; parity)
    bad = Dict{Symbol,Int}()
    for n in outs
        m = count(reinterpret(UInt64, F[n][k1:k2]) .!= reinterpret(UInt64, J[n][k1:k2]))
        m > 0 && (bad[n] = m)
    end
    (; npoints = k2 - k1 + 1, bad)
end
Libdl_dlsym(s) = TF.Libdl.dlsym(TF.load(), s)

function main()
    TF.init(8; rundir = joinpath(TF.ROOT, "lib", "jl_runs", "SUBC_harness"))
    TF.run!()
    println("SUBC run to jstep = $(TF.getscalar(:jstep)), t = $(TF.getscalar(:timen)); leos = ",
            TF.getscalar(:leos), ", lopac = ", TF.getscalar(:lopac))
    dn = TF.var(:dn); tn = TF.var(:tn)
    T0 = copy(tn)
    gam, gmu = TF.getscalar(:gam), TF.getscalar(:gmu)
    ck, cm0, csige = TF.getscalar(:ck), TF.getscalar(:cm0), TF.getscalar(:csige)
    k1 = Int(TF.getscalar(:ngrs)) - 1; k2 = Int(TF.getscalar(:ngre)) + 2
    perturb!() = for k in k1:k2; tn[k] = T0[k] * (1.0 + 0.37 * sin(1.7k)); end

    for (label, prep!) in (("snapshot", () -> (tn .= T0)), ("perturbed T", perturb!))
        for parity in (false, true)
            prep!()
            r = compare_routine(:eos_, EOS_OUT,
                    (o, a, b; parity) -> eos_ideal!(o, dn, tn, a, b, gam, gmu, ck, cm0; parity); parity)
            println("eos  (leos=2) $label, parity=$parity: $(r.npoints) zones x $(length(EOS_OUT)) outputs: ",
                    isempty(r.bad) ? "BITWISE IDENTICAL" : "mismatches $(r.bad)")
        end
    end
    lop = TF.getscalar(:lopac); TF.setscalar!(:lopac, Int32(6))
    for (label, prep!) in (("snapshot", () -> (tn .= T0)), ("perturbed T", perturb!))
        for parity in (false, true)
            prep!()
            r = compare_routine(:opac_, OPAC_OUT,
                    (o, a, b; parity) -> opac_kramers!(o, dn, tn, a, b, gmu, cm0, csige; parity); parity)
            println("opac (lopac=6) $label, parity=$parity: $(r.npoints) zones x $(length(OPAC_OUT)) outputs: ",
                    isempty(r.bad) ? "BITWISE IDENTICAL" : "mismatches $(r.bad)")
        end
    end
    TF.setscalar!(:lopac, lop); tn .= T0
    TF.finalize()

    # ------------------------------------------- micro-checks of libm/libgcc
    N = 1_000_000
    xs = [exp(-30 + 60 * ((i * 0.6180339887498949) % 1)) for i in 1:N]   # 1e-13 .. 1e13
    ys = [-50 + 100 * ((i * 0.7548776662466927) % 1) for i in 1:N]
    cnt(f, g, v) = count(x -> reinterpret(UInt64, f(x)) != reinterpret(UInt64, g(x)), v)
    println("\nmicro-checks over $N arguments (count of results differing in any bit):")
    println("  exp: Julia vs glibc            ", cnt(exp, c_exp, ys))
    println("  log: Julia vs glibc            ", cnt(log, c_log, xs))
    println("  pow(x,-3.5): Julia vs glibc    ", cnt(x -> x^-3.5, x -> c_pow(x, -3.5), xs))
    println("  pow(x,0.25): Julia vs glibc    ", cnt(x -> x^0.25, x -> c_pow(x, 0.25), xs))
    println("  sqrt: Julia vs glibc           ", cnt(sqrt, c_sqrt, xs))
    n3 = Ref(3)                                   # runtime integer exponent (like r**mup1)
    println("  x^n (n=3 runtime Int): Julia vs __powidf2   ", cnt(x -> x^n3[], x -> c_powi(x, 3), xs))
    println("  x*x*x vs __powidf2(x,3)                     ", cnt(x -> x * x * x, x -> c_powi(x, 3), xs))
    n2 = Ref(2)
    println("  x^n (n=2 runtime Int): Julia vs __powidf2   ", cnt(x -> x^n2[], x -> c_powi(x, 2), xs))
    println("  x^4 (literal): Julia vs (x*x)*(x*x) [gfortran x**4]  ", cnt(x -> x^4, x -> (x * x) * (x * x), xs))
    println("  x^3 (literal): Julia vs x*(x*x)    [gfortran x**3]  ", cnt(x -> x^3, x -> x * (x * x), xs))
end
main()
