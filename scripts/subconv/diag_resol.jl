# layer width and zones per H_P at the given steps: julia diag_resol.jl <snaps.jls> step...   (layer = diag_layer.jl definition)
using Serialization, Printf
shockif(s) = findfirst(<(-2e5), s.u)
function layer(s)
    cd = s.cd; js = shockif(s); js === nothing && return (js, 0, 0)
    bad(j) = cd.mach_i[j] < 0.5 && (cd.N2_L[j] < 0 || cd.N2_eq[j] < 0)
    jo = 0; for j in js-1:-1:max(2, js-4); bad(j) && (jo = j; break); end
    jo == 0 && return (js, 0, 0)
    ji = jo; while ji - 1 >= 2 && bad(ji - 1); ji -= 1; end
    (js, ji, jo)
end
function main()
    v = Any[]; open(ARGS[1]) do io; while !eof(io); push!(v, deserialize(io)); end; end
    for st in parse.(Int, ARGS[2:end])
        s = v[findlast(x -> x.step == st, v)]; cd = s.cd; js, ji, jo = layer(s)
        ji == 0 && (@printf("step %d: no layer\n", st); continue)
        w = sum(0.5 * (s.r[j+1] - s.r[j-1]) / cd.HP[j] for j in ji:jo); z = [cd.HP[j] / (0.5 * (s.r[j+1] - s.r[j-1])) for j in ji:jo]
        @printf("step %d  n_zones %d  R*=%.4g Rsun  layer mass %.3g Msun  width = sum dr/HP = %.2f HP  zones per HP at the layer interfaces = %s  total zones %d\n", st, length(s.ρ), s.r[js]/6.957e10, cd.layer[3], w, join((@sprintf("%.2f", x) for x in z), ","), length(s.ρ))
    end
end
main()
