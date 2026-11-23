using Serialization, Printf
d = ARGS[1]; f = filter(x -> occursin(r"^snaps_\d+_\d+\.jls$", x), readdir(d))[1]
open(joinpath(d, f)) do io
  while !eof(io)
    s = deserialize(io); c = s.cd
    c === nothing && continue
    fc = c.fconv; ft = c.Ftot; j = argmax(filter(!isnan, abs.(c.Fconv_over_Ftot)) |> x -> isempty(x) ? [0.0] : x)
    fo = [isnan(v) ? 0.0 : v for v in c.Fconv_over_Ftot]; jm = argmax(fo)
    sup = [isnan(v) ? -Inf : v for v in c.sup]
    @printf("step %d dt %.3e  max Fconv/Ftot %.3g at int %d; max sup %.3g; layer %s; ncv(sup>0,wτ>0) %d\n", s.step, s.dt, fo[jm], jm, maximum(sup), join((@sprintf("%.3g", v) for v in c.layer), ","), count(>(0), [isnan(v) ? 0.0 : v for v in c.fconv]))
  end
end
