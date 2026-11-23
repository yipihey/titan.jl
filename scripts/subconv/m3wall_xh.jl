using Serialization, Printf
open(ARGS[1]) do io
  while !eof(io)
    s = deserialize(io)
    @printf("%4d dt %.2e  Hp36 %.3e Hp37 %.3e Hp38 %.3e | H2 36 %.2e 37 %.2e 38 %.2e | T36 %.0f T37 %.0f w37 %.2e\n", s.step, s.dt, s.Y[8,36], s.Y[8,37], s.Y[8,38], s.Y[2,36], s.Y[2,37], s.Y[2,38], s.T[36], s.T[37], s.cd.tdc_w[37])
  end
end
