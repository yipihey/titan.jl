using WasmTarget, ForwardDiff
include(joinpath(@__DIR__, "..", "qsw.jl")); const Q = QSW; const RT = Roots
A = (Q.table_init, (Int32, Int32, Float64, Float64, Float64, Float64)); B = (Q.table_set, (Int32, Int32, Float64, Float64))
C = (Q.kap_R, (Float64, Float64)); D = (Q.shoot, (Float64, Float64, Float64, Float64, Float64))
E = (RT.all_roots, (Float64, Float64, Float64, Float64)); F = (RT.root, (Int32,))
for (name, fs) in (("AB", [A, B]), ("ABC", [A, B, C]), ("BC", [B, C]), ("AD", [A, D]), ("E", [E]), ("CA", [C, A]), ("DEFAB", [D, E, F, A, B]))
    try; b = compile_multi(fs); println(name, " OK ", length(b)); write("out/m_$name.wasm", b)
    catch e; println(name, " FAIL ", first(sprint(showerror, e), 160)); end
end
