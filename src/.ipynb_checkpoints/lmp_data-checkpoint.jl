module lmp_data

include("useful_funcs.jl")
export dist, str2float, diag, copy_array

include("read_dump.jl")
export read_dump, read_dump_prim, split_dump! 


end # module





