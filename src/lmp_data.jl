module lmp_data

include("useful_funcs.jl")
export dist, str2float, diag, copy_array

include("read_dump.jl")
export read_dump, read_dump_prim, split_dump! 

include("box.jl")
export genr_box_diag, genr_box_inv, PBC!

end # module





