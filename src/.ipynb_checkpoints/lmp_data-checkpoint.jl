module lmp_data

using Statistics
using Printf

include("Constants.jl")
#=
export Const_k_b, Const_n_a, Const_density_wat,
    Const_kcal2j, Const_kcalm2j, Const_kcalm2t,
    Const_g2kg, Const_kg2g, Const_gm2g, Const_gm2kg,
    Const_an2m, Const_an2nm, Const_nm2m, Const_cm2an, Const_dm2m, Const_m2dm,
    Const_fs2s, Const_ps2s, Const_ns2s
=#

include("useful_funcs.jl")
export dist, str2float, diag, copy_array, leastsq

include("read_dump.jl")
export read_dump, read_dump_prim, split_info!, split_atom, momentum!, mass!, time_step!, coord_scl!

include("box.jl")
export genr_box_diag, genr_box_inv, PBC!

include("diffusion.jl")
export MSD!, diffusion

include("RDF.jl")
export RDF


end # module





