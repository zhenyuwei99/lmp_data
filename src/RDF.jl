"""
    mutable struct RDF
This contains the information of the radiual distribution function (RDF)

# Elements List:
- `r::Vector{Float64}`: the radius vector which coorespond the bins coordinates.
- `g::Vector{Float64}`: the amplitude vector of each bin corresponds to `r`.
- `rho::Float64`: the number density of neighbour atoms.
- `num_bins`: # of bins in calculation of RDF.
- `dr`: length of each bins in calculation of RDF, with the unit of length.
"""
mutable struct RDF
    r::Vector{Float64}
    g::Vector{Float64}
    rho::Float64
    num_bins::Int64
    dr::Float64
end

"""
    function RDF(data_01::Dict, data_02::Dict; r_cut=10, num_bins=200)
This will retrun a varibale in RDF type, which includes information of RDF between atom_01 and atom_02.

# `Kwg`
- `data_01`: A Dict varibale information of atom_01.
- `data_02`: A Dict varibale information of atom_02.
- `num_bins`: # of bins that separate `r_cut` into different parts.
- `r_cut`: cut off radius.

# Notice
- `data_01` and `data_02` should be Dict with information of one kind of atom, accoplising by `split_atom()`.
"""
function RDF(data_01::Dict, data_02::Dict; r_cut=10, num_bins=200)
    # Reading Input
    
    try
        data_01["coord_scl"]
    catch
        coord_scl!(data_01)
    end
    
    try
        data_02["coord_scl"]
    catch
        coord_scl!(data_02)
    end
    
    num_atoms_01 = data_01["num_atoms"]
    num_atoms_02 = data_02["num_atoms"]
    num_steps = data_01["num_steps"]
    coord_scl_01 = copy_array(data_01["coord_scl"])
    coord_scl_02 = copy_array(data_02["coord_scl"])
    
    box_vol = genr_box_vol(data_02)
    box_diag = genr_box_diag(data_02)
    dr = r_cut / num_bins
    rho_02 = num_atoms_02 / box_vol
    
    # 

    g_x = [i for i in range(dr, r_cut, step=dr)]
    g_raw = zeros(num_steps, num_atoms_01, num_bins)
    g_time_avg = zeros(num_atoms_01, num_bins)
    d_vol = [4π * (n*dr)^2 * dr for n = 1:num_bins]
    for atom_01 = 1:num_atoms_01
        for atom_02 = 1:num_atoms_02
            r_scl = coord_scl_02[:, atom_02, :] - coord_scl_01[:, atom_01, :]
            r_scl .-= round.(r_scl)
            r_diff = r_scl * box_diag
            r_diff = sqrt.(sum(r_diff.^2, dims=2))
            step = [i[1] for i in findall(x->0<x<=r_cut, r_diff)]
            bin = convert.(Int64, ceil.(r_diff ./ dr))
            for i in step
                g_raw[i, atom_01, bin[i]] += 1
            end
        end
        g_time_avg[atom_01, :] = mean(g_raw[:, atom_01, :], dims=1) ./ d_vol'
    end
    
    g = vec(mean(g_time_avg, dims=1)/rho_02)
    return RDF(g_x, g, rho_02, num_bins, dr)
    
end

"""
    function RDF_double(data_01::Dict, data_02::Dict; r_cut=10, num_bins=200, dim=1)
This will retrun a varibale in RDF type, which includes information of RDF between atom_01 and atom_02 symetry with axies `dim`.

# `Kwg`
- `data_01`: A Dict varibale information of atom_01.
- `data_02`: A Dict varibale information of atom_02.
- `num_bins`: # of bins that separate `r_cut` into different parts.
- `r_cut`: cut off radius.
- `dim`: symetry axies. 1, 2, 3 for x, y, z axies respectively.
# Notice
- `data_01` and `data_02` should be Dict with information of one kind of atom, accoplising by `split_atom()`.
"""
function RDF_double(data_01::Dict, data_02::Dict; r_cut=10, num_bins=200, dim=1)
    # Reading Input
    
    try
        data_01["coord_scl"]
    catch
        coord_scl!(data_01)
    end
    
    try
        data_02["coord_scl"]
    catch
        coord_scl!(data_02)
    end
    
    num_atoms_01 = data_01["num_atoms"]
    num_atoms_02 = data_02["num_atoms"]
    num_steps = data_01["num_steps"]
    coord_scl_01 = copy_array(data_01["coord_scl"])
    coord_scl_02 = copy_array(data_02["coord_scl"])
    
    box_vol = genr_box_vol(data_02)
    box_diag = genr_box_diag(data_02)
    dr = 2 * r_cut / num_bins
    rho_02 = num_atoms_02 / box_vol
    
    # 

    g_x = [i for i in range(-r_cut + dr, r_cut, step=dr)]
    g_raw = zeros(num_steps, num_atoms_01, num_bins)
    g_time_avg = zeros(num_atoms_01, num_bins)
    d_vol = [4π * ( (n - num_bins/2)*dr)^2 * dr for n = 1:num_bins]
    for atom_01 = 1:num_atoms_01
        for atom_02 = 1:num_atoms_02
            r_scl = coord_scl_02[:, atom_02, :] - coord_scl_01[:, atom_01, :]
            r_scl .-= round.(r_scl)
            r_diff = r_scl * box_diag
            r_diff = r_diff[:, dim] ./ abs.(r_diff[:, dim] .+ 1e-6) .* sqrt.(sum(r_diff.^2, dims=2))
            step = [i[1] for i in findall(x->-r_cut<x<=r_cut, r_diff)]
            bin = convert.(Int64, ceil.(r_diff ./ dr).+ num_bins / 2) 
            for i in step
                g_raw[i, atom_01, bin[i]] += 1
            end
        end
        g_time_avg[atom_01, :] = mean(g_raw[:, atom_01, :], dims=1) ./ d_vol'
    end
    
    g = vec(mean(g_time_avg, dims=1)/rho_02)
    return RDF(g_x, g, rho_02, num_bins, dr)
    
end

"""
    function num_coordination(g::RDF, r_start, r_end)
This will return a value in Float64, which express the coordination number between range of `r_statr` to `r_end`.

# Paramters List:
- `g`: the RDF information in type `RDF` which is created by function RDF.
- `r_start`: the start radius of calculation of coordination number.
- `r_end`: the end radius of calculation of coordination number.
"""
function num_coordination(g::RDF, r_start, r_end)
	pos = findall(x->r_start<x<r_end, g.r)
	g_range = g.g[pos]
	r_range = g.r[pos]
	res = sum(g.rho * 4π * r_range.^2 * g.dr .* g_range)
	return res
end