"""
    function RDF(data_01::Dict, data_02::Dict, r_cut, num_bins)

This will retrun an array, which includes information of RDF between atom_01 and atom_02.

# `Kwg`
- `data_01`: A Dict varibale information of atom_01.
- `data_02`: A Dict varibale information of atom_02.
- `num_bins`: # of bins that separate `r_cut` into different parts.
- `r_cut`: cut off radius.

# Notice
- `data_01` and `data_02` should be Dict with information of one kind of atom, accoplising by `split_atom()`.
"""
function RDF(data_01::Dict, data_02::Dict, r_cut, num_bins)
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
    
    box_vol = genr_box_vol(data)
    box_diag = genr_box_diag(data)
    r_delta = r_cut / num_bins
    rho_02 = num_atoms_02 / box_vol
    
    # 

    g_x = [i for i in range(r_delta, r_cut, step=r_delta)]
    g_raw = zeros(num_steps, num_atoms_01, num_bins)
    g_time_avg = zeros(num_atoms_01, num_bins)
    d_vol = [4π * (n*r_delta)^2 * r_delta for n = 1:num_bins]
    for atom_01 = 1:num_atoms_01
        for atom_02 = 1:num_atoms_02
            r_scl = coord_scl_02[:, atom_02, :] - coord_scl_01[:, atom_01, :]
            r_scl .-= round.(r_scl)
            r_diff = r_scl * box_diag
            r_diff = sqrt.(sum(r_diff.^2, dims=2))
            step = [i[1] for i in findall(x->0<x<=r_cut, r_diff)]
            bin = convert.(Int64, ceil.(r_diff ./ r_delta))
            for i in step
                g_raw[i, atom_01, bin[i]] += 1
            end
        end
        g_time_avg[atom_01, :] = mean(g_raw[:, atom_01, :], dims=1) ./ d_vol'
    end
    
    res = hcat(g_x, mean(g_time_avg, dims=1)' / rho_02)
    return res
end