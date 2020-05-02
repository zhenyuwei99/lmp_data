
"""
    genr_box_diag(data::Dict)
This will return the matrix of box size information.
"""
function genr_box_diag(data::Dict)
    box_diag = copy_array(data["box_info"])
    box_diag = mean(box_diag, dims=1)[1, :, :]
    box_diag = box_diag[:, 2] - box_diag[:, 1]
    box_diag = diag(box_diag)
    return box_diag
end

"""
    genr_box_vol(data::Dict)
This will return the volume of box in unit of your chossen.
"""
function genr_box_vol(data::Dict)
    box_vol = copy_array(data["box_info"])
    box_vol = mean(box_vol, dims=1)[1, :, :]
    box_vol = box_vol[:, 2] - box_vol[:, 1]
    return prod(box_vol)
end

"""
    genr_box_inv(data::Dict)
This will return the inverse matrix of box size information.
"""
function genr_box_inv(data::Dict)
    box_diag = genr_box_diag(data)
    return inv(box_diag)
end

"""
      PBC!(data::Dict)
This will add a variable "coord_corr" to `data` which is created by `read_dump`. `coord_corr` has handled the PBC issues.
"""
function PBC!(data::Dict)
    coord = data["coord"]
    box_diag = genr_box_diag(data)
    box_inv = genr_box_inv(data)
    coord_corr = zeros(size(coord))
    for step = 1 : data["num_steps"]
        coord_corr[step, :, :] = coord[step, :, :] * box_inv
    end
    coord_diff = coord_corr[2:end, :, :] .- coord_corr[1:end-1, :, :]
    coord_diff .-= round.(coord_diff)
    coord_corr[2:end, :, :] .= 0
    coord_corr[1, :, :] *= box_diag
    for step = 2 : data["num_steps"]
        coord_corr[step, :, :] = coord_corr[step-1, :, :] .+ coord_diff[step-1, :, :] * box_diag
    end
    data["coord_corr"] = coord_corr;
    
    println("\"coord_corr\"\thas been added to `data`")
end