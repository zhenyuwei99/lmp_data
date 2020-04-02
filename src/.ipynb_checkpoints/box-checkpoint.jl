```
    genr_box_diag(data::Dict)
This will return the matrix of box size information
```
function genr_box_diag(data::Dict)
    box_diag = copy_array(data["box_info"])
    box_diag = mean(box_diag, dims=1)[1, :, :]
    box_diag = box_diag[:, 2] - box_diag[:, 1]
    box_diag = diag(box_diag)
    return box_diag
end

```
    genr_box_inv(data::Dict)
This will return the inverse matrix of box size information
```
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
    box_inv = inv(box_diag)
    coord_corr = zeros(size(coord))
    for step = 1 : data["num_steps"]
        coord_corr[step, :, :] = coord[step, :, :] * box_inv
    end
    coord_corr .-= round.(coord_corr)
    for step = 1 : data["num_steps"]
        coord_corr[step, :, :] *= box_diag
    end
    data["coord_corr"] = coord_corr;
end