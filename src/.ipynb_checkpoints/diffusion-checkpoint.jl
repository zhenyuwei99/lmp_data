"""
    MSD!(data::Dict; alpha=10)
This will add a variable "MSD" to `data` which is created by `read_dump`. "MSD" is the abbreviation of Mean Square Displacement (MSD), which could be used to calculate diffusion coefficients with the Einstein relationship. `alpha` is the ration between simulation time and the maximum of MSD calculating times.
"""
function MSD!(data::Dict; alpha=10)
    try
        data["coord_corr"]
    catch
        PBC!(data)
    end
    # Reading Input
    num_steps = data["num_steps"]
    num_atoms = data["num_atoms"]
    coord_corr = copy_array(data["coord_corr"])
    
    # MSD
    num_msd_steps = convert(Int64, round(num_steps / alpha))
    msd = zeros(num_msd_steps, num_atoms, 4)
    for step = 0 : num_msd_steps-1
        msd[step+1, :, 1:3] = mean((coord_corr[1+step:end, :, :] .- coord_corr[1:end-step, :, :]).^2, dims=1)
        msd[step+1, :, 4] = sum(msd[step+1, :, 1:3], dims=2)
    end
    data["MSD"] = msd
    return " MSD has been added to `data`"
end
    
"""
    diffusion(data::Dict, alpha=10)
This will return the diffusion coefficients of each atoms in `data`. 

Notes: time_step shoulb be added by `time_step!()` before calling this function.
"""
function diffusion(data::Dict, alpha=10)
    try 
        data["MSD"]
    catch
        MSD!(data, alpha)
    end
    
    try
        data["time_step"]
    catch
        error("\"time_step\" is not contained in `data`. Please call `time_step!()` before calling this function")
    end
    
    # Reading Input
    MSD = copy_array(data["MSD"])
    time_step = prod(data["time_step"])
    num_atoms = data["num_atoms"]
    
    # Diffusion
    num_step_msd = size(MSD, 1)
    t = [i*time_step for i=0:num_step_msd-1]
    res = zeros(num_atoms)
    
    for atom = 1:num_atoms
        res[atom] = leastsq(t, MSD[:, atom, 4])[2]
    end
    return res
end