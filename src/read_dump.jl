################
# read_dump
################
"""
    read_dump_prim(dump_name::String)

This will read data in `dump_name` file directly, returning a Dict variable, without spliting them into different array.
"""
function read_dump_prim(dump_name::String)
    io = open(dump_name, "r")
    str = read(io, String)
    
    # Reading all digital data
    reg = r"\b[\d\se\+ -\.]+\n"
    res = findall(reg, str)
    res = [str[n] for n in res]
    num_steps = Int64(length(res) / 4)
    
    # Reading Simulation Info
    step_vec = parse(Int64, res[1][2:end])
    num_atoms = parse(Int64, res[2][2:end])
    box_info = [str2float(x) for x in split(res[3], "\n")[2:end-1]]
    box_info = Array(hcat(box_info...)')
    box_temp = copy_array(box_info)
    atom_info = [str2float(x) for x in split(res[4], "\n")[2:end-1]]
    atom_info = Array(hcat(atom_info...)')
    atom_temp = copy_array(atom_info)
    
    box_info = Array{Float64}(undef, (num_steps, size(box_info, 1), size(box_info, 2)))
    box_info[1, :, :] .= box_temp
    atom_info = Array{Float64}(undef, (num_steps, size(atom_info, 1), size(atom_info, 2)))
    atom_info[1, :, :] .= sortslices(atom_temp, dims=1)
    
    # Rearanging Data
    for step = 2:num_steps
        step_vec = vcat(step_vec, parse(Int64, res[(step-1)*4+1][2:end]))
        box_info[step, :, :] = Array(hcat([str2float(x) for x in split(res[(step-1)*4+3], "\n")[2:end-1]]...)')
        atom_temp = Array(hcat([str2float(x) for x in split(res[(step-1)*4+4], "\n")[2:end-1]]...)')
        atom_info[step, :, :] = sortslices(atom_temp, dims=1)
    end
    
    # Output
    res = Dict(
        "step_vec" => step_vec,
        "num_atoms" => num_atoms,
        "num_steps" => num_steps,
        "box_info" => box_info,
        "atom_info" => atom_info,
    )
end

"""
    split_dump!(data::Dict; id=false, atom_type=false, mol=false, element=false, mass=false, coord=false, coord_scl=false, vel=false, force=false, acc=false)

This will split data readed from dump file, adding variable named same as `kwg` before to `data`.
The value of `kwg` can be `false` or column range of property of interest.

# Example

```julia-repl
dump_name = "xxxx" # name of dump file
data = read_dump_prim(dump_name)
split_dump!(data, id=1, atom_type=2, coord=3:5)
```
Code above will add three elements, named as "id", "atom_type", and "coord", to data which contains variable in column 1, 2, 3:5 of `data["atom_info"]`
"""
function split_dump!(data::Dict; id=false, atom_type=false, mol=false, element=false, mass=false, 
        coord=false, coord_scl=false, vel=false, force=false, acc=false)
    atom_info = data["atom_info"]
    if id != false
        data["id"] = convert.(Int64, atom_info[:, :, id])
    end
    if atom_type != false
        data["atom_type"] = convert.(Int64, atom_info[:, :, atom_type])
    end
    if mol != false
        data["mol"] = convert.(Int64, atom_info[:, :, mol])
    end
    if element != false
        data["element"] = convert.(Int64, atom_info[:, :, element])
    end
    if mass != false
        data["mass"] = atom_info[:, :, mass]
    end
    if coord != false
        data["coord"] = atom_info[:, :, coord]
    end
    if coord_scl != false
        data["coord_scl"] = atom_info[:, :, coord_scl]
    end
    if vel != false
        data["vel"] = atom_info[:, :, vel]
    end
    if acc != false
        data["acc"] = atom_info[:, :, acc]
    end
    if force != false
        data["force"] = atom_info[:, :, force]
    end
    return 0
end
        
"""
    read_dump(dump_name::String; id=false, atom_type=false, mol=false, element=false, mass=false, coord=false, coord_scl=false, vel=false, force=false, acc=false)

This will read data in `dump_name` file and return a Dict which contains elements
- `step_vec`
- `num_atoms`
- `num_steps`
- `box_info`
- `atom_info`

# `kwg`

Each `kwg` can be `false` or column range of property of interest in `atom_info`. If a non-`false` `kwg` is activated, a element with the same name as `kwg` will be add to Dict that will return.

# Example

```julia-repl
dump_name = "xxxx"  # name of dump file
data = read_dump(dump_name, id=1, atom_type=2, coord=3:5)
```
Code above will return a Dict contains elements:
- `step_vec`
- `num_atoms`
- `num_steps`
- `box_info`
- `atom_info`
- `id`
- `atom_type`
- `coord`

Where `id`, `atom_type`, and `coord` are values of `atom_info` in column 1, 2, 3:5 respectively.
"""
function read_dump(dump_name::String; id=false, atom_type=false, mol=false, element=false, mass=false, 
        coord=false, coord_scl=false, vel=false, force=false, acc=false)
    data = read_dump_prim(dump_name)
    split_dump!(data, id=id, atom_type=atom_type, mol=mol, element=element, mass=mass, 
        coord=coord, vel=vel, acc=acc, force=force)
    return data
end

"""
    function momentum!(data::Dict)

This will add a component called "momentum" into `data` while `vel` and `mass` is contained in `data`
"""
function momentum!(data::Dict)
    try
        data["mass"]
    catch
        error("Info of \"mass\" is not contained in `data` ")
    end
    try
        data["vel"]
    catch
        error("Info of \"vel\" is not contained in `data` ")
    end
    
    momentum = zeros(size(data["vel"]))
    mass_diag = diag(data["mass"])
    for dim = 1 : 3
        momentum[:, :, dim] .= data["vel"][:, :, dim] * mass_diag
    end
    data["momentum"] = momentum
    return "momentum has been added to `data`"
end

"""
    function mass!(data::Dict, mass_vec, id_vec=false)

This will add a component called "mass" into `data`. Elements in `mass_vec` are the mass of atoms with id correspond to `id_vec`. If `id_vec` is set to default, it means `id_vec = [1, 2, 3, ... length(mass_vec)]'
"""
function mass!(data::Dict, mass_vec, id_vec=false)
    if id_vec == false
        id_vec = [i for i = 1:length(mass_vec)]
    end
    id = copy_array(data["atom_type"][1, :])
    mass = zeros(size(id))
    for i = 1:length(mass_vec)
        mass[findall(x->x==id_vec[i], id)] .= mass_vec[i]
    end
    data["mass"] = mass
    return "mass has been added to `data`"
end


"""
    function time_step!(data::Dict, time_len, converter, dump_len)
This will add a component, which is an array contain all information needed to calulate the time in MD_simulation, called "time_step" into `data`.

# Parameters:
- `data`: Dict variable create by `read_dump()`.
- `time_len`: length of time step of simulation.
- `converter`: unit converter that convert time_len to second.
- `dump_len`: length of dump step.
"""
function time_step!(data::Dict, time_len, converter, dump_len)
    data["time_step"] = [time_len, converter, dump_len]
    return "time_step has been added to `data`"
end