function str_cmp(str_cmp::AbstractString, str_goal::String, len)
    str_cmp[1:len] == str_goal[1:len]
end

function str_cmp(str_cmp::AbstractString, str_goal::String)
    len = length(str_goal)
    if length(str_cmp) <= len
        return false
    end
    str_cmp[1:len] == str_goal[1:len]
end

function dist(coord_01, coord_02)
    num_dims = length(coord_01)
    dist = 0
    for dim = 1 : num_dims
        dist += (coord_01[dim] - coord_02[dim])^2
    end
    sqrt(dist)
end
function str2float(str::AbstractString)
    str = split(str)
    [parse(Float64, str[n]) for n = 1:length(str)]
end


function diag(vec::AbstractArray)
    len = length(vec)
    result = zeros(len, len)
    for i = 1 : len
        result[i, i] = vec[i]
    end
    result
end

function copy_array(goal)
    res = zeros(typeof(goal).parameters[1], size(goal))
    res[:] = goal[:]
    return res
end