"""
Contains methods to implement distributed calculations
"""
module distributedcalculate


export sample_ntimes, calculate_points

using ..calculators, ..sample, ..clusters
using Distributed


"""
    _sample_and_calculate(inputs)

Internal helper function to help implement distributed calculations
"""
function _sample_and_calculate(inputs::InputAdaptiveSampler)
    start_dir = pwd()
    if inputs.cal.calculator.tmp_dir != start_dir
        cd(inputs.cal.calculator.tmp_dir)
    end
    return sample_multiple_adaptive_lines(inputs, basename="base-$(myid())", id="Pid $(myid())")
end


"""
    sample_ntimes(input, n=1)

Calls sample_multiple_adaptive_lines `n` times using distributed pmap
"""
function sample_ntimes(input::InputAdaptiveSampler, n=1)
    inp = [input for i in 1:n]
    tmp = pmap(_sample_and_calculate, inp)
    energy = hcat(map( x -> x["Energy"], tmp)...)
    points = hcat(map( x -> x["Points"], tmp)...)
    mindis = vcat(map( x -> x["Mindis"], tmp)...)
    return Dict("Energy" => energy, "Points" => points, "Mindis" => mindis)
end



"""
    _calculate_points(cal, c1_points, c2_points)

Internal helper function to ease use of distributed calculations.
Users should call "calculate_points" instead.
"""
function _calculate_points(cal, c1_points, c2_points)
    start_dir = pwd()
    if cal.calculator.tmp_dir != start_dir
        cd(cal.calculator.tmp_dir)
    end
    return bsse_corrected_energy(cal, c1_points, c2_points, basename="base-$(myid())",
                                 id="Pid $(myid())")
end


"""
    calculate_points(calculator, c1_points, c2_points; batch_size=16)

Calculates using distributed pmap energy between two given clusters arrays

# Atributes
- `calculator` : calculator used in energy calculations
- `c1_points`  : array of clusters
- `c2_points`  : array of clusters
- `batch_size` : bactch size for distributed calculations
"""
function calculate_points(calculator, c1_points, c2_points; batch_size=16)
    l = length(c1_points)
    if length(c2_points) != l
        @error "cluster1 array has different dimensions than cluster2 array"
        throw(error("cluster1 array has different dimensions than cluster2 array"))
    end
    if batch_size > l
        c1_in = [c1_points]
        c2_in = [c2_points]
    else
        c1_in=[ c1_points[(i-batch_size+1):i]  for i in batch_size:batch_size:l ]
        c2_in=[ c2_points[(i-batch_size+1):i]  for i in batch_size:batch_size:l ]
        if mod(l,batch_size) != 0
            push!(c1_in, c1_points[end-mod(l,batch_size):end])
            push!(c2_in, c2_points[end-mod(l,batch_size):end])
        end
    end
    ca_in = fill(ca ,size(c1_in))
    energy = pmap(_calculate_points, ca_in, c1_in, c2_in)
    energy_out = vcat(energy...)
    return reshape(energy_out, size(c1_points))
end


end  # module calculate
