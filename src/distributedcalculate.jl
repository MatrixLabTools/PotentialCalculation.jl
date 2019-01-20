module distributedcalculate


export sample_and_calculate, sample_ntimes, calculate_points

using ..calculators, ..sample, ..clusters
using Distributed



function sample_and_calculate(inputs)
    start_dir = pwd()
    if inputs.cal.calculator.tmp_dir != start_dir
        cd(inputs.cal.calculator.tmp_dir)
    end
    return sample_multiple_adaptive_lines(inputs, basename="base-$(myid())", id="Pid $(myid())")
end

function sample_ntimes(input, n)
    inp = [input for i in 1:n]
    tmp = pmap(sample_and_calculate, inp)
    energy = hcat(map( x -> x["Energy"], tmp)...)
    points = hcat(map( x -> x["Points"], tmp)...)
    mindis = vcat(map( x -> x["Mindis"], tmp)...)
    return Dict("Energy" => energy, "Points" => points, "Mindis" => mindis)
end

function _calculate_points(cal, c1_points, c2_points)
    start_dir = pwd()
    if cal.calculator.tmp_dir != start_dir
        cd(cal.calculator.tmp_dir)
    end
    return bsse_corrected_energy(cal, c1_points, c2_points, basename="base-$(myid())",
                                 id="Pid $(myid())")
end

function calculate_points(ca, c1_points, c2_points; batch_size=16)
    l = length(c1_points)
    if length(c2_points) != l
        @error "cluster1 has different dimensions than cluster2"
        throw(error("cluster1 has different dimensions than cluster2"))
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
