module distributedcalculate


export sample_and_calculate

using ..calculators, ..sample, ..clusters
using Distributed



function sample_and_calculate(inputs)
    start_dir = pwd()
    if inputs.cal.calculator.tmp_dir != start_dir
        cd(inputs.cal.calculator.tmp_dir)
    end
    return sample_multiple_adaptive_lines(inputs, basename="base-$(myid())")
end

function sample_nlines(cal, cl1, cl2, n)

end

end  # module calculate
