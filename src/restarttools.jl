module restarttools
# methods that use fileaccess and other parts
export write_save_file,
       load_restart_file,
       continue_calculation,
       calculate_with_different_method,
       load_clusters_and_make_input,
       load_clusters_and_sample_input,
       calculate_adaptive_sample_inputs

using ..fileaccess, ..calculators, ..sample, ..distributedcalculate
using Distributed


"""
    write_save_file(fname, calculator, points, energy, cluster1, cluster2)

Saves restart/final information for energy calculation.

# Arguments
- `fname` : name of restartfile
- `calculator::Calculator` : calculator used in calculations
- `points` : 2d array of point where to calculate energy
- `energy` : energy for points that have been caculeted
- `cluster1` : cluster1
- `cluster2` : cluster2
"""
function write_save_file(fname, calculator::Calculator, points, energy, cluster1, cluster2)
    data = Dict("Energy" => energy, "Points" => points, "cluster1" => cluster1, "cluster2" => cluster2,
                "Basis" => calculator.basis, "Method" => calculator.method )
    save_jld_data(fname, data)
end


"""
    load_restart_file(fname)

Loads restart information from file `fname` and returns Dict with
"not_calculated", "calculeted" and "settings".

"settings" option has information of used computational method,
in form of Dict that as keys "Method", "Basis", "cluster1" and "cluster2"

"calculated" is a Dict with content "Energy" and "Points", that are 2-d arrays each,
and "clusters" that is tuple of 2 Clusters

"not_calculated" is a Dict with contenct "c1_points" that "c2_points" that are 2-d arrays
"""
function load_restart_file(fname)
    data = load_jld_data(fname)

    settings = Dict("Basis" => data["Basis"], "Method" => data["Method"],
                    "cluster1" => data["cluster1"], "cluster2" => data["cluster2"])

    sp = size(data["Points"])
    se = size(data["Energy"])
    needs_to_calculate = data["Points"][1:sp[1],sp[2]+1:se[2]]
    lc1 = length(data["cluster1"])
    c1_points = map(x -> x[1:lc1], needs_to_calculate)
    c2_points = map(x -> x[lc1+1:end], needs_to_calculate)
    not_calculated = Dict("c1_points" => c1_points, "c2_points" => c2_points)

    cp = data["Points"][1:se[1],1:se[2]]
    calculated = Dict("Points" => cp, "Energy" => data["Energy"],
                       "clusters" => (data["cluster1"], data["cluster2"] ) )
    return Dict("settings" => settings, "calculated" => calculated,
                "not_calculated" => not_calculated)
end



"""
    continue_calculation(fname, calculator::Calculator; batch_size=16)

This funtion is ways to use restart files to continue calculation

# Arguments
- `fname` : name of restart file
- `calculatror::Calculator` : calculator used for calculations (basis and method is read from file)
- `batch_size` : batch_size used in calculation
"""
function continue_calculation(fname, calculator::Calculator; batch_size=16)
    data = load_restart_file(fname)
    calculator.basis = data["setting"]["Basis"]
    calculator.method = data["setting"]["Method"]
    calculator.calculation_type = Energy()
    new_data = calculate_points(calculator, data["not_calculated"]["c1_points"],
                                data["not_calculated"]["c1_points"], batch_size=batch_size)
    energy = hcat(data["calculated"]["Energy"], new_data["Energy"])
    points = hcat(data["calculated"]["Points"], new_data["Points"])
    return Dict("Energy"=>energy, "Points"=>points, "Method"=>data["setting"]["Basis"],
                "Basis"=>data["setting"]["Basis"], "cluster1"=>data["calculated"]["cluster1"],
                "cluster2"=>data["calculated"]["cluster2"])
end


"""
    continue_calculation(fname, calculator::Calculator; batch_size=16, save_file="")

With this function you can use already chosen points on to which to do energy calculation

# Arguments
- `fname` : name of save/restart file where points are read
- `calculatror::Calculator` : calculator used for calculations
- `batch_size` : batch_size used in calculation
"""
function calculate_with_different_method(fname, calculator::Calculator; batch_size=16,
                                         save_file="")
    data = load_jld_data(fname)
    @info "File $(fname) loaded"
    lc1 = length(data["cluster1"])
    c1_points = map(x -> x[1:lc1], data["Points"])
    c2_points = map(x -> x[lc1+1:end], data["Points"])
    @info "Starting calculation"
    energy = calculate_points(calculator, c1_points, c2_points, batch_size=batch_size)
    out = Dict("Energy" => energy, "Points"=> data["Points"],
               "cluster1"=>data["cluster1"], "cluster2"=>data["cluster2"],
               "Method"=>calculator.method, "Basis"=>calculator.basis)
    if save_file != ""
        write_save_file(save_file, calculator, out["Points"], out["Energy"],
                        out["cluster1"], out["cluster2"])
    end
    return out
end


"""
load_clusters_and_sample_input(fname_cluster1, cluster2, calculator, nsamples;
                                      nlines=1, max_e=0, unit="cm-1", npoints=10,
                                      maxdis=9.0, sstep=0.1, startdistance=2.5)

Loads cluster1 from xyz-file and returns them all as an Array
that can then be used with calculate_adaptive_sample_inputs to calculate energy data.
Differs from load_clusters_and_sample_input by taking every point from file

# Arguments
- `fname_cluster1` : file where cluster1 is sampled
- `cluster2` :  cluster2
- `nlines` : number of lines to be sampled in calculation
- `max_e` : Point that is closest and has less energy than this will be starting point for a line
- `unit` : Unit in which `max_e` is given
- `npoints` : Number of points in potential
- `maxdis` : Maximum distance in potential calculation
- `sstep` : Search step size for adaptive search
- `startdistance` : Distance from which adaptive seach is started
"""
function load_clusters_and_make_input(fname_cluster1, cluster2, calculator;
                                      nlines=1, max_e=0, unit="cm-1", npoints=10,
                                      maxdis=9.0, sstep=0.1, startdistance=2.5)
    cluster1 = read_xyz(fname_cluster1)

    return  [InputAdaptiveSampler(calculator, c, cluster2,
                nlines, max_e, unit=unit, npoints=npoints, maxdis=maxdis,
                sstep=sstep, startdistance=startdistance) for c in cluster1 ]
end


"""
load_clusters_and_sample_input(fname_cluster1, cluster2, calculator, nsamples;
                                      nlines=1, max_e=0, unit="cm-1", npoints=10,
                                      maxdis=9.0, sstep=0.1, startdistance=2.5)

Loads cluster1 from xyz-file and takes `nsamples` samples of it and returns them as an Array
that can then be used with calculate_adaptive_sample_inputs to calculate energy data

# Arguments
- `fname_cluster1` : file where cluster1 is sampled
- `cluster2` :  cluster2
- `nsamples` :  number of samples taken from `fname_cluster1`
- `nlines` : number of lines to be sampled in calculation
- `max_e` : Point that is closest and has less energy than this will be starting point for a line
- `unit` : Unit in which `max_e` is given
- `npoints` : Number of points in potential
- `maxdis` : Maximum distance in potential calculation
- `sstep` : Search step size for adaptive search
- `startdistance` : Distance from which adaptive seach is started
"""
function load_clusters_and_sample_input(fname_cluster1, cluster2, calculator, nsamples;
                                      nlines=1, max_e=0, unit="cm-1", npoints=10,
                                      maxdis=9.0, sstep=0.1, startdistance=2.5)
    cluster1 = read_xyz(fname_cluster1)

    return  [InputAdaptiveSampler(calculator, rand(cluster1), cluster2,
                nlines, max_e, unit=unit, npoints=npoints, maxdis=maxdis,
                sstep=sstep, startdistance=startdistance) for _ in 1:nsamples ]
end


"""
    calculate_adaptive_sample_inputs(inputs; save_file_name="", save_step=nworkers())

Uses adaptive line sampler to `inputs` in distributed fashion and save data too

# Arguments
- `inputs` : calculation inputs array (eg. from load_clusters_and_sample_input)
- `save_file_name` : file where data saved
- `save_step=nworkers()` : number of calculated items before data is saved
"""
function calculate_adaptive_sample_inputs(inputs; save_file_name="", save_step=nworkers())
    t= collect(1:save_step:length(inputs))
    i_range = [ t[i-1]:t[i]-1  for i in 2:length(t) ]
    if t[end] < length(inputs)
        push!(i_range, t[end]:length(inputs))
    end
    @info "i_range $(i_range)"

    c = Channel(length(i_range))

    for r in i_range
        @async put!(c, pmap(distributedcalculate._sample_and_calculate,inputs[r]) )
    end

    tmp = take!(c)
    energy = hcat(map( x -> x["Energy"], tmp)...)
    points = hcat(map( x -> x["Points"], tmp)...)
    mindis = vcat(map( x -> x["Mindis"], tmp)...)
    if save_file_name != ""
        @info "Saving data to file $(save_file_name)"
        write_save_file(save_file_name, inputs[1].cal, points, energy,
                        inputs[1].cl1, inputs[2].cl1)
    end

    for i in i_range[1+1:end]
        tmp = take!(c)
        energy = hcat(energy, hcat(map( x -> x["Energy"], tmp)...) )
        points = hcat(points, hcat(map( x -> x["Points"], tmp)...) )
        mindis = vcat(mindis, vcat(map( x -> x["Mindis"], tmp)...) )
        if save_file_name != ""
            @info "Saving data to file $(save_file_name)"
            write_save_file(save_file_name, inputs[1].cal, points, energy,
                            inputs[1].cl1, inputs[2].cl1)
        end
    end
    return Dict("Energy" => energy, "Points" => points, "Mindis" => mindis)
end


end  # module restarttools
