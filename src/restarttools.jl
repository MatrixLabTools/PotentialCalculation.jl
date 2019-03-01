"""
module restarttools

Primary tools to do calculations.
Also contains all methods to restart calculations.
"""
module restarttools

export calculate_adaptive_sample_inputs,
       calculate_with_different_method,
       continue_calculation,
       load_clusters_and_make_input,
       load_clusters_and_sample_input,
       load_data_file,
       load_restart_file,
       write_restart_file,
       write_save_file

using ..fileaccess, ..calculators, ..sample, ..clusters
using Distributed, ProgressMeter


"""
    write_save_file(fname, calculator, points, energy, cluster1, cluster2)

Saves final information for energy calculation.

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
    write_restart_file(fname, calculator, points, restart_energy, cluster1, cluster2)

Saves restart information for energy calculation.

# Arguments
- `fname` : name of restartfile
- `calculator::Calculator` : calculator used in calculations
- `points` : 2d array of point where to calculate energy
- `energy` : energy for points that have been caculeted
- `cluster1` : cluster1
- `cluster2` : cluster2
"""
function write_restart_file(fname, calculator::Calculator, points, restart_energy, cluster1, cluster2)
    @info "Writing restart file $(fname)"
    data = Dict("restart_energy" => restart_energy, "Points" => points, "cluster1" => cluster1, "cluster2" => cluster2,
                "Basis" => calculator.basis, "Method" => calculator.method )
    save_jld_data(fname, data)
end

"""
    load_data_file(fname)

Loads saved data
"""
function load_data_file(fname)
    return load_jld_data(fname)
end

"""
    load_restart_file(fname)

Loads restart information from file `fname` and adds to it key `not_calculated`,
which holds information of which collumns of `Points` have not been calculated.
"""
function load_restart_file(fname)
    data = load_jld_data(fname)
    not_calculated = Set(1:size(data["Points"])[2])
    if haskey(data, "restart_energy")
        setdiff!(not_calculated, x[1] for x in data["restart_energy"])
    else
        @warn "File is not a restart file" fname
    end
    push!(data, "not_calculated"=>not_calculated)
    return data
end



"""
    continue_calculation(fname, calculator::Calculator; save_file="", restart_file="", pbar=true)

Restarts calculation from given file

# Arguments
- `fname` : file from which calculation is restarted
- `calculator::Calculator` : calculator used in calculations
- `save_file` : file where final results are saved, if given
- `restart_file` : file where new restart information is saved, if given
- `pbar=true` : show progress bar
"""
function continue_calculation(fname, calculator::Calculator; save_file="", restart_file="", pbar=true)
    data = load_restart_file(fname)
    @info "File $(fname) loaded - continua calculation"
    flush(stdout)
    calculator.basis = data["Basis"]
    calculator.method = data["Method"]


    lc1 = length(data["cluster1"])
    lc2 = length(data["cluster2"])
    rc1 =  1:lc1           #UnitRange for cluster1
    rc2 =  lc1+1:lc1+lc2   #UnitRange for cluster2

    ncol = length(data["not_calculated"])
    t = max(trunc(Int,ceil(2*nworkers()/size(data["Points"])[1])), 2)
    l = min(t, ncol)
    @debug "Using $(l) Channels"
    c = Channel(l)

    if pbar
        prog = Progress(size(data["Points"])[1]*ncol*3 , dt=1, desc="Calculating points:")
        pchannel = RemoteChannel(()->Channel{Bool}(2*nworkers()), myid())
        @async while take!(pchannel)
            ProgressMeter.next!(prog)
            flush(stderr)
            flush(stdout)
        end
    else
        pchannel = undef
    end

    @async for collumn in data["not_calculated"]
        c1_points = map( x -> x[1:lc1], data["Points"][:,collumn])
        c2_points = map( x -> x[lc1+1:end], data["Points"][:,collumn])
        @async put!(c , (collumn ,pmap((x,y)->_calculate_points(calculator, x,y ,pchannel=pchannel),
                                      c1_points, c2_points ) ))
        sleep(0.1)  # make sure that FIFO queue is filled in right order
    end



    for collumn in 2:ncol
        push!(data["restart_energy"], take!(c))
        if restart_file != ""
            write_restart_file(restart_file, calculator, data["Points"], data["restart_energy"],
                               data["cluster1"], data["cluster2"])
            flush(stdout)
        end
    end

    push!(data["restart_energy"], take!(c))
    energy = similar(data["Points"], Float64)
    for col in data["restart_energy"]
        energy[:, col[1]] = col[2]
    end
    if save_file != ""
        write_save_file(save_file, calculator, data["Points"], energy,
                           data["cluster1"], data["cluster2"])
        flush(stdout)
    end
    return Dict("Energy" => energy, "Points"=> data["Points"],
               "cluster1"=>data["cluster1"], "cluster2"=>data["cluster2"],
               "Method"=>calculator.method, "Basis"=>calculator.basis)

end



"""
    calculate_with_different_method(fname, calculator::Calculator;
                                    save_file="", restart_file="", pbar=true)

With this function you can use already chosen points on to which to do energy calculation.

# Arguments
- `fname` : name of save/restart file where points are read
- `calculatror::Calculator` : calculator used for calculations
- `save_file=""` : save final results here, if not empty
- `restart_file=""` : save restarts here, if not empty
- `pbar=true` : show progress bar
"""
function calculate_with_different_method(fname, calculator::Calculator;
                                         save_file="", restart_file="", pbar=true)
    data = load_jld_data(fname)
    @info "File $(fname) loaded - calculating with differnt method"
    flush(stdout)
    lc1 = length(data["cluster1"])
    c1_points = map(x -> x[1:lc1], data["Points"])
    c2_points = map(x -> x[lc1+1:end], data["Points"])
    #inputs = fill(calculator, size(c1_points)[1] )

    ncol = length(c1_points[1,:])
    t = max(trunc(Int,ceil(2*nworkers()/size(c1_points)[1])), 2)
    l = min(t, ncol)
    @debug "Using $(l) Channels"
    c = Channel(l)

    if pbar
        prog = Progress(length(c1_points)*3 ,dt=1, desc="Calculating points:")
        pchannel = RemoteChannel(()->Channel{Bool}(2*nworkers()), myid())
        @async while take!(pchannel)
            ProgressMeter.next!(prog)
            flush(stderr)
            flush(stdout)
        end
    else
        pchannel = undef
    end

    @async for collumn in 1:ncol
        @async put!(c , (collumn ,pmap( (x,y) -> _calculate_points(calculator,x,y, pchannel=pchannel),
                                        c1_points[:,collumn], c2_points[:,collumn] ) ))
        sleep(0.1)  # make sure that FIFO queue is filled in right order
    end

    tmp_energy = []
    for collumn in 2:ncol
        push!(tmp_energy, take!(c))
        if restart_file != ""
            write_restart_file(restart_file, calculator, data["Points"], tmp_energy,
                               data["cluster1"], data["cluster2"])
            flush(stdout)
        end
    end

    push!(tmp_energy, take!(c))
    energy = similar(c1_points, Float64)
    for col in tmp_energy
        energy[:, col[1]] = col[2]
    end
    if save_file != ""
        write_save_file(save_file, calculator, data["Points"], energy,
                           data["cluster1"], data["cluster2"])
        flush(stdout)
    end
    return Dict("Energy" => energy, "Points"=> data["Points"],
               "cluster1"=>data["cluster1"], "cluster2"=>data["cluster2"],
               "Method"=>calculator.method, "Basis"=>calculator.basis)
end


"""
load_clusters_and_make_input(cluster1::String, cluster2::Cluster, calculator, nsamples;
                                      max_e=0, unit="cm-1", npoints=10,
                                      maxdis=9.0, sstep=0.1, startdistance=2.5)

Loads cluster1 from xyz-file and returns them all as an Array
that can then be used with [`calculate_adaptive_sample_inputs`](@ref) to calculate energy data.
Differs from [`load_clusters_and_sample_input`](@ref) by taking every point from file

# Arguments
- `cluster1::String` : file from where cluster1 is sampled
- `cluster2::Cluster` :  cluster2
- `max_e` : Point that is closest and has less energy than this will be starting point for a line
- `unit` : Unit in which `max_e` is given
- `npoints` : Number of points in potential
- `maxdis` : Maximum distance in potential calculation
- `sstep` : Search step size for adaptive search
- `startdistance` : Distance from which adaptive seach is started
"""
function load_clusters_and_make_input(cluster1::String, cluster2::Cluster, calculator;
                                      max_e=0, unit="cm-1", npoints=10,
                                      maxdis=9.0, sstep=0.1, startdistance=2.5)
    cluster1 = read_xyz(cluster1)

    return  [InputAdaptiveSampler(calculator, c, cluster2,
                1, max_e, unit=unit, npoints=npoints, maxdis=maxdis,
                sstep=sstep, startdistance=startdistance) for c in cluster1 ]
end


"""
load_clusters_and_make_input(cluster1::Cluster, cluster2::Cluster, calculator;
                                      nlines=1, max_e=0, unit="cm-1", npoints=10,
                                      maxdis=9.0, sstep=0.1, startdistance=2.5)

Loads cluster1 from xyz-file and returns them all as an Array
that can then be used with [`calculate_adaptive_sample_inputs`](@ref) to calculate energy data.
Differs from [`load_clusters_and_sample_input`](@ref) by taking every point from file

# Arguments
- `cluster1::Cluster` : cluster1
- `cluster2::Cluster` : cluster2
- `nlines` : number of lines to be sampled in calculation
- `max_e` : Point that is closest and has less energy than this will be starting point for a line
- `unit` : Unit in which `max_e` is given
- `npoints` : Number of points in potential
- `maxdis` : Maximum distance in potential calculation
- `sstep` : Search step size for adaptive search
- `startdistance` : Distance from which adaptive seach is started
"""
function load_clusters_and_make_input(cluster1::Cluster, cluster2::Cluster, calculator;
                                      nlines=1, max_e=0, unit="cm-1", npoints=10,
                                      maxdis=9.0, sstep=0.1, startdistance=2.5)

    return  [InputAdaptiveSampler(calculator, cluster1, cluster2,
                1, max_e, unit=unit, npoints=npoints, maxdis=maxdis,
                sstep=sstep, startdistance=startdistance) for _ in 1:nlines ]
end



"""
load_clusters_and_sample_input(fname_cluster1, cluster2, calculator, nsamples;
                                      max_e=0, unit="cm-1", npoints=10,
                                      maxdis=9.0, sstep=0.1, startdistance=2.5)

Loads cluster1 from xyz-file and takes `nsamples` samples of it and returns them as an Array
that can then be used with [`calculate_adaptive_sample_inputs`](@ref) to calculate energy data

# Arguments
- `cluster1` : file from where cluster1 is sampled
- `cluster2` : cluster2
- `nsamples` : number of samples taken from `fname_cluster1`
- `max_e` : Point that is closest and has less energy than this will be starting point for a line
- `unit` : Unit in which `max_e` is given
- `npoints` : Number of points in potential
- `maxdis` : Maximum distance in potential calculation
- `sstep` : Search step size for adaptive search
- `startdistance` : Distance from which adaptive seach is started
"""
function load_clusters_and_sample_input(cluster1::String, cluster2, calculator, nsamples;
                                      max_e=0, unit="cm-1", npoints=10,
                                      maxdis=9.0, sstep=0.1, startdistance=2.5)
    cluster1 = read_xyz(cluster1)

    return  [InputAdaptiveSampler(calculator, rand(cluster1), cluster2,
                1, max_e, unit=unit, npoints=npoints, maxdis=maxdis,
                sstep=sstep, startdistance=startdistance) for _ in 1:nsamples ]
end


"""
function load_clusters_and_sample_input(cluster1::String, cluster2::String, calculator, nsamples;
                                      nlines=1, max_e=0, unit="cm-1", npoints=10,
                                      maxdis=9.0, sstep=0.1, startdistance=2.5)

# Arguments
- `cluster1` : file from where cluster1 is sampled
- `cluster2` : file from where cluster2 is sampled
- `nsamples` : number of samples taken from `fname_cluster1`
- `nlines` : number of lines to be sampled in calculation
- `max_e` : Point that is closest and has less energy than this will be starting point for a line
- `unit` : Unit in which `max_e` is given
- `npoints` : Number of points in potential
- `maxdis` : Maximum distance in potential calculation
- `sstep` : Search step size for adaptive search
- `startdistance` : Distance from which adaptive seach is started
"""
function load_clusters_and_sample_input(cluster1::String, cluster2::String, calculator, nsamples;
                                      nlines=1, max_e=0, unit="cm-1", npoints=10,
                                      maxdis=9.0, sstep=0.1, startdistance=2.5)
    cluster1 = read_xyz(cluster1)
    cluster2 = read_xyz(cluster2)

    return  [InputAdaptiveSampler(calculator, rand(cluster1), rand(cluster2),
                nlines, max_e, unit=unit, npoints=npoints, maxdis=maxdis,
                sstep=sstep, startdistance=startdistance) for _ in 1:nsamples ]
end



"""
    calculate_adaptive_sample_inputs(inputs; save_file_name="", save_step=nworkers(), pbar=true)

Uses [`adaptive_line_sampler`](@ref) to `inputs` in distributed fashion

# Arguments
- `inputs` : calculation inputs array (eg. from [`load_clusters_and_sample_input`](@ref))
- `save_file_name` : file where data saved
- `save_step=nworkers()` : number of calculated items before data is saved
- `pbar=true` : show progress bar
"""
function calculate_adaptive_sample_inputs(inputs; save_file_name="", save_step=nworkers(), pbar=true)
    @info "Starting adaptive sampling"
    flush(stdout)
    t= collect(1:save_step:length(inputs))
    i_range = [ t[i-1]:t[i]-1  for i in 2:length(t) ]
    if t[end] < length(inputs)
        push!(i_range, t[end]:length(inputs))
    end

    lc = min(length(i_range), trunc(Int, ceil(2*nworkers()/save_step)) )
    c = Channel(lc)

    if pbar
        prog = ProgressUnknown(dt=1, desc="Adaptive Sampling Points calculated:")
        pchannel = RemoteChannel(()->Channel{Bool}(2*nworkers()), myid())
        @async while take!(pchannel)
            ProgressMeter.next!(prog)
            flush(stderr)
            flush(stdout)
        end
    else
        pchannel = undef
    end

    @async for r in i_range
        @async put!(c, pmap( x->_sample_and_calculate(x, pchannel=pchannel), inputs[r]) )
    end

    tmp = take!(c)
    energy = hcat(map( x -> x["Energy"], tmp)...)
    points = hcat(map( x -> x["Points"], tmp)...)
    mindis = vcat(map( x -> x["Mindis"], tmp)...)
    if save_file_name != ""
        @info "Saving data to file $(save_file_name)"
        write_save_file(save_file_name, inputs[1].cal, points, energy,
                        inputs[1].cl1, inputs[1].cl2)
        flush(stdout)
    end

    for i in i_range[2:end]
        tmp = take!(c)
        energy = hcat(energy, hcat(map( x -> x["Energy"], tmp)...) )
        points = hcat(points, hcat(map( x -> x["Points"], tmp)...) )
        mindis = vcat(mindis, vcat(map( x -> x["Mindis"], tmp)...) )
        if save_file_name != ""
            @info "Saving data to file $(save_file_name)"
            write_save_file(save_file_name, inputs[1].cal, points, energy,
                            inputs[1].cl1, inputs[1].cl2)
            flush(stdout)
        end
    end
    #pbar && ProgressMeter.finish!(prog)
    return Dict("Energy" => energy, "Points" => points, "Mindis" => mindis,
                "cluster1"=>inputs[1].cl1, "cluster2"=>inputs[1].cl2,
                "Method"=>inputs[1].cal.method, "Basis"=>inputs[1].cal.basis)
end



"""
    _sample_and_calculate(inputs; pchannel=undef )

Internal helper function to help implement distributed calculations
"""
function _sample_and_calculate(inputs::InputAdaptiveSampler; pchannel=undef)
    start_dir = pwd()
    if inputs.cal.calculator.tmp_dir != start_dir
        cd(inputs.cal.calculator.tmp_dir)
    end
    return sample_multiple_adaptive_lines(inputs, basename="base-$(myid())",
                                         id="Pid $(myid())", pchannel=pchannel)
end


"""
    _calculate_points(cal, c1_points, c2_points; pchannel=undef)

Internal helper function to ease use of distributed calculations.
Users should call "calculate_points" instead.
"""
function _calculate_points(cal, c1_points, c2_points; pchannel=undef)
    start_dir = pwd()
    if cal.calculator.tmp_dir != start_dir
        cd(cal.calculator.tmp_dir)
    end
    return bsse_corrected_energy(cal, c1_points, c2_points, basename="base-$(myid())",
                                 id="Pid $(myid())", pchannel=pchannel)
end

end  # module restarttools
