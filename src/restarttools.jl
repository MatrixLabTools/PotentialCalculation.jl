"""
module restarttools

Primary tools to do calculations.
Also contains all methods to restart calculations.
"""
module restarttools

export calculate_adaptive_sample_inputs,
       calculate_energy_for_xyzfile,
       calculate_with_different_method,
       calculate_potential,
       continue_calculation,
       create_inputs,
       load_clusters_and_make_input,
       load_clusters_and_sample_input,
       load_data_file,
       load_restart_file,
       test_work,
       write_restart_file,
       write_save_file

using ..fileaccess, ..calculators, ..sample, ..clusters
using Distributed, ProgressMeter


"""
    write_save_file(fname, calculator, points, energy, cluster1, cluster2)

Saves final information for energy calculation.

# Arguments
- `fname` : name of restartfileusing PotentialCalculation
- `calculator::Calculator` : calculator used in calculations
- `points` : 2d array of point where to calculate energy
- `energy` : energy for points that have been caculeted
- `cluster1` : cluster1
- `cluster2` : cluster2
"""
function write_save_file(
    fname::AbstractString,
    calculator::Calculator,
    points,
    energy,
    cluster1,
    cluster2
)
    data = Dict("Energy" => energy, "Points" => points, "cluster1" => cluster1, "cluster2" => cluster2,
                "Basis" => calculator.basis, "Method" => calculator.method )
    flush(stdout)
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
function write_restart_file(
    fname::AbstractString,
    calculator::Calculator,
    points,
    restart_energy,
    cluster1,
    cluster2
)
    @info "Writing restart file $(fname)"
    data = Dict("restart_energy" => restart_energy, "Points" => points, "cluster1" => cluster1, "cluster2" => cluster2,
                "Basis" => calculator.basis, "Method" => calculator.method )
    flush(stdout)
    save_jld_data(fname, data)
end

"""
    load_data_file(fname)

Loads saved data
"""
function load_data_file(fname::AbstractString)
    return load_jld_data(fname)
end

"""
    load_restart_file(fname)

Loads restart information from file `fname` and adds to it key `not_calculated`,
which holds information of which collumns of `Points` have not been calculated.
"""
function load_restart_file(fname::AbstractString)
    data = load_jld_data(fname)
    if haskey(data, "restart_energy")
        if ! ( typeof(data["restart_energy"]) <: Tuple)
            #NOTE Old style restart file
            not_calculated = Set(1:size(data["Points"])[2])
            setdiff!(not_calculated, x[1] for x in data["restart_energy"])
            push!(data, "not_calculated"=>not_calculated)
            iscalculated = trues(size(data["Points"]))
            for i in not_calculated
                iscalculated[:,i] .= false
            end
            energy = similar(data["Points"],Float64)
            for (i,e) in data["restart_energy"]
                energy[:,i] .= e
            end
            data["restart_energy"] = (energy, iscalculated)
        end
    else
        @warn "File is not a restart file" fname
    end
    return data
end



"""
continue_calculation(fname, calculator::Calculator; save_file="", restart_file="",
                              save_after=nworkers(),  pbar=true)

Restarts calculation from given file

# Arguments
- `fname` : file from which calculation is restarted
- `calculator::Calculator` : calculator used in calculations
- `save_file` : file where final results are saved, if given
- `restart_file` : file where new restart information is saved, if given
- `save_after=nworkers()` : make restart file when given ammount of points is calculated
- `pbar=true` : show progress bar
"""
function continue_calculation(
    fname::AbstractString,
    calculator::Calculator;
    save_file="",
    restart_file="",
    save_after=nworkers(),
    pbar=true
)
    data = load_restart_file(fname)
    @info "File $(fname) loaded - continuing calculation"
    calculator.basis = data["Basis"]
    calculator.method = data["Method"]

    energy,iscal = data["restart_energy"]

    lc1 = length(data["cluster1"])
    lc2 = length(data["cluster2"])
    rc1 =  1:lc1           #UnitRange for cluster1
    rc2 =  lc1+1:lc1+lc2   #UnitRange for cluster2

    lcal = length(data["Points"])-length(data["Points"][iscal])

    if pbar
        prog = Progress(lcal  ,dt=1, desc="Calculating points:")
    end

    jobs = RemoteChannel(()->Channel(2*nworkers()))
    results = RemoteChannel(()->Channel(2*nworkers()))

    iscalculated = deepcopy(iscal)
    @async for i in eachindex(data["Points"])
        if ! iscal[i]
            x = (calculator,data["Points"][i][rc1], data["Points"][i][rc2])
            put!(jobs, (i,x) )
        end
    end

    for p in workers()
        remote_do(parallel_work, p, jobs, results, _calculate_points)
    end

    n = 0
    for j in 1:lcal
        i,tmp = take!(results)
        energy[i] = tmp
        iscalculated[i] = true
        pbar && ProgressMeter.next!(prog)
        n += 1
        if restart_file != "" && n >= save_after && j < lcal
            write_restart_file(restart_file, calculator, data["Points"], (energy,iscalculated),
                               data["cluster1"], data["cluster2"])
            n = 0
        end
    end

    if save_file != ""
        write_save_file(save_file, calculator, data["Points"], energy,
                           data["cluster1"], data["cluster2"])
    end
    return Dict("Energy" => energy, "Points"=> data["Points"],
               "cluster1"=>data["cluster1"], "cluster2"=>data["cluster2"],
               "Method"=>calculator.method, "Basis"=>calculator.basis)

end



"""
    calculate_potential(fname::AbstractString, calculator::Calculator;
                        save_file="", restart_file="", pbar=true, save_after=nworkers()
                       )
With this function you can use already chosen points on to which to do energy calculation.

# Arguments
- `fname` : name of save/restart file where points are read
- `calculatror::Calculator` : calculator used for calculations

# Keywords
- `save_file=""` : save final results here, if not empty
- `restart_file=""` : save restarts here, if not empty
- `pbar=true` : show progress bar
- `save_after=nworkers()` : make restart file when given ammount of points is calculated
"""
function calculate_potential(fname::AbstractString, calculator::Calculator;
                             save_file="", restart_file="", pbar=true, save_after=nworkers()
                            )
    data = load_jld_data(fname)
    @info "File $(fname) loaded - calculating with different method"
    flush(stdout)
    lc1 = length(data["cluster1"])
    c1_points = map( x-> x[1:lc1], data["Points"])
    c2_points = map( x-> x[lc1+1:end], data["Points"])

    if pbar
        prog = Progress(length(c1_points) ,dt=1, desc="Calculating points:")
    end

    jobs = RemoteChannel(()->Channel(2*nworkers()))
    results = RemoteChannel(()->Channel(2*nworkers()))

    @async for i in eachindex(c1_points)
        x = (calculator,c1_points[i], c2_points[i])
        put!(jobs, (i,x) )
    end

    for p in workers()
        remote_do(parallel_work, p, jobs, results, _calculate_points)
    end

    energy = similar(c1_points, Float64)
    iscalculated = falses(size(energy))
    n = 0
    lcal=length(c1_points)
    for j in 1:lcal
        i,tmp = take!(results)
        energy[i] = tmp
        iscalculated[i] = true
        pbar && ProgressMeter.next!(prog)
        n += 1
        if restart_file != "" && n >= save_after && j < lcal
            write_restart_file(restart_file, calculator, data["Points"], (energy,iscalculated),
                               data["cluster1"], data["cluster2"])
            n = 0
        end
    end

    if save_file != ""
        write_save_file(save_file, calculator, data["Points"], energy,
                           data["cluster1"], data["cluster2"])
    end
    return Dict("Energy" => energy, "Points"=> data["Points"],
               "cluster1"=>data["cluster1"], "cluster2"=>data["cluster2"],
               "Method"=>calculator.method, "Basis"=>calculator.basis)
end


function calculate_with_different_method(fname::AbstractString, calculator::Calculator;
                             save_file="", restart_file="", pbar=true, save_after=nworkers()
                            )
    @warn "calculate_with_different_method is deprecated use calculate_potential instead"
    return calculate_potential(fname, calculator; save_file=save_file, restart_file=restart_file,
                               pbar=pbar, save_after=save_after)
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
    @warn "load_clusters_and_make_input deprecated use create_inputs instead"
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

    @warn "load_clusters_and_make_input deprecated use create_inputs instead"
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
    @warn "load_clusters_and_sample_input deprecated use create_inputs instead"
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
    @warn "load_clusters_and_sample_input is deprecated use create_inputs instead"
    cluster1 = read_xyz(cluster1)
    cluster2 = read_xyz(cluster2)

    return  [InputAdaptiveSampler(calculator, rand(cluster1), rand(cluster2),
                nlines, max_e, unit=unit, npoints=npoints, maxdis=maxdis,
                sstep=sstep, startdistance=startdistance) for _ in 1:nsamples ]
end



"""
    calculate_adaptive_sample_inputs(inputs; save_after="", save_after=nworkers(), pbar=true)

Uses [`adaptive_line_sampler`](@ref) to `inputs` in distributed fashion

# Arguments
- `inputs` : calculation inputs array (eg. from [`create_inputs`](@ref))
- `save_after` : file where data saved
- `save_after=nworkers()` : number of calculated items before data is saved
- `pbar=true` : show progress bar
"""
function calculate_potential(inputs::AbstractArray{InputAdaptiveSampler};
                             save_file="", save_after=nworkers(), pbar=true)
    @info "Starting adaptive sampling"
    flush(stdout)
    t= collect(1:save_after:length(inputs))
    i_range = [ t[i-1]:t[i]-1  for i in 2:length(t) ]
    if t[end] < length(inputs)
        push!(i_range, t[end]:length(inputs))
    end

    if pbar
        prog = Progress(length(inputs) ,dt=1, desc="Calculating points:")
    end

    jobs = RemoteChannel(()->Channel(2*nworkers()))
    results = RemoteChannel(()->Channel(2*nworkers()))

    @async for x in inputs
        put!(jobs, (1, (x,) ) )
    end

    for p in workers()
        remote_do(parallel_work, p, jobs, results, _sample_and_calculate)
    end

    _,tmp = take!(results)
    energy = tmp["Energy"]
    points = tmp["Points"]
    mindis = tmp["Mindis"]
    pbar && ProgressMeter.next!(prog)
    if save_file != ""
        @info "Saving data to file $(save_file)"
        write_save_file(save_file, inputs[1].cal, points, energy,
                        inputs[1].cl1, inputs[1].cl2)
    end
    if length(inputs) > 1
        for i in 2:length(inputs)
            _,tmp = take!(results)
            energy = hcat(energy, tmp["Energy"] )
            points = hcat(points, tmp["Points"] )
            mindis = vcat(mindis, tmp["Mindis"] )
            pbar && ProgressMeter.next!(prog)
            if save_file != ""
                @info "Saving data to file $(save_file)"
                write_save_file(save_file, inputs[1].cal, points, energy,
                                inputs[1].cl1, inputs[1].cl2)
            end
        end
    end
    #pbar && ProgressMeter.finish!(prog)
    return Dict("Energy" => energy, "Points" => points, "Mindis" => mindis,
                "cluster1"=>inputs[1].cl1, "cluster2"=>inputs[1].cl2,
                "Method"=>inputs[1].cal.method, "Basis"=>inputs[1].cal.basis)
end


function calculate_adaptive_sample_inputs(inputs::AbstractArray{InputAdaptiveSampler};
                                         save_file_name="", save_step=nworkers(), pbar=true)
    @warn "calculate_adaptive_sample_inputs is deprecated use calculate_potential instead"
    return calculate_potential(inputs, save_file=save_file_name, save_after=save_step, pbar=pbar)
 end


"""
    create_inputs(cluster1, cluster2, calculator::Calculator;  kwargs...)



Creates inputs for the calculations. Use [`calculate_adaptive_sample_inputs`](@ref) on
results to do the calculation itself.

Inputs for the molecules (named cluster1/2) can be loaded from xyz-file or given by as
[`Cluster`](@ref) structures. If given as file names, then `nsamples` points are randomly
picked from the trajectories for calculation. Total number of lines sampled is thus
`nlines` times `nsamples`.

Parallelisation is done over `nsamples`, thus it is recommended for it to be multiple of
number workers processes.

# Arguments
- `cluster1` : something that can be interpreted as [`Cluster`](@ref) or an Array of it
- `cluster2` : something that can be interpreted as [`Cluster`](@ref) or an Array of it
- `calculator::Calculator` : [`Calculator`](@ref) used in sampling

# Keywords
- max_e=15000 : energy treshold in sampling. Enegy is lower than this value.
- maxdis=9.0 : maximun distance in calculation
- nline=1  : number of lines sampled for given structure pair
- npoints=10 : number of points in each sampled line
- nsamples=2 : number of points picked randomly from trajectories for line sampling
- sstep=0.1 : step size used in sampling
- startdistance= 3.5 : starting distance in sampling
- unit="cm-1" : unit for `max_e`

"""
function create_inputs(cluster1::AbstractString,
                cluster2::AbstractString,
                calculator::Calculator;
                nlines=1,
                nsamples=2,
                max_e=15000,
                unit="cm-1",
                npoints=10,
                maxdis=9.0,
                sstep=0.1,
                startdistance=3.0
               )

    c1 = read_xyz(cluster1)
    c2 = read_xyz(cluster2)

    return  create_inputs(c1, c2, calculator;
                nlines=nlines,
                nsamples=nsamples,
                max_e=max_e,
                unit=unit,
                npoints=npoints,
                maxdis=maxdis,
                sstep=sstep,
                startdistance=startdistance
            )
end


function create_inputs(cluster1::AbstractString,
                cluster2,
                calculator::Calculator;
                nlines=1,
                nsamples=2,
                max_e=15000,
                unit="cm-1",
                npoints=10,
                maxdis=9.0,
                sstep=0.1,
                startdistance=3.0
               )

    c1 = read_xyz(cluster1)

    return  create_inputs(c1, cluster2, calculator;
                nlines=nlines,
                nsamples=nsamples,
                max_e=max_e,
                unit=unit,
                npoints=npoints,
                maxdis=maxdis,
                sstep=sstep,
                startdistance=startdistance
            )
end


function create_inputs(cluster1,
                cluster2::AbstractString,
                calculator::Calculator;
                nlines=1,
                nsamples=2,
                max_e=15000,
                unit="cm-1",
                npoints=10,
                maxdis=9.0,
                sstep=0.1,
                startdistance=3.0
               )

    c2 = read_xyz(cluster2)

    return  create_inputs(cluster1, c2, calculator;
                nlines=nlines,
                nsamples=nsamples,
                max_e=max_e,
                unit=unit,
                npoints=npoints,
                maxdis=maxdis,
                sstep=sstep,
                startdistance=startdistance
            )
end

function create_inputs(cluster1::AbstractCluster,
                cluster2::AbstractCluster,
                calculator::Calculator;
                nlines=1,
                nsamples=2,
                max_e=15000,
                unit="cm-1",
                npoints=10,
                maxdis=9.0,
                sstep=0.1,
                startdistance=3.0
               )


    return  [InputAdaptiveSampler(calculator, cluster1, cluster2, nlines, max_e;
                                  unit=unit, npoints=npoints, maxdis=maxdis,
                                  sstep=sstep, startdistance=startdistance) for _ in 1:nsamples ]
end

function create_inputs(cluster1::AbstractCluster,
                cluster2::AbstractArray{<:AbstractCluster},
                calculator::Calculator;
                nlines=1,
                nsamples=2,
                max_e=15000,
                unit="cm-1",
                npoints=10,
                maxdis=9.0,
                sstep=0.1,
                startdistance=3.0
               )


    return  [InputAdaptiveSampler(calculator, cluster1, rand(cluster2), nlines, max_e;
                                  unit=unit, npoints=npoints, maxdis=maxdis,
                                  sstep=sstep, startdistance=startdistance) for _ in 1:nsamples ]
end


function create_inputs(cluster1::AbstractArray{<:AbstractCluster},
                cluster2::AbstractCluster,
                calculator::Calculator;
                nlines=1,
                nsamples=2,
                max_e=15000,
                unit="cm-1",
                npoints=10,
                maxdis=9.0,
                sstep=0.1,
                startdistance=3.0
               )


    return  [InputAdaptiveSampler(calculator, rand(cluster1), cluster2, nlines, max_e;
                                  unit=unit, npoints=npoints, maxdis=maxdis,
                                  sstep=sstep, startdistance=startdistance) for _ in 1:nsamples ]
end

function create_inputs(cluster1::AbstractArray{<:AbstractCluster},
                cluster2::AbstractArray{<:AbstractCluster},
                calculator::Calculator;
                nlines=1,
                nsamples=2,
                max_e=15000,
                unit="cm-1",
                npoints=10,
                maxdis=9.0,
                sstep=0.1,
                startdistance=3.0
               )


    return  [InputAdaptiveSampler(calculator, rand(cluster1), rand(cluster2), nlines, max_e;
                                  unit=unit, npoints=npoints, maxdis=maxdis,
                                  sstep=sstep, startdistance=startdistance) for _ in 1:nsamples ]
end


"""
    _sample_and_calculate(inputs; pchannel=undef )

Internal helper function to help implement distributed calculations
"""
function _sample_and_calculate(inputs::InputAdaptiveSampler; pchannel=undef)
    start_dir = pwd()
    if typeof(inputs.cal.calculator) == Orca
        if inputs.cal.calculator.tmp_dir != start_dir
            cd(inputs.cal.calculator.tmp_dir)
        end
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
    return bsse_corrected_energy(cal, c1_points, c2_points, basename="base-$(myid())",
                                 id="Pid $(myid())", pchannel=pchannel)
end


function _calculate_points(cal::Calculator{Orca}, c1_points, c2_points; pchannel=undef)
    start_dir = pwd()
    if cal.calculator.tmp_dir != start_dir
        cd(cal.calculator.tmp_dir)
    end
    return bsse_corrected_energy(cal, c1_points, c2_points, basename="base-$(myid())",
                                 id="Pid $(myid())", pchannel=pchannel)
end



"""
    calculate_energy_for_xyzfile(fname, cal; pbar=true)

Reads xyz-file and calculates energy of each point on it
"""
function calculate_energy_for_xyzfile(fname, cal; pbar=true)
    println("Calculating points for a file $(fname)")
    points = read_xyz(fname)
    out=@showprogress pmap(points) do x
        _calculate_energy(cal,x)
    end
    return out
end




function _calculate_energy(cal, points; pchannel=undef)
    return calculate_energy(cal, points, basename="base-$(myid())",
                                 id="Pid $(myid())", pchannel=pchannel)
end

function _calculate_energy(cal::Calculator{Orca}, points; pchannel=undef)
    start_dir = pwd()
    if cal.calculator.tmp_dir != start_dir
        cd(cal.calculator.tmp_dir)
    end
    return calculate_energy(cal, points, basename="base-$(myid())",
                                 id="Pid $(myid())", pchannel=pchannel)
end


function parallel_work(jobs, results, f)
    while true
        i,t = take!(jobs)
        x = f(t...)
        put!(results, (i,x))
    end
end

end  # module restarttools
