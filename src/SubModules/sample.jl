module Sample

export line_sampler, adaptive_line_sampler, sample_multiple_adaptive_lines,
       InputAdaptiveSampler

using LinearAlgebra
using ProgressMeter
using Rotations
using ..Atoms
using ..Calculators
using ..Clusters
using ..UnitConversions



"""
    _line_sampler(cl1::Cluster, cl2::Cluster, u, step, npoints)

Local helper function.

Returns clusters tuple of cluster arrays that hold positions where clusters
are on different distances from each other in a line like fashion.

# Arguments
- `cl1::Cluster` : first cluster and it's initial location
- `cl2::Cluster` : second cluster and it's initial location
- `u`  : (unit) vector descriping  direction on which `cl2` is moved
- `step` : distance between points
- `npoints` : nunber of points sampled
"""
function _line_sampler(cl1::Cluster, cl2::Cluster, u, step, npoints)
    c1 = deepcopy(cl1)
    c2 = deepcopy(cl2)
    out1 = Array{typeof(c1)}(undef, npoints)
    out2 = Array{typeof(c2)}(undef, npoints)
    for i in eachindex(out1)
        out1[i] = deepcopy(c1)
        out2[i] = deepcopy(c2)
        move!(c2, step*u)
    end
    return out1, out2
end


"""
    line_sampler(cl1::Cluster, cl2::Cluster; npoints=10, mindis=2.0, maxdis=9.0)

Samples cluster position in a line like fashion with even intervals.

Returns clusters tuple of cluster arrays that hold positions where clusters
are on different distances from each other in a line like fashion.

Sampling is done by first generating random direction and placing `cl2` in that
direction from `cl2` at `mindis` and then picking `npoints` in even distances
to `maxdis`.
"""
function line_sampler(cl1::Cluster, cl2::Cluster; npoints=10, mindis=2.0, maxdis=9.0)
    c1 = deepcopy(cl1)
    c2 = deepcopy(cl2)

    center_cluster!(c1)
    center_cluster!(c2)

    rotate_randomly!(c1)
    rotate_randomly!(c2)

    step = (maxdis - mindis) / npoints
    u = [1.0, 0.0, 0.0]
    move!(c2, mindis * u )
    while minimum(distances(c1,c2)) < mindis
        move!(c2, step * u)
    end
     return _line_sampler(c1,c2,u,step,npoints)
end


"""
  adaptive_line_sampler(cal::Calculator, cl1::Cluster, cl2::Cluster, max_e=0; unit="cm-1",
               npoints=10, maxdis=9.0, sstep=0.1, startdistance=3.0, pchannel=undef)

Calculates potential on a line type distances that is suitable for visualization.
Closest distance used in calculations is looked by searching closest distance that
has energy lower than `max_e`

# Arguments
- `max_e`           : Point that is closest and has less energy than this will be starting point for a line

# Keywords
- `unit`            : Unit in which `max_e` is given
- `npoints`         : Number of points in potential
- `maxdis`          : Maximum distance in potential calculation
- `sstep`           : Search step size for adaptive search, also defines accuracy of search
- `startdistance`   : Distance from which adaptive seach is started
- `basename="base"` : prefix for temporary files in calculations
- `id=""`           : extra identificaltion to help to do parallel computations
- `pchannel=undef`  : (Remote)Channel where progess information is added
"""
function adaptive_line_sampler(cal::Calculator, cl1::Cluster, cl2::Cluster, max_e=0;
                               unit="cm-1", npoints=10, maxdis=9.0, sstep=0.1, startdistance=3.0,
                               basename="base", id="", pchannel=undef)
    @debug "Starting adaptive_line_sampler"

    # take deepcopys so that originals are not moved
    c1 = deepcopy(cl1)
    c2 = deepcopy(cl2)

    center_cluster!(c1)
    center_cluster!(c2)

    # Rotate so that we get random orientation
    rotate_randomly!(c1)
    rotate_randomly!(c2)

    # We move the clusters in x-direction
    u = [1.0, 0.0, 0.0]

    # Find good initial point for the search
    # TODO make this more robust
    r = 0.5*startdistance
    move!(c2, 0.5*startdistance*u)
    while minimum(distances(c1,c2)) < startdistance
        move!(c2, sstep * u)
        r += sstep
    end

    # Search for the intial point
    # We need to find a point where potential energy is emax
    # To do this we find the most distant point where energy is more than emax
    # and the closest point where energy is less than e
    emax = energy_from(max_e, unit)
    Δe = bsse_corrected_energy(cal, [c1], [c2], basename=basename, id=id, pchannel=pchannel) .- emax
    ctemp = [deepcopy(c2)] # History for points we tried
    fless = false # We have a point where Δe < 0
    fmore = false # We have a point where Δe > 0
    if Δe[1] > 0
        fmore = true
    else
        fless = true
    end
    while fmore == false || fless == false
        if Δe[end] < 0
            move!(c2, -1*sstep*u)
            r -= sstep
        else
            move!(c2, sstep*u)
            r += sstep
        end
        tmp = bsse_corrected_energy(cal, [c1], [c2], basename=basename, id=id, pchannel=pchannel) .- emax
        push!(Δe,  tmp...)
        push!(ctemp, deepcopy(c2))
        if Δe[end] > 0
            fmore = true
        else
            fless = true
        end
    end
    if length(Δe) == 1
        error("Search fail!")
    end

    # Pretare output and perform some checks
    @debug "adaptive_line_sampler found initial point in $(length(e)) steps"
    out2 = [ctemp[Δe.<0][end]] # We take the closest point where e < emax
    eout = [Δe[Δe.<0][end] + emax] # Calculated potential energy
    out1 = [c1]  # cluster1 coordinates
    c2 = deepcopy(out2[1])
    rr = minimum(distances(c1,c2))
    @debug "Minimum distance at start" r rr
    @assert maxdis > rr "You need to increase maximum distance in search"

    # Lets now find the points
    step = (maxdis - rr ) / npoints
    move!(c2,step*u) # We already calculated the first one, so we need to move
    tmp = _line_sampler(c1, c2, u, step, npoints-1)

    # Now we can calculate the potential
    et = bsse_corrected_energy(cal, tmp[1], tmp[2], basename=basename, id=id, pchannel=pchannel)
    push!(eout, et...)
    push!(out1, tmp[1]...)
    push!(out2, tmp[2]...)
    @debug "Minimum distance at end $(minimum(distances(out1[end],out2[end])))"
    return Dict("Energy" => eout, "Points" =>  out1 .+ out2,  "Mindis" =>  minimum(distances(out1[1],out2[1])))
end


"""
    InputAdaptiveSampler

Structure used to help use of adaptive line samplers

# Fields
- `cal::Calculator`   : [`calculator`](@ref) used in calculations
- `cl1::Cluster`      : first cluster
- `cl2::Cluster`      : second cluster
- `nlines`            : number of lines calculated
- `max_e`             : maximum energy for configuration
- `unit=`             : unit for maximum energy ( see [`energy_from`](@ref) for supported units)
- `npoints`           : number of points in line
- `maxdis`            : maximum cluster distance
- `sstep`             : search step used by [`adaptive_line_sampler`](@ref)
- `startdistance`     : distance where search is started
"""
mutable struct InputAdaptiveSampler
    # This is needed for parallelization.
    # pmap needs that all data is send to othe process
    # wich is intended to be done with this struct

    # TODO clean this - it is way too messy now
    cal
    cl1
    cl2
    nlines
    max_e
    unit
    npoints
    maxdis
    sstep
    startdistance
    function InputAdaptiveSampler(cal::Calculator, cl1::Cluster, cl2::Cluster,
                                  nlines, max_e=0;
                                  unit="cm-1", npoints=10, maxdis=9.0, sstep=0.1,
                                  startdistance=3.0)
        new(cal,cl1,cl2,nlines,max_e,unit,npoints,maxdis,sstep,startdistance)
    end
end


"""
sample_multiple_adaptive_lines(cal::Calculator, cl1::Cluster, cl2::Cluster, nlines, max_e=0;
                               unit="cm-1", npoints=10, maxdis=9.0, sstep=0.1, startdistance=3.0,
                               basename="base", id="")

Calculates energy of two given clusters in different distances of each other.
Distances are sampled in line like fashion in random directions with even distances.
Sampling of line is done with [`adaptive_line_sampler`](@ref).

# Arguments
- `cal::Calculator`   : [`calculator`](@ref) used in calculations
- `cl1::Cluster`      : first cluster
- `cl2::Cluster`      : second cluster
- `nlines`            : number of lines calculated
- `max_e=0`           : maximum energy for configuration
- `unit="cm-1"`       : unit for maximum energy ( see [`energy_from`](@ref) for supported units)
- `npoints=10`        : number of points in line
- `maxdis=9.0`        : maximum cluster distance
- `sstep=0.1`         : search step used by [`adaptive_line_sampler`](@ref)
- `startdistance=3.0` : distance where search is started
- `basename="base"` : prefix for temporary files in calculations
- `id=""`           : extra identificaltion to help to do parallel computations
- `pchannel=undef`                : (Remote)Channel where progess information is added

# Returns
[`Dict`](@ref) with keys
* "Energy" : array holding energy
* "Points" : array of [`Cluster`](@ref) representing points where energy was calculated
* "Mindis" : array of minimum distances in each line
"""
function sample_multiple_adaptive_lines(cal::Calculator, cl1::Cluster, cl2::Cluster, nlines, max_e=0;
                               unit="cm-1", npoints=10, maxdis=9.0, sstep=0.1, startdistance=3.0,
                               basename="base", id="", pchannel=undef)
    c1 = deepcopy(cl1)
    c2 = deepcopy(cl2)
    out = Vector{Dict}(undef, nlines)
    rtmp = Float64[]
    sr = startdistance
    for i in 1:nlines
        tmp = adaptive_line_sampler(cal, c1, c2, max_e; unit=unit, npoints=npoints,
                                    maxdis=maxdis, startdistance=sr, basename=basename,
                                     id=id, pchannel=pchannel)
        push!(rtmp, tmp["Mindis"])
        sr = sum(rtmp) / length(rtmp)
        out[i] = tmp
    end
    energy = hcat(map( x -> x["Energy"], out)...)
    points = hcat(map( x -> x["Points"], out)...)
    mindis = map(x -> x["Mindis"], out)
    return Dict("Energy" => energy, "Points" => points, "Mindis" => mindis)
end


"""
sample_multiple_adaptive_lines(inputs::InputAdaptiveSampler; basename="base", id="", pchannel=undef)

Calculates energy of two given clusters in different distances of each other.
Distances are sampled in line like fashion in random directions with even distances.
Sampling of line is done with [`adaptive_line_sampler`](@ref).

# Arguments
- `inputs::InputAdaptiveSampler`  : input information in clean package
- `basename="base"`               : prefix for temporary files in calculations
- `id=""`                         : extra identificaltion to help to do parallel computations
- `pchannel=undef`                : (Remote)Channel where progess information is added

# Returns
[`Dict`](@ref) with keys
* "Energy" : array holding energy
* "Points" : array of [`Cluster`](@ref) representing points where energy was calculated
* "Mindis" : array of minimum distances in each line
"""
function sample_multiple_adaptive_lines(inputs::InputAdaptiveSampler; basename="base", id="", pchannel=undef)
    sample_multiple_adaptive_lines(inputs.cal, inputs.cl1, inputs.cl2,
                 inputs.nlines, inputs.max_e,
                 unit=inputs.unit, npoints=inputs.npoints, maxdis=inputs.maxdis,
                 startdistance=inputs.startdistance, basename=basename, id=id, pchannel=pchannel)
end

end  # module sample
