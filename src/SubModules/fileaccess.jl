"""
module Fileaccess

This module contains low level fileaccess methods
"""
module Fileaccess



export load_jld_data
export read_xyz
export save_jld_data


using FileIO

using ..Atoms
using ..Clusters
using ..Molecules



"""
    save_jld_data(fname::AbstractString, data::Dict)

Saves given data to file in jld format.
"""
function save_jld_data(fname::AbstractString, data::Dict)
    # Make saved data independed of Cluster type
    points = data["Points"]
    xyz = [ p.xyz for p in points ]
    symbols = map(x->x.id, points[1].atoms)
    cluster1 = 1:length(data["cluster1"])
    cluster2 = last(cluster1)+1 : last(cluster1) + length(data["cluster2"])

    new_data = filter( data ) do (key,val)
        ! ( key in ["Points", "cluster1", "cluster2"] )
    end
    new_data["cluster1"] = cluster1
    new_data["cluster2"] = cluster2
    new_data["symbols"] = symbols
    new_data["xyz"] = xyz
    save(fname, data)
    @info "Data writing to file \"$(fname)\" done"
end


"""
    load_jld_data(fname::AbstractString)

Loads data from file and returns it as a [`Dict`](@ref)
"""
function load_jld_data(fname::AbstractString)
    data = load(fname)
    if haskey(data, "Points") # Old format
        # Has Cluster data in the save file
        return data
    else
        # Construct clusters
        atoms = AtomOnlySymbol.(data["symbols"])
        points = [ Cluster(xyz, atoms) for xyz in data["xyz"] ]
        cluster1 = points[1][data["cluster1"]]
        cluster2 = points[1][data["cluster2"]]

        new_data = filter( data ) do (key,val)
            ! ( key in ["xyz", "cluster1", "cluster2", "symbols"] )
        end
        new_data["Points"] = points
        new_data["cluster1"] = cluster1
        new_data["cluster2"] = cluster2
        return new_data
    end
end


"""
read_xyz(fname::AbstractString)

Reads xyz file and returns it as an Array of [`Cluster`](@ref)
"""
function read_xyz(fname::AbstractString)
    lines = Vector{String}()
    open(fname, "r") do file
        lines = readlines(file)
    end

    # How many atoms
    natoms = parse(Int, lines[1])

    # How many clusters - use of floor allows extra empty lines at end
    nclusters = Int(floor(length(lines)/(natoms+2)))
    @debug "Type of nclusters $(typeof(nclusters))"

    xyz = zeros(Float64, 3, natoms)
    atoms = Vector{AtomOnlySymbol}(undef, natoms)
    clusters = Vector{Cluster{AtomOnlySymbol}}()

    for nc in 1:nclusters
        for na in 1:natoms
            cont = split(lines[(nc-1)*(natoms+2)+na+2])
            atoms[na] = AtomOnlySymbol(cont[1])
            xyz[1,na] = parse(Float64, cont[2])
            xyz[2,na] = parse(Float64, cont[3])
            xyz[3,na] = parse(Float64, cont[4])
        end
        push!(clusters, deepcopy(Cluster(xyz, atoms)))
    end
    return clusters
end




end #module
