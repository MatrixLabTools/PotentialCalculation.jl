"""
module Fileaccess

This module contains low level fileaccess methods
"""
module Fileaccess



export load_jld_data
export read_xyz
export save_jld_data


#using JLD
#using JLD2
using FileIO

using ..Atoms
using ..Clusters
using ..Molecules



"""
    save_jld_data(fname::AbstractString, data::Dict)

Saves given data to file in jld format.
"""
function save_jld_data(fname::AbstractString, data::Dict)
    save(fname, data)
    @info "Data writing to file \"$(fname)\" done"
end


"""
    load_jld_data(fname::AbstractString)

Loads data from file and returns it as a [`Dict`](@ref)
"""
function load_jld_data(fname::AbstractString)
    data = load(fname)
    return data
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
