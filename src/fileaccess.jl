"""
module fileaccess

This module contains low level fileaccess methods
"""
module fileaccess



export load_jld_data,
       read_xyz,
       save_jld_data


using JLD

using ..atoms
using ..clusters
using ..molecules



"""
save_jld_data(fname::AbstractString, data::Dict)

Saves given data to file in jld format.

Following keys are scanned from `data` and saved if present
- Method
- Basis
- cluster1
- cluster2
- Points
- Energy
- restart_energy
"""
function save_jld_data(fname::AbstractString, data::Dict)
    k = keys(data)
    jldopen(fname, "w") do file
        if haskey(data,"Method")
            file["Method"] = data["Method"]
            @info "Method information saved"
        end
        if haskey(data,"Basis")
            file["Basis"] = data["Basis"]
            @info "Basis informaiton saved"
        end
        if haskey(data,"cluster1")
            file["cluster1"] = data["cluster1"]
            @info "cluster1 information saved"
        end
        if haskey(data,"cluster2")
            file["cluster2"] = data["cluster2"]
            @info "cluster2 information saved"
        end
        if haskey(data,"Points")
            file["Points"] = data["Points"]
            @info "Points information saved"
        end
        if haskey(data,"Energy")
            file["Energy"] = data["Energy"]
            @info "Energy information saved"
        end
        if haskey(data,"restart_energy")
            file["restart_energy"] = data["restart_energy"]
            @info "Restart information saved"
        end
    end
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
