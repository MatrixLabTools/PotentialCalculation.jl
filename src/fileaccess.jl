"""
module fileaccess

This module contains low level fileaccess methods
"""
module fileaccess



export load_jld_data,
       read_xyz,
       save_jld_data


using JLD, FileIO
using ..clusters, ..atoms, ..molecules



"""
save_jld_data(fname, data::Dict)

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
function save_jld_data(fname, data::Dict)
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
load_jld_data(fname)

Loads data from file and returns it as a [`Dict`](@ref)
"""
function load_jld_data(fname)
    data = load(fname)
    function test_old(data)
        if haskey(data,"cluster1") && length(data["cluster1"]) != length(data["cluster1"].atoms)
            return true
        elseif haskey(data,"cluster2") && length(data["cluster2"]) != length(data["cluster2"].atoms)
            return true
        elseif haskey(data,"Points") && length(data["Points"][1]) != length(data["Points"][1].atoms)
            return true
        else
            return false
        end
    end
    if test_old(data)
        @warn "Loaded file has old type data"
        _newform(x) = typeof(x)(x.xyz',x.atoms)
        old_data = deepcopy(data)
        if haskey(data,"cluster1")
            @info "changing cluster1"
            data["cluster1"] = _newform(data["cluster1"])
        end
        if haskey(data,"cluster2")
            @info "changing cluster2"
            data["cluster2"] = _newform(data["cluster2"])
        end
        if haskey(data,"Points")
            @info "changing Points"
            data["Points"] = map(x -> _newform(x), data["Points"])
        end
        if ! test_old(data)
            @warn "Data changed to new form"
            @warn "You should consider saving data to new form using 'save_jld_data'"
            return data
        else
            return old_data
        end
    end
    return data
end



"""
read_xyz(fname)

Reads xyz file and returns it as an Array of [`Cluster`](@ref)
"""
function read_xyz(fname)
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
        push!(clusters, deepcopy(Cluster{AtomOnlySymbol}(xyz, atoms)))
    end
    return clusters
end




end #module
