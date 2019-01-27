module fileaccess
# low level fileaccesss


export read_h5file, give_radius, give_energy, save_jld_data, load_jld_data, read_xyz,
       make_savable_data
using HDF5, JLD, FileIO
using ..clusters, ..atoms, ..molecules


# Old form function
function convert_to_array(atoms::String)
    tmp = split(replace(replace(replace(atoms,"'"=>""),
                                              "]"=>""),
                                              "["=>""),  ", ")
    return convert.(String,tmp)
end

# Old form function
function toclusters(data::AbstractArray)
    s = size(data)
    tmp = reshape(data,s[1],s[2],prod(s[3:end]))
    out = [ ClusterNoSymbols(tmp[:,:,i]) for i in axes(tmp)[3] ]
    return reshape(out,s[3:end])
end

# This is old form function
function read_h5file(fname::String)
    h5open(fname,"r") do fdata
        @info "Opening file $(fname) for reading"
        g = fdata["geometry"]
        e = fdata["energy"]
        c1 = toclusters(read(g,"c1"))
        c2 = toclusters(read(g,"c2"))
        points = toclusters(read(g,"points"))
        typeofcal = read(attrs(fdata),"type")
        @info "File data type is: \"$(typeofcal)\""
        c1_atoms = Molecule{AtomOnlySymbol}( convert_to_array(read(attrs(fdata),"c1_atoms")) )
        c2_atoms = Molecule{AtomOnlySymbol}( convert_to_array(read(attrs(fdata),"c2_atoms")) )
        energy = read(e,"points")
        c1_energy = read(e,"c1")
        c2_energy = read(e,"c2")
        out = Dict("c1"=>c1,
                "c2"=>c2,
                "points"=>points,
                "typeofcal"=>typeofcal,
                "c1_molecule"=>c1_atoms,
                "c2_molecule"=>c2_atoms,
                "energy"=>energy)
        if exists(e, "c1_energy")
            @info "File has energy data for each fragment"
            c1_energy = read(e,"c1")
            c2_energy = read(e,"c2")
            push!(out,"c1_energy"=>c1_energy, "c2_energy"=>c2_energy)
        end
        if exists(e,"bsse1")
            @info "File has BSSE correction data for fragments"
            bsse1 = read(e,"bsse1")
            bsse2 = read(e,"bsse2")
            push!(out, "bsse1"=>bsse1, "bsse2"=>bsse2)
        end
        return out
    end
end

# This is old form
function give_radius(data, flat=false, T=Array)
    l1 = length(data["c1_molecule"])
    l2 = l1+length(data["c2_molecule"])
    s2 = l1+1
    @debug "Bounds are 1:$(l1) and $(s2):$(l2)"
    out = similar(data["points"],T{Float64,1})
    for i in eachindex(out)
        out[i] = vec(distances(data["points"][i],1:l1, s2:l2))
    end
    if flat
        return T(hcat(out...)')
    else
        return T(out)
    end
end


function give_energy(data, flat=false)
    if haskey(data, "bsse1") && haskey(data,"bsse2")
        out = data["energy"] .- data["bsse1"] .- data["bsse2"]
    else
        out = data["energy"]
    end
    if flat
        return vcat(out...)
    else
        return out
    end
end


function save_jld_data(fname, data)
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


function load_jld_data(fname)
    return load(fname)
end



# This is old form
function make_savable_data(inp, data)
    e = hcat( [data[i][j]["Energy"] for i in 1:length(data) for j in 1:length(data[1])  ]...  )
    p = hcat( [data[i][j]["Points"] for i in 1:length(data) for j in 1:length(data[1])  ]...  )
    return Dict("Points"=>p, "Energy"=>e, "c1"=>inp.cl1, "c2"=>inp.cl2,
                "Basis" => inp.cal.basis, "Method" => inp.cal.method)
end

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

    xyz = zeros(Float64, natoms,3)
    atoms = Vector{AtomOnlySymbol}(undef, natoms)
    clusters = Vector{Cluster{AtomOnlySymbol}}()

    for nc in 1:nclusters
        for na in 1:natoms
            cont = split(lines[(nc-1)*(natoms+2)+na+2])
            atoms[na] = AtomOnlySymbol(cont[1])
            xyz[na,1] = parse(Float64, cont[2])
            xyz[na,2] = parse(Float64, cont[3])
            xyz[na,3] = parse(Float64, cont[4])
        end
        push!(clusters, Cluster{AtomOnlySymbol}(xyz, atoms))
    end
    return clusters
end




end #module
