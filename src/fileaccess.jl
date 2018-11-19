module fileaccess

export read_h5file, give_radius, give_energy

using HDF5
using ..clusters, ..atoms, ..molecules


function convert_to_array(atoms::String)
    tmp = split(replace(replace(replace(atoms,"'"=>""),
                                              "]"=>""),
                                              "["=>""),  ", ")
    return convert.(String,tmp)
end

function toclusters(data::AbstractArray)
    s = size(data)
    tmp = reshape(data,s[1],s[2],prod(s[3:end]))
    out = [ ClusterNoSymbols(tmp[:,:,i]) for i in axes(tmp)[3] ]
    return reshape(out,s[3:end])
end


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



end #module
