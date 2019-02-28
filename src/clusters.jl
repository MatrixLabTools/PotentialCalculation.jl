module clusters

export AbstractCluster,
       AbstractClusterWithSymbols,
       center_cluster!,
       center_coordinates,
       Cluster,
       cluster_angle,
       ClusterNoSymbols,
       dihedral_angle,
       distances,
       move!,
       print_xyz,
       rotate_x!,
       rotate_y!,
       rotate_z!


using ..atoms
using Distances: Euclidean, pairwise, euclidean
using LinearAlgebra

import Base.==, Base.+



abstract type AbstractCluster end
abstract type AbstractClusterWithSymbols <: AbstractCluster end

"""
    ClusterNoSymbols <: AbstractCluster

Basic cluster has only location of atoms.

# Fields
- `xyz::Array{Float64,2}` : location of atoms in 3d space, first index is x, y, z coordinate
"""
mutable struct ClusterNoSymbols <: AbstractCluster
    xyz::Array{Float64,2}
    function ClusterNoSymbols(xyz::AbstractArray{<:Real,2})
        if size(xyz,1) != 3
            throw(DimensionMismatch("ClusterNoSymbols - xyz has wrong dimensions"))
        end
        new(xyz)
    end
    function ClusterNoSymbols(xyz::AbstractArray{<:Real,1})
        if length(xyz) != 3
            throw(DimensionMismatch("ClusterNoSymbols - atoms has wrong dimensions length=$(size(xyz))"))
        end
        new(reshape(xyz,3,1))
    end
end


"""
    Cluster{T<:AbstractAtom} <: AbstractClusterWithSymbols

Structure to hold location data of clusters/molecules

# Fields
- `xyz::Array{Float64,2}` : location of atoms in 3d space, first index is x, y, z coordinate
- `atoms::Vector{T}` : atom type information
"""
mutable struct Cluster{T<:AbstractAtom} <: AbstractClusterWithSymbols
    "Location of atoms"
    xyz::Array{Float64,2}
    "Symbols for atoms"
    atoms::Vector{T}
    function Cluster{T}(xyz::AbstractArray{<:Real,2}, atoms:: Vector{<:T}) where T<:AbstractAtom
        if size(xyz,2) != length(atoms)
            throw(DimensionMismatch("Cluster has different sizes for atoms $(size(atoms)) and xyz $(size(xyz))"))
        elseif size(xyz,1) != 3
            throw(DimensionMismatch("Cluster - xyz has wrong dimensions"))
        end
        new(xyz,atoms)
    end
    function Cluster{T}(xyz::AbstractArray{<:Real,1}, atoms:: Vector{<:T}) where T<:AbstractAtom
        if length(atoms) != 1
            throw(DimensionMismatch("Cluster has different number for atoms $(length(atoms)) and xyz 1"))
        elseif length(xyz) != 3
            throw(DimensionMismatch("Cluster - xyz has wrong dimensions"))
        end
        new(reshape(xyz,3,1),atoms)
    end
    function Cluster{T}(xyz::AbstractArray{<:Real,1}, atom::T) where T<:AbstractAtom
        if length(xyz) != 3
            throw(DimensionMismatch("Cluster - xyz has wrong dimensions size=$(size(xyz))"))
        end
        new(reshape(xyz,3,1), [atom])
    end
end



Base.show(io::IO, C::AbstractCluster) = print(io,"$(typeof(C)) of ", size(C), " atoms")

Base.size(a::AbstractCluster) = size(a.xyz,2)

Base.length(a::AbstractCluster) = size(a)

Base.lastindex(a::AbstractCluster) = length(a)


function (+)(c1::T, c2::T) where T <: AbstractClusterWithSymbols
    return T(hcat(c1.xyz,c2.xyz),vcat(c1.atoms,c2.atoms))
end

function (+)(c1::ClusterNoSymbols, c2::ClusterNoSymbols)
    return ClusterNoSymbols(hcat(c1.xyz,c2.xyz))
end


function Base.getindex(C::T, i::Int) where T <: AbstractClusterWithSymbols
    return T(C.xyz[:,i], [C.atoms[i]])
end

function Base.getindex(C::ClusterNoSymbols, i::Int)
    return ClusterNoSymbols(C.xyz[:,i])
end


function Base.getindex(C::T, i::AbstractUnitRange) where T <: AbstractClusterWithSymbols
    return T(C.xyz[:,i], C.atoms[i])
end

function Base.getindex(C::ClusterNoSymbols, i::AbstractUnitRange)
    return ClusterNoSymbols(C.xyz[:,i])
end


function Base.print(io::IO, C::AbstractClusterWithSymbols)
    for i in 1:size(C)
        println(io, C.atoms[i].id, "    ",C.xyz[1,i], "  ", C.xyz[2,i], "  ",C.xyz[3,i])
    end
end


"""
    distances(c1::AbstractCluster, c2::AbstractCluster)

Return distances between atoms of given clusters
"""
function distances(c1::AbstractCluster, c2::AbstractCluster)
    return pairwise(Euclidean(),c1.xyz,c2.xyz,dims=2)
end

"""
    distances(c::AbstractCluster, ar1::AbstractRange, ar2::AbstractRange)

Returns distances between atoms in given unit ranges
"""
function distances(c::AbstractCluster, ur1::AbstractUnitRange, ur2::AbstractUnitRange)
    return pairwise(Euclidean(),c.xyz[:,ur1],c.xyz[:,ur2],dims=2)
end

"""
    distances(c::AbstractCluster, i, j)

Returns distance of atoms `i` and `j`
"""
function distances(c::AbstractCluster, i, j)
    return euclidean(c.xyz[:,i], c.xyz[:,j])
end

"""
    distances(c1::AbstractCluster, i, c2::AbstractCluster, j)

Return distance between atom `i` in `c1` and atom `j` in `c2`
"""
function distances(c1::AbstractCluster, i, c2::AbstractCluster, j)
    return euclidean(c1.xyz[:,i], c2.xyz[:,j])
end

"""
    distances(c1::AbstractCluster, ur1::AbstractUnitRange,
              c2::AbstractCluster, ur2::AbstractUnitRange)

Return distance between atoms `ur1` in `c1` and atom `ur2` in `c2`
"""
function distances(c1::AbstractCluster, ur1::AbstractUnitRange,
                   c2::AbstractCluster, ur2::AbstractUnitRange)
    return pairwise(Euclidean(),c1.xyz[:,ur1],c2.xyz[:,ur2],dims=2)
end

"""
    center_coordinates(c::AbstractCluster)

Gives coordinates to aricmetric mean of clusters atoms
"""
function center_coordinates(c::AbstractCluster)
    return sum(c.xyz, dims=2) / length(c)
end

"""
    move!(c::AbstractCluster,r)

Moves cluster by `r`
"""
function move!(c::AbstractCluster,r)
    c.xyz .+= r
end

"""
    center_cluster!(c::AbstractCluster)

Centers cluster to origin of coordinates
"""
function center_cluster!(c::AbstractCluster)
    move!(c, -center_coordinates(c))
end


"""
    rotate_x!(c::AbstractCluster, θ)

Rotates cluster around x-axis by angle `θ`
"""
function rotate_x!(c::AbstractCluster, θ)
    R = Matrix{Float64}(I, 3,3)
    R[2,2] = cos(θ)
    R[3,3] = R[2,2]
    R[3,2] = -sin(θ)
    R[2,3] = -R[3,2]
    c.xyz = R * c.xyz
end


"""
    rotate_y!(c::AbstractCluster, θ)

Rotates cluster around y-axis by angle `θ`
"""
function rotate_y!(c::AbstractCluster, θ)
    R = Matrix{Float64}(I, 3,3)
    R[1,1] = cos(θ)
    R[3,3] = R[1,1]
    R[3,1] = sin(θ)
    R[1,3] = -R[3,1]
    c.xyz = R * c.xyz
end


"""
    rotate_z!(c::AbstractCluster, θ)

Rotates cluster around z-axis by angle `θ`
"""
function rotate_z!(c::AbstractCluster, θ)
    R = Matrix{Float64}(I, 3,3)
    R[1,1] = cos(θ)
    R[2,2] = R[1,1]
    R[2,1] = -sin(θ)
    R[1,2] = -R[2,1]
    c.xyz = R * c.xyz
end


"""
    print_xyz(io::IO, c::AbstractClusterWithSymbols, note=""; printheader=true)

Prints cluster in xyz file format

# Arguments
- `io::IO` : stream where writing is done
- `c::AbstractClusterWithSymbols` : cluster that is writen
- `note=""` : message writen on note line
- `printheader=true` : wheather or not header is writen (number of atoms and note)
"""
function print_xyz(io::IO, c::AbstractClusterWithSymbols, note=""; printheader=true)
    if printheader
        println(io, "    ",length(c))
        println(io, note)
    end
    for i in 1:length(c)
        println(io, c.atoms[i].id, "   ", c.xyz[1,i], "  ", c.xyz[2,i], "  ", c.xyz[3,i])
    end
end



"""
cluster_angle(c::AbstractCluster, i, j, k)

Calculates angle (radians) between atoms i,j,k in cluster
"""
function cluster_angle(c::AbstractCluster, i, j, k)
    r1 = c.xyz[:,i] - c.xyz[:,j]
    r2 = c.xyz[:,k] - c.xyz[:,j]
    return acos(dot(r1,r2)/sqrt(dot(r1,r1)*dot(r2,r2)))
end


"""
cluster_angle(c1::AbstractCluster, i, j, c2::AbstractCluster, k)

Calculates angle (radians) between atons in different clusters

# Arguments
- `c1::AbstractCluster` : first cluster
- `i` : index in `c1`
- `j` : index in `c1`
- `c2::AbstractCluster` : second cluster
- `k` : index in `c2`
"""
function cluster_angle(c1::AbstractCluster, i, j, c2::AbstractCluster, k)
    r1 = c1.xyz[:,i] - c1.xyz[:,j]
    r2 = c2.xyz[:,k] - c1.xyz[:,j]
    return acos(dot(r1,r2)/sqrt(dot(r1,r1)*dot(r2,r2)))
end


"""
cluster_dihedral_angle(c::AbstractCluster, i, j, k, m)
"""
function dihedral_angle(c::AbstractCluster, i, j, k, m)
    r1 = c.xyz[:,j] - c.xyz[:,i]
    r2 = c.xyz[:,k] - c.xyz[:,j]
    r3 = c.xyz[:,m] - c.xyz[:,k]
    t1 = cross(cross(r1,r2), cross(r2,r3))
    t2 = dot(cross(r1,r2), cross(r2,r3))
    return atan( dot(t1, r2./norm(r2)), t2 )
end

end #module
