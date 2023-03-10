module Clusters

export AbstractCluster
export AbstractClusterWithSymbols
export center_cluster!
export center_coordinates
export Cluster
export cluster_angle
export ClusterNoSymbols
export dihedral_angle
export distances
export move!
export print_xyz
export rotate_x!
export rotate_y!
export rotate_z!
export rotate_randomly!


using AtomsBase
using ..Atoms
using Distances: Euclidean, pairwise, euclidean
using ExtXYZ
using LinearAlgebra
using Rotations
using StaticArrays
using Unitful
using UnitfulAtomic

import Base.==, Base.+



abstract type AbstractCluster end
abstract type AbstractClusterWithSymbols{T} <: AbstractCluster where T<:AbstractAtom end

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
    Cluster{T} <: AbstractClusterWithSymbols{T} where T<:AbstractAtom

Structure to hold location data of clusters/molecules

# Fields
- `xyz::Array{Float64,2}` : location of atoms in 3d space, first index is x, y, z coordinate
- `atoms::Vector{T}` : atom type information
"""
mutable struct Cluster{T} <: AbstractClusterWithSymbols{T}
    "Location of atoms"
    xyz::Array{Float64,2}
    "Symbols for atoms"
    atoms::Vector{T}
    function Cluster(xyz::AbstractArray{<:Real,2}, atoms::AbstractVector{T}) where T<:AbstractAtom
        if size(xyz,2) != length(atoms)
            throw(DimensionMismatch("Cluster has different sizes for atoms $(size(atoms)) and xyz $(size(xyz))"))
        elseif size(xyz,1) != 3
            throw(DimensionMismatch("Cluster - xyz has wrong dimensions"))
        end
        new{T}(xyz,atoms)
    end
    function Cluster(xyz::AbstractArray{<:Real,1}, atoms::AbstractVector{T}) where T<:AbstractAtom
        if length(atoms) != 1
            throw(DimensionMismatch("Cluster has different number for atoms $(length(atoms)) and xyz 1"))
        elseif length(xyz) != 3
            throw(DimensionMismatch("Cluster - xyz has wrong dimensions"))
        end
        new{T}(reshape(xyz,3,1),atoms)
    end
    function Cluster(xyz::AbstractVector{<:Real}, atom::T) where T<:AbstractAtom
        if length(xyz) != 3
            throw(DimensionMismatch("Cluster - xyz has wrong dimensions size=$(size(xyz))"))
        end
        new{T}(reshape(xyz,3,1), [atom])
    end
end



Base.show(io::IO, C::AbstractCluster) = print(io,"$(typeof(C)) of ", size(C), " atoms")

Base.size(a::AbstractCluster) = size(a.xyz,2)

Base.length(a::AbstractCluster) = size(a)

Base.lastindex(a::AbstractCluster) = length(a)


function (+)(c1::Cluster{T}, c2::Cluster{T}) where T
    return Cluster(hcat(c1.xyz,c2.xyz),vcat(c1.atoms,c2.atoms))
end

function (+)(c1::ClusterNoSymbols, c2::ClusterNoSymbols)
    return ClusterNoSymbols(hcat(c1.xyz,c2.xyz))
end


function Base.getindex(C::Cluster, i::Int)
    return Cluster(C.xyz[:,i], [C.atoms[i]])
end

function Base.getindex(C::ClusterNoSymbols, i::Int)
    return ClusterNoSymbols(C.xyz[:,i])
end


function Base.getindex(C::Cluster, i::AbstractUnitRange)
    return Cluster(C.xyz[:,i], C.atoms[i])
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
    return pairwise(Euclidean(),c1.xyz,c2.xyz,dims=2) .* u"Å"
end

"""
    distances(c::AbstractCluster, ar1::AbstractRange, ar2::AbstractRange)

Returns distances between atoms in given unit ranges
"""
function distances(c::AbstractCluster, ur1::AbstractUnitRange, ur2::AbstractUnitRange)
    return pairwise(Euclidean(),c.xyz[:,ur1],c.xyz[:,ur2],dims=2) .* u"Å"
end

"""
    distances(c::AbstractCluster, i, j)

Returns distance of atoms `i` and `j`
"""
function distances(c::AbstractCluster, i, j)
    return euclidean(c.xyz[:,i], c.xyz[:,j]) .* u"Å"
end

"""
    distances(c1::AbstractCluster, i, c2::AbstractCluster, j)

Return distance between atom `i` in `c1` and atom `j` in `c2`
"""
function distances(c1::AbstractCluster, i, c2::AbstractCluster, j)
    return euclidean(c1.xyz[:,i], c2.xyz[:,j]) .* u"Å"
end

"""
    distances(c1::AbstractCluster, ur1::AbstractUnitRange,
              c2::AbstractCluster, ur2::AbstractUnitRange)

Return distance between atoms `ur1` in `c1` and atoms `ur2` in `c2`
"""
function distances(c1::AbstractCluster, ur1::AbstractUnitRange,
                   c2::AbstractCluster, ur2::AbstractUnitRange)
    return pairwise(Euclidean(),c1.xyz[:,ur1],c2.xyz[:,ur2],dims=2) .* u"Å"
end

"""
    center_coordinates(c::AbstractCluster)

Gives coordinates to aricmetric mean of clusters atoms
"""
function center_coordinates(c::AbstractCluster)
    return sum(c.xyz, dims=2) ./ length(c) .* u"Å"
end


Unitful.uconvert(::typeof(u"Å"), r::Real) = r*u"Å"

"""
    move!(c::AbstractCluster,r)

Moves cluster by `r`
"""
function move!(c::AbstractCluster,r)
    tmp = @. uconvert(u"Å", r) |> ustrip
    c.xyz .+= tmp
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
    c.xyz = RotX(θ) * c.xyz
end


"""
    rotate_y!(c::AbstractCluster, θ)

Rotates cluster around y-axis by angle `θ`
"""
function rotate_y!(c::AbstractCluster, θ)
    c.xyz = RotY(θ) * c.xyz
end


"""
    rotate_z!(c::AbstractCluster, θ)

Rotates cluster around z-axis by angle `θ`
"""
function rotate_z!(c::AbstractCluster, θ)
    c.xyz = RotZ(θ) * c.xyz
end


"""
    rotate_randomly!(c::AbstractCluster)

Rotate cluster by random angle and axis
"""
function rotate_randomly!(c::AbstractCluster)
    c.xyz = rand(RotMatrix{3}) * c.xyz
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



## AtomsBase support

function AtomsBase.bounding_box(::AbstractCluster)
    a = SVector{3}( [Inf, 0., 0.] .* u"bohr" )
    b = SVector{3}( [0., Inf, 0.] .* u"bohr" )
    c = SVector{3}( [0., 0., Inf] .* u"bohr" )
    return SVector(a, b, c)
end


function AtomsBase.boundary_conditions(::AbstractCluster)
    return SVector{3, BoundaryCondition}(DirichletZero(), DirichletZero(), DirichletZero())
end


function Cluster(sys::AtomsBase.FlexibleSystem)
    a = AtomOnlySymbol.( atomic_symbol(sys) )
    pos = map( position(sys) ) do r
        ustrip.(u"Å", r)
    end
    return Cluster(hcat(pos...), a)
end


function AtomsBase.FlexibleSystem(c::Cluster; kwargs...)
    a = [
        Symbol(c.atoms[i].id) => SVector{3}( c.xyz[:,i] .*u"Å")
        for i in 1:length(c)
    ]
    return isolated_system(a; kwargs...)
end

end #module Clusters
