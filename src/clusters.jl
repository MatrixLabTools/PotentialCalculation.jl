module clusters

export AbstractCluster, AbstractClusterWithSymbols,
       ClusterNoSymbols, Cluster, ClusterWithMass, distances,
       center_coordinates, move!, center_cluster!,
       rotate_x!, rotate_y!, rotate_z!, print_xyz

using ..atoms
using Distances: Euclidean, pairwise
using LinearAlgebra

import Base.==, Base.+



abstract type AbstractCluster end
abstract type AbstractClusterWithSymbols <: AbstractCluster end


mutable struct ClusterNoSymbols{} <: AbstractCluster
    xyz::Array{Float64,2}
    function ClusterNoSymbols(xyz::AbstractArray{<:Real,2})
        if size(xyz,1) != 3
            error("ClusterNoSymbols - xyz has wrong dimensions")
        end
        new(xyz)
    end
    function ClusterNoSymbols(xyz::AbstractArray{<:Real,1})
        if length(xyz) != 3
            error("ClusterNoSymbols - atoms has wrong dimensions length=$(size(xyz))")
        end
        new(reshape(xyz,3,1))
    end
end


"""
    Cluster{T<:AbstractAtom} <: AbstractClusterWithSymbols

Structure to hold location data of clusters/molecules
"""
mutable struct Cluster{T<:AbstractAtom} <: AbstractClusterWithSymbols
    "Location of atoms"
    xyz::Array{Float64,2}
    "Symbols for atoms"
    atoms::Vector{T}
    function Cluster{T}(xyz::AbstractArray{<:Real,2}, atoms:: Vector{<:T}) where T<:AbstractAtom
        if size(xyz,2) != length(atoms)
            error("Cluster has different sizes for atoms $(size(atoms)) and xyz $(size(xyz))")
        elseif size(xyz,1) != 3
            error("Cluster - xyz has wrong dimensions")
        end
        new(xyz,atoms)
    end
    function Cluster{T}(xyz::AbstractArray{<:Real,1}, atoms:: Vector{<:T}) where T<:AbstractAtom
        if length(xyz) != 3
            error("Cluster - xyz has wrong dimensions size=$(size(xyz))")
        end
        if length(atoms) != 1
            error("Cluster - atoms has wrong dimensions length=$(size(xyz))")
        end
        new(reshape(xyz,3,1), atoms)
    end
end


mutable struct ClusterWithMass{T<:AbstractAtom} <: AbstractClusterWithSymbols
    xyz::Array{Float64,2}
    atoms::Vector{T}
    mass::Float64
    function ClusterWithMass{T}(xyz::AbstractArray{<:Real,2}, atoms::Vector{<:AbstractAtomWithMass}) where T<:AbstractAtom
        if size(xyz,2) != length(atoms)
            error("Cluster has different sizes for atoms $(size(atoms)) and xyz $(size(xyz))")
        elseif size(xyz,1) != 3
            error("Cluster - xyz has wrong dimensions")
        end
        m = 0.0
        for x in atoms
            m += x.mass
        end
        new(xyz,atoms,m)
    end
    function ClusterWithMass{T}(xyz::AbstractArray{<:Real,1}, atoms::Vector{<:AbstractAtomWithMass}) where T<:AbstractAtom
        if length(xyz) != 3
            error("ClusterWithMass - xyz has wrong dimensions size=$(size(xyz))")
        end
        if length(atoms) != 1
            error("ClusterWithMass - atoms has wrong dimensions length=$(size(xyz))")
        end
        m = 0.0
        for x in atoms
            m += x.mass
        end
        new(reshape(xyz,3,1),atoms,m)
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


function Base.getindex(C::T, i::AbstractRange) where T <: AbstractClusterWithSymbols
    return T(C.xyz[:,i], C.atoms[i])
end

function Base.getindex(C::ClusterNoSymbols, i::AbstractRange)
    return ClusterNoSymbols(C.xyz[:,i])
end


function Base.print(io::IO, C::AbstractClusterWithSymbols)
    for i in 1:size(C)
        println(io, C.atoms[i].id, "    ",C.xyz[1,i], "  ", C.xyz[2,i], "  ",C.xyz[3,i])
    end
end

function distances(c1::AbstractCluster, c2::AbstractCluster)
    return pairwise(Euclidean(),c1.xyz,c2.xyz)
end

function distances(c::AbstractCluster, ur1::UnitRange, ur2::UnitRange)
    return pairwise(Euclidean(),c.xyz[:,ur1],c.xyz[:,ur2])
end


function center_coordinates(c::AbstractCluster)
    return sum(c.xyz, dims=2) / length(c)
end

function move!(c::AbstractCluster,r)
    c.xyz = c.xyz .+ r
end

function center_cluster!(c::AbstractCluster)
    move!(c, -center_coordinates(c))
end

function rotate_x!(c::AbstractCluster, theta)
    R = Matrix{Float64}(I, 3,3)
    R[2,2] = cos(theta)
    R[3,3] = R[2,2]
    R[3,2] = -sin(theta)
    R[2,3] = -R[3,2]
    c.xyz = R * c.xyz
end

function rotate_y!(c::AbstractCluster, theta)
    R = Matrix{Float64}(I, 3,3)
    R[1,1] = cos(theta)
    R[3,3] = R[1,1]
    R[3,1] = sin(theta)
    R[1,3] = -R[3,1]
    c.xyz = R * c.xyz
end


function rotate_z!(c::AbstractCluster, theta)
    R = Matrix{Float64}(I, 3,3)
    R[1,1] = cos(theta)
    R[2,2] = R[1,1]
    R[2,1] = -sin(theta)
    R[1,2] = -R[2,1]
    c.xyz = R * c.xyz
end



function print_xyz(io::IO, c::AbstractClusterWithSymbols, note=""; printheader=true)
    if printheader
        println(io, "    ",length(c))
        println(io, note)
    end
    for i in 1:length(c)
        println(io, c.atoms[i].id, "   ", c.xyz[1,i], "  ", c.xyz[2,i], "  ", c.xyz[3,i])
    end
end

end #module
