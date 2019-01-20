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


mutable struct ClusterNoSymbols <: AbstractCluster
    xyz::Array{Float64,2}
    function ClusterNoSymbols(xyz::Array{<:Real,2})
        if size(xyz,2) != 3
            error("Cluster - xyz has wrong dimensions")
        end
        new(xyz)
    end
end


mutable struct Cluster{T<:AbstractAtom} <: AbstractClusterWithSymbols
    xyz::Array{Float64,2}
    atoms::Vector{T}
    function Cluster{T}(xyz::Array{<:Real,2}, atoms:: Vector{<:T}) where T<:AbstractAtom
        if size(xyz,1) != length(atoms)
            error("Cluster has different sizes for atoms $(size(atoms)) and xyz $(size(xyz))")
        elseif size(xyz,2) != 3
            error("Cluster - xyz has wrong dimensions")
        end
        new(xyz,atoms)
    end
end


mutable struct ClusterWithMass{T<:AbstractAtom} <: AbstractClusterWithSymbols
    xyz::Array{Float64,2}
    atoms::Vector{T}
    mass::Float64
    function ClusterWithMass{T}(xyz::Array{<:Real,2}, atoms::Vector{<:AbstractAtomWithMass}) where T<:AbstractAtom
        if size(xyz,1) != length(atoms)
            error("Cluster has different sizes for atoms $(size(atoms)) and xyz $(size(xyz))")
        elseif size(xyz,2) != 3
            error("Cluster - xyz has wrong dimensions")
        end
        m = 0.0
        for x in atoms
            m += x.mass
        end
        new(xyz,atoms,m)
    end
end


Base.show(io::IO, C::AbstractCluster) = print(io,"$(typeof(C)) of ", size(C), " atoms")

Base.size(a::AbstractCluster) = size(a.xyz,1)

Base.length(a::AbstractCluster) = size(a)

Base.lastindex(a::AbstractCluster) = length(a)


function (+)(c1::T, c2::T) where T <: AbstractClusterWithSymbols
    return T(vcat(c1.xyz,c2.xyz),vcat(c1.atoms,c2.atoms))
end

function (+)(c1::ClusterNoSymbols, c2::ClusterNoSymbols)
    return ClusterNoSymbols(vcat(c1.xyz,c2.xyz))
end


function Base.getindex(C::T, i::Int) where T <: AbstractClusterWithSymbols
    return T(reshape(C.xyz[i,:],1,3), [C.atoms[i]])
end

function Base.getindex(C::ClusterNoSymbols, i::Int)
    return ClusterNoSymbols(reshape(C.xyz[i,:],1,3))
end


function Base.getindex(C::T, i::AbstractRange) where T <: AbstractClusterWithSymbols
    return T(C.xyz[i,:], C.atoms[i])
end

function Base.getindex(C::ClusterNoSymbols, i::AbstractRange)
    return ClusterNoSymbols(C.xyz[i,:])
end


function Base.print(io::IO, C::AbstractClusterWithSymbols)
    for i in 1:size(C)
        println(io, C.atoms[i].id, "    ",C.xyz[i,1], "  ", C.xyz[i,2], "  ",C.xyz[i,3])
    end
end

function distances(c1::AbstractCluster, c2::AbstractCluster)
    return pairwise(Euclidean(),transpose(c1.xyz),transpose(c2.xyz))
end

function distances(c::AbstractCluster, ur1::UnitRange, ur2::UnitRange)
    return pairwise(Euclidean(),transpose(c.xyz[ur1,:]),transpose(c.xyz[ur2,:]))
end


function center_coordinates(c::AbstractCluster)
    return sum(c.xyz, dims=1) / length(c)
end

function move!(c::AbstractCluster,r)
    c.xyz = c.xyz .+ reshape(r, (1,3))
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
    c.xyz = c.xyz * R
end

function rotate_y!(c::AbstractCluster, theta)
    R = Matrix{Float64}(I, 3,3)
    R[1,1] = cos(theta)
    R[3,3] = R[1,1]
    R[3,1] = sin(theta)
    R[1,3] = -R[3,1]
    c.xyz = c.xyz * R
end


function rotate_z!(c::AbstractCluster, theta)
    R = Matrix{Float64}(I, 3,3)
    R[1,1] = cos(theta)
    R[2,2] = R[1,1]
    R[2,1] = -sin(theta)
    R[1,2] = -R[2,1]
    c.xyz = c.xyz * R
end



function print_xyz(io::IO, c::AbstractClusterWithSymbols, note=""; printheader=true)
    if printheader
        println(io, "    ",length(c))
        println(io, note)
    end
    for i in 1:length(c)
        println(io, c.atoms[i].id, "   ", c.xyz[i,1], "  ", c.xyz[i,2], "  ", c.xyz[i,3])
    end
end

end #module
