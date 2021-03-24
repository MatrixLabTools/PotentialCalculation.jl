"""
    module Molecules

Holds information of molecules.

Main reason form this module is to implement identical atom information when fitting potential
"""
module Molecules

export AbstractMolecule,
       Molecule,
       MoleculeIdenticalInformation,
       makeidentical!


using ..IdenticalTools
using ..Atoms


abstract type AbstractMolecule end


"""
    Molecule{T<:AbstractAtom} <: AbstractMolecule

Molecule representation

# Field
- `atoms::Vector{T}` : Contains information of atoms
"""
struct Molecule{T<:AbstractAtom} <: AbstractMolecule
    atoms::Vector{T}
    function Molecule{T}(atoms::AbstractVector{T}) where T <:AbstractAtom
         new(atoms)
    end
    function Molecule{T}(atomnames) where T <:AbstractAtom
         new(convert.(T,atomnames))
    end
end

"""
    Molecule{T<:AbstractAtom} <: AbstractMolecule

Molecule representation with information of identical atoms

# Field
- `atoms::Vector{T}` : Contains information of atoms
- `identical::Identical` : identical atom information
"""
struct MoleculeIdenticalInformation{T<:AbstractAtom} <: AbstractMolecule
    atoms::Vector{T}
    identical::Identical
    function MoleculeIdenticalInformation{T}(atoms::Vector{T}) where T <:AbstractAtom
         new(atoms,Identical(length(atoms)))
    end
    function MoleculeIdenticalInformation{T}(atomnames) where T <:AbstractAtom
         new(convert.(T,atomnames),Identical(length(atomnames)))
    end
end


"""
    makeidentical!(Mol::MoleculeIdenticalInformation,x)

Adds information for identical atoms to molecule.

# Arguments
- `Mol::MoleculeIdenticalInformation` : moleculet to which identical information is added
- `x` : collection of indices of identical atoms

## Throws
Error if `x` is out of bounds of `Mol`
"""
function makeidentical!(Mol::MoleculeIdenticalInformation,x)
    if maximum(x) <= length(Mol.atoms) && minimum(x) >= 1
        push!(Mol.identical,x)
    else
        error("Molecule has only $(length(Mol.atoms)) atoms")
    end
end


"""
    areidentical(mol::MoleculeIdenticalInformation,x)

Checks if atoms with indices given in collection `x` are identical
"""
IdenticalTools.areidentical(mol::MoleculeIdenticalInformation,x) = IdenticalTools.areidentical(mol.identical,x)

Base.convert(t::Type{<:AbstractMolecule}, m::AbstractMolecule) = t(m.atoms)

Base.length(mol::AbstractMolecule) = length(mol.atoms)

Base.show(io::IO, mol::AbstractMolecule) = print(io,"Molecule of $(length(mol)) atoms")

end #module
