module molecules

export AbstractMolecule, Molecule, MoleculeIdenticalInformation


using ..identical, ..atoms


abstract type AbstractMolecule end

struct Molecule{T<:AbstractAtom} <: AbstractMolecule
    atoms::Vector{T}
    function Molecule{T}(atoms::Vector{T}) where T <:AbstractAtom
         new(atoms)
    end
    function Molecule{T}(atomnames) where T <:AbstractAtom
         new(convert.(T,atomnames))
    end
end

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

function Base.push!(Mol::MoleculeIdenticalInformation,x)
    if maximum(x) <= length(Mol.atoms)
        push!(Mol.identical,x)
    else
        error("Molecule has only $(length(Mol.atoms)) atoms")
    end
 end

areidentical(mol::MoleculeIdenticalInformation,x) = areidentical(mol.identical,x)

Base.convert(t::Type{<:AbstractMolecule}, m::AbstractMolecule) = t(m.atoms)

Base.length(mol::AbstractMolecule) = length(mol.atoms)

Base.show(io::IO, mol::AbstractMolecule) = print(io,"Molecule of $(length(mol)) atoms")

end #module
