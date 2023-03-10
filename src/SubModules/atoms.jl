module Atoms

using AtomsBase
using Unitful
using UnitfulAtomic

export proton_mass, masses

export AbstractAtom
export AbstractAtomWithMass
export AtomOnlySymbol
export AtomWithMass


using ..IdenticalTools


const proton_mass = 1836.0u"me_au"
const masses =  Dict("H"  =>   1.00782503223u"u",
                     "O"  =>  15.99491461956u"u",
                     "C"  =>  12.0u"u",
                     "Ar" =>  39.9623831225u"u",
                     "N"  =>  14.0030740048u"u",
                     "D"  =>   2.01410177812u"u",
                     "Ne" =>  19.9924401762u"u")



abstract type AbstractAtom end
abstract type AbstractAtomWithMass <: AbstractAtom end

"""
    AtomOnlySymbol <: AbstractAtom

Atom that only has symbol

# Fields
- `id::String` : hold description of atom
"""
struct AtomOnlySymbol <: AbstractAtom
    id::String
end

function AtomOnlySymbol(a::AtomsBase.Atom)
    return AtomOnlySymbol(String(atomic_symbol(a)))
end

function AtomOnlySymbol(s::Symbol)
    return AtomOnlySymbol(String(s))
end

"""
    AtomWithMass <: AbstractAtomWithMass

Atom that has symbol and mass

# Fields
- `id::String`    : hold description of atom
- `mass::Float64` : mass of atom
"""
struct AtomWithMass <: AbstractAtomWithMass
    id::String
    mass::Float64
    AtomWithMass(id::AbstractString) = new(id, austrip(masses[id]) )
    function AtomWithMass(id::AbstractString, mass)
        @assert dimension(mass) == dimension(u"kg") || dimension(mass) == NoDims
        new(id,austrip(mass))
    end
end


Base.convert(t::Type{<:AbstractAtom}, x::AbstractString) = t(x)
Base.convert(t::Type{<:AbstractAtom}, x::AbstractAtom) = t(x.id)

end #module

