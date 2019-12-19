module atoms


export m_electron, m_u, m_au, proton_mass, masses

export AbstractAtom, AbstractAtomWithMass,
       AtomOnlySymbol, AtomWithMass


using ..identical


const m_electron = 9.10938291E-31
const m_u = 1.660539040E-27
const m_au = m_u/m_electron
const proton_mass = 1836.0
const masses =  Dict("H"  =>   1.00782503223*m_au,
                     "O"  =>  15.99491461956*m_au,
                     "C"  =>  12.0*m_au,
                     "Ar" =>  39.9623831225*m_au,
                     "N"  =>  14.0030740048*m_au,
                     "D"  =>   2.01410177812*m_au,
                     "Ne" =>  19.9924401762*m_au)



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
    AtomWithMass(id::String) = new(id,masses[id])
    AtomWithMass(id::String, mass::Real) = new(id,mass)
end


Base.convert(t::Type{<:AbstractAtom}, x::String) = t(x)
Base.convert(t::Type{<:AbstractAtom}, x::AbstractAtom) = t(x.id)

end #module
