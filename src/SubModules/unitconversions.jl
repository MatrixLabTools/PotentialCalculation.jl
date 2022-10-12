module UnitConversions

using Unitful
using UnitfulAtomic
using UnitfulEquivalences



export change_energy_unit
export energy_from
export energy_to



"""
    energy_to(e, unit)

Changes energy from hartrees to given unit

Suported units include:
  * cm^-1/cm⁻¹/cm-1
  * eV
  * hartree
  * kcal/mol
  * kJ/mol
  * K

If unit is not recognized defaults to hartree.

Source is [Wikipedia](https://en.wikipedia.org/wiki/Hartree)
"""
function energy_to(e, unit::AbstractString)
    if unit in ["cm-1", "cm^-1", "cm⁻¹"]
        return e*219474.6313702
    elseif unit in ["ev", "eV"]
        return e*27.21138602
    elseif unit in ["kcal/mol"]
        return e*627.509474
    elseif unit in ["kJ/mol", "kj/mol"]
        return e*2625.499639
    elseif unit in ["K"]
        return e*315775.13
    else
        return e
    end
end


function energy_to(e, u::Unitful.EnergyUnits)
    return (ustrip ∘ uconvert)(u, e*u"hartree")
end

function energy_to(e, u)
    return (ustrip ∘ uconvert)(u, e*u"hartree", Spectral())
end


"""
    energy_from(e, unit)

Changes given energy unit to hartree

Suported units include:
  * cm^-1/cm⁻¹/cm-1
  * eV
  * hartree
  * kcal/mol
  * kJ/mol
  * K

If unit is not recognized defaults to hartree.

Source is [Wikipedia](https://en.wikipedia.org/wiki/Hartree)
"""
function energy_from(e, unit)
    if unit in ["cm-1", "cm^-1", "cm⁻¹"]
        return e/219474.6313702
    elseif unit in ["ev", "eV"]
        return e/27.21138602
    elseif unit in ["kcal/mol"]
        return e/627.509474
    elseif unit in ["kJ/mol", "kj/mol"]
        return e/2625.499639
    elseif unit in ["K"]
        return e/315775.13
    else
        return e
    end
end

function energy_from(e::Unitful.Energy)
    return austrip(e)
end

function energy_from(e::Unitful.Wavenumber)
    return (ustrip ∘ uconvert)(u"hartree", e, Spectral())
end

function energy_from(e::Unitful.Frequency)
    return (ustrip ∘ uconvert)(u"hartree", e, Spectral())
end


"""
    change_energy(e, from, to)

Changes energy unit

Suported units include:
  * cm^-1/cm⁻¹/cm-1
  * eV
  * hartree
  * kcal/mol
  * kJ/mol
  * K

If unit is not recognized defaults to hartree.

Source is [Wikipedia](https://en.wikipedia.org/wiki/Hartree)
"""
function change_energy_unit(e, from, to)
    tmp = energy_from(e, from)
    return energy_to(tmp,to)
end


end #module
