module unitconversions

export change_energy_unit,
       energy_from,
       energy_to



"""
    energy_to(e, unit)

Changes energy from hartrees to given unit

Suported units include:
  * cm^-1
  * eV
  * hartree
  * kcal/mol
  * kJ/mol
  * K

If unit is not recognized defaults to hartree.

Source is [Wikipedia](https://en.wikipedia.org/wiki/Hartree)
"""
function energy_to(e, unit)
    if unit in ["cm-1", "cm^-1"]
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

"""
    energy_from(e, unit)

Changes given energy unit to hartree

Suported units include:
  * cm^-1
  * eV
  * hartree
  * kcal/mol
  * kJ/mol
  * K

If unit is not recognized defaults to hartree.

Source is [Wikipedia](https://en.wikipedia.org/wiki/Hartree)
"""
function energy_from(e, unit)
    if unit in ["cm-1", "cm^-1"]
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


"""
    change_energy(e, from, to)

Changes energy unit

Suported units include:
  * cm^-1
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
