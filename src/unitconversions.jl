module unitconversions

export energy_to, energy_from, change_energy

"""
    energy_to(e, unit)

Changes energy from hartrees to given unit

Suported units include:
  * cm^-1
  * eV

If unit is not recognized defaults to hartree
"""
function energy_to(e, unit)
    if unit in ["cm-1", "cm^-1"]
        return e*27.2113834*8065.54477
    elseif unit in ["ev", "eV"]
        return e*27.2113834
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

If unit is not recognized defaults to hartree
"""
function energy_from(e, unit)
    if unit in ["cm-1", "cm^-1"]
        return e/(27.2113834*8065.54477)
    elseif unit in ["ev", "eV"]
        return e/27.2113834
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
"""
function change_energy(e, from, to)
    tmp = energy_from(e, from)
    return energy_to(tmp,to)
end


end #module
