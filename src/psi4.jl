module psi4

using PyCall
using ..calculators
using ..clusters

export Psi4

global gpsi4init=false
global gPsi4 = undef


function initpsi4(;memory="500 MiB", quiet=true)
    global gPsi4 = pyimport("psi4")
    global gpsi4init = true
    gPsi4.set_memory(memory)
    quiet && gPsi4.core.be_quiet()
end

function setpsi4memory(memory)
    !
end

mutable struct Psi4 <: AbstractCalculationProgram
    memory
    tmpdir
    function Psi4(;memory="500MiB", tmpdir=undef)
        initpsi4(memory=memory)
        @debug gPsi4
        new(memory, tmpdir)
    end
end


function calculators.calculate_energy(cal::Calculator{Psi4}, point::Cluster; basename="base", ghost=undef, id="", pchannel=undef)
    ! gpsi4init && initpsi4(memory=cal.calculator.memory)
    s=sprint( (io, x) -> print_xyz(io,x, printheader=false), point)
    c = gPsi4.geometry(s)
    out = gPsi4.energy(cal.method*"/"*cal.basis, molecule=c)
    pchannel != undef && put!(pchannel,true)
    return out
end


function calculators.calculate_energy(cal::Calculator{Psi4}, points; basename="base", ghost=undef, id="", pchannel=undef)
    return map( x -> calculate_energy(cal, x, basename=basename, ghost=ghost, id=id, pchannel=pchannel), points )
end


function calculators.bsse_corrected_energy(cal::Calculator{Psi4}, c1::Cluster, c2::Cluster;
                               basename="base", id="", pchannel=undef)
    ! gpsi4init && initpsi4(memory=cal.calculator.memory)
    s1=sprint( (io, x) -> print_xyz(io,x, printheader=false), c1)
    s2=sprint( (io, x) -> print_xyz(io,x, printheader=false), c2)
    c = gPsi4.geometry(s1*"\n--\n"*s2)
    out = gPsi4.energy(cal.method*"/"*cal.basis, molecule=c, bsse_type="cp")
    pchannel != undef && put!(pchannel,true)
    return out
end

function calculators.bsse_corrected_energy(cal::Calculator{Psi4}, c1, c2; basename="base", id="",
                                pchannel=undef)
    return map( (x,y) -> bsse_corrected_energy(cal, x, y, basename=basename, id=id, pchannel=pchannel ), c1, c2  )
end


end  # module psi4
