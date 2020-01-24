module psi4

using PyCall
using ..calculators
using ..clusters

export Psi4

"""
    gpsi4init=false

Global variable to see if Psi4 was initiated for current process
"""
global gpsi4init=false

"""
    gPsi4 = undef

Hold Psi4 object given by PyCall. Used to do calculation in current process
"""
global gPsi4 = undef

"""
    initpsi4(;memory="500 MiB", quiet=true, nthreads=1)

Used to intialize Psi4 environment
"""
function initpsi4(;memory="500 MiB", quiet=true, nthreads=1)
    global gPsi4 = pyimport("psi4")
    global gpsi4init = true
    quiet && gPsi4.core.be_quiet()
    gPsi4.set_memory(memory)
    nthreads > 1 && gPsi4.set_num_threads(nthreads)
end

"""
    mutable struct Psi4 <: AbstractCalculationProgram

Holds information that calculations are to be done with Psi4
"""
mutable struct Psi4 <: AbstractCalculationProgram
    memory::String
    nthreads::UInt
    function Psi4(;memory="500MiB", nthreads=1)
        initpsi4(memory=memory, nthreads=nthreads)
        @debug gPsi4
        new(memory, nthreads)
    end
end


function calculators.calculate_energy(cal::Calculator{Psi4}, point::Cluster;
                                      basename="base", ghost=undef, id="", pchannel=undef)
    ! gpsi4init && initpsi4(memory=cal.calculator.memory)
    s=sprint( (io, x) -> print_xyz(io,x, printheader=false), point)
    c = gPsi4.geometry(s)
    out = gPsi4.energy(cal.method*"/"*cal.basis, molecule=c)
    pchannel != undef && put!(pchannel,true)
    return out
end


function calculators.calculate_energy(cal::Calculator{Psi4}, points;
                                      basename="base", ghost=undef, id="", pchannel=undef)
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

function calculators.bsse_corrected_energy(cal::Calculator{Psi4}, c1, c2;
                                           basename="base", id="", pchannel=undef)
    return map( (x,y) -> bsse_corrected_energy(cal, x, y, basename=basename, id=id, pchannel=pchannel ), c1, c2  )
end


end  # module psi4
