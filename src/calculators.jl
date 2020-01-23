module calculators

export AbstractCalculator,
       AbstractCalculationProgram,
       bsse_corrected_energy,
       calculate_energy,
       Calculator,
       clean_calculation_files,
       getBSSEsteps,
       Orca,
       read_energy,
       write_input






using ..clusters
using Distributed


"""
print_orca_xyz(io::IO, c::AbstractClusterWithSymbols; ghost=undef)

Prints ORCA input style xyz input. Idexes in `ghost` are marked as ghost atoms.
"""
function print_orca_xyz(io::IO, c::AbstractClusterWithSymbols; ghost=undef)
    #TODO remove expicit referenzes to Cluster.xyz
    if ghost == undef
        for i in 1:length(c)
            println(io, c.atoms[i].id, "   ", c.xyz[1,i], "  ", c.xyz[2,i], "  ", c.xyz[3,i])
        end
    else
        for i in 1:length(c)
            if i in ghost
                println(io, c.atoms[i].id, ":   ", c.xyz[1,i], "  ", c.xyz[2,i], "  ", c.xyz[3,i])
            else
                println(io, c.atoms[i].id, "   ", c.xyz[1,i], "  ", c.xyz[2,i], "  ", c.xyz[3,i])
            end
        end
    end
end





abstract type AbstractCalculationProgram end

"""
    mutable struct Orca <: AbstractCalculationProgram

Holds information of how to use ORCA program.

Note that if you are using more than 1 core you need to specify exact
path to orca executable.

# Fields
- `executable` : path to "orca" executable - default = `"orca"`
- `ncore` : number of cores used by orca - default = `1`
- `maxmem` : maximum memory per core for orca - default = `1000`
- `tmp_dir` : directory where calculations are done - default = `mktempdir()`
"""
mutable struct Orca <: AbstractCalculationProgram
    "path for orca excecutable"
    executable
    "number of cores in calculation"
    ncore::UInt
    "maximum memory per core"
    memcore::UInt
    "directory where calculation is done"
    tmp_dir
    function Orca(;executable="orca",
                   ncore=1, maxmem=1000, tmp_dir=mktempdir())
        new(executable, ncore, maxmem, tmp_dir)
    end
end

"""
    getBSSEsteps(cal::AbstractCalculationProgram)

Used to set up progressbar. Returns 1, if calculator program has a method to calculate
counter poise correction with a single command and 3, if program needs to be called
multiple times to make counter poise correction.
"""
getBSSEsteps(cal::AbstractCalculationProgram) = 1
getBSSEsteps(cal::Orca) = 3


"""
    mutable struct Calculator

Struct used to hold together different bits of information needed in calculations.

# Fields
- `method` : calculation method informaton
- `basis` : basis information
- `calculator` : calculation program information

# Examples
```julia
julia> Calculator{Orca}("method", "basis")
Calculator{Orca}("method", "basis", Orca("orca", 0x0000000000000001, 0x00000000000003e8, "/tmp/tmpxEoJYW"))

julia> Calculator{Orca}("method", "basis", Orca())
Calculator{Orca}("method", "basis", Orca("orca", 0x0000000000000001, 0x00000000000003e8, "/tmp/tmpVg959k"))
```
"""
mutable struct Calculator{T<:AbstractCalculationProgram}
    method
    basis
    calculator::T
    function Calculator{T}(method::AbstractString, basis::AbstractString) where T <: AbstractCalculationProgram
        new(method, basis, T())
    end
    function Calculator{T}(method::AbstractString, basis::AbstractString, calculator::T) where T <: AbstractCalculationProgram
        new(method, basis, calculator)
    end
end



"""
    write_input(io::IO, cal::Calculator,
                     c::AbstractClusterWithSymbols; ghost=undef)

Writes input files for calculation program. Only ORCA input is supported now.

# Arguments
- `io::IO` : stream where writing is done
- `cal::Calculator` : calculation information
- `c::AbstractClusterWithSymbols` : molecule cluster that is calculated
- `ghost=undef` : collection for atom idexes that are considered ghost atoms

"""
function write_input(io::IO, cal::Calculator{Orca},
                     c::AbstractClusterWithSymbols; ghost=undef)
    if typeof(cal.calculator) == Orca
        println(io, "! ", cal.method)
        println(io, "! ", cal.basis)
        println(io, "! MINIPRINT")
        if cal.calculator.ncore > 1
            println(io, "%pal nprocs $(cal.calculator.ncore) end")
        end
        println(io, "%maxcore $(cal.calculator.memcore)")
        println(io, "* xyz 0 1")
        print_orca_xyz(io, c, ghost=ghost)
        println(io,"*")
    else
        error("Calculator type not recogniced")
    end
end



"""
    read_energy(fname)

Read ORCA output for final energy. Throws an error if calculation failed.
"""
function read_energy(fname)
    lines = ""
    open(fname, "r") do f
        lines = readlines(f)
    end
    if ! occursin("****ORCA TERMINATED NORMALLY****", lines[end-1])
        error("Orca failed : $(lines[end-1])")
    end
    for l in reverse(lines)
        if occursin("FINAL SINGLE POINT ENERGY", l)
            return parse(Float64, split(l)[end])
        end
    end
end



"""
    calculate_energy(cal::Calculator, points; basename="base", ghost=undef, id="", pchannel=pchannel)

Calculates energy for given clusters.

# Arguments
- `cal::Calculator` : calculation information
- `points` : collection of [`Cluster`](@ref) which energy is calculated

# Keywords
- `basename="base"` : base name for input/ouput files
- `ghost=undef` : indices for atoms considered as ghost
- `id=""` : additional information for calculator - needed for multiprocessing in same folder
- `pchannel=undef`  : (Remote)Channel where progess information is added
"""
function calculate_energy(cal::Calculator{Orca}, points; basename="base", ghost=undef, id="", pchannel=undef)
    clean_calculation_files(basename=basename)
    inname = "$(basename).inp"
    outname= "$(basename).out"
    cmd = pipeline(`$(cal.calculator.executable) $(inname)`, outname)
    function do_calculation(point)
        ts = time()
        open(inname,"w") do io
            write_input(io, cal, point, ghost=ghost)
        end
        run(cmd)
        out = read_energy(outname)
        te = time()
        @debug "Calculation done in $(round(te-ts, digits=1)) seconds"
        pchannel != undef && put!(pchannel,true)
        return out
    end
    return map(x -> do_calculation(x), points)
end


"""
    calculate_energy(cal::Calculator, point::Cluster; basename="base", ghost=undef, id="", pchannel=pchannel)

Calculates energy for given cluster.

# Arguments
- `cal::Calculator` : calculation information
- `point::Cluster` : cluster which energy is calculated

# Keywords
- `basename="base"` : base name for input/ouput files
- `ghost=undef` : indices for atoms considered as ghost
- `id=""` : additional information for calculator - needed for multiprocessing in same folder
- `pchannel=undef`  : (Remote)Channel where progess information is added
"""
function calculate_energy(cal::Calculator{Orca}, point::Cluster; basename="base", ghost=undef, id="", pchannel=undef)
    clean_calculation_files(basename=basename)
    inname = "$(basename).inp"
    outname= "$(basename).out"
    cmd = pipeline(`$(cal.calculator.executable) $(inname)`, outname)

    ts = time()
    open(inname,"w") do io
        write_input(io, cal, point, ghost=ghost)
    end
    run(cmd)
    out = read_energy(outname)
    te = time()
    @debug "Calculation done in $(round(te-ts, digits=1)) seconds"
    pchannel != undef && put!(pchannel,true)
    return out
end


"""
    bsse_corrected_energy(cal::Calculator, c1, c2; basename="base", id="", pchannel=pchannel)

Calculates energy of combined clusters taking into account basis set superposition error

# Arguments
- `cal::Calculator` : calcualation information
- `c1` : collection of clusters
- `c2` : collection of clusters

# Keywords
- `basename="base"` : base name for input/ouput files
- `id=""` : additional information for calculator - needed for multiprocessing in same folder
- `pchannel=undef`  : (Remote)Channel where progess information is added
"""
function bsse_corrected_energy(cal::Calculator{Orca}, c1, c2; basename="base", id="",
                                pchannel=undef)
    # expects ORCA calculator
    points = c1 .+ c2
    e = calculate_energy(cal, points, basename=basename, id=id, pchannel=pchannel)
    l1 = length(c1[1])
    l2 = length(c2[1]) + l1
    bsse1 = calculate_energy(cal, points, basename=basename, ghost=1:l1, id=id, pchannel=pchannel)
    bsse2 = calculate_energy(cal, points, basename=basename, ghost=(l1+1):l2, id=id, pchannel=pchannel)
    return e .- bsse1 .- bsse2
end

"""
    bsse_corrected_energy(cal::Calculator, c1::Cluster, c2::Cluster; basename="base", id="", pchannel=pchannel)

Calculates energy of combined cluster taking into account basis set superposition error

# Arguments
- `cal::Calculator` : calcualation information
- `c1::Cluster` : cluster
- `c2::Cluster` : cluster

# Keywords
- `basename="base"` : base name for input/ouput files
- `id=""` : additional information for calculator - needed for multiprocessing in same folder
- `pchannel=undef`  : (Remote)Channel where progess information is added
"""
function bsse_corrected_energy(cal::Calculator{Orca}, c1::Cluster, c2::Cluster;
                               basename="base", id="", pchannel=undef)
    # expects ORCA calculator
    points = c1 + c2
    e = calculate_energy(cal, points, basename=basename, id=id, pchannel=pchannel)
    l1 = length(c1)
    l2 = length(c2) + l1
    bsse1 = calculate_energy(cal, points, basename=basename, ghost=1:l1, id=id, pchannel=pchannel)
    bsse2 = calculate_energy(cal, points, basename=basename, ghost=(l1+1):l2, id=id, pchannel=pchannel)
    return e .- bsse1 .- bsse2
end

"""
    clean_calculation_files(;dir=".", basename="base")

Deletes files that calculations produced

# Arguments
- `dir="."` : directory to be cleaned
- `basename="base"` : base name used in calculations
"""
function clean_calculation_files(;dir=".", basename="base")
    filenames=readdir(dir)
    # There is a  possible exploit here related to basename
    i = map( x -> occursin(basename*".", x) || occursin(basename*"_", x), filenames)
    rm.(filenames[i])
end

end #module
