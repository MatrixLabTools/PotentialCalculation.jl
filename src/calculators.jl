module calculators

export AbstractCalculator,
       Orca,
       Calculator,
       write_input,
       read_energy,
       calculate_energy,
       clean_calculation_files,
       bsse_corrected_energy



using ..clusters



#TODO remove expicit referenzes to Cluster.xyz
function print_orca_xyz(io::IO, c::AbstractClusterWithSymbols; ghost=undef)
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


mutable struct Calculator
    method
    basis
    calculator
end


function write_input(io::IO, cal::Calculator,
                     c::AbstractClusterWithSymbols; ghost=undef)
    if typeof(cal.calculator) == Orca
        println(io, "! ", cal.method)
        println(io, "! ", cal.basis)
        println(io, "! MINIPRINT")
        if cal.calculator.ncore > 1
            println(io, "%pal nprogs $(cal.calculator.ncore)")
        end
        println(io, "%maxcore $(cal.calculator.memcore)")
        println(io, "* xyz 0 1")
        print_orca_xyz(io, c, ghost=ghost)
        println(io,"*")
    else
        error("Calculator type not recogniced")
    end
end


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


function calculate_energy(cal::Calculator, points; basename="base", ghost=undef, id="")
    clean_calculation_files(basename=basename)
    inname = "$(basename).inp"
    outname= "$(basename).out"
    cmd = pipeline(`$(cal.calculator.executable) $(inname)`, outname)
    out = Float64[]
    for p in points
        ts = time()
        open(inname,"w") do io
            write_input(io, cal, p, ghost=ghost)
        end
        run(cmd)
        push!(out, read_energy(outname) )
        te = time()
        @info "$(id) : Calculation done in $(round(te-ts, digits=1)) seconds"
    end
    return out
end


function calculate_energy(cal::Calculator, point::Cluster; basename="base", ghost=undef, id="")
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
    @info "$(id) : Calculation done in $(round(te-ts, digits=1)) seconds"
    return out
end



function bsse_corrected_energy(cal::Calculator, c1, c2; basename="base", id="")
    # expects ORCA calculator
    points = c1 .+ c2
    e = calculate_energy(cal, points, basename=basename, id=id)
    l1 = length(c1[1])
    l2 = length(c2[1]) + l1
    bsse1 = calculate_energy(cal, points, basename=basename, ghost=1:l1, id=id)
    bsse2 = calculate_energy(cal, points, basename=basename, ghost=(l1+1):l2, id=id)
    return e .- bsse1 .- bsse2
end

function bsse_corrected_energy(cal::Calculator, c1::Cluster, c2::Cluster; basename="base", id="")
    # expects ORCA calculator
    points = c1 + c2
    e = calculate_energy(cal, points, basename=basename, id=id)
    l1 = length(c1)
    l2 = length(c2) + l1
    bsse1 = calculate_energy(cal, points, basename=basename, ghost=1:l1, id=id)
    bsse2 = calculate_energy(cal, points, basename=basename, ghost=(l1+1):l2, id=id)
    return e .- bsse1 .- bsse2
end


function clean_calculation_files(;dir=".", basename="base")
    filenames=readdir(dir)
    # There is a  possible exploit here related to basename
    i = map( x -> occursin(basename, x), filenames)
    rm.(filenames[i])
end

end #module
