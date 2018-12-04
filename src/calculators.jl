module calculators

export AbstactCalculationType, Energy, Gradient, AbstractCalculator, Orca,
       Calculator, write_input, read_energy, calculate_energy, clean_calculation_files,
       bsse_corrected_energy



using ..clusters
using Distributed



function print_orca_xyz(io::IO, c::AbstractClusterWithSymbols; ghost=undef)
    if ghost == undef
        for i in 1:length(c)
            println(io, c.atoms[i].id, "   ", c.xyz[i,1], "  ", c.xyz[i,2], "  ", c.xyz[i,3])
        end
    else
        for i in 1:length(c)
            if i in ghost
                println(io, c.atoms[i].id, ":   ", c.xyz[i,1], "  ", c.xyz[i,2], "  ", c.xyz[i,3])
            else
                println(io, c.atoms[i].id, "   ", c.xyz[i,1], "  ", c.xyz[i,2], "  ", c.xyz[i,3])
            end
        end
    end
end



abstract type AbstactCalculationType end

struct Energy <: AbstactCalculationType
end

struct Gradient <: AbstactCalculationType
    analytic::Bool
    gradient() = new(true)
    gradient(analytic::Bool) = new(analytic)
end



abstract type AbstractCalculationProgram end

mutable struct Orca <: AbstractCalculationProgram
    excecutable
    ncore::UInt
    memcore::UInt
    tmp_dir
    function Orca(;excecutable="/home/teemu/mpiapps/apps/ORCA/4.0.1/orca",
                   ncore=1, maxmem=1000, tmp_dir=mktempdir())
        cd(tmp_dir)
        @info "Changed working directory to $(tmp_dir)"
        new(excecutable, ncore, maxmem, tmp_dir)
    end
end


mutable struct Calculator
    method
    basis
    calculator
    calculation_type
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


function calculate_energy(cal::Calculator, points; basename="base", ghost=undef)
    clean_calculation_files(basename=basename)
    inname = "$(basename).inp"
    outname= "$(basename).out"
    cmd = pipeline(`$(cal.calculator.excecutable) $(inname)`, outname)
    out = similar(points,Float64)
    for i in eachindex(points)
        ts = time()
        open(inname,"w") do io
            write_input(io, cal, points[i], ghost=ghost)
        end
        run(cmd)
        out[i] = read_energy(outname)
        te = time()
        @info "Pid $(myid()) : Calculation done in $(round(te-ts, digits=1)) seconds"
    end
    return out
end



function bsse_corrected_energy(cal::Calculator, c1, c2; basename="base")
    # expects ORCA calculator
    points = c1 .+ c2
    e = calculate_energy(cal, points, basename=basename)
    l1 = length(c1[1])
    l2 = length(c2[1]) + l1
    bsse1 = calculate_energy(cal, points, basename=basename, ghost=1:l1)
    bsse2 = calculate_energy(cal, points, basename=basename, ghost=(l1+1):l2)
    return e .- bsse1 .- bsse2
end


function clean_calculation_files(;dir=".", basename="base")
    #TODO dir does not do anything atm
    filenames=readdir(dir)
    i = map( x -> occursin(basename, x), filenames)
    rm.(filenames[i])
end

end #module
