using Test
using PotentialCalculation.identical
using PotentialCalculation.atoms
using PotentialCalculation.molecules



@testset "molecules" begin
    m = Molecule{AtomOnlySymbol}(["A", "B", "Ca"])
    mi = MoleculeIdenticalInformation{AtomOnlySymbol}(["A", "B", "Ca"])

    @test length(m) == 3
    @test length(m) == length(mi)

    push!(mi, (1,2))

    @test ! areidentical(mi, (3,2))
    @test areidentical(mi, (1,2))

    @test try
        push!(mi, (3,4))
    catch
        true
    end

    @test try
        push!(mi, (0,2))
    catch
        true
    end

end
