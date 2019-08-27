using Test
using PotentialCalculation.identical
using PotentialCalculation.atoms
using PotentialCalculation.molecules



@testset "molecules" begin
    m = Molecule{AtomOnlySymbol}(["A", "B", "Ca"])
    mm = Molecule{AtomOnlySymbol}(AtomOnlySymbol.(["A", "B", "Ca"]))
    mi = MoleculeIdenticalInformation{AtomOnlySymbol}(["A", "B", "Ca"])
    mmi = MoleculeIdenticalInformation{AtomOnlySymbol}(AtomOnlySymbol.(["A", "B", "Ca"]))

    mt = convert(typeof(mi), m)
    @test mt.atoms == m.atoms

    show(devnull, mm)

    @test length(m) == 3
    @test length(m) == length(mi)

    makeidentical!(mi, (1,2))

    @test ! areidentical(mi, (3,2))
    @test areidentical(mi, (1,2))

    @test try
        makeidentical!(mi, (3,4))
    catch
        true
    end

    @test try
        makeidentical!(mi, (0,2))
    catch
        true
    end

end
