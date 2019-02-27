using Test
using PotentialCalculation.atoms

@testset "atoms" begin

    a = AtomOnlySymbol("H")

    am = convert(AtomWithMass, a)

    @test a.id == am.id

    x = AtomWithMass("x", 3.2)

    O = AtomWithMass("O")

    @test am.mass == masses["H"]
    @test O.mass == masses["O"]
end
