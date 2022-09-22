using Test
using PotentialCalculation

@testset "atoms" begin

    a = AtomOnlySymbol("H")

    am = convert(AtomWithMass, a)
    aa = convert(AtomOnlySymbol, "H")

    @test a.id == am.id

    x = AtomWithMass("x", 3.2u"u")

    O = AtomWithMass("O")

    @test am.mass == austrip( masses["H"] )
    @test O.mass == austrip( masses["O"] )
end
