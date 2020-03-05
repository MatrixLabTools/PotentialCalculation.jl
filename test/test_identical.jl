using Test
using PotentialCalculation.IdenticalTools

@testset "IdenticalTools" begin
    ii = Identical()
    push!(ii, (2,3))
    @test ! areidentical(ii, (1,2))
    @test areidentical(ii, (3,2))

    jj = Identical(10)
    @test length(jj.identical) == 10

    push!(jj, (3,6,8))
    @test length(jj.identical) == 8

    @test ! areidentical(jj, (2,3))
    @test areidentical(jj, (3,6))
    @test areidentical(jj, (3,8))
    @test areidentical(jj, (3,6,8))
end
