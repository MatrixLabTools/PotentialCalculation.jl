using Test
using PotentialCalculation.Atoms
using PotentialCalculation.Clusters
using PotentialCalculation.Sample



@testset "Sample" begin
    formic_acid=Cluster(
     [-6.7041359778      1.3501192944      0.0102209137
    -5.3688853815      1.2229556023      0.0440598937
    -7.2470157373      2.4374213225      0.0651311769
     -5.0398812618      2.1435406993      0.1155201154
     -7.2001330967      0.3718768293     -0.0703451879]',
      AtomOnlySymbol.(["C", "O", "O", "H", "H"])
    )
    Ar = Cluster( Atom(:Ar, rand(3)u"Å"  )  )

    tmp = deepcopy(formic_acid)

    Clusters.rotate_randomly!(tmp)

    @test tmp.xyz != formic_acid.xyz

    t = line_sampler(formic_acid, tmp)

    @test 2.0u"Å" <= minimum(distances(t[1][1], t[2][1]))

    if Sys.which("orca") !== nothing && Sys.which("orca_scf") !== nothing
        @info "Orca binary found. Testing ORCA."
        ca = Calculator("blyp d3bj TIGHTSCF", "def2-svp", Orca())
        a = PotentialCalculation.random_sampler(ca, [formic_acid], [Ar], 2; npoints=3, move_step=0.2u"Å")
        @test haskey(a, "Points")
        @test haskey(a, "cluster1")
        @test haskey(a, "cluster2")
        @test haskey(a, "Energy")
        @test size(a["Points"]) == size(a["Energy"]) == (3,2)
        @test length(a["cluster1"]) == length(formic_acid)
        @test length(a["cluster2"]) == length(Ar)
    end
end
