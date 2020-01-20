using Test
using PotentialCalculation.atoms
using PotentialCalculation.clusters
using PotentialCalculation.sample



@testset "sample" begin
    formic_acid=Cluster(
     [-6.7041359778      1.3501192944      0.0102209137
    -5.3688853815      1.2229556023      0.0440598937
    -7.2470157373      2.4374213225      0.0651311769
     -5.0398812618      2.1435406993      0.1155201154
     -7.2001330967      0.3718768293     -0.0703451879]',
      AtomOnlySymbol.(["C", "O", "O", "H", "H"]) )

      tmp = deepcopy(formic_acid)

      sample.random_rotation!(tmp)

      @test tmp.xyz != formic_acid.xyz

      t = line_sampler(formic_acid, tmp)

      @test 2.0 <= minimum(distances(t[1][1], t[2][1]))
end
