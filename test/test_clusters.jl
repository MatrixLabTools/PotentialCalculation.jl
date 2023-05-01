using Test
using PotentialCalculation.IdenticalTools
using PotentialCalculation.Atoms
using PotentialCalculation.Clusters
using AtomsBase



@testset "Clusters" begin
    formic_acid=Cluster(
     [-6.7041359778      1.3501192944      0.0102209137
    -5.3688853815      1.2229556023      0.0440598937
    -7.2470157373      2.4374213225      0.0651311769
     -5.0398812618      2.1435406993      0.1155201154
     -7.2001330967      0.3718768293     -0.0703451879]',
      AtomOnlySymbol.(["C", "O", "O", "H", "H"]) )

     Ar = Cluster( Atom(:Ar, rand(3)u"Å"  )  )

     t1 = formic_acid[5]
     t2 = formic_acid[end]

     @test position(t1) == position(t2)


     @test length(formic_acid(1:2)) == 2

     @test length(formic_acid + Ar) == 6

     distances(formic_acid, Ar)
     distances(formic_acid, 1:2, 3:4)

     distances(formic_acid, 1:2, formic_acid, 2:4)

     print(devnull, formic_acid)

    
    cc = center_coordinates(formic_acid)
    center_cluster!(formic_acid)
    @test all(abs.(center_coordinates(formic_acid)) .< 1E-12u"Å" .* ones(size(cc)))

     l = distances(formic_acid, 2,4)
     a = cluster_angle(formic_acid, 1,2,4)
     d = dihedral_angle(formic_acid, 3,1,2,4)

     cluster_angle(formic_acid, 1,2, formic_acid, 3)

     @test l != distances(formic_acid, 1,2)
     @test a != cluster_angle(formic_acid, 1,2,3)
     @test d != dihedral_angle(formic_acid, 1,2,3,5)

     for f in [rotate_x!, rotate_y!, rotate_z!]
         f(formic_acid, rand())
         tl = distances(formic_acid, 2,4)
         ta = cluster_angle(formic_acid, 1,2,4)
         td = dihedral_angle(formic_acid, 3,1,2,4)
         @test l ≈ tl && a ≈ ta && d ≈ td
     end

     @test_throws DimensionMismatch Cluster(rand(2), [AtomOnlySymbol("H")])
     @test_throws DimensionMismatch Cluster(rand(2), AtomOnlySymbol("H"))
     @test_throws DimensionMismatch Cluster(rand(3),AtomOnlySymbol.(["H", "O"]))
     @test_throws DimensionMismatch Cluster(rand(3,3),AtomOnlySymbol.(["H", "O"]))
     @test_throws DimensionMismatch Cluster(rand(4,2),AtomOnlySymbol.(["H", "O"]))

     # AtomsBase 
     fa = FlexibleSystem(formic_acid; e=1)
     pf = Cluster(collect(formic_acid))
     @test haskey(fa, :e)
     cfa = Cluster(fa)
     @test all( atomic_symbol(cfa) .== atomic_symbol(formic_acid) )
     @test all( position(cfa) .== position(formic_acid) )
end
