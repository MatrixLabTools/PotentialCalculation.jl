using PotentialCalculation

ca = Calculator("RI-MP2 RIJK", "def2-svp def2-svp/C def2/JK", Orca(), Energy())
c1 = Cluster{AtomOnlySymbol}([0.0 0.0 0.0; 1.2 0.0 0.0 ], AtomOnlySymbol.(["N", "N"]))
c2 = Cluster{AtomOnlySymbol}(rand(1,3), AtomOnlySymbol.(["Ar"]))

a = adaptive_line_sampler(ca, c1, c2, 5000)
