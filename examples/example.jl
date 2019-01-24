using Distributed

addprocs(8)
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using PotentialCalculation



ca = Calculator("blyp d3bj", "ma-def2-tzvp", Orca())
c1=Cluster{AtomOnlySymbol}([-6.7041359778      1.3501192944      0.0102209137;
-5.3688853815      1.2229556023      0.0440598937;
-7.2470157373      2.4374213225      0.0651311769;
 -5.0398812618      2.1435406993      0.1155201154;
 -7.2001330967      0.3718768293     -0.0703451879], AtomOnlySymbol.(["C", "O", "O", "H", "H"]) )


#c2 = Cluster{AtomOnlySymbol}([0.0 0.0 0.0; 1.2 0.0 0.0 ], AtomOnlySymbol.(["N", "N"]))
c2 = Cluster{AtomOnlySymbol}(rand(1,3), AtomOnlySymbol.(["Ar"]))


inp = InputAdaptiveSampler(ca,c1,c2, 2, 5000, startdistance=2.4, npoints=10)

data = sample_ntimes(inp, 8)
