using Test
using PotentialCalculation

@testset "restarttools" begin
# Tests basic calculations

fname = tempname()
rname = tempname()
sname = tempname()
xyzname = tempname()

ca = Calculator("blyp d3bj", "ma-def2-tzvp", Orca())

formic_acid=Cluster{AtomOnlySymbol}(
 [-6.7041359778      1.3501192944      0.0102209137
-5.3688853815      1.2229556023      0.0440598937
-7.2470157373      2.4374213225      0.0651311769
 -5.0398812618      2.1435406993      0.1155201154
 -7.2001330967      0.3718768293     -0.0703451879]',
  AtomOnlySymbol.(["C", "O", "O", "H", "H"]) )

 Ar = Cluster{AtomOnlySymbol}(rand(3), AtomOnlySymbol.(["Ar"]))

open(xyzname,"w") do io
    print_xyz(io, formic_acid)
end

input1=load_clusters_and_make_input(xyzname, Ar, ca)
inputs=load_clusters_and_sample_input(xyzname, Ar, ca, 2)

data1=calculate_adaptive_sample_inputs(inputs, save_file_name=fname)
data2=calculate_with_different_method(fname,ca,save_file=sname, restart_file=rname)
data3=continue_calculation(rname,ca, save_file=sname, restart_file=rname)

@test all(isapprox.(data1["Energy"], data2["Energy"], atol=2E-6))
@test all(isapprox.(data1["Energy"], data3["Energy"], atol=2E-6))
@test all(isapprox.(data2["Energy"], data3["Energy"], atol=2E-6))

ldata=load_jld_data(sname)
for fn in [fname, sname]
    ldata=load_jld_data(fn)
    @test length(ldata["cluster1"]) + length(ldata["cluster2"]) == length(ldata["Points"][1])
end

rm(fname)
rm(rname)
rm(sname)
rm(xyzname)
end
