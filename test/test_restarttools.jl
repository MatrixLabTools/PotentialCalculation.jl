using Test

using Distributed

addprocs(2)
@everywhere using PotentialCalculation
@everywhere using PotentialCalculation.psi4

fname = tempname()
rname = tempname()
sname = tempname()
xyzname = tempname()

formic_acid=Cluster{AtomOnlySymbol}(
 [-6.7041359778      1.3501192944      0.0102209137
-5.3688853815      1.2229556023      0.0440598937
-7.2470157373      2.4374213225      0.0651311769
 -5.0398812618      2.1435406993      0.1155201154
 -7.2001330967      0.3718768293     -0.0703451879]',
  AtomOnlySymbol.(["C", "O", "O", "H", "H"]) )

Ar = Cluster{AtomOnlySymbol}(rand(3), AtomOnlySymbol.(["Ar"]))
N2 = Cluster{AtomOnlySymbol}([[0.0 1.0]; [0.0 0.0]; [0.0 0.0]], AtomOnlySymbol.(["N", "N"]))

open(xyzname,"w") do io
    print_xyz(io, formic_acid)
end

pbar=true



@testset "Orca" begin
    ca = Calculator{Orca}("blyp d3bj", "def2-svp", Orca())

    input1=load_clusters_and_make_input(xyzname, Ar, ca)
    inputs=load_clusters_and_sample_input(xyzname, N2, ca, 2, npoints=5)
    inputss=load_clusters_and_sample_input(xyzname, xyzname, ca, 2)

    data1=calculate_adaptive_sample_inputs(inputs, save_file_name=fname, pbar=pbar)
    data2=calculate_with_different_method(fname,ca,save_file=sname, restart_file=rname, pbar=pbar)
    data3=continue_calculation(rname,ca, save_file=sname, restart_file=rname, pbar=pbar)

    @test all(isapprox.(data1["Energy"], data2["Energy"], atol=2E-6))
    @test all(isapprox.(data1["Energy"], data3["Energy"], atol=2E-6))
    @test all(isapprox.(data2["Energy"], data3["Energy"], atol=2E-6))
end

@testset "Psi4" begin
    ca = Calculator{Psi4}("blyp-d3bj", "def2-svp",Psi4(memory="1000MiB", nthreads=2))

    input1=load_clusters_and_make_input(xyzname, Ar, ca)
    inputs=load_clusters_and_sample_input(xyzname, N2, ca, 2, npoints=5)
    inputss=load_clusters_and_sample_input(xyzname, xyzname, ca, 2)

    data1=calculate_adaptive_sample_inputs(inputs, save_file_name=fname, pbar=pbar)
    data2=calculate_with_different_method(fname,ca,save_file=sname, restart_file=rname, pbar=pbar)
    data3=continue_calculation(rname,ca, save_file=sname, restart_file=rname, pbar=pbar)

    @test all(isapprox.(data1["Energy"], data2["Energy"], atol=2E-6))
    @test all(isapprox.(data1["Energy"], data3["Energy"], atol=2E-6))
    @test all(isapprox.(data2["Energy"], data3["Energy"], atol=2E-6))
end

@testset "restarttools" begin
    savedata = load_data_file(sname)

    ldata=load_jld_data(sname)
    for fn in [fname, sname]
        ldata=load_jld_data(fn)
        @test length(ldata["cluster1"]) + length(ldata["cluster2"]) == length(ldata["Points"][1])
    end
end

rm(fname)
rm(rname)
rm(sname)
rm(xyzname)
