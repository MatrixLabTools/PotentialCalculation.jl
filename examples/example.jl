using Distributed

addprocs(8)


#Adopt or drop these depending your path
@everywhere using Pkg
@everywhere Pkg.activate("..")

@everywhere using PotentialCalculation



mp2 = Calculator("RI-MP2 RIJK", "aug-cc-pVTZ aug-cc-pVTZ/C def2/JK TIGHTSCF", Orca())

Ar = Cluster{AtomOnlySymbol}(rand(1,3), AtomOnlySymbol.(["Ar"]))
inputs = load_clusters_and_sample_input("some trajectory.xyz",
                      Ar, mp2, 32, max_e=5000, npoints=50)

fname="file to save results"
data = calculate_adaptive_sample_inputs(inputs, save_file_name=fname)

#Calculate with different method using same points
ccf12 = Calculator("CCSD(T)-F12/RI", "cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/C TIGHTSCF", Orca(maxmem=3500))
data2 = calculate_with_different_method(fname, ccf12, save_file="final results file", restart_file="restart file")

#If you need to restart calculation
data3=continue_calculation("restart file", ccf12, save_file="final results file", restart_file="new restart file")
