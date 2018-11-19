module PotentialCalculation

export m_electron, m_u, m_au, proton_mass, masses,
       AbstractAtom, AbstractAtomWithMass,
       AtomOnlySymbol, AtomWithMass,
       AbstractMolecule, Molecule, MoleculeIdenticalInformation,
       AbstractCluster, AbstractClusterWithSymbols,
         ClusterNoSymbols, Cluster, ClusterWithMass, distances,
         center_coordinates, move!, center_cluster!,
         rotate_x!, rotate_y!, rotate_z!, print_xyz,
       AbstractIdentical, Identical, areidentical,

       read_h5file, give_radius, give_energy,
       energy_from, energy_to, change_energy,

       AbstactCalculationType, Energy, Gradient, AbstractCalculator, Orca,
       Calculator, write_input, read_energy, calculate_energy, clean_calculation_files,
       bsse_corrected_energy,
       line_sampler, adaptive_line_sampler, sample_multiple_adaptive_lines


include("identical.jl")
include("atoms.jl")
include("molecules.jl")
include("clusters.jl")

include("unitconversions.jl")
include("fileaccess.jl")

include("calculators.jl")
include("sample.jl")


using .atoms
using .molecules
using .clusters
using .identical
using .unitconversions
using .fileaccess
using .calculators
using .sample


greet() = print("Hello World!")


#ca = Calculator("RI-MP2 RIJK", "def2-svp def2-svp/C def2/JK", Orca(), Energy())
#c1 = Cluster{AtomOnlySymbol}([0.0 0.0 0.0; 1.2 0.0 0.0 ], AtomOnlySymbol.(["N", "N"]))
#c2 = Cluster{AtomOnlySymbol}(rand(1,3), AtomOnlySymbol.(["Ar"]))
# log1 = Logging.SimpleLogger(stdout,Logging.Debug)
# Logging.global_logger(log1)

end # module
