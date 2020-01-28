var documenterSearchIndex = {"docs":
[{"location":"install/#Install-1","page":"Install","title":"Install","text":"","category":"section"},{"location":"install/#","page":"Install","title":"Install","text":"Start Julia and hit \"]\" to enter pkg REPL then type","category":"page"},{"location":"install/#","page":"Install","title":"Install","text":"pkg> add registry add https://github.com/MatrixLabTools/PackageRegistry\npkg> add PotentialCalculation","category":"page"},{"location":"install/#","page":"Install","title":"Install","text":"Currently there are two backends, ORCA and Psi4. To do any calculation you need to have at least one of these installed.","category":"page"},{"location":"install/#Testing-installation-1","page":"Install","title":"Testing installation","text":"","category":"section"},{"location":"install/#","page":"Install","title":"Install","text":"pkg> test PotentialCalculation","category":"page"},{"location":"install/#","page":"Install","title":"Install","text":"checking for backends is made:","category":"page"},{"location":"install/#","page":"Install","title":"Install","text":"Search PATH for orca-executable, for ORCA backend\nTrying to import psi4 python module","category":"page"},{"location":"install/#","page":"Install","title":"Install","text":"If either of these succeeds then k.o. backend is tested.","category":"page"},{"location":"install/#","page":"Install","title":"Install","text":"note: Note\nWhen doing calculations with ORCA you can specify orca-binary location and thus it is not necessary for it to be in PATH. Only testing need it to be in PATH. Also using ORCA to do parallel calculations (using ORCA with more than one core) requires explicit specification of orca-binary location.","category":"page"},{"location":"reference/#Reference-1","page":"References","title":"Reference","text":"","category":"section"},{"location":"reference/#","page":"References","title":"References","text":"Modules = [PotentialCalculation.clusters,\n           PotentialCalculation.calculators,\n           PotentialCalculation.psi4,\n           PotentialCalculation.unitconversions,\n           PotentialCalculation.atoms,\n           PotentialCalculation.restarttools]","category":"page"},{"location":"reference/#PotentialCalculation.clusters.Cluster","page":"References","title":"PotentialCalculation.clusters.Cluster","text":"Cluster{T} <: AbstractClusterWithSymbols{T} where T<:AbstractAtom\n\nStructure to hold location data of clusters/molecules\n\nFields\n\nxyz::Array{Float64,2} : location of atoms in 3d space, first index is x, y, z coordinate\natoms::Vector{T} : atom type information\n\n\n\n\n\n","category":"type"},{"location":"reference/#PotentialCalculation.clusters.ClusterNoSymbols","page":"References","title":"PotentialCalculation.clusters.ClusterNoSymbols","text":"ClusterNoSymbols <: AbstractCluster\n\nBasic cluster has only location of atoms.\n\nFields\n\nxyz::Array{Float64,2} : location of atoms in 3d space, first index is x, y, z coordinate\n\n\n\n\n\n","category":"type"},{"location":"reference/#PotentialCalculation.clusters.center_cluster!-Tuple{AbstractCluster}","page":"References","title":"PotentialCalculation.clusters.center_cluster!","text":"center_cluster!(c::AbstractCluster)\n\nCenters cluster to origin of coordinates\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.center_coordinates-Tuple{AbstractCluster}","page":"References","title":"PotentialCalculation.clusters.center_coordinates","text":"center_coordinates(c::AbstractCluster)\n\nGives coordinates to aricmetric mean of clusters atoms\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.cluster_angle-Tuple{AbstractCluster,Any,Any,AbstractCluster,Any}","page":"References","title":"PotentialCalculation.clusters.cluster_angle","text":"cluster_angle(c1::AbstractCluster, i, j, c2::AbstractCluster, k)\n\nCalculates angle (radians) between atons in different clusters\n\nArguments\n\nc1::AbstractCluster : first cluster\ni : index in c1\nj : index in c1\nc2::AbstractCluster : second cluster\nk : index in c2\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.cluster_angle-Tuple{AbstractCluster,Any,Any,Any}","page":"References","title":"PotentialCalculation.clusters.cluster_angle","text":"cluster_angle(c::AbstractCluster, i, j, k)\n\nCalculates angle (radians) between atoms i,j,k in cluster\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.dihedral_angle-Tuple{AbstractCluster,Any,Any,Any,Any}","page":"References","title":"PotentialCalculation.clusters.dihedral_angle","text":"clusterdihedralangle(c::AbstractCluster, i, j, k, m)\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.distances-Tuple{AbstractCluster,AbstractCluster}","page":"References","title":"PotentialCalculation.clusters.distances","text":"distances(c1::AbstractCluster, c2::AbstractCluster)\n\nReturn distances between atoms of given clusters\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.distances-Tuple{AbstractCluster,AbstractUnitRange,AbstractCluster,AbstractUnitRange}","page":"References","title":"PotentialCalculation.clusters.distances","text":"distances(c1::AbstractCluster, ur1::AbstractUnitRange,\n          c2::AbstractCluster, ur2::AbstractUnitRange)\n\nReturn distance between atoms ur1 in c1 and atoms ur2 in c2\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.distances-Tuple{AbstractCluster,AbstractUnitRange,AbstractUnitRange}","page":"References","title":"PotentialCalculation.clusters.distances","text":"distances(c::AbstractCluster, ar1::AbstractRange, ar2::AbstractRange)\n\nReturns distances between atoms in given unit ranges\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.distances-Tuple{AbstractCluster,Any,AbstractCluster,Any}","page":"References","title":"PotentialCalculation.clusters.distances","text":"distances(c1::AbstractCluster, i, c2::AbstractCluster, j)\n\nReturn distance between atom i in c1 and atom j in c2\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.distances-Tuple{AbstractCluster,Any,Any}","page":"References","title":"PotentialCalculation.clusters.distances","text":"distances(c::AbstractCluster, i, j)\n\nReturns distance of atoms i and j\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.move!-Tuple{AbstractCluster,Any}","page":"References","title":"PotentialCalculation.clusters.move!","text":"move!(c::AbstractCluster,r)\n\nMoves cluster by r\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.print_xyz","page":"References","title":"PotentialCalculation.clusters.print_xyz","text":"print_xyz(io::IO, c::AbstractClusterWithSymbols, note=\"\"; printheader=true)\n\nPrints cluster in xyz file format\n\nArguments\n\nio::IO : stream where writing is done\nc::AbstractClusterWithSymbols : cluster that is writen\nnote=\"\" : message writen on note line\nprintheader=true : wheather or not header is writen (number of atoms and note)\n\n\n\n\n\n","category":"function"},{"location":"reference/#PotentialCalculation.clusters.rotate_randomly!-Tuple{AbstractCluster}","page":"References","title":"PotentialCalculation.clusters.rotate_randomly!","text":"rotate_randomly!(c::AbstractCluster)\n\nRotate cluster by random angle and axis\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.rotate_x!-Tuple{AbstractCluster,Any}","page":"References","title":"PotentialCalculation.clusters.rotate_x!","text":"rotate_x!(c::AbstractCluster, θ)\n\nRotates cluster around x-axis by angle θ\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.rotate_y!-Tuple{AbstractCluster,Any}","page":"References","title":"PotentialCalculation.clusters.rotate_y!","text":"rotate_y!(c::AbstractCluster, θ)\n\nRotates cluster around y-axis by angle θ\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.clusters.rotate_z!-Tuple{AbstractCluster,Any}","page":"References","title":"PotentialCalculation.clusters.rotate_z!","text":"rotate_z!(c::AbstractCluster, θ)\n\nRotates cluster around z-axis by angle θ\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.calculators.Calculator","page":"References","title":"PotentialCalculation.calculators.Calculator","text":"mutable struct Calculator\n\nStruct used to hold together different bits of information needed in calculations.\n\nFields\n\nmethod : calculation method informaton\nbasis : basis information\ncalculator : calculation program information\n\nExamples\n\njulia> Calculator(\"method\", \"basis\", Orca())\nCalculator{Orca}(\"method\", \"basis\", Orca(\"orca\", 0x0000000000000001, 0x00000000000003e8, \"/tmp/tmpVg959k\"))\n\n\n\n\n\n","category":"type"},{"location":"reference/#PotentialCalculation.calculators.Orca","page":"References","title":"PotentialCalculation.calculators.Orca","text":"mutable struct Orca <: AbstractCalculationProgram\n\nHolds information of how to use ORCA program.\n\nNote that if you are using more than 1 core you need to specify exact path to orca executable.\n\nFields\n\nexecutable : path to \"orca\" executable - default = \"orca\"\nncore::Uint : number of cores used by orca - default = 1\nmaxmem:.Uint : maximum memory per core for orca in mega bytes - default = 1000\ntmp_dir : directory where calculations are done - default = mktempdir()\n\nExamples\n\njulia> Orca()\nOrca(\"orca\", 0x0000000000000001, 0x00000000000003e8, \"/tmp/jl_k4hnkY\")\n\njulia> Orca(executable=\"/opt/Orca/bin/orca\", maxmem=4000, tmp_dir=\"/tmp/calculations/\")\nOrca(\"/opt/Orca/bin/orca\", 0x0000000000000001, 0x0000000000000fa0, \"/tmp/calculations/\")\n\n\n\n\n\n","category":"type"},{"location":"reference/#PotentialCalculation.calculators.bsse_corrected_energy-Tuple{Calculator{Orca},Any,Any}","page":"References","title":"PotentialCalculation.calculators.bsse_corrected_energy","text":"bsse_corrected_energy(cal::Calculator, c1, c2; basename=\"base\", id=\"\", pchannel=pchannel)\n\nCalculates energy of combined clusters taking into account basis set superposition error\n\nArguments\n\ncal::Calculator : calcualation information\nc1 : collection of clusters\nc2 : collection of clusters\n\nKeywords\n\nbasename=\"base\" : base name for input/ouput files\nid=\"\" : additional information for calculator - needed for multiprocessing in same folder\npchannel=undef  : (Remote)Channel where progess information is added\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.calculators.bsse_corrected_energy-Tuple{Calculator{Orca},Cluster,Cluster}","page":"References","title":"PotentialCalculation.calculators.bsse_corrected_energy","text":"bsse_corrected_energy(cal::Calculator, c1::Cluster, c2::Cluster; basename=\"base\", id=\"\", pchannel=pchannel)\n\nCalculates energy of combined cluster taking into account basis set superposition error\n\nArguments\n\ncal::Calculator : calcualation information\nc1::Cluster : cluster\nc2::Cluster : cluster\n\nKeywords\n\nbasename=\"base\" : base name for input/ouput files\nid=\"\" : additional information for calculator - needed for multiprocessing in same folder\npchannel=undef  : (Remote)Channel where progess information is added\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.calculators.calculate_energy-Tuple{Calculator{Orca},Any}","page":"References","title":"PotentialCalculation.calculators.calculate_energy","text":"calculate_energy(cal::Calculator, points; basename=\"base\", ghost=undef, id=\"\", pchannel=pchannel)\n\nCalculates energy for given clusters.\n\nArguments\n\ncal::Calculator : calculation information\npoints : collection of Cluster which energy is calculated\n\nKeywords\n\nbasename=\"base\" : base name for input/ouput files\nghost=undef : indices for atoms considered as ghost\nid=\"\" : additional information for calculator - needed for multiprocessing in same folder\npchannel=undef  : (Remote)Channel where progess information is added\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.calculators.calculate_energy-Tuple{Calculator{Orca},Cluster}","page":"References","title":"PotentialCalculation.calculators.calculate_energy","text":"calculate_energy(cal::Calculator, point::Cluster; basename=\"base\", ghost=undef, id=\"\", pchannel=pchannel)\n\nCalculates energy for given cluster.\n\nArguments\n\ncal::Calculator : calculation information\npoint::Cluster : cluster which energy is calculated\n\nKeywords\n\nbasename=\"base\" : base name for input/ouput files\nghost=undef : indices for atoms considered as ghost\nid=\"\" : additional information for calculator - needed for multiprocessing in same folder\npchannel=undef  : (Remote)Channel where progess information is added\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.calculators.clean_calculation_files-Tuple{}","page":"References","title":"PotentialCalculation.calculators.clean_calculation_files","text":"clean_calculation_files(;dir=\".\", basename=\"base\")\n\nDeletes files that calculations produced\n\nArguments\n\ndir=\".\" : directory to be cleaned\nbasename=\"base\" : base name used in calculations\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.calculators.getBSSEsteps-Tuple{AbstractCalculationProgram}","page":"References","title":"PotentialCalculation.calculators.getBSSEsteps","text":"getBSSEsteps(cal::AbstractCalculationProgram)\n\nUsed to set up progressbar. Returns 1, if calculator program has a method to calculate counter poise correction with a single command and 3, if program needs to be called multiple times to make counter poise correction.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.calculators.read_energy-Tuple{Any}","page":"References","title":"PotentialCalculation.calculators.read_energy","text":"read_energy(fname)\n\nRead ORCA output for final energy. Throws an error if calculation failed.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.calculators.write_input-Tuple{IO,Calculator{Orca},AbstractClusterWithSymbols}","page":"References","title":"PotentialCalculation.calculators.write_input","text":"write_input(io::IO, cal::Calculator,\n                 c::AbstractClusterWithSymbols; ghost=undef)\n\nWrites input files for calculation program. Only ORCA input is supported now.\n\nArguments\n\nio::IO : stream where writing is done\ncal::Calculator : calculation information\nc::AbstractClusterWithSymbols : molecule cluster that is calculated\nghost=undef : collection for atom idexes that are considered ghost atoms\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.calculators.print_orca_xyz-Tuple{IO,AbstractClusterWithSymbols}","page":"References","title":"PotentialCalculation.calculators.print_orca_xyz","text":"printorcaxyz(io::IO, c::AbstractClusterWithSymbols; ghost=undef)\n\nPrints ORCA input style xyz input. Idexes in ghost are marked as ghost atoms.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.psi4.Psi4","page":"References","title":"PotentialCalculation.psi4.Psi4","text":"mutable struct Psi4 <: AbstractCalculationProgram\n\nHolds information that calculations are to be done with Psi4\n\nFields\n\nmemory=\"500MiB\" : memory used by Psi4\nnthreads=1  : number of threads used by Psi4\n\n\n\n\n\n","category":"type"},{"location":"reference/#PotentialCalculation.psi4.gPsi4","page":"References","title":"PotentialCalculation.psi4.gPsi4","text":"gPsi4 = undef\n\nHold Psi4 object given by PyCall. Used to do calculation in current process\n\n\n\n\n\n","category":"constant"},{"location":"reference/#PotentialCalculation.psi4.gpsi4init","page":"References","title":"PotentialCalculation.psi4.gpsi4init","text":"gpsi4init=false\n\nGlobal variable to see if Psi4 was initiated for current process\n\n\n\n\n\n","category":"constant"},{"location":"reference/#PotentialCalculation.psi4.initpsi4-Tuple{}","page":"References","title":"PotentialCalculation.psi4.initpsi4","text":"initpsi4(;memory=\"500 MiB\", quiet=true, nthreads=1)\n\nUsed to intialize Psi4 environment\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.unitconversions.change_energy_unit-Tuple{Any,Any,Any}","page":"References","title":"PotentialCalculation.unitconversions.change_energy_unit","text":"change_energy(e, from, to)\n\nChanges energy unit\n\nSuported units include:\n\ncm^-1/cm⁻¹/cm-1\neV\nhartree\nkcal/mol\nkJ/mol\nK\n\nIf unit is not recognized defaults to hartree.\n\nSource is Wikipedia\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.unitconversions.energy_from-Tuple{Any,Any}","page":"References","title":"PotentialCalculation.unitconversions.energy_from","text":"energy_from(e, unit)\n\nChanges given energy unit to hartree\n\nSuported units include:\n\ncm^-1/cm⁻¹/cm-1\neV\nhartree\nkcal/mol\nkJ/mol\nK\n\nIf unit is not recognized defaults to hartree.\n\nSource is Wikipedia\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.unitconversions.energy_to-Tuple{Any,Any}","page":"References","title":"PotentialCalculation.unitconversions.energy_to","text":"energy_to(e, unit)\n\nChanges energy from hartrees to given unit\n\nSuported units include:\n\ncm^-1/cm⁻¹/cm-1\neV\nhartree\nkcal/mol\nkJ/mol\nK\n\nIf unit is not recognized defaults to hartree.\n\nSource is Wikipedia\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.atoms.AtomOnlySymbol","page":"References","title":"PotentialCalculation.atoms.AtomOnlySymbol","text":"AtomOnlySymbol <: AbstractAtom\n\nAtom that only has symbol\n\nFields\n\nid::String : hold description of atom\n\n\n\n\n\n","category":"type"},{"location":"reference/#PotentialCalculation.atoms.AtomWithMass","page":"References","title":"PotentialCalculation.atoms.AtomWithMass","text":"AtomWithMass <: AbstractAtomWithMass\n\nAtom that has symbol and mass\n\nFields\n\nid::String    : hold description of atom\nmass::Float64 : mass of atom\n\n\n\n\n\n","category":"type"},{"location":"reference/#PotentialCalculation.restarttools","page":"References","title":"PotentialCalculation.restarttools","text":"module restarttools\n\nPrimary tools to do calculations. Also contains all methods to restart calculations.\n\n\n\n\n\n","category":"module"},{"location":"reference/#PotentialCalculation.restarttools.calculate_energy_for_xyzfile-Tuple{Any,Any}","page":"References","title":"PotentialCalculation.restarttools.calculate_energy_for_xyzfile","text":"calculate_energy_for_xyzfile(fname, cal; pbar=true)\n\nReads xyz-file and calculates energy of each point on it\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools.calculate_potential-Tuple{AbstractArray{InputAdaptiveSampler,N} where N}","page":"References","title":"PotentialCalculation.restarttools.calculate_potential","text":"calculate_adaptive_sample_inputs(inputs; save_after=\"\", save_after=nworkers(), pbar=true)\n\nUses adaptive_line_sampler to inputs in distributed fashion\n\nArguments\n\ninputs : calculation inputs array (eg. from create_inputs)\nsave_after : file where data saved\nsave_after=nworkers() : number of calculated items before data is saved\npbar=true : show progress bar\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools.calculate_potential-Tuple{AbstractString,Calculator}","page":"References","title":"PotentialCalculation.restarttools.calculate_potential","text":"calculate_potential(fname::AbstractString, calculator::Calculator;\n                    save_file=\"\", restart_file=\"\", pbar=true, save_after=nworkers()\n                   )\n\nWith this function you can use already chosen points on to which to do energy calculation.\n\nArguments\n\nfname : name of save/restart file where points are read\ncalculatror::Calculator : calculator used for calculations\n\nKeywords\n\nsave_file=\"\" : save final results here, if not empty\nrestart_file=\"\" : save restarts here, if not empty\npbar=true : show progress bar\nsave_after=nworkers() : make restart file when given ammount of points is calculated\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools.continue_calculation-Tuple{AbstractString,Calculator}","page":"References","title":"PotentialCalculation.restarttools.continue_calculation","text":"continuecalculation(fname, calculator::Calculator; savefile=\"\", restartfile=\"\",                               saveafter=nworkers(),  pbar=true)\n\nRestarts calculation from given file\n\nArguments\n\nfname : file from which calculation is restarted\ncalculator::Calculator : calculator used in calculations\nsave_file : file where final results are saved, if given\nrestart_file : file where new restart information is saved, if given\nsave_after=nworkers() : make restart file when given ammount of points is calculated\npbar=true : show progress bar\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools.create_inputs-Tuple{AbstractString,AbstractString,Calculator}","page":"References","title":"PotentialCalculation.restarttools.create_inputs","text":"create_inputs(cluster1, cluster2, calculator::Calculator;  kwargs...)\n\nCreates inputs for the calculations. Use calculate_adaptive_sample_inputs on results to do the calculation itself.\n\nInputs for the molecules (named cluster1/2) can be loaded from xyz-file or given by as Cluster structures. If given as file names, then nsamples points are randomly picked from the trajectories for calculation. Total number of lines sampled is thus nlines times nsamples.\n\nParallelisation is done over nsamples, thus it is recommended for it to be multiple of number workers processes.\n\nArguments\n\ncluster1 : something that can be interpreted as Cluster or an Array of it\ncluster2 : something that can be interpreted as Cluster or an Array of it\ncalculator::Calculator : Calculator used in sampling\n\nKeywords\n\nmax_e=15000 : energy treshold in sampling. Enegy is lower than this value.\nmaxdis=9.0 : maximun distance in calculation\nnline=1  : number of lines sampled for given structure pair\nnpoints=10 : number of points in each sampled line\nnsamples=2 : number of points picked randomly from trajectories for line sampling\nsstep=0.1 : step size used in sampling\nstartdistance= 3.5 : starting distance in sampling\nunit=\"cm-1\" : unit for max_e\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools.load_clusters_and_make_input-Tuple{Cluster,Cluster,Any}","page":"References","title":"PotentialCalculation.restarttools.load_clusters_and_make_input","text":"loadclustersandmakeinput(cluster1::Cluster, cluster2::Cluster, calculator;                                       nlines=1, max_e=0, unit=\"cm-1\", npoints=10,                                       maxdis=9.0, sstep=0.1, startdistance=2.5)\n\nLoads cluster1 from xyz-file and returns them all as an Array that can then be used with calculate_adaptive_sample_inputs to calculate energy data. Differs from load_clusters_and_sample_input by taking every point from file\n\nArguments\n\ncluster1::Cluster : cluster1\ncluster2::Cluster : cluster2\nnlines : number of lines to be sampled in calculation\nmax_e : Point that is closest and has less energy than this will be starting point for a line\nunit : Unit in which max_e is given\nnpoints : Number of points in potential\nmaxdis : Maximum distance in potential calculation\nsstep : Search step size for adaptive search\nstartdistance : Distance from which adaptive seach is started\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools.load_clusters_and_make_input-Tuple{String,Cluster,Any}","page":"References","title":"PotentialCalculation.restarttools.load_clusters_and_make_input","text":"loadclustersandmakeinput(cluster1::String, cluster2::Cluster, calculator, nsamples;                                       max_e=0, unit=\"cm-1\", npoints=10,                                       maxdis=9.0, sstep=0.1, startdistance=2.5)\n\nLoads cluster1 from xyz-file and returns them all as an Array that can then be used with calculate_adaptive_sample_inputs to calculate energy data. Differs from load_clusters_and_sample_input by taking every point from file\n\nArguments\n\ncluster1::String : file from where cluster1 is sampled\ncluster2::Cluster :  cluster2\nmax_e : Point that is closest and has less energy than this will be starting point for a line\nunit : Unit in which max_e is given\nnpoints : Number of points in potential\nmaxdis : Maximum distance in potential calculation\nsstep : Search step size for adaptive search\nstartdistance : Distance from which adaptive seach is started\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools.load_clusters_and_sample_input-Tuple{String,Any,Any,Any}","page":"References","title":"PotentialCalculation.restarttools.load_clusters_and_sample_input","text":"loadclustersandsampleinput(fnamecluster1, cluster2, calculator, nsamples;                                       maxe=0, unit=\"cm-1\", npoints=10,                                       maxdis=9.0, sstep=0.1, startdistance=2.5)\n\nLoads cluster1 from xyz-file and takes nsamples samples of it and returns them as an Array that can then be used with calculate_adaptive_sample_inputs to calculate energy data\n\nArguments\n\ncluster1 : file from where cluster1 is sampled\ncluster2 : cluster2\nnsamples : number of samples taken from fname_cluster1\nmax_e : Point that is closest and has less energy than this will be starting point for a line\nunit : Unit in which max_e is given\nnpoints : Number of points in potential\nmaxdis : Maximum distance in potential calculation\nsstep : Search step size for adaptive search\nstartdistance : Distance from which adaptive seach is started\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools.load_clusters_and_sample_input-Tuple{String,String,Any,Any}","page":"References","title":"PotentialCalculation.restarttools.load_clusters_and_sample_input","text":"function loadclustersandsampleinput(cluster1::String, cluster2::String, calculator, nsamples;                                       nlines=1, max_e=0, unit=\"cm-1\", npoints=10,                                       maxdis=9.0, sstep=0.1, startdistance=2.5)\n\nArguments\n\ncluster1 : file from where cluster1 is sampled\ncluster2 : file from where cluster2 is sampled\nnsamples : number of samples taken from fname_cluster1\nnlines : number of lines to be sampled in calculation\nmax_e : Point that is closest and has less energy than this will be starting point for a line\nunit : Unit in which max_e is given\nnpoints : Number of points in potential\nmaxdis : Maximum distance in potential calculation\nsstep : Search step size for adaptive search\nstartdistance : Distance from which adaptive seach is started\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools.load_data_file-Tuple{AbstractString}","page":"References","title":"PotentialCalculation.restarttools.load_data_file","text":"load_data_file(fname)\n\nLoads saved data\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools.load_restart_file-Tuple{AbstractString}","page":"References","title":"PotentialCalculation.restarttools.load_restart_file","text":"load_restart_file(fname)\n\nLoads restart information from file fname and adds to it key not_calculated, which holds information of which collumns of Points have not been calculated.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools.write_restart_file-Tuple{AbstractString,Calculator,Any,Any,Any,Any}","page":"References","title":"PotentialCalculation.restarttools.write_restart_file","text":"write_restart_file(fname, calculator, points, restart_energy, cluster1, cluster2)\n\nSaves restart information for energy calculation.\n\nArguments\n\nfname : name of restartfile\ncalculator::Calculator : calculator used in calculations\npoints : 2d array of point where to calculate energy\nenergy : energy for points that have been caculeted\ncluster1 : cluster1\ncluster2 : cluster2\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools.write_save_file-Tuple{AbstractString,Calculator,Any,Any,Any,Any}","page":"References","title":"PotentialCalculation.restarttools.write_save_file","text":"write_save_file(fname, calculator, points, energy, cluster1, cluster2)\n\nSaves final information for energy calculation.\n\nArguments\n\nfname : name of restartfileusing PotentialCalculation\ncalculator::Calculator : calculator used in calculations\npoints : 2d array of point where to calculate energy\nenergy : energy for points that have been caculeted\ncluster1 : cluster1\ncluster2 : cluster2\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools._calculate_points-Tuple{Any,Any,Any}","page":"References","title":"PotentialCalculation.restarttools._calculate_points","text":"_calculate_points(cal, c1_points, c2_points; pchannel=undef)\n\nInternal helper function to ease use of distributed calculations. Users should call \"calculate_points\" instead.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PotentialCalculation.restarttools._sample_and_calculate-Tuple{InputAdaptiveSampler}","page":"References","title":"PotentialCalculation.restarttools._sample_and_calculate","text":"_sample_and_calculate(inputs; pchannel=undef )\n\nInternal helper function to help implement distributed calculations\n\n\n\n\n\n","category":"method"},{"location":"#PotentialCalculation.jl-1","page":"Home","title":"PotentialCalculation.jl Documentation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"PotentialCalculation, PotentialFitting and PotentialDB together form a group of packages that can be used to calculate, fit and store molecular potentials that can be used on QM-MM simulations of noble gas matrix isolation experiments.","category":"page"},{"location":"#[PotentialCalculation](https://github.com/MatrixLabTools/PotentialCalculation.jl)-Features-1","page":"Home","title":"PotentialCalculation Features","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Calculate potentials with ORCA or Psi4\nAutomatic sampling of calculation points\nSupports parallelisation of calculation across compute nodes","category":"page"},{"location":"#[PotentialFitting](https://github.com/MatrixLabTools/PotentialFitting.jl)-Features-1","page":"Home","title":"PotentialFitting Features","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Define new potentials easily\nFit potential using methods provided by ScikitLearn","category":"page"},{"location":"#[PotentialDB](https://github.com/MatrixLabTools/PotentialDB.jl)-Features-1","page":"Home","title":"PotentialDB Features","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Store computed data\nPublic storage to computed data\nEasily manage calculated potentials","category":"page"},{"location":"#Contents-1","page":"Home","title":"Contents","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Pages = [\"install.md\", \"use.md\", \"reference.md\"]","category":"page"},{"location":"use/#Using-PotentialCalculation-1","page":"Usage","title":"Using PotentialCalculation","text":"","category":"section"},{"location":"use/#Basic-usage-1","page":"Usage","title":"Basic usage","text":"","category":"section"},{"location":"use/#","page":"Usage","title":"Usage","text":"Load package","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"using PotentialCalculation","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"or if using parallelization (recommended).","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"using Distributed\n@everywhere using PotentialCalculation","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"note: Note\nDo not forget to use @everywhere macro when importing PotentialCalculation!","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"Creating inputs and doing basic calculation, where two molecules are calculated with various distances and orientations from each other.","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"\n# Creating calculation method\nmp2 = Calculator(\n        \"RI-MP2 RIJK TIGHTSCF\",\n        \"aug-cc-pVTZ aug-cc-pVTZ/C def2/JK\",\n         Orca()\n      )\n\n# Creating argon atom\nAr = Cluster(zeros(3), AtomOnlySymbol(\"Ar\"))\n\n# File where other molecule is taken\ntrajfile=\"Some trajectory.xyz\"\n\n# Create input for calculations\ninputs = create_inputs(\n            trajfile,    # First molecule will be picked from here\n            Ar,          # Second molecule will be picked from here\n            mp2;         # Calculator to do electronic structure calculations\n            nsaples=32,  # How many lines are calculated\n            max_e=10000, # Maximum energy in cm⁻¹ -\n                         #   limit to areas where energy is less than this\n            npoints=50   # Number of points per line\n        )  \n\n# Do calculation\ndata = calculate_potential(inputs, save_file=\"save file\")","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"The molecules used in the input can be created by hand or read from xyz-trajectory (recommended). If trajectory file (or array of Cluster) are used, then random points of them is used.","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"New calculations with different method can be done on previous points.","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"#New method\nccf12 = Calculator(\n           \"CCSD(T)-F12/RI TIGHTSCF\",\n           \"cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/C\",\n           Orca(maxmem=4000)\n        )\n\ndata2 = calculate_potential(\n          \"previous results file\",\n          ccf12,  # Calculate with this method now\n          save_file=\"new save file\",\n          restart_file=\"restart file\"\n       )","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"To restart calculations from restart file.","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"data3 = continue_calculation(\n          \"restart file\",\n          Orca(maxmem=4000),\n          save_file=\"final save file\",\n          restart_file=\"new restart file\"\n        )","category":"page"},{"location":"use/#Calculators-1","page":"Usage","title":"Calculators","text":"","category":"section"},{"location":"use/#","page":"Usage","title":"Usage","text":"PotentilaCalculation can use either Orca or Psi4 as a backend for calculations.","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"To create ORCA calculator you can use","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"Orca(\n   executable=\"path to orca binary\",\n   maxmem=4000,\n   ncore=1,\n   tmp_dir=\"directory where calculations are done\"\n)","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"If you are using ncore>1 then you need specify executable, in other cases PATH is searched for orca-binary.  ","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"For Psi4 use","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"using PotentialCalculation.psi4\n\nPsi4(\n   memory=\"1GiB\",\n   nthreads=1,\n)","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"You need to import Psi4 explicitly with using PotentialCalculation.psi4. All Psi4 global environmental variables are present. To access them you need to use PotentialCalculation.psi4.gPsi4-handel after using PotentialCalculation.psi4.initpsi()-function to initialize Psi4 environment.","category":"page"},{"location":"use/#Using-with-SLURM/PBS-etc.-1","page":"Usage","title":"Using with SLURM/PBS etc.","text":"","category":"section"},{"location":"use/#","page":"Usage","title":"Usage","text":"Julia support use of varios scheduling software like Slurm or PBS etc. trough ClusterManagers.jl package. It is recommended that you use any these then you should use ClusterMangers.","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"To PotentialCalculation with Slurm simply start with","category":"page"},{"location":"use/#","page":"Usage","title":"Usage","text":"using Distributed\nusing ClusterManagers\n\naddprocs_slurm(number_of_processes) # ncore option in Slurm job file\n@everywhere using PotentialCalculation","category":"page"}]
}
