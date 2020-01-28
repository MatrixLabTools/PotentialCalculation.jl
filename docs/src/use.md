# Using PotentialCalculation

## Basic usage

Load package

```julia
using PotentialCalculation
```

or if using parallelization (recommended).

```julia
using Distributed
@everywhere using PotentialCalculation
```

!!! note "Note"
    Do not forget to use `@everywhere` macro when importing PotentialCalculation!


Creating inputs and doing basic calculation, where two molecules are calculated
with various distances and orientations from each other.

```julia
# Creating calculation method
mp2 = Calculator(
        "RI-MP2 RIJK TIGHTSCF",
        "aug-cc-pVTZ aug-cc-pVTZ/C def2/JK",
         Orca()
      )

# Creating argon atom
Ar = Cluster(zeros(3), AtomOnlySymbol("Ar"))

# File where other molecule is taken
trajfile="Some trajectory.xyz"

# Create input for calculations
inputs = create_inputs(
            trajfile,    # First molecule will be picked from here
            Ar,          # Second molecule will be picked from here
            mp2;         # Calculator to do electronic structure calculations
            nsaples=32,  # How many lines are calculated
            max_e=10000, # Maximum energy in cm⁻¹ -
                         #   limit to areas where energy is less than this
            npoints=50   # Number of points per line
        )  

# Do calculation
data = calculate_potential(inputs, save_file="save file")
```

The molecules used in the input can be created by hand or read from
xyz-trajectory (recommended). If trajectory file (or array of [`Cluster`](@ref))
are used, then random points of them are picked for the calculation.



New calculations with different method can be done on previous points.

```julia
#New method
ccf12 = Calculator(
           "CCSD(T)-F12/RI TIGHTSCF",
           "cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/C",
           Orca(maxmem=4000)
        )

data2 = calculate_potential(
          "previous results file",
          ccf12,  # Calculate with this method now
          save_file="new save file",
          restart_file="restart file"
       )
```

To restart calculations from restart file.

```julia
data3 = continue_calculation(
          "restart file",
          Orca(maxmem=4000),
          save_file="final save file",
          restart_file="new restart file"
        )
```

## Calculators

PotentilaCalculation can use either [Orca](@ref)
or [Psi4](@ref) as a backend for calculations.

To create ORCA calculator you can use

```julia
Orca(
   executable="path to orca binary",
   maxmem=4000,
   ncore=1,
   tmp_dir="directory where calculations are done"
)
```

If you are using `ncore>1` then you need specify `executable`, in other cases
PATH is searched for orca-binary.  

For Psi4 use

```julia
using PotentialCalculation.psi4

Psi4(
   memory="1GiB",
   nthreads=1,
)
```
You need to import Psi4 explicitly with `using PotentialCalculation.psi4`. All
Psi4 global environmental variables are present. To access them you need to use
`PotentialCalculation.psi4.gPsi4`-handel after using
`PotentialCalculation.psi4.initpsi()`-function to initialize Psi4 environment.

### Adding Calculation Method and Basis Set

Backend calculators are then combined with calculation method and basis set,
to create Calculator structure that has all the information needed to perform
electronic structure calculations.

This is done with
```julia
Calculator(
   "Method",
   "Basis",
   Orca()    # Backend program
)
```

## General Structure of Calculations

### 1-Sampling

Fist phase of calculations is sampling of calculations points. It is recommended
that is is done with some cheap method, like RI-MP2 or DFT, then check that the
sampling is acceptable, use visualization tools from
[PotentialFitting](https://github.com/MatrixLabTools/PotentialFitting.jl) to
do this.

General points are:
- `npoints` needs to be more than 50 for potential to catch rapid changes around
  minimum.
- `max_e` defines usability of the potential, if you just need simple spectrum
  5000 cm⁻¹ is probably enough, but for other cases you might need more.
  So it is recommended to use at least the default 15000 cm⁻¹ value.
- Generate first trajectory with some variability in molecule geometry.
  PotentialCalculation can then sample points out of it, so that the calculated
  potential energy surface has greater applicability.
- Default search parameters have been tested with several molecules and are
  probably fine.

### 2-Calculation

In second phase of calculation final potential is generated. To do this, it is
recommended to use good method. For small molecules CCSD(T)-F12, with
basis set limit approximation, is doable and thus recommended. That is,
when using ORCA, you should use something like

```julia
ccf12 = Calculator(
           "CCSD(T)-F12/RI TIGHTSCF",
           "cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/C",
           Orca(maxmem=4000)
        )
```

For larger molecules ORCA has DLPNO-CCSD(T)-F12, but in the limited testing that
method had some jumps in the potential and was thus deemed to not be a good one.
But it might be ok, with some tuning of DLPNO-parameters.

### 3-Fitting

Third phase is fitting of potential and that is done with
[PotentialFitting](https://github.com/MatrixLabTools/PotentialFitting.jl).

### 4-Storing/Sharing potential

To store or share potential you can use [PotentialDB](https://github.com/MatrixLabTools/PotentialDB.jl).


## Using with SLURM/PBS etc.

Julia support use of varios scheduling software like [Slurm](https://www.schedmd.com/)
or [PBS](https://www.pbspro.org/) etc. trough
[ClusterManagers.jl](https://github.com/JuliaParallel/ClusterManagers.jl)
package. It is recommended that you use any these then you should use
ClusterMangers.

To PotentialCalculation with Slurm simply start with

```julia
using Distributed
using ClusterManagers

addprocs_slurm(number_of_processes) # ntasks option in Slurm job file
@everywhere using PotentialCalculation
```
