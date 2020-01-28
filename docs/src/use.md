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
are used, then random points of them is used.



New calculations with different method can be done on previous points.

```julia
#New method
ccf12 = Calculator(
           "CCSD(T)-F12/RI TIGHTSCF",
           "cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/C",
           Orca(maxmem=3500)
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

### Calculators

PotentilaCalculation can use either [ORCA](https://orcaforum.kofo.mpg.de)
or [Psi4](http://www.psicode.org/) as a backend for calculations.

To create ORCA calculator you can use

```@docs
Orca
```

For Psi4 use

```@docs
Psi4
```
