# PotentialCalculation.jl
[![codecov](https://codecov.io/gh/tjjarvinen/PotentialCalculation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/tjjarvinen/PotentialCalculation.jl)

A Julia package to calculate potential energy between two molecules.

Fitting of potential energy is done by separate [package](https://github.com/tjjarvinen/PotentialFitting.jl).

Currently only supported backend is [ORCA](https://orcaforum.kofo.mpg.de).

### Installation

```julia
julia> Pkg.add("https://github.com/tjjarvinen/PotentialCalculation.jl")
```


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

#Creating calculation method
mp2 = Calculator("RI-MP2 RIJK",
                 "aug-cc-pVTZ aug-cc-pVTZ/C def2/JK TIGHTSCF",
                  Orca())

# Creating argon atom
Ar = Cluster{AtomOnlySymbol}(rand(1,3), AtomOnlySymbol.(["Ar"]))

# File where other molecule is taken
trajfile="Some trajectory.xyz"

# Create input for calculations
inputs = load_clusters_and_sample_input(trajfile, Ar, mp2,
               32,         #How many lines are calculated
               max_e=5000, #Maximum energy in cm⁻¹ -
                           #   limit to areas where energy is less than this
               npoints=50) #Number of points per line

# Do calculation
data = calculate_adaptive_sample_inputs(inputs, save_file_name="save file")
```

New calculations with different method can be done on previous points.

```julia
#New method
ccf12 = Calculator("CCSD(T)-F12/RI",
                   "cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/C TIGHTSCF",
                    Orca(maxmem=3500))

data2 = calculate_with_different_method("previous results file",
                                        ccf12,
                                        save_file="new save file",
                                        restart_file="restart file")
```

To restart calculations from restart file.

```julia
# Create calculator.
# Method and basis does not matter as they are read from restart file
cal = Calculator("", "", Orca(maxmem=3500))

data3 = continue_calculation("restart file", cal,
                             save_file="final save file",
                             restart_file="new restart file")
```
