# PotentialCalculation.jl


| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url] [![][codecov-img]][codecov-url] |



A Julia package to calculate potential energy between two molecules.

Fitting of potential energy is done by separate [package](https://github.com/MatrixLabTools/PotentialFitting.jl).

Currently supported backendes are [ORCA](https://orcaforum.kofo.mpg.de)
and [Psi4](http://www.psicode.org/).

### Installation

Hit "]" to enter "pkg>"
```julia
pkg> add registry add https://github.com/MatrixLabTools/PackageRegistry
pkg> add PotentialCalculation
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


[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://MatrixLabTools.github.io/PotentialCalculation.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://MatrixLabTools.github.io/PotentialCalculation.jl/stable

[travis-img]: https://travis-ci.org/MatrixLabTools/PotentialCalculation.jl.svg?branch=master
[travis-url]: https://travis-ci.org/MatrixLabTools/PotentialCalculation.jl

[codecov-img]: https://codecov.io/gh/MatrixLabTools/PotentialCalculation.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/MatrixLabTools/PotentialCalculation.jl
