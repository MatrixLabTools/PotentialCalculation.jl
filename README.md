# PotentialCalculation.jl
[![Build Status](https://travis-ci.org/MatrixLabTools/PotentialCalculation.jl.svg?branch=master)](https://travis-ci.org/MatrixLabTools/PotentialCalculation.jl)
[![codecov](https://codecov.io/gh/MatrixLabTools/PotentialCalculation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MatrixLabTools/PotentialCalculation.jl)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://MatrixLabTools.github.io/PotentialCalculation.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://MatrixLabTools.github.io/PotentialCalculation.jl/dev/)

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
mp2 = Calculator("RI-MP2 RIJK TIGHTSCF",
                 "aug-cc-pVTZ aug-cc-pVTZ/C def2/JK",
                  Orca()
                )

# Creating argon atom
Ar = Cluster(rand(3), AtomOnlySymbol("Ar"))

# File where other molecule is taken
trajfile="Some trajectory.xyz"

# Create input for calculations
inputs = createinputs(trajfile, Ar, mp2;
               nsaples=32,  # How many lines are calculated
               max_e=10000, # Maximum energy in cm⁻¹ -
                            #   limit to areas where energy is less than this
               npoints=50)  # Number of points per line

# Do calculation
data = calculate_potential(inputs, save_file="save file")
```

New calculations with different method can be done on previous points.

```julia
#New method
ccf12 = Calculator("CCSD(T)-F12/RI TIGHTSCF",
                   "cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/C",
                   Orca(maxmem=3500)
                  )

data2 = calculate_potential("previous results file",
                            ccf12,  # Calculate with this method now
                            save_file="new save file",
                            restart_file="restart file"
                          )
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
