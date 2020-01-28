# [PotentialCalculation.jl Documentation](@id PotentialCalculation.jl)


[PotentialCalculation](https://github.com/MatrixLabTools/PotentialCalculation.jl),
[PotentialFitting](https://github.com/MatrixLabTools/PotentialFitting.jl)
and [PotentialDB](https://github.com/MatrixLabTools/PotentialDB.jl) together
form a group of packages that can be used to calculate, fit and store molecular
potentials, used on QM-MM simulations of noble gas matrix isolation
experiments.

## Features

### [PotentialCalculation](https://github.com/MatrixLabTools/PotentialCalculation.jl)

- Calculate potentials with [ORCA](https://orcaforum.kofo.mpg.de) or [Psi4](http://www.psicode.org/)
- Automatic sampling of calculation points
- Supports parallelisation of calculation across compute nodes

### [PotentialFitting](https://github.com/MatrixLabTools/PotentialFitting.jl)

- Define new potentials easily
- Fit potential using methods provided by [ScikitLearn](https://github.com/cstjean/ScikitLearn.jl/)

### [PotentialDB](https://github.com/MatrixLabTools/PotentialDB.jl)

- Store computed data
- Public storage to computed data
- Easily manage calculated potentials


## Contents

```@contents
Pages = ["install.md", "use.md", "reference.md"]
```
