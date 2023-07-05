# Install

Start Julia and hit "]" to enter into pkg REPL, after that type

```julia
pkg> add registry add https://github.com/MatrixLabTools/PackageRegistry
pkg> add PotentialCalculation
```

Currently there are two backends, [ORCA](https://orcaforum.kofo.mpg.de)
and [Psi4](http://www.psicode.org/). To do any calculation, you need to have
at least one of these installed. With Psi4 you need to also load [Psi4Calculator](https://github.com/MatrixLabTools/Psi4Calculator.jl)

## Testing installation

```julia
pkg> test PotentialCalculation
```

checking for backends is made:

- Search PATH for `orca`-executable, for ORCA backend
- Trying to import psi4 python module
If either of these succeeds then that backend is tested.

!!! note "Note"
    When doing calculations with ORCA you can specify `orca`-binary location and
    thus it is not necessary for it to be in PATH. Only testing need it to be in
    PATH. Also using ORCA to do parallel calculations (using ORCA with more than
    one core) requires explicit specification of `orca`-binary location.
