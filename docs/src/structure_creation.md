# Building Molecules

Input structures are defined by [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl),
but internally calculations are done with a structure called [`Cluster`](@ref).
It is [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl) compatible and
you can convert back and forth and use it where AtomsBase structures can be used.
Atoms base is exported by PotentialCalculation, so there is no need to load it.


While you can create [`Cluster`](@ref) directly, it is recommended that you use
[AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl) to create structures.
You can e.g. use [AtomsIO](https://github.com/mfherbst/AtomsIO.jl)
to load data molecule or just create molecules by hand on the fly.

## Example

```julia
using PotentialCalculation
using AtomsIO

# Load from file with AtomsIO
system = load_system("molecule.xyz")

# Convert to Cluster
csystem = Cluster(system)

# Convert back to AtomsBase standard structure
fsystem = FlexibleSystem(csystem)

# Create structure by hand
N2 = isolated_system( [Atom(:N, [1., 0., 0.].*u"Å"), Atom(:N, [0., 0., 0.].*u"Å")] )
```

