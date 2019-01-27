using Distributed

addprocs(2)
@everywhere using PotentialCalculation


include("tests.jl")
