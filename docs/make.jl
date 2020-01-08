using Documenter
using PotentialCalculation

makedocs(sitename="PotentialCalculation.jl",
         pages=["Home" => "index.md"]

)

deploydocs(
    repo = "github.com/MatrixLabTools/PotentialCalculation.jl.git",
)
