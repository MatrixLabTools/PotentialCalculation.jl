using Documenter
using PotentialCalculation, PotentialCalculation.psi4

makedocs(sitename="PotentialCalculation.jl",
         pages=["Home" => "index.md",
                "Install" => "install.md",
                "Usage" => "use.md",
                "References" => "reference.md"]

)

deploydocs(
    repo = "github.com/MatrixLabTools/PotentialCalculation.jl.git",
)
