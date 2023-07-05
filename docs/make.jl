using Documenter
using PotentialCalculation

makedocs(sitename="PotentialCalculation.jl",
         pages=["Home" => "index.md",
                "Install" => "install.md",
                "Building Molecules" => "structure_creation.md",
                "Usage" => "use.md",
                "References" => "reference.md"]

)

deploydocs(
    repo = "github.com/MatrixLabTools/PotentialCalculation.jl.git",
)
