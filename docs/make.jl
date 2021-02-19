using Documenter
using WaveProp

makedocs(
    sitename = "WaveProp",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [WaveProp],
    pages = [
        "Home"      => ["index.md",]
        # "Tutorials" => ["Tutorials/soundsoft_scattering.md",]
        "Modules"   => ["Utils/Utils.md",
                        "Geometry/Geometry.md",
                        "Interpolation/Interpolation.md",
                        "Integration/Integration.md",
                        "Mesh/Mesh.md",
                        "PhysicalProblem/PhysicalProblem.md",
                        "FEM/FEM.md",
                        "Nystrom/Nystrom.md",
                        "IO/IO.md"
                        ]
        # "Examples" => "examples.md"
        # "Documentation" => "devdocs.md"
    ]
)

# Deploy to gh-pages. See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/WaveProp/WaveProp.jl.git",
    devbranch = "main"
)
