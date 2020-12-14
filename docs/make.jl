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
        "Tutorials" => ["Tutorials/soundsoft_scattering.md",]
        "Modules"   => ["Utils/Utils.md",
                        "Geometry/Geometry.md"]
        # "Examples" => "examples.md"
        # "Documentation" => "devdocs.md"
    ]
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
# deploydocs(
#     repo = "github.com/WaveProp/WaveProp.jl.git"
# )
