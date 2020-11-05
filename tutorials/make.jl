using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

using Literate

flist = ["sphere_scattering.jl", "singular_integrals.jl"]

for file in flist
    fname =  joinpath(@__DIR__,file)
    Literate.markdown(fname,@__DIR__)    
    # Literate.notebook(fname,@__DIR__)    
end    