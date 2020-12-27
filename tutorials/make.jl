using Pkg

Pkg.activate(@__DIR__)
Pkg.update()
Pkg.instantiate()

using Literate

if isempty(ARGS)
    flist = ["sphere_scattering.jl", "singular_integrals.jl"]
else
    flist = ARGS
end    

for file in flist
    fname =  joinpath(@__DIR__,file)
    Literate.markdown(fname,@__DIR__,documenter=true)    
    # Literate.notebook(fname,@__DIR__;documenter=true,execute=false)    
end    