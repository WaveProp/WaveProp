using Test
using WaveProp
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Circle" begin
    geo = Circle()
    bnd = boundary(geo)
    els = Mesh._meshgen(bnd[1];gridsize=0.1)
    el  = els[1]
    Geometry._global_add_entity(geo)
    Ω,mesh =  Mesh.meshgen(geo;gridsize=0.1)
    # fig = plot()
    # for el in els   
    #     plot!(fig,el;gridsize=0.01)    
    # end
    # display(fig)
end