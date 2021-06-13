using Test
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Potential test" begin
    pde  = Helmholtz(;dim=3,k=1)
    Ω,M  = WaveProp.IO.gmsh_sphere(dim=2)
    Γ    = boundary(Ω)
    mesh = NystromMesh(view(M,Γ),order=1)
    𝓢    = SingleLayerPotential(pde,mesh)
    𝓓    = DoubleLayerPotential(pde,mesh)
    @test Nystrom.kernel_type(𝓢) == Nystrom.SingleLayer()
    @test Nystrom.kernel_type(𝓓) == Nystrom.DoubleLayer()
    σ    = γ₀(x->0.,mesh)
    u(x) = 𝓢[σ](x)
    x₀ = SVector(1.,1.,1.)
    @test u(x₀) ≈ 0
end
