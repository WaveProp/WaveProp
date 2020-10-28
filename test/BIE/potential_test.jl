using Test
using WaveProp
using WaveProp.BIE
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Potential test" begin
    pde  = Helmholtz(;dim=3,k=1)
    Ω,mesh  = WaveProp.IO.gmsh_sphere(dim=2)
    compute_quadrature!(mesh,order=1,dim=2,need_normal=true)
    𝓢    = SingleLayerPotential(pde,mesh)
    𝓓    = DoubleLayerPotential(pde,mesh)
    @test BIE.kernel_type(𝓢) == BIE.SingleLayer()
    @test BIE.kernel_type(𝓓) == BIE.DoubleLayer()
    σ    = γ₀(x->0.,mesh)
    u(x) = 𝓢[σ](x)
    x₀ = Point(1.,1.,1.)
    @test u(x₀) ≈ 0
end