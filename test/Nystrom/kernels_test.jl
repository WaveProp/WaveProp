using Test
using WaveProp.Nystrom

@testset "Kernels" begin
    pde  = Helmholtz(;dim=3,k=1)
    G    = SingleLayerKernel(pde)
    dG    = DoubleLayerKernel(pde)
    @test Nystrom.kernel_type(G) == Nystrom.SingleLayer()
    @test Nystrom.kernel_type(dG) == Nystrom.DoubleLayer()
    @test Nystrom.return_type(G) == ComplexF64
    @test Nystrom.return_type(dG) == ComplexF64
    @test Nystrom.combined_field_coefficients(G) == (0,-1)
    @test Nystrom.combined_field_coefficients(dG) == (1,0)

    pde  = Laplace(;dim=3)
    G    = SingleLayerKernel(pde)
    dG    = DoubleLayerKernel(pde)
    @test Nystrom.kernel_type(G) == Nystrom.SingleLayer()
    @test Nystrom.kernel_type(dG) == Nystrom.DoubleLayer()
    @test Nystrom.return_type(G) == Float64
    @test Nystrom.return_type(dG) == Float64
    @test Nystrom.combined_field_coefficients(G) == (0,-1)
    @test Nystrom.combined_field_coefficients(dG) == (1,0)    
end