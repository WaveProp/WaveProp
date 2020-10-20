using Test
using WaveProp.BIE

@testset "Kernels" begin
    pde  = Helmholtz(;dim=3,k=1)
    G    = SingleLayerKernel(pde)
    dG    = DoubleLayerKernel(pde)
    @test BIE.kernel_type(G) == BIE.SingleLayer()
    @test BIE.kernel_type(dG) == BIE.DoubleLayer()
    @test BIE.return_type(G) == ComplexF64
    @test BIE.return_type(dG) == ComplexF64
    @test BIE.combined_field_coefficients(G) == (0,-1)
    @test BIE.combined_field_coefficients(dG) == (1,0)

    pde  = Laplace(;dim=3)
    G    = SingleLayerKernel(pde)
    dG    = DoubleLayerKernel(pde)
    @test BIE.kernel_type(G) == BIE.SingleLayer()
    @test BIE.kernel_type(dG) == BIE.DoubleLayer()
    @test BIE.return_type(G) == Float64
    @test BIE.return_type(dG) == Float64
    @test BIE.combined_field_coefficients(G) == (0,-1)
    @test BIE.combined_field_coefficients(dG) == (1,0)    
end