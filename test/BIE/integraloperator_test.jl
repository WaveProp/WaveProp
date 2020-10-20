using Test
using WaveProp
using WaveProp.BIE
using WaveProp.Geometry
using WaveProp.Integration

@testset "Integral operator test" begin
    pde  = Helmholtz(;dim=3,k=1)
    Ω,M  = WaveProp.IO.gmsh_sphere(dim=2)
    qrule = Gauss{ReferenceTriangle,1}()
    mesh  = NystromMesh(M,qrule)
    𝐒    = SingleLayerOperator(pde,mesh)
    𝐃    = DoubleLayerOperator(pde,mesh)
    @test BIE.kernel_type(𝐒) == BIE.SingleLayer()
    @test BIE.kernel_type(𝐃) == BIE.DoubleLayer()
end