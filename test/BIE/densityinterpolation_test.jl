using Test, LinearAlgebra
using WaveProp
using WaveProp.BIE
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Greens interpolant test" begin
    # construct interior solution
    pde  = Laplace(dim=3)
    xout = Point(3,3,3)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    Ω,M   = WaveProp.IO.gmsh_sphere(dim=2,h=0.1)
    qrule = Gauss{ReferenceTriangle,1}()
    mesh  = NystromMesh(M,qrule)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    𝐒     = SingleLayerOperator(pde,mesh)
    𝐃     = DoubleLayerOperator(pde,mesh)
    e0    = WaveProp.Utils.error_interior_green_identity(𝐒,𝐃,γ₀u,γ₁u)
    δ𝐒    = GreensCorrection(𝐒) |> Matrix
    δ𝐃    = GreensCorrection(𝐃) |> Matrix
    e1 = WaveProp.Utils.error_interior_green_identity(𝐒+δ𝐒,𝐃+δ𝐃,γ₀u,γ₁u)
    @test norm(e1,Inf) < norm(e0,Inf)
end