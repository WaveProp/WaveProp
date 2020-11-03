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
    Ω,mesh   = WaveProp.IO.gmsh_sphere(dim=2,h=0.2)
    compute_quadrature!(mesh,dim=2,order=1,need_normal=true)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    𝐒     = SingleLayerOperator(pde,mesh) 
    𝐃     = DoubleLayerOperator(pde,mesh) 
    e0    = WaveProp.Utils.error_interior_green_identity(𝐒,𝐃,γ₀u,γ₁u)
    δS    = GreensCorrection(𝐒) 
    δD    = GreensCorrection(𝐃) 
    Sfull = Matrix(𝐒) + δS
    Dfull = Matrix(𝐃) + δD
    e1 = WaveProp.Utils.error_interior_green_identity(Sfull,Dfull,γ₀u,γ₁u)
    @test norm(e1,Inf) < norm(e0,Inf)
end
