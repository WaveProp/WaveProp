using Test, LinearAlgebra
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Greens interpolant test" begin
    # construct interior solution
    pde  = Helmholtz(dim=3,k=1)
    xout = Point(3,3,3)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    Ω,M  = WaveProp.IO.gmsh_sphere(dim=2,h=0.2)
    Γ    = boundary(Ω)
    mesh = NystromMesh(view(M,Γ))
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    S     = SingleLayerOperator(pde,mesh) 
    D     = DoubleLayerOperator(pde,mesh) 
    e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)
    δS    = singular_weights_dim(S) 
    δD    = singular_weights_dim(D) 
    Sfull = Matrix(S) + δS
    Dfull = Matrix(D) + δD
    e1 = WaveProp.Nystrom.error_interior_green_identity(Sfull,Dfull,γ₀u,γ₁u)
    @test 10*norm(e1,Inf) < norm(e0,Inf)
    S,D = Nystrom.single_double_layer(pde,mesh)
    e2 = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)
    @test e1 ≈ e2
end