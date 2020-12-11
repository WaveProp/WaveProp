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
    mesh = NystromMesh(view(M,Γ);order=1)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    S     = SingleLayerOperator(pde,mesh) 
    D     = DoubleLayerOperator(pde,mesh) 
    e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)
    Sdim  = Nystrom.assemble_dim(S)
    Ddim  = Nystrom.assemble_dim(D)
    e1    = WaveProp.Nystrom.error_interior_green_identity(Sdim,Ddim,γ₀u,γ₁u)
    @test 10*norm(e1,Inf) < norm(e0,Inf)
    S,D   = Nystrom.single_double_layer(pde,mesh)
    e2    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)
    @test e1 ≈ e2
end

@testset "Multiple obstacles" begin
    # construct interior solution
    pde  = Helmholtz(dim=2,k=1)
    xout = Point(-10,0)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    geo1  = Geometry.Kite()
    geo2  = Circle(center=(10,0))
    Ω,M   = meshgen([geo1,geo2],gridsize=0.1)
    Γ     = boundary(Ω)
    # generate a Nystrom mesh with Gauss-Legendre quadrature
    qrule   = GaussLegendre(5)
    e2qrule = Dict(E=>qrule for E in etypes(M))
    mesh    = NystromMesh(M,Γ,e2qrule)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    S     = SingleLayerOperator(pde,mesh) 
    D     = DoubleLayerOperator(pde,mesh) 
    Smk   = Nystrom.assemble_dim(S)
    Dmk   = Nystrom.assemble_dim(D)
    e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)    
    e1 = WaveProp.Nystrom.error_interior_green_identity(Smk,Dmk,γ₀u,γ₁u)
    @test norm(e0,Inf) > 1e-5
    @test norm(e1,Inf) < 1e-5
end