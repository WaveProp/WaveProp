using Test, LinearAlgebra
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Single scattering" begin
    # construct interior solution
    pde  = Helmholtz(dim=2,k=1)
    xout = SVector(3,3)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    geo  = Geometry.Kite()
    Ω,M  = meshgen(geo)
    Γ    = boundary(Ω)
    # generate a Nystrom mesh with trapezoidal quadrature
    mesh = NystromMesh(M,Γ;quad_rule=TrapezoidalP(100))
    γ₀u  = γ₀(u,mesh)
    γ₁u  = γ₁(dudn,mesh)
    S    = SingleLayerOperator(pde,mesh) 
    D    = DoubleLayerOperator(pde,mesh) 
    Smk  = Nystrom.assemble_mk(S)
    Dmk  = Nystrom.assemble_mk(D)
    e0   = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)    
    e1   = WaveProp.Nystrom.error_interior_green_identity(Smk,Dmk,γ₀u,γ₁u)
    @test norm(e0,Inf) > 1e-5
    @test norm(e1,Inf) < 1e-5
end

@testset "Multiple scattering" begin
    # construct interior solution
    pde  = Helmholtz(dim=2,k=1)
    xout = SVector(-10,0)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    geo1  = Geometry.Kite()
    geo2  = Circle(center=(10,0))
    Ω,M  = meshgen([geo1,geo2],gridsize=100)
    Γ    = boundary(Ω)
    # generate a Nystrom mesh with trapezoidal quadrature
    mesh  = NystromMesh(M,Γ;quad_rule=TrapezoidalP(100))
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    S     = SingleLayerOperator(pde,mesh) 
    D     = DoubleLayerOperator(pde,mesh) 
    Smk   = Nystrom.assemble_mk(S)
    Dmk   = Nystrom.assemble_mk(D)
    e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)    
    e1 = WaveProp.Nystrom.error_interior_green_identity(Smk,Dmk,γ₀u,γ₁u)
    @test norm(e0,Inf) > 1e-5
    @test norm(e1,Inf) < 1e-5
end

@testset "Change of variables" begin
    # construct interior solution
    pde  = Helmholtz(dim=2,k=1)
    xout = SVector(-10,0)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    geo  = ParametricEntity(ReferenceLine()) do u
        SVector(cos(2π*u[1]),sin(2π*u[1]))
    end  
    bdy  = ParametricBody{2,1,Float64}([geo])
    Ω,M  = meshgen(bdy)
    Γ    = boundary(Ω)
    # generate a Nystrom mesh with trapezoidal quadrature
    mesh  = NystromMesh(M,Γ;quad_rule=TrapezoidalP(40))
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    S     = SingleLayerOperator(pde,mesh) 
    D     = DoubleLayerOperator(pde,mesh) 
    Smk   = Nystrom.assemble_mk(S)
    Dmk   = Nystrom.assemble_mk(D)
    e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)    
    e1    = WaveProp.Nystrom.error_interior_green_identity(Smk,Dmk,γ₀u,γ₁u)
    @test norm(e0,Inf) > 1e-5
    @test norm(e1,Inf) < 1e-5
    norm(e0,Inf)
    norm(e1,Inf)
end

@testset "Composite surface" begin
    # construct interior solution
    cov  = KressP()
    qstd = TrapezoidalP(10)
    qrule = SingularQuadratureRule(qstd,cov)
    pde  = Helmholtz(dim=2,k=1)
    xout = SVector(-2,0)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    bdy  = Circle()
    Ω,M  = meshgen(bdy;gridsize=0.05)
    Γ    = boundary(Ω)
    # generate a Nystrom mesh with trapezoidal quadrature
    mesh  = NystromMesh(M,Γ;quad_rule=qrule)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    S     = SingleLayerOperator(pde,mesh) 
    D     = DoubleLayerOperator(pde,mesh) 
    Smk   = Nystrom.assemble_mk(S)
    Dmk   = Nystrom.assemble_mk(D)
    e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)    
    e1    = WaveProp.Nystrom.error_interior_green_identity(Smk,Dmk,γ₀u,γ₁u)
    @test norm(e0,Inf) > 1e-5
    @test norm(e1,Inf) < 1e-5
    norm(e0,Inf)
    norm(e1,Inf)
end
