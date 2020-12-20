using Test, LinearAlgebra
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Single element" begin
    pde  = op = Helmholtz(dim=2,k=1)
    el          = LagrangeLine((1,0),(1.1,0))
    q   = GaussLegendre(10)
    x̂,ŵ = q()
    yi  = map(el,x̂)
    νi  = map(x->normal(el,x),x̂)
    # x   = el(0.51)
    x  = SVector(0.5,-0.1)
    nb = 3*length(yi)
    G   = SingleLayerKernel(op)
    dG  = DoubleLayerKernel(op) 
    w  = Nystrom.singular_weights_ldim(G,el,yi,νi,x,:offsurface)
    ϕ  = y->cos(y[1])    
    I  = integrate(y->G(x,y)*ϕ(y),el)
    Ia = w*ϕ.(yi)
    @test norm(I - Ia,Inf) < 1e-5
end

@testset "Single obstacle" begin
    pde  = Helmholtz(dim=2,k=1)
    xout = SVector(-10,0)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    # geo  = Circle()
    Ω,M   = WaveProp.IO.gmsh_disk()
    M = GenericMesh{2}(M)
    Γ     = boundary(Ω)
    # generate a Nystrom mesh with Gauss-Legendre quadrature
    qrule   = GaussLegendre(10)
    e2qrule = Dict(E=>qrule for E in etypes(M))
    mesh    = NystromMesh(M,Γ,e2qrule)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    S     = SingleLayerOperator(pde,mesh) 
    D     = DoubleLayerOperator(pde,mesh) 
    e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)    
    Sldim = Nystrom.assemble_ldim(S)
    Dldim = Nystrom.assemble_ldim(D)
    e1    = WaveProp.Nystrom.error_interior_green_identity(Sldim,Dldim,γ₀u,γ₁u)
    @test norm(e0,Inf) > 1e-5
    @test norm(e1,Inf) < 1e-5    
end

@testset "Single obstacle" begin
    pde  = Helmholtz(dim=2,k=1)
    xout = SVector(-10,0)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    geo  = Circle()
    Ω,M   = meshgen(geo,gridsize=0.1)
    Γ     = boundary(Ω)
    # generate a Nystrom mesh with Gauss-Legendre quadrature
    qrule   = GaussLegendre(5)
    e2qrule = Dict(E=>qrule for E in etypes(M))
    mesh    = NystromMesh(M,Γ,e2qrule)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    S     = SingleLayerOperator(pde,mesh) 
    D     = DoubleLayerOperator(pde,mesh) 
    e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)    
    Sldim = Nystrom.assemble_ldim(S)
    Dldim = Nystrom.assemble_ldim(D)
    e1    = WaveProp.Nystrom.error_interior_green_identity(Sldim,Dldim,γ₀u,γ₁u)
    @test norm(e0,Inf) > 1e-5
    @test norm(e1,Inf) < 1e-5    
end

