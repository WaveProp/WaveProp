using Test, LinearAlgebra, OrderedCollections
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Single element" begin
    pde  = op = Helmholtz(dim=2,k=1)
    el          = LagrangeLine((1,0),(1.1,0))
    q   = GaussLegendre(20)
    x̂,ŵ = q()
    yi  = map(el,x̂)
    νi  = map(x->normal(el,x),x̂)
    nb = 3*length(yi)
    G   = SingleLayerKernel(op)
    dG  = DoubleLayerKernel(op)
    ##
    ϕ  = y->cos(y[1])  # regular density
    x  = el(0.51)
    ν  = normal(el,0.51)
    w  = Nystrom.singular_weights_ldim(G,el,yi,νi,x,:onsurface)
    Ie  = integrate(y->G(x,y)*ϕ(y),el)
    Ia = w*ϕ.(yi)
    @test norm(Ie - Ia,Inf)/norm(Ie) < 1e-5
    y₀  = el(0.51)
    x   = x + 1e-2*ν # outside
    w  = Nystrom.singular_weights_ldim(G,el,yi,νi,x,:outside)
    Ie  = integrate(y->G(x,y)*ϕ(y),el)
    Ia = w*ϕ.(yi)
    @test norm(Ie - Ia,Inf)/norm(Ie) < 1e-5
    x   = x - 1e-2*ν # inside
    w  = Nystrom.singular_weights_ldim(G,el,yi,νi,x,:inside)
    Ie  = integrate(y->G(x,y)*ϕ(y),el)
    Ia = w*ϕ.(yi)
    @test norm(Ie - Ia,Inf)/norm(Ie) < 1e-5
end

@testset "Single obstacle" begin
    pde  = Helmholtz(dim=2,k=1)
    xout = SVector(-10,0)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    # geo  = Circle()
    Geometry.clear_entities!()
    Ω,M   = WaveProp.IO.gmsh_disk()
    Γ     = boundary(Ω)
    # generate a Nystrom mesh with Gauss-Legendre quadrature
    qrule   = GaussLegendre(3)
    e2qrule = OrderedDict(E=>qrule for E in etypes(M))
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
    K     = AdjointDoubleLayerOperator(pde,mesh)
    H     = HyperSingularOperator(pde,mesh)
    e0    = WaveProp.Nystrom.error_interior_derivative_green_identity(K,H,γ₀u,γ₁u)
    Kldim  = Nystrom.assemble_ldim(K)
    Hldim  = Nystrom.assemble_ldim(H)
    e1    = WaveProp.Nystrom.error_interior_derivative_green_identity(Kldim,Hldim,γ₀u,γ₁u)
    @test 10*norm(e1,Inf) < norm(e0,Inf)
    @test norm(e1,Inf) < 1e-4
end

@testset "Single obstacle" begin
    pde  = Helmholtz(dim=2,k=1)
    xout = SVector(-10,0)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    Geometry.clear_entities!()
    geo  = Kite()
    Ω    = Domain(geo)
    M   = meshgen(Ω,h=0.05)
    Γ     = boundary(Ω)
    # generate a Nystrom mesh with Gauss-Legendre quadrature
    qrule   = GaussLegendre(5)
    e2qrule = OrderedDict(E=>qrule for E in etypes(M))
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
