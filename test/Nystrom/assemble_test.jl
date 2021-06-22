using Test, LinearAlgebra, OrderedCollections
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Single element" begin
    pde  = op = Helmholtz(dim=2,k=1)
    el   = LagrangeLine((1,0),(1.1,0))
    q    = GaussLegendre(20)
    shand = Kress()
    qs  = SingularQuadratureRule(q,shand)
    x̂,ŵ = q()
    yi  = map(el,x̂)
    νi  = map(x->normal(el,x),x̂)
    G   = SingleLayerKernel(op)
    ##
    ϕ  = y->cos(y[1])  # regular density
    s  = 0.51 # location of singularity in parameter space
    x  = el(s)
    ν  = normal(el,s)
    k  = let el = el, G=G
        v -> G(x,el(v))*measure(el,v)
    end
    w  = Nystrom.singular_weights(k,q,qs,s)
    Ie  = integrate(y->G(x,y)*ϕ(y),el)
    Ia  = transpose(w)*ϕ.(yi)
    @test norm(Ie - Ia,Inf)/norm(Ie) < 1e-5
    y₀  = el(s)
    x   = x + 1e-2*ν # outside
    w  = Nystrom.singular_weights(k,q,qs,s)
    Ie  = integrate(y->G(x,y)*ϕ(y),el)
    Ia = transpose(w)*ϕ.(yi)
    @test norm(Ie - Ia,Inf)/norm(Ie) < 1e-5
    x   = x - 1e-2*ν # inside
    w  = Nystrom.singular_weights(k,q,qs,s)
    Ie  = integrate(y->G(x,y)*ϕ(y),el)
    Ia = transpose(w)*ϕ.(yi)
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
    mesh  = NystromMesh(view(M,Γ);order=5)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    S     = SingleLayerOperator(pde,mesh)
    D     = DoubleLayerOperator(pde,mesh)
    e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)
    Sa    = Nystrom.assemble(S)
    Da    = Nystrom.assemble(D)
    e1    = WaveProp.Nystrom.error_interior_green_identity(Sa,Da,γ₀u,γ₁u)
    @test norm(e0,Inf) > 1e-5
    @test norm(e1,Inf) < 1e-5
    # K     = AdjointDoubleLayerOperator(pde,mesh)
    # H     = HyperSingularOperator(pde,mesh)
    # e0    = WaveProp.Nystrom.error_interior_derivative_green_identity(K,H,γ₀u,γ₁u)
    # Kldim  = Nystrom.assemble_ldim(K)
    # Hldim  = Nystrom.assemble_ldim(H)
    # e1    = WaveProp.Nystrom.error_interior_derivative_green_identity(Kldim,Hldim,γ₀u,γ₁u)
    # @test 10*norm(e1,Inf) < norm(e0,Inf)
    # @test norm(e1,Inf) < 1e-4
end

@testset "Single obstacle" begin
    pde  = Helmholtz(dim=2,k=1)
    xout = SVector(-10,0)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    Geometry.clear_entities!()
    geo  = Kite()
    Ω    = Domain(geo)
    M    = meshgen(Ω,n=100)
    Γ    = boundary(Ω)
    # generate a Nystrom mesh with Gauss-Legendre quadrature
    mesh    = NystromMesh(view(M,Γ);order=5)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    S     = SingleLayerOperator(pde,mesh)
    D     = DoubleLayerOperator(pde,mesh)
    e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)
    Sa    = Nystrom.assemble(S)
    Da    = Nystrom.assemble(D)
    e1    = WaveProp.Nystrom.error_interior_green_identity(Sa,Da,γ₀u,γ₁u)
    @test norm(e0,Inf) > 1e-5
    @test norm(e1,Inf) < 1e-5
end
