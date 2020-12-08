using Test, LinearAlgebra
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Greens interpolant test" begin
    # construct interior solution
    pde  = Helmholtz(dim=2,k=1)
    xout = Point(3,3)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    geo  = Geometry.Kite()
    Ω,M  = meshgen(geo,gridsize=100)
    Γ    = boundary(Ω)
    # generate a Nystrom mesh with trapezoidal quadrature
    E    = etypes(M)[1] # element type
    qrule   = TrapezoidalP(40)
    e2qrule = Dict(E=>qrule)
    mesh = NystromMesh(view(M,Γ),e2qrule)
    els = mesh.elements[E]
    el = first(els)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    S     = SingleLayerOperator(pde,mesh) 
    D        = DoubleLayerOperator(pde,mesh) 
    Smk,Dmk  = Nystrom.self_interaction(S,mesh)
    norm(Smk - S,Inf)
    e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)    
    e1 = WaveProp.Nystrom.error_interior_green_identity(Smk,Dmk,γ₀u,γ₁u)
    @test norm(e0,Inf) > 1e-5
    @test norm(e1,Inf) < 1e-5
end