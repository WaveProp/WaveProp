using Test, LinearAlgebra
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Greens identity test" begin
    # construct interior solution
    Geometry.clear!()
    pde  = Helmholtz(dim=2,k=1)
    xout = SVector(3,3)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    geo = Circle()
    Ω   = Domain(geo)
    M   = meshgen(Ω,h=0.1)
    mesh = NystromMesh(view(M,boundary(Ω));order=5)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    𝐒     = SingleLayerOperator(pde,mesh) 
    𝐃     = DoubleLayerOperator(pde,mesh) 
    e0    = WaveProp.Nystrom.error_interior_green_identity(𝐒,𝐃,γ₀u,γ₁u)
    norm(e0,Inf)/norm(γ₀u,Inf)
    shand = Kress()
    q     = SingularQuadratureRule(GaussLegendre(8),shand)
    δ𝐒    = singular_weights(𝐒,q)
    δ𝐃    = singular_weights(𝐃,q)
    SS    = 𝐒 + δ𝐒
    DD    = 𝐃 + δ𝐃 
    e0    = WaveProp.Nystrom.error_interior_green_identity(SS,DD,γ₀u,γ₁u)
    norm(e0,Inf)/norm(γ₀u,Inf)
    # δS    = singular_weights_dim(𝐒) 
    # δD    = singular_weights_dim(𝐃) 
    # Sfull = Matrix(𝐒) + δS
    # Dfull = Matrix(𝐃) + δD
    # e1 = WaveProp.Nystrom.error_interior_green_identity(Sfull,Dfull,γ₀u,γ₁u)
    # @test norm(e1,Inf) < norm(e0,Inf)
    # S,D = Nystrom.single_double_layer(pde,mesh)
    # e2 = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)
    # @test e1 ≈ e2
end
