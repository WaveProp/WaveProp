using Test, LinearAlgebra
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.SingularIntegration
using WaveProp.Mesh

@testset "Greens identity test" begin
    # construct interior solution
    pde  = Helmholtz(dim=2,k=1)
    xout = Point(3,3)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    Ω,mesh = WaveProp.IO.gmsh_disk(dim=1,h=0.025,order=2)
    mesh   = GenericMesh{2}(mesh)
    compute_quadrature!(mesh,dim=1,order=3,need_normal=true)
    sum(mesh.qweights)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    𝐒     = SingleLayerOperator(pde,mesh) 
    𝐃     = DoubleLayerOperator(pde,mesh) 
    e0    = WaveProp.Nystrom.error_interior_green_identity(𝐒,𝐃,γ₀u,γ₁u)
    norm(e0,Inf)/norm(γ₀u,Inf)
    qstd  = GaussLegendre(3)
    shand = Kress()
    q     = SingularQuadratureRule(GaussLegendre(8),shand)
    δ𝐒    = singular_weights(𝐒,qstd,q)
    δ𝐃    = singular_weights(𝐃,qstd,q)
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
