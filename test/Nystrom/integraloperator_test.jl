using Test, LinearAlgebra
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Basic tests" begin
    pde   = Helmholtz(;dim=3,k=1)
    Ω,M   = WaveProp.IO.gmsh_sphere(dim=2)
    Γ     = boundary(Ω)
    mesh  = NystromMesh(view(M,Γ))
    𝐒     = SingleLayerOperator(pde,mesh)
    𝐃     = DoubleLayerOperator(pde,mesh)
    @test Nystrom.kernel_type(𝐒) == Nystrom.SingleLayer()
    @test Nystrom.kernel_type(𝐃) == Nystrom.DoubleLayer()
end

# this is a rought test that the Greens identities are satisfied. Note that
# since we do not correct for the singular character of the matrix, the errors
# are rather large. This only works, of course, because the kernels are weakly
# singular (integrable), and therefore the silly strategy of ignoring the
# diagonal will converge (albeit slowly).
@testset "Greens identity test" begin
    # construct interior solution
    pde  = Laplace(dim=3)
    xout = Point(3,3,3)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    Ω,M   = WaveProp.IO.gmsh_sphere(dim=2,h=0.1)
    Γ     = boundary(Ω)
    mesh  = NystromMesh(view(M,Γ))
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    𝐒     = SingleLayerOperator(pde,mesh) |> Matrix
    𝐃     = DoubleLayerOperator(pde,mesh) |> Matrix
    ee    = WaveProp.Nystrom.error_interior_green_identity(𝐒,𝐃,γ₀u,γ₁u) / norm(γ₀u,Inf)  
    @test norm(ee,Inf) < 5e-2
end

@testset "Greens identity test in 2d" begin
    # construct interior solution
    pde  = Laplace(dim=2)
    xout = Point(3,3)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    Ω,M   = WaveProp.IO.gmsh_disk(dim=1,h=0.01)
    M    = GenericMesh{2}(M)
    Γ     = boundary(Ω)
    mesh = NystromMesh(view(M,Γ))
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    𝐒     = SingleLayerOperator(pde,mesh) 
    𝐃     = DoubleLayerOperator(pde,mesh) 
    ee = WaveProp.Nystrom.error_interior_green_identity(𝐒,𝐃,γ₀u,γ₁u) / norm(γ₀u,Inf)  
    @test norm(ee,Inf) < 5e-2
end