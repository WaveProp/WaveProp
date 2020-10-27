using Test, LinearAlgebra
using WaveProp
using WaveProp.BIE
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Basic tests" begin
    pde   = Helmholtz(;dim=3,k=1)
    Ω,mesh   = WaveProp.IO.gmsh_sphere(dim=2)
    compute_quadrature!(mesh;order=1,dim=2,need_normal=true)
    𝐒     = SingleLayerOperator(pde,mesh)
    𝐃     = DoubleLayerOperator(pde,mesh)
    @test BIE.kernel_type(𝐒) == BIE.SingleLayer()
    @test BIE.kernel_type(𝐃) == BIE.DoubleLayer()
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
    Ω,mesh   = WaveProp.IO.gmsh_sphere(dim=2,h=0.1)
    compute_quadrature!(mesh,order=1,dim=2,need_normal=true)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    𝐒     = SingleLayerOperator(pde,mesh) |> Matrix
    𝐃     = DoubleLayerOperator(pde,mesh) |> Matrix
    𝐒[diagind(𝐒)] .= 0
    𝐃[diagind(𝐃)] .= 0
    ee = WaveProp.Utils.error_interior_green_identity(𝐒,𝐃,γ₀u,γ₁u) / norm(γ₀u,Inf)  
    @test norm(ee,Inf) < 5e-2
end

@testset "Greens identity test in 2d" begin
    # construct interior solution
    pde  = Laplace(dim=2)
    xout = Point(3,3)
    u    = (x)   -> SingleLayerKernel(pde)(xout,x)
    dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
    Ω,mesh   = WaveProp.IO.gmsh_disk(dim=1,h=0.01)
    mesh    = GenericMesh{2}(mesh)
    compute_quadrature!(mesh,order=1,dim=1,need_normal=true)
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    𝐒     = SingleLayerOperator(pde,mesh) |> Matrix
    𝐃     = DoubleLayerOperator(pde,mesh) |> Matrix
    ee = WaveProp.Utils.error_interior_green_identity(𝐒,𝐃,γ₀u,γ₁u) / norm(γ₀u,Inf)  
    @test norm(ee,Inf) < 5e-2
end