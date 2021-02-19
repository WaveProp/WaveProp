using Test
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Density test" begin
    pde  = Helmholtz(;dim=3,k=1)
    Ω,M  = WaveProp.IO.gmsh_sphere(dim=2)
    Γ    = boundary(Ω)
    mesh = NystromMesh(view(M,Γ),order=1)
    σ     = γ₀(x->x[1],mesh)
    @test eltype(σ) == Float64
    σ     = γ₀(x->exp(im*x[3]),mesh)
    @test eltype(σ) == ComplexF64
    σ     = γ₁((x,n)->sum(x.*n),mesh)
    @test eltype(σ) == Float64
end