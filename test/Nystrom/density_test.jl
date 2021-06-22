using Test
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh
using LinearAlgebra

@testset "Density test" begin
    pde  = Helmholtz(;dim=3,k=1)
    Ω,M  = WaveProp.IO.gmsh_sphere(dim=2)
    Γ    = boundary(Ω)
    mesh = NystromMesh(view(M,Γ),order=1)
    σ     = γ₀(target->norm(coords(target)),mesh)
    @test eltype(σ) == Float64
    σ     = γ₀(mesh) do target
        x = coords(target)
        exp(im*2*norm(x))
    end
    @test eltype(σ) == ComplexF64
end
