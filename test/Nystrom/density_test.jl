using Test
using StaticArrays
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.ParametricSurfaces
using WaveProp.Mesh
using LinearAlgebra
using Random
import WaveProp.Nystrom: vals
Random.seed!(1)

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

@testset "TangentialDensity test" begin
    clear_entities!()
    Ω   = ParametricSurfaces.Sphere(;radius=1) |> Domain
    Γ   = boundary(Ω)
    M   = meshgen(Γ,(2,2))
    mesh = NystromMesh(view(M,Γ),order=2)
    pde = Elastostatic(;dim=3,μ=2,λ=3)
    T = Nystrom.default_density_eltype(pde)
    xout = SVector(3,3,3)
    c = rand(T)

    σ = trace(mesh) do target  # density defined with a tangential field
        x = coords(target)
        cross(normal(target), SingleLayerKernel(pde)(xout,target)*c)
    end
    tan_σ = TangentialDensity(σ)
    @test vals(Density(tan_σ)) ≈ vals(σ)

    ncross_σ = ncross(σ)
    ncross_tan_σ = ncross(tan_σ)
    @test vals(TangentialDensity(ncross_σ)) ≈ vals(ncross_tan_σ)
    @test vals(Density(ncross_tan_σ)) ≈ vals(ncross_σ)
end
