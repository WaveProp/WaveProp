using Test
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh
using WaveProp.ParametricSurfaces
using StaticArrays


@testset "Interior representation" begin
    # test the interior representation formula u(r) = 𝒮[γ₁u](r) - 𝒟[γ₀u](r)
    @testset "2d" begin
        clear_entities!()
        Ω  = ParametricSurfaces.Circle() |> Domain
        Γ    = boundary(Ω)
        M    = meshgen(Γ,(10,))
        mesh = NystromMesh(view(M,Γ),order=5)
        x₀   = SVector(0.1,-0.2)
        xout = SVector(3,3)
        ops = (
            Laplace(;dim=2),
            Helmholtz(;dim=2,k=1.2),
            Elastostatic(;dim=2,μ=2,λ=3)
        )
        for pde in ops
            S    = SingleLayerPotential(pde,mesh)
            D    = DoubleLayerPotential(pde,mesh)
            T    = Nystrom.default_density_eltype(pde)
            c    = T<:Number ? one(T) : ones(T)
            u    = (qnode) -> SingleLayerKernel(pde)(xout,qnode)*c
            dudn = (qnode) -> transpose(DoubleLayerKernel(pde)(xout,qnode))*c
            γ₀u   = γ₀(u,mesh)
            γ₁u   = γ₁(dudn,mesh)
            uₐ(x) = S[γ₁u](x) - D[γ₀u](x)
            @test isapprox(u(x₀),uₐ(x₀),rtol=1e-3)
        end
    end
    @testset "3d" begin
        clear_entities!()
        Ω  = ParametricSurfaces.Sphere() |> Domain
        Γ    = boundary(Ω)
        M    = meshgen(Γ,(4,4))
        mesh = NystromMesh(view(M,Γ),order=5)
        x₀   = SVector(0.1,-0.2,0.1)
        xout = SVector(3,3,3)
        ops = (
            Laplace(;dim=3),
            Helmholtz(;dim=3,k=1.2),
            Elastostatic(;dim=3,μ=2,λ=3),
            Maxwell(;k=1)
        )
        for pde in ops
            S    = SingleLayerPotential(pde,mesh)
            D    = DoubleLayerPotential(pde,mesh)
            T    = Nystrom.default_density_eltype(pde)
            c    = T<:Number ? one(T) : ones(T)
            u    = (qnode) -> SingleLayerKernel(pde)(xout,qnode)*c
            dudn = (qnode) -> transpose(DoubleLayerKernel(pde)(xout,qnode))*c
            γ₀u   = γ₀(u,mesh)
            γ₁u   = γ₁(dudn,mesh)
            uₐ(x) = S[γ₁u](x) - D[γ₀u](x)
            @test isapprox(u(x₀),uₐ(x₀),rtol=1e-3)
        end
    end
end
