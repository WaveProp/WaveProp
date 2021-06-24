using Test, LinearAlgebra
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh
using WaveProp.ParametricSurfaces
using WaveProp.Utils
using StaticArrays
using Random
Random.seed!(1)

@testset "DIM" begin
    # construct interior solution
    @testset "Greens identity (interior) 2d" begin
        clear_entities!()
        Ω   = ParametricSurfaces.Circle() |> Domain
        Γ   = boundary(Ω)
        M   = meshgen(Γ,(5,))
        mesh = NystromMesh(view(M,Γ),order=5)
        xout = SVector(3,3)
        ops = (
            Laplace(;dim=2),
            Helmholtz(;dim=2,k=1.2),
            Elastostatic(;dim=2,μ=2,λ=3)
        )
        for pde in ops
            T    = Nystrom.default_density_eltype(pde)
            c    = rand(T)
            u    = (qnode) -> SingleLayerKernel(pde)(xout,qnode)*c
            dudn = (qnode) -> transpose(DoubleLayerKernel(pde)(xout,qnode))*c
            γ₀u   = γ₀(u,mesh)
            γ₁u   = γ₁(dudn,mesh)
            γ₀u_norm = norm(norm.(γ₀u,Inf),Inf)
            γ₁u_norm = norm(norm.(γ₁u,Inf),Inf)
            # single and double layer
            S     = SingleLayerOperator(pde,mesh)
            Smat  = S |> Matrix
            fill_zero_diagonal!(Smat)
            D     = DoubleLayerOperator(pde,mesh)
            Dmat = D |> Matrix
            fill_zero_diagonal!(Dmat)
            e0   = Nystrom.error_interior_green_identity(Smat,Dmat,γ₀u,γ₁u)/γ₀u_norm
            Sdim  = Nystrom.assemble_dim(S)
            Ddim  = Nystrom.assemble_dim(D)
            e1    = Nystrom.error_interior_green_identity(Sdim,Ddim,γ₀u,γ₁u)/γ₀u_norm
            @test norm(e0,Inf) > 1e-2
            @test norm(e1,Inf) < 1e-2
            # adjoint double-layer and hypersingular
            K     = AdjointDoubleLayerOperator(pde,mesh)
            Kmat     = K |> Matrix
            fill_zero_diagonal!(Kmat)
            H     = HyperSingularOperator(pde,mesh)
            Hmat = H |> Matrix
            fill_zero_diagonal!(Hmat)
            e0   = WaveProp.Nystrom.error_interior_derivative_green_identity(Kmat,Hmat,γ₀u,γ₁u)/γ₁u_norm
            Kdim  = Nystrom.assemble_dim(K)
            Hdim  = Nystrom.assemble_dim(H)
            e1   = WaveProp.Nystrom.error_interior_derivative_green_identity(Kdim,Hdim,γ₀u,γ₁u)/γ₁u_norm
            @test norm(e0,Inf) > 1e-2
            @test norm(e1,Inf) < 1e-2
        end
    end

    @testset "Greens identity (interior) 3d" begin
        clear_entities!()
        Ω   = ParametricSurfaces.Sphere() |> Domain
        Γ   = boundary(Ω)
        M   = meshgen(Γ,(2,2))
        mesh = NystromMesh(view(M,Γ),order=5)
        xout = SVector(3,3,3)
        ops = (
            Laplace(;dim=3),
            Helmholtz(;dim=3,k=1.2),
            Elastostatic(;dim=3,μ=2,λ=3),
            Maxwell(;k=1)
        )
        for pde in ops
            T    = Nystrom.default_density_eltype(pde)
            c    = rand(T)
            u    = (qnode) -> SingleLayerKernel(pde)(xout,qnode)*c
            dudn = (qnode) -> transpose(DoubleLayerKernel(pde)(xout,qnode))*c
            γ₀u   = γ₀(u,mesh)
            γ₁u   = γ₁(dudn,mesh)
            γ₀u_norm = norm(norm.(γ₀u,Inf),Inf)
            γ₁u_norm = norm(norm.(γ₁u,Inf),Inf)
            # single and double layer
            S     = SingleLayerOperator(pde,mesh)
            Smat  = S |> Matrix
            fill_zero_diagonal!(Smat)
            D     = DoubleLayerOperator(pde,mesh)
            Dmat = D |> Matrix
            fill_zero_diagonal!(Dmat)
            e0   = Nystrom.error_interior_green_identity(Smat,Dmat,γ₀u,γ₁u)/γ₀u_norm
            Sdim  = Nystrom.assemble_dim(S)
            Ddim  = Nystrom.assemble_dim(D)
            e1    = Nystrom.error_interior_green_identity(Sdim,Ddim,γ₀u,γ₁u)/γ₀u_norm
            @test norm(e0,Inf) > 1e-2
            @test norm(e1,Inf) < 1e-2
            # adjoint double-layer and hypersingular
            pde isa Maxwell && continue
            K     = AdjointDoubleLayerOperator(pde,mesh)
            Kmat     = K |> Matrix
            fill_zero_diagonal!(Kmat)
            H     = HyperSingularOperator(pde,mesh)
            Hmat = H |> Matrix
            fill_zero_diagonal!(Hmat)
            e0   = WaveProp.Nystrom.error_interior_derivative_green_identity(Kmat,Hmat,γ₀u,γ₁u)/γ₁u_norm
            Kdim  = Nystrom.assemble_dim(K)
            Hdim  = Nystrom.assemble_dim(H)
            e1   = WaveProp.Nystrom.error_interior_derivative_green_identity(Kdim,Hdim,γ₀u,γ₁u)/γ₁u_norm
            @test norm(e0,Inf) > 1e-2
            @test norm(e1,Inf) < 1e-2
        end
    end

    @testset "Greens identity (exterior) 2d" begin
        clear_entities!()
        Ω   = ParametricSurfaces.Circle() |> Domain
        Γ   = boundary(Ω)
        M   = meshgen(Γ,(7,))
        mesh = NystromMesh(view(M,Γ),order=5)
        xin = SVector(0.1,0.2)
        ops = (
            Laplace(;dim=2),
            Helmholtz(;dim=2,k=1.2),
            Elastostatic(;dim=2,μ=2,λ=3)
        )
        for pde in ops
            T    = Nystrom.default_density_eltype(pde)
            c    = rand(T)
            u    = (qnode) -> SingleLayerKernel(pde)(xin,qnode)*c
            dudn = (qnode) -> transpose(DoubleLayerKernel(pde)(xin,qnode))*c
            γ₀u   = γ₀(u,mesh)
            γ₁u   = γ₁(dudn,mesh)
            γ₀u_norm = norm(norm.(γ₀u,Inf),Inf)
            γ₁u_norm = norm(norm.(γ₁u,Inf),Inf)
            # single and double layer
            S     = SingleLayerOperator(pde,mesh)
            Smat  = S |> Matrix
            fill_zero_diagonal!(Smat)
            D     = DoubleLayerOperator(pde,mesh)
            Dmat = D |> Matrix
            fill_zero_diagonal!(Dmat)
            e0   = Nystrom.error_exterior_green_identity(Smat,Dmat,γ₀u,γ₁u)/γ₀u_norm
            Sdim  = Nystrom.assemble_dim(S)
            Ddim  = Nystrom.assemble_dim(D)
            e1    = Nystrom.error_exterior_green_identity(Sdim,Ddim,γ₀u,γ₁u)/γ₀u_norm
            @test norm(e0,Inf) > 3e-2
            @test norm(e1,Inf) < 3e-2
            # adjoint double-layer and hypersingular
            K     = AdjointDoubleLayerOperator(pde,mesh)
            Kmat     = K |> Matrix
            fill_zero_diagonal!(Kmat)
            H     = HyperSingularOperator(pde,mesh)
            Hmat = H |> Matrix
            fill_zero_diagonal!(Hmat)
            e0   = WaveProp.Nystrom.error_exterior_derivative_green_identity(Kmat,Hmat,γ₀u,γ₁u)/γ₁u_norm
            Kdim  = Nystrom.assemble_dim(K)
            Hdim  = Nystrom.assemble_dim(H)
            e1   = WaveProp.Nystrom.error_exterior_derivative_green_identity(Kdim,Hdim,γ₀u,γ₁u)/γ₁u_norm
            @test norm(e0,Inf) > 3e-2
            @test norm(e1,Inf) < 3e-2
        end
    end

    @testset "Greens identity (exterior) 3d" begin
        clear_entities!()
        Ω   = ParametricSurfaces.Sphere() |> Domain
        Γ   = boundary(Ω)
        M   = meshgen(Γ,(3, 3))
        mesh = NystromMesh(view(M,Γ),order=5)
        xin = SVector(0.1,-0.2,0)
        ops = (
            Laplace(;dim=3),
            Helmholtz(;dim=3,k=1.2),
            Elastostatic(;dim=3,μ=2,λ=3)
            #Maxwell(;k=1)
        )
        for pde in ops
            T    = Nystrom.default_density_eltype(pde)
            c    = rand(T)
            u    = (qnode) -> SingleLayerKernel(pde)(xin,qnode)*c
            dudn = (qnode) -> transpose(DoubleLayerKernel(pde)(xin,qnode))*c
            γ₀u   = γ₀(u,mesh)
            γ₁u   = γ₁(dudn,mesh)
            γ₀u_norm = norm(norm.(γ₀u,Inf),Inf)
            γ₁u_norm = norm(norm.(γ₁u,Inf),Inf)
            # single and double layer
            S     = SingleLayerOperator(pde,mesh)
            Smat  = S |> Matrix
            fill_zero_diagonal!(Smat)
            D     = DoubleLayerOperator(pde,mesh)
            Dmat = D |> Matrix
            fill_zero_diagonal!(Dmat)
            e0   = Nystrom.error_exterior_green_identity(Smat,Dmat,γ₀u,γ₁u)/γ₀u_norm
            Sdim  = Nystrom.assemble_dim(S)
            Ddim  = Nystrom.assemble_dim(D)
            e1    = Nystrom.error_exterior_green_identity(Sdim,Ddim,γ₀u,γ₁u)/γ₀u_norm
            @test norm(e0,Inf) > 2e-2
            @test norm(e1,Inf) < 2e-2
            # adjoint double-layer and hypersingular
            pde isa Maxwell && continue
            K     = AdjointDoubleLayerOperator(pde,mesh)
            Kmat     = K |> Matrix
            fill_zero_diagonal!(Kmat)
            H     = HyperSingularOperator(pde,mesh)
            Hmat = H |> Matrix
            fill_zero_diagonal!(Hmat)
            e0   = WaveProp.Nystrom.error_exterior_derivative_green_identity(Kmat,Hmat,γ₀u,γ₁u)/γ₁u_norm
            Kdim  = Nystrom.assemble_dim(K)
            Hdim  = Nystrom.assemble_dim(H)
            e1   = WaveProp.Nystrom.error_exterior_derivative_green_identity(Kdim,Hdim,γ₀u,γ₁u)/γ₁u_norm
            @test norm(e0,Inf) > 6e-2
            @test norm(e1,Inf) < 6e-2
        end
    end
end
