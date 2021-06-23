using Test, LinearAlgebra, OrderedCollections
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh
using StaticArrays

@testset "Greens interpolant test" begin
    # construct interior solution
    @testset "Helmholtz" begin
        pde  = Helmholtz(dim=3,k=1)
        xout = SVector(3,3,3)
        u    = (qnode) -> SingleLayerKernel(pde)(xout,qnode)
        dudn = (qnode) -> DoubleLayerKernel(pde)(xout,qnode)
        clear_entities!()
        Ω,M  = WaveProp.IO.gmsh_sphere(dim=2,h=0.4)
        Γ    = boundary(Ω)
        mesh = NystromMesh(M,Γ;order=4)
        γ₀u   = γ₀(u,mesh)
        γ₁u   = γ₁(dudn,mesh)
        S     = SingleLayerOperator(pde,mesh)
        Smat  = S |> Matrix
        Smat[diagind(Smat)] .= 0
        D     = DoubleLayerOperator(pde,mesh)
        Dmat = D |> Matrix
        Dmat[diagind(Dmat)] .= 0
        e0    = Nystrom.error_interior_green_identity(Smat,Dmat,γ₀u,γ₁u) / norm(γ₀u,Inf)
        Sdim  = Nystrom.assemble_dim(S)
        Ddim  = Nystrom.assemble_dim(D)
        e1    = Nystrom.error_interior_green_identity(Sdim,Ddim,γ₀u,γ₁u) /  norm(γ₀u,Inf)
        @test 10*norm(e1,Inf) < norm(e0,Inf)
        S,D   = Nystrom.single_doublelayer_dim(pde,mesh)
        e2    = Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u) / norm(γ₀u,Inf)
        @test e1 ≈ e2
        @test norm(e1,Inf) < 1e-2
        K     = AdjointDoubleLayerOperator(pde,mesh)
        H     = HyperSingularOperator(pde,mesh)
        e0    = Nystrom.error_interior_derivative_green_identity(K,H,γ₀u,γ₁u) / norm(γ₁u,Inf)
        Kdim  = Nystrom.assemble_dim(K)
        Hdim  = Nystrom.assemble_dim(H)
        e1    = Nystrom.error_interior_derivative_green_identity(Kdim,Hdim,γ₀u,γ₁u) / norm(γ₁u,Inf)
        @test 10*norm(e1,Inf) < norm(e0,Inf)
        @test norm(e1,Inf) < 1e-2
    end
    @testset "Elastostatic" begin
        pde  = Elastostatic(dim=3,μ=1,λ=2.)
        xout = SVector(3,3,3)
        c    = ones(Nystrom.default_density_eltype(pde))
        u    = (dof)   -> SingleLayerKernel(pde)(xout,dof)*c
        dudn = (dof) -> transpose(DoubleLayerKernel(pde)(xout,dof))*c
        clear_entities!()
        Ω,M  = WaveProp.IO.gmsh_sphere(dim=2,h=0.4)
        Γ    = boundary(Ω)
        mesh = NystromMesh(view(M,Γ);order=4)
        γ₀u   = γ₀(u,mesh)
        γ₁u   = γ₁(dudn,mesh)
        Sop     = SingleLayerOperator(pde,mesh)
        Dop     = DoubleLayerOperator(pde,mesh)
        e0    = WaveProp.Nystrom.error_interior_green_identity(Sop,Dop,γ₀u,γ₁u) / norm(norm.(γ₀u,Inf),Inf)
        Sdim  = Nystrom.assemble_dim(Sop)
        Ddim  = Nystrom.assemble_dim(Dop)
        e1    = WaveProp.Nystrom.error_interior_green_identity(Sdim,Ddim,γ₀u,γ₁u) / norm(norm.(γ₀u,Inf),Inf)
        @test 10*norm(norm.(e1,Inf),Inf) < norm(norm.(e0,Inf),Inf)
        Sdim2,Ddim2   = Nystrom.single_doublelayer_dim(pde,mesh)
        e2    = WaveProp.Nystrom.error_interior_green_identity(Sdim2,Ddim2,γ₀u,γ₁u) /norm(norm.(γ₀u,Inf),Inf)
        @test e1 ≈ e2
        Kop     = AdjointDoubleLayerOperator(pde,mesh)
        Hop     = HyperSingularOperator(pde,mesh)
        e0    = WaveProp.Nystrom.error_interior_derivative_green_identity(Kop,Hop,γ₀u,γ₁u) / norm(norm.(γ₁u,Inf),Inf)
        Kdim  = Nystrom.assemble_dim(Kop)
        Hdim  = Nystrom.assemble_dim(Hop)
        e1    = WaveProp.Nystrom.error_interior_derivative_green_identity(Kdim,Hdim,γ₀u,γ₁u) / norm(norm.(γ₁u,Inf),Inf)
        @test 10*norm(norm(e1,Inf),Inf) < norm(norm(e0,Inf),Inf)
        @test norm(norm(e1,Inf),Inf) < 1e-2
    end
#     @testset "Maxwell" begin
#         pde  = Maxwell(dim=3,k=1.0)
#         xout = SVector(3,3,3)
#         c    = ones(Nystrom.default_density_eltype(pde))
#         u    = (x)   -> SingleLayerKernel(pde)(xout,x)*c
#         dudn = (x,n) -> transpose(DoubleLayerKernel(pde)(xout,x,n))*c
#         Geometry.clear!()
#         geo    = WaveProp.Geometry.Sphere()
#         Ω   = Domain(geo)
#         Γ   = boundary(Ω)
#         M = meshgen(Ω;h=10,n=16)
#         mesh = NystromMesh(view(M,Γ);order=3)
#         γ₀u   = γ₀(u,mesh)
#         normals = mesh.qnormals
#         γ₁u   = γ₁(dudn,mesh)
#         S     = SingleLayerOperator(pde,mesh)
#         D     = DoubleLayerOperator(pde,mesh)
#         e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)
#         Sdim  = Nystrom.assemble_dim(S)
#         Ddim  = Nystrom.assemble_dim(D)
#         e1    = WaveProp.Nystrom.error_interior_green_identity(Sdim,Ddim,γ₀u,γ₁u)
#         @test 10*norm(norm.(e1,Inf),Inf) < norm(norm.(e0,Inf),Inf)
#         S,D   = Nystrom.single_doublelayer_dim(pde,mesh)
#         e2    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)
#         @test e1 ≈ e2
#         norm(norm.(e1,Inf),Inf) / norm(norm(γ₀u,Inf),Inf)
#     end
end

# @testset "Multiple obstacles" begin
#     # construct interior solution
#     pde  = Helmholtz(dim=2,k=1)
#     xout = SVector(-10,0)
#     u    = (x)   -> SingleLayerKernel(pde)(xout,x)
#     dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
#     Geometry.clear!()
#     geo1  = Geometry.Kite()
#     geo2  = Circle(center=(10,0))
#     Ω     = Domain([geo1,geo2])
#     M     = meshgen(Ω,h=0.1)
#     Γ     = boundary(Ω)
#     # generate a Nystrom mesh with Gauss-Legendre quadrature
#     qrule   = GaussLegendre(5)
#     e2qrule = OrderedDict(E=>qrule for E in etypes(M))
#     mesh    = NystromMesh(M,Γ,e2qrule)
#     γ₀u   = γ₀(u,mesh)
#     γ₁u   = γ₁(dudn,mesh)
#     S     = SingleLayerOperator(pde,mesh)
#     D     = DoubleLayerOperator(pde,mesh)
#     Smk   = Nystrom.assemble_dim(S)
#     Dmk   = Nystrom.assemble_dim(D)
#     e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)
#     e1 = WaveProp.Nystrom.error_interior_green_identity(Smk,Dmk,γ₀u,γ₁u)
#     @test norm(e0,Inf) > 1e-5
#     @test norm(e1,Inf) < 1e-5
#     @test sum(qweights(mesh[Γ[2]])) ≈ 2π
# end


# @testset "Extract domain mesh" begin
#     # construct interior solution
#     pde  = Helmholtz(dim=2,k=1)
#     xout = SVector(-10,0)
#     u    = (x)   -> SingleLayerKernel(pde)(xout,x)
#     dudn = (x,n) -> DoubleLayerKernel(pde)(xout,x,n)
#     Geometry.clear!()
#     geo1  = Geometry.Kite()
#     geo2  = Circle(center=(10,0))
#     Ω     = Domain([geo1,geo2])
#     M   = meshgen(Ω,h=0.1)
#     Γ     = boundary(Ω)
#     # generate a Nystrom mesh with Gauss-Legendre quadrature
#     qrule   = GaussLegendre(5)
#     e2qrule = OrderedDict(E=>qrule for E in etypes(M))
#     Γ_mesh    = NystromMesh(M,Γ,e2qrule)
#     Γ1_mesh   = Γ_mesh[Γ[1]]
#     Γ2_mesh   = Γ_mesh[Γ[2]]
#     γ₀u   = γ₀(u,Γ2_mesh)
#     γ₁u   = γ₁(dudn,Γ2_mesh)
#     S     = SingleLayerOperator(pde,Γ2_mesh)
#     D     = DoubleLayerOperator(pde,Γ2_mesh)
#     Smk   = Nystrom.assemble_dim(S)
#     Dmk   = Nystrom.assemble_dim(D)
#     e0    = WaveProp.Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)
#     e1 = WaveProp.Nystrom.error_interior_green_identity(Smk,Dmk,γ₀u,γ₁u)
#     @test norm(e0,Inf) > 1e-5
#     @test norm(e1,Inf) < 1e-5
# end
