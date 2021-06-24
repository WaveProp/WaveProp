using Test, LinearAlgebra
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh
using WaveProp.ParametricSurfaces

@testset "Basic tests" begin
    clear_entities!()
    Ω,M   = WaveProp.IO.gmsh_sphere(dim=3,h=0.05,order=2)
    Γ     = boundary(Ω)
    # mesh  = view(M,boundary(Ω))
    nmesh = NystromMesh(M,Γ,order=4)
    area   = sum(qweights(nmesh))
    @test isapprox(area,π,atol=1e-5)
    clear_entities!()
    Ω     = ParametricSurfaces.Sphere() |> Domain
    Γ     = boundary(Ω)
    M     = meshgen(Γ,(4,4))
    nmesh = NystromMesh(view(M,Γ),order=5)
    area   = sum(qweights(nmesh))
    @test isapprox(area,π,atol=1e-5)
    clear_entities!()
    geo   = ParametricSurfaces.Circle(radius=0.5)
    Ω     = Domain(geo)
    Γ     = boundary(Ω)
    M     = meshgen(Γ,(100,))
    nmesh = NystromMesh(M,Γ,order=5)
    @test isapprox(sum(qweights(nmesh)),π,atol=1e-10)
end

@testset "Mesh quadrature" begin
    @testset "Cube" begin
        # generate a mesh
        clear_entities!()
        (lx,ly,lz) = widths = (1.,1.,2.)
        Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
        ∂Ω = boundary(Ω)
        mesh  = NystromMesh(view(M,∂Ω),order=1)
        A     = 2*(lx*ly + lx*lz + ly*lz)
        @test A ≈ sum(qweights(mesh))
        # generate a Nystrom mesh for volume. Note that since we are creating a
        # volume Nystrom mesh, we pass the kwarg compute_normal=false since
        # normals do not make sense in this case
        mesh  = NystromMesh(view(M,Ω),order=1)
        V     = prod(widths)
        # sum only weights corresponding to tetras
        @test V ≈ sum(qweights(mesh))
    end
    @testset "Sphere" begin
        r = 0.5
        Geometry.clear_entities!()
        Ω, M = WaveProp.IO.gmsh_sphere(;radius=r,order=1)
        Γ = boundary(Ω)
        A = 4π*r^2
        mesh = NystromMesh(view(M,Γ),order=2)
        # the coarse tolerance below is because we use flat elements to
        # approximate the surface area and volume of a sphere
        @test isapprox(A,sum(qweights(mesh));atol=1e-2)
        # test volume
        mesh = NystromMesh(view(M,Ω),order=2)
        A = 4π*r^2
        V = (4/3)*π*r^3
        @test isapprox(V,sum(qweights(mesh));atol=1e-2)
    end
    @testset "Disk" begin
        r = rx = ry = 0.5
        Geometry.clear_entities!()
        Ω, M = WaveProp.IO.gmsh_disk(;rx,ry)
        Γ = boundary(Ω)
        mesh = NystromMesh(view(M,Ω),order=2)
        A = π*r^2
        # the coarse tolerance below is because we use flat elements to
        # approximate the surface area and volume of a sphere
        @test isapprox(A,sum(qweights(mesh));atol=1e-2)
        # test perimeter
        mesh = NystromMesh(view(M,Γ),order=2)
        P = 2π*r
        @test isapprox(P,sum(qweights(mesh));atol=1e-2)
    end
end
