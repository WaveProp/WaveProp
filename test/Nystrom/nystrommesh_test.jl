using Test, LinearAlgebra
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Basic tests" begin
    Geometry.clear!()
    Ω,M   = WaveProp.IO.gmsh_sphere(dim=2,h=0.05)
    Γ     = boundary(Ω)
    mesh  = view(M,boundary(Ω))
    nmesh = NystromMesh(mesh,order=2)
    @test isapprox(sum(qweights(nmesh)),π,atol=0.1)
    Geometry.clear!()
    geo   = Circle(radius=0.5)
    Ω     = Domain(geo)
    M     = meshgen(Ω;h=(0.1))
    mesh  = view(M,boundary(Ω))
    nmesh = NystromMesh(mesh,order=2)
    @test isapprox(sum(qweights(nmesh)),π,atol=0.1)
end

@testset "Mesh quadrature" begin
    @testset "Cube" begin
        # generate a mesh
        Geometry.clear!()
        (lx,ly,lz) = widths = (1.,1.,2.)
        Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
        ∂Ω = boundary(Ω)
        mesh  = NystromMesh(view(M,∂Ω),order=1)
        A     = 2*(lx*ly + lx*lz + ly*lz)
        @test A ≈ sum(qweights(mesh))
        # generate a Nystrom mesh for volume. Note that since we are creating a
        # volume Nystrom mesh, we pass the kwarg compute_normal=false since
        # normals do not make sense in this case
        mesh  = NystromMesh(view(M,Ω),order=1,compute_normal=false)
        V     = prod(widths)
        # sum only weights corresponding to tetras
        @test V ≈ sum(qweights(mesh))
        # finally generate a Nystrom mesh for the surface AND volume        
        Ω₁ = union(Ω,∂Ω)
        mesh  = NystromMesh(view(M,Ω₁),order=1,compute_normal=false)
        @test V+A ≈ sum(qweights(mesh))
    end
    @testset "Sphere" begin
        r = 0.5
        Geometry.clear!()
        Ω, M = WaveProp.IO.gmsh_sphere(;radius=r)
        Γ = boundary(Ω)
        A = 4π*r^2
        mesh = NystromMesh(view(M,Γ),order=2)
        # the coarse tolerance below is because we use flat elements to
        # approximate the surface area and volume of a sphere
        @test isapprox(A,sum(qweights(mesh));atol=0.1)        
        # test volume
        mesh = NystromMesh(view(M,Ω),order=2;compute_normal=false)
        A = 4π*r^2
        V = (4/3)*π*r^3
        @test isapprox(V,sum(qweights(mesh));atol=0.1)
    end
    @testset "Disk" begin
        r = rx = ry = 0.5
        Geometry.clear!()
        Ω, M = WaveProp.IO.gmsh_disk(;rx,ry)
        Γ = boundary(Ω)
        mesh = NystromMesh(view(M,Ω),order=2,compute_normal=false)
        A = π*r^2
        # the coarse tolerance below is because we use flat elements to
        # approximate the surface area and volume of a sphere
        @test isapprox(A,sum(qweights(mesh));atol=0.1)        
        # test perimeter
        mesh = NystromMesh(view(M,Γ),order=2)
        P = 2π*r 
        @test isapprox(P,sum(qweights(mesh));atol=0.1)
    end
    @testset "Append" begin
        Geometry.clear!()
        Ω, M = WaveProp.IO.gmsh_disk()
        Γ = boundary(Ω)
        mesh = NystromMesh(view(M,Ω),order=2,compute_normal=false)
        Σ_mesh = NystromMesh{2,Float64}()
        append!(Σ_mesh,mesh)
        append!(Σ_mesh,mesh)
    end    
end
