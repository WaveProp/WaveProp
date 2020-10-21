using Test
using WaveProp
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Element iterator" begin
    (lx,ly,lz) = widths = (1.,1.,2.)
    Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
    idx  = 2
    E    = etypes(M)[idx]
    iter = ElementIterator{E}(M)
    @test eltype(iter) == E
    i = WaveProp.Mesh._compute_etype_index(iter)
    @test i == idx
    @test length(iter) == size(M.el2vtx[idx],2)
end

@testset "Mesh quadrature" begin
    @testset "Cube" begin
        (lx,ly,lz) = widths = (1.,1.,2.)
        Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
        qrule = Gauss{ReferenceTriangle,1}()
        Q = quadgen(M,qrule;dim=2,need_normal=true)
        A = 2*(lx*ly + lx*lz + ly*lz)
        @test A ≈ sum(Q.weights)
        qrule = Gauss{ReferenceTetrahedron,1}()
        Q = quadgen(M,qrule;dim=3)
        V = prod(widths)
        @test V ≈ sum(Q.weights)
    end
    @testset "Sphere" begin
        r = 0.5
        Ω, M = WaveProp.IO.gmsh_sphere(;radius=r)
        qrule = Gauss{ReferenceTriangle,1}()
        Q = quadgen(M,qrule;dim=2)
        A = 4π*r^2
        # the coarse tolerance below is because we use flat elements to
        # approximate the surface area and volume of a sphere
        @test isapprox(A,sum(Q.weights);atol=0.1)        
        qrule = Gauss{ReferenceTetrahedron,1}()
        Q = quadgen(M,qrule;dim=3)
        V = (4/3)*π*r^3
        @test isapprox(V,sum(Q.weights);atol=0.1)
    end
end

@testset "Nystrom mesh" begin
    @testset "Cube" begin
        (lx,ly,lz) = widths = (1.,1.,2.)
        Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
        qrule = Gauss{ReferenceTriangle,3}()
        mesh  = NystromMesh(M,qrule)
        @test length(mesh) == 3 * length(elements(mesh))
    end
end