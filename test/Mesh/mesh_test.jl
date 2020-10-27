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
    @test length(iter) == size(M.el2nodes[idx],2)
end

@testset "Mesh quadrature" begin
    @testset "Cube" begin
        # generate a mesh (no quadrature by default)
        (lx,ly,lz) = widths = (1.,1.,2.)
        Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
        # generate a quadrature for triangles
        qrule = Gauss{ReferenceTriangle,1}()
        E     = etypes(M)[2] # triangles
        M     = Mesh._compute_quadrature!(M,E,qrule)
        A     = 2*(lx*ly + lx*lz + ly*lz)
        @test A ≈ sum(M.qweights)
        # generate a quadrature for tetras
        E     = etypes(M)[3] # tetras
        qrule = Gauss{ReferenceTetrahedron,1}()
        M     = Mesh._compute_quadrature!(M,E,qrule)
        V     = prod(widths)
        # sum only weights corresponding to tetras
        qtags = M.el2qnodes[3]
        @test V ≈ sum(M.qweights[qtags])
        # sum all weigths (tetras and triangles)
        @test V+A ≈ sum(M.qweights)
    end
    @testset "Cube" begin
        (lx,ly,lz) = widths = (1.,1.,2.)
        Ω, mesh  = WaveProp.IO.gmsh_box(;widths=widths)
        compute_quadrature!(mesh;dim=2,order=1)
        area = 2*(lx*ly + lx*lz + ly*lz)
        @test area ≈ sum(qweights(mesh))
        # volume
        Ω, mesh  = WaveProp.IO.gmsh_box(;widths=widths)
        compute_quadrature!(mesh;dim=3,order=1)
        V = prod(widths)
        @test V ≈ sum(qweights(mesh))
    end
    @testset "Sphere" begin
        r = 0.5
        Ω, M = WaveProp.IO.gmsh_sphere(;radius=r)
        A = 4π*r^2
        compute_quadrature!(M;dim=2,order=1)
        # the coarse tolerance below is because we use flat elements to
        # approximate the surface area and volume of a sphere
        @test isapprox(A,sum(qweights(M));atol=0.1)        
        # test volume
        Ω, M = WaveProp.IO.gmsh_sphere(;radius=r)
        A = 4π*r^2
        compute_quadrature!(M;dim=3,order=1)
        V = (4/3)*π*r^3
        @test isapprox(V,sum(qweights(M));atol=0.1)
    end
    @testset "Disk" begin
        r = rx = ry = 0.5
        Ω, M = WaveProp.IO.gmsh_disk(;rx,ry)
        A = π*r^2
        compute_quadrature!(M;dim=2,order=1)
        # the coarse tolerance below is because we use flat elements to
        # approximate the surface area and volume of a sphere
        @test isapprox(A,sum(qweights(M));atol=0.1)        
        # test perimeter
        Ω, M = WaveProp.IO.gmsh_disk(;rx,ry)
        P = 2π*r
        compute_quadrature!(M;dim=1,order=1)
        @test isapprox(P,sum(qweights(M));atol=0.1)
    end
end
