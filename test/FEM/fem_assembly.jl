using Test
using LinearAlgebra
using SparseArrays
using WriteVTK
using WaveProp
using WaveProp.Mesh
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Interpolation
using WaveProp.FEM
using WaveProp.IO

## Assembly test
@testset "Mass matrix 1-2-3D Lagrange P0 P1" begin
    ## 2D
    Ω, mesh3d = WaveProp.IO.gmsh_rectangle(;h=0.1);
    mesh2d = GenericMesh{2}(mesh3d);
    for ord in [1,2,]
        ## matrix assembly
        for (pu, pv) in [(1, 1), (1, 1), (0, 1), (1, 0),]
            u = LagrangeBasis{ReferenceTriangle,pu}()
            v = LagrangeBasis{ReferenceTriangle,pv}()
            @time A = assembly(mesh2d, Ω, u, v; order=ord)
            m, n = size(A)
            surface = ones(1, m) * A * ones(n)
            @test abs(surface[1] - 1) < 1.e-12
            @show m,n
            u = LagrangeBasis{ReferenceLine,pu}()
            v = LagrangeBasis{ReferenceLine,pv}()
            @time A = assembly(mesh2d, boundary(Ω), u, v; order=ord)
            m, n = size(A)
            perimeter = ones(1, m) * A * ones(n)
            @test abs(perimeter[1] - 4) < 1.e-12
            @show m,n
        end
    end
    ## 3D
    Ω, mesh3d = WaveProp.IO.gmsh_box(;h=0.5);
    for ord in [1,2,]
        ## matrix assembly
        for (pu, pv) in [(0, 0), (1, 1), (0, 1), (1, 0),]
            u = LagrangeBasis{ReferenceTetrahedron,pu}()
            v = LagrangeBasis{ReferenceTetrahedron,pv}()
            @time A = assembly(mesh3d, Ω, u, v; order=ord)
            m, n = size(A)
            volume = ones(1, m) * A * ones(n)
            @test abs(volume[1] - 1) < 1.e-12
            @show m,n
            u = LagrangeBasis{ReferenceTriangle,pu}()
            v = LagrangeBasis{ReferenceTriangle,pv}()
            @time A = assembly(mesh3d, boundary(Ω), u, v; order=ord)
            m, n = size(A)
            surface = ones(1, m) * A * ones(n)
            @test abs(surface[1] - 6) < 1.e-12
            @show m,n
            #u = LagrangeBasis{ReferenceLine,pu}()
            #v = LagrangeBasis{ReferenceLine,pv}()
            #@show skeleton(boundary(Ω)) |> tags |> length
            #@time A = assembly(mesh3d, skeleton(boundary(Ω)), u, v; order=ord)
            #m, n = size(A)
            #edgesLength = ones(1, m) * A * ones(n)
            #@test abs(edgesLength[1] - 12) < 1.e-12
            #@show m,n
        end
    end
end


## Helmholtz
k = 1
a_helmholtz(u, v) = (i,j,x̂,el,x) -> (inv(transpose(jacobian(el,x̂))) * grad(u)(x̂)[j]) ⋅ (inv(transpose(jacobian(el,x̂))) *  grad(v)(x̂)[i]) - k^2 * u(x̂)[j] * v(x̂)[i]
@testset "Helmholtz matrix 1-2-3D Lagrange P1" begin
    ## 2D
    Ω, mesh3d = WaveProp.IO.gmsh_rectangle(;h=0.1);
    mesh2d = GenericMesh{2}(mesh3d);
    for ord in [1,2,]
        ## matrix assembly
        for (pu, pv) in [(1, 1),]
            u = LagrangeBasis{ReferenceTriangle,pu}()
            v = LagrangeBasis{ReferenceTriangle,pv}()
            @time A = assembly(mesh2d, Ω, u, v; order=ord, f=a_helmholtz)
            m, n = size(A)
            @show m,n
        end
    end
    ## 3D
    Ω, mesh3d = WaveProp.IO.gmsh_box(;h=0.5);
    for ord in [1,2,]
        ## matrix assembly
        for (pu, pv) in [(1, 1),]
            u = LagrangeBasis{ReferenceTetrahedron,pu}()
            v = LagrangeBasis{ReferenceTetrahedron,pv}()
            @time A = assembly(mesh3d, Ω, u, v; order=ord, f=a_helmholtz)
            m, n = size(A)
            @show m,n
        end
    end
end


## Poisson test
p, q = 2, 3
sol(x,y) = sin(p*π*x) * sin(q*π*y)
f(x,y) = (p*π)^2 * sin(p*π*x) * sin(q*π*y) + (q*π)^2 * sin(p*π*x) * sin(q*π*y)
a_poisson(u, v) = (i,j,x̂,el,x) -> (inv(transpose(jacobian(el,x̂))) * grad(u)(x̂)[j]) ⋅ (inv(transpose(jacobian(el,x̂))) *  grad(v)(x̂)[i])
l_poisson(v) = (i,x̂,el,x) -> f(x[1], x[2]) * v(x̂)[i]
Ω, mesh3d = WaveProp.IO.gmsh_rectangle(;h=0.02);
mesh2d = GenericMesh{2}(mesh3d);
u = LagrangeBasis{ReferenceTriangle,1}()
v = LagrangeBasis{ReferenceTriangle,1}()
A = assembly(mesh2d, Ω, u, v; order=2, f=a_poisson);
b = assembly(mesh2d, Ω, v; order=2, f=l_poisson);
# Dirichlet BC
vΩdofnb = local_dof_numbering(mesh2d, Ω, v);
vΓdofnb = local_dof_numbering(mesh2d, boundary(Ω), trace(v));
glob2loc = SparseVector(length(nodes(mesh2d)),dofs(vΩdofnb),collect(1:length(dofs(vΩdofnb))));
ϵ = 1e-16
for dof in dofs(vΓdofnb)
    i = glob2loc[dof]
    A[i,i] = 1 / ϵ
    b[i] = 0
end
u = factorize(A) \ b;
vtkfile = vtk_mesh_file(mesh3d, Ω, joinpath(@__DIR__,"poisson"))
vtkfile["u", VTKPointData()] = u;
vtkfile |> vtk_save
# Test
uref = [sol(x...) for x in nodes(mesh2d)]
norm(u .- uref)
@test norm(u .- uref) < 4e-3 # for h = 0.02


## Modes for 2D Neumann Laplacian
ak(u, v) = (i,j,x̂,el,x) -> (inv(transpose(jacobian(el,x̂))) * grad(u)(x̂)[j]) ⋅ (inv(transpose(jacobian(el,x̂))) *  grad(v)(x̂)[i])
am(u, v) = (i,j,x̂,el,x) -> u(x̂)[j] * v(x̂)[i]
Ω, mesh3d = WaveProp.IO.gmsh_disk(;h=0.01);
mesh2d = GenericMesh{2}(mesh3d);
u = LagrangeBasis{ReferenceTriangle,1}()
v = LagrangeBasis{ReferenceTriangle,1}()
K = assembly(mesh2d, Ω, u, v; order=2, f=ak);
M = assembly(mesh2d, Ω, u, v; order=2, f=am);
F = eigen(Matrix(K), Matrix(M))
@show F.values |> length
for km in 1:min(100,length(F.values))
    vtkfile = vtk_mesh_file(mesh3d, Ω, joinpath(@__DIR__,"modes_$(km)"))
    vtkfile["u", VTKPointData()] = F.vectors[:,km];
    vtkfile |> vtk_save
end