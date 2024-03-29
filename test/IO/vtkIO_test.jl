using WaveProp
using WaveProp.Geometry
using WaveProp.IO
using WaveProp.Mesh
using WriteVTK

@testset "VTK export" begin
    # This test should simply not throw an error
    Geometry.clear!()
    Ω, M = WaveProp.IO.gmsh_sphere()
    vtk_mesh_file(M, joinpath(@__DIR__,"ball")) |> vtk_save
    rm(joinpath(@__DIR__,"ball.vtu"))
    vtk_mesh_file(M, Ω, joinpath(@__DIR__,"ball")) |> vtk_save
    rm(joinpath(@__DIR__,"ball.vtu"))
    vtk_mesh_file(M, external_boundary(Ω), joinpath(@__DIR__,"sphere")) |> vtk_save
    rm(joinpath(@__DIR__,"sphere.vtu"))
end
