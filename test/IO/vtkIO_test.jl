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


@testset "Partitioning" begin
    # This test should simply not throw an error
    Geometry.clear!()
    Ω, M = WaveProp.IO.gmsh_disk(;nΩ=2, h=0.1)
    vtmfile = vtk_multiblock(joinpath(@__DIR__,"disk"))
    for ω in Ω
        Ωi = Domain(ω)
        vtkfile = vtk_mesh_file(M, Ωi, vtmfile)
    end
    vtk_save(vtmfile)
    rm(joinpath(@__DIR__,"disk.vtm"))
    rm(joinpath(@__DIR__,"disk_1.vtu"))
    rm(joinpath(@__DIR__,"disk_2.vtu"))
    rm(joinpath(@__DIR__,"disk_3.vtu"))
end