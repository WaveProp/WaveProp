# # Sphere sound-soft scattering

# In this example we solve the scattering by a sphere. First, we need to load
# the modules
using LinearAlgebra
using WaveProp
using WaveProp.BIE
using WaveProp.Geometry
using WriteVTK

# We now specify which *pde* we want to solve, as well as the geometry:
dim = 3; k=2π
pde = Helmholtz(dim=dim,k=k)
Ω,mesh   = WaveProp.IO.gmsh_sphere(dim=2,h=0.1)
WaveProp.Mesh.compute_quadrature!(mesh,dim=2,order=1,need_normal=true)

# We are interested in solving a sound-soft scattering problem. We will take an
# incident field `uᵢ` as a plane wave, and solve for the scattered field using a
# combined field integral equation:
uᵢ       = (x)   -> exp(im*k*dot((1,0,0),x))
rhs      = γ₀((x)-> -uᵢ(x),mesh)
S,D      = single_double_layer(pde,mesh)
L        = I/2 + D - im*k*S
ϕ        = L\rhs    

# Now that the density ϕ has been found, we can represent the solution 
𝒮        = SingleLayerPotential(pde,mesh)
𝒟        = DoubleLayerPotential(pde,mesh)
u(x)     = uᵢ(x) + 𝒟[ϕ](x) - im*k*𝒮[ϕ](x)

# and output it on a plane for visualization

## output the mesh
fname = joinpath(@__DIR__,"mesh")
vtkfile = WaveProp.IO.vtk_mesh_file(mesh,Ω,fname) |> vtk_save

## compute the solution on a plane
rect,out_mesh  =  WaveProp.IO.gmsh_rectangle(origin=(-5,-5,0),dx=10,dy=10,h=0.1,dim=2)
vals           = [u(x) for x in WaveProp.Mesh.nodes(out_mesh)]
fname = joinpath(@__DIR__,"output")
vtkfile = WaveProp.IO.vtk_mesh_file(out_mesh,rect,fname) 
vtkfile["total field (real part)",VTKPointData()] = real(vals)
vtkfile["total field (imag part)",VTKPointData()] = imag(vals)
vtk_save(vtkfile)






