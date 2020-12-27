using Test, BenchmarkTools
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

vtx = ((-1,0),(1,0),(0,1))
el  = LagrangeTriangle(vtx)
u   = LagrangeBasis{ReferenceTriangle,1}()
v   = LagrangeBasis{ReferenceTriangle,1}()
q   = Gauss{ReferenceTriangle,3}()

M   = zeros(3,3)
FEM.elementary_matrix!(M,el,u,v,q)
M2  = FEM.mass_matrix(el,u,v,q)
@test M ≈ M2

@btime FEM.elementary_matrix!($M,$el,$u,$v,$q)
@btime FEM.mass_matrix($el,$u,$v,$q)





