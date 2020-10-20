using Test
using WaveProp
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.BIE

r    = 0.5
Ω, M = WaveProp.IO.gmsh_sphere(radius=r,dim=2)
qrule = Gauss{ReferenceTriangle,1}()

mesh = NystromMesh(M,qrule)

