using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Mesh
using WaveProp.Integration

Ω,M = WaveProp.IO.gmsh_disk()
E   = etypes(view(M,boundary(Ω)))[1]
etype2qrule = Dict(E => GaussLegendre(10))
Γ = NystromSpace(M,boundary(Ω),etype2qrule)