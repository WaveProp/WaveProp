using Test, LinearAlgebra, OrderedCollections
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh
using WaveProp.Utils
using WaveProp.Nystrom: assemble_dim

# Formulation to find J = i k n × H given M = γ₁H / ik

pde  = Maxwell(dim=3,k=1.0)

c    = ones(Nystrom.default_density_eltype(pde))

xout = SVector(3,3,3)

# exact interior solution
H    = (x)   -> SingleLayerKernel(pde)(xout,x)*c

# neumann trace = n × ∇ × H on the boundary
dH = (x,n) -> transpose(DoubleLayerKernel(pde)(xout,x,n))*c

Geometry.clear!()
geo    = WaveProp.Geometry.Sphere()
Ω   = Domain(geo)
Γ   = boundary(Ω)
# Ω,M  = WaveProp.IO.gmsh_sphere(dim=2,h=0.1)
# Γ    = boundary(Ω)
M = meshgen(Ω;h=10,n=32)
mesh = NystromMesh(view(M,Γ);order=2)
# mesh = NystromMesh(view(M,Γ);order=2)

# traces on the mesh
γ₀u   = γ₀(H,mesh) # exact dirichlet trace (to be found)
γ₁u   = γ₁(dH,mesh) # exact neumann trace (boundary condition)

S     = SingleLayerOperator(pde,mesh) |> assemble_dim
D     = DoubleLayerOperator(pde,mesh) |> assemble_dim

ee    = γ₀u/2 + D*γ₀u - S*γ₁u
umax = norm(norm.(γ₀u,Inf),Inf)
norm(norm.(ee,Inf),Inf) / umax

# S     = SingleLayerOperator(pde,mesh)
# D     = DoubleLayerOperator(pde,mesh)

# rhs_vec = S*γ₁u
# L_mat   = I/2 + D
# T = eltype(rhs_vec)

# L   = blockmatrix_to_matrix(L_mat)
# rhs = blockvector_to_vector(rhs_vec)

# sol = L \ rhs

# γ₀u_approx = vector_to_blockvector(sol,T)

# ee = γ₀u_approx - γ₀u

# norm(norm.(ee,Inf),Inf) / norm(norm.(γ₀u,Inf),Inf)
