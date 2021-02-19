# # Spectrally accurate Nyström method for the Helmholtz equation in 2D
# ## Carlos Pérez Arancibia (cperez@mat.uc.cl)

using WaveProp, WaveProp.Geometry, WaveProp.Integration, WaveProp.Nystrom, WaveProp.Mesh
using LinearAlgebra
using Plots

pyplot()

#=

## Validation

=#
k   = 20
pde = Helmholtz(dim=2,k=k)
# geo = Kite()
δ = 1e-4
ref = ReferenceLine() # the [0,1] segment
bnd = ParametricEntity(ref) do u
     x = cos(2π*u[1])
     y = (δ+x^2)*sin(2π*u[1])
     SVector(x,y)
end     
geo = ParametricBody(boundary=bnd)
G   = SingleLayerKernel(pde)
xₛ  = (0.5,0)
u   = (x) -> G(xₛ,x)
xₒ  = (3,-1)
subfig = plot(geo,label="Γ")
plot!(subfig,[xₛ],m=:o,ms=10,label="source location")
plot!(subfig,[xₒ],m=:o,ms=10,label="test location")
η     = k
Ω,M   = meshgen(geo)
Γ     = boundary(Ω)
niter = 7
er    = zeros(niter)
nn    = 40*collect(1:niter)

for i in 1:niter
     quad_rule = TrapezoidalP(nn[i])
     mesh      = NystromMesh(M,Γ;quad_rule)
     S = SingleLayerOperator(pde,mesh)
     D = DoubleLayerOperator(pde,mesh)
     Smk = Nystrom.assemble_mk(S)
     Dmk = Nystrom.assemble_mk(D)
     L     = I/2 + Dmk - im*η*Smk # create a dense matrix
     rhs = γ₀(u,mesh)
     σ   = L\rhs
     𝒮 = SingleLayerPotential(pde,mesh)
     𝒟 = DoubleLayerPotential(pde,mesh)
     uₙ = (x) -> 𝒟[σ](x) - im*η*𝒮[σ](x)
     xx = (5,5)
     er[i] = abs(u(xx) - uₙ(xx))/abs(u(xx))
end
fig = plot(nn,er,m=:x,yscale=:log10,ylabel="log₁₀(error)",xlabel="n",label="(uₙ - u)/u")

niter = 9
er_dim    = []
er_ldim    = []
nn    = []
quad_rule = GaussLegendre(6)
for i in 1:niter
     Ω,M   = meshgen(geo;gridsize=(1/2)^i)
     Γ     = boundary(Ω)     
     mesh      = NystromMesh(M,Γ;quad_rule)
     S = SingleLayerOperator(pde,mesh)
     D = DoubleLayerOperator(pde,mesh)
     Sdim = Nystrom.assemble_dim(S)
     Ddim = Nystrom.assemble_dim(D)
     Sldim = Nystrom.assemble_ldim(S)
     Dldim = Nystrom.assemble_ldim(D)
     Ldim     = I/2 + Ddim - im*η*Sdim # create a dense matrix
     Lldim     = I/2 + Dldim - im*η*Sldim # create a dense matrix
     rhs = γ₀(u,mesh)
     σdim    = Ldim\rhs
     σldim   = Lldim\rhs
     𝒮 = SingleLayerPotential(pde,mesh)
     𝒟 = DoubleLayerPotential(pde,mesh)
     xx = (5,5)
     uₙ = (x) -> 𝒟[σdim](x) - im*η*𝒮[σdim](x)
     push!(er_dim,abs(u(xx) - uₙ(xx))/abs(u(xx)))
     uₙ = (x) -> 𝒟[σldim](x) - im*η*𝒮[σldim](x)
     push!(er_ldim,abs(u(xx) - uₙ(xx))/abs(u(xx)))
     push!(nn,length(mesh))
end
fig = plot(nn,er_dim,m=:x,xscale=:log2,yscale=:log10,ylabel="log₁₀(error)",xlabel="n",label="gdim")
plot!(fig,nn,er_ldim,m=:x,xscale=:log2,yscale=:log10,ylabel="log₁₀(error)",xlabel="n",label="ldim")


plot(fig,subfig,layout=(2,1))

#===============================================
We consider a kite-shaped curve given by: 
```math
     \Gamma = \{(\cos(t)+0.65\cos(2t)-0.65, 1.5\sin(t))\in\mathbb R^2,t\in (0,2\pi)\}
```
and we let $\Omega\subset\mathbb{R}^2$ be the (open) domain enclosed by
$\Gamma$, i.e., $\Gamma=\partial\Omega$.

We can construct this **parametric curve** and plot it through the `Geometry` module:
===============================================#

geo = Kite()
plot(geo)

#===============================================
The boundary value problem that we consider corresponds to an incident planewave 
$u^{\rm inc}(x_1,x_2) = \exp(ik(x_1\cos\theta+x_2\sin\theta))$ where 
$\theta$ is the angle of incidence and $k>0$ is the wavenumber, that impinges on 
$\Gamma$ giving rise to a scattered field $u$ which satisfies:
```math
\begin{aligned}
\Delta u+k^2 u &=0\quad\text{in}\quad\mathbb R^2\setminus\overline\Omega,\\
u &=-u^{\rm inc}\quad\text{on}\quad\Gamma,\\
\lim_{|x|\to\infty}|x|^{1/2}\left\{\frac{\partial u}{\partial |x|}-iku\right\} &= 0
\end{aligned}
```

We proceed to define the partial differential equation (**PDE**) that we wish so solve,
as well as the incident wavefield
===============================================#
θ   = 0.
d   = (cos(θ),sin(θ))
uᵢ  = (x) -> exp.(im*k*(x⋅d))
∂uᵢ(x,ν) = im*k*(ν⋅d)*uᵢ(x)

#=
In order to solve the exterior boundary value problem, we resort to the combined field integral equation formulation, whereby we seek the scattered field in the form:

$u(x) = (\mathcal D-i\eta S)[\varphi](x),\quad x\in\mathbb R^2\setminus\Gamma,$

where $\mathcal D$ and $\mathcal S$ are the double- and single-layer potentials defined as:

$(\mathcal D\varphi)(x) = \int_{\Gamma}\frac{\partial G(x,y)}{\partial n(y)}\varphi(y)ds(y)\quad\text{and}\quad (\mathcal S\varphi)(x) = \int_{\Gamma}G(x,y)\varphi(y)ds(y).$

Here $G(x,y) = \frac{i}4 H_0^{(1)}(k|x-y|)$ is the free-space Green's function for Helmholtz equation and $\varphi:\Gamma\to\mathbb C$ is an unknown density function that will be determined by solving an integral equation on curve $\Gamma$. 

The constant $\eta>0$ is intruduced so as to ensure uniqueness of solutions of the boundary integral equation. We now define $\eta$ and the tuple _coeff_ that will be used below to construct the combined field potential and operator:

Imposing the Dirichlet boundary condition on $\Gamma$, we arrive at the following second-kind integral equation:

$\left(\frac{I}{2}+K-i\eta S\right)\varphi = -u^{\rm inc}\quad\text{on}\quad\Gamma,$

where operators above are the so-called double- and single-layer integral operators defined as:

$(K\varphi)(x) = \int_{\Gamma}\frac{\partial G(x,y)}{\partial n(y)}\varphi(y)ds(y)\quad\text{and}\quad(S\varphi)(x) = \int_{\Gamma}G(x,y)\varphi(y)ds(y),\qquad x\in\Gamma.$

where we have used the jump relations:

$$\lim_{\epsilon\to 0\pm}(\mathcal D\varphi)(x+\epsilon\nu(x))=\pm\frac{\varphi(x)}{2}+(K\varphi)(x)$$

and 

$$\lim_{\epsilon\to 0\pm}(\mathcal S\psi)(x+\epsilon\nu(x))=(S\psi)(x)$$

for $x\in\Gamma.$

We will approximate the integral operators $S$ and $K$ using a Nyström method, for which we must specify (i) a 
mesh of $\Gamma$ and (ii) a quadrature rule to be used in each element of the mesh. This is done as follows
=#

n     = 10*k
quad_rule = TrapezoidalP(n)
mesh  = NystromMesh(M,Γ;quad_rule)

# Now that an appropriate mesh is available, we can create the integral operators:
S = SingleLayerOperator(pde,mesh)
D = DoubleLayerOperator(pde,mesh)

#=====
The `S` and `D` variables are `IntegralOperator`s, which are in fact a light wrapper around the 
information needed to defined these operators mathematically. They are subtypes of `AbstractMatrix`, 
and behave as such. For instance, we can now *solve* for the unknown density $\sigma$ as follows:
=====#

η   = k
L   = I/2 + D - im*η*S # create a dense matrix
rhs = -γ₀(uᵢ,mesh)
σ   = L\rhs

#=

There are two possible issues with the solution method we just employed: **efficiency** and **accuracy**. We will return to those in a bit, but to get to the end of things, let us now visualize the solution obtained:

=#

𝒮 = SingleLayerPotential(pde,mesh)
𝒟 = DoubleLayerPotential(pde,mesh)
uₛ = (x) -> 𝒟[σ](x) - im*η*𝒮[σ](x)
u  = x -> uᵢ(x) + uₛ(x) # the total field
xgrid = ygrid = -5:(2π/k/10):5
U  = [(x,y) ∈ mesh ? NaN+im*NaN : u((x,y)) for y in ygrid, x in xgrid]
fig1 = heatmap(xgrid,ygrid,real.(U),clims=(-2,2))
plot!(fig1,mesh,lw=6,lc=:white)

#=

Because the integral operators we just used are (weakly) singular, special care is needed to compute the 
entries of `S` and `D` if high-order convergence is to be obtained. There are several techniques for doing 
so, and deciding which to use is highly application dependent. For two-dimensional smooth obstacles admiting a global parametrization, the most accurate method is probably due to Martensen and Kussmaul. To *assemble* the operators `S` and `K` usign this method, we can simply write

=#

Smk = Nystrom.assemble_mk(S)
Dmk = Nystrom.assemble_mk(D)

# The rest follows exactly as before:

L   = I/2 + Dmk - im*η*Smk # create a dense matrix
σ   = L\rhs
poly = copy(qnodes(mesh))
push!(poly,qnodes(mesh)[1])
U    = [ (x,y) in mesh ? NaN+im*NaN : u((x,y)) for y in ygrid, x in xgrid]
fig2 = heatmap(xgrid,ygrid,real.(U),clims=(-2,2))
plot!(fig2,mesh,lw=6,lc=:white)

# The difference between the naive method and the Martensen-Kusmall method is
# evident. In particular the shadow:
plot(fig1,fig2,layout=(1,2),colorbar=:false,ylims=(-5,5))

# Adding more scatterers is simple:
geo   = [Kite(radius=0.1,center=(xc,yc)) for yc in -1:1:1, xc in -1:0.5:1] |> vec
Ω,M   = meshgen(geo)
Γ     = boundary(Ω)
mesh  = NystromMesh(M,Γ;quad_rule)
S     = SingleLayerOperator(pde,mesh)
D     = DoubleLayerOperator(pde,mesh)
Smk   = Nystrom.assemble_mk(S)
Dmk   = Nystrom.assemble_mk(D)
L     = I/2 + Dmk - im*η*Smk # create a dense matrix
rhs   = -γ₀(uᵢ,mesh)
σ     = L\rhs
## now output the solution. 
𝒮     = SingleLayerPotential(pde,mesh)
𝒟     = DoubleLayerPotential(pde,mesh)
uₛ = (x) -> 𝒟[σ](x) - im*η*𝒮[σ](x)
u  = x -> uᵢ(x) .+ uₛ(x) # the total field
xgrid = ygrid = -2:(2π/k/10):2
U  = [(x,y) in mesh ? NaN+im*NaN : u((x,y)) for y in ygrid, x in xgrid]
# U  = [u(SVector(x,y)) for y in ygrid, x in xgrid]
heatmap(xgrid,ygrid,real.(U),clims=(-2,2))
plot!(mesh,lw=2,lc=:black)
plot!(size=(2000,2000))

#=

## Adding your own shapes

Not everybody wants to solve scattering by an arrays of kites, so lets create our own shape and solve scattering problem using that. 
=#
## define the shape throught its boundary
ref = ReferenceLine() # the [0,1] segment
bnd = ParametricEntity(ref) do u
     r = 1 + 0.25*cos(10π*u[1])     
     SVector(r*cos(2π*u[1]),r*sin(2π*u[1]))
end     
geo = ParametricBody(boundary=bnd)
## solve the scattering problem
Ω,M   = meshgen(geo)
Γ     = boundary(Ω)
mesh  = NystromMesh(M,Γ;quad_rule=TrapezoidalP(1000))
S     = SingleLayerOperator(pde,mesh)
D     = DoubleLayerOperator(pde,mesh)
Smk   = Nystrom.assemble_mk(S)
Dmk   = Nystrom.assemble_mk(D)
L     = I/2 + Dmk - im*η*Smk # create a dense matrix
rhs   = -γ₀(uᵢ,mesh)
σ     = L\rhs
## now output the solution. 
𝒮     = SingleLayerPotential(pde,mesh)
𝒟     = DoubleLayerPotential(pde,mesh)
uₛ = (x) -> 𝒟[σ](x) - im*η*𝒮[σ](x)
u  = x -> uᵢ(x) .+ uₛ(x) # the total field
xgrid = ygrid = -3:(2π/k/10):3
U  = [(x,y) ∈ mesh ? NaN+im*NaN : u(SVector(x,y)) for y in ygrid, x in xgrid]
heatmap(xgrid,ygrid,real.(U),clims=(-2,2))
plot!(mesh,lw=4,lc=:white)
plot!(size=(2000,2000))

# # Composite quadrature rules 
k = 40
pde = Helmholtz(k=k,dim=2)
Ω,M   = meshgen(geo;gridsize=0.01)
Γ     = boundary(Ω)
mesh  = NystromMesh(M,Γ;quad_rule=TrapezoidalP(20))
S     = SingleLayerOperator(pde,mesh)
D     = DoubleLayerOperator(pde,mesh)
Smk   = Nystrom.assemble_dim(S)
Dmk   = Nystrom.assemble_dim(D)
L     = I/2 + Dmk - im*k*Smk # create a dense matrix
rhs   = -γ₀(uᵢ,mesh)
σ     = L\rhs
## now output the solution. 
𝒮     = SingleLayerPotential(pde,mesh)
𝒟     = DoubleLayerPotential(pde,mesh)
uₛ = (x) -> 𝒟[σ](x) - im*k*𝒮[σ](x)
u  = x -> uᵢ(x) .+ uₛ(x) # the total field
xgrid = ygrid = -3:(2π/k/10):3
U  = [(x,y) ∈ mesh ? NaN+im*NaN : u(SVector(x,y)) for y in ygrid, x in xgrid]
heatmap(xgrid,ygrid,real.(U),clims=(-2,2))
plot!(mesh,lw=1,lc=:white)
