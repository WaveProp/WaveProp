# Sound-soft scattering

## Problem formulation

In this example we will consider a sound-soft acoustic scattering problem. This
means we will be solving `Helmholtz` equation outside of the region enclosed by
a closed curve curve ``\Gamma``. Letting ``\Omega \subset \mathbb{R}^2`` be the
region enclosed by ``\Gamma``, we must solve

```math
\begin{aligned}
    \Delta u  + k^2 u &= 0, \quad x \in \mathbb{R}^2 \setminus \Omega \\
    u &= 0, \quad x \in \Gamma
\end{aligned}
```
where ``u`` denotes the total field. 

In order to solve the exterior boundary value problem, we resort to the combined
field integral equation formulation, whereby we split the total field ``u`` as an
incident and scattered field, ``u = u_i + u_s``, and seek for the scattered field
in the form:
```math
    u_s(x) = (\mathcal D-i\eta S)[\varphi](x),\quad x\in\mathbb R^2\setminus\Gamma,
```
where $\mathcal D$ and $\mathcal S$ are the double- and single-layer potentials defined as:
```math
\mathcal D[\varphi](x) = \int_{\Gamma}\frac{\partial G(x,y)}{\partial n(y)}\varphi(y)ds(y)\quad\text{and}\quad \mathcal S[\varphi](x) = \int_{\Gamma}G(x,y)\varphi(y)ds(y).
```

Here $G(x,y) = \frac{i}4 H_0^{(1)}(k|x-y|)$ is the free-space Green's function
for Helmholtz equation and $\varphi:\Gamma\to\mathbb C$ is an unknown density
function that will be determined by solving an integral equation on curve
$\Gamma$. The constant $\eta>0$ is introduced so as to ensure uniqueness of
solutions of
the boundary integral equation. We will take ``\eta = k`` as customary. 

Imposing the sound-soft (i.e. homogeneous Dirichlet) boundary condition on $\Gamma$, we arrive at the following second-kind integral equation:
```math
    \left(\frac{I}{2}+D-i\eta S\right)\varphi = -u^{\rm inc}\quad\text{on}\quad\Gamma,
```

where operators above are the so-called double- and single-layer integral
operators defined as:

```math
(D\varphi)(x) = \int_{\Gamma}\frac{\partial G(x,y)}{\partial n(y)}\varphi(y)ds(y)\quad\text{and}\quad(S\varphi)(x) = \int_{\Gamma}G(x,y)\varphi(y)ds(y),\qquad x\in\Gamma.
```
where we have used the jump relations:
```math
\lim_{\epsilon\to 0\pm}(\mathcal D\varphi)(x+\epsilon\nu(x))=\pm\frac{\varphi(x)}{2}+(K\varphi)(x)
```
and
```math
lim_{\epsilon\to 0\pm}(\mathcal S\psi)(x+\epsilon\nu(x))=(S\psi)(x)
```
for $x\in\Gamma.$

The numerical approximation of these operators will give rise to dense matrices,
and the singularity of the Greens function will require the use of special
quadrature rules which are implemented in this package. 

## Scattering by a single obstacle

We now proceed to show how to solve such a problem using `WaveProp`. We will
begin by defining the **PDE**, as well as the incident field ``u_i``:
```@example soundsoft-scattering sound-soft-scattering
    using WaveProp
    using LinearAlgebra
    k   = 30
    pde = WaveProp.PhysicalProblem.Helmholtz(dim=2,k=k)   
    θ   = 0                         # angle of incident wave
    d   = (cos(θ),sin(θ))           # direction of incide wave
    uᵢ  = (x) -> exp.(im*k*(x⋅d))   # incident field
```

Next we need to define a closed curve ``\Gamma``. In our first example we will
take the star-shaped curve given by 
```math
     \Gamma = \{r(u) \cos(2\pi u), r(u)\sin(2\pi u))\in\mathbb R^2,u\in (0,1)\} \quad \text{where} \quad
     r(u) = 1 + 0.25 \cos(10\pi u)
```
To implement this, we create a `ParametricBody` which has boundary given by the
aforementioned curve (represented as a `ParametricEntity`):
```@example soundsoft-scattering sound-soft-scattering
    using WaveProp.Geometry # import various geometry utility functions into namespace
    using Plots
    pyplot() # hide
    ref = ReferenceLine() # the [0,1] segment, i.e. the domain of the parametrization
    bnd = ParametricEntity(ref) do u
        r = 1 + 0.25*cos(10π*u[1])     
        SVector(r*cos(2π*u[1]),r*sin(2π*u[1]))
    end     
    geo = ParametricBody(boundary=bnd)
    plot(geo)
```

!!! tip

    For many structures defined in `WaveProp`, `Plots` recipes are provided so that
    you can always try to just call `plot` on it.

With the domain ``\Omega`` defined through ``\Gamma``, we must now create a mesh
for it. Meshes for `ParametricBody` objects are surface meshes, meaning that
only `\Gamma` will be discretized. In this very simple example where the initial
geometry is known in closed parametric form, the *mesh* will be composed of a
single mesh element: a `ParametricElement` representing the curve itself. We
generate a simple mesh for our surface using:
```@example soundsoft-scattering sound-soft-scattering
    using WaveProp.Mesh # meshing utilities
    Ω,M   = meshgen(geo) 
    bnd   = boundary(Ω)
```
The first argument returned by the `meshgen` function is a `Domain` object which
helps navigate the topology of a mesh.

!!! info
    If this seems overly complicated, it is because more realistic problems
    require to distinguish an entity (of geometrical nature) from the 
    mesh elements comprising it. For instance, the geometric entity ``\Gamma``
    could have been meshed using many mesh segments.

Since we are going to employ a Nystrom method, the last step before actually
assembling the integral operators $S$ and $D$ is to generate a quadrature for
the surface. This quadrature + mesh information is stored in a `NystromMesh`
object, which is the equivalent of appending to a mesh a finite element space.
The `NystromMesh` is created through
```@example soundsoft-scattering sound-soft-scattering
    using WaveProp.Integration # integration utilities
    n         = 10*k # number of dof, proportional to k
    quad_rule = TrapezoidalP(n)
    Γ         = NystromMesh(M,bnd;quad_rule)
```

We are now in position to assemble the integral operators over $\Gamma$ and
solve for $\varphi$:
```@example soundsoft-scattering sound-soft-scattering
    using WaveProp.Nystrom
    S     = SingleLayerOperator(pde,Γ)
    D     = DoubleLayerOperator(pde,Γ)
    Smk   = Nystrom.assemble_mk(S)
    Dmk   = Nystrom.assemble_mk(D)
    L     = I/2 + Dmk - im*k*Smk # create a dense matrix
    rhs   = -γ₀(uᵢ,Γ)
    φ     = L\rhs
```

That is it, we have "solved" the scattering problem! Recalling the combined
field representation, what we have now is a numerical approximation to 
```math
    u(x) = u_i(x) + u_s(x) = u_i(x) + (\mathcal D-i\eta S)[\varphi](x),\quad x\in\mathbb R^2\setminus\Gamma,
```
To plot the solution, we will construct this representation and evaluate at some
points near $\Gamma$:
```@example soundsoft-scattering sound-soft-scattering
    𝒮     = SingleLayerPotential(pde,Γ)
    𝒟     = DoubleLayerPotential(pde,Γ)
    uₛ     = (x) -> 𝒟[φ](x) - im*k*𝒮[φ](x)
    u      = x -> uᵢ(x) .+ uₛ(x) # the total field
    xgrid = ygrid = -3:(2π/k/10):3
    # evalue solution on points outside
    U  = [(x,y) ∈ Γ ? NaN+im*NaN : u(SVector(x,y)) for y in ygrid, x in xgrid]
    heatmap(xgrid,ygrid,real.(U),clims=(-2,2))
    plot!(Γ,lw=4,lc=:white)
```

## Multiple scattering

Now that we have seen how to solve a scattering problem, let us make $\Gamma$
slightly more complicated by adding more obstacles. This can be achieved easily
by creating a vector `ParametricBody` as follows:
```@example soundsoft-scattering sound-soft-scattering
r     = 0.25
geo   = [Circle(radius=r,center=(-0.5,-sqrt(3)/4)), Circle(radius=r,center=(0.5,-sqrt(3)/4)), Circle(radius=r,center=(0,sqrt(3)/4))]
Ω,M   = meshgen(geo)
Γ     = boundary(Ω)
mesh  = NystromMesh(M,Γ;quad_rule)
S     = SingleLayerOperator(pde,mesh)
D     = DoubleLayerOperator(pde,mesh)
Smk   = Nystrom.assemble_mk(S)
Dmk   = Nystrom.assemble_mk(D)
L     = I/2 + Dmk - im*k*Smk # create a dense matrix
rhs   = -γ₀(uᵢ,mesh)
σ     = L\rhs
## now output the solution. 
𝒮     = SingleLayerPotential(pde,mesh)
𝒟     = DoubleLayerPotential(pde,mesh)
uₛ = (x) -> 𝒟[σ](x) - im*k*𝒮[σ](x)
u  = x -> uᵢ(x) .+ uₛ(x) # the total field
xgrid = ygrid = -2:(2π/k/10):2
U  = [(x,y) in mesh ? NaN+im*NaN : u((x,y)) for y in ygrid, x in xgrid]
heatmap(xgrid,ygrid,real.(U),clims=(-2,2))
plot!(mesh,lw=2,lc=:black)
```





