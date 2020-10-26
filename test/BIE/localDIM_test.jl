using Test
using WaveProp.Utils
using WaveProp.BIE
using WaveProp.Geometry
using WaveProp.Integration
using QuadGK
using LinearAlgebra
using Plots

pde = Laplace(;dim=2)
n    = 40
h    = 0.1
τ    = line(Point(0,0),Point(h,0))
d⃗    = -20*h/n*normal(τ,0.5)
Ω    = extrude(τ,d⃗)
Γ    = boundary(Ω)

Σ    = Geometry.circle(;radius=10*h,center=Point(h/2,0))


qrule = GaussLegendre{n}()

x̂, ŵ = qrule()
mesh = GenericQuadrature{2,Float64}()
cov = IMT{10,1}()
x̃,w̃ = Integration._push_forward_quad(cov,x̂,ŵ)
for i in 1:4
    if i == 2
        x,w,n⃗ = Integration._push_forward_quad_with_normal(Γ[i],x̃,w̃)
    elseif i == 4
        xr = -x̃ + 1 |> reverse
        wr = reverse(w̃)
        x,w,n⃗ = Integration._push_forward_quad_with_normal(Γ[i],xr,wr)
    else
        x,w,n⃗ = Integration._push_forward_quad_with_normal(Γ[i],x̂,ŵ)    
    end
    append!(mesh.nodes,x)
    append!(mesh.weights,w)
    append!(mesh.normals,n⃗)
end

plot(mesh.nodes,m=:x,label="")

# mesh = GenericQuadrature{2,Float64}(Γ,qrule)

##
# plot(mesh.nodes)
# side = Γ[2]
# x,w = quadgen(side,qrule)
# cov = IMT()
# x̂,ŵ = qrule()
# x,w = push_forward_map(cov,x̂,ŵ)
# sum(w)
# scatter(x,zero(x))
# scatter!(x)
##

# tmp  = Geometry.circle(;radius=h,center=Point(0,0))
# tmp2  = Geometry.circle(;radius=h,center=Point(3*h,0))
# mesh = GenericQuadrature{2,Float64}((tmp,tmp2),qrule)

xs = Point(10,10)

G = SingleLayerKernel(pde)
dG = DoubleLayerKernel(pde)
γ₀Φ = [G(xs,y) for y in qnodes(mesh)]
γ₁Φ = [dG(xs,y,ny) for (y,ny) in zip(qnodes(mesh),qnormals(mesh))]

S = SingleLayerOperator(pde,mesh)
D = DoubleLayerOperator(pde,mesh)

ee = D*γ₀Φ - S*γ₁Φ + γ₀Φ/2
@info norm(ee[1:n],Inf)

xx = [pt[1] for pt in mesh.nodes[1:n]]
plot(xx,ee[1:n],m=:x)

ee[1:n]

# ee = Utils.error_interior_green_identity(S,D,γ₀Φ,γ₁Φ) / norm(γ₀Φ,Inf)  
# norm(ee,Inf)


# @testset "Kernels" begin
#     pde  = Laplace(;dim=2)
#     h    = 0.01
#     k    = SingleLayerKernel(pde)
#     Γ    = ParametricLine() do u
#         Point(h*u,0.)    
#     end    
#     N = 1
#     qrule    = GaussLegendre{2*N}()
#     x̂, ŵ     = qrule()
#     xi,wi,ni    = quadgen_with_normal(Γ,qrule)
#     l        = sum(wi) # approximate length of Γ    
#     d⃗        = 5*l*normal(Γ,0.5)    
#     bnd     = extrude(Γ,d⃗)    
#     els     = bnd[2:end]
#     xq,wq,nq = quadgen_with_normal(els,qrule)
#     xs = xi[N]
#     w  = singular_quadrature(k,Γ,xs,qrule)
#     f  = u-> k(xs,Γ(u))*sqrt(det(transpose(jacobian(Γ,u))*jacobian(Γ,u)))
#     ff = (x) -> k(xs,x)
#     val = quadgk(f,0,1)[1]    
#     val1 = ff.(xi) .* wi |> sum
#     val2 = w  |> sum
#     e1 = val - val1 |> abs
#     e1 = e1 / abs(val)
#     e2 = val - val2 |> abs
#     e2 = e2 / abs(val)
#     @show val, val1, val2
#     @show e1, e2
# end

# xx = [x[1] for x in xi]; yy = [x[2] for x in xi]; scatter(xx,yy,m=:x);p
# xx = [x[1] for x in xq]; yy = [x[2] for x in xq]; scatter!(xx,yy,m=:o);
# xx = [x[1] for x in xs]; yy = [x[2] for x in xs]; scatter!(xx,yy,m=:+)