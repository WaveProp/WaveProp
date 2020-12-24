function complement(el::AbstractElement{ReferenceLine,2},h)
    d      = -normal(el,0.5)
    bnd    = extrude(el,h*d) |> boundary 
    return bnd[2:end]
end

function complement_quadrature(el::AbstractElement{ReferenceLine,2},h,qstd::AbstractQuadratureRule,qsing::SingularQuadratureRule)
    # create the complement of the element
    rside,opside,lside = complement(el,h)
    # compute the singular quadrature (sides)
    x̂,ŵ   = qsing() # singular quadrature on [0,1] with singularity at x=0
    xr    = map(rside,x̂)
    wr    = map(zip(x̂,ŵ)) do (x̂,ŵ)
        μ = measure(rside,x̂)
        μ*prod(ŵ)
    end
    νr    = map(x̂->normal(rside,x̂),x̂)
    # left boundary. Reverse nodes since singularity is on endpoint
    x̂     = reverse(1-x̂)
    ŵ     = reverse(ŵ)
    xl    = map(lside,x̂)
    wl    = map(zip(x̂,ŵ)) do (x̂,ŵ)
        μ = measure(lside,x̂)
        μ*prod(ŵ)
    end
    νl    = map(x̂->normal(lside,x̂),x̂) # normal has to be flipped to point outside
    # now compute regular quadrature on opposite side
    x̂,ŵ   = qstd() # regular quadrature
    xop   = map(opside,x̂)
    wop   = map(zip(x̂,ŵ)) do (x̂,ŵ)
            μ = measure(opside,x̂)
        μ*prod(ŵ)
    end
    νop    = map(x̂->normal(opside,x̂),x̂)
    vcat(xr,xop,xl), vcat(wr,wop,wl), vcat(νr,νop,νl)    
end    

function source_locations(xi,xq)
    num_sources = 3*length(xi)
    bbox  = bounding_box(vcat(xi,xq)) # bounding box for interpolation and quadrature nodes
    N = xi |> first |> length
    if N == 2
        Σ = Geometry.circle(;radius=2*radius(bbox),center=center(bbox))
        ss = [(i-1)/num_sources for i in 1:num_sources]
        xs = map(Σ,ss)
        return xs
    else
        notimplemented()    
    end        
end    

function interpolant_trace(op,(xq,wq,nq),xs)
    G,dG = SingleLayerKernel(op), DoubleLayerKernel(op)
    γ₀Φ  = [G(x,y)    for (y,w) in zip(xq,wq), x in xs]
    γ₁Φ  = [dG(x,y,n) for (y,w,n) in zip(xq,wq,nq), x in xs]    
    return γ₀Φ, γ₁Φ
end  

function interpolation_matrix(op,(xi,νi),xs)
    G,dG = SingleLayerKernel(op), DoubleLayerKernel(op)
    T    = return_type(G)
    M    = Matrix{T}(undef,2*length(xi),length(xs))
    for j in 1:length(xs)
        for i in 1:length(xi)
            M[i,j]            = G(xs[j],xi[i])
            M[length(xi)+i,j] = dG(xs[j],xi[i],νi[i])
        end
    end    
    return M
end    

# function singular_weights(k::AbstractKernel,τ,x,qrule,qinner=qrule)    
#     op       =  k.op
#     # interpolation surface
#     xi,wi,νi = push_forward_quad_with_normal(τ,qrule)
#     # auxiliary integration surface
#     xq,wq,νq = _auxiliary_quadrature(τ,qinner)
#     # find appropriate source points xₛ given the points xq
#     Σ           = _auxiliary_expansion_surface(xq)
#     num_sources = 3*length(xi)
#     ss          = [(i-1)/num_sources for i in 1:num_sources]
#     xs          = [Σ(s) for s in ss]
#     # precompute integration matrices of size `length(xq)×length(xs)`
#     γ₀Φ, γ₁Φ = _auxiliary_integration_matrix(op,(xq,wq,νq),xs)
#     M        = _auxiliary_interpolation_matrix(op,(xi,νi),xs)
#     # compute things which depend on location `x`
#     G,dG     = SingleLayerKernel(op), DoubleLayerKernel(op)
#     # G,dG    = _dim_auxiliary_kernels(k)    
#     𝐊₀       = [G(x,y)*w for (y,w) in zip(xq,wq)]
#     𝐊₁       = [dG(x,y,ny)*w for (y,ny,w) in zip(xq,νq,wq)]
#     tmp      = transpose(𝐊₀)*γ₁Φ - transpose(𝐊₁)*γ₀Φ
#     for i in 1:length(xs)
#         tmp[i] += -G(xs[i],x)/2
#     end
#     α,β = combined_field_coefficients(k)
#     D  = [diagm(α*ones(length(xi))) ; diagm(β*ones(length(xi)))]
#     w  = (tmp/M)*D
#     return w
# end 

# function _auxiliary_quadrature(τ,qrule,cov=Kress(order=5))
#     x̂, ŵ  = qrule()   
#     x̃, w̃  = Integration._push_forward_quad(cov,x̂,ŵ)
#     # build the opposite face quadrature
#     (x,w,ν),(h,t⃗) = _auxiliary_quadrature_opposite_face(τ,x̂,ŵ)
#     # knowing the extrusion parameter d⃗, build the sides quadratures
#     a,b = boundary(τ)
#     xl  = map(x->a+x*(h*t⃗),x̃)
#     xr  = map(x->b+x*(h*t⃗),x̃)
#     wl  = h .* w
#     wr  = h .* w
#     νl  = svector(i->Point(t⃗[2],-t⃗[1]),length(ŵ))
#     νr  = svector(i->Point(-t⃗[2],t⃗[1]),length(ŵ))
#     vcat(xl,xr,x), vcat(wl,wr,w), vcat(νl,νr,ν)
# end    

# function _auxiliary_quadrature_opposite_face(τ,x̂,ŵ)
#     x,w,n = Integration._push_forward_quad_with_normal(τ,x̂,ŵ)
#     l     = sum(w) # length of τ
#     h     = 1*l/length(w)
#     # average normal direction
#     n̄ = mapreduce(+,zip(n,w)) do (n,w)
#         w*n
#     end 
#     d⃗ = -normalize(n̄)
#     # translate nodes, flip normal
#     x = map(x->x+h*d⃗,x)
#     n = -n    
#     return (x,w,n),(h,d⃗)
# end    

# function _auxiliary_expansion_surface(x)
#     bbox  = bounding_box(x)
#     return Geometry.circle(;radius=2*radius(bbox),center=center(bbox))
# end    

# function _auxiliary_integration_matrix(op,(xq,wq,nq),xs)
#     G,dG = SingleLayerKernel(op), DoubleLayerKernel(op)
#     γ₀Φ  = [G(x,y)    for (y,w) in zip(xq,wq), x in xs]
#     γ₁Φ  = [dG(x,y,n) for (y,w,n) in zip(xq,wq,nq), x in xs]    
#     return γ₀Φ, γ₁Φ
# end    

# function _auxiliary_interpolation_matrix(op,(xi,νi),xs)
#     G,dG = SingleLayerKernel(op), DoubleLayerKernel(op)
#     T    = return_type(G)
#     M    = Matrix{T}(undef,2*length(xi),length(xs))
#     for j in 1:length(xs)
#         for i in 1:length(xi)
#             M[i,j]            = G(xs[j],xi[i])
#             M[length(xi)+i,j] = dG(xs[j],xi[i],νi[i])
#         end
#     end    
#     return M
# end    

# function _precompute_interpolation_matrix(op,(xi,ni),xs)
#     G,dG = SingleLayerKernel(op), DoubleLayerKernel(op)
#     T    = return_type(G)
#     M    = Matrix{T}(undef,2*length(xi),length(xs))
#     for i in 1:length(xi)
#         for j in 1:length(xs)
#             M[i,j]            = G(xs[j],xi[i])
#             M[length(xi)+i,j] = dG(xs[j],xi[i],ni[i])
#         end
#     end    
#     return M
# end    
