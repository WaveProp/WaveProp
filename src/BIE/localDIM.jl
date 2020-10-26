struct AuxiliaryElement{Op,T,S}
    Ω::T
    Σ::S
end    

function _auxiliary_quadrature(τ,qrule,cov=IMT{10,1}())
    x̂, ŵ  = qrule()
    x̃,w̃   = _push_forward_quad(cov,x̂,ŵ)
    # quadrature for τ and opposite side
    x,w,n = _push_forward_quad_with_normal(τ,x̂,ŵ)
    l = sum(w)
    h = 5*l/length(w)
    d⃗ = mapreduce(+,zip(n,w)) do (n,w)
        w*n
    end 
    d⃗ = normalize(d⃗)
    a,b = boundary(τ)
    x  = map(x->x+d⃗,x) #translate quadrature nodes
    xl = map(x->a+h*x,x̃)
    xr = map(x->b+h*x,x̃)
    wl = h .* w
    wr = h .* w
    vcat(xr,x,xl), vcat(wr,w,wl), vcat(nr,n,nl)
end    




function singular_quadrature(k::AbstractKernel,Γ,x,qrule)    
    op       =  k.op
    # generate the quadrature for Γ and Γᶜ
    x̂,ŵ      = qrule()
    xi,wi,ni = quadgen_with_normal(Γ,qrule)
    l        = sum(wi) # approximate length of Γ    
    d⃗        = -5*l*normal(Γ,0.5)   
    bnd      = extrude(Γ,d⃗) 
    Γᶜ       = bnd[2:end]
    xq,wq,nq = quadgen_with_normal(Γᶜ,qrule)
    # find appropriate source points xₛ given the points xi and xq
    num_sources = 3*length(xi)
    bbox        = bounding_box(vcat(xi,xq))
    xc,r        = center(bbox),radius(bbox)
    xs          = _circle_sources(num_sources,xc,20*r)
    # precompute integration matrices of size `length(xq)×length(xs)`
    γ₀Φ, γ₁Φ = _precompute_integration_matrices(op,(xq,wq,nq),xs)
    # precompute interpolation matrix of size `2*length(xi)×length(xs)`
    M        = _precompute_interpolation_matrix(op,(xi,ni),xs)
    # compute things which depend on location `x`
    G,dG = SingleLayerKernel(op), DoubleLayerKernel(op)
    # G,dG    = _dim_auxiliary_kernels(k)    
    𝐊₀       = [G(x,y)*w for (y,w) in zip(xq,wq)]
    𝐊₁       = [dG(x,y,ny)*w for (y,ny,w) in zip(xq,nq,wq)]
    tmp      = transpose(𝐊₀)*γ₁Φ - transpose(𝐊₁)*γ₀Φ
    @show tmp      
    if x ∈ xi
        for i in 1:length(xs)
            tmp[i] += -G(xs[i],x)/2
        end
    end  
    @show tmp      
    α,β = combined_field_coefficients(k)
    D  = [diagm(α*ones(length(xi))) ; diagm(β*ones(length(xi)))]
    @info cond(M)
    w  = (tmp/M)*D
    return w
end 

function _precompute_integration_matrices(op,(xq,wq,nq),xs)
    G,dG = SingleLayerKernel(op), DoubleLayerKernel(op)
    γ₀Φ    = [G(x,y)    for (y,w) in zip(xq,wq), x in xs]
    γ₁Φ   = [dG(x,y,n) for (y,w,n) in zip(xq,wq,nq), x in xs]    
    return γ₀Φ, γ₁Φ
end    

function _precompute_interpolation_matrix(op,(xi,ni),xs)
    G,dG = SingleLayerKernel(op), DoubleLayerKernel(op)
    T    = return_type(G)
    M    = Matrix{T}(undef,2*length(xi),length(xs))
    for i in 1:length(xi)
        for j in 1:length(xs)
            M[i,j]            = G(xs[j],xi[i])
            M[length(xi)+i,j] = dG(xs[j],xi[i],ni[i])
        end
    end    
    return M
end    

function _compute_source_locations(num_sources,pts)
    N    = length(first(pts))    
    bbox = bounding_box(pts)
    xc,r = center(bbox),radius(bbox)
    if N == 2
        xs = _circle_sources(n,xc,r)
    elseif N == 3
        @notimplemented    
    else
        @notimplemented
    end            
end    

function _circle_sources(n,xc,r)
    dtheta = 2π / n
    theta = dtheta.*collect(0:1:n-1)
    map(theta) do θ
        x = r*cos(θ) + xc[1]
        y = r*sin(θ) + xc[2]
        Point(x,y)
    end    
end    

# function _precompute_dim_matrices(k,Γ,Γᶜ)
#     op   = pde(k)    
#     α, β = combined_field_coefficients(k)
#     G,dG = SingleLayerKernel(op), DoubleLayerKernel(op)
#     xs  = _compute_source_locations(Γ,Γᶜ)
#     𝐀    = [G(x,y) for y in xq, x in xs]
#     𝐁    = [dG(x,y,ny) for (y,ny) in zip(xq,nq), x in xs]
#     𝐂    = _interpolation_matrix((xi,ni),G,dG)
# end    



# function singular_weights(k::AbstractKernel,x,el::AbstractElement,q,on_surface=true)
#     # generate a quadrature on the extruded element    
#     d⃗    = normal(el,0.5)
#     Γ    = extrude(el,d⃗)  
#     quad = quadgen(Γ,q)
#     xq,nq,wq = nodes(quad), normals(quad), weights(quad)
#     # compute the kernels evaluated on the  extended surface
#     G,dG = _get_dim_kernels(k)
#     γ₀G = [G(x,y)*w     for (y,w) in zip(xq,wq)]
#     γ₁G = [dG(x,y,ny)*w for (y,ny,w) in zip(xq,nq,wq)]
#     # compute the basis function 
#     xs  = _compute_source_locations(τ)
#     ns  = length(xs)
#     γ₀Ψ = [G(x,y) for y in xq, x in xs]
#     γ₁Ψ = [dG(x,y,ny) for (y,ny) in zip(xq,nq), x in xs]
#     tmp = transpose(γ₁G)*γ₀Ψ - transpose(γ₀G)*γ₁Ψ # 1×ns
#     if on_surface
#         for i in 1:ns
#             tmp[i] += -G(x,y[i])/2 
#         end
#     end
#     ## get the interpolation matrix
#     P   = interpolation_matrix(k,xi,xs)
#     # compute weights
#     w = tmp * P
# end

# function _green_interpolation_matrix(pde,x,n⃗,d⃗)    
#     _bbox = bounding_box(x)
#     xl,xh = low_corner(_bbox), high_corner(_bbox)
#     bbox = bounding_box(SVector(xl,xh,xl+d⃗,xh+d⃗))    
#     xc,r = center(bbox), radius(bbox)
#     ns   = 3*length(x)
#     xs   = _circle_sources(ns,xc,5*r)
#     M    = _green_interpolation_matrix(pde,x,n⃗,xs)
# end

# function _green_interpolation_matrix(pde,xi,n⃗,xs,T=default_kernel_eltype(pde))
#     G  = SingleLayerKernel{T}(pde)
#     dG = DoubleLayerKernel{T}(pde)
#     ni, ns  = length(xi), length(xs)
#     M       = Matrix{T}(undef,2*ni,ns)
#     for i in 1:ni
#         for j in 1:ns
#             M[i,j]    = G(xs[j],xi[i])
#             M[ni+i,j] = dG(xs[j],xi[i],n⃗[i])
#         end
#     end
#     return M
# end

# struct GreenBasis{Op,N,T}
#     xs::Vector{Point{N,T}}
# end    