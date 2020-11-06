"""
    SingularQuadratureRule{D,Q,S} <: AbstractQuadratureRule{D}

A quadrature rule over `D` intended to integrate functions which are singular at a point
`s ∈ D`.

A singular quadrature is rule is composed of a *regular* quadrature rule (e.g.
`GaussLegendre`) and a [`SingularityHandler`](@ref) to transform the regular
quadrature. The regular quadrature rule generates nodes and weights on the
`domain(sing_handler)`, and those are mapped into an appropriate quadrature over
the `D = range(sing_handler)` using the singularity handler. 

Besides the methods `(q::AbstractQuadratureRule)()` and `(q::AbstractQuadratureRule)(el::AbstractElement)` described in
[`AbstractQuadratureRule`](@ref), 
singular quadrature objects also support
`(q::SingularQuadratureRule)(el::AbstractElement,s)`, which returns nodes `x` and
weights `w` for integrating a function over `el` with a possible (integrable)
singularity at parametric location `s ∈ D`. 
"""
struct SingularQuadratureRule{D,Q,S} <: AbstractQuadratureRule{D}
    qrule::Q
    singularity_handler::S
    function SingularQuadratureRule(q::Q,s::S) where {Q,S}
        @assert domain(s) == domain(q) "domain of quadrature must coincide with domain of singularity handler"
        D = range(s) |> typeof
        new{D,Q,S}(q,s)
    end    
end 

# getters
qrule(q::SingularQuadratureRule) = q.qrule
singularity_handler(q::SingularQuadratureRule) = q.singularity_handler

function (qs::SingularQuadratureRule)()
    qstd = qs |> qrule
    cov  = qs |> singularity_handler    
    x,w  = qstd(cov)
    return x,w
end

function (q::SingularQuadratureRule)(el::AbstractElement,s)
    x̂,ŵ   = singular_quadrature(q,s)
    x   = map(x->el(x),x̂)
    w   = map(zip(x̂,ŵ)) do (x̂,ŵ)
        μ = measure(el,x̂)
        μ*prod(ŵ)
    end 
    return x,w
end    

"""
    integrate(f,q::SingularQuadratureRule,el,s)

Integrate `f` over `el` assuming a singularity of `f` at `el(s)`.
"""
function integrate(f,q::SingularQuadratureRule,el,s)
    x,w = q(el,s)
    integrate(f,x,w)
end    

"""
    singular_quadrature(q::SingularQuadratureRule,s)

Return the nodes and weights to integrate a function over `domain(q)`. The
function can be (weakly) singular at the location `s`.
"""
function singular_quadrature(q::SingularQuadratureRule{ReferenceLine},s)
    @assert 0 < s < 1    
    # split the domain into two
    l1    = line(s,0)
    l2    = line(s,1)
    # apply the quadrature to each segment
    x1,w1 = q(l1)
    x2,w2 = q(l2)
    # combine the nodes and weights
    return vcat(x1,x2), vcat(w1,w2)
end 

function singular_quadrature(q::SingularQuadratureRule{ReferenceTriangle,<:Any,Duffy},s)
    # split the domain
    t1    = triangle((0,0),s,(1,0))
    t2    = triangle((0,0),s,(0,1))
    t3    = triangle((1,0),s,(0,1))
    # apply the quadrature
    x1,w1 = q(t1)
    x2,w2 = q(t2)
    x3,w3 = q(t3)
    return vcat(x1,x2,x3), vcat(w1,w2,w3)
end 

function singular_quadrature(qrule::SingularQuadratureRule{ReferenceSquare,<:Any,Duffy},s)
    x,w   = qrule()    
    # split the domain
    t1    = triangle((0,0),s,(1,0))
    t2    = triangle((0,0),s,(0,1))
    t3    = triangle((1,0),s,(1,1))
    t4    = triangle((0,1),s,(1,1))
    # apply the quadrature
    x1,w1 = qrule(t1)
    x2,w2 = qrule(t2)
    x3,w3 = qrule(t3)
    x4,w4 = qrule(t4)
    return vcat(x1,x2,x3,x4), vcat(w1,w2,w3,w4)
end 

# """
#     singular_weights(q::SingularQuadratureRule,xi,k,s)

# Return the weights to integrate `∫k(x)f(x)dx ≈ sum(f.(xi) .* w)` where a
# singularity can be efficiently handled at location `s`.
# """
# function singular_weights(qsin::SingularQuadratureRule,xi,k,s)
#     x,w = qsin(k,s)
#     wlag = barycentric_lagrange_weights(xi)
#     L   = barycentric_lagrange_matrix(xi,x,wlag)
#     return transpose(L)*w
# end 

# function singular_weights(qsin::SingularQuadratureRule,xi,K,s,el)
#     k   = (v) -> begin
#         jac = jacobian(el,v)
#         μ   = transpose(jac)*jac |> det |> sqrt
#         K(el(v))*μ
#     end    
#     singular_weights(qsin,xi,k,s)
# end 

# function (qrule::SingularQuadratureRule{<:Any,<:TensorProductSingularityHandler})(::ReferenceSquare,s)
#     sx,sy = s[1],s[2]
#     # split the domain
#     s1    = rectangle(s,(0,sy),(0,0),(sx,0))
#     s2    = rectangle(s,(sx,0),(1,0),(1,sy))
#     s3    = rectangle(s,(1,sy),(1,1),(sx,1))
#     s4    = rectangle(s,(sx,1),(0,1),(0,sy))
#     # apply the quadrature
#     x1,w1 = qrule(s1)
#     x2,w2 = qrule(s2)
#     x3,w3 = qrule(s3)
#     x4,w4 = qrule(s4)
#     return vcat(x1,x2,x3,x4), vcat(w1,w2,w3,w4)
# end 


