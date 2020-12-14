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

"""
    (qs::SingularQuadratureRule)()

Return the nodes and weights for integrating over `domain(qs)`
"""
function (qs::SingularQuadratureRule)()
    qstd = qs |> qrule
    x̂,ŵ  = qstd() # reference quadrature
    cov  = qs |> singularity_handler   
    # use the change of variables cov to modify reference quadrature
    x    = map(x->cov(x),x̂)
    w    = map(zip(x̂,ŵ)) do (x̂,ŵ)
        jac = jacobian(cov,x̂)    
        μ   = abs(det(jac))
        μ*prod(ŵ)
    end 
    return x,w
end

"""
    singular_quadrature(q::SingularQuadratureRule,s)

Use `q` to produce nodes `x` and weigths `w` for integrating a function over
`domain(q)` with a possible singularity at `s ∈ domain(q)`.

The actual implementation depends closely on the type of `q`.
"""
function singular_quadrature end

"""
    singular_quadrature(k,q::SingularQuadratureRule,s)

Return nodes and weights to integrate a function over `domain(q)` with a
factored weight `k`. 
"""
function singular_quadrature(k,q::SingularQuadratureRule,s)
    x,w = singular_quadrature(q,s)
    w   = map(zip(x,w)) do (x,w)
        k(x)*w
    end
    return x,w
end    

"""
    singular_weights(k,ui,q::SingularQuadratureRule,s)

Return weights to integrate a function over `domain(q)` with value known only at
nodes `ui` and a factored weight `k`. 
"""
function singular_weights(k,ui,q::SingularQuadratureRule,s)
    x,w = singular_quadrature(k,q,s)
    x   = map(x->x[1],x)
    wlag = barycentric_lagrange_weights(ui)
    L    = barycentric_lagrange_matrix(ui,x,wlag)
    return transpose(L)*w
end    

"""
    singular_weights(q::SingularQuadratureRule,xi,k,s)

Return the weights to integrate `∫k(x)f(x)dx ≈ sum(f.(xi) .* w)` where a
singularity can be efficiently handled at location `s`.
"""
function singular_weights(qsin::SingularQuadratureRule,xi,k,s)
    x,w = qsin(k,s)
    wlag = barycentric_lagrange_weights(xi)
    L   = barycentric_lagrange_matrix(xi,x,wlag)
    return transpose(L)*w
end 

# function singular_weights(qsin::SingularQuadratureRule,xi,K,s,el)
#     k   = (v) -> begin
#         jac = jacobian(el,v)
#         μ   = transpose(jac)*jac |> det |> sqrt
#         K(el(v))*μ
#     end    
#     singular_weights(qsin,xi,k,s)
# end 
