"""
    SingularQuadratureRule{D,Q,S} <: AbstractQuadratureRule{D}

A quadrature rule over `D` intended to integrate functions which are singular at a point
`s ∈ D`.

A singular quadrature is rule is composed of a *regular* quadrature rule (e.g.
`GaussLegendre`) and a [`AbstractSingularityHandler`](@ref) to transform the regular
quadrature. The regular quadrature rule generates nodes and weights on the
`domain(sing_handler)`, and those are mapped into an appropriate quadrature over
`D = range(sing_handler)` using the singularity handler. 

Besides the methods `(q::AbstractQuadratureRule)()` and
`(q::AbstractQuadratureRule)(el::AbstractElement)` described in
[`AbstractQuadratureRule`](@ref), singular quadrature objects also support
`(q::SingularQuadratureRule)(el::AbstractElement,s)`, which returns nodes `x`
and weights `w` for integrating a function over `el` with a possible
(integrable) singularity at parametric location `s ∈ D`. In practice this
typically means that the generated singular quadrature will accumulate points
near `s`, and the exact distribution of how the points are accumulated depends
on the `singularity_handler` employed. 
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

@generated function (qs::SingularQuadratureRule{D,Q,S})() where {D,Q,S}
    qstd  = Q()    
    shand = S()
    x̂,ŵ   = qstd() # reference quadrature
    # use the change of variables cov to modify reference quadrature
    x    = map(x->shand(x),x̂)
    w    = map(x̂,ŵ) do x̂,ŵ
        jac = jacobian(shand,x̂)    
        μ   = abs(det(jac))
        μ*prod(ŵ)
    end 
    return :($x,$w)
end

"""
    singular_quadrature(q::SingularQuadratureRule,s)

Use `q` to produce nodes `x` and weigths `w` for integrating a function over
`domain(q)` with a possible singularity at `s ∈ domain(q)`.

The actual implementation depends on the type of `q`.
"""
function singular_quadrature end

function singular_quadrature(q::SingularQuadratureRule{ReferenceLine},s)
    # remember that s is typically a static vector of length 1    
    @assert 0 < s[1] < 1    
    # qstd  = qrule(q)
    # shand = singularity_handler(q)
    # x̂,ŵ   = qstd()
    # x1    = map(x̂) do x
    #     x*(1-s) + s
    # end    
    # x2    = map(x̂) do x
    #     x*(1-s) + s
    # end    
    x̂,ŵ = q() 
    # left domain
    x1  = @. s[1] * (1-x̂)
    w1  = @. ŵ*s[1]
    # right domain
    x2  = @. s[1]  + x̂*(1-s[1])
    w2  = @. ŵ*(1-s[1])
    # split the domain into two
    # l1    = line(s[1],0)
    # l2    = line(s[1],1)
    # # apply the quadrature to each segment
    # x1,w1 = q(l1)
    # x2,w2 = q(l2)
    # combine the nodes and weights
    return vcat(x1,x2), vcat(w1,w2)
end 


"""
    singular_quadrature(k,q::SingularQuadratureRule,s)

Return nodes and weights to integrate a function over `domain(q)` with a
factored weight `k`. 
"""
function singular_quadrature(k,q::SingularQuadratureRule,s)
    x,w = singular_quadrature(q,s)
    T = Base.promote_op(k,eltype(x))
    assert_concrete_type(T)
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

function singular_weights(k,qstd::AbstractQuadratureRule,q::SingularQuadratureRule,s)
    ui  = qnodes(qstd)
    ui  = map(x->x[1],ui)
    x,w = singular_quadrature(k,q,s)
    x   = map(x->x[1],x)
    wlag = barycentric_lagrange_weights(qstd)
    L    = barycentric_lagrange_matrix(ui,x,wlag)
    return transpose(L)*w
end    

@generated function barycentric_lagrange_weights(Q::AbstractQuadratureRule)
    q  = Q()
    xs = qnodes(q) # static vector
    N  = length(xs)
    T  = eltype(xs) |> eltype
    x  = map(x->x[1],xs)
    # x  = Vector{T}(reinterpret(T,xs))
    w  =  barycentric_lagrange_weights(x) 
    ws = SVector{N}(w)
    return :($ws)
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
