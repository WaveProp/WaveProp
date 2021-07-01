"""
    SingularQuadratureRule{D,Q,S} <: AbstractQuadratureRule{D}

A quadrature rule over `D` intended to integrate functions which are singular at
a known point `s ∈ D`.

A singular quadrature is rule is composed of a *regular* quadrature rule (e.g.
`GaussLegendre`) and a [`AbstractSingularityHandler`](@ref) to transform the regular
quadrature. The regular quadrature rule generates nodes and weights on the
`domain(sing_handler)`, and those are mapped into an appropriate quadrature over
`D = range(sing_handler)` using the singularity handler.
"""
struct SingularQuadratureRule{D,Q,S} <: AbstractQuadratureRule{D}
    qrule::Q
    singularity_handler::S
    function SingularQuadratureRule(q::Q,s::S) where {Q,S}
        @assert domain(s) == domain(q) "domain of quadrature must coincide with domain of singularity handler"
        D = image(s) |> typeof
        new{D,Q,S}(q,s)
    end
end

# getters
qrule(q::SingularQuadratureRule) = q.qrule
singularity_handler(q::SingularQuadratureRule) = q.singularity_handler

domain(qs::SingularQuadratureRule) = domain(singularity_handler(qs))
image(qs::SingularQuadratureRule)  = image(singularity_handler(qs))

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

function (qs::SingularQuadratureRule{ReferenceLine})(s)
    @assert s ∈ image(qs)
    x̂,ŵ = qs()
    # left domain
    x1  = @. s[1] * (1-x̂)
    w1  = @. ŵ*s[1]
    # right domain
    x2  = @. s[1]  + x̂*(1-s[1])
    w2  = @. ŵ*(1-s[1])
    # combine the nodes and weights
    x = vcat(x1,x2)
    w = vcat(w1,w2)
    return x,w
end

"""
    singular_quadrature(k,q::SingularQuadratureRule,s)

Return nodes and weights to integrate a function over `domain(q)` with a
factored weight `k`.
"""
function singular_quadrature(k,q::SingularQuadratureRule,s)
    x,w = q(s)
    T = Base.promote_op(k,eltype(x))
    assert_concrete_type(T)
    w   = map(zip(x,w)) do (x,w)
        k(x)*w
    end
    return x,w
end

"""
    singular_weights(k,xi,q::SingularQuadratureRule,s)
"""
function singular_weights(k,xi,q::SingularQuadratureRule,s)
    x,w = singular_quadrature(k,q,s)
    x   = map(x->x[1],x)
    ws = map(lag_basis) do li
        integrate(x,w) do x
            f(x)*li(x)
        end
    end
    # wlag = barycentric_lagrange_weights(xi)
    # L    = barycentric_lagrange_matrix(xi,x,wlag)
    # return transpose(L)*w
    return ws
end

function singular_weights(k,qreg::AbstractQuadratureRule,qsin::SingularQuadratureRule,s)
    xq,wq = qsin(s)
    xi    = qnodes(qreg)
    lag_basis = lagrange_basis(xi)
    map(lag_basis) do l
        integrate(xq,wq) do x
            k(x)*l(x)
        end
    end
    # return transpose(L)*w
end
