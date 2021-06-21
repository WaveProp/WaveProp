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
