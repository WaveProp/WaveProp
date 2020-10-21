"""
    quadgen(el::AbstractElement,qrule::AbstractQuadratureRule) -> (x,w)

Generate a quadrature for `el` using `qrule`. Return the quadrature nodes `x` and weights `w`.

The quadrature is computed by generating a quadrature on the reference
element of `el`, followed by the push-forward map. This requires:

    - el(x)
    - jacobian(el,x)
"""
function quadgen(el::AbstractElement,qrule)
    # generate a quadrature on the reference element    
    msg = "the domain of `qrule` must coincide with the domain `el`"
    @assert Geometry.domain(el) == domain(qrule) msg
    x̂,ŵ = qrule()
    # modify the quadrature using the push-forward map
    x   = map(el,x̂)
    μ   = map(x̂) do x̂
        jac = jacobian(el,x̂)
        g   = transpose(jac)*jac |> det
        sqrt(g)
    end 
    w   = ŵ .* μ 
    return x,w
end    
# FIXME: the function above is somewhat inneficient when the ambient and
# geometric dimensions of the element are the same. In that case `μ` simplifies
# to the usual `|det(jac)|`. This should be easy to fix by checking e.g. whether
# `jac` is a square matrix. Since these are static arrays there should be no
# runtime overhead compared to the hand-written version
