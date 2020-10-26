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
    @assert domain(el) == domain(qrule) msg
    # modify the quadrature using the push-forward map
    x,w = push_forward_quad(el,qrule)
    return x,w
end    
