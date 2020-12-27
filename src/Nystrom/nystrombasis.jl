"""
    struct NystromBasis{D,Q}
    
Polynomial basis for an reference elemennt of type `D` on the quadrature nodes
of an `AbstractQuadratureRule` of type `Q`. 

This is useful when the value of a function defined on the quadrature nodes is
needed on points on other points for e.g. integrating singular integrals. 
"""
struct NystromBasis{D,Q} <: PolynomialBasis
end    