abstract type AbstractQuadratureRule{N} end

struct Gauss{N} <: AbstractQuadratureRule{N} end

struct ClenshawCurtis{N} <: AbstractQuadratureRule{N} end

"""
    quadgen(el::AbstractReferenceShape,qrule::AbstractQuadratureRule) -> (x,w)

Generate a quadrature for `el` using `qrule`. Return the quadrature nodes `x` and weights `w`.

For performance reasons, it is best for the `(x,w)` to be computable from the
types of `el` and `qrule` alone.
"""
function quadgen end    

@generated function quadgen(el::ReferenceLine,qrule::Gauss{N}) where {N}
    x,w = gausslegendre(N) # gives integral in [-1,1]. Converted to [0,1] below
    x = svector(i->Point{1,Float64}(0.5*(x[i]+1)),N) 
    w = svector(i->w[i]/2,N)
    return :($x,$w)
end

@generated function quadgen(el::ReferenceTriangle,qrule::Gauss{N}) where {N}
    if N == 1
        x = svector(i->Point(1/3,1/3),1)
        w = svector(i->1/2,1)
    elseif N == 2
        x = SVector(Point(1/6,1/6),
                    Point(2/3,1/6),
                    Point(1/6,2/3))
        w = svector(i->1/6,3)
    else
        @notimplemented
    end            
    return :($x,$w)
end

@generated function quadgen(el::ReferenceTetrahedron,qrule::Gauss{N}) where {N}
    if N == 1
        x = SVector((Point(1/4,1/4,1/4),))
        w = SVector(1/6)
    elseif N == 2
        a = (5-√5)/20
        b = (5+3*√5)/20
        x = SVector(Point(a,a,a),
                    Point(a,a,b),
                    Point(a,b,a),
                    Point(b,a,a)
                )
        w = svector(i->1/24,4)
    else
        @notimplemented
    end            
    return :($x,$w)
end

"""
    quadgen(el::AbstractElement,qrule::AbstractQuadratureRule) -> (x,w)

Generate a quadrature for `el` using `qrule`. Return the quadrature nodes `x` and weights `w`.

The quadrature is computed by generating a quadrature on the reference
element of `el`, followed by the push-forward map. This requires `typeof(el)` to
implement: 

    - el(x)
    - jacobian(el,x)
"""
function quadgen(el::AbstractElement,qrule)
    # generate a quadrature on the reference element    
    F̂   = reference_element(el)
    x̂,ŵ = quadgen(F̂,qrule)
    # modify the quadrature using the push-forward map
    x   = el.(x̂)
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

function quadgen(M::GenericMesh,qrule;dim=ambient_dimension(M))
    N,T = ambient_dimension(M), eltype(M)
    Q = GenericQuadrature{N,T}()
    for (i,etype) in enumerate(etypes(M))
        geometric_dimension(etype) == dim || continue 
        tags     = M.el2vtx[i]  # Np × Nel matrix
        Np, Nel  = size(tags)   # num of pts per element, num. of elements
        for n in 1:Nel
            el_vtx = M.vtx[tags[:,n]] # get the coordinates of nodes in this element
            el  = etype(el_vtx)       # construct the element
            x,w = quadgen(el,qrule)
            append!(Q.nodes,x)
            append!(Q.weights,w)
        end    
    end     
    return Q   
end    