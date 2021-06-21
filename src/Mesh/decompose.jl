"""
    decompose(s::AbstractReferenceShape,x,[target_shape])

Decompose an [`AbstractReferenceShape`](@ref) into [`LagrangeElement`](@ref)s so
that `x` is a vertex of the children elements.

# Examples
s = ReferenceLine()
el1, el2 = decompose(s,0.3)
el1(1) == el2(0) == 0.3 # true
"""
function decompose(ln::ReferenceLine,x::Float64)
    @assert x ∈ ln
    a,b = vertices(ln)
    return LagrangeLine(x,a),LagrangeLine(x,b)
end

function decompose(tri::ReferenceTriangle,x)
    @assert x ∈ tri
    a,b,c = vertices(tri)
    return LagrangeTriangle(a,x,b), LagrangeTriangle(b,x,c), LagrangeTriangle(c,x,a)
end

function decompose(sq::ReferenceSquare,x)
    @assert x ∈ sq
    a,b,c,d = vertices(sq)
    return LagrangeTriangle(a,x,b), LagrangeTriangle(b,x,c), LagrangeTriangle(c,x,d), LagrangeTriangle(d,x,a)
end

function decompose(sq::ReferenceSquare,x,::ReferenceSquare)
    @assert x ∈ sq
    a,b,c,d = vertices(sq)
    return  LagrangeSquare(x,(0,x[2]),(0,0),(x[1],0)),
            LagrangeSquare(x,(x[1],0),(1,0),(1,x[2])),
            LagrangeSquare(x,(1,x[2]),(1,1),(x[1],1)),
            LagrangeSquare(x,(x[1],1),(0,1),(0,x[2]))
end
