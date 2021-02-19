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

function decompose(sq::ReferenceSquare,x,target::ReferenceSquare)
    @assert x ∈ sq
    a,b,c,d = vertices(sq)
    return  LagrangeRectangle(x,(0,x[2]),(0,0),(x[1],0)),
            LagrangeRectangle(x,(x[1],0),(1,0),(1,x[2])),
            LagrangeRectangle(x,(1,x[2]),(1,1),(x[1],1)),
            LagrangeRectangle(x,(x[1],1),(0,1),(0,x[2]))
end    