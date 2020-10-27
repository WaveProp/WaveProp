quadgk(f,l::AbstractReferenceShape) = quadgk(f,typeof(l))

quadgk(f,::Type{ReferenceLine}) = quadgk(f,0,1)[1]

function quadgk(f,::Type{ReferenceSquare})
    I    = x-> quadgk(y->f(Point(x,y)),0,1)[1]
    out  = quadgk(I,0,1)[1]
end    

function quadgk(f,::Type{ReferenceTriangle})
    I    = x -> quadgk(y->f(Point(x,y)),0,1-x)[1]
    out  = quadgk(I,0,1)[1]
end    