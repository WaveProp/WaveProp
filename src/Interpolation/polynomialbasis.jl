"""
    abstract type PolynomialBasis{D,P}
    
A polynomial basis over `D` of order `≤P`.
"""
abstract type PolynomialBasis{D,P} end

struct MonomialBasis{D,P} <: PolynomialBasis{D,P} end
MonomialBasis(d::AbstractReferenceShape,p::Int) = MonomialBasis{typeof(d),p}()

# monomial basis over reference segment
Base.length(::MonomialBasis{ReferenceLine,P}) where {P} = P+1
Base.eltype(::MonomialBasis{ReferenceLine}) = Monomial{1}
function Base.iterate(basis::MonomialBasis{ReferenceLine,P},state=0) where {P}
    if state > P
        return nothing        
    else
        return Monomial(state),state+1
    end
end  

# monomial basis over reference triangle
Base.length(::MonomialBasis{ReferenceTriangle,P}) where {P} = Int((P+1)*(P+2)/2)
Base.eltype(::MonomialBasis{ReferenceTriangle}) = Monomial{2}
function Base.iterate(basis::MonomialBasis{ReferenceTriangle,P},state=(0,0)) where {P}
    i,j = state    
    if j > P
        return nothing        
    else
        if i == P - j
            next = (0,j+1)    
        else
            next = (i+1,j)
        end
        return Monomial(state),next
    end
end  

# monomial basis over reference square
Base.length(::MonomialBasis{ReferenceSquare,P}) where {P} = (P)^2 + 1
Base.eltype(::MonomialBasis{ReferenceSquare}) = Monomial{2}
function Base.iterate(basis::MonomialBasis{ReferenceSquare,P},state=(0,0)) where {P}
    if state[1] > P || state[2] > P
        return nothing        
    else
        if state[1]==state[2]
            return Monomial(state),state .+ (1,0)
        elseif state[1]>state[2]
            return Monomial(state),state .+ (-1,1)    
        elseif state[2]>state[1]
            return Monomial(state),state .+ (1,0)        
        end
    end
end  

# vandermond matrix
function vandermond(xi,basis::PolynomialBasis)
    return [p(x) for x in xi, p in basis]
end

function lagrange_basis(xi,m::MonomialBasis)
    L = vandermond(xi,m)
    C = L\I
    return C
end    


struct LagrangeBasis{D,P} <: PolynomialBasis{D,P} 
end    

length(b::LagrangeBasis{ReferenceTriangle,P}) where {P} = (P+1)*(P+2)/2 |> Int

"""
    (b::LagrangeBasis)(x)

Evaluate all base elements in `b` at the point `x`. Return a `StaticVector` of
length `length(b)`.
"""
function (b::LagrangeBasis{ReferenceTriangle,1})(x)
    SVector(1 - x[1] - x[2], x[1], x[2])
end    

function (b::LagrangeBasis{ReferenceTriangle,0})(x)
    SVector(1)
end    

function (b::LagrangeBasis{ReferenceTriangle,0})(x̂::SVector{<:Any,<:SVector})
    mapreduce(x->b(x),hcat,x̂)    
end    

function (b::LagrangeBasis{ReferenceTriangle,1})(x̂::SVector{<:Any,<:SVector})
    mapreduce(x->b(x),hcat,x̂)    
end    
