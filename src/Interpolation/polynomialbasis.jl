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

# monomial basis over reference square
Base.length(::MonomialBasis{ReferenceSquare,P}) where {P} = (P)^2 + 1
Base.eltype(::MonomialBasis{ReferenceSquare}) = Monomial{2}
function Base.iterate(basis::MonomialBasis{ReferenceSquare,P},state=SVector(0,0)) where {P}
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
