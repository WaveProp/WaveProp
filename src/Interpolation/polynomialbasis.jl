"""
    abstract type PolynomialBasis{D,P}

A polynomial basis over `D` of order `≤P`.
"""
abstract type PolynomialBasis{D,P} end

"""
        struct MonomialBasis{D,P}

A basis of monomials spanning all polynomials of *order* `P` over `D`.

For multi-dimensional domains `D`, the precise definion of *order* depends in
fact on `D`. When `D` is a unit hypercube, the order of ``xᵅ`` where `α` is a
multi-index (e.g. an `NTuple`) is infinity norm. When `D` is a reference
simplex, the norm is the *L¹-norm*.
"""
struct MonomialBasis{D,P} <: PolynomialBasis{D,P} end
MonomialBasis(d::AbstractReferenceShape,p::Int) = MonomialBasis{typeof(d),p}()

Base.getindex(b::MonomialBasis,::Colon) = SVector{length(b)}(collect(b))

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

"""
    struct LagrangeBasis{D,P}

A Lagrange basis of order `P` over the reference shape
`D<:AbstractReferenceShape`. Spans the same space [`MonomialBasis{D,P}`](@ref)

Can be used as a function through `(b::LagrangeBasis)(x)`, where it Evaluate all
base elements in `b` at the point `x`. Returns a `StaticVector` of length
`length(b)`.
"""
struct LagrangeBasis{D,P} <: PolynomialBasis{D,P}
end

Base.getindex(b::LagrangeBasis,I) = b[:][I]

Base.length(::LagrangeBasis{ReferenceLine,P}) where {P} = P+1 |> Int
Base.length(::LagrangeBasis{ReferenceTriangle,P}) where {P} = (P+1)*(P+2)/2 |> Int
Base.length(::LagrangeBasis{ReferenceTetrahedron,P}) where {P} = (P+1)*(P+2)*(P+3)/6 |> Int


function (b::LagrangeBasis{D,0})(x::SVector{<:Any,<:Number}) where {D}
    SVector(1)
end
Base.getindex(b::LagrangeBasis{D,0},::Colon) where {D} = (x->SVector{1},)

function (b::LagrangeBasis{ReferenceLine,1})(x::SVector{<:Any,<:Number})
    map(f->f(x),b[:])
end
Base.getindex(b::LagrangeBasis{ReferenceLine,1},::Colon) = (x->1-x[1],x->x[1])

function (b::LagrangeBasis{ReferenceTriangle,1})(x::SVector{<:Any,<:Number})
    SVector(1 - x[1] - x[2], x[1], x[2])
end
(b::LagrangeBasis{ReferenceTriangle,1})()              = (x -> 1 - x[1] - x[2], x -> x[1], x -> x[2])
Base.getindex(b::LagrangeBasis{ReferenceTriangle,1},I) = b()[I]

function (b::LagrangeBasis{ReferenceTetrahedron,1})(x::SVector{<:Any,<:Number})
    SVector(1 - x[1] - x[2] - x[3], x[1], x[2], x[3])
end
Base.getindex(b::LagrangeBasis{ReferenceTetrahedron,1},::Colon) = (x -> 1 - x[1] - x[2] - x[3], x -> x[1], x -> x[2], x -> x[3])

function (b::LagrangeBasis{D,p})(x̂::SVector{<:Any,<:SVector}) where {D,p}
    mapreduce(x->b(x),hcat,x̂)
end

trace(::LagrangeBasis{ReferenceTriangle,p}) where {p}    = LagrangeBasis{ReferenceLine,p}()
trace(::LagrangeBasis{ReferenceTetrahedron,p}) where {p} = LagrangeBasis{ReferenceTriangle,p}()

struct gradLagrangeBasis{D,P} <: PolynomialBasis{D,P}
end

Base.length(::gradLagrangeBasis{ReferenceLine,P}) where {P} = P+1 |> Int
Base.length(::gradLagrangeBasis{ReferenceTriangle,P}) where {P} = (P+1)*(P+2)/2 |> Int
Base.length(::gradLagrangeBasis{ReferenceTetrahedron,P}) where {P} = (P+1)*(P+2)*(P+3)/6 |> Int

function (b::gradLagrangeBasis{ReferenceLine,1})(x::SVector{<:Any,<:Number})
    SVector(SVector(-1),
            SVector(1))
end

function (b::gradLagrangeBasis{ReferenceTriangle,1})(x::SVector{<:Any,<:Number})
    SVector(SVector(-1, -1),
            SVector(1, 0),
            SVector(0, 1))
end
(b::gradLagrangeBasis{ReferenceTriangle,1})()               = ( x -> SVector(-1, -1),
                                                                x -> SVector(1, 0),
                                                                x-> SVector(0, 1)
                                                                )
Base.getindex(b::gradLagrangeBasis{ReferenceTriangle,1},I)  = b()[I]

function (b::gradLagrangeBasis{ReferenceTetrahedron,1})(x::SVector{<:Any,<:Number})
    SVector(SVector(-1, -1, -1),
            SVector(1, 0, 0),
            SVector(0, 1, 0),
            SVector(0, 0, 1))
end

function (b::gradLagrangeBasis{D,p})(x̂::SVector{<:Any,<:SVector}) where {D,p}
    mapreduce(x->b(x),hcat,x̂)
end


grad(::LagrangeBasis{D,p}) where {D,p} = gradLagrangeBasis{D,p}()


# FIXME: need a better way to related grad to functions in a basis
blist = (LagrangeBasis{ReferenceTriangle,1}(),)
for b in blist
    for i in 1:length(b)
        f = b[i]
        @eval grad(::typeof($f)) = grad($b)[$i]
    end
end
