struct Monomial{N}
    θ::SVector{N,Int}
end 
Monomial(p::Int) = Monomial(SVector(p))

dimension(m::Monomial{N}) where {N} = N

Base.exponent(m::Monomial) = m.θ

function (m::Monomial)(x::SVector) 
    @assert length(x) == dimension(m)
    x .^ exponent(m) |> prod
end
(m::Monomial)(x::Number)  = m(SVector(x))
(m::Monomial)(x::Vector)  = m(SVector(x))