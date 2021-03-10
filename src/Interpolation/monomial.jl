struct Monomial{N}
    θ::SVector{N,Int}
end
Monomial(p::Int) = Monomial(SVector(p))
Monomial(p::Tuple) = Monomial(SVector(p))

dimension(m::Monomial{N}) where {N} = N

Base.exponent(m::Monomial) = m.θ

function (m::Monomial)(x::SVector)
    @assert length(x) == dimension(m)
    x .^ exponent(m) |> prod
end
(m::Monomial)(x)  = m(SVector(x))
