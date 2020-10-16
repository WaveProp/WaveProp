abstract type AbstractQuadrature{N,T} end

struct GenericQuadrature{N,T}
    nodes::Vector{Point{N,T}}
    weights::Vector{T}
end

function GenericQuadrature{N,T}() where {N,T}
    GenericQuadrature{N,T}([],[])
end