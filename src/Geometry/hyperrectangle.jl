"""
    HyperRectangle{N,T}

Hyperrectangle in `N` dimensions described by a `low_corner::SVector{N,T}` and a
`high_corner::SVector{N,T}`
"""
struct HyperRectangle{N,T}
    low_corner::SVector{N,T}
    high_corner::SVector{N,T}
end

low_corner(r::HyperRectangle)  = r.low_corner
high_corner(r::HyperRectangle) = r.high_corner

# promote
HyperRectangle(l::Tuple,h::Tuple)     = HyperRectangle(SVector(l),SVector(h))
HyperRectangle(l::SVector,h::SVector) = HyperRectangle(promote(l,h)...)
# 1d case
HyperRectangle(a::Number,b::Number) = HyperRectangle(SVector(a),SVector(b))

Base.:(==)(h1::HyperRectangle, h2::HyperRectangle) = (h1.low_corner == h2.low_corner) && (h1.high_corner == h2.high_corner)

Base.isapprox(h1::HyperRectangle,h2::HyperRectangle;kwargs...) = isapprox(h1.low_corner,h2.low_corner;kwargs...) && isapprox(h1.high_corner,h2.high_corner;kwargs...)

Base.in(point,h::HyperRectangle) = all(h.low_corner .<= point .<= h.high_corner)

Base.eltype(h::HyperRectangle{N,T}) where {N,T}     = T

ambient_dimension(h::HyperRectangle{N}) where {N}   = N

geometric_dimension(h::HyperRectangle{N}) where {N} = N

diameter(cub::HyperRectangle) = norm(high_corner(cub) .- low_corner(cub),2)

function bounding_box(data)
    isempty(data)  && (error("data cannot be empty") )
    low_corner  = first(data)
    high_corner = first(data)
    for pt in data
        low_corner  = min.(low_corner,pt)
        high_corner = max.(high_corner,pt)
    end
    return HyperRectangle(low_corner,high_corner)
end

center(rec::HyperRectangle) = (rec.low_corner + rec.high_corner) / 2

radius(rec::HyperRectangle) = diameter(rec) / 2
