"""
    HyperRectangle{N,T}

Hyperrectangle in `N` dimensions given by `low_corner::SVector{N,T}` and `high_corner::SVector{N,T}`
"""
struct HyperRectangle{N,T}
    low_corner::SVector{N,T}
    high_corner::SVector{N,T}
end

# promote 
HyperRectangle(l::Tuple,h::Tuple) = HyperRectangle(SVector(l),SVector(h))
# 1d case
HyperRectangle(a::Number,b::Number) = HyperRectangle(SVector(a),SVector(b))

Base.:(==)(h1::HyperRectangle, h2::HyperRectangle) = (h1.low_corner == h2.low_corner) && (h1.high_corner == h2.high_corner)
Base.in(point,h::HyperRectangle) = all(h.low_corner .<= point .<= h.high_corner)

Base.eltype(h::HyperRectangle{N,T}) where {N,T}    = T
ambient_dimension(h::HyperRectangle{N}) where {N}          = N
geometric_dimension(h::HyperRectangle{N}) where {N}          = N

diameter(cub::HyperRectangle) = norm(cub.high_corner .- cub.low_corner,2)

function bounding_box(data)
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
