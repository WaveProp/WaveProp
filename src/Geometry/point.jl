const Point{N,T} = SVector{N,T}

geometric_dimension(::Type{<:Point}) = 0
geometric_dimension(pt::Point) = geometric_dimension(typeof(pt))