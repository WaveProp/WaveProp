const Point{N,T} = SVector{N,T}

geometric_dimension(::Type{<:Point}) = 0
geometric_dimension(pt::Point) = geometric_dimension(typeof(pt))
ambient_dimension(::Point{N}) where N = N
ambient_dimension(::Type{<:Point{N}}) where N = N