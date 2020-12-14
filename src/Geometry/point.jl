# a point is just an SVector for us
geometric_dimension(::Type{<:SVector}) = 0
geometric_dimension(pt::SVector) = geometric_dimension(typeof(pt))
ambient_dimension(::SVector{N}) where N = N
ambient_dimension(::Type{<:SVector{N}}) where N = N