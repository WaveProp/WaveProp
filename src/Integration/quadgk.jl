"""
    quadgk(f,s::AbstractReferenceShape;kwargs...)

Integrate `f` over the reference shape `s`.
"""
quadgk(f,::Type{ReferenceLine};kwargs...)     = quadgk(f,0,1;kwargs...)

quadgk(f,l::AbstractReferenceShape;kwargs...) = quadgk(f,typeof(l);kwargs...)
