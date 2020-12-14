"""
    hcubature(f,s::AbstractReferenceShape;kwargs...)

Integrate `f` over the reference shape `s`.
"""
hcubature(f,l::AbstractReferenceShape;kwargs...) = hcubature(f,typeof(l);kwargs...)

hcubature(f,::Type{ReferenceLine};kwargs...)     = hcubature(f,(0,),(1,);kwargs...)

hcubature(f,::Type{ReferenceSquare};kwargs...)   = hcubature(f,(0,0),(1,1);kwargs...)
