"""
    abstract type AbstractReferenceShape{N}
    
A reference polygon in `ℜᴺ`.

Reference shapes are used mostly for defining [`AbstractElement`](@ref)s as
transformations mapping an `AbstractReferenceShape` into some region of `ℜᴹ`. 

See e.g. [`ReferenceLine`](@ref) or [`ReferenceTriangle`](@ref) for some
examples of concrete subtypes.
"""
abstract type AbstractReferenceShape{N} end

ambient_dimension(::Type{<:AbstractReferenceShape{N}}) where {N} = N
geometric_dimension(::Type{<:AbstractReferenceShape{N}}) where {N} = N
ambient_dimension(::AbstractReferenceShape{N}) where {N} = N
geometric_dimension(::AbstractReferenceShape{N}) where {N} = N

"""
    struct ReferenceLine
    
Singleton type representing the `[0,1]` segment.
"""
struct ReferenceLine <: AbstractReferenceShape{1}
end
Base.in(x,::ReferenceLine) = 0 ≤ x[1] ≤ 1
center(::Type{ReferenceLine}) = 0.5
center(::ReferenceLine)       = 0.5

"""
    struct ReferenceTriangle
    
Singleton type representing the triangle with vertices `(0,0),(0,1),(1,0)`
"""
struct ReferenceTriangle <: AbstractReferenceShape{2} 
end    
Base.in(x,::ReferenceTriangle) = 0 ≤ x[1] ≤ 1 && 0 ≤ x[2] ≤ 1 - x[1]

"""
    struct ReferenceSquare
    
Singleton type representing the square with vertices `(0,0),(0,1),(1,1),(1,0)`
"""
struct ReferenceSquare <: AbstractReferenceShape{2} 
end    
Base.in(x,::ReferenceSquare) = 0 ≤ x[1] ≤ 1 && 0 ≤ x[2] ≤ 1
center(::Type{ReferenceSquare}) = SVector(0.5,0.5)
center(::ReferenceSquare)       = SVector(0.5,0.5)

"""
    struct ReferenceTetrahedron
    
Singleton type representing the tetrahedron with vertices `(0,0,0),(0,0,1),(0,1,0),(1,0,0)`
"""
struct ReferenceTetrahedron <: AbstractReferenceShape{3} 
end    
Base.in(x,::ReferenceTetrahedron) = 0 ≤ x[1] ≤ 1 && 0 ≤ x[2] ≤ 1 - x[1] && 0 ≤ x[3] ≤ 1 - x[1] - x[2]

# TODO: generalize structs above to `ReferenceSimplex{N}` and `ReferenceCuboid{N}`