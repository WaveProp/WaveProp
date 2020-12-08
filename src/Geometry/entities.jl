"""
    abstract type AbstractEntity
    
Entity of geometrical nature. 
"""
abstract type AbstractEntity end

"""
    struct ElementaryEntity <: AbstractEntity
    
The basic building block of geometrical objects. 

# Fields:
- `dim::UInt8`: the geometrical dimension of the entity (e.g. line has `dim=1`, surface has `dim=2`, etc)
- `tag::Int64`: an integer tag associated to the entity
- `boundary::Vector{ElementaryEntity}`: the entities of dimension `dim-1` forming the entity's boundary  
"""
struct ElementaryEntity <: AbstractEntity
    dim::UInt8
    tag::Int64
    boundary::Vector{ElementaryEntity}
    function ElementaryEntity(d::Integer, t::Integer, boundary::Vector{ElementaryEntity})
        msg = "an elementary entities in the boundary has a wrong dimension"
        for b in boundary
            @assert geometric_dimension(b) == d-1 msg
        end
        # modify global variable TAGS by adding the new (d,t) for the
        # entity. It shows a warning if (d,t) already exists.
        _add_tag!(d,t) 
        new(d, t, boundary)
    end
end

"""
    ElementaryEntity(dim,tag)

Construct an [`ElementaryEntity`](@ref) with an empty boundary .
"""
function ElementaryEntity(dim,tag)
    ElementaryEntity(dim,tag,ElementaryEntity[])
end
ElementaryEntity(dim) = ElementaryEntity(dim,_new_tag(dim))

geometric_dimension(ω::ElementaryEntity) = ω.dim

"""
    tag(ω::ElementaryEntity)

Return the unique tag (for a given dimension) of the elementary entity.
"""
tag(ω::ElementaryEntity) = (geometric_dimension(ω), ω.tag)

"""
    boundary(ω::ElementaryEntity)

Return the vector of elementary entities making the boundary.
"""
boundary(ω::ElementaryEntity) = ω.boundary

"""
    ==(Ω1::ElementaryEntity,Ω2::ElementaryEntity) 

Two elementary entities are considered equal if their `dim` and `tag` fields
match.

Notice that this implies `dim` and `tag` of an elementary entity should uniquely
define it, and therefore global variables like [`TAGS`](@ref) are needed to make
sure newly created elementary entities have a new `(dim,tag)` identifier. 
"""
function Base.:(==)(Ω1::ElementaryEntity, Ω2::ElementaryEntity)
    Ω1.dim == Ω2.dim || return false
    Ω1.tag == Ω2.tag || return false
    Ω1.boundary == Ω2.boundary || return false
    return true
end