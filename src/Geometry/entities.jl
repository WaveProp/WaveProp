"""
    abstract type AbstractEntity

Entity of geometrical nature. Identifiable throught its `(dim,tag)` key.
"""
abstract type AbstractEntity end

"""
    key(e::AbstractEntity)

The `(dim,tag)` pair used as a key to identify various abstract entities.
"""
key(e::AbstractEntity) = geometric_dimension(e),tag(e)

# reasonable defaults which assume the filds `tag` and `dim` and `boundary`
# fields exist. Some
# `AbstractEntities` need to override this method.
tag(e::AbstractEntity) = e.tag
geometric_dimension(e::AbstractEntity) = e.dim
boundary(e::AbstractEntity) = e.boundary

function Base.show(io::IO,ent::AbstractEntity)
    T = typeof(ent)
    d = geometric_dimension(ent)
    t = tag(ent)
    print(io,"$T with (dim,tag)=($d,$t)")
end

"""
    ==(Ω1::AbstractEntity,Ω2::AbstractEntity)

Two entities are considered equal
`geometric_dimension(Ω1)==geometric_dimension(Ω2)` and
`abs(tag(Ω1))=abs(tag(Ω2))`. For entities of co-dimension one, the sign of
`tag(Ω)` is used to determine the orientation of the normal vector.

Notice that this implies `dim` and `tag` of an elementary entity should uniquely
define it (up to the sign of `tag`), and therefore global variables like
[`TAGS`](@ref) are needed to make sure newly created [`AbstractEntities`](@ref)
have a new `(dim,tag)` identifier.
"""
function Base.:(==)(Ω1::AbstractEntity, Ω2::AbstractEntity)
    d1,t1 = geometric_dimension(Ω1),tag(Ω1)
    d2,t2 = geometric_dimension(Ω2),tag(Ω2)
    d1 == d2  || (return false)
    abs(t1) == abs(t2) || (return false)
    # boundary(Ω1) == boundary(Ω2) || return false # this should not be needed
    return true
end
Base.hash(ent::AbstractEntity,h::UInt)= hash((geometric_dimension(ent),abs(tag(ent))),h)

function normal(ent::AbstractEntity, u)
    s = tag(ent) |> sign
    jac::SMatrix = jacobian(ent, u)
    s*normal(jac)
end
function normal(jac::SMatrix{N,M}) where {N,M}
    msg = "computing the normal vector requires the element to be of co-dimension one."
    @assert (N - M == 1) msg
    if M == 1 # a line in 2d
        t = jac[:,1] # tangent vector
        n = SVector(t[2], -t[1])
        return n / norm(n)
    elseif M == 2 # a surface in 3d
        t₁ = jac[:,1]
        t₂ = jac[:,2]
        n  = cross(t₁, t₂)
        return n / norm(n)
    else
        notimplemented()
    end
end

"""
    struct ElementaryEntity <: AbstractEntity

The most basic representation of an [`AbstractEntity`](@ref).

# Fields:
- `dim::UInt8`: the geometrical dimension of the entity (e.g. line has `dim=1`,
  surface has `dim=2`, etc)
- `tag::Int64`: an integer tag associated to the entity
- `boundary::Vector{AbstractEntity}`: the entities of dimension `dim-1`
  forming the entity's boundary
"""
struct ElementaryEntity <: AbstractEntity
    dim::UInt8
    tag::Int64
    boundary::Vector{<:AbstractEntity}
    function ElementaryEntity(d::Integer, t::Integer, boundary::Vector{<:AbstractEntity})
        msg = "an elementary entities in the boundary has the wrong dimension"
        for b in boundary
            @assert geometric_dimension(b) == d-1 msg
        end
        ent = new(d, t, boundary)
        # every entity gets added to a global variable ENTITIES so that we can
        # ensure the (d,t) pair is a UUID for an entity, and to easily retrieve
        # different entities.
        global_add_entity!(ent)
        return ent
    end
end

geometric_dimension(ω::ElementaryEntity) = ω.dim
tag(ω::ElementaryEntity) = ω.tag
boundary(ω::ElementaryEntity) = ω.boundary

"""
    ElementaryEntity(dim,tag)

Construct an [`ElementaryEntity`](@ref) with an empty boundary .
"""
function ElementaryEntity(dim,tag)
    ElementaryEntity(dim,tag,ElementaryEntity[])
end

ElementaryEntity(dim) = ElementaryEntity(dim,new_tag(dim))

function ElementaryEntity(;boundary,dim::Int=_compute_dim_from_boundary(boundary))
    t = new_tag(dim)
    ElementaryEntity(UInt8(dim),t,boundary)
end

function _compute_dim_from_boundary(boundary)
    dmin,dmax = extrema(geometric_dimension,boundary)
    @assert dmin == dmax "all entities in `boundary` must have the same dimension"
    dim = dmin+1
end

"""
    PointEntity{N,T} <: AbstractEntity

Zero-dimension geometrical entity. As a subtype of [`AbstractEntity`],(@ref) the
`(dim,tag)` of all created point entities get added to the global `ENTITIES`.
Intended usage is to build higher dimensionional entities, and *not* to
represent regular points such as grid points.
"""
struct PointEntity <: AbstractEntity
    tag::Int
    coords::SVector
end

coords(p::PointEntity) = p.coords
tag(p::PointEntity)    = p.tag

geometric_dimension(::PointEntity) = 0
ambient_dimension(p::PointEntity)   = length(coords(p))

#####################################################################

# Variables and functions to globally keep track of entities

#####################################################################

"""
    const TAGS::Dict{Int,Vector{Int}}

Global dictionary storing the used entity tags (the value) for a given dimension
(the key).
"""
const TAGS = Dict{Int,Vector{Int}}()

"""
    const ENTITIES

Global dictionary storing the used entity tags (the value) for a given dimension
(the key).
"""
const ENTITIES = Dict{Tuple{Int,Int},AbstractEntity}()

function global_add_entity!(ent::AbstractEntity)
    d,t = geometric_dimension(ent), tag(ent)
    _add_tag!(d,t) # add this tag to global list to make sure it is not used again
    msg = "overwriting ENTITIES: value in key ($d,$t) will be replaced"
    haskey(ENTITIES,(d,t)) && (@warn msg)
    ENTITIES[(d,t)] = ent
    return d,t
end

"""
    new_tag(dim)

Generate a unique tag for an `AbstractEntity` of dimension `dim`.

The implementation consists of adding one to the maximum value of `TAGS[dim]`

# See also: [`TAGS`](@ref).

"""
function new_tag(dim)
    if !haskey(TAGS,dim)
        return 1
    else
        tnew = maximum(TAGS[dim]) + 1
        return tnew
    end
end

function _add_tag!(dim,tag)
    if is_new_tag(dim,tag)
        # now add key
        if haskey(TAGS,dim)
            push!(TAGS[dim],tag)
        else
            TAGS[dim] = [tag,]
        end
    else
        # print warning but don't add duplicate tag
        msg  = "entity of dimension $dim and tag $tag already exists in TAGS.
        Creating a possibly duplicate entity."
        @warn msg
    end
    return TAGS
end

function is_new_tag(dim,tag)
    if haskey(TAGS,dim)
        existing_tags = TAGS[dim]
        if in(tag,existing_tags)
            return false
        end
    end
    return true
end

"""
    clear_entities!()

Empty the global variables used to keep track of the various entities
created.

# See also: [`ENTITIES`](@ref), [`TAGS`](@ref)
"""
function clear_entities!()
    empty!(TAGS)
    empty!(ENTITIES)
    nothing
end
