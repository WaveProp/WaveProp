"""
    abstract type AbstractEntity
    
Entitty of geometrical nature. 
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
        msg = "an elementaty entities in the boundary has a wrong dimension"
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

"""Return geometric the dimension of the elementary entity."""
geometric_dimension(ω::ElementaryEntity) = ω.dim

"""Return the unique tag (for a given dimension) of the elementary entity."""
tag(ω::ElementaryEntity) = (geometric_dimension(ω), ω.tag)

"""Return the vector of elementary entities making the boundary."""
boundary(ω::ElementaryEntity) = ω.boundary

"""Test equality between elementary entities."""
function Base.:(==)(Ω1::ElementaryEntity, Ω2::ElementaryEntity)
    Ω1.dim == Ω2.dim || return false
    Ω1.tag == Ω2.tag || return false
    Ω1.boundary == Ω2.boundary || return false
    return true
end

"""
Represent a physical domain (union of elementary entities).
"""
struct Domain
    entities::Vector{ElementaryEntity}
end
Domain(ω::ElementaryEntity) = Domain([ω,])
Domain() = Domain(ElementaryEntity[])


"""Return a vector of all elementary entities making up a domain."""
entities(Ω::Domain) = Ω.entities

"""Return the dimension of the domain."""
function geometric_dimension(Ω::Domain)
    dims = [geometric_dimension(ent) for ent in entities(Ω)] |> unique!
    msg = "currently not able to handle a `Domain` with entities of different
           geometrical dimension. Got entities of dimensions ($dims)"
    @assert length(dims) == 1 msg
    return first(dims)
end

"""
The length of a domain corresponds to the number of elementary entities that make it.
"""
Base.length(Ω::Domain) = length(entities(Ω))

"""Return all the boundaries of the domain."""
skeleton(Ω::Domain) = union(Domain.(boundary.(entities(Ω)))...)

"""Test equality between two sub-domains."""
function Base.:(==)(Ω1::Domain, Ω2::Domain)
    return entities(Ω1) == entities(Ω2)
end


"""
Check whether an `ElementaryEntity` belongs to a `Domain`.
"""
Base.in(ω::ElementaryEntity, Ω::Domain) = in(ω, entities(Ω))

Base.getindex(Ω::Domain, i) = entities(Ω)[i]
Base.lastindex(Ω::Domain) = length(Ω)

"""
Domain can be viewed as an iterable collection.
"""
function Base.iterate(Ω::Domain, state=1)
    # Check if we are done iterating
    if state > length(Ω)
        return nothing
    end
    # Return (result, state)
    return (Ω[state], state+1)
end

"""
Check whether a Domain is empty.
"""
Base.isempty(Ω::Domain) = length(entities(Ω)) == 0

"""
Domain which is not in the intersection of domains Ω1 and Ω2.
"""
function Base.setdiff(Ω1::Domain, Ω2::Domain)
    Domain(setdiff(entities(Ω1), entities(Ω2)))
end

"""
Union of domains.
"""
function Base.union(Ωs::Domain...)
    Domain(Vector{ElementaryEntity}(unique(vcat(entities.(Ωs)...))))
end

"""
Check that two domains have same dimension.

If one of the domain (or both) are empty, the assertion is assumed to be true.
"""
function assertequaldim(Ω1::Domain, Ω2::Domain)
    if isempty(Ω1) || isempty(Ω2) return nothing end
    msg = "The dimension of the first domain should be equal to the dimension
    of the second domain."
    @assert geometric_dimension(Ω1) == geometric_dimension(Ω2) msg
end

"""
Intersection between domain Ω1 and domain Ω2.
"""
function Base.intersect(Ω1::Domain, Ω2::Domain)
    assertequaldim(Ω1, Ω2)
    Ωinter = Domain(intersect(entities(Ω1), entities(Ω2)))
    if isempty(Ωinter)
        if isempty(Ω1) || isempty(Ω2)
            return Domain()
        else
            return intersect(skeleton(Ω1), skeleton(Ω2))
        end
    else
        return Ωinter
    end
end

function Base.push!(Ω::Domain,ent::ElementaryEntity)
    push!(entities(Ω),ent)
end    

"""
Determine whether every element of domain Ω1 is also in domain Ω2.
"""
function Base.issubset(Ω1::Domain, Ω2::Domain)
    assertequaldim(Ω1, Ω2)
    return issubset(entities(Ω1), entities(Ω2))
end

"""
Remove domain Ω1 from domain Ω2.
"""
function remove(Ω1::Domain, Ω2::Domain)
    assertequaldim(Ω1, Ω2)
    return Domain(setdiff(entities(Ω2), intersect(entities(Ω1), entities(Ω2))))
end

"""Return the internal boundaries inside a domain."""
function internal_boundary(Ω::Domain)
    Ω1 = Domain(Ω[1])
    γ = Domain()
    for ω2 in Ω[2:end]
        Ω2 = Domain(ω2)
        γ1 = Domain(vcat(boundary.(entities(Ω1))...))
        γ2 = Domain(vcat(boundary.(entities(Ω2))...))
        γ = union(γ, intersect(γ1, γ2))
        Ω1 = union(Ω1, Ω2)
    end
    return γ
end

"""Return the external boundaries inside a domain."""
function external_boundary(Ω::Domain)
    return remove(internal_boundary(Ω), skeleton(Ω))
end

"""
    boundary(Ω)

Return a domain comprising the external boundary of Ω. 

See also: [`external_boundary`](@ref)
"""
boundary(Ω) = external_boundary(Ω)

"""
Return all tags of the elementary entities in the domain `Ω` corresponding to the dimension `d`.
"""
function tags(Ω::Domain, d::Integer)
    if isempty(Ω)
        return Tuple{Int64,Int64}[]
    elseif d == geometric_dimension(Ω)
        return Vector{Tuple{Int64,Int64}}(vcat(tag.(entities(Ω))))
    elseif d < geometric_dimension(Ω)
        return unique(tags(skeleton(Ω),d))
    else
        error("Asking for tags with dimension > dimension of domain")
    end
end
function tags(Ω::Domain)
    isempty(Ω) ? Tuple{Int64, Int64}[] : tags(Ω, geometric_dimension(Ω))
end

"""
Return all tags of the elementary entities in the domain `Ω` corresponding to the dimensions contained in `dims`.
"""
function tags(Ω::Domain, dims::Vector{T}) where T <: Integer
    tgs = Vector{Tuple{Int64, Int64}}(undef, 0)
    for d in dims push!(tgs, tags(Ω, d)...) end
    return tgs
end