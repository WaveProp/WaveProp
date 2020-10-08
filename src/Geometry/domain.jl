"""
Indivisible geometric entity.
"""
struct ElementaryEntity
    dim::UInt8
    tag::Int64
    boundary::Vector{ElementaryEntity}
    function ElementaryEntity(d::Integer, t::Integer, boundary::Vector{ElementaryEntity})
        msg = "An elementaty entities in the boundary has a wrong dimension"
        for b in boundary
            @assert dim(b) == d-1 msg
        end
        new(d, t, boundary)
    end
end


"""Return the dimension of the elementary entity."""
dim(ω::ElementaryEntity) = ω.dim


"""Return the unique tag (for a given dimension) of the elementary entity."""
tag(ω::ElementaryEntity) = (dim(ω), ω.tag)


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
dim(Ω::Domain) = dim(entities(Ω)[1])


"""
The length of a domain corresponds to the number of elementary entities that make it.
"""
Base.length(Ω::Domain) = length(entities(Ω))


"""Return all the boundaryies of the domain."""
skeleton(Ω::Domain) = union(Domain.(boundary.(entities(Ω)))...)


"""Test equality between two sub-domains."""
function Base.:(==)(Ω1::Domain, Ω2::Domain)
    return entities(Ω1) == entities(Ω2)
end


"""
Check whether an `ElementaryEntity` belongs to a `Domain`.
"""
Base.in(ω::ElementaryEntity, Ω::Domain) = Base.in(ω, entities(Ω))


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
    @assert dim(Ω1) == dim(Ω2) msg
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


"""Return the exterior boundaries inside a domain."""
function exterior_boundary(Ω::Domain)
    return remove(internal_boundary(Ω), skeleton(Ω))
end


"""
Return all tags of the elementary entities in the domain `Ω` corresponding to the dimension `d`.
"""
function tags(Ω::Domain, d::Integer)
    if isempty(Ω)
        return Tuple{Int64,Int64}[]
    elseif d == dim(Ω)
        return Vector{Tuple{Int64,Int64}}(vcat(tag.(entities(Ω))))
    elseif d < dim(Ω)
        return unique(tags(skeleton(Ω),d))
    else
        error("Asking for tags with dimension > dimension of domain")
    end
end
function tags(Ω::Domain)
    isempty(Ω) ? Tuple{Int64, Int64}[] : tags(Ω, dim(Ω))
end


"""
Return all tags of the elementary entities in the domain `Ω` corresponding to the dimensions contained in `dims`.
"""
function tags(Ω::Domain, dims::Vector{T}) where T <: Integer
    tgs = Vector{Tuple{Int64, Int64}}(undef, 0)
    for d in dims push!(tgs, tags(Ω, d)...) end
    return tgs
end