"""
    struct Domain

Represent a physical domain as a union of elementary entities.
"""
struct Domain
    entities::Vector{ElementaryEntity}
end
Domain(ω::ElementaryEntity) = Domain([ω,])
Domain() = Domain(ElementaryEntity[])

"""
    entities(Ω::Domain)

Return a vector of all elementary entities making up a domain.
"""
entities(Ω::Domain) = Ω.entities

"""
    geometric_dimension(Ω::Domain)

If all entities forming the domain have the same `geometric_dimension`, return
that value; otherwise throw an assertion error. 
"""
function geometric_dimension(Ω::Domain)
    dims = [geometric_dimension(ent) for ent in entities(Ω)] |> unique!
    msg = "currently not able to handle a `Domain` with entities of different
           geometrical dimension. Got entities of dimensions ($dims)"
    @assert length(dims) == 1 msg
    return first(dims)
end

"""
    length(Ω:::Domain)

The length of a domain corresponds to the number of elementary entities that make it.
"""
Base.length(Ω::Domain) = length(entities(Ω))

"""
    skeleton(Ω::Domain)

Return all the boundaries of the domain, i.e. the domain's skeleton.
"""
skeleton(Ω::Domain) = union(Domain.(boundary.(entities(Ω)))...)

"""
    ===(Ω1::Domain,Ω2::Domain)

Two `Domain`s are equal if all their entities are equal.
"""
function Base.:(==)(Ω1::Domain, Ω2::Domain)
    return entities(Ω1) == entities(Ω2)
end


"""
    in(ω::ElementaryEntity,Ω::Domain)

Check whether an `ElementaryEntity` belongs to a `Domain`. Not that his only
does a *shallow* comparisson, meaning it only checks that `ω` is not one of the
entities in `entities(Ω)`; thus `ω` could still belong equal one of the
boundaries of an entity in `entities(Ω)`.
"""
Base.in(ω::ElementaryEntity, Ω::Domain) = in(ω, entities(Ω))

Base.getindex(Ω::Domain, i) = entities(Ω)[i]
Base.lastindex(Ω::Domain) = length(Ω)

"""
    iterate()

Iterating over domain means iterating over its entities. 
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
    union(Ωs::Domain...)

Union of domains. This done perform a union on the true geometric objects in
`entities(Ω)`, but simply on their identification through `(dim,tag)`. 
"""
function Base.union(Ω1::Domain,Ωs::Domain...)
    Domain(Vector{ElementaryEntity}(union(entities(Ω1),entities.(Ωs)...)))
end    
Base.union(Ω::Domain) = Domain(unique(entities(Ω)))

# function Base.union(Ωs::Domain...)
#     Domain(Vector{ElementaryEntity}(unique(vcat(entities.(Ωs)...))))
# end

"""
    assertequaldim(Ω1::Domain,Ω2::Domain)

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