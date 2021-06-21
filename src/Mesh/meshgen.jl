# Basic utilities for meshing parametric surfaces. Produces a `GenericMesh`

"""
    meshgen(Ω::Domain;h::Tuple)

Generate a `GenericMesh` for the domain `Ω`. Requires the entities forming `Ω`
to `ParametricEntity`.
"""
function meshgen(Ω::Domain,sz)
    # extract the ambient dimension for these entities (i.e. are we in 2d or
    # 3d). Only makes sense if all entities have the same ambient dimension.
    N  = ambient_dimension(first(Ω))
    @assert all(p->ambient_dimension(p)==N,entities(Ω))
    mesh = GenericMesh{N,Float64}()
    meshgen!(mesh,Ω,sz) # fill in
    return mesh
end

"""
    meshgen!(mesh,Ω;h)

Similar to [`meshgen`](@ref), but append entries to `mesh`.
"""
function meshgen!(mesh::GenericMesh,Ω::Domain,sz)
    for ent in Ω
        ent isa ParametricEntity || error("meshgen! only works on parametric entites")
        _meshgen!(mesh,ent,sz)
    end
    return mesh
end

function _meshgen!(mesh::GenericMesh,ent::ParametricEntity,sz)
    # extract relevant fields and mesh the entity
    f = parametrization(ent)
    d = domain(ent)
    els    = _meshgen(f,d,sz)
    # push related information to mesh
    E      = eltype(els)
    vals   = get!(mesh.elements,E,Vector{E}())
    istart = length(vals) + 1
    append!(vals,els)
    iend   = length(vals)
    haskey(mesh.ent2tags,ent) && @debug "skipping entity $(key(ent)): already present in mesh"
    mesh.ent2tags[ent] = Dict(E=>collect(istart:iend)) # add key
    return mesh
end

"""
    _meshgen(f,d,sz)

Create a UniformCartesianMesh` of `d` push-forward map. The cartesian mesh has size `sz`, and is uniform in parameter
space.

!!! warning
    For parametrizations with large `jacobian`, the mesh size in physical space can
    be far from uniform.
"""
function _meshgen(f,d,sz)
    grid =UniformCartesianMesh(d,sz)
    iter = ElementIterator(grid)
    els = [ParametricElement(f,d) for d in iter]
    return els
end

# Element iterator interface
function Base.size(iter::ElementIterator{<:ParametricElement,<:GenericMesh})
    E = eltype(iter)
    M = mesh(iter)
    els::Vector{E} = M.elements[E]
    return (length(els),)
end

function Base.getindex(iter::ElementIterator{<:ParametricElement,<:GenericMesh},i::Int)
    E = eltype(iter)
    M = mesh(iter)
    els::Vector{E} = M.elements[E]
    return els[i]
end

function Base.iterate(iter::ElementIterator{<:ParametricElement,<:GenericMesh}, state=1)
    state > length(iter) && (return nothing)
    iter[state], state + 1
end
