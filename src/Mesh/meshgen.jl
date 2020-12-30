"""
    meshgen(Ω::Domain;gridsize::Tuple) 

Generate a `GenericMesh` for the domain `Ω`. Requires the entities forming `Ω` to
`ParametricEntity` or `ClosedEntity`. 
"""
function meshgen(Ω::Domain;gridsize=floatmax())
    N       = ambient_dimension(first(Ω))
    @assert all(p->ambient_dimension(p)==N,Ω) # check that ambient dimensions are consistent
    mesh = GenericMesh{N,Float64}()
    meshgen!(mesh,Ω;gridsize)
end

"""
    meshgen!(mesh,Ω;gridsize)

Similar to [`meshgen`](@ref), but append entries to `mesh`. 
"""
function meshgen!(mesh::GenericMesh,Ω::Domain;gridsize)
    for ent in Ω
        ent isa ParametricEntity || ent isa AbstractParametricBody || error("unable to generate mesh for entity of type $(typeof(ent))")
        _meshgen!(mesh,ent;gridsize)
    end 
    return mesh       
end    

function _meshgen!(mesh::GenericMesh,ent::ParametricEntity;gridsize)
    # mesh the entity
    els = _meshgen(ent;gridsize)
    # push related information to mesh
    E   = eltype(els)
    E in etypes(mesh) || push!(etypes(mesh),E)
    vals   = get!(mesh.elements,E,Vector{E}())
    istart = length(vals) + 1
    append!(vals,els)
    iend   = length(vals)
    haskey(mesh.ent2tags,ent) && @warn "skipping entity $(key(ent)): already present in mesh"
    mesh.ent2tags[ent] = OrderedDict(E=>collect(istart:iend)) # add key
    return mesh
end    

function _meshgen!(mesh::GenericMesh,p::AbstractParametricBody;gridsize)
    mesh.ent2tags[p] = OrderedDict{DataType,Vector{Int}}() # no mesh is generated for ClosedEntity, only its boundary, so empty entry
    for ent in boundary(p)
        _meshgen!(mesh,ent;gridsize)
    end
    return mesh
end    

"""
    _meshgen(p::ParametricEntity;gridsize)

Mesh `p` by creating a `CartesianMesh` for its parameter space followed by the
push-forward map. Note that `gridsize` specifies the size in parameter space;
for parametrizations with large `jacobian` this is poor indication of the
mesh-size in real space. 
"""
function _meshgen(p::ParametricEntity;gridsize)
    # mesh the domain of p
    f    = parametrization(p)
    d    = domain(p)
    grid = CartesianMesh(d;gridsize)
    els = [ParametricElement(f,d) for d in ElementIterator(grid)]
    return els
end