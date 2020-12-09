"""
    meshgen(p::AbstractParametricBody;gridsize::Tuple) --> Ω,M

Return a `Domain` Ω and a `GenericMesh` M for the parametric body.
"""
function meshgen(p::AbstractParametricBody;gridsize)
    N       = ambient_dimension(p)    
    mesh    = GenericMesh{N,Float64}()
    dim,tag = Geometry._global_add_entity!(p) # add entity to global list to make it easily retriavable using only (dim,tag)
    ent     = ElementaryEntity(dim,tag)
    Ω       = Domain(ent) 
    mesh.ent2tags[ent] = Dict{DataType,Vector{Int}}() # no mesh is generated for domain, only its boundary, empty entry
    for bnd in boundary(p)
        # add elements for each boundary segment    
        bnd_dim,bnd_tag = Geometry._global_add_entity!(bnd)
        bnd_ent         = ElementaryEntity(bnd_dim,bnd_tag)
        els             = _meshgen(bnd;gridsize)
        E               = eltype(els)
        push!(etypes(mesh),E)
        vals            = get!(mesh.elements,E,Vector{E}())
        istart          = length(vals) + 1
        append!(vals,els)
        iend            = length(vals)
        # add mapping from entity to elements
        push!(boundary(ent),bnd_ent)
        mesh.ent2tags[bnd_ent] = Dict(E=>collect(istart:iend)) # add key
    end 
    unique!(etypes(mesh))   
    return Ω,mesh
end    

function _meshgen(p::ParametricEntity;gridsize)
    # mesh the domain of p
    f    = Geometry.parametrization(p)
    d    = domain(p)
    grid = CartesianMesh(d;gridsize)
    els = [Geometry.ParametricElement(f,d) for d in elements(grid)]
    return els
end