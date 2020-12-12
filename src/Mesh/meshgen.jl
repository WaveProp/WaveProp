"""
    meshgen(p::AbstractParametricBody;gridsize::Tuple) --> Ω,M

Return a `Domain` Ω and a `GenericMesh` M for the parametric body.
"""
function meshgen!(mesh::GenericMesh,Ω::Domain,p::AbstractParametricBody;gridsize)
    dim,tag = Geometry._global_add_entity!(p) # add entity to global list to make it easily retriavable using only (dim,tag)
    ent     = ElementaryEntity(dim,tag)
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
    push!(Ω.entities,ent)
    return Ω,mesh
end    

function meshgen(p::AbstractParametricBody;gridsize=floatmax())
    N       = ambient_dimension(p)    
    mesh    = GenericMesh{N,Float64}()
    Ω       = Domain()
    meshgen!(mesh,Ω,p;gridsize)
end   

function meshgen(ps::Vector{<:AbstractParametricBody};gridsize=floatmax())
    N       = ambient_dimension(first(ps))
    @assert all(p->ambient_dimension(p)==N,ps)
    mesh    = GenericMesh{N,Float64}()
    Ω       = Domain()    
    for p in ps
        meshgen!(mesh,Ω,p;gridsize)
    end
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