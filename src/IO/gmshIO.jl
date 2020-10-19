"""
    read_geo(fname::String;dim=3)

Read a `.geo` file and generate a domain [`Ω::Domain`](@ref) of dimension `d`.
"""
function read_geo(fname, dim=3)
    assert_extension(fname, ".geo")    
    occursin(r".geo$", fname) || error(msg)
    gmsh.initialize()
    gmsh.open(fname)    
    Ω = _initialize_domain(dim)
    gmsh.finalize()
    return Ω
end

"""
    _initialize_domain(d)

Construct a `Domain` from the current `gmsh` model, starting from entities of dimension `d`. 

This is a helper function, and should not be called by itself since it assumes that `gmsh` has been initialized.
"""
function _initialize_domain(dim)
    Ω = Domain() # Create empty domain
    dim_tags = gmsh.model.getEntities(dim)
    for (_, tag) in dim_tags
        ent = ElementaryEntity(dim, tag)
        _fill_entity_boundary!(ent)
        push!(Ω, ent)
    end        
    return Ω
end

"""
    _fill_entity_boundary!

Use the `gmsh` API to add the boundary of an `ElementaryEntity`. 

This is a helper function, and should not be called by itself since it assumes that `gmsh` has been initialized.
"""
function _fill_entity_boundary!(ent)
    combine = true # FIXME: what should we use here?
    dim_tags = gmsh.model.getBoundary((ent.dim, ent.tag), combine)
    for (d, t) in dim_tags
        bnd = ElementaryEntity(d, t)
        _fill_entity_boundary!(bnd)
        push!(ent.boundary, bnd)
    end    
    return ent
end    

"""
    read_msh(fname::String)

Read a `.msh` file and return a domain [`Ω::Domain`](@ref) together with a mesh [`M::GenericMesh`](@ref).
"""
function read_msh(fname;dim=3)
    assert_extension(fname, ".msh")    
    gmsh.initialize()    
    gmsh.open(fname)
    Ω = _initialize_domain(dim)
    M = _initialize_mesh(Ω)
    gmsh.finalize()    
    return Ω, M
end    

function _initialize_mesh(Ω::Domain)
    tags, coord, _ = gmsh.model.mesh.getNodes()
    pts    = reinterpret(SVector{3,Float64}, coord) |> collect
    etypes = gmsh.model.mesh.getElementTypes()
    el2vtx = Vector{Matrix{Int}}(undef, length(etypes))
    for (i, etype) in enumerate(etypes)
        etags, ntags = gmsh.model.mesh.getElementsByType(etype)
        _, _, _, Np, rest = gmsh.model.mesh.getElementProperties(etype)
        ntags     = reshape(ntags, Int(Np), :)
        el2vtx[i] = ntags
    end    
    dict    = Dict{ElementaryEntity,Dict{Int32,Vector{Int}}}()
    ent2tag = _domain_to_mesh!(dict, Ω)
    return GenericMesh(pts, etypes, el2vtx, ent2tag)
end    

"""
    _domain_to_mesh!(dict,Ω)

For each entity in Ω, push into `Dict` that entity as a key and a `Dict{Int,Vector{Int}}` as value. The value gives an element type, as well as the tags associated to that entity and element type.
"""
function _domain_to_mesh!(dict, Ω)
    isempty(Ω) && (return dict)
    for ent in entities(Ω)
        _ent_to_mesh!(dict, ent)
    end    
    Γ = skeleton(Ω)
    _domain_to_mesh!(dict, Γ)
    return dict
end

"""
    _ent_to_mesh!(dict,ent)

For each element type used to mesh `ent`, push into `dict::Dictionary` the pair `etype=>tags`, 
where `etype::Int32` determines the type of the element (see [`type_tag_to_etype`](@ref)), and
`tags::Vector{Int}` gives the tags of those elements.
"""
function _ent_to_mesh!(dict, ent)
    etypes, etags, ntags = gmsh.model.mesh.getElements(ent.dim, ent.tag)
    etypes_to_etags = Dict(etypes[i] => etags[i] for i in 1:length(etypes))
    push!(dict, ent => etypes_to_etags)
end    

"""
    gmsh_sphere(;radius=0.5,center=(0,0,0)) -> Ω, M

Use `gmsh` API to generate a sphere and return `Ω::Domain` and `M::GenericMesh`.
"""
function gmsh_sphere(;radius=0.5,center=(0.,0.,0.))
    gmsh.initialize()
    gmsh.model.occ.addSphere(center...,radius)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate()
    Ω = _initialize_domain(3)
    M = _initialize_mesh(Ω)
    gmsh.finalize()
    return Ω,M
end    

"""
    gmsh_box(;origin=(0,0,0),widths=(0,0,0)) -> Ω, M

Use `gmsh` API to generate an axis aligned box. Return `Ω::Domain` and `M::GenericMesh`.
"""
function gmsh_box(;origin=(0.,0.,0.),widths=(1.,1.,1.))
    gmsh.initialize()
    gmsh.model.occ.addBox(origin...,widths...)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate()
    Ω = _initialize_domain(3)
    M = _initialize_mesh(Ω)
    gmsh.finalize()
    return Ω,M
end    