using Artifacts

const gmsh_path = artifact"gmsh4.8.4"

const gmsh_dirs = readdir(gmsh_path,join=true)

@assert length(gmsh_dirs)==1 "there should be only one directory for the untared gmsh file. Got $dirs"

include(joinpath(first(gmsh_dirs),"lib","gmsh.jl"))

export gmsh # gmsh module

"""
macro gmsh(ex)

Initialize `gmsh` through `gmsh.initilize(), execute `ex`, the close `gmsh`
through `gmsh.finalize()`.
"""
macro gmsh(ex)
    return quote
        gmsh.initialize()
        try
            _gmsh_set_verbosity(0)
            $(esc(ex))
        finally
            # make sure we finalize gmsh if something goes wrong
            gmsh.finalize()
        end
    end
end

"""
    read_geo(fname::String;dim=3)

Read a `.geo` file and generate a [`GenericMesh`](@ref) together with a
[`Domain`](@ref) for it using `gmsh`.
"""
function read_geo(fname;dim=3,h=nothing,order=nothing)
    assert_extension(fname, ".geo")
    Ω,M = @gmsh begin
        gmsh.open(fname)
        h     === nothing || _gmsh_set_meshsize(h)
        order === nothing || _gmsh_set_meshorder(order)
        gmsh.model.mesh.generate(dim)
        Ω = _initialize_domain(dim)
        M = _initialize_mesh(Ω)
        if dim == 2
            M = convert_to_2d(M)
        end
        Ω,M
    end
    return Ω,M
end

"""
    read_msh(fname::String)

Similar to [`read_geo`](@ref), but the mesh is simply read from the input file
instead of generated.
"""
function read_msh(fname;dim=3)
    assert_extension(fname, ".msh")
    Ω,M = @gmsh begin
        gmsh.open(fname)
        Ω = _initialize_domain(dim)
        M = _initialize_mesh(Ω)
        Ω, M
    end
    return Ω,M
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
        if haskey(ENTITIES,(dim,tag))
            ent = ENTITIES[(dim,tag)]
        else
            ent = ElementaryEntity(dim, tag)
            _fill_entity_boundary!(ent)
        end
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
    combine  = true # FIXME: what should we use here?
    oriented = false
    dim_tags = gmsh.model.getBoundary((ent.dim, ent.tag),combine,oriented)
    for (d, t) in dim_tags
        if haskey(ENTITIES,(d,t))
            bnd = ENTITIES[(d,t)]
        else
            bnd = ElementaryEntity(d, t)
            _fill_entity_boundary!(bnd)
        end
        push!(ent.boundary, bnd)
    end
    return ent
end

"""
    _initialize_mesh(Ω::Domain)

Performs all the GMSH API calls to extract the information necessary to
construct the mesh.
"""
function _initialize_mesh(Ω::Domain)
    tags, coord, _ = gmsh.model.mesh.getNodes()
    nodes = reinterpret(SVector{3,Float64}, coord) |> collect
    # map gmsh type tags to actual internal types
    etypes = [_type_tag_to_etype(e) for e in gmsh.model.mesh.getElementTypes()]
    # Recursively populating the dictionaries
    elements = OrderedDict{DataType,Matrix{Int}}()
    ent2tags = OrderedDict{AbstractEntity,OrderedDict{DataType,Vector{Int}}}()
    elements, ent2tag = _domain_to_mesh!(elements, ent2tags, Ω)
    return GenericMesh{3,Float64}(;nodes, elements, ent2tags)
end

"""
    _domain_to_mesh!(elements, ent2tag, Ω::Domain)

Recursively populating the dictionaries `elements` and `ent2tag`.
"""
function _domain_to_mesh!(elements, ent2tag, Ω::Domain)
    isempty(Ω) && (return elements, ent2tag)
    for ω in Ω
        _ent_to_mesh!(elements, ent2tag, ω)
    end
    Γ = skeleton(Ω)
    _domain_to_mesh!(elements, ent2tag, Γ)
end

"""
    _ent_to_mesh!(elements, ent2tag, ω::ElementaryEntity)

For each element type used to mesh `ω`:
- push into `elements::OrderedDict` the pair `etype=>ntags`;
- push into `ent2tag::OrderedDict` the pair `etype=>etags`;

where:
- `etype::DataType` determines the type of the element (see
    [`_type_tag_to_etype`](@ref));
- `ntags::Matrix{Int}` gives the indices of the nodes defining those
    elements;
- `etags::Vector{Int}` gives the indices of those elements in `elements`.
"""
function _ent_to_mesh!(elements, ent2tag, ω::ElementaryEntity)
    ω in keys(ent2tag) && (return elements, ent2tag)
    etypes_to_etags = OrderedDict{DataType,Vector{Int}}()
    # Loop on GMSH element types (integer)
    type_tags, _, ntagss = gmsh.model.mesh.getElements(key(ω)...)
    for (type_tag, ntags) in zip(type_tags, ntagss)
        _, _, _, Np, _ = gmsh.model.mesh.getElementProperties(type_tag)
        ntags = reshape(ntags, Int(Np), :)
        etype = _type_tag_to_etype(type_tag)
        if etype in keys(elements)
            etag = size(elements[etype], 2) .+ collect(1:size(ntags,2))
            ntags = hcat(elements[etype], ntags)
        else
            etag = collect(1:size(ntags, 2))
        end
        push!(elements, etype => ntags)
        push!(etypes_to_etags, etype => etag)
    end
    push!(ent2tag, ω => etypes_to_etags)
    return elements, ent2tag
end

"""
    gmsh_disk(;rx=0.5,ry=0.5,center=(0,0,0)) -> Ω, M

Use `gmsh` API to generate a disk and return a [`Domain`](@ref) `Ω` and a
[`GenericMesh`](@ref) of the disk.
"""
function gmsh_disk(;rx=0.5,ry=0.5,center=(0.,0.,0.),dim=2,h=min(rx,ry)/10,order=1)
    @gmsh begin
        _gmsh_set_meshsize(h)
        _gmsh_set_meshorder(order)
        gmsh.model.occ.addDisk(center...,rx,ry)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(dim)
        Ω = _initialize_domain(2)
        M = _initialize_mesh(Ω)
        M = convert_to_2d(M)
        return Ω,M
    end
end

"""
    gmsh_sphere(;radius=0.5,center=(0,0,0),dim=3,h=radius/10,order=1) -> Ω, M

Use `gmsh` API to generate a sphere and return `Ω::Domain` and `M::GenericMesh`.
Only entities of dimension `≤ dim` are meshed.
"""
function gmsh_sphere(;radius=0.5,center=(0., 0., 0.),dim=3,h=radius/10,order=1,recombine=false)
    Ω,M = @gmsh begin
        _gmsh_set_meshsize(h)
        _gmsh_set_meshorder(order)
        gmsh.model.occ.addSphere(center..., radius)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(dim)
        recombine && gmsh.model.mesh.recombine()
        Ω = _initialize_domain(3)
        M = _initialize_mesh(Ω)
        Ω, M
    end
    return Ω,M
end

"""
    gmsh_box(;origin=(0,0,0),widths=(0,0,0)) -> Ω, M

Use `gmsh` API to generate an axis aligned box. Return `Ω::Domain` and `M::GenericMesh`.
"""
function gmsh_box(;origin=(0., 0., 0.),widths=(1., 1., 1.),h=0.1,dim=3)
    Ω, M = @gmsh begin
        _gmsh_set_meshsize(h)
        gmsh.model.occ.addBox(origin..., widths...)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(dim)
        Ω = _initialize_domain(3)
        M = _initialize_mesh(Ω)
        Ω,M
    end
    return Ω, M
end

"""
    gmsh_rectangle(;origin,widths)
"""
function gmsh_rectangle(;origin=(0.,0.,0.),dx=1,dy=1,dim=2,h=0.1)
    Ω,M = @gmsh begin
        _gmsh_set_meshsize(h)
        gmsh.model.occ.addRectangle(origin...,dx,dy)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(dim)
        Ω = _initialize_domain(dim)
        M = _initialize_mesh(Ω)
        if dim == 2
            M = convert_to_2d(M)
        end
        Ω,M
    end
    return Ω,M
end

function _gmsh_set_meshsize(hmax, hmin=hmax)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", hmin)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", hmax)
end

function _gmsh_set_meshorder(order)
    gmsh.option.setNumber("Mesh.ElementOrder", order)
end

function _gmsh_set_verbosity(i)
    gmsh.option.setNumber("General.Verbosity",i)
end

function _gmsh_summary(model)
    gmsh.model.setCurrent(model)
    @printf("List of entities in model `%s`: \n", model)
    @printf("|%10s|%10s|%10s|\n","name","dimension","tag")
    ents = gmsh.model.getEntities()
    # pgroups = gmsh.model.getPhysicalGroups()
    for ent in ents
        name = gmsh.model.getEntityName(ent...)
        dim, tag = ent
        @printf("|%10s|%10d|%10d|\n", name, dim, tag)
    end
    println()
end

function _gmsh_summary()
    models = gmsh.model.list()
    for model in models
        _gmsh_summary(model)
    end
end

"""
    _type_tag_to_etype(tag)

Mapping of `gmsh` element types, encoded as an integer, to the internal
equivalent of those. This function assumes `gmsh` has been initilized.
"""
function _type_tag_to_etype(tag)
    T = SVector{3,Float64} # point type
    name,dim,order,num_nodes,ref_nodes,num_primary_nodes  = gmsh.model.mesh.getElementProperties(tag)
    num_nodes = Int(num_nodes) #convert to Int64
    if occursin("Point",name)
        etype = LagrangePoint{1,T}
    elseif occursin("Line",name)
    etype = LagrangeLine{num_nodes,T}
    elseif occursin("Triangle",name)
        etype = LagrangeTriangle{num_nodes,T}
    elseif occursin("Quadrilateral",name)
        etype = LagrangeSquare{num_nodes,T}
    elseif occursin("Tetrahedron",name)
        etype = LagrangeTetrahedron{num_nodes,T}
    else
        error("gmsh element of family $name does not an internal equivalent")
    end
    return etype
end

"""
    _etype_to_type_tag(etype)

Mapping of internal element types, to the integer tag of `gmsh` elements. This
function assumes `gmsh` has been initialized.
"""
function _etype_to_type_tag(el::LagrangeElement)
    etype = typeof(el)
    tag = 1
    while true
        E   = _type_tag_to_etype(tag)
        E === etype && (return tag)
        tag = tag + 1
    end
end

# use gmsh API to extract the reference nodes
# TODO: we should just code these internally at some point instead of relying on
# gmsh.
@generated function reference_nodes(::LagrangeElement{D,K}) where {D,K}
    if D == ReferenceLine
        family_name = "line"
        @gmsh begin
            tag = gmsh.model.mesh.getElementType(family_name,K)
            name,dim,order,num_nodes,ref_nodes,num_primary_nodes  = gmsh.model.mesh.getElementProperties(tag)
        end
        num_nodes = Int(num_nodes)
        ref_nodes = SVector{num_nodes}(ref_nodes)
        return :($ref_nodes)
    elseif D == ReferenceTriangle
        family_name = "triangle"
        @gmsh begin
            tag = gmsh.model.mesh.getElementType(family_name,K)
            name,dim,order,num_nodes,ref_nodes,num_primary_nodes  = gmsh.model.mesh.getElementProperties(tag)
        end
        num_nodes = Int(num_nodes)
        ref_nodes = reshape(ref_nodes,2,num_nodes)
        # convert nodes to 2d static vectors
        sref_nodes = svector(num_nodes) do i
            SVector{2}(ref_nodes[:,i])
        end
        return :($sref_nodes)
    else
        notimplemented()
    end
end
