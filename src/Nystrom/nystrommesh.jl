Base.@kwdef struct NystromMesh{N,T} <: AbstractMesh{N,T}
    elements::Dict{DataType,Any} = Dict{DataType,Any}()
    # quadrature info
    qnodes::Vector{Point{N,T}} = Vector{Point{N,T}}()
    qweights::Vector{T} = Vector{T}()
    qnormals::Vector{Point{N,T}} = Vector{Point{N,T}}()
    el2qnodes::Dict{DataType,Matrix{Int}} = Dict{DataType,Matrix{Int}}()
    etype2qrule::Dict{DataType,AbstractQuadratureRule}= Dict{DataType,AbstractQuadratureRule}()
    ent2tags::Dict{ElementaryEntity,Dict{DataType,Vector{Int}}} = Dict{ElementaryEntity,Dict{DataType,Vector{Int}}}()
end

# getters
qnodes(m::NystromMesh) = m.qnodes
qweights(m::NystromMesh) = m.qweights
qnormals(m::NystromMesh) = m.qnormals
el2qnodes(m::NystromMesh) = m.el2qnodes
el2qnodes(m::NystromMesh,E::DataType) = m.el2qnodes[E]

Geometry.entities(mesh::NystromMesh) = collect(keys(mesh.ent2tags))

function dof(mesh::NystromMesh,domain::Domain)
    idxs = Int[]
    for ent in entities(domain)
        dict = mesh.ent2tags[ent]
        for (E,tags) in dict
            append!(idxs,view(mesh.el2qnodes[E],:,tags))    
        end    
    end    
    return idxs
end    
dof(mesh,ent::ElementaryEntity) = dof(mesh,Domain(ent))

Base.length(m::NystromMesh) = length(qnodes(m))

function NystromMesh(mesh::GenericMesh,Ω::Domain,e2qrule,compute_normal::Bool=true)
    N,T        = ambient_dimension(mesh), eltype(mesh)
    # initialize empty fields
    nys_mesh = NystromMesh{N,T}(;etype2qrule=e2qrule)
    # loop over entities
    for ent in entities(Ω) 
        dict    = Dict{DataType,Vector{Int}}()    
        submesh = view(mesh,ent)
        for E in etypes(submesh)
            haskey(e2qrule,E) || error("no quadrature rule found for element of type $E")
            qrule = e2qrule[E]
            if haskey(nys_mesh.elements,E)
                istart = length(nys_mesh.elements[E])+1
            else
                istart = 1
            end
            _build_nystrom_mesh!(nys_mesh,submesh,E,qrule,compute_normal)
            iend   = length(nys_mesh.elements[E])       
            dict[E] = collect(istart:iend)
        end  
        nys_mesh.ent2tags[ent] = dict # new entry 
    end   
    return nys_mesh 
end    

@noinline function _build_nystrom_mesh!(nys_mesh,submesh,E,qrule,compute_normal)
    els::Vector{E}         = get!(nys_mesh.elements,E,Vector{E}())
    el2qnodes   = get!(nys_mesh.el2qnodes,E,Matrix{Int}(undef,0,0))
    el2qnodes   = vec(el2qnodes)
    x̂,ŵ         = qrule() #nodes and weights on reference element
    for el in elements(submesh,E)
        # compute quadrature on element
        x = map(el,x̂)
        w = map(zip(x̂,ŵ)) do (x̂,ŵ)
            ŵ*measure(el,x̂)    
        end 
        # add element
        push!(els,el)
        # append quadrature information
        istart = length(qnodes(nys_mesh)) + 1
        append!(qweights(nys_mesh),w)
        append!(qnodes(nys_mesh),x)
        if compute_normal
            ν = map(x->normal(el,x),x̂)   
            append!(qnormals(nys_mesh),ν)
        end
        iend   = length(qnodes(nys_mesh))
        # add information to recover qnodes from index of el
        append!(el2qnodes,collect(istart:iend))
    end
    el2qnodes = reshape(el2qnodes,length(x̂),:)
    nys_mesh.el2qnodes[E] = el2qnodes
    return nys_mesh
end

NystromMesh(submesh::SubMesh,args...) = NystromMesh(submesh.mesh,submesh.domain,args...)

function NystromMesh(mesh::SubMesh;order=1,compute_normal::Bool=true)
    etype2qrule = Integration._qrule_for_mesh(mesh,order)
    NystromMesh(mesh,etype2qrule,compute_normal)
end    

# interface methods
etypes(m::NystromMesh) = keys(m.elements) |> collect

"""
    near_interaction_list(X,Y;dim,atol)

Given a target mesh `X` and a source mesh `Y`, return a dictionary with keys given by element types of the source mesh `Y`. To each key, which represents an element type, we associate a vector encoding whose `i`-th entry encodes information about the points in `X` which are close to the a given element of the key type. 
"""
function near_interaction_list(X,Y;dim,atol)
    dict = Dict{DataType,Vector{Vector{Tuple{Int,Int}}}}()    
    for E in etypes(Y)
        geometric_dimension(E) == dim || continue
        ilist = _etype_near_interaction_list(X,Y,E,atol)
        push!(dict,E=>ilist)
    end
    return dict
end    

function _etype_near_interaction_list(X,Y,E,atol)
    ilist = Vector{Vector{Tuple{Int,Int}}}()
    e2n   = el2qnodes(Y,E)
    npts,nel = size(e2n)
    for n in 1:nel
        ynodes = qnodes(Y)[e2n[:,n]]
        inear  = _near_interaction_list(X,ynodes,atol)
        push!(ilist,inear)
    end 
    ilist       
end    

function _near_interaction_list(X,ynodes,atol)
    ilist    = Tuple{Int,Int}[]
    for (i,x) in enumerate(qnodes(X))            
        d,j = findmin([norm(x-y) for y in ynodes])        
        if d ≤ atol
            push!(ilist,(i,j))    
        end
    end            
    return ilist
end   