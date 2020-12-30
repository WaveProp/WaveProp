"""
    struct NystromMesh{N,T} <: AbstractMesh{N,T}
    
A mesh data structure for solving boundary integral equation using Nyström
methods.

A `NystromMesh` can be constructed from a `mesh::AbstractMesh` and a dictionary
`etype2qrule` mapping element types in `mesh` (the keys) to an appropriate
quadrature rule for integration over elements of that type (the value).

The degrees of freedom of a `NystromMesh` are associated to its quadrature nodes
`qnodes::Vector{SVector{N,T}}`. 
"""
Base.@kwdef struct NystromMesh{N,T} <: AbstractMesh{N,T}
    elements::OrderedDict{DataType,Any} = OrderedDict{DataType,Any}()
    etype2qrule::OrderedDict{DataType,AbstractQuadratureRule}= OrderedDict{DataType,AbstractQuadratureRule}()
    # quadrature info
    qnodes::Vector{SVector{N,T}} = Vector{SVector{N,T}}()
    qweights::Vector{T} = Vector{T}()
    qnormals::Vector{SVector{N,T}} = Vector{SVector{N,T}}()
    # mappings from elements --> dofs and entities --> elements
    elt2dof::OrderedDict{DataType,Matrix{Int}} = OrderedDict{DataType,Matrix{Int}}()
    ent2elt::OrderedDict{AbstractEntity,OrderedDict{DataType,Vector{Int}}} = OrderedDict{AbstractEntity,OrderedDict{DataType,Vector{Int}}}()
end

function Base.append!(msh1::NystromMesh,msh2::NystromMesh)
    Ω = domain(msh2)    
    @assert ambient_dimension(msh1) == ambient_dimension(msh2)
    @assert eltype(msh1) == eltype(msh2)
    # extract relevant quadrature rules
    etype2dof = OrderedDict{DataType,Vector{Int}}()
    for E in etypes(msh2,Ω)
        if haskey(msh1.elements,E)
            etype2dof[E]        = msh1.elt2dof[E] |> vec
            @assert msh1.etype2qrule[E] == msh2.etype2qrule[E]
        else    
            etype2dof[E]        = Int[]    
            msh1.elements[E]    = E[]
            msh1.etype2qrule[E] = msh2.etype2qrule[E]
        end
    end
    # loop over entities and extract the information
    for ent in entities(msh2)
        haskey(msh1.ent2elt,ent) && continue # entity alreay meshed
        msh1.ent2elt[ent] = OrderedDict{DataType,Vector{Int}}()
        for (E,v) in ent2elt(msh2,ent)
            # add dofs information    
            msh2_dofs  = msh2.elt2dof[E][:,v]
            idx_start = length(msh1.qnodes) + 1 
            append!(msh1.qnodes,msh2.qnodes[msh2_dofs])
            append!(msh1.qweights,msh2.qweights[msh2_dofs])
            if !isempty(qnormals(msh2)) # otherwise assume it did not have normal stored
                append!(msh1.qnormals,msh2.qnormals[msh2_dofs])
            end     
            idx_end = length(msh1.qnodes) # end index of dofs added
            append!(etype2dof[E],collect(idx_start:idx_end))
            # append elements
            idx_start = length(msh1.elements[E]) + 1
            append!(msh1.elements[E],msh2.elements[E][v])
            idx_end   = length(msh1.elements[E])
            # and mapping ent -> elt
            msh1.ent2elt[ent][E] = collect(idx_start:idx_end)
        end
    end
    # finally do some rehshaping
    for (E,v) in etype2dof
        dof_per_el = size(msh2.elt2dof[E],1)
        msh1.elt2dof[E] = reshape(v,dof_per_el,:)
    end    
    return msh1
end    

# getters
qnodes(m::NystromMesh)   = m.qnodes
qweights(m::NystromMesh) = m.qweights
qnormals(m::NystromMesh) = m.qnormals
elements(m::NystromMesh) = m.elements

# various mappings
elt2dof(m::NystromMesh)  = m.elt2dof
elt2dof(m::NystromMesh,E::DataType) = m.elt2dof[E]
ent2elt(m::NystromMesh) = m.ent2elt
ent2elt(m::NystromMesh,ent::AbstractEntity) = m.ent2elt[ent]

Geometry.domain(mesh::NystromMesh)   = Domain(entities(mesh))
Geometry.entities(mesh::NystromMesh) = collect(keys(mesh.ent2elt))

function dom2elt(mesh::NystromMesh,domain::Domain)
    dict = OrderedDict{DataType,Vector{Int}}()
    for ent in entities(domain)
        other = ent2elt(mesh,ent)    
        mergewith!(append!,dict,other)
    end
    return dict
end    

function dom2dof(mesh::NystromMesh,domain::Domain)
    idxs = Int[]
    for ent in entities(domain)
        dict = ent2elt(mesh,ent)
        for (E,tags) in dict
            append!(idxs,view(elt2dof(mesh,E),:,tags))    
        end    
    end    
    return idxs
end    
dom2dof(mesh,ent::AbstractEntity) = dom2dof(mesh,Domain(ent))
ent2dof(mesh,ent::AbstractEntity) = dom2dof(mesh,ent)

Base.length(m::NystromMesh) = length(qnodes(m))

function etypes(m::NystromMesh,Ω::Domain=domain(m)) 
    ee = DataType[]
    for ent in entities(Ω)
        dict = ent2elt(m,ent)
        union!(ee, keys(dict))
    end    
    @assert allunique(ee)
    return ee
end

function NystromMesh(mesh::GenericMesh,Ω::Domain,e2qrule,compute_normal::Bool=true)
    N,T        = ambient_dimension(mesh), eltype(mesh)
    # initialize empty fields
    nys_mesh = NystromMesh{N,T}(;etype2qrule=e2qrule)
    # loop over entities
    for ent in entities(Ω) 
        dict    = OrderedDict{DataType,Vector{Int}}()    
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
        nys_mesh.ent2elt[ent] = dict # new entry 
    end   
    return nys_mesh 
end    

@noinline function _build_nystrom_mesh!(nys_mesh,submesh,E,qrule,compute_normal)
    els::Vector{E}         = get!(nys_mesh.elements,E,Vector{E}())
    el2qnodes   = get!(nys_mesh.elt2dof,E,Matrix{Int}(undef,0,0))
    el2qnodes   = vec(el2qnodes)
    x̂,ŵ         = qrule() #nodes and weights on reference element
    for el in ElementIterator(submesh,E)
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
    nys_mesh.elt2dof[E] = el2qnodes
    return nys_mesh
end

NystromMesh(submesh::SubMesh,args...) = NystromMesh(submesh.mesh,submesh.domain,args...)

function NystromMesh(mesh::SubMesh;order,compute_normal::Bool=true)
    etype2qrule = Mesh._qrule_for_mesh(mesh,order)
    NystromMesh(mesh,etype2qrule,compute_normal)
end    

function NystromMesh(mesh::GenericMesh,Ω::Domain;quad_rule,compute_normal::Bool=true)
    etype2qrule = OrderedDict(E=>quad_rule for E in etypes(view(mesh,Ω)))
    NystromMesh(mesh,Ω,etype2qrule,compute_normal)
end    

# construct a mesh as a restriction to a domain
function Base.getindex(msh2::NystromMesh,Ω::Domain)
    N,T      = ambient_dimension(msh2), eltype(msh2)
    # initialize empty fields
    msh1 = NystromMesh{N,T}()    
    # extract relevant quadrature rules
    etype2dof = OrderedDict{DataType,Vector{Int}}()
    for E in etypes(msh2,Ω)
        msh1.etype2qrule[E] = msh2.etype2qrule[E]
        etype2dof[E] = Int[]
        msh1.elements[E] = E[]
    end
    # loop over entities and extract the information
    for ent in Ω
        @assert ent in domain(msh2)
        msh1.ent2elt[ent] = OrderedDict{DataType,Vector{Int}}()
        for (E,v) in msh2.ent2elt[ent]
            # add dofs information    
            msh2_dofs  = msh2.elt2dof[E][:,v]
            idx_start = length(msh1.qnodes) + 1
            append!(msh1.qnodes,msh2.qnodes[msh2_dofs])
            append!(msh1.qweights,msh2.qweights[msh2_dofs])
            if !isempty(qnormals(msh2)) 
                append!(msh1.qnormals,msh2.qnormals[msh2_dofs])
            end     
            idx_end = length(msh1.qnodes)
            append!(etype2dof[E],collect(idx_start:idx_end))
            # append elements
            idx_start = length(msh1.elements[E]) + 1
            append!(msh1.elements[E],msh2.elements[E][v])
            idx_end   = length(msh1.elements[E])
            # and mapping ent -> elt
            msh1.ent2elt[ent][E] = collect(idx_start:idx_end)
        end
    end
    # finally do some rehshaping
    for (E,v) in etype2dof
        dof_per_el = size(msh2.elt2dof[E],1)
        msh1.elt2dof[E] = reshape(v,dof_per_el,:)
    end    
    return msh1
end   
Base.getindex(m::NystromMesh,ent::AbstractEntity) = m[Domain(ent)]

function _elt2dof(m::NystromMesh,Ω::Domain)
    out = OrderedDict{DataType,Matrix{Int}}()
    for ent in Ω
        dict = ent2elt(m,ent)
        for (E,v) in dict
            dofs = get!(out,E,Matrix{Int}(undef,0,0))  
            # make into vec to append  
            dofs = vec(dofs)
            dof_per_el = size(elt2dof(m,E),1)
            append!(dofs,elt2dof(m,E)[:,v])
            # convert back to matrix shape
            dofs = reshape(dofs,dof_per_el,:)
            out[E] = dofs
        end    
    end
    return out        
end    

"""
    near_interaction_list(pts,Y;dim,atol)

Given a target mesh `X` and a source mesh `Y`, return a dictionary with keys
given by element types of the source mesh `Y`. To each key, which represents an
element type, we associate a vector whose `i`-th entry encodes
information about the points in `X` which are close to a given element of
the key type. 
"""
function near_interaction_list(pts,Y;atol)
    dict = OrderedDict{DataType,Vector{Vector{Tuple{Int,Int}}}}()    
    for E in etypes(Y)
        ilist = _etype_near_interaction_list(pts,Y,E,atol)
        push!(dict,E=>ilist)
    end
    return dict
end    

function _etype_near_interaction_list(pts,Y,E,atol)
    ilist = Vector{Vector{Tuple{Int,Int}}}()
    e2n   = elt2dof(Y,E)
    npts,nel = size(e2n)
    for n in 1:nel
        ynodes = qnodes(Y)[e2n[:,n]]
        inear  = _near_interaction_list(pts,ynodes,atol)
        push!(ilist,inear)
    end 
    ilist       
end    

function _near_interaction_list(pts,ynodes,atol)
    ilist    = Tuple{Int,Int}[]
    for (i,x) in enumerate(pts)            
        d,j = findmin([norm(x-y) for y in ynodes])        
        if d ≤ atol
            push!(ilist,(i,j))    
        end
    end            
    return ilist
end   

elements(m::NystromMesh,E::DataType) = m.elements[E]