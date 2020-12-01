struct NystromMesh{N,T} <: AbstractMesh{N,T}
    elements::Dict{DataType,ElementIterator}
    # quadrature info
    qnodes::Vector{Point{N,T}}
    qweights::Vector{T}
    qnormals::Vector{Point{N,T}}
    el2qnodes::Dict{DataType,Matrix{Int}}
    etype2qrule::Dict{DataType,AbstractQuadratureRule}
end

# getters
qnodes(m::NystromMesh) = m.qnodes
qweights(m::NystromMesh) = m.qweights
qnormals(m::NystromMesh) = m.qnormals
el2qnodes(m::NystromMesh) = m.el2qnodes
el2qnodes(m::NystromMesh,E::DataType) = m.el2qnodes[E]

Base.length(m::NystromMesh) = length(qnodes(m))

function NystromMesh{N,T}() where {N,T}
    elements    = Dict{DataType,ElementIterator}()
    qnodes      = Vector{Point{N,T}}()
    qweights    = Vector{T}()
    qnormals    = Vector{Point{N,T}}()
    el2qnodes   = Dict{DataType,Matrix{Int}}()
    etype2qrule = Dict{DataType,AbstractQuadratureRule}()
    NystromMesh{N,T}(elements,qnodes,qweights,qnormals,el2qnodes,etype2qrule)
end    

function NystromMesh(mesh::AbstractMesh,e2qrule,compute_normal::Bool)
    N,T        = ambient_dimension(mesh), eltype(mesh)
    # initialize empty fields
    nmesh = NystromMesh{N,T}()
    # loop over element types, then call inner function. This allows for the
    # heavy lifting to be type-stable.
    for E in etypes(mesh)
        haskey(e2qrule,E) || error("no quadrature rule found for element of type $E")
        qrule = e2qrule[E]
        _build_nystrom_mesh!(nmesh,mesh,E,qrule,compute_normal)
    end
    return nmesh
end 

@noinline function _build_nystrom_mesh!(nmesh,mesh,E,qrule,compute_normal)
    iter        = ElementIterator{E}(mesh)
    nmesh.elements[E] = iter    
    x̂,ŵ         = qrule() #nodes and weights on reference element
    el2qnodes   = Int[]
    for el in iter
        # push forward to element    
        x = map(el,x̂)
        w = map(zip(x̂,ŵ)) do (x̂,ŵ)
            ŵ*measure(el,x̂)    
        end 
        # append quadrature information and keep track of which nodes belong to
        # which element 
        istart = length(qnodes(nmesh)) + 1
        append!(qweights(nmesh),w)
        append!(qnodes(nmesh),x)
        if compute_normal
            ν = map(x->normal(el,x),x̂)   
            append!(qnormals(nmesh),ν)
        end
        iend   = length(qnodes(nmesh))
        append!(el2qnodes,collect(istart:iend))
    end
    nmesh.el2qnodes[E] = reshape(el2qnodes,length(x̂),:)
    return nmesh
end

function NystromMesh(mesh::AbstractMesh;order=1,compute_normal::Bool=true)
    etype2qrule = Integration._qrule_for_mesh(mesh,order)
    NystromMesh(mesh,etype2qrule,compute_normal)
end    

function NystromMesh(mesh::AbstractMesh,Ω::Domain;kwargs...)
    submesh = SubMesh(mesh,Ω)
    NystromMesh(submesh;kwargs...)
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