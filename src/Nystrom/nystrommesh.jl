struct NystromMesh{N,T} <: AbstractMesh{N,T}
    elements::Dict{DataType,ElementIterator}
    # quadrature info
    qnodes::Vector{Point{N,T}}
    qweights::Vector{T}
    qnormals::Vector{Point{N,T}}
    el2qnodes::Dict{DataType,Matrix{Int}}
    etype2qrule::Dict{DataType,AbstractQuadratureRule}
end

function NystromMesh{N,T}() where {N,T}
    elements    = Dict{DataType,ElementIterator}()
    qnodes      = Vector{Point{N,T}}()
    qweights    = Vector{T}()
    qnormals    = Vector{Point{N,T}}()
    el2qnodes   = Dict{DataType,Matrix{Int}}()
    etype2qrule = Dict{DataType,AbstractQuadratureRule}()
    NystromMesh{N,T}(elements,qnodes,qweights,qnormals,el2qnodes,etype2qrule)
end    

qnodes(m::NystromMesh) = m.qnodes
qweights(m::NystromMesh) = m.qweights
qnormals(m::NystromMesh) = m.qnormals

el2qnodes(m::NystromMesh) = m.el2qnodes
el2qnodes(m::NystromMesh,E::DataType) = m.el2qnodes[E]

function NystromMesh(mesh::AbstractMesh,e2qrule)
    N,T        = ambient_dimension(mesh), eltype(mesh)
    # initialize empty fields
    nmesh = NystromMesh{N,T}()
    # loop over element types, then call inner function. This allows for the
    # heavy lifting to be type-stable.
    for E in etypes(mesh)
        haskey(e2qrule,E) || error("no quadrature rule found for element of type $E")
        qrule = e2qrule[E]
        _build_nystrom_mesh!(nmesh,mesh,E,qrule)
    end
    return nmesh
end 

@noinline function _build_nystrom_mesh!(nmesh,mesh,E,qrule)
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
        ν = map(x->normal(el,x),x̂)   
        # append quadrature information and keep track of which nodes belong to
        # which element 
        istart = length(qnodes(nmesh)) + 1
        append!(qweights(nmesh),w)
        append!(qnodes(nmesh),x)
        append!(qnormals(nmesh),ν)
        iend   = length(qnodes(nmesh))
        append!(el2qnodes,collect(istart:iend))
    end
    nmesh.el2qnodes[E] = reshape(el2qnodes,length(x̂),:)
    return nmesh
end

function NystromMesh(mesh::AbstractMesh;order=1)
    etype2qrule = Integration._qrule_for_mesh(mesh,order)
    NystromMesh(mesh,etype2qrule)
end    

function NystromMesh(mesh::AbstractMesh,Ω::Domain;kwargs...)
    submesh = SubMesh(mesh,Ω)
    NystromMesh(submesh;kwargs...)
end    

# interface methods
etypes(m::NystromMesh) = keys(m.elements) |> collect