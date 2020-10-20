struct NystromMesh{N,T} 
    elements::Vector{AbstractElement}
    quadrature::GenericQuadrature{N,T}
end    

function NystromMesh(M::GenericMesh,qrule;dim=ambient_dimension(M)-1)
    N,T = ambient_dimension(M), eltype(M)
    Q = GenericQuadrature{N,T}()
    els = Vector{AbstractElement}()
    for (i,etype) in enumerate(etypes(M))
        geometric_dimension(etype) == dim || continue 
        tags     = M.el2vtx[i]  # Np × Nel matrix
        Np, Nel  = size(tags)   # num of pts per element, num. of elements
        for n in 1:Nel
            el_vtx = M.vtx[tags[:,n]] # get the coordinates of nodes in this element
            el  = etype(el_vtx)       # construct the element
            push!(els,el)
            @assert domain(el) == domain(qrule)
            x̂,ŵ = qrule()
            x,w = push_forward_map(el,x̂,ŵ)
            append!(Q.nodes,x)
            append!(Q.weights,w)
            n⃗ = map(u->normal(el,u),x̂) 
            append!(Q.normals,n⃗)
        end    
    end 
    return NystromMesh{N,T}(els,Q)
end    

nodes(m::NystromMesh)    = m.quadrature.nodes
normals(m::NystromMesh)  = m.quadrature.normals
weights(m::NystromMesh)  = m.quadrature.weights
elements(m::NystromMesh) = m.elements


# function NystromMesh(mesh::GenericMesh{N,T}) where {N,T}
#     qrule = Gauss{ReferenceTriangle,1}()
#     quad  = quadgen(mesh,qrule,dim=2,need_normal=true)
# end 