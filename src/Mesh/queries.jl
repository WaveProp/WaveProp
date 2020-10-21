"""
    nearest_element_list(X,Y;[tol])

Return a vector of integers, where the `i` entry of the vector gives the index of the nearest element in `Y` to the *ith*-node.

An optional keywork argument `tol` can be passed so that only elements which are closer than `tol` are considered. If a node `x ∈ X` with index `i` has no element in `Y` closer than `tol`, the value -1 is stored indicating such a case.
"""
function nearest_element_list(X::NystromMesh,Y::NystromMesh; tol=0)
    n,m  = length(X),length(Y)
    list = fill(-1,n) # idxel of nearest element for each x ∈ X. -1 means there is not element in `Y` closer than `tol`
    xnodes = nodes(X)
    ynodes = nodes(Y)
    in2e   = _idx_nodes_to_elements(Y)
    if X == Y
        # special case (which arises often for integral operators) where X==Y.
        # No distance computation is needed then
        for i=1:n
            @assert length(in2e[i]) == 1
            list[i] = in2e[i] |>  first
        end
    else
        @notimplemented
    end
    return list
end

"""
    _idx_nodes_to_elements(q::NystromMesh)

For each node in `q`, return the indices of the elements to which it belongs.

Depending on the quadrature type, more efficient methods can be defined and overloaded if needed.
"""
function _idx_nodes_to_elements(mesh::NystromMesh)
    list = [Int[] for _ in 1:length(mesh)]
    for n in 1:length(mesh.el2quad)
        for i in mesh.el2quad[n]
            push!(list[i],n)
        end
    end
    return list
end