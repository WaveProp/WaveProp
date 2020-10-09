"""
$(TYPEDEF)

It includes the following fields:
$(FIELDS)

`etypes` stores the different type tags (e.g. `0 => node`, 
`1 => 2-point line`, `3 => 3-point triangle`, and so on). 

`el2vtx::Dict{Int,Matrix{Int}}` stores, for each element type in `etypes`, 
a matrix of indices describing all elements of that type. The `j` column of the matrix stores the indices in `vtx` which comprise the `j` element.
"""
struct GenericMesh{N,T}
    vtx::Vector{Point{N,T}}
    # element types
    etypes::Vector{Int}
    el2vtx::Dict{Int,Matrix{Int}}
end