"""
    abstract type AbstractMesh{N,T}
    
An abstract mesh structure in dimension `N` with primite data of type `T`. 

By default `T=Float64`.

An instance of an `AbstractMesh` is expected to implement the following methods:

#TODO: decide on the interface
"""
abstract type AbstractMesh{N,T} end

"""
    struct GenericMesh{N,T,E}

A simple data structure representing a generic mesh in an ambient space of dimension `N`, with data of type `T`, and composed of elements of type `E`. 
"""
struct GenericMesh{N,T,E}
    vtx::Vector{Point{N,T}}
    # element types
    etypes::E
    # mapping from element type to indices of vtx in each element
    el2vtx::Vector{Matrix{Int}}
end


