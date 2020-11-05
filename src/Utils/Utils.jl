"""
    Utils

Module containing various utility functions for `WaveProp`.
"""
module Utils

using DocStringExtensions
using StaticArrays

using WaveProp

export svector, notimplemented, @abstractmethod, assert_extension

"""
    svector(f,n)

Just like `ntuple`, but convert output to a `StaticVector`.
"""
svector(f,n) = ntuple(f,n) |> SVector

"""
    notimplemented()

Things which should probably be implemented at some point.
"""
function notimplemented()
    error("not (yet) implemented")
end 

"""
    abstractmethod

A method of an `abstract type` for which concrete subtypes are expected
to provide an implementation.
"""
function abstractmethod()
    error("this method needs to be implemented by the concrete subtype.")
end 

"""
    assert_extension(fname,ext,[msg])

Check that `fname` is of extension `ext`. Print the message `msg` as an assertion error otherwise.
"""
function assert_extension(fname::String,ext::String,msg="file extension must be $(ext)")
    r = Regex("$(ext)\$")    
    @assert occursin(r,fname) msg
end

function debug(mod)
    @eval ENV["JULIA_DEBUG"] = $(mod)
end

end # module