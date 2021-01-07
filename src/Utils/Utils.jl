"""
    Utils

Module containing various utility functions for `WaveProp`.
"""
module Utils

using DocStringExtensions
using StaticArrays

using WaveProp

export 
    svector, 
    notimplemented, 
    abstractmethod, 
    assert_extension,
    assert_concrete_type

"""
    svector(f,n)

Just like [`Base.ntuple`](https://docs.julialang.org/en/v1/base/base/#Base.ntuple), but convert output to an `SVector`.
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
function abstractmethod(T)
    error("this method needs to be implemented by the concrete subtype $T.")
end 

"""
    assert_extension(fname,ext,[msg])

Check that `fname` is of extension `ext`. Print the message `msg` as an assertion error otherwise.
"""
function assert_extension(fname::String,ext::String,msg="file extension must be $(ext)")
    r = Regex("$(ext)\$")    
    @assert occursin(r,fname) msg
end

function assert_concrete_type(T::DataType)
    isconcretetype(T) || throw(ConcreteInferenceError(T)) 
end

function debug(mod="WaveProp")
    @eval ENV["JULIA_DEBUG"] = $(mod)
end

struct ConcreteInferenceError <: Exception
    T
end
Base.showerror(io::IO, e::ConcreteInferenceError) = print(io, "unable to infer concrete type from function signature: T=$(e.T)" )

end # module