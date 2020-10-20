"""
    Utils

Various utilities functions for `WaveProp`.

Exports:
$EXPORTS    
"""
module Utils

using DocStringExtensions
using StaticArrays

using WaveProp

export svector, @notimplemented, assert_extension

svector(f,n) = ntuple(f,n) |> SVector

macro notimplemented()
    quote
        error("not (yet) implemented")
    end
end 

function assert_extension(fname::String,ext::String,msg="file extension must be $(ext)")
    r = Regex("$(ext)\$")    
    @assert occursin(r,fname) msg
end

function debug(mod)
    @eval ENV["JULIA_DEBUG"] = $(mod)
end

end # module