"""
    Utils

Various utilities functions for `WaveProp`.

Exports:
$EXPORTS    
"""
module Utils

using DocStringExtensions
using StaticArrays

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

function debug(flag)
    if flag
        @eval ENV["JULIA_DEBUG"] = "WaveProp"
    else
        @eval ENV["JULIA_DEBUG"] = ""    
    end
end

debug(true)

end # module