"""
    Utils

Various utilities functions for `WaveProp`.

Exports:
$EXPORTS    
"""
module Utils

using DocStringExtensions
using StaticArrays

export svector

svector(f,n) = ntuple(f,n) |> SVector

end