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
    matrix_to_blockmatrix,
    blockmatrix_to_matrix,
    blockvector_to_vector,
    vector_to_blockvector,
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
    blockmatrix_to_matrix(A::Matrix{B}) where {B<:SMatrix}

Convert a `Matrix{B}`, where `B<:SMatrix`, to the equivalent `Matrix{T}`, where `T = eltype(B)`
"""
function blockmatrix_to_matrix(A::Matrix{B})  where B <: SMatrix
    T = eltype(B) 
    sblock = size(B)
    ss     = size(A).*sblock # matrix size when viewed as matrix over T
    Afull = Matrix{T}(undef,ss)
    for i=1:ss[1], j=1:ss[2]
        bi, ind_i = divrem(i-1,sblock[1]) .+ (1,1)
        bj, ind_j = divrem(j-1,sblock[2]) .+ (1,1)
        Afull[i,j] = A[bi,bj][ind_i,ind_j]
    end
    return Afull
end

"""
    blockvector_to_vector(A::Vector{B}) where {B<:SVector}

Convert a `Vector{B}`, where `B<:SVector`, to the equivalent `Vector{T}`, where `T = eltype(B)`
"""
function blockvector_to_vector(A::Vector{B})  where B <: SVector
    T = eltype(B)     
    reinterpret(T,A) |> collect
end

"""
    matrix_to_blockmatrix(A::Matrix,B)

Convert a `Matrix{T}` to a `Matrix{B}`, where `B<:Type{SMatrix}`. The element
type of `B` must match that of `A`, and the size of `A` must be divisible by the
size of `B` along each dimension. 
"""
function matrix_to_blockmatrix(A::Matrix,B::Type{<:SMatrix})
    @assert eltype(A) == eltype(B)
    @assert sum(size(A) .% size(B)) == 0 "block size $(size(B)) not compatible with size of A=$(size(A))"
    sblock = size(B)
    nblock = div.(size(A),sblock)
    Ablock = Matrix{B}(undef,nblock)
    for i in 1:nblock[1]
        istart = (i-1)*sblock[1] + 1
        iend = i*sblock[1]
        for j in 1:nblock[2]
            jstart = (j-1)*sblock[2] + 1
            jend   = j*sblock[2]
            Ablock[i,j] = A[istart:iend,jstart:jend]
        end
    end
    return Ablock
end

"""
    vector_to_blockvector(A::Vector,B)

Convert a `Vector{T}` to a `Vector{B}`, where `B<:Type{SVector}`. The element
type of `B` must match that of `A`, and the size of `A` must be divisible by the
size of `B` along each dimension. 
"""
function vector_to_blockvector(A::Vector,B::Type{<:SVector})
    @assert eltype(A) == eltype(B)
    @assert sum(size(A) .% size(B)) == 0 "block size $(size(B)) not compatible with size of A=$(size(A))"
    T = eltype(B)     
    reinterpret(B,A) |> collect
end

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