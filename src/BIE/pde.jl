abstract type AbstractPDE{N} end

ambient_dimension(::Type{<:AbstractPDE{N}}) where {N} = N
ambient_dimension(op) = ambient_dimension(typeof(op))

"""
    struct Laplace{N}
    
Laplace equation in `N` dimension: Δu = 0.
"""
struct Laplace{N} <: AbstractPDE{N} end    

Laplace(;dim=3) = Laplace{dim}()

getname(::Laplace) = "Laplace"

default_kernel_eltype(::Laplace)  = Float64
default_density_eltype(::Laplace) = Float64

"""
    struct Helmholtz{N,T}
    
Helmholtz equation in `N` dimensions: Δu + k²u = 0.
"""
struct Helmholtz{N,K} <: AbstractPDE{N}
    k::K
end    

Helmholtz(;k,dim=3) = Helmholtz{dim,typeof(k)}(k)

getname(::Helmholtz) = "Helmholtz"

default_kernel_eltype(::Helmholtz)  = ComplexF64
default_density_eltype(::Helmholtz) = ComplexF64