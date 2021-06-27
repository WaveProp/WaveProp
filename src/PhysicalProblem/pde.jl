abstract type AbstractPDE{N} end

ambient_dimension(::Type{<:AbstractPDE{N}}) where {N} = N
ambient_dimension(op::AbstractPDE) = ambient_dimension(typeof(op))

"""
    struct Laplace{N}

Laplace equation in `N` dimension: Δu = 0.
"""
struct Laplace{N} <: AbstractPDE{N} end

Laplace(;dim=3) = Laplace{dim}()

function Base.show(io::IO,pde::Laplace)
    print(io,"Δu = 0")
end

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

function Base.show(io::IO,pde::Helmholtz)
    # k = parameters(pde)
    print(io,"Δu + k u = 0")
end

parameters(pde::Helmholtz) = pde.k

default_kernel_eltype(::Helmholtz)  = ComplexF64
default_density_eltype(::Helmholtz) = ComplexF64

"""
    struct Elastostatic{N,T} <: AbstractPDE{N}

Elastostatic equation in `N` dimensions: μΔu + (μ+λ)∇(∇⋅u) = 0. Note that the
displacement u is a vector of length `N` since this is a vectorial problem.
"""
struct Elastostatic{N,T} <: AbstractPDE{N}
    μ::T
    λ::T
end
Elastostatic(;μ,λ,dim=3)               = Elastostatic{dim}(promote(μ,λ)...)
Elastostatic{N}(μ::T,λ::T) where {N,T} = Elastostatic{N,T}(μ,λ)

function Base.show(io::IO,pde::Elastostatic)
    # μ,λ = parameters(pde)
    print(io,"μΔu + (μ+λ)∇(∇⋅u) = 0")
end

parameters(pde::Elastostatic) = pde.μ, pde.λ

default_kernel_eltype(::Elastostatic{N}) where {N}  = SMatrix{N,N,Float64,N*N}
default_density_eltype(::Elastostatic{N}) where {N} = SVector{N,Float64}

"""
    Maxwell{T} <: AbstractPDE{3}

Normalized Maxwell's equation ∇ × ∇ × E - k² E = 0, where
k = ω √ϵμ.
"""
struct Maxwell{T} <: AbstractPDE{3}
    k::T
end

Maxwell(;dim=3,k::T) where {T}        = Maxwell{T}(k)

parameters(pde::Maxwell) = pde.k

function Base.show(io::IO,pde::Maxwell)
    # k = parameters(pde)
    print(io,"∇ × ∇ × E - k² E = 0")
end

default_kernel_eltype(::Maxwell)   = SMatrix{3,3,ComplexF64,9}
default_density_eltype(::Maxwell)  = SVector{3,ComplexF64}

"""
    MaxwellCFIE{T} <: AbstractPDE{3}

Normalized Maxwell's equation ∇ × ∇ × E - k² E = 0, where
k = ω √ϵμ.
"""
struct MaxwellCFIE{T} <: AbstractPDE{3}
    k::T   # wavenumber
end
MaxwellCFIE(;k) = MaxwellCFIE{ComplexF64}(k)
parameters(pde::MaxwellCFIE) = pde.k

function Base.show(io::IO,pde::MaxwellCFIE)
    # k = parameters(pde)
    print(io,"∇ × ∇ × E - k² E = 0")
end

default_kernel_eltype(::MaxwellCFIE)   = SMatrix{3,3,ComplexF64,9}
default_density_eltype(::MaxwellCFIE)  = SVector{3,ComplexF64}
