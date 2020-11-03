"""
    Utils

Various utilities functions for `WaveProp`.

"""
module Utils

using DocStringExtensions
using StaticArrays

using WaveProp

export svector, @notimplemented, assert_extension

"""
    svector(f,n)

Just like [`ntuple`](@ref), but convert output to a `StaticVector`.
"""
svector(f,n) = ntuple(f,n) |> SVector
# FIXME: how to cross-reference the function `ntuple` from base in the docstring?

"""
    @notimplemented

Things which should probably be implemented.
"""
macro notimplemented()
    quote
        error("not (yet) implemented")
    end
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

error_green_formula(SL,DL,γ₀u,γ₁u,u,σ)                      = σ*u + SL*γ₁u  - DL*γ₀u
error_derivative_green_formula(ADL,H,γ₀u,γ₁u,un,σ)          = σ*un + ADL*γ₁u - H*γ₀u
error_interior_green_identity(SL,DL,γ₀u,γ₁u)                = error_green_formula(SL,DL,γ₀u,γ₁u,γ₀u,-1/2)
error_interior_derivative_green_identity(ADL,H,γ₀u,γ₁u)     = error_derivative_green_formula(ADL,H,γ₀u,γ₁u,γ₁u,-1/2)
error_exterior_green_identity(SL,DL,γ₀u,γ₁u)                = error_green_formula(SL,DL,γ₀u,γ₁u,γ₀u,1/2)
error_exterior_derivative_green_identity(ADL,H,γ₀u,γ₁u)     = error_derivative_green_formula(ADL,H,γ₀u,γ₁u,γ₁u,1/2)

end # module