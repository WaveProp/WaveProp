"""
    abstract type KernelType
    
Trait used classify a given kernel function `k`. 
"""
abstract type KernelType end

struct UnknownKernelType <: KernelType end
struct CombinedField{a,b} <: KernelType end
struct DerivativeCombinedField{a,b} <: KernelType end

combined_field_coefficients(::CombinedField{a,b}) where {a,b} = (a,b)
combined_field_coefficients(::DerivativeCombinedField{a,b}) where {a,b} = (a,b)

const SingleLayer         = CombinedField{1,0}
const DoubleLayer         = CombinedField{0,1}
const AdjointDoubleLayer  = DerivativeCombinedField{1,0}
const HyperSingular       = DerivativeCombinedField{0,1}

kernel_type(::Type) = UnknownKernelType()
kernel_type(x)      = kerneltype(typeof(x))
