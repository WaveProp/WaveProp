abstract type AbstractQuadratureRule{N} end

struct Trapezoidal{N} <: AbstractQuadratureRule{N} end

"""
    trapezoidal(n::Integer) -> (x,w)

Compute nodes `x` and weights `w` of the trapezoidal rule for integrating a function over `[-1,1]`.

Note that for analytic periodic functions, the trapezoidal rule converges exponentially fast. 
"""
function trapezoidal(n::Integer)
    h = 2/(n-1)
    x = [-1.0 + k*h for k in 0:n-1]
    w = [h for k in 1:n]
    w[1]   = h/2
    w[end] = h/2
    return x,w
end    

struct Fejer1{N} <: AbstractQuadratureRule{N} end

"""
    fejer1(n::Integer) -> (x,w)

Compute nodes `x` and weights `w` of Fejer's first quadrature rule for integrating a function over `[-1,1]`.
"""
function fejer1(n::Integer)
    theta = [(2j-1)*π/(2n) for j=1:n]
    x = -cos.(theta)
    w = zero(x)
    for j in 1:n
        tmp = 0.0
        for l in 1:floor(n/2)
            tmp += 1/(4*l^2-1) *cos(2*l*theta[j])
        end
        w[j] = 2/n * (1 - 2*tmp)
    end
    return x,w
end

struct Gauss{N} <: AbstractQuadratureRule{N} end

struct DoubleExponential{N} <: AbstractQuadratureRule{N} end

"""
    double_exponential(n::Integer) -> (x,w)

Compute nodes `x` and weights `w` of the double expoenential quadrature rule for integrating a function over `[-1,1]`.
"""
function double_exponential(n::Integer,h=_estimate_h(n))
    @assert isodd(n) "double exponential requires odd number of points."
    x̂ = (t) -> tanh(π/2*sinh(t))
    ŵ = (t) -> h*π/2*cosh(t)/cosh(π/2*sinh(t))^2
    krange = -(n-1)÷2 : (n-1)÷2
    x = [x̂(k*h) for k in krange]
    w = [ŵ(k*h) for k in krange]
    return x,w
end

# FIXME: find out a how to estimate the stepsize h in the double exponential
# formula. I think what you want is to find `h` s.t. n*h = constant, where the
# constant should depend on the desired precision. 
function _estimate_h(n;tol=1e-16)
    x̂ = (t) -> tanh(π/2*sinh(t))
    ŵ = (t) -> h*π/2*cosh(t)/cosh(π/2*sinh(t))^2
    h = 2^(-8)
    while true
        if ŵ(n*h) < tol || 1-x̂(n*h) < tol
            return h    
        else
            h = 2*h
        end    
    end
end    

