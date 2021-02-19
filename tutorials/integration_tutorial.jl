using WaveProp.Geometry
using WaveProp.Integration
using Test

# # Integration

# The `WaveProp.Integration` module implements the various quadrature rules and
# integration routines required throughout the package. The interface of
# the `Integration` module is given through `AbstractQuadratureRule{D}`, where the
# type parameter `D` is a subtype of `AbstractReferenceShape`.

# For a given a domain `D`, an $N$-point quadrature rule over `D` is a set of nodes $\{x_i\}$
# and weights $\{w_i\}$, $1 \leq i \leq N$, which can be used to integrate a function over `D`. That is
# ```math
# \int_D f(x) dx \approx \sum_{i=1}^{N} f(x_i)w_i
# ```

# There are several concrete subtypes of `AbstractQuadratureRule`. To see a list
# of available implementations you can type:
subtypes(Integration.AbstractQuadratureRule)

# Different quadrature rules may require different input arguments in order to
# be instantiate. For example, the 1d quadrature rules of `GaussLegendre` type
# require as argument `n`, the number of quadrature nodes. You can instantiate a
# 10-point `GaussLegendre` rule as
q = GaussLegendre(10)

# To get the nodes and weights of `q` you can simply type:
x̂,ŵ = q()

# Similarly, you can `integrate` a callable object usign the syntax
I = integrate(x->cos(x)*exp(x),q)

# By default, calling an `AbstractQuadrature` object with no arguments returns
# the underlying quadrature rule for integration over `domain(q)`. The basic
# strategy for integration over an arbitrarily shaped element `el` is simple:
# first, develop a quadrature rules on reference shape `D == domain(el)`, and
# then *lift* the quadrature nodes `xᵢ` and weights `wᵢ` to the element `el`
# using its parametrization. Notice that this requires `el` to implement `el(u)`
# and `jacobian(el,u)` for `u ∈ D`. 

# We demonstrate how to do that in a simple element: a line segment
# connecting `a=(0,0)` to `b=(2,2)`:
a,b = (0,0), (2,2)
el = line(a,b)
x,w = q(el)

# It is easy to check that `sum(w) == √8`, i.e. the length of the line
@test sum(w) ≈  √8

