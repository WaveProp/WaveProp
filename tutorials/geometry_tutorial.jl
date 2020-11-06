# # Geometry 

# The `WaveProp.Geometry` module handles everything related to the geometrical
# descriptions of curves/surfaces. Its main interface is exposed through
# `AbstractElement`, which represents an `M`-dimensional geometical object
# embedded in ℝᴺ. In abstract terms, an *element* τ is a subset of ℝᴺ described
# by 
# * a reference element ``\tau \subset \mathbb{R}^N``
# * a *push-forward map* χ : τ̂ → τ ⊂ ℝᴺ.

# Concrete implementations of `AbstractElement`s are responsible for
# implementing a few crucial methods used to derive almost all functionality of
# an `el::AbstractElement` used throughout this package. In the [next section](#interface) we
# describe this interface.

# ## Interface<a id='interface'></a>

# The interface is composed of the following methods:
# - `domain(el)` : return the reference element `R::AbstractReferenceElement`
# - `el(u)` : evaluate the parametrization of the element at (parametric).
#   coordinate `u ∈ domain(el)`.
# - `jacobian(u)` : evaluate the jacobian of the parametrization at (parametric)
#   coordinate `u ∈ domain(el)`.

# ## Lagrangian elements<a id='lagrange-element'></a>

# An concrete implementation of `AbstractElement` is the parametric type
# `LagrangeElement{D,Np,N,T}`. The type parameters represent:
# - `D` : the domain of the element. Thus, `domain(el)` return an instance of
#   the singleton type `R`.
# - `Np` : the number of points/DOF which describe the element. For example, for
#   a flat triangle, `Np=3`, while for a flat quadrilateral element `Np=4`.
# - `N` : the dimension of the ambient/embedding space. 
# - `T` : the data type used in the presentation of `el`. Defaults to `Float64`.




