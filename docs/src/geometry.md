## Geometry

For the purpose of this package, a *geometrical element* ``\tau``
is fully described by a map ``\chi_\tau`` taking a *reference shape* ``\hat{\tau} \subset \mathbb{R}^n`` and mapping it into ... 

### Interface

The `Geometry` module defines and exports various types of `AbstractElements`. Given an `el <: AbstractElement`, the following methods are important:
- `el(x)` : the push-forward of a point `x` defined on the reference element ``\hat{\tau}``
- `jacobian(el,x)` : the jacobian matrix at the point `x`.


```@autodocs
Modules = [WaveProp.Geometry]
```




