# Geometry module

The `WaveProp.Geometry` module is of fundamental importance as it defines the
various geometrical objects used to perform the numerical simulations. At a
high-level, the usual workflow is:
```math
    \fbox{Geometry} \rightarrow \fbox{Mesh} \rightarrow \fbox{PDE} \rightarrow \fbox{Solution}
```
Here the `Geometry` handles everything related to defining geometrical entities,
as well as **geometrical elements**. 

