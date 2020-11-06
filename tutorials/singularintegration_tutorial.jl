# # Singular Integrals

# This example shows how to use the `SingularIntegration` module for computing
# quadrature rules for functions with point singularities (such as those
# appearing in the numerical discretization of boundary integral equations).

# ## Change of variables

# The first set of *tricks* revolve a round a simple change of variables. We
# focus first on the one-dimensional case, where we wish to integrate 
```math
    \int_0^1 f(x) dx,
```
# and where the function ``f`` can have an integrable singularity at 

using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.SingularIntegration
using QuadGK

f        = (x) -> x==0 ? 0.0 : log(abs(x))*cos(x)
I,_        = quadgk(f,0,1,rtol=1e-16)

rows = GaussLegendre.([5,10,20,40,80])
cols = [identity,IMT{1,2}(),Kress{8}(), Window{1,1,7}()]    
ee   = []
for qstd in rows
    for shandler in cols
        q        = SingularQuadratureRule(qstd,shandler)
        Is       = integrate(f,q)
        er       = abs(I-Is)
        push!(ee,er)
    end    
end

using NamedArrays
ee = reshape(ee,length(cols),length(rows)) |> transpose |> NamedArray
setnames!(ee,string.(rows),1)
setnames!(ee,string.(cols),2)
setdimnames!(ee,["Base quadrature","Singularity handler"])
show(ee)

# ## 
# using Plots
# qstd = GaussLegendre(10)
# x̂,ŵ  = qstd()
# cols = [IMT{1,2}(),Kress{8}(), Window{0.5,1,7}()]    
# fig = plot()
# for shandler in cols
#     phi       = shandler.(x̂)
#     phip      = [jacobian(shandler,x)[1] for x in x̂]
#     plot!(fig,x̂,f.(phi) .* phip,label=string(shandler),m=:x)
# end 
# display(fig)   

