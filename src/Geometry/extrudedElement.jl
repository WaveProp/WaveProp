struct ExtrudedElement{B,N}
    base::B
    path::SVector{N,Float64}
end 

extrude(el::AbstractElement,d⃗) = ExtrudedElement(el,d⃗)

geometric_dimension(el::ExtrudedElement) = (el |> base |> geometric_dimension) + 1
ambient_dimension(el::ExtrudedElement)   = el |> base |> ambient_dimension

base(el::ExtrudedElement) = el.base
path(el::ExtrudedElement) = el.path

function (el::ExtrudedElement)(u)
    @assert length(u) == geometric_dimension(el)
    d⃗ = path(el)
    return base(el)(u[1:end-1]) + u[end]*d⃗
end    

function jacobian(el::ExtrudedElement,u)
    @assert length(u) == geometric_dimension(el)
    j = jacobian(base(el),u)
    d⃗ = path(el)
    return hcat(j,d⃗)    
end    

function boundary(el::ExtrudedElement)
    a,b = el |> base |> boundary
    d⃗   = path(el)
    τ   = base(el)
    return τ,line(b,b+d⃗),reverse(translate(τ,d⃗)),line(a+d⃗,a)
end    