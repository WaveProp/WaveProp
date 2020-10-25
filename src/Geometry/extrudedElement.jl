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

function boundary(el::ExtrudedElement)
    a,b = el |> base |> boundary
    d⃗   = path(el)
    return el,line(b,b+d⃗),reverse(translate(el,d⃗)),line(a+d⃗,a)
end    