function assemble_ldim(iop::IntegralOperator;compress=Matrix)
    X = target_surface(iop)
    Y = source_surface(iop)
    out        = compress(iop)
    correction = assemble_ldim_correction(iop)
    axpy!(1,correction,out) # out <-- out + correction
    return out
end    

function assemble_ldim_correction(iop::IntegralOperator;compress=Matrix)
    X    = target_surface(iop)    
    Y    = source_surface(iop)    
    pde  = iop.kernel.op
    dict_near = near_interaction_list(qnodes(X),Y;atol=1e-1)
    T = eltype(iop)
    # vectors to store nonzero indices (i,j) and value
    Is = Int[]
    Js = Int[]
    Vs = T[]
    for (E,list_near) in dict_near
        _assemble_ldim_correction!(Is,Js,Vs,iop,E,list_near)
    end    
    Sp = sparse(Is,Js,Vs,size(iop)...)
    return Sp
end  

@noinline function _assemble_ldim_correction!(Is,Js,Vs,iop,E,list_near)
    X    = target_surface(iop)    
    Y    = source_surface(iop)    
    T    = eltype(iop)  
    k    = iop.kernel  
    op   = iop.kernel.op
    el2qnodes = Y.el2qnodes[E]
    els  = elements(Y,E)
    @assert length(els) == length(list_near)
    for nel in 1:length(els)
        el    = els[nel]    
        tags  = el2qnodes[:,nel]     
        # interpolation nodes on source element    
        yi = Y.qnodes[tags]     # interpolation/quadrature nodes 
        νi = Y.qnormals[tags]
        # now loop over target points in X which need correction when integrated
        # over el
        for (iglob,jloc) in list_near[nel]
            x   = X.qnodes[iglob]
            y₀  = yi[jloc]
            d   = norm(x-y₀)
            location = d == 0 ? :onsurface : :offsurface
            w   = singular_weights_ldim(k,el,yi,νi,x,location)
            v   = transpose(w) - iop[iglob,tags] # remove the "naive" approximation
            @assert length(w) == length(yi)
            append!(Is,fill(iglob,length(w)))
            append!(Js,tags)
            append!(Vs,v)
        end
    end    
end

"""
    singular_weights_ldim(k,τ,yi,x)

Given a kernel `k`, generate quadrature weights `w` such that
```math
\\int_\\tau k(x,y) \\varphi(y) ds\\_y \\approx \\sum_q \\varphi(y_q) w_q
```
"""
@inline function singular_weights_ldim(k::AbstractKernel,el,yi,νi,x,location=:onsurface)    
    T = return_type(k)
    h  = 0.1
    nb = 3*length(yi)
    op =  k.op
    # interpolation surface
    Γ = auxiliary_domain_ldim(el,-h)
    # find appropriate source points xₛ given the points xq
    γ₀B, γ₁B    = basis_ldim(op,el,h,nb)
    G   = SingleLayerKernel(op)
    dG  = DoubleLayerKernel(op)
    r   = zeros(T,nb)
    for n in 1:nb
        r[n] = sum(Γ) do b
            integrate(ReferenceLine) do v
                y  = b(v)    
                νy = normal(b,v)    
                μ  = measure(b,v)
                (G(x,y)*γ₁B[n](y,νy) - dG(x,y,νy)*γ₀B[n](y))*μ
            end    
        end    
        if location == :onsurface
            r[n] = r[n] - γ₀B[n](x)/2    
        end
    end    
    T    = return_type(G)
    M    = Matrix{T}(undef,2*length(yi),nb)
    for j in 1:nb
        for i in 1:length(yi)
            M[i,j]            = γ₀B[j](yi[i])
            M[length(yi)+i,j] = γ₁B[j](yi[i],νi[i])
        end
    end    
    α,β = combined_field_coefficients(k)
    D  = [diagm(α*ones(length(yi))) ; diagm(β*ones(length(yi)))]
    w  = (transpose(r)/M)*D
    return w
end 

function basis_ldim(pde,el,h,nb)
    ν  = normal(el,0.5)    
    xc = el(0.5) + h/2*ν
    r  = 10*h
    xs = _circle_sources(nsources=nb, radius=r, center=xc)
    G   = SingleLayerKernel(pde)
    dG  = DoubleLayerKernel(pde)
    γ₀B = map(xs) do x
        y -> G(x,y)
    end
    γ₁B = map(xs) do x
        (y,ν) -> dG(x,y,ν)
    end
    return γ₀B, γ₁B    
end    

function auxiliary_domain_ldim(el::AbstractElement,h)
    @assert domain(el) == ReferenceLine()
    a,b = el(0),el(1) # the endpoints
    ν   = normal(el,0.5) 
    Γ₁  = LagrangeLine(b,b+h*ν)
    Γ₂  = LagrangeLine(b+h*ν,a+h*ν)
    Γ₃  = LagrangeLine(a+h*ν,a)
    return Γ₁,Γ₂,Γ₃
end    