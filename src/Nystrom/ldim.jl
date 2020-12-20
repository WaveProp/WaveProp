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
            ν₀  = νi[jloc]
            r   = x-y₀
            d   = norm(r)
            side = dot(r,ν₀)
            location = d == 0 ? :onsurface : side > 0 ? :outside : :inside
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
@inline function singular_weights_ldim(k::AbstractKernel,el,yi,νi,x,side=:onsurface)    
    # pre-compute things which do not depend on x
    op  = k.op
    T   = return_type(k)
    K0  = SingleLayerKernel(op)
    K1  = DoubleLayerKernel(op)
    Γ,σ,γ₀B,γ₁B,F = _auxiliary_quantities_ldim(k,el,yi,νi,side)
    # compute things which depend on target point x
    nb  = length(γ₀B)
    I   = zeros(T,nb)
    _integrate_basis_ldim!(I,Γ,K0,K1,γ₀B,γ₁B,x,σ) # modify r
    α,β = combined_field_coefficients(k)
    D   = [diagm(α*ones(length(yi))) ; diagm(β*ones(length(yi)))]
    w   = ((transpose(I)/F.R)*adjoint(F.Q))*D
    return w
end 

function _auxiliary_quantities_ldim(k::AbstractKernel,el,yi,νi,side)
    T = return_type(k)
    h  = 0.1
    nb = 3*length(yi)
    op =  k.op
    # interpolation surface
    if side == :onsurface
        Γ = auxiliary_domain_ldim(el,-h)
        σ = -0.5
    elseif side == :outside
        Γ = auxiliary_domain_ldim(el,-h)
        σ = 0.0
    elseif side == :inside
        Γ = auxiliary_domain_ldim(el,h)
        σ = 0.0
    end    
    # find a "basis" for the interpolation
    γ₀B, γ₁B    = basis_ldim(op,el,h,nb)
    # build and factor interpolation matrix
    M    = Matrix{T}(undef,2*length(yi),nb)
    for j in 1:nb
        for i in 1:length(yi)
            M[i,j]            = γ₀B[j](yi[i])
            M[length(yi)+i,j] = γ₁B[j](yi[i],νi[i])
        end
    end   
    F   = qr(M)
    return Γ,σ,γ₀B,γ₁B,F
end    

function _integrate_basis_ldim!(r,Γ,K0,K1,γ₀B,γ₁B,x,σ)
    nb = length(r)
    for n in 1:nb
        r[n] = sum(Γ) do b
            integrate(ReferenceLine) do v
                y  = b(v)    
                νy = normal(b,v)    
                μ  = measure(b,v)
                (K0(x,y)*γ₁B[n](y,νy) - K1(x,y,νy)*γ₀B[n](y))*μ
            end    
        end    
        r[n] = r[n] + σ*γ₀B[n](x)
    end    
    return r
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