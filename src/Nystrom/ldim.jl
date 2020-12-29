Base.@kwdef struct LocalDIMParameters
    compression=Matrix
    atol_near::Float64=0.1
    rtol_near::Float64=5
    nb_multiplier::Int=3
    r_multiplier::Float64=10
    aux_surf_multiplier::Float64=5
end    

function assemble_ldim(iop::IntegralOperator,p::LocalDIMParameters=LocalDIMParameters())
    X = target_surface(iop)
    Y = source_surface(iop)
    out        = p.compression(iop)
    correction = assemble_ldim_correction(iop,p)
    axpy!(1,correction,out) # out <-- out + correction
end    

function assemble_ldim_correction(iop::IntegralOperator,p::LocalDIMParameters)
    X    = target_surface(iop)    
    Y    = source_surface(iop)    
    pde  = iop.kernel.op
    dict_near = near_interaction_list(qnodes(X),Y;atol=p.atol_near)
    T = eltype(iop)
    # vectors to store nonzero indices (i,j) and value
    Is = Int[]
    Js = Int[]
    Vs = T[]
    for (E,list_near) in dict_near
        _assemble_ldim_correction!(Is,Js,Vs,iop,E,list_near,p)
    end    
    Sp = sparse(Is,Js,Vs,size(iop)...)
    return Sp
end  

@noinline function _assemble_ldim_correction!(Is,Js,Vs,iop,E,list_near,p::LocalDIMParameters)
    X    = target_surface(iop)    
    Y    = source_surface(iop)    
    T    = eltype(iop)  
    k    = iop.kernel  
    α,β = combined_field_coefficients(k)    
    op   = iop.kernel.op
    K0 = SingleLayerKernel(op)
    K1 = DoubleLayerKernel(op)
    el2qnodes = elt2dof(Y,E)
    els  = elements(Y,E)
    @assert length(els) == length(list_near)
    ni      = size(el2qnodes,1)
    nb      = p.nb_multiplier*ni
    M       = Matrix{T}(undef,2*ni,nb)
    I       = zeros(T,nb)
    for nel in 1:length(els)
        el    = els[nel]    
        tags  = @view el2qnodes[:,nel]     
        # interpolation nodes on source element    
        yi = @view Y.qnodes[tags]   
        νi = @view Y.qnormals[tags]
        wi = @view Y.qweights[tags]
        # basis
        l       = sum(wi)
        h       = sum(wi)/length(wi) # typical space of a grid
        γ₀B,γ₁B = basis_ldim(op,el,p.r_multiplier*h,nb)
        # interpolation matrix
        interpolation_matrix!(M,yi,νi,γ₀B,γ₁B)
        # F   = qr!(M)
        D   = [diagm(α*ones(length(yi))) ; diagm(β*ones(length(yi)))]
        # now loop over target points in X which need correction when integrated
        # over el
        for (iglob,jloc) in list_near[nel]
            x   = X.qnodes[iglob]
            y₀  = yi[jloc]
            ν₀  = νi[jloc]
            r   = x-y₀
            d   = norm(r)
            d > p.rtol_near*h && continue
            side = dot(r,ν₀)
            location = d == 0 ? :onsurface : side > 0 ? :outside : :inside
            if location == :onsurface
                Γ = auxiliary_domain_ldim(el,-p.aux_surf_multiplier*h)
                σ = -0.5
            elseif location == :outside
                Γ = auxiliary_domain_ldim(el,-h)
                σ = 0.0
            elseif location == :inside
                Γ = auxiliary_domain_ldim(el,h)
                σ = 0.0
            end   
            _integrate_basis_ldim!(I,Γ,K0,K1,γ₀B,γ₁B,x,σ) # modify r
            w   = (transpose(I)/M)*D
            v   = transpose(w) - iop[iglob,tags] # remove the "naive" approximation
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
    h = 0.1
    # pre-compute things which do not depend on x
    op  = k.op
    T   = return_type(k)
    K0  = SingleLayerKernel(op)
    K1  = DoubleLayerKernel(op)
    nb  = 3*length(yi)
    γ₀B,γ₁B = basis_ldim(op,el,10*h,nb)
    M       = Matrix{T}(undef,2*length(yi),nb)
    interpolation_matrix!(M,yi,νi,γ₀B,γ₁B)
    F   = qr(M)
    # compute things which depend on target point x
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
    I   = zeros(T,nb)
    _integrate_basis_ldim!(I,Γ,K0,K1,γ₀B,γ₁B,x,σ) # modify r
    α,β = combined_field_coefficients(k)
    D   = [diagm(α*ones(length(yi))) ; diagm(β*ones(length(yi)))]
    w   = ((transpose(I)/F.R)*adjoint(F.Q))*D
    return w
end 

function _singular_weights_ldim_kernel(el,K0,K1,γ₀B,γ₁B,x,side)
    nb = length(γ₀B)    
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
    I   = zeros(T,nb)
    _integrate_basis_ldim!(I,Γ,K0,K1,γ₀B,γ₁B,x,σ) # modify r
    α,β = combined_field_coefficients(k)
    D   = [diagm(α*ones(length(yi))) ; diagm(β*ones(length(yi)))]
    w   = ((transpose(I)/F.R)*adjoint(F.Q))*D
    return w
end

function interpolation_matrix!(M,yi,νi,γ₀B,γ₁B)
    nb = length(γ₀B)
    for j in 1:nb
        for i in 1:length(yi)
            M[i,j]            = γ₀B[j](yi[i])
            M[length(yi)+i,j] = γ₁B[j](yi[i],νi[i])
        end
    end   
    return M
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

function _integrate_basis_ldim_static!(r,Γ,K0,K1,γ₀B,γ₁B,x,σ,q,qs)
    nb = length(r)
    for n in 1:nb
        x̂,ŵ = q()    
        b = Γ[2]    
        r[n] = sum(zip(x̂,ŵ)) do (x̂,ŵ)
                y  = b(x̂)    
                νy = normal(b,x̂)    
                μ  = measure(b,x̂)
                (K0(x,y)*γ₁B[n](y,νy) - K1(x,y,νy)*γ₀B[n](y))*μ*ŵ
        end    
        b = Γ[1]    
        x̂,ŵ = qs()
        r[n] += sum(zip(x̂,ŵ)) do (x̂,ŵ)
                y  = b(x̂)    
                νy = normal(b,x̂)    
                μ  = measure(b,x̂)
                (K0(x,y)*γ₁B[n](y,νy) - K1(x,y,νy)*γ₀B[n](y))*μ*ŵ
        end    
        b = Γ[3]    
        x̂,ŵ = qs()
        x̂ = -x̂ .+ 1 |> reverse
        ŵ = reverse(ŵ)
        r[n] += sum(zip(x̂,ŵ)) do (x̂,ŵ)
                y  = b(x̂)    
                νy = normal(b,x̂)    
                μ  = measure(b,x̂)
                (K0(x,y)*γ₁B[n](y,νy) - K1(x,y,νy)*γ₀B[n](y))*μ*ŵ
        end    
        r[n] = r[n] + σ*γ₀B[n](x)
    end    
    return r
end

function basis_ldim(pde,el,r,nb)
    ν  = normal(el,0.5)    
    xc = el(0.5) # midpoint of the element
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