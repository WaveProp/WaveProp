function assemble_alpert(iop::IntegralOperator;compress=Matrix)
    X = target_surface(iop)
    Y = source_surface(iop)
    @assert X === Y "Alpert method should only be used for computing integral operators"
    N   = length(X) # number of dof
    out        = compress(iop)
    correction = assemble_alpert_correction!(out,iop)
    axpy!(1,correction,out) # out <-- out + correction
    return 
end    

function assemble_alpert_correction!(out,iop)
    X = target_surface(iop)
    Y = source_surface(iop)
    @assert X === Y "Martensen-Kussmaul method should only be used for computing integral operators"
    N   = length(X) # number of dof
    out = zero(iop) # output matrix to be filled
    for target_ent in entities(X)
        for source_ent in entities(Y)    
            if source_ent == target_ent
                _assemble_alpert_correction!(out,iop,target_ent)    
            else 
                # nothing special to do here, just fill up the matrix using the
                # regular quadrature    
                target_dof  = dof(X,target_ent)
                source_dof  = dof(Y,source_ent)        
                @. out[target_dof,source_dof] = iop[target_dof,source_dof]
            end        
        end    
    end   
    return out     
end    

# Assemble the self-interaction part for a given entity
function _assemble_alpert_correction!(out,iop,ent::AbstractEntity)
    X         = target_surface(iop)
    etype2tag = X.ent2elt[ent]
    # loop over individual elements composing ent. Note that because we must
    # handle meshes composed of elements of possiby different types, we must
    # loop over element types, and then over the elements of that type
    # (contained in the entity)
    for (target_E,target_tags) in etype2tag, target_t in target_tags
        for (source_E,source_tags) in etype2tag, source_t in source_tags
            if source_E === target_E && target_t === source_t # same element
                dof = X.elt2dof[source_E][:,source_t] # dof for filling out in global matrix
                el   = X.elements[source_E][source_t]                 # retrieve the element (i.e. access to its parametrization)     
                _assemble_alpert_correction!(out,iop,el,dof)
            else
                target_dof  = X.elt2dof[target_E][:,target_t] # dof for filling out in global matrix
                source_dof  = X.elt2dof[source_E][:,source_t] # dof for filling out in global matrix
                out[target_dof,source_dof] = iop[target_dof,source_dof]
            end        
        end
    end 
    return out   
end    

function _assemble_alpert_correction!(out,iop,el::AbstractElement,idxs)
    X     = target_surface(iop)
    qrule = X.etype2qrule[typeof(el)] # quadrature used on elements of type el
    x̂,ŵ   = qrule()                   # nodes on reference segment [0,1]
    # compute all quantities needed for mk rule, which includes higher order
    # derivatives of the parametrization. 
    N     = length(idxs)
    dx    =  [derivative(el,u)/(2π)     for u in x̂]
    d2x   =  [derivative2(el,u)/(4*π^2) for u in x̂]
    τ     = [norm(τ) for τ in dx]
    Δs    = 2π / N
    psi   = -0.577215664901532
    x     = qnodes(X)
    ν     = qnormals(X)
    w     = qweights(X)
    R     = nystrom_weights(N)
    k     = iop.kernel.op.k
    K     = kernel(iop)
    φ     = mk_split(K)
    for (iloc,iglob) in enumerate(idxs)
        for (jloc,jglob) in enumerate(idxs)
            lSin = log(4*(sin(Δs*(iloc-jloc)/2))^2) # singular part factored
            if kernel_type(iop) === SingleLayer()
                K1   = φ(x[iglob],x[jglob])*w[jglob]
                if iloc != jloc
                    K2   = K(x[iglob],x[jglob])*w[jglob]-K1*lSin # what is left after factoring out lSin and φ
                else    
                    K2   = (im/4+psi/2/pi-1/4/pi*log(k^2/4*τ[jloc]^2))*w[jglob] + 2*log(w[jglob]/(τ[jloc]*Δs))*K1
                end    
            elseif kernel_type(iop) == DoubleLayer()
                K1   = φ(x[iglob],x[jglob],ν[jglob])*w[jglob]
                if iloc != jloc
                    K2   = K(x[iglob],x[jglob],ν[jglob])*w[jglob]-K1*lSin
                else    
                    K2   = -1/(4*pi)*(dx[jloc][1]*d2x[jloc][2]-dx[jloc][2]*d2x[jloc][1])/τ[jloc]^3*w[jglob] + 2*log(w[jglob]/(τ[jloc]*Δs))*K1
                end    
            end
            out[iglob,jglob] = (R[iloc,jloc]*K1 + K2)
        end
    end   
    return out
end    

# analytical splitting of kernel in 2d
function mk_split(SL::SingleLayerKernel{T,S}) where {T,S<:Helmholtz}
    k = SL.op.k
    ϕ = (x,y) -> begin
        d = norm(x-y)
        (-1/(4*pi))*besselj0(k*d) |> T
    end
    return ϕ
end    

function mk_split(DL::DoubleLayerKernel{T,S}) where {T,S<:Helmholtz}
    k = DL.op.k
    ϕ = (x,y,ν) -> begin
        x == y && (return zero(T)) # analitic limit taken by hand to avoid 0/0
        r = x-y
        d = norm(r)
        (-k/(4*pi))*besselj1(k*d)/d*dot(r,ν) |> T
    end
    return ϕ
end    

function nystrom_weights(M)
    @assert M%2 == 0  "number of points `M` must be even"
    Md2 = div(M,2)    
    R   = zeros(M,M)
    for p=1:M            
        tp = pi/Md2*(p-1)
        for j=1:M
            tj = pi/Md2*(j-1)
            R[p,j]=-2*sum((ones(Md2-1)./collect(1:Md2-1)).*cos.(collect(1:Md2-1)*(tp-tj)))-pi/(Md2^2)*cos(Md2*(tp-tj))
        end            
    end
    return R
end

function  _alpert_params(order)
    if order == 2    
        a = 1
        chi_p = [1.591549430918953e-01]
        w_p = [5.000000000000000e-01]
    elseif order == 6        
        a = 3
        chi_p= [4.004884194926570e-03;
                7.745655373336686e-02;
                3.972849993523248e-01;
                1.075673352915104e+00;
                2.003796927111872e+00]

        w_p = [ 1.671879691147102e-02;
                1.636958371447360e-01;
                4.981856569770637e-01;
                8.372266245578912e+00;
                9.841730844088381e+00]
    elseif order == 10        
        a = 6
        chi_p = [1.175089381227308e-03;
             1.877034129831289e-02;
             9.686468391426860e-02;
             3.004818668002884e-01;
             6.901331557173356e-01;
             1.293695738083659e+00;
             2.090187729798780e+00;
             3.016719313149212e+00;
             4.001369747872486e+00;
             5.000025661793423e+00]
        w_p = [4.560746882084207e-03;
           3.810606322384757e-02;
           1.293864997289512e-01;
           2.884360381408835e-01;
           4.958111914344961e-01;
           7.077154600594529e-01;
           8.741924365285083e-01;
           9.661361986515218e-01;
           9.957887866078700e-01;
           9.998665787423845e-01]
    elseif order == 16
        a = 10;
        chi_p = [8.371529832014113e-04
                  1.239382725542637e-02 
                  6.009290785739468e-02
                  1.805991249601928e-01
                  4.142832599028031e-01
                  7.964747731112430e-01
                  1.348993882467059e+00
                  2.073471660264395e+00 
                  2.947904939031494e+00
                  3.928129252248612e+00 
                  4.957203086563112e+00
                  5.986360113977494e+00 
                  6.997957704791519e+00
                  7.999888757524622e+00  
                  8.999998754306120e+00]
        w_p = [3.190919086626234e-03
               2.423621380426338e-02
               7.740135521653088e-02
               1.704889420286369e-01
               3.029123478511309e-01
               4.652220834914617e-01
               6.401489637096768e-01
               8.051212946181061e-01
               9.362411945698647e-01
               1.014359775369075e+00
               1.035167721053657e+00
               1.020308624984610e+00
               1.004798397441514e+00
               1.000395017352309e+00
               1.000007149422537e+00]
    end
    chi_p = [chi_p' -chi_p[end:-1:1]']';
    w_p   = [w_p' w_p[end:-1:1]']';
    
    return (a, chi_p, w_p)

end