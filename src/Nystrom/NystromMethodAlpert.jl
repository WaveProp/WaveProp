using SpecialFunctions; # need to evaluate the Hankel function (Green's function)

function matricesAlpert(k,prt,order)   
    
    (nskip, chi, weights) = Alpert_params(order)
    # nskip = a
    ninterp = order + 3;

    τ = prt.τ
    ν = prt.ν
    s = prt.s
    ds= prt.ds
        
    N = prt.N
    Δs = prt.Δs
       

    x1 = prt.x[:,1];     x2 = prt.x[:,2]
    dx1 = prt.dx[:,1];   dx2 = prt.dx[:,2]
    
    S  = zeros(N,N)+im*zeros(N,N)
    K  = zeros(N,N)+im*zeros(N,N)
    Kp = zeros(N,N)+im*zeros(N,N)
    T  = zeros(N,N)+im*zeros(N,N)
        
    for i=1:N
        
        for j=1:N

            r = sqrt((x1[i]-x1[j])^2+(x2[i]-x2[j])^2);            

            H0 = hankelh1(0,k*r); H1 = hankelh1(1,k*r)

            # if i!=j
            
                S2  = im/4*H0*τ[j]*Δs*ds[j]
                K2  = im*k/4*H1/r*((x1[i]-x1[j])*dx2[j]-(x2[i]-x2[j])*dx1[j])*Δs*ds[j]
                Kp2 = im*k/4*H1/r*((x1[j]-x1[i])*dx2[i]-(x2[j]-x2[i])*dx1[i])*(τ[j]/τ[i])*Δs*ds[j]
                
                T2 = ((x1[i]-x1[j])*dx1[j]+(x2[i]-x2[j])*dx2[j])*((x1[i]-x1[j])*dx1[i]+(x2[i]-x2[j])*dx2[i])/r^2                
                T2 = im/2*T2*(k^2*H0-2*k*H1/r)+im*k/2*(dx1[i]*dx1[j]+dx2[i]*dx2[j])*H1/r
                T2 = -0.5*T2*Δs*ds[j]/τ[j]

            # end 

            S[i,j]  = S2
            K[i,j]  = K2
            Kp[i,j] = Kp2            
            T[i,j]  = T2
            
        end
        
        # @. S[i, mod(i+(-nskip+1:nskip-1)-1, N)+1] = 0

        S[i, mod.(i.+(-nskip+1:nskip-1).-1, N).+1] .= 0
        K[i, mod.(i.+(-nskip+1:nskip-1).-1, N).+1] .= 0  
        Kp[i,mod.(i.+(-nskip+1:nskip-1).-1, N).+1] .= 0
        T[i, mod.(i.+(-nskip+1:nskip-1).-1, N).+1] .= 0  
            
    end
    

    for p = 1:size(chi, 1)
        
        if prt.PD == 0
            t = Δs*collect(0:N-1)             
            w = mod.(t .+ Δs*chi[p], 2*pi)
            dw = ones(size(w))            
        else
            t = Δs*collect(1:N) 
            (w,dw) = polySigmoid(mod.(t .+ Δs*chi[p], 2*pi),prt.PD)
        end
                        
        (aux_x, aux_dx) = prt.param(w)
        # @show aux_x 
        # return
        
        aux_x1  = aux_x[:,1]
        aux_x2  = aux_x[:,2]

        aux_dx1 = aux_dx[:,1]
        aux_dx2 = aux_dx[:,2]

        tau_aux = sqrt.(aux_dx1.^2+aux_dx2.^2)
        
        aux_n1 = aux_dx2./tau_aux 
        aux_n2 =-aux_dx1./tau_aux
                
        aux_r = sqrt.((x1-aux_x1).^2+(x2-aux_x2).^2)

        aux_H0 = hankelh1.(0,k*aux_r)
        aux_H1 = hankelh1.(1,k*aux_r)

        ind = floor.(chi[p] - ninterp/2 + 1) .+ (0 : ninterp - 1)
        
        L   = lagrangeInterp(ind', chi[p])
        
        kerS = im/4*aux_H0.*tau_aux.*dw
        kerK = im*k/4*aux_H1./aux_r.*((x1-aux_x1).*aux_dx2-(x2-aux_x2).*aux_dx1).*dw
        kerKp= im*k/4*aux_H1./aux_r.*((aux_x1-x1).*ν[:,1]+(aux_x2-x2).*ν[:,2]).*(tau_aux.*dw)

        kerT = ((x1-aux_x1).*aux_dx1+(x2-aux_x2).*aux_dx2).*((x1-aux_x1).*dx1+(x2-aux_x2).*dx2)./aux_r.^2                
        kerT = im/2*kerT.*(k^2*aux_H0-2*k*aux_H1./aux_r)+
                im*k/2*(dx1.*aux_dx1+dx2.*aux_dx2)./aux_r.*aux_H1 
        kerT = -0.5*kerT.*dw./tau_aux
        
        for i0=1:N        
           j0  = i0 .+ ind
           j0  = floor.(Int,mod.(j0.-1, N) .+ 1)
           S[i0,j0] .= S[i0,j0]  + Δs*weights[p]*kerS[i0]*L[:,1]
           K[i0,j0] .= K[i0,j0]  + Δs*weights[p]*kerK[i0]*L[:,1]
           Kp[i0,j0].= Kp[i0,j0] + Δs*weights[p]*kerKp[i0]*L[:,1]
           T[i0,j0] .= T[i0,j0]  + Δs*weights[p]*kerT[i0]*L[:,1]
        end

    end
        
    T = k^2*S.*(ν*ν') + T
    
    mat = (S=S,K=K,Kp=Kp,T=T);

end

########################################################################

function  Alpert_params(order)
    
    if order == 2    
        a = 1
        chi_p = [.591549430918953e-01]
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
    w_p = [w_p' w_p[end:-1:1]']';
    
    return (a, chi_p, w_p)

end

##############################################################################

function lagrangeInterp(xx, t)
    #  L is a row vector calculating the coefficients of lagrange interpolation
    #  at t xx are data 
    n = length(xx) # here xx = index = [-2, -1, 0, 1, 2]
    denom = xx'*ones(1, n) - ones(n, 1)*xx
    numer = t*ones(n, n) - ones(n, 1)*xx
    for i = 1:n
      denom[i, i] = 1.0
      numer[i, i] = 1.0
    end
    temp = numer./denom
    L = prod(temp,dims=2)    
    return L
end