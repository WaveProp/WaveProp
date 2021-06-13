function assemble(iop,qdict;compress=Matrix)
    out        = compress(iop)
    correction = singular_weights(iop,qdict)
    axpy!(1,correction,out) # out <-- out + correction
end

function assemble(iop;compress=Matrix,k=2)
    Y = target_surface(iop)
    dict = Dict{DataType,SingularQuadratureRule}()
    for (E,q) in etypes(Y)
        # use a finer quadrature to compute singular integrals. Oversemple by a
        # factor of k
        q     = Integration.refine(q,k)
        # combine the oversample quadrature with a singularity handler to obtain
        # an appropriate SingularQuadratureRule
        shand = Kress(;order=3)
        qsing = SingularQuadratureRule(q,shand)
        # associate to each element type E a singular quadrature rule
        dict[E] = qsing
    end
    return assemble(iop,dict;compress)
end

function singular_weights(iop::IntegralOperator,qdict)
    X,Y = iop.X, iop.Y
    T = eltype(iop)
    K = iop.kernel
    dict_near = near_interaction_list(qnodes(X),Y;atol=0.1)
    Is    = Int[]
    Js    = Int[]
    Vs    = T[]
    for (E,list_near) in dict_near # for each element type
        qstd      = Y.etype2qrule[E]
        q         = qdict[E] # singular quadrature rule
        _singular_weights!(Is,Js,Vs,iop,E,qstd,q,list_near)
    end
    return Sp = sparse(Is,Js,Vs,size(iop)...)
end

function _singular_weights!(Is,Js,Vs,iop,E,qstd,q,list_near)
    X,Y = iop.X, iop.Y
    K   = iop.kernel
    ui        = qnodes(qstd)
    elt2dof   = Y.elt2dof[E]
    for (n,el) in enumerate(elements(Y)[E]) # for each element of type E in the mesh Y
        j_glob  = elt2dof[:,n]
        for (i,jloc) in list_near[n] # for each near observation point
            x   = X.qnodes[i]   # observation point
            vs  = ui[jloc]    # regularization center in parameter space
            if kernel_type(K) == SingleLayer()
                k    = (v) -> K(x,el(v))*measure(el,v)
            elseif kernel_type(K) == DoubleLayer()
                k    = (v) -> K(x,el(v),normal(el,v)) * measure(el,v)
            elseif kernel_type(K) == AdjointDoubleLayer()
                νx   = X.qnormals[i] # point x must also have a normal
                k    = (v) -> K(x,el(v),νx) * measure(el,v)
            elseif kernel_type(K) == HyperSingular()
                νx  = X.qnormals[i] # point x must also have a normal
                k   = (v) -> K(x,el(v),νx,normal(el,v)) * measure(el,v)
            else
                error("unkown kernel type")
            end
            w  = singular_weights(k,qstd,q,vs)
            append!(Is,fill(i,length(j_glob)))
            append!(Js,j_glob)
            w  = w - iop[i,j_glob]
            append!(Vs,w)
        end
    end
end
