"""
    near_interaction_list(X,Y;dim,atol)

Given a target mesh `X` and a source mesh `Y`, return a dictionary with keys given by element types of the source mesh `Y`. To each key, which represents an element type, we associate a vector encoding whose `i`-th entry encodes information about the points in `X` which are close to the a given element of the key type. 
"""
function near_interaction_list(X,Y;dim,atol)
    dict = Dict{DataType,Vector{Vector{Tuple{Int,Int}}}}()    
    for E in etypes(Y)
        geometric_dimension(E) == dim || continue
        ilist = _etype_near_interaction_list(X,Y,E,atol)
        push!(dict,E=>ilist)
    end
    return dict
end    

function _etype_near_interaction_list(X,Y,E,atol)
    ilist = Vector{Vector{Tuple{Int,Int}}}()
    e2n   = el2qnodes(Y,E)
    npts,nel = size(e2n)
    for n in 1:nel
        ynodes = qnodes(Y)[e2n[:,n]]
        inear  = _near_interaction_list(X,ynodes,atol)
        push!(ilist,inear)
    end 
    ilist       
end    

function _near_interaction_list(X,ynodes,atol)
    ilist    = Tuple{Int,Int}[]
    for (i,x) in enumerate(qnodes(X))            
        d,j = findmin([norm(x-y) for y in ynodes])        
        if d ≤ atol
            push!(ilist,(i,j))    
        end
    end            
    return ilist
end   