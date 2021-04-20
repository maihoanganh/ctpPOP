function get_indcons(m::Int64,supp_g::Vector{Vector{Vector{UInt64}}},I::Vector{Vector{UInt16}},p::Int64,lI::Vector{Int64};assign::String="first")
    J=[UInt64[] for i in 1:p]
    ind=UInt64(0)
    indvec=zeros(Int64,1)
    l_indvec=Int64(1)
    rind=zeros(Int64,1)
    temp=UInt64[]
    ncc=UInt64[]
    for i in 1:m
        rind=unique(vcat(supp_g[i]...))
        if assign=="first"
            ind=findfirst(k->issubset(rind,I[k]),1:p)
            if ind!=nothing
                push!(J[ind], i)
            else
                push!(ncc, i)
            end
        elseif assign=="all"
            indvec=findall(k->issubset(rind, I[k]), 1:p)
            l_indvec=length(indvec)
            if l_indvec!=0
                for j in 1:l_indvec
                    push!(J[indvec[j]], i)
                end
            else
                push!(ncc, i)
            end
        elseif assign=="min"
            for j in 1:p
                if issubset(rind, I[j])
                    push!(temp,j)
                end
            end
            if temp!=[]              
                push!(J[temp[argmin(lI[temp])]], i)
            else
                push!(ncc, i)
            end
            temp=UInt64[]   
        else
            println("No assign!!!")
            
        end
    end
    lJ=[length(J[i]) for i in 1:p]
    return J,lJ,ncc
end


function decomp_obj(supp_f,lmon_f,I,p,n)
    Indf=[UInt64[] for i in 1:p]
    
    ind=zeros(UInt64,1)
    for t in 1:lmon_f
        ind=unique(supp_f[t])#findall(j-> supp_f[j,t]!=0,1:n)
        append!(Indf[findfirst(i->issubset(ind,I[i]),1:p)],t)
    end   
    lIndf=[length(Indf[i]) for i in 1:p]
    return Indf,lIndf
end


function get_ncgraph(tsupp,basis;obj="eigen")
    lb=length(basis)
    G=SimpleGraph(lb)
    ltsupp=length(tsupp)
    for i = 1:lb, j = i+1:lb
        bi = [basis[i][end:-1:1]; basis[j]]
        bi=_sym_canon(bi)
        if obj=="trace"
            bi=_cyclic_canon(bi)
        end
        # if nx>0
        #     bi=comm(bi, nx)
        #     proj!(bi)
        # end
        if ncbfind(tsupp, ltsupp, bi)!=0
           add_edge!(G, i, j)
        end
    end
    return G
end


function clique_decomp(n::Int,m::Int,dg::Vector{Int},supp;order="min",alg="MD",minimize=false)
    if alg==false
        cliques=[UInt16[i for i=1:n]]
        cql=1
        cliquesize=[n]
    else
        G=SimpleGraph(n)
        for i=1:m+1
            if order=="min"||i==1||order==ceil(Int, dg[i-1]/2)
                for j = 1:supp[i].n
                    add_clique!(G,supp[i].rowval[supp[i].colptr[j]:(supp[i].colptr[j+1]-1)])
                end
            else
                add_clique!(G,unique(supp[i].rowval))
            end
        end
        if alg=="NC"
            cliques,cql,cliquesize=max_cliques(G)
        else
            cliques,cql,cliquesize=chordal_cliques!(G, method=alg, minimize=minimize)
        end
    end
    uc=unique(cliquesize)
    sizes=[sum(cliquesize.== i) for i in uc]
    println("------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("------------------------------------------------------")
    return cliques,cql,cliquesize
end