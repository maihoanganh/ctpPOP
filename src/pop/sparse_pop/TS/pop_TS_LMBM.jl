function get_blocks_TS_LMBM(t::Int64,n::Int64,m::Int64,l::Int64,Usupp::Matrix{UInt64},sk::Int64,sk_g::Vector{Int64},sk_h::Vector{Int64},lmon_g::Vector{UInt64},lmon_h::Vector{UInt64},supp_g::Vector{Matrix{UInt64}},supp_h::Vector{Matrix{UInt64}},coe_g::Vector{Vector{Float64}},coe_h::Vector{Vector{Float64}},v::Matrix{UInt64})
    
    
    Usupp=sortslices(Usupp,dims=2)
    Usupp=unique(Usupp,dims=2)
    lUsupp=size(Usupp,2)

    block_g0=Vector{Vector{UInt64}}(undef,1)
    block_g=Vector{Vector{Vector{UInt64}}}(undef,m)
    block_h=Vector{Vector{Vector{UInt64}}}(undef,l)

    lblock_g0=UInt64(0)
    lblock_g=Vector{UInt64}(undef,m)
    lblock_h=Vector{UInt64}(undef,l)

    lt_block_g0=Vector{UInt64}(undef,1)
    lt_block_g=Vector{Vector{UInt64}}(undef,m)
    lt_block_h=Vector{Vector{UInt64}}(undef,l)

    old_block_g0=Vector{Vector{UInt64}}(undef,1)
    old_block_g=Vector{Vector{Vector{UInt64}}}(undef,m)
    old_block_h=Vector{Vector{Vector{UInt64}}}(undef,l)
    
    y=1
    graph=LightGraphs.SimpleGraphs.SimpleGraph{Int64}

    iter=1
    t_supp=1
    while iter<=t
        graph=SimpleGraph(sk)
        for p in 1:sk, q in p:sk
            if bfind(Usupp,lUsupp,v[:,p]+v[:,q],n)!=0
               add_edge!(graph,p,q)
            end
        end
        block_g0=connected_components(graph)
        lblock_g0=length(block_g0)
        lt_block_g0=[length(block_g0[j]) for j in 1:lblock_g0]

        
        for i in 1:m
            graph=SimpleGraph(sk_g[i])
            for p in 1:sk_g[i]
                for q in p:sk_g[i] 
                    while y<=lmon_g[i]
                        if bfind(Usupp,lUsupp,v[:,p]+v[:,q]+supp_g[i][:,y],n)!=0
                            break
                        else
                            y+=1
                        end
                    end
                    if y<=lmon_g[i]
                       add_edge!(graph,p,q)
                    end
                    y=1
                end
            end
            block_g[i]=connected_components(graph)
            lblock_g[i]=length(block_g[i])
            lt_block_g[i]=[length(block_g[i][j]) for j in 1:lblock_g[i]]
        end


       
        for i in 1:l
            graph=SimpleGraph(sk_h[i])
            for p in 1:sk_h[i]
                for q in p:sk_h[i]
                    while y<=lmon_h[i]
                        if bfind(Usupp,lUsupp,v[:,p]+v[:,q]+supp_h[i][:,y],n)!=0
                            break
                        else
                            y+=1
                        end
                    end
                    if y<=lmon_h[i]
                       add_edge!(graph,p,q)
                    end
                    y=1
                end
            end
            block_h[i]=connected_components(graph)
            lblock_h[i]=length(block_h[i])
            lt_block_h[i]=[length(block_h[i][j]) for j in 1:lblock_h[i]]
        end

        Usupp=zeros(UInt64,n,Int64(0.5*sum(lt_block_g0[j]*(lt_block_g0[j]+1) for j in 1:lblock_g0)))
        t_supp=1
        for j in 1:lblock_g0
            for p in 1:lt_block_g0[j], q in p:lt_block_g0[j]
                Usupp[:,t_supp]=v[:,block_g0[j][p]]+v[:,block_g0[j][q]]
                t_supp+=1
            end
        end
        
        #=
       for i in 1:m
            for j in 1:lblock_g[i]
                for p in 1:lt_block_g[i][j], q in p:lt_block_g[i][j]
                    for z in 1:lmon_g[i]
                        Usupp=[Usupp v[:,block_g[i][j][p]]+v[:,block_g[i][j][q]]+supp_g[i][:,z]]
                    end
                end
            end
        end
        for i in 1:l
            for j in 1:lblock_h[i]
                for p in 1:lt_block_h[i][j], q in p:lt_block_h[i][j]
                    for z in 1:lmon_h[i]
                        Usupp=[Usupp v[:,block_h[i][j][p]]+v[:,block_h[i][j][q]]+supp_h[i][:,z]]
                    end
                end
            end
        end
        =#
        Usupp=sortslices(Usupp,dims=2)
        Usupp=unique(Usupp,dims=2)
        lUsupp=size(Usupp,2)


        if iter==1
            old_block_g0=block_g0
            old_block_g=block_g
            old_block_h=block_h
        else
            if old_block_g0==block_g0 && old_block_g==block_g && old_block_h==block_h
                println("  Stable sparse order: t=",iter-1)
                break
            else
                old_block_g0=block_g0
                old_block_g=block_g
                old_block_h=block_h
            end
        end
        iter+=1
    end
    #=
    println("  -----------------------------")
    println("  #block_g0 = ",[lt_block_g0[i] for i in 1:lblock_g0])
    println("  -----------------------------")
    println("  #block_g = ",[[lt_block_g[i][j] for j in 1:lblock_g[i]] for i in 1:m])
    println("  -----------------------------")
    println("  #block_h = ",[[lt_block_h[i][j] for j in 1:lblock_h[i]] for i in 1:l])
    println("  -----------------------------")
    =#
    return Usupp,lUsupp,block_g0,block_g,block_h,lblock_g0,lblock_g,lblock_h,lt_block_g0,lt_block_g,lt_block_h

end


function get_constant_trace_TS_LMBM(x,g,h,k;UseEq=true)
    n=length(x)
    m=length(g)
    l=length(h)
            
    v=get_basis(n,2*k)
    
    sk=binomial(k+n,n)
    sk_g=Vector{Int64}(undef,m)
    sk_h=Vector{Int64}(undef,l)
    s2k_g=Vector{Int64}(undef,m)
    s2k_h=Vector{Int64}(undef,l)
    
    lmon_g=Vector{UInt64}(undef,m)
    supp_g=Vector{Matrix{UInt64}}(undef,m)
    coe_g=Vector{Vector{Float64}}(undef,m)
    
    lmon_h=Vector{UInt64}(undef,l)
    supp_h=Vector{Matrix{UInt64}}(undef,l)
    coe_h=Vector{Vector{Float64}}(undef,l)
        
    
    ceil0=Int64(0)
    for i in 1:m
        ceil0=ceil(Int64,maxdegree(g[i])/2)
        sk_g[i]=binomial(k-ceil0+n,n)
        s2k_g[i]=binomial(2*(k-ceil0)+n,n)
        lmon_g[i],supp_g[i],coe_g[i]=info(g[i],x,n)
    end
                    
    for i in 1:l
        ceil0=ceil(Int64,maxdegree(h[i])/2)
        sk_h[i]=binomial(k-ceil0+n,n)
        s2k_h[i]=binomial(2*(k-ceil0)+n,n)
        lmon_h[i],supp_h[i],coe_h[i]=info(h[i],x,n)
    end
    
    lV=sk
    
    if m!=0
        lV+=sum(sk_g[i]*lmon_g[i] for i in 1:m)
    end
    
    if UseEq
        if l!=0
            lV+=sum(lmon_h[i]*s2k_h[i] for i in 1:l)
        end
    end
     
    V=Matrix{UInt64}(undef,n,lV)
    V[:,1:sk]=2*v[:,1:sk]
    t=sk+1
    for i in 1:m
        for r in 1:lmon_g[i]
            for j in 1:sk_g[i]
                @inbounds V[:,t]=2*v[:,j]+supp_g[i][:,r]
                @inbounds t+=1
            end
        end
    end
     
    if UseEq
        for i in 1:l
            for r in 1:lmon_h[i]
                for j in 1:s2k_h[i]
                    @inbounds V[:,t]=v[:,j]+supp_h[i][:,r]
                    @inbounds t+=1
                end
            end
        end    
    end  
    V=unique(V,dims=2)
    V=sortslices(V,dims=2)
    
    
    lV=size(V,2)
                                    
    model=JuMP.Model(with_optimizer(Mosek.Optimizer, QUIET=true))
    cons=[AffExpr(0) for i in 1:lV]

    p0=@variable(model, [1:sk], lower_bound=1)
    p=Vector{Vector{VariableRef}}(undef, m)
    

    
    for j in 1:sk
        @inbounds add_to_expression!(cons[bfind(V,lV,2*v[:,j],n)],p0[j])
    end


    for i in 1:m
        p[i]=@variable(model, [1:sk_g[i]], lower_bound=1)
        for j in 1:sk_g[i]
            for r in 1:lmon_g[i]
                @inbounds add_to_expression!(cons[bfind(V,lV,2*v[:,j]+supp_g[i][:,r],n)],p[i][j]*coe_g[i][r])
            end
        end
    end

    
    if UseEq
        q=Vector{Vector{VariableRef}}(undef, l)   
        for i in 1:l
            q[i]=@variable(model, [1:s2k_h[i]])
            for j in 1:s2k_h[i]
                for r in 1:lmon_h[i]
                    @inbounds add_to_expression!(cons[bfind(V,lV,v[:,j]+supp_h[i][:,r],n)],q[i][j]*coe_h[i][r])
                end
            end
        end
    end
    @constraint(model, cons[2:end].==0)
    @variable(model, lambda)
    @constraint(model, cons[1]==lambda)
    @objective(model, Min, lambda)
    optimize!(model)
                                    
    println("  Computing constant trace status: ", termination_status(model))
                                    
    ak=value(lambda)
    println("  Constant trace: ak = ",ak)
    
    P=sqrt.(value.(p0))
    Pg=[sqrt.(value.(p[i])) for i in 1:m]
    
    return n,m,l,v,sk,sk_g,sk_h,s2k_g,s2k_h,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,ak,P,Pg
end



function model_POP_TS_LMBM(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64,t::Int64;UseEq::Bool=true)

    n,m,l,v,sk,sk_g,sk_h,s2k_g,s2k_h,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,ak,P,Pg=get_constant_trace_TS_LMBM(x,g,h,k,UseEq=UseEq)

    

    lmon_f,supp_f,coe_f=info(f,x,n)

    
    Usupp=[2*v[:,1:sk] supp_f]
    
    if m>0
        Usupp=[Usupp hcat(supp_g...)]
    end
        
    if l>0
        Usupp=[Usupp hcat(supp_h...)]
    end
    Usupp,lUsupp,block_g0,block_g,block_h,lblock_g0,lblock_g,lblock_h,lt_block_g0,lt_block_g,lt_block_h=get_blocks_TS_CGAL(t,n,m,l,Usupp,sk,sk_g,sk_h,lmon_g,lmon_h,supp_g,supp_h,coe_g,coe_h,v[:,1:sk])
    
    omega=lblock_g0+sum(lblock_g)
    
    println("  Number of blocks: omega=",omega)
    println("  Size of the largest block: s^max=",maximum([maximum(lt_block_g0);maximum(maximum(lt_block_g[i]) for i in 1:m)]))
    
   
    
    Order(alpha)=bfind(Usupp,lUsupp,alpha,n)
    
    w=[Int64(0.5*lt_block_g0[j]*(lt_block_g0[j]+1)) for j in 1:lblock_g0]
    w_g=[[@inbounds Int64(0.5*lt_block_g[i][j]*(lt_block_g[i][j]+1)) for j in 1:lblock_g[i]] for i in 1:m]
    
    u=sum(w)
    u_g=[sum(w_g[i]) for i in 1:m]

    
    
    
    d=u+sum(u_g)
   
    
    
    zeta=d-lUsupp
    
    

    IndM=[@inbounds Vector{Vector{Int64}}([]) for j in 1:lUsupp]
    invIndeM=spzeros(UInt64,sk,sk)
    r=UInt64(0)
    t_iter=UInt64(1)
    
    a_ind1=Vector{UInt64}([])
    a_ind2=Vector{UInt64}([])
    a_val=Vector{Float64}([])
    a_len=UInt64(0)
    
    t_blo=0
    
    for j in 1:lblock_g0
        for p in 1:lt_block_g0[j], q in 1:p
            @inbounds r=Order(v[:,block_g0[j][p]]+v[:,block_g0[j][q]])
            @inbounds push!(IndM[r],[block_g0[j][p],block_g0[j][q]])
            @inbounds invIndeM[block_g0[j][p],block_g0[j][q]]=t_iter
            @inbounds t_iter+=1
        end
        t_blo+=lt_block_g0[j]
    end
    l_IndM=[length(IndM[r]) for r in 1:lUsupp]
   
    
    
    t_a=1
    
    
    I=zeros(UInt64,2)
            
    for r in 1:lUsupp
        if l_IndM[r]>1
            for i in 2:l_IndM[r]
                I=IndM[r][1]
                push!(a_ind1,invIndeM[I[1],I[2]])
                push!(a_ind2,t_a)
                push!(a_val,ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
                a_len+=1
                
                I=IndM[r][i]
                push!(a_ind1,invIndeM[I[1],I[2]])
                push!(a_ind2,t_a)
                push!(a_val,-ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
               
                a_len+=1
                t_a+=1
             end
             IndM[r]=Vector{Int64}[IndM[r][end]]
         end
    end
    
    
    Usupp_g=Vector{Matrix{UInt64}}(undef,m)
    lUsupp_g=zeros(UInt64,m)

    
    for i in 1:m     
        Usupp_g[i]=zeros(UInt64,n,u_g[i])
        t_iter=1
        for j in 1:lblock_g[i]
            for p in 1:lt_block_g[i][j], q in p:lt_block_g[i][j]
                @inbounds Usupp_g[i][:,t_iter]=v[:,block_g[i][j][p]]+v[:,block_g[i][j][q]]
                t_iter+=1
            end
        end
        
        Usupp_g[i]=sortslices(Usupp_g[i],dims=2)
        Usupp_g[i]=unique(Usupp_g[i],dims=2)
        lUsupp_g[i]=size(Usupp_g[i],2) 
    end
   
    IndMg=[[Vector{Vector{Int64}}([]) for j in 1:lUsupp_g[i]] for i in 1:m]
    invIndeMg=[spzeros(UInt64,sk_g[i],sk_g[i]) for i in 1:m]
    t_Blo=u
 
 
    for i in 1:m
        r=UInt64(0)
        t_iter=UInt64(1)

        t_blo=0
        for j in 1:lblock_g[i]
            for p in 1:lt_block_g[i][j], q in 1:p
                @inbounds r=bfind(Usupp_g[i],lUsupp_g[i],v[:,block_g[i][j][p]]+v[:,block_g[i][j][q]],n)
                @inbounds push!(IndMg[i][r],[block_g[i][j][p];block_g[i][j][q]])
                @inbounds invIndeMg[i][p+t_blo,q+t_blo]=t_iter
                t_iter+=1
            end
            t_blo+=lt_block_g[i][j]
        end

                

        l_IndM=[length(IndMg[i][q]) for q in 1:lUsupp_g[i]]
        
        for rr in 1:lUsupp_g[i]
            if l_IndM[rr]>1
                for j in 2:l_IndM[rr]
                    I=IndMg[i][rr][1]
                    push!(a_ind1,invIndeMg[i][I[1],I[2]]+t_Blo)
                    push!(a_ind2,t_a)
                    push!(a_val,ak/Pg[i][I[1]]/Pg[i][I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
                    a_len+=1
                    
                    
                    
                    I=IndMg[i][rr][j]
                    push!(a_ind1,invIndeMg[i][I[1],I[2]]+t_Blo)
                    push!(a_ind2,t_a)
                    push!(a_val,-ak/Pg[i][I[1]]/Pg[i][I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
                    a_len+=1
                    
                    t_a+=1
                end
                IndMg[i][rr]=Vector{Int64}[IndMg[i][rr][end]]
            end
            
            
           
            I=IndMg[i][rr][1]
            push!(a_ind1,invIndeMg[i][I[1],I[2]]+t_Blo)
            push!(a_ind2,t_a)
            push!(a_val,-ak/Pg[i][I[1]]/Pg[i][I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
            a_len+=1
            
            for p in 1:lmon_g[i]  
                I=IndM[Order(supp_g[i][:,p]+Usupp_g[i][:,rr])][end]
                push!(a_ind1,invIndeM[I[1],I[2]])
                push!(a_ind2,t_a)
                push!(a_val,coe_g[i][p]*ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
                a_len+=1
            end
            t_a+=1
            
         end
         t_Blo+=u_g[i]
    end
    
    Usupp_h=Vector{Matrix{UInt64}}(undef,l)
    lUsupp_h=zeros(UInt64,l)
    
    
    
    for i in 1:l
        Usupp_h[i]=zeros(UInt64,n,Int(0.5*sum(lt_block_h[i][j]*(lt_block_h[i][j]+1) for j in 1:lblock_h[i])))
        t_iter=1
        for j in 1:lblock_h[i]
            for p in 1:lt_block_h[i][j], q in p:lt_block_h[i][j]
                @inbounds Usupp_h[i][:,t_iter]=v[:,block_h[i][j][p]]+v[:,block_h[i][j][q]]
                t_iter+=1
            end
        end
        Usupp_h[i]=sortslices(Usupp_h[i],dims=2)
        Usupp_h[i]=unique(Usupp_h[i],dims=2)
        lUsupp_h[i]=size(Usupp_h[i],2)
    end
    

    zeta+=sum(lUsupp_h)+1
    
    println("  Number of equality trace constraints: zeta=",zeta)
    
    
    @simd for j in 1:l
        @simd for r in 1:lUsupp_h[j]
            @simd for p in 1:lmon_h[j]
                      I=IndM[Order(supp_h[j][:,p]+Usupp_h[j][:,r])][1]
                      push!(a_ind1,invIndeM[I[1],I[2]])
                      push!(a_ind2,t_a)
                      push!(a_val,coe_h[j][p]*ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
                      a_len+=1
                
                   end
                    @inbounds t_a+=1
               end       
    end  
    
    
    
    push!(a_ind1,1)
    push!(a_ind2,t_a)
    push!(a_val,ak)
    a_len+=1
   
            


    a0=zeros(Float64,d)
    
    
    @simd for p in 1:lmon_f
        I=IndM[Order(supp_f[:,p])][1]
        a0[invIndeM[I[1],I[2]]]=coe_f[p]*ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1)
   end
    
    
    a_val,a0,norm_a0,opnorm_a=rescale_TS_LMBM(a_ind1,a_ind2,a_val,a_len,a0,zeta)
    
    
    s=Vector{UInt64}(undef,omega)
    r_blo=Vector{UInt64}(undef,omega)
    u_vec=Vector{UInt64}(undef,omega)
    
    t_s=1
    t_Blo=0
    
    for j in 1:lblock_g0
        s[t_s]=lt_block_g0[j]
        r_blo[t_s]=t_Blo
        u_vec[t_s]=w[j]
        
        t_Blo+=w[j]
        t_s+=1
    end
    
    for i in 1:m
        for j in 1:lblock_g[i]
            s[t_s]=lt_block_g[i][j]
            r_blo[t_s]=t_Blo
            u_vec[t_s]=w_g[i][j]
            
            t_Blo+=w_g[i][j]
            t_s+=1
        end 
    end
    
    
    #test_TS_Mosek(a0,a,s,zeta,d,omega)
    
    
    Ind=Vector{Vector{UInt64}}(undef,omega)
    Index=Vector{UInt64}(1:a_len)
    pos=zeros(UInt64,1)
    
    @fastmath @inbounds @simd  for j in 1:omega
        pos=findall(t->(a_ind1[t]>=(r_blo[j]+1)) && (a_ind1[t]<=(r_blo[j]+u_vec[j])),Index)
        Ind[j]=Index[pos]
        deleteat!(Index, pos)
    end
    
    a_block=[@fastmath @inbounds sparse(a_ind2[Ind[i]],a_ind1[Ind[i]].-r_blo[i],a_val[Ind[i]],zeta,u_vec[i]) for i in 1:omega]
    
    #@time a_block=[SparseMatrixCSC{Float64}(a[r_blo[i]+1:r_blo[i]+u_vec[i],:]') for i in 1:omega]
    a0_block=[@fastmath @inbounds a0[r_blo[i]+1:r_blo[i]+u_vec[i]]' for i in 1:omega]
 
    
  
   
    
    println("Modeling time:")
 
    return omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak,P,v[:,1:sk],n,m,l
        
end

function rescale_TS_LMBM(a_ind1::Vector{UInt64},a_ind2::Vector{UInt64},a_val::Vector{Float64},a_len::UInt64,a0::Vector{Float64},zeta::UInt64)
    

    norm_a=zeros(Float64,zeta)
    
    @fastmath @inbounds @simd  for t in 1:a_len
        norm_a[a_ind2[t]]+=a_val[t]^2
    end
        
    norm_a=sqrt.(norm_a)
   
    
   
    @fastmath @inbounds @simd  for t in 1:a_len
        a_val[t]/=norm_a[a_ind2[t]]
    end
    
    
    opnorm_a=svds(sparse(a_ind1,a_ind2,a_val), nsv = 1)[1].S[1]
    
    a_val=a_val./opnorm_a

    
    norm_a0=norm(a0)
    a0=a0./norm_a0
    
    
    return a_val,a0,norm_a0,opnorm_a
end




function SmallEig_TS_LMBM(mat::Matrix{Float64},s::UInt64;EigAlg::String="Arpack")
    if s==1
        return mat[1,1], ones(Float64,1)
    #=elseif norm(mat)==0
        return 0.0, [1;zeros(Float64,s-1)]=#
    elseif EigAlg=="Arpack"
       @fastmath @inbounds E=eigs(mat,nev = 1,which=:SR,tol=1e-2) 
       return E[1][1],E[2][:,1]
    elseif EigAlg=="ArnoldiMethod"      
       @fastmath @inbounds E=partialeigen(partialschur(mat, nev=1,tol=1e-2, which=SR())[1])
       return E[1][1],E[2][:,1]
    elseif EigAlg=="Normal"
       @fastmath @inbounds E=eigen(Symmetric(mat),1:1)
       return E.values[1],E.vectors[:,1]
    elseif EigAlg=="Mix"
       try 
           @fastmath @inbounds E=partialeigen(partialschur(mat, nev=1,tol=1e-2, which=SR())[1])
           return E[1][1],E[2][:,1]
       catch
           @fastmath @inbounds E=eigs(mat,nev = 1,which=:SR,tol=1e-2) 
           return E[1][1],E[2][:,1]
       end
    else
       println("No eigenvalue algorithm!!!")
    end  
    
end

                                        
function SmallEig_block_TS_LMBM(vec::Adjoint{Float64,Array{Float64,1}},sk::UInt64;EigAlg::String="Arpack")
    
    B=zeros(Float64,sk,sk)
    r=1
    @fastmath @inbounds for i in 1:sk, j in 1:i
        B[i,j]=vec[r]/sqrt(1+(j<i))
        B[j,i]= copy(B[i,j])
        r+=1
    end
 
    return SmallEig_TS_LMBM(B,sk,EigAlg=EigAlg)

end




function SmallEigBlocks_TS_LMBM(a0_block::Vector{Adjoint{Float64,Array{Float64,1}}},a_block::Vector{SparseMatrixCSC{Float64,UInt64}},y::Vector{Float64},s::Vector{UInt64},omega::UInt64;EigAlg::String="Arpack")
    
    
    eigval,eigvec=Float64(Inf),zeros(Float64,1)

    smalleigval=eigval
    smalleigvec=eigvec
    ind=0

 
    @fastmath @inbounds @simd for p in 1:omega
        eigval,eigvec=SmallEig_block_TS_LMBM(a0_block[p]+y'*a_block[p],s[p],EigAlg=EigAlg)
        if smalleigval>eigval
            smalleigval=eigval
            smalleigvec=eigvec
            ind=p
        end
    end

    return smalleigval,smalleigvec,ind
end       



function solve_POP_TS_LMBM(omega::UInt64,a0_block::Vector{Adjoint{Float64,Array{Float64,1}}},a_block::Vector{SparseMatrixCSC{Float64,UInt64}},s::Vector{UInt64},zeta::UInt64,norm_a0::Float64,opnorm_a::Float64,ak::Float64;EigAlg::String="Arpack",tol::Float64=1e-4)
                
    function phi(nvar::Cint,xp::Ptr{Cdouble},gp::Ptr{Cdouble})
        zvar=unsafe_wrap(Array, xp, (convert(Int, nvar),))
        grad=unsafe_wrap(Array, gp, (convert(Int, nvar),))  
        eigval,eigvec,ind=SmallEigBlocks_TS_LMBM(a0_block,a_block,zvar,s,omega,EigAlg=EigAlg)
        grad[:]=-a_block[ind]*[@fastmath @inbounds eigvec[i]*eigvec[j]*sqrt(1+(j<i)) for i in 1:s[ind] for j in 1:i]                
        grad[nvar]+=1/opnorm_a/ak
        return(convert(Cdouble,-eigval+zvar[nvar]/opnorm_a/ak))
    end

                        
                        
    opt_val,z= lmbm(phi,zeros(Float64,zeta);printinfo=true,tol=tol)
                     
                        
    opt_val*=-norm_a0
    println()
    println("####################################")
    println("opt_val = ",opt_val)
    println("####################################")
    println("Solving time:")
    return opt_val,z            
                
end
    


function POP_TS_LMBM(x::Vector{PolyVar{true}},f::Polynomial{true},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64,t::Int64;EigAlg="Arpack",tol=1e-4,UseEq=true)

    @time begin
    
    
    @time omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak,P,vsk,n,m,l=model_POP_TS_LMBM(x,f,g,h,k,t,UseEq=UseEq)
    
 
    @time opt_val,z=solve_POP_TS_LMBM(omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak,EigAlg=EigAlg,tol=tol)
    
    @time opt_sol=extract_optimizer_TS_LMBM(a0_block[1]+z'*a_block[1],P,s[1],vsk,n,m,l,opt_val,f,g,h,x,EigAlg=EigAlg) 
       println("Total time:")
                    end
    return opt_val,opt_sol
end
                    
function get_Grammat_TS_LMBM(vec::Adjoint{Float64,Array{Float64,1}},sk::UInt64,P::Vector{Float64};EigAlg="Arpack")
    t_iter=1
    Gr=zeros(Float64,sk,sk)
    @fastmath @inbounds for i in 1:sk, j in i:sk
        Gr[i,j]=vec[t_iter]/sqrt(1+(j<i))
        Gr[j,i]= copy(Gr[i,j])
        t_iter+=1
    end
    eigval,~=SmallEig_TS_LMBM(-Gr,sk,EigAlg=EigAlg)
    @fastmath @inbounds @simd for i in 1:sk
        Gr[i,i]+=eigval
    end
    t_iter=1               
    @fastmath @inbounds for i in 1:sk, j in i:sk
        Gr[i,j]*=P[i]*P[j]
        Gr[j,i]= copy(Gr[i,j])
        t_iter+=1
    end                   
    return Gr
end
                
                
                
function extract_optimizer_TS_LMBM(vec::Adjoint{Float64,Array{Float64,1}},P::Vector{Float64},lu0::UInt64,basis_sigma0::Matrix{UInt64},n::Int64,l_g::Int64,l_h::Int64,opt_val::Float64,f::Polynomial{true},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},x::Vector{PolyVar{true}};EigAlg="Arpack")
     
    Gr=get_Grammat_LMBM(vec,lu0,P,EigAlg=EigAlg)                
                    
    #extraction of optimizers
    V=nullspace(Gr,atol=1e-5)
    rk=size(V,2)

    println("Dimension of the null space of Gram matrix = ", rk)
    sol=Vector{Float64}[]
    if rk >0
        # extraction of Henrion and Lasserre
        V= rref_with_pivots!(Matrix(V'),1e-2)
        U=Matrix(V[1]')

        w=basis_sigma0[:,V[2]]
        N=Vector{Array{Float64}}(undef,n)
        flag=UInt8(1)
 
        
        for i in 1:n
            @inbounds kk=UInt16[]
            for j=1:size(w)[2]
                @inbounds xwj=w[:,j]
                @inbounds xwj[i]+=1
                for t=1:lu0
                    if xwj==basis_sigma0[:,t]
                        if kk==UInt16[]
                            @inbounds kk=t
                        else
                            @inbounds kk=[kk;t]
                        end
                        break
                    end
                end
            end
            
            if rk!=length(kk)
                @inbounds flag=0
                break
            else
                @inbounds N[i]=U[kk,:]
            end
        end

        if flag==1

            # Create random convex combination
            rands = rand(n);rands = rands/sum(rands)
            M = zeros(length(V[2]),rk)
            for i=1:n
                @inbounds M+=rands[i]*N[i]
            end

            F= schur(M);
            L=F.Z
            # Extract solution
            for i=1:rk
                @inbounds atom=Float64[]
                for j = 1:n
                    @inbounds coordinatej=L[:,i]'*N[j]*L[:,i]
                    @inbounds push!(atom,coordinatej[1])
                end
                println("------------------------------------")
                println("atom ",i," = ",atom)
                @inbounds flag=1
                
                @inbounds check=polynomial(f)(x => atom)-opt_val
                
                println("  check gap of lower bound  = ",check)
                if abs(check)>1e-1
                    @inbounds flag=0
                end


                for i=1:l_g
                    @inbounds check=polynomial(g[i])(x => atom)
                    println("  check inequality constraint ",i," = ",check)
                    if check<-1e-1
                        @inbounds flag=0
                    end
                end
                
                for i=1:l_h
                    @inbounds check=polynomial(h[i])(x => atom)
                    println("  check equality constraint ",i," = ",check)
                    if abs(check)>1e-1
                        @inbounds flag=0
                    end
                end

                if flag ==1
                    @inbounds sol=atom
                    println("####################################")
                    println("Optimal solution: opt_sol = ",atom)
                    println("####################################")
                end

            end
        end
    end
    
    println("Extracting solutuion time:")

    return sol

end