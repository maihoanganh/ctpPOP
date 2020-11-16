function get_blocks_TS(t::Int64,n::Int64,m::Int64,l::Int64,Usupp::Matrix{UInt64},sk::Int64,sk_g::Vector{UInt64},sk_h::Vector{UInt64},lmon_g::Vector{UInt64},lmon_h::Vector{UInt64},supp_g::Vector{Matrix{UInt64}},supp_h::Vector{Matrix{UInt64}},coe_g::Vector{Vector{Float64}},coe_h::Vector{Vector{Float64}},v::Matrix{UInt64})
    
    
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






function model_POP_TS(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Matrix{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Matrix{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::Matrix{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64,t::Int64;use_eqcons_to_get_constant_trace::Bool=true)

    v,sk,sk_g,s2k_g,sk_h,s2k_h,ak,P,Pg=get_constant_trace_dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)

    
    Usupp=[2*v[:,1:sk] supp_f]
    
    if m>0
        Usupp=[Usupp hcat(supp_g...)]
    end
        
    if l>0
        Usupp=[Usupp hcat(supp_h...)]
    end
   
    Usupp,lUsupp,block_g0,block_g,block_h,lblock_g0,lblock_g,lblock_h,lt_block_g0,lt_block_g,lt_block_h=get_blocks_TS(t,n,m,l,Usupp,sk,sk_g,sk_h,lmon_g,lmon_h,supp_g,supp_h,coe_g,coe_h,v[:,1:sk])
    
    omega=Int64(lblock_g0+sum(lblock_g))
    
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
    
    
    a_val,a0,norm_a0,opnorm_a=rescale_dense(a_ind1,a_ind2,a_val,a_len,a0,zeta)
    
    
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
    
    a_block=Vector{SparseMatrixCSC{Float64}}(undef,omega)
       
    a0_block=Vector{Adjoint{Float64,Array{Float64,1}}}(undef,omega)
    
    @fastmath @inbounds @simd  for i in 1:omega
        a_block[i]=sparse(a_ind2[Ind[i]],a_ind1[Ind[i]].-r_blo[i],a_val[Ind[i]],zeta,u_vec[i])
        a0_block[i]=a0[r_blo[i]+1:r_blo[i]+u_vec[i]]'
    end
    
  
   
    
    println("Modeling time:")
 
    return omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak
end



function POP_TS_CGAL(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Matrix{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Matrix{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::Matrix{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64,t::Int64;maxit::Int64=Int64(1e5),tol::Float64=1e-4,use_eqcons_to_get_constant_trace::Bool=true,check_tol_each_iter::Bool=true)

    @time begin
    
    
    @time omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak=model_POP_TS(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,t,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)
    

      
    @time opt_val=solve_POP_dense_CGAL(omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak,maxit=maxit,tol=tol)
    
              
       println("Total time:")
                    end
    return opt_val
end