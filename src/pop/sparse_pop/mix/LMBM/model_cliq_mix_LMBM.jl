function initial_mix_LMBM(n::Int64,m::Int64,l::Int64,k::Int64,dg::Vector{Int64},dh::Vector{Int64})

    v=get_basis(n,2*k)
    s2k=size(v,2)
    
    sk=binomial(k+n,n)
    sk_g=Vector{UInt64}(undef,m)
    sk_h=Vector{UInt64}(undef,l)
    s2k_g=Vector{UInt64}(undef,m)
    s2k_h=Vector{UInt64}(undef,l)
    
        
    
    ceil0=Int64(0)
    for i in 1:m
        ceil0=ceil(Int64,dg[i]/2)
        sk_g[i]=binomial(k-ceil0+n,n)
        s2k_g[i]=binomial(2*(k-ceil0)+n,n)
    end
                    
    for i in 1:l
        ceil0=ceil(Int64,dh[i]/2)
        sk_h[i]=binomial(k-ceil0+n,n)
        s2k_h[i]=binomial(2*(k-ceil(Int64,dh[i]/2))+n,n)
    end
    return v,s2k,sk,sk_g,sk_h,s2k_g,s2k_h
end


function get_blocks_Cliq_mix_LMBM(k::Int64,n::Int64,m::Int64,l::Int64,Usupp::SparseMatrixCSC{UInt64,Int64},lmon_g::Vector{UInt64},lmon_h::Vector{UInt64},supp_g,supp_h,coe_g::Vector{Vector{Float64}},coe_h::Vector{Vector{Float64}},v::Matrix{UInt64},s2k::UInt64,sk::UInt64,sk_g::Vector{UInt64},sk_h::Vector{UInt64},s2k_g::Vector{UInt64},s2k_h::Vector{UInt64})
    
    #println(typeof(Usupp))
    
    Usupp=sortslices(Usupp,dims=2)
    Usupp=unique(Usupp,dims=2)
    lUsupp=size(Usupp,2)


    block_g=Vector{Vector{Vector{UInt64}}}(undef,m)
    block_h=Vector{Vector{Vector{UInt64}}}(undef,l)


    lblock_g=Vector{UInt64}(undef,m)
    lblock_h=Vector{UInt64}(undef,l)


    lt_block_g=Vector{Vector{UInt64}}(undef,m)
    lt_block_h=Vector{Vector{UInt64}}(undef,l)

    
    


    
    
    graph=SimpleGraph(sk)
    for p in 1:sk, q in 1:p
        if bfind(Usupp,lUsupp,v[:,p]+v[:,q],n)!=0
           add_edge!(graph,p,q)
        end
    end
    block_g0=connected_components(graph)
    lblock_g0=length(block_g0)
    lt_block_g0=[length(block_g0[j]) for j in 1:lblock_g0]

    y=1
    for i in 1:m
        graph=SimpleGraph(sk_g[i])
        for p in 1:sk_g[i]
            for q in 1:p
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
            for q in 1:p
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
        for p in 1:lt_block_g0[j], q in 1:p
            Usupp[:,t_supp]=v[:,block_g0[j][p]]+v[:,block_g0[j][q]]
            t_supp+=1
        end
    end
    
    #=
   for i in 1:m
        for j in 1:lblock_g[i]
            for p in 1:lt_block_g[i][j], q in 1:p
                for z in 1:lmon_g[i]
                    Usupp=[Usupp v[:,block_g[i][j][p]]+v[:,block_g[i][j][q]]+supp_g[i][:,z]]
                end
            end
        end
    end
    for i in 1:l
        for j in 1:lblock_h[i]
            for p in 1:lt_block_h[i][j], q in 1:p
                for z in 1:lmon_h[i]
                    Usupp=[Usupp v[:,block_h[i][j][p]]+v[:,block_h[i][j][q]]+supp_h[i][:,z]]
                end
            end
        end
    end
      =#
    
    Usupp=unique(Usupp,dims=2)
    lUsupp=size(Usupp,2)

    
    return Usupp,lUsupp,block_g0,block_g,block_h,lblock_g0,lblock_g,lblock_h,lt_block_g0,lt_block_g,lt_block_h

end





function get_constant_trace_Cliq_mix_LMBM(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g,coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h,coe_h::Vector{Vector{Float64}},dg::Vector{Int64},dh::Vector{Int64},v::Matrix{UInt64},sk::UInt64,sk_g::Vector{UInt64},s2k_g::Vector{UInt64},s2k_h::Vector{UInt64};UseEq::Bool=true)
    
    lV=sk
    if m>0
        lV+=sum(@fastmath @inbounds sk_g[i]*lmon_g[i] for i in 1:m)
    end
    if UseEq
        if l>0
            lV+=sum(@fastmath @inbounds lmon_h[i]*s2k_h[i] for i in 1:l)
        end
    end


    V=Matrix{UInt64}(undef,n,lV)
    V[:,1:sk]=2*v[:,1:sk]
    t=sk+1
    @fastmath @inbounds for i in 1:m
        @fastmath @inbounds for r in 1:lmon_g[i]
            @fastmath @inbounds for j in 1:sk_g[i]
                V[:,t]=2*v[:,j]+supp_g[i][:,r]
                t+=1
            end
        end
    end
    if UseEq
        @fastmath @inbounds for i in 1:l
            @fastmath @inbounds for r in 1:lmon_h[i]
                @fastmath @inbounds for j in 1:s2k_h[i]
                    V[:,t]=v[:,j]+supp_h[i][:,r]
                    t+=1
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
                                    
    
                                    
    ak=value(lambda)
    
    #=println("  Computing constant trace status: ", termination_status(model))
    println("  Constant trace: ak = ",ak)=#
    
    P=sqrt.(value.(p0))
    Pg=[sqrt.(value.(p[i])) for i in 1:m]
    
    return ak,P,Pg
end


function model_Cliq_mix_LMBM(n::Int64,m::Int64,l::Int64,lmon_f::Int64,supp_f::SparseMatrixCSC{UInt64,Int64},coe_f::Vector{Float64},lmon_g::Vector{UInt64},supp_g,coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h,coe_h::Vector{Vector{Float64}},dg::Vector{Int64},dh::Vector{Int64},mex_cliq::Vector{UInt64},lmex_cliq::Int64,v::Matrix{UInt64},Usupp::Matrix{UInt64},lUsupp::UInt64,block_g0::Vector{Vector{UInt64}},block_g::Vector{Vector{Vector{UInt64}}},block_h::Vector{Vector{Vector{UInt64}}},lblock_g0::UInt64,lblock_g::Vector{UInt64},lblock_h::Vector{UInt64},lt_block_g0::Vector{UInt64},lt_block_g::Vector{Vector{UInt64}},lt_block_h::Vector{Vector{UInt64}},sk::UInt64,sk_g::Vector{UInt64},sk_h::Vector{UInt64},s2k::UInt64,s2k_g::Vector{UInt64},s2k_h::Vector{UInt64};get_constant_trace::Bool=false,UseEq::Bool=true)
    
   
    ak,P,Pg=get_constant_trace_Cliq_mix_LMBM(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,dg,dh,v,sk,sk_g,s2k_g,s2k_h,UseEq=UseEq)
                                                          
    
    omega=lblock_g0+sum(lblock_g)
    
    #=println("  Number of blocks: omega=",omega)
    if m>0
        println("  Size of the largest block: s^max=",maximum([maximum(lt_block_g0);maximum(maximum(lt_block_g[i]) for i in 1:m)]))
    else
        println("  Size of the largest block: s^max=",maximum(lt_block_g0))
    end=#
    
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
    
    sum_lUsupp_h=sum(lUsupp_h)
    zeta+=sum_lUsupp_h
    
    #println("  Number of equality trace constraints: zeta=",zeta)
    
    
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
   
            


    a0=zeros(Float64,d)
    
    
    @simd for p in 1:lmon_f
        I=IndM[Order(supp_f[:,p])][1]
        a0[invIndeM[I[1],I[2]]]=coe_f[p]*ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1)
   end
    
    
    a_mex_ind=Vector{UInt64}([])
    a_mex_val=Vector{Float64}([])
    
    
    
    
    
    for j in 1:lmex_cliq
        I=IndM[mex_cliq[j]][1]
        push!(a_mex_ind,invIndeM[I[1],I[2]])
        push!(a_mex_val,ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))  
    end
    
    
    
    s=Vector{Int64}(undef,omega)
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
    
   
    
    if get_constant_trace
        return omega,a0,a_ind1,a_ind2,a_val,a_len,a_mex_ind,a_mex_val,r_blo,u_vec,s,zeta,d,ak
    else
        return omega,a0,a_ind1,a_ind2,a_val,a_len,a_mex_ind,a_mex_val,r_blo,u_vec,s,zeta,d
    end
end
