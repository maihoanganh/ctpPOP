function get_constant_trace_Cliq_CS_LMBM(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,dg,dh,v,k;UseEq=true)
    
    sk=binomial(k+n,n)
    sk_g=Vector{Int64}(undef,m)
    s2k_g=Vector{Int64}(undef,m)
    s2k_h=Vector{Int64}(undef,l)
    
        
    
    ceil0=Int64(0)
    for i in 1:m
        ceil0=ceil(Int64,dg[i]/2)
        sk_g[i]=binomial(k-ceil0+n,n)
        s2k_g[i]=binomial(2*(k-ceil0)+n,n)
    end
                    
    for i in 1:l
        s2k_h[i]=binomial(2*(k-ceil(Int64,dh[i]/2))+n,n)
    end
    
    
    lV=sk
    if m>0
        lV+=sum(sk_g[i]*lmon_g[i] for i in 1:m)
    end
    if UseEq
        if l>0
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
    
    

    
    for j in 1:sk
        @inbounds add_to_expression!(cons[bfind(V,lV,2*v[:,j],n)],p0[j])
    end

    p=Vector{Vector{VariableRef}}(undef, m)
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
                                    
    #println("  Computing constant trace status: ", termination_status(model))
                                    
    ak=value(lambda)
    #println("  Constant trace: ak = ",ak)
    
    P=sqrt.(value.(p0))
    Pg=[sqrt.(value.(p[i])) for i in 1:m]
    
    return sk,sk_g,s2k_g,s2k_h,ak,P,Pg
end


function model_Cliq_CS_LMBM(n,m,l,lmon_f,supp_f,coe_f,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,dg,dh,v,s2k,sort_v,re_ind,mex_cliq,lmex_cliq,k;get_constant_trace=false,UseEq=true)
    
    sk,sk_g,s2k_g,s2k_h,ak,P,Pg=get_constant_trace_Cliq_CS_LMBM(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,dg,dh,v,k,UseEq=UseEq)
                                                          
        
    Order(alpha)=re_ind[bfind(sort_v,s2k,alpha,n)]        
    #println("  Number of blocks: omega=",m+1)
    #println("  Size of the largest block: s^max=",sk)
    
    u=Int64(0.5*sk*(sk+1))
    u_g=[@inbounds Int64(0.5*sk_g[i]*(sk_g[i]+1)) for i in 1:m]
    d=u+sum(u_g)
    zeta=d-s2k+sum(s2k_h)
    #println("  Number of equality trace constraints: zeta=",zeta) 
            

    IndM=[@inbounds Vector{Vector{UInt64}}([]) for j in 1:s2k]
    invIndeM=[@inbounds Vector{UInt64}([]) for j in 1:sk]
    r=UInt64(0)
    t_iter=UInt64(1)
    
    

    a_ind1=Vector{UInt64}([])
    a_ind2=Vector{UInt64}([])
    a_val=Vector{Float64}([])
    a_len=UInt64(0)
    
            
    for i in 1:sk, j in 1:i
        @inbounds r=Order(v[:,i]+v[:,j])
        @inbounds push!(IndM[r],[i;j])
        @inbounds push!(invIndeM[i],t_iter)
        @inbounds t_iter+=1
    end
    l_IndM=[length(IndM[r]) for r in 1:s2k]
    
         
    t_a=UInt64(1)
    I=zeros(UInt64,2)
            
    for r in 1:s2k
        if l_IndM[r]>1
            for i in 2:l_IndM[r]
                I=IndM[r][1]
                push!(a_ind1,invIndeM[I[1]][I[2]])
                push!(a_ind2,t_a)
                push!(a_val,ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
                a_len+=1
                
                
                I=IndM[r][i]
                push!(a_ind1,invIndeM[I[1]][I[2]])
                push!(a_ind2,t_a)
                push!(a_val,-ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
               
                a_len+=1
                t_a+=1
             end
             IndM[r]=Vector{Int64}[IndM[r][1]]
         end
    end
   
    
    
    
    
    IndMg=[[Vector{Vector{UInt64}}([]) for j in 1:s2k_g[i]] for i in 1:m]
    invIndeMg=[[Vector{UInt64}([]) for j in 1:sk_g[i]] for i in 1:m]
    t_Blo=u
 
    for j in 1:m
        r=UInt64(0)
        t_iter=UInt64(1)
        for p in 1:sk_g[j], q in 1:p
            @inbounds r=Order(v[:,p]+v[:,q])
            @inbounds push!(IndMg[j][r],[p;q])
            @inbounds push!(invIndeMg[j][p],t_iter)
            t_iter+=1
        end
                
        
        l_IndM=[length(IndMg[j][q]) for q in 1:s2k_g[j]]

        for r in 1:s2k_g[j]
            if l_IndM[r]>1
                for i in 2:l_IndM[r]
                    I=IndMg[j][r][1]
                    push!(a_ind1,invIndeMg[j][I[1]][I[2]]+t_Blo)
                    push!(a_ind2,t_a)
                    push!(a_val,ak/Pg[j][I[1]]/Pg[j][I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
                    a_len+=1
                    
                    I=IndMg[j][r][i]
                    push!(a_ind1,invIndeMg[j][I[1]][I[2]]+t_Blo)
                    push!(a_ind2,t_a)
                    push!(a_val,-ak/Pg[j][I[1]]/Pg[j][I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
                    a_len+=1
                    
                    t_a+=1
                end
                IndMg[j][r]=Vector{Int64}[IndMg[j][r][1]]
            end
            
            @inbounds I=IndMg[j][r][1]
                    push!(a_ind1,invIndeMg[j][I[1]][I[2]]+t_Blo)
                    push!(a_ind2,t_a)
                    push!(a_val,-ak/Pg[j][I[1]]/Pg[j][I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
                    a_len+=1
                  
            @simd for p in 1:lmon_g[j]  
                @inbounds I=IndM[Order(supp_g[j][:,p]+v[:,r])][end]
                    push!(a_ind1,invIndeM[I[1]][I[2]])
                    push!(a_ind2,t_a)
                    push!(a_val,coe_g[j][p]*ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
                    a_len+=1
                     
                   end
                   t_a+=1
            
         end
        
         
        
         t_Blo+=u_g[j]
    end
   
    
           

    @simd for j in 1:l
        @simd for r in 1:s2k_h[j]
            @simd for p in 1:lmon_h[j]
                      @inbounds I=IndM[Order(supp_h[j][:,p]+v[:,r])][1]
                      push!(a_ind1,invIndeM[I[1]][I[2]])
                      push!(a_ind2,t_a)
                      push!(a_val,coe_h[j][p]*ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
                      a_len+=1
                     
                   end
                    @inbounds t_a+=1
               end       
    end 
    


    a0=zeros(Float64,d)
    
    
    
    
    @simd for p in 1:lmon_f
        @inbounds I=IndM[Order(supp_f[:,p])][1]
        @inbounds a0[invIndeM[I[1]][I[2]]]=coe_f[p]*ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1)
      
   end
    
    
    
    a_mex_ind=Vector{UInt64}([])
    a_mex_val=Vector{Float64}([])
    
    
    
    
    for j in 1:lmex_cliq
        @inbounds I=IndM[mex_cliq[j]][1]
        @inbounds push!(a_mex_ind,invIndeM[I[1]][I[2]])
        @inbounds push!(a_mex_val,ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
    end
        

    
    
    
    omega=m+1

    r_blo=Vector{UInt64}(undef,omega)
    uvec=Vector{UInt64}(undef,omega)
    s=Vector{UInt64}(undef,omega)
    
    
    t=1
    r_blo[t]=0
    s[t]=sk
    uvec[t]=u
    t+=1
    t_Blo=u
    for j in 1:m
        uvec[t]=u_g[j]
        r_blo[t]=t_Blo
        s[t]=sk_g[j]
        t+=1
        t_Blo+=u_g[j]
    end
    if get_constant_trace
        return omega,a0,a_ind1,a_ind2,a_val,a_len,a_mex_ind,a_mex_val,r_blo,uvec,s,zeta,d,ak
    else
        return omega,a0,a_ind1,a_ind2,a_val,a_len,a_mex_ind,a_mex_val,r_blo,uvec,s,zeta,d
    end
end
