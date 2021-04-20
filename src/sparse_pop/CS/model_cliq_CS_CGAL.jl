function get_constant_trace_Cliq_CS_CGAL(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Vector{Vector{UInt64}}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Vector{Vector{UInt64}}},coe_h::Vector{Vector{Float64}},dg::Vector{Int64},dh::Vector{Int64},v::Vector{Vector{UInt64}},k::Int64;use_eqcons_to_get_constant_trace::Bool=true)

    sk=numlet(n,k)
    sk_g=Vector{Int64}(undef,m)
    s2k_g=Vector{Int64}(undef,m)
    sk_h=Vector{Int64}(undef,l)
    s2k_h=Vector{Int64}(undef,l)
    
    
    
    ceil0=Int64(0)
    @fastmath @inbounds @simd for i in 1:m
        ceil0=ceil(Int64,dg[i]/2)
        sk_g[i]=numlet(n,k-ceil0)
        s2k_g[i]=numlet(n,2*(k-ceil0))
    end
                    
    @fastmath @inbounds @simd for i in 1:l
        ceil0=ceil(Int64,dh[i]/2)
        sk_h[i]=numlet(n,k-ceil0)
        s2k_h[i]=numlet(n,2*(k-ceil0))
    end
    
    
    lV=sk
    
    if m!=0
        lV+=sum(sk_g[i]*lmon_g[i] for i in 1:m)
    end
    if use_eqcons_to_get_constant_trace
        if l!=0
            lV+=sum(lmon_h[i]*s2k_h[i] for i in 1:l)
        end
    end
    V=[Vector{UInt64}([]) for j in 1:lV]
    V[1:sk]=[_sym_canon([v[i][end:-1:1]; v[i]]) for i in 1:sk]
    t=sk+1
    @fastmath @inbounds @simd for i in 1:m
        @fastmath @inbounds @simd for r in 1:lmon_g[i]
            @fastmath @inbounds @simd for j in 1:sk_g[i]
                V[t]=_sym_canon([v[j][end:-1:1]; supp_g[i][r]; v[j]])
                t+=1
            end
        end
    end
    if use_eqcons_to_get_constant_trace               
        @fastmath @inbounds @simd for i in 1:l
            @fastmath @inbounds @simd for r in 1:lmon_h[i]
                @fastmath @inbounds for j in 1:sk_h[i], w in j:sk_h[i]
                    V[t]=_sym_canon([v[j][end:-1:1]; supp_h[i][r]; v[w]])
                    t+=1
                end
            end
        end 
    end
    V=unique(V)
    V=sort(V)
    lV=length(V)
                                    
    model=JuMP.Model(with_optimizer(Mosek.Optimizer, QUIET=true))
    cons=[AffExpr(0) for i in 1:lV]

    p0=@variable(model, [1:sk], lower_bound=1)
    p=[@variable(model, [1:sk_g[i]], lower_bound=1) for i in 1:m]
    

    
    @fastmath @inbounds @simd for j in 1:sk
        add_to_expression!(cons[ncbfind(V,lV,_sym_canon([v[j][end:-1:1]; v[j]]))],p0[j])
    end


    @fastmath @inbounds @simd for i in 1:m
        @fastmath @inbounds @simd for j in 1:sk_g[i]
            @fastmath @inbounds @simd for r in 1:lmon_g[i]
                add_to_expression!(cons[ncbfind(V,lV,_sym_canon([v[j][end:-1:1]; supp_g[i][r]; v[j]]))],p[i][j]*coe_g[i][r])
            end
        end
    end
    
    if use_eqcons_to_get_constant_trace
        q=[@variable(model, [1:Int64(sk_h[i]*(sk_h[i]+1)/2)]) for i in 1:l]
        t_q=1
        @fastmath @inbounds @simd for i in 1:l
            @fastmath @inbounds for j in 1:sk_h[i], w in j:sk_h[i]
                @fastmath @inbounds @simd for r in 1:lmon_h[i]
                    add_to_expression!(cons[ncbfind(V,lV,_sym_canon([v[j][end:-1:1]; supp_h[i][r]; v[w]]))],p[i][t_q]*coe_h[i][r])
                    t_q+=1
                end
            end
            t_q=1
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


function model_Cliq_CS_CGAL(n::Int64,m::Int64,l::Int64,lmon_f::Int64,supp_f::Vector{Vector{UInt64}},coe_f::Vector{Float64},lmon_g::Vector{UInt64},supp_g::Vector{Vector{Vector{UInt64}}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Vector{Vector{UInt64}}},coe_h::Vector{Vector{Float64}},dg::Vector{Int64},dh::Vector{Int64},v::Vector{Vector{UInt64}},s2k::Int64,sort_v::Vector{Vector{UInt64}},lsort_v::UInt64,mex_cliq::Vector{UInt64},lmex_cliq::Int64,k::Int64;use_eqcons_to_get_constant_trace::Bool=true)
    
    sk,sk_g,s2k_g,s2k_h,ak,P,Pg=get_constant_trace_Cliq_CS_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,dg,dh,v,k,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)

                                                          
        
    Order(alpha)=ncbfind(sort_v,lsort_v,_sym_canon(alpha))
    #println("  Number of blocks: omega=",m+1)
    #println("  Size of the largest block: s^max=",sk)
    
    u=Int64(0.5*sk*(sk+1))
    u_g=[@fastmath @inbounds Int64(0.5*sk_g[i]*(sk_g[i]+1)) for i in 1:m]
    d=u+sum(u_g)
    
    #println("  Number of equality trace constraints: zeta=",zeta) 
            

    IndM=[@fastmath @inbounds Vector{Vector{UInt64}}([]) for j in 1:lsort_v]
    invIndeM=[@fastmath @inbounds Vector{UInt64}([]) for j in 1:lsort_v]
    r=UInt64(0)
    t_iter=UInt64(1)
    
    a_ind1=Vector{UInt64}([])
    a_ind2=Vector{UInt64}([])
    a_val=Vector{Float64}([])
    a_len=UInt64(0)
    
    a_ind1_=Vector{UInt64}([])
    a_val_=Vector{Float64}([])
    a_len_=UInt64(0)
    
            
    @fastmath @inbounds for i in 1:sk, j in 1:i
        @inbounds r=Order([v[i][end:-1:1];v[j]])
        @inbounds push!(IndM[r],[i;j])
        @inbounds push!(invIndeM[i],t_iter)
        @inbounds t_iter+=1
    end
    l_IndM=[@fastmath @inbounds length(IndM[r]) for r in 1:lsort_v]
    
         
    t_a=UInt64(1)
    I=zeros(UInt64,2)
    I_=zeros(UInt64,2)
            
    @fastmath @inbounds @simd for r in 1:lsort_v
        if l_IndM[r]>1
            @fastmath @inbounds @simd for i in 2:l_IndM[r]
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
   
    
    sort_v_g=Vector{Vector{Vector{UInt64}}}(undef,m)
    lsort_v_g=zeros(UInt64,m)
    
    for j in 1:m
        sort_v_g[j]=[_sym_canon(v[i]) for i in 1:s2k_g[j]]
        sort_v_g[j]=unique(sort_v_g[j])  
        sort_v_g[j]=sort(sort_v_g[j])
        lsort_v_g[j]=length(sort_v_g[j])
    end
    
    
    IndMg=[@fastmath @inbounds [@fastmath @inbounds Vector{Vector{UInt64}}([]) for j in 1:lsort_v_g[i]] for i in 1:m]
    invIndeMg=[@fastmath @inbounds [@fastmath @inbounds Vector{UInt64}([]) for j in 1:lsort_v_g[i]] for i in 1:m]
    t_Blo=u
    
    vec=spzeros(Float64,d)
 
    @fastmath @inbounds @simd for j in 1:m
        r=UInt64(0)
        t_iter=UInt64(1)
        @fastmath @inbounds for p in 1:sk_g[j], q in 1:p
            r=ncbfind(sort_v_g[j],lsort_v_g[j],_sym_canon([v[p][end:-1:1];v[q]]))
            push!(IndMg[j][r],[p;q])
            push!(invIndeMg[j][p],t_iter)
            t_iter+=1
        end
                
        
        l_IndM=[@fastmath @inbounds length(IndMg[j][q]) for q in 1:lsort_v_g[j]]
       
        @fastmath @inbounds @simd for r in 1:lsort_v_g[j]
            if l_IndM[r]>1
                @fastmath @inbounds @simd for i in 2:l_IndM[r]
                    I1,I2=IndMg[j][r][1]
                    push!(a_ind1,invIndeMg[j][I1][I2]+t_Blo)
                    push!(a_ind2,t_a)
                    push!(a_val,ak/Pg[j][I1]/Pg[j][I2]*((0.5*sqrt(2)-1)*(I2<I1)+1))
                    a_len+=1
                  
             
                    I1,I2=IndMg[j][r][i]
                    push!(a_ind1,invIndeMg[j][I1][I2]+t_Blo)
                    push!(a_ind2,t_a)
                    push!(a_val,-ak/Pg[j][I1]/Pg[j][I2]*((0.5*sqrt(2)-1)*(I2<I1)+1))
                    a_len+=1 
                    t_a+=1
                end
                IndMg[j][r]=Vector{Int64}[IndMg[j][r][1]]

            end
            I1_,I2_=IndMg[j][r][1]
            push!(a_ind1,invIndeMg[j][I1_][I2_]+t_Blo)
            push!(a_ind2,t_a)
            push!(a_val,-ak/Pg[j][I1_]/Pg[j][I2_]*((0.5*sqrt(2)-1)*(I2_<I1_)+1))
            a_len+=1

            
           @fastmath @inbounds @simd  for p in 1:lmon_g[j]  
                I1,I2=IndM[Order([v[I1_][end:-1:1];supp_g[j][p];v[I2_]])][end]
                vec[invIndeM[I1][I2]]+=coe_g[j][p]*ak/P[I1]/P[I2]*((0.5*sqrt(2)-1)*(I2<I1)+1)
           end
           a_ind1_,a_val_=findnz(vec)
            a_len_=length(a_ind1_)
            a_ind1=[a_ind1;a_ind1_]
            a_val=[a_val;a_val_]
            a_ind2=[a_ind2;t_a*ones(UInt64,a_len_)]
            a_len+=a_len_
            @inbounds t_a+=1
            vec=spzeros(Float64,d)
         end
        
         t_Blo+=u_g[j]
    end
   
    
    sort_v_h=Vector{Vector{Vector{UInt64}}}(undef,l)
    lsort_v_h=zeros(UInt64,l)
    
    for j in 1:l
        sort_v_h[j]=[_sym_canon(v[i]) for i in 1:s2k_h[j]]
        sort_v_h[j]=unique(sort_v_h[j])            
        sort_v_h[j]=sort(sort_v_h[j])
        lsort_v_h[j]=length(sort_v_h[j])
    end       

    @fastmath @inbounds @simd for j in 1:l
        @fastmath @inbounds for r in 1:lsort_v_h[j]
            I1_,I2_=IndM[Order(sort_v_h[j][r])][1]
            @fastmath @inbounds @simd  for p in 1:lmon_h[j]
            I1,I2=IndM[Order([v[I1_][end:-1:1];supp_h[j][p];v[I2_]])][1]
                vec[invIndeM[I1][I2]]+=coe_h[j][p]*ak/P[I1]/P[I2]*((0.5*sqrt(2)-1)*(I2<I1)+1)
                   end
                    a_ind1_,a_val_=findnz(vec)
                    a_len_=length(a_ind1_)
                    a_ind1=[a_ind1;a_ind1_]
                    a_val=[a_val;a_val_]
                    a_ind2=[a_ind2;t_a*ones(UInt64,a_len_)]
                    a_len+=a_len_
                    @inbounds t_a+=1
                    vec=spzeros(Float64,d)
               end       
    end
    
    
    zeta=t_a-1
    #println("  Number of equality trace constraints: zeta=",zeta)
    


    a0=zeros(Float64,d)
    
    
    
    
    @fastmath @inbounds @simd  for p in 1:lmon_f
        I1,I2=IndM[Order(supp_f[p])][1]
        a0[invIndeM[I1][I2]]+=coe_f[p]*ak/P[I1]/P[I2]*((0.5*sqrt(2)-1)*(I2<I1)+1)
   end
    
    
    
    a_mex_ind=Vector{UInt64}([])
    a_mex_val=Vector{Float64}([])
    
    
    
    
    @fastmath @inbounds @simd for j in 1:lmex_cliq
        I=IndM[mex_cliq[j]][1]
        push!(a_mex_ind,invIndeM[I[1]][I[2]])
        push!(a_mex_val,ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1))
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
    @fastmath @inbounds @simd for j in 1:m
        uvec[t]=u_g[j]
        r_blo[t]=t_Blo
        s[t]=sk_g[j]
        t+=1
        t_Blo+=u_g[j]
    end
    return omega,a0,a_ind1,a_ind2,a_val,a_len,a_mex_ind,a_mex_val,r_blo,uvec,s,zeta,d,ak
   
end
