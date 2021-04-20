function initial_mix_CGAL(n::Int64,m::Int64,l::Int64,k::Int64,dg::Vector{Int64},dh::Vector{Int64})

    v=get_ncbasis(n,2*k)
    s2k=length(v)
    
    sk=numlet(n,k)
    sk_g=Vector{UInt64}(undef,m)
    sk_h=Vector{UInt64}(undef,l)
    s2k_g=Vector{UInt64}(undef,m)
    s2k_h=Vector{UInt64}(undef,l)
    
        
    
    ceil0=Int64(0)
    for i in 1:m
        ceil0=ceil(Int64,dg[i]/2)
        sk_g[i]=numlet(n,k-ceil0)
        s2k_g[i]=numlet(n,2*(k-ceil0))
    end
                    
    for i in 1:l
        ceil0=ceil(Int64,dh[i]/2)
        sk_h[i]=numlet(n,k-ceil0)
        s2k_h[i]=numlet(n,2*(k-ceil(Int64,dh[i]/2)))
    end
    return v,s2k,sk,sk_g,sk_h,s2k_g,s2k_h
end


function get_blocks_Cliq_mix_CGAL(k::Int64,n::Int64,m::Int64,l::Int64,Usupp::Vector{Vector{UInt64}},lmon_g::Vector{UInt64},lmon_h::Vector{UInt64},supp_g::Vector{Vector{Vector{UInt64}}},supp_h::Vector{Vector{Vector{UInt64}}},coe_g::Vector{Vector{Float64}},coe_h::Vector{Vector{Float64}},v::Vector{Vector{UInt64}},s2k::UInt64,sk::UInt64,sk_g::Vector{UInt64},sk_h::Vector{UInt64},s2k_g::Vector{UInt64},s2k_h::Vector{UInt64})
    
    #println(typeof(Usupp))
    Usupp=_sym_canon.(Usupp)
    Usupp=unique(Usupp)
    Usupp=sort(Usupp)
    lUsupp=length(Usupp)


    block_g=Vector{Vector{Vector{UInt64}}}(undef,m)
    block_h=Vector{Vector{Vector{UInt64}}}(undef,l)


    lblock_g=Vector{UInt64}(undef,m)
    lblock_h=Vector{UInt64}(undef,l)


    lt_block_g=Vector{Vector{UInt64}}(undef,m)
    lt_block_h=Vector{Vector{UInt64}}(undef,l)

    
    


    
    
    graph=SimpleGraph(sk)
    for p in 1:sk, q in 1:p
        if ncbfind(Usupp,lUsupp,_sym_canon([v[p][end:-1:1];v[q]]))!=0
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
                    if ncbfind(Usupp,lUsupp,_sym_canon([v[p][end:-1:1];supp_g[i][y];v[q]]))!=0
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
                    if ncbfind(Usupp,lUsupp,_sym_canon([v[p][end:-1:1];supp_h[i][y];v[q]]))!=0
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
    
    lUsupp=Int64(0.5*sum(lt_block_g0[j]*(lt_block_g0[j]+1) for j in 1:lblock_g0))
    if m>0
        lUsupp+=Int64(0.5*sum(lmon_g[i]*lt_block_g[i][j]*(lt_block_g[i][j]+1) for i in 1:m for j in 1:lblock_g[i]))
    end
    if l>0
        lUsupp+=Int64(0.5*sum(lmon_h[i]*lt_block_h[i][j]*(lt_block_h[i][j]+1) for i in 1:l for j in 1:lblock_h[i]))
    end
    Usupp=Vector{Vector{UInt64}}(undef,lUsupp)

    t_supp=1
    for j in 1:lblock_g0
        for p in 1:lt_block_g0[j], q in 1:p
            Usupp[t_supp]=_sym_canon([v[block_g0[j][p]][end:-1:1];v[block_g0[j][q]]])
            t_supp+=1
        end
    end
    
    
   for i in 1:m
        for j in 1:lblock_g[i]
            for p in 1:lt_block_g[i][j], q in 1:p
                for z in 1:lmon_g[i]
                    Usupp[t_supp]=_sym_canon([v[block_g[i][j][p]][end:-1:1];supp_g[i][z];v[block_g[i][j][q]]])
                    t_supp+=1
                end
            end
        end
    end
    for i in 1:l
        for j in 1:lblock_h[i]
            for p in 1:lt_block_h[i][j], q in 1:p
                for z in 1:lmon_h[i]
                    Usupp[t_supp]=_sym_canon([v[block_h[i][j][p]][end:-1:1];supp_h[i][z];v[block_h[i][j][q]]])
                    t_supp+=1
                end
            end
        end
    end
     
    
    Usupp=unique(Usupp)
    lUsupp=length(Usupp)

    
    return Usupp,lUsupp,block_g0,block_g,block_h,lblock_g0,lblock_g,lblock_h,lt_block_g0,lt_block_g,lt_block_h

end





function get_constant_trace_Cliq_mix_CGAL(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Vector{Vector{UInt64}}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Vector{Vector{UInt64}}},coe_h::Vector{Vector{Float64}},dg::Vector{Int64},dh::Vector{Int64},v::Vector{Vector{UInt64}},sk::UInt64,sk_g::Vector{UInt64},s2k_g::Vector{UInt64},s2k_h::Vector{UInt64};use_eqcons_to_get_constant_trace::Bool=true)
    
    
    
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
    
    return ak,P,Pg
end


function model_Cliq_mix_CGAL(n::Int64,m::Int64,l::Int64,lmon_f::Int64,supp_f::Vector{Vector{UInt64}},coe_f::Vector{Float64},lmon_g::Vector{UInt64},supp_g::Vector{Vector{Vector{UInt64}}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Vector{Vector{UInt64}}},coe_h::Vector{Vector{Float64}},dg::Vector{Int64},dh::Vector{Int64},mex_cliq::Vector{UInt64},lmex_cliq::Int64,v::Vector{Vector{UInt64}},Usupp::Vector{Vector{UInt64}},lUsupp::UInt64,block_g0::Vector{Vector{UInt64}},block_g::Vector{Vector{Vector{UInt64}}},block_h::Vector{Vector{Vector{UInt64}}},lblock_g0::UInt64,lblock_g::Vector{UInt64},lblock_h::Vector{UInt64},lt_block_g0::Vector{UInt64},lt_block_g::Vector{Vector{UInt64}},lt_block_h::Vector{Vector{UInt64}},sk::UInt64,sk_g::Vector{UInt64},sk_h::Vector{UInt64},s2k::UInt64,s2k_g::Vector{UInt64},s2k_h::Vector{UInt64};use_eqcons_to_get_constant_trace::Bool=true)
    
   
    ak,P,Pg=get_constant_trace_Cliq_mix_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,dg,dh,v,sk,sk_g,s2k_g,s2k_h,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)
                      
  

    omega=lblock_g0+sum(lblock_g)
    
    #=println("  Number of blocks: omega=",omega)
    if m>0
        println("  Size of the largest block: s^max=",maximum([maximum(lt_block_g0);maximum(maximum(lt_block_g[i]) for i in 1:m)]))
    else
        println("  Size of the largest block: s^max=",maximum(lt_block_g0))
    end=#
    
    
    
    
    
    w=[Int64(0.5*lt_block_g0[j]*(lt_block_g0[j]+1)) for j in 1:lblock_g0]
    w_g=[[@inbounds Int64(0.5*lt_block_g[i][j]*(lt_block_g[i][j]+1)) for j in 1:lblock_g[i]] for i in 1:m]
    
    u=sum(w)
    u_g=[sum(w_g[i]) for i in 1:m]

    
    
    
    d=u+sum(u_g)
    
    
    
    Usupp=Vector{Vector{UInt64}}(undef,u)
    t_iter=1
    for j in 1:lblock_g0
        for p in 1:lt_block_g0[j], q in 1:p
            @inbounds Usupp[t_iter]=_sym_canon([v[block_g0[j][p]][end:-1:1];v[block_g0[j][q]]])
            t_iter+=1
        end
    end
    Usupp=unique(Usupp)
    Usupp=sort(Usupp)
    lUsupp=length(Usupp) 
   
    Order(alpha)=ncbfind(Usupp,lUsupp,_sym_canon(alpha))
    

    IndM=[@inbounds Vector{Vector{Int64}}([]) for j in 1:lUsupp]
    invIndeM=spzeros(UInt64,sk,sk)
    r=UInt64(0)
    t_iter=UInt64(1)
    
    a_ind1=Vector{UInt64}([])
    a_ind2=Vector{UInt64}([])
    a_val=Vector{Float64}([])
    a_len=UInt64(0)
                                    
    a_ind1_=Vector{UInt64}([])
    a_val_=Vector{Float64}([])
    a_len_=UInt64(0)                                
    
    t_blo=0
    
    for j in 1:lblock_g0
        for p in 1:lt_block_g0[j], q in 1:p
            @inbounds r=Order([v[block_g0[j][p]][end:-1:1];v[block_g0[j][q]]])
            @inbounds push!(IndM[r],[block_g0[j][p],block_g0[j][q]])
            @inbounds invIndeM[block_g0[j][p],block_g0[j][q]]=t_iter
            @inbounds t_iter+=1
        end
        t_blo+=lt_block_g0[j]
    end
    l_IndM=[length(IndM[r]) for r in 1:lUsupp]
   
    
    
    t_a=1
    
    
    I=zeros(UInt64,2)
    I_=zeros(UInt64,2)
            
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
    
    
    Usupp_g=Vector{Vector{Vector{UInt64}}}(undef,m)
    lUsupp_g=zeros(UInt64,m)

    
    for i in 1:m     
        Usupp_g[i]=Vector{Vector{UInt64}}(undef,u_g[i])
        t_iter=1
        for j in 1:lblock_g[i]
            for p in 1:lt_block_g[i][j], q in p:lt_block_g[i][j]
                @inbounds Usupp_g[i][t_iter]=_sym_canon([v[block_g[i][j][p]];v[block_g[i][j][q]]])
                t_iter+=1
            end
        end
        Usupp_g[i]=unique(Usupp_g[i])
        Usupp_g[i]=sort(Usupp_g[i])
        
        lUsupp_g[i]=length(Usupp_g[i]) 
    end
   
    IndMg=[[Vector{Vector{Int64}}([]) for j in 1:lUsupp_g[i]] for i in 1:m]
    invIndeMg=[spzeros(UInt64,sk_g[i],sk_g[i]) for i in 1:m]
    t_Blo=u
    vec=spzeros(Float64,d)
 
    for i in 1:m
        r=UInt64(0)
        t_iter=UInt64(1)

        t_blo=0
        for j in 1:lblock_g[i]
            for p in 1:lt_block_g[i][j], q in 1:p
                @inbounds r=ncbfind(Usupp_g[i],lUsupp_g[i],_sym_canon([v[block_g[i][j][p]][end:-1:1];v[block_g[i][j][q]]]))
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
            
            
           
            I_=IndMg[i][rr][1]
            push!(a_ind1,invIndeMg[i][I_[1],I_[2]]+t_Blo)
            push!(a_ind2,t_a)
            push!(a_val,-ak/Pg[i][I_[1]]/Pg[i][I_[2]]*((0.5*sqrt(2)-1)*(I_[2]<I_[1])+1))
            a_len+=1
            
            for p in 1:lmon_g[i]  
                I=IndM[Order([v[I_[1]][end:-1:1];supp_g[i][p];v[I_[2]]])][end]
                vec[invIndeM[I[1],I[2]]]+=coe_g[i][p]*ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1)
                
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
         t_Blo+=u_g[i]
    end
    
    Usupp_h=Vector{Vector{Vector{UInt64}}}(undef,l)
    lUsupp_h=zeros(UInt64,l)
    
    
    
    for i in 1:l
        Usupp_h[i]=Vector{Vector{UInt64}}(undef,Int(0.5*sum(lt_block_h[i][j]*(lt_block_h[i][j]+1) for j in 1:lblock_h[i])))
        t_iter=1
        for j in 1:lblock_h[i]
            for p in 1:lt_block_h[i][j], q in p:lt_block_h[i][j]
                @inbounds Usupp_h[i][t_iter]=_sym_canon([v[block_h[i][j][p]][end:-1:1];v[block_h[i][j][q]]])
                t_iter+=1
            end
        end
        Usupp_h[i]=unique(Usupp_h[i])
        Usupp_h[i]=sort(Usupp_h[i])
        lUsupp_h[i]=length(Usupp_h[i])
    end
    
    sum_lUsupp_h=sum(lUsupp_h)
    
    
    #println("  Number of equality trace constraints: zeta=",zeta)
    
    
    @simd for j in 1:l
        @simd for r in 1:lUsupp_h[j]
            I_=IndM[Order(Usupp_h[j][r])][1]
            @simd for p in 1:lmon_h[j]
                      I=IndM[Order([v[I_[1]][end:-1:1];supp_h[j][p];v[I_[2]]])][1]
                      vec[invIndeM[I[1],I[2]]]+=coe_h[j][p]*ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1)
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


    a0=zeros(Float64,d)
    
    
    @simd for p in 1:lmon_f
        I=IndM[Order(supp_f[p])][1]
        a0[invIndeM[I[1],I[2]]]+=coe_f[p]*ak/P[I[1]]/P[I[2]]*((0.5*sqrt(2)-1)*(I[2]<I[1])+1)
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
    
   
    
    return omega,a0,a_ind1,a_ind2,a_val,a_len,a_mex_ind,a_mex_val,r_blo,u_vec,s,zeta,d,ak
    
end
