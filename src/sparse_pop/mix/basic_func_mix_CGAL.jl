function run_get_blocks_Cliq_mix_CGAL(n::Int64,p::Int64,t::Int64,k::Int64,lI::Vector{Int64},lJ::Vector{Int64},lW::Vector{Int64},I::Vector{Vector{UInt16}},J::Vector{Vector{UInt64}},W::Vector{Vector{UInt64}},Supp::Vector{Vector{UInt64}},IndA::Vector{Vector{UInt64}},lmon_g::Vector{UInt64},lmon_h::Vector{UInt64},supp_g::Vector{Vector{Vector{UInt64}}},supp_h::Vector{Vector{Vector{UInt64}}},coe_g::Vector{Vector{Float64}},coe_h::Vector{Vector{Float64}},dg::Vector{Int64},dh::Vector{Int64})
    
    
    Usupp=Vector{Vector{Vector{UInt64}}}(undef,p)
    lUsupp=Vector{UInt64}(undef,p)
    
    block_g0=Vector{Vector{Vector{UInt64}}}(undef,p)
    block_g=Vector{Vector{Vector{Vector{UInt64}}}}(undef,p)
    block_h=Vector{Vector{Vector{Vector{UInt64}}}}(undef,p)

    lblock_g0=Vector{UInt64}(undef,p)
    lblock_g=Vector{Vector{UInt64}}(undef,p)
    lblock_h=Vector{Vector{UInt64}}(undef,p)

    lt_block_g0=Vector{Vector{UInt64}}(undef,p)
    lt_block_g=Vector{Vector{Vector{UInt64}}}(undef,p)
    lt_block_h=Vector{Vector{Vector{UInt64}}}(undef,p)
    
    v=Vector{Vector{Vector{UInt64}}}(undef,p)
    s2k=Vector{UInt64}(undef,p)
    sk=Vector{UInt64}(undef,p)
    sk_g=Vector{Vector{UInt64}}(undef,p)
    sk_h=Vector{Vector{UInt64}}(undef,p)
    s2k_g=Vector{Vector{UInt64}}(undef,p)
    s2k_h=Vector{Vector{UInt64}}(undef,p)
    
    
    
    for j in 1:p
        v[j],s2k[j],sk[j],sk_g[j],sk_h[j],s2k_g[j],s2k_h[j]=initial_mix_CGAL(lI[j],lJ[j],lW[j],k,dg[J[j]],dh[W[j]])
        v[j]=[I[j][v[j][i]] for i in 1:s2k[j]]
         

        Usupp[j],lUsupp[j],block_g0[j],block_g[j],block_h[j],lblock_g0[j],lblock_g[j],lblock_h[j],lt_block_g0[j],lt_block_g[j],lt_block_h[j]=get_blocks_Cliq_mix_CGAL(k,lI[j],lJ[j],lW[j],[Supp[IndA[j]];[[v[j][i][end:-1:1]; v[j][i]] for i in 1:sk[j]]],lmon_g[J[j]],lmon_h[W[j]],[supp_g[i] for i in J[j]],[supp_h[i] for i in W[j]],coe_g[J[j]],coe_h[W[j]],v[j],s2k[j],sk[j],sk_g[j],sk_h[j],s2k_g[j],s2k_h[j])
        Usupp[j]=unique(Usupp[j])
        Usupp[j]=sort(Usupp[j])
        lUsupp[j]=length(Usupp[j])
    end
    
    
    lSupp=sum(lUsupp)
    Supp=Vector{Vector{UInt64}}(undef,lSupp)
    t_blo=0
    for j in 1:p
        Supp[1+t_blo:lUsupp[j]+t_blo]=Usupp[j]
        t_blo+=lUsupp[j]
    end
    
    Supp=unique(Supp)
    lSupp=length(Supp)
    
    IndA=[UInt64[] for i in 1:p]
    ind=zeros(UInt64,1)
    for j in 1:lSupp
        ind=unique(Supp[j])
        for i in 1:p
            if issubset(ind,I[i])
                push!(IndA[i],j)
            end
        end
    end
    
    
    
    for iter in 2:t
    
        for j in 1:p

            #println("Clique ",j,"th: ==================")
            Usupp[j],lUsupp[j],block_g0[j],block_g[j],block_h[j],lblock_g0[j],lblock_g[j],lblock_h[j],lt_block_g0[j],lt_block_g[j],lt_block_h[j]=get_blocks_Cliq_mix_CGAL(k,lI[j],lJ[j],lW[j],Supp[IndA[j]],lmon_g[J[j]],lmon_h[W[j]],[supp_g[i] for i in J[j]],[supp_h[i] for i in W[j]],coe_g[J[j]],coe_h[W[j]],v[j],s2k[j],sk[j],sk_g[j],sk_h[j],s2k_g[j],s2k_h[j])
            if iter==t
                Usupp[j]=unique(Usupp[j])
                Usupp[j]=sort(Usupp[j])
                lUsupp[j]=length(Usupp[j])
            end
        end
        
        if iter<t
            
            lSupp=sum(lUsupp)
            Supp=Vector{Vector{UInt64}}(undef,lSupp)
            t_blo=0
            for j in 1:p
                Supp[1+t_blo:lUsupp[j]+t_blo]=Usupp[j]
                t_blo+=lUsupp[j]
            end
        
            Supp=unique(Supp)
            lSupp=length(Supp)

            IndA=[UInt64[] for i in 1:p]
            for j in 1:lSupp
                ind=unique(Supp[j])
                for i in 1:p
                    if issubset(ind,I[i])
                        push!(IndA[i],j)
                    end
                end
            end
    #println("==================")
      
        end
    end
    
    return v,Usupp,lUsupp,block_g0,block_g,block_h,lblock_g0,lblock_g,lblock_h,lt_block_g0,lt_block_g,lt_block_h,sk,sk_g,sk_h,s2k,s2k_g,s2k_h
end


function run_model_cliq_mix_CGAL(p::Int64,lI::Vector{Int64},lJ::Vector{Int64},lW::Vector{Int64},I::Vector{Vector{UInt16}},J::Vector{Vector{UInt64}},W::Vector{Vector{UInt64}},lIndf::Vector{Int64},Indf::Vector{Vector{UInt64}},supp_f::Vector{Vector{UInt64}},coe_f,lmon_g::Vector{UInt64},supp_g::Vector{Vector{Vector{UInt64}}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Vector{Vector{UInt64}}},coe_h::Vector{Vector{Float64}},dg::Vector{Int64},dh::Vector{Int64},mex_cliq::Vector{Vector{UInt64}},lmex_cliq::Vector{Int64},v::Vector{Vector{Vector{UInt64}}},Usupp::Vector{Vector{Vector{UInt64}}},lUsupp::Vector{UInt64},block_g0::Vector{Vector{Vector{UInt64}}},block_g::Vector{Vector{Vector{Vector{UInt64}}}},block_h::Vector{Vector{Vector{Vector{UInt64}}}},lblock_g0::Vector{UInt64},lblock_g::Vector{Vector{UInt64}},lblock_h::Vector{Vector{UInt64}},lt_block_g0::Vector{Vector{UInt64}},lt_block_g::Vector{Vector{Vector{UInt64}}},lt_block_h::Vector{Vector{Vector{UInt64}}},sk::Vector{UInt64},sk_g::Vector{Vector{UInt64}},sk_h::Vector{Vector{UInt64}},s2k::Vector{UInt64},s2k_g::Vector{Vector{UInt64}},s2k_h::Vector{Vector{UInt64}};use_eqcons_to_get_constant_trace::Bool=true)
    
   
    
    omega_cliq=Vector{Int64}(undef,p)
    a0_cliq=Vector{Vector{Float64}}(undef,p)
    a_ind1_cliq=Vector{Vector{UInt64}}(undef,p)
    a_ind2_cliq=Vector{Vector{UInt64}}(undef,p)
    a_val_cliq=Vector{Vector{Float64}}(undef,p)
    a_len_cliq=Vector{UInt64}(undef,p)
    ak=Vector{Float64}(undef,p)
    
    a_mex_ind=Vector{Vector{UInt64}}(undef,p)
    a_mex_val=Vector{Vector{Float64}}(undef,p)
    r_cliq=Vector{Vector{UInt64}}(undef,p)
    u_cliq=Vector{Vector{UInt64}}(undef,p)
    s_cliq=Vector{Vector{UInt64}}(undef,p)
    
    zeta_cliq=Vector{UInt64}(undef,p)
    d_cliq=Vector{UInt64}(undef,p)
    
    

    for j in 1:p
        #println("Clique ",j,"th: ==================")

        omega_cliq[j],a0_cliq[j],a_ind1_cliq[j],a_ind2_cliq[j],a_val_cliq[j],a_len_cliq[j],a_mex_ind[j],a_mex_val[j],r_cliq[j],u_cliq[j],s_cliq[j],zeta_cliq[j],d_cliq[j],ak[j]=model_Cliq_mix_CGAL(lI[j],lJ[j],lW[j],lIndf[j],supp_f[Indf[j]],coe_f[Indf[j]],lmon_g[J[j]],[supp_g[i] for i in J[j]],coe_g[J[j]],lmon_h[W[j]],[supp_h[i] for i in W[j]],coe_h[W[j]],dg[J[j]],dh[W[j]],mex_cliq[j],lmex_cliq[j],v[j],Usupp[j],lUsupp[j],block_g0[j],block_g[j],block_h[j],lblock_g0[j],lblock_g[j],lblock_h[j],lt_block_g0[j],lt_block_g[j],lt_block_h[j],sk[j],sk_g[j],sk_h[j],s2k[j],s2k_g[j],s2k_h[j],use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)
       
    end
    
    println("  Largest constant trace: a=", maximum(ak))
    return omega_cliq,a0_cliq,a_ind1_cliq,a_ind2_cliq,a_val_cliq,a_len_cliq,a_mex_ind,a_mex_val,r_cliq,u_cliq,s_cliq,zeta_cliq,d_cliq,ak[1]


end


function get_mex_mix_CGAL(n::Int64,lUsupp::Vector{UInt64},Usupp::Vector{Vector{Vector{UInt64}}},p::Int64,I::Vector{Vector{UInt16}},lI::Vector{Int64})
    
    Uv=Vector{Vector{UInt64}}(undef,sum(lUsupp))
    
    t_U=0
    for j in 1:p
        Uv[1+t_U:lUsupp[j]+t_U]=Usupp[j]
        t_U+=lUsupp[j]
    end
    
    Uv=unique(Uv)
    lUv=length(Uv)
    
    mex=[Vector{UInt64}[] for j in 1:lUv]
    vec=zeros(UInt64,n)
    w=0
    
    for r in 1:lUv
        vec=Uv[r]
        for j in 1:p
            if issubset(unique(vec),I[j])
                w=ncbfind(Usupp[j],lUsupp[j],vec)
                if w>0
                    append!(mex[r],[[j,w]])
                end
            end
        end
    end
    
    
    l_mex=[length(mex[j]) for j in 1:lUv]
   
    
    
    Indmex=findall(j->l_mex[j]>1,1:lUv)
    mex=mex[Indmex]
    lt_mex=length(Indmex)
    l_mex=l_mex[Indmex]
    
    
    
    mex_cliq=[UInt64[] for j in 1:p]
    Index=zeros(UInt64,2)
    
    for j in 1:lt_mex
        for i in 1:l_mex[j]
            Index=mex[j][i]
            append!(mex_cliq[Index[1]],[Index[2]])
        end
    end

    lmex_cliq=[length(mex_cliq[j]) for j in 1:p]
    return mex,lt_mex,l_mex,mex_cliq,lmex_cliq
end

