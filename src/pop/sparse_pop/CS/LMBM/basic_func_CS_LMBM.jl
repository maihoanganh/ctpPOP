function run_model_cliq_CS_LMBM(p,lI,lJ,lW,I,J,W,lIndf,Indf,supp_f,coe_f,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,dg,dh,v,s2k,sort_v,re_ind,mex_cliq,lmex_cliq,k;UseEq=true)
    omega_cliq=Vector{Int64}(undef,p)
    a0_cliq=Vector{Vector{Float64}}(undef,p)
    a_ind1_cliq=Vector{Vector{UInt64}}(undef,p)
    a_ind2_cliq=Vector{Vector{UInt64}}(undef,p)
    a_val_cliq=Vector{Vector{Float64}}(undef,p)
    a_len_cliq=Vector{UInt64}(undef,p)
    
    a_mex_ind=Vector{Vector{UInt64}}(undef,p)
    a_mex_val=Vector{Vector{Float64}}(undef,p)
    r_cliq=Vector{Vector{UInt64}}(undef,p)
    u_cliq=Vector{Vector{UInt64}}(undef,p)
    s_cliq=Vector{Vector{UInt64}}(undef,p)
    
    zeta_cliq=Vector{UInt64}(undef,p)
    d_cliq=Vector{UInt64}(undef,p)
    
    
    j=1
    
    #println("Clique ",j,"th: ==================")
       omega_cliq[j],a0_cliq[j],a_ind1_cliq[j],a_ind2_cliq[j],a_val_cliq[j],a_len_cliq[j],a_mex_ind[j],a_mex_val[j],r_cliq[j],u_cliq[j],s_cliq[j],zeta_cliq[j],d_cliq[j],ak=model_Cliq_CS_LMBM(lI[j],lJ[j],lW[j],lIndf[j],supp_f[I[j],Indf[j]],coe_f[Indf[j]],lmon_g[J[j]],[supp_g[i][I[j],:] for i in J[j]],coe_g[J[j]],lmon_h[W[j]],[supp_h[i][I[j],:] for i in W[j]],coe_h[W[j]],dg[J[j]],dh[W[j]],v[j],s2k[j],sort_v[j],re_ind[j],mex_cliq[j],lmex_cliq[j],k,get_constant_trace=true,UseEq=UseEq)
    j+=1
    
    while j<= p
        #println("Clique ",j,"th: ==================")
           omega_cliq[j],a0_cliq[j],a_ind1_cliq[j],a_ind2_cliq[j],a_val_cliq[j],a_len_cliq[j],a_mex_ind[j],a_mex_val[j],r_cliq[j],u_cliq[j],s_cliq[j],zeta_cliq[j],d_cliq[j]=model_Cliq_CS_LMBM(lI[j],lJ[j],lW[j],lIndf[j],supp_f[I[j],Indf[j]],coe_f[Indf[j]],lmon_g[J[j]],[supp_g[i][I[j],:] for i in J[j]],coe_g[J[j]],lmon_h[W[j]],[supp_h[i][I[j],:] for i in W[j]],coe_h[W[j]],dg[J[j]],dh[W[j]],v[j],s2k[j],sort_v[j],re_ind[j],mex_cliq[j],lmex_cliq[j],k,UseEq=UseEq)
        j+=1
    end
    
    #println("==================")
    return omega_cliq,a0_cliq,a_ind1_cliq,a_ind2_cliq,a_val_cliq,a_len_cliq,a_mex_ind,a_mex_val,r_cliq,u_cliq,s_cliq,zeta_cliq,d_cliq,ak

end


function get_mex_CS_LMBM(k,n,lI,p,I)
    
    
    v=Vector{Matrix{UInt64}}(undef,p)
    sort_v=Vector{Matrix{UInt64}}(undef,p)
    re_ind=Vector{Vector{UInt64}}(undef,p)
    s2k=[binomial(2*k+lI[j],lI[j]) for j in 1:p]
    
    Uv=spzeros(UInt64,n,sum(s2k))
    
    t_U=0
    for j in 1:p
        v[j]=get_basis(lI[j],2*k)
        sort_v[j]=sortslices(v[j],dims=2)
        re_ind[j]=Vector{UInt64}(undef,s2k[j])
        @simd for i in 1:s2k[j]
            @inbounds re_ind[j][bfind(sort_v[j],s2k[j],v[j][:,i],lI[j])]=i
        end
        Uv[I[j],1+t_U:s2k[j]+t_U]=v[j]
        t_U+=s2k[j]
    end
    
    
    Uv=unique(Uv,dims=2)
    lUv=size(Uv,2)
    
    
    
    mex=[Vector{UInt64}[] for j in 1:lUv]
    vec=spzeros(UInt64,n)
    for r in 1:lUv
        vec=Uv[:,r]
        for j in 1:p
            if issubset(findall(y->y>0,vec),I[j])
                append!(mex[r],[[j,re_ind[j][bfind(sort_v[j],s2k[j],vec[I[j]],lI[j])]]])
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
    
    return mex,lt_mex,l_mex,mex_cliq,lmex_cliq,v,s2k,sort_v,re_ind
end