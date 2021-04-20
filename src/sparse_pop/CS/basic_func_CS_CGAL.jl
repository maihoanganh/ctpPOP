function run_model_cliq_CS_CGAL(p::Int64,lI::Vector{Int64},lJ::Vector{Int64},lW::Vector{Int64},I::Vector{Vector{UInt16}},J::Vector{Vector{UInt64}},W::Vector{Vector{UInt64}},lIndf::Vector{Int64},Indf::Vector{Vector{UInt64}},supp_f::Vector{Vector{UInt64}},coe_f::Vector{Float64},lmon_g::Vector{UInt64},supp_g::Vector{Vector{Vector{UInt64}}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Vector{Vector{UInt64}}},coe_h::Vector{Vector{Float64}},dg::Vector{Int64},dh::Vector{Int64},v::Vector{Vector{Vector{UInt64}}},s2k::Vector{Int64},sort_v::Vector{Vector{Vector{UInt64}}},lsort_v::Vector{UInt64},mex_cliq::Vector{Vector{UInt64}},lmex_cliq::Vector{Int64},k::Int64;use_eqcons_to_get_constant_trace::Bool=true)
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
    
    
    @fastmath @inbounds @simd for j in 1:p
omega_cliq[j],a0_cliq[j],a_ind1_cliq[j],a_ind2_cliq[j],a_val_cliq[j],a_len_cliq[j],a_mex_ind[j],a_mex_val[j],r_cliq[j],u_cliq[j],s_cliq[j],zeta_cliq[j],d_cliq[j],ak[j]=model_Cliq_CS_CGAL(lI[j],lJ[j],lW[j],lIndf[j],supp_f[Indf[j]],coe_f[Indf[j]],lmon_g[J[j]],[supp_g[i] for i in J[j]],coe_g[J[j]],lmon_h[W[j]],[supp_h[i] for i in W[j]],coe_h[W[j]],dg[J[j]],dh[W[j]],v[j],s2k[j],sort_v[j],lsort_v[j],mex_cliq[j],lmex_cliq[j],k,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)
        
    end
    println("  Largest constant trace: a=", maximum(ak))
    #println("==================")
    return omega_cliq,a0_cliq,a_ind1_cliq,a_ind2_cliq,a_val_cliq,a_len_cliq,a_mex_ind,a_mex_val,r_cliq,u_cliq,s_cliq,zeta_cliq,d_cliq,ak[1]

end


function get_mex_CS_CGAL(k::Int64,n::Int64,lI::Vector{Int64},p::Int64,I::Vector{Vector{UInt16}})
    
    
    v=Vector{Vector{Vector{UInt64}}}(undef,p)
    sort_v=Vector{Vector{Vector{UInt64}}}(undef,p)
    lsort_v=Vector{UInt64}(undef,p)
    s2k=[numlet(lI[j],2*k) for j in 1:p]
    
    Uv=Vector{Vector{UInt64}}(undef,sum(s2k))
    
    t_U=0
    @fastmath @inbounds @simd for j in 1:p
        v[j]=get_ncbasis(lI[j],2*k)
        v[j]=[I[j][v[j][i]] for i in 1:s2k[j]]
        sort_v[j]=_sym_canon.(v[j])
        sort_v[j]=unique(sort_v[j])
        sort_v[j]=sort(sort_v[j])
        lsort_v[j]=length(sort_v[j])

        Uv[1+t_U:s2k[j]+t_U]=v[j]
        t_U+=s2k[j]
    end
    
    Uv=_sym_canon.(Uv)
    Uv=unique(Uv)
    lUv=length(Uv)
    
    mex=[Vector{UInt64}[] for j in 1:lUv]
    vec=zeros(UInt64,n)
    @fastmath @inbounds for r in 1:lUv
        vec=Uv[r]
        @fastmath @inbounds for j in 1:p
            if issubset(unique(vec),I[j])
                append!(mex[r],[[j,ncbfind(sort_v[j],lsort_v[j],vec)]])
            end
        end
    end
   
    
    l_mex=[@fastmath @inbounds length(mex[j]) for j in 1:lUv]
   
    
    
    Indmex=findall(j->l_mex[j]>1,1:lUv)
    mex=mex[Indmex]
    lt_mex=length(Indmex)
    l_mex=l_mex[Indmex]
    
    
    
    mex_cliq=[@fastmath @inbounds UInt64[] for j in 1:p]
    Index=zeros(UInt64,2)
    
    @fastmath @inbounds @simd for j in 1:lt_mex
        @fastmath @inbounds @simd for i in 1:l_mex[j]
            Index=mex[j][i]
            push!(mex_cliq[Index[1]],Index[2])
        end
    end

    lmex_cliq=[@fastmath @inbounds length(mex_cliq[j]) for j in 1:p]
    
    return mex,lt_mex,l_mex,mex_cliq,lmex_cliq,v,s2k,sort_v,lsort_v
end