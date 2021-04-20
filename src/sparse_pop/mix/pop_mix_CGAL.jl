function model_mix_CGAL(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Vector{Vector{UInt64}}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Vector{Vector{UInt64}}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::Vector{Vector{UInt64}},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64,t::Int64;use_eqcons_to_get_constant_trace::Bool=true)
    
    
    
    function matsupp(supp,lmon)
        mat=spzeros(UInt64,n,lmon)
        for j in 1:lmon
            mat[unique(supp[j]),j].=1
        end
        return mat
    end
    
    I,p,lI=clique_decomp(n,m+l,[dg;dh],[[matsupp(supp_f,lmon_f)];[matsupp(supp_g[i],lmon_g[i]) for i in 1:m];[matsupp(supp_h[i],lmon_h[i]) for i in 1:l]],order="max",alg="MD",minimize=true)
    
    indJ_out=Vector{UInt64}([])
    indJ_in=Vector{UInt64}([])
    for j in 1:m
        if coe_g[j][end]<0
            push!(indJ_out,j)
        else
            push!(indJ_in,j)
        end
    end
    
    
    posJ_in,lposJ_in,~=get_indcons(length(indJ_in),supp_g[indJ_in],I,p,lI,assign="all")
    posJ_out,lposJ_out,~=get_indcons(length(indJ_out),supp_g[indJ_out],I,p,lI,assign="min")
    
    J=[union(indJ_in[posJ_in[j]],indJ_out[posJ_out[j]]) for j in 1:p]
    lJ=[lposJ_in[j]+lposJ_out[j] for j in 1:p]
    
    W,lW,~=get_indcons(l,supp_h,I,p,lI,assign="min")
    
    println("  Number of cliques: p=", p)
    println("  Largest clique size: u=", maximum(lI))
    
    lSupp=lmon_f+sum(lmon_g)+sum(lmon_h)
    Supp=Vector{Vector{UInt64}}(undef,lSupp)
    Supp[1:lmon_f]=supp_f
    t_blo=lmon_f
    for j in 1:m
        Supp[1+t_blo:lmon_g[j]+t_blo]=supp_g[j]
        t_blo+=lmon_g[j]
    end
    for j in 1:l
        Supp[1+t_blo:lmon_h[j]+t_blo]=supp_h[j]
        t_blo+=lmon_h[j]
    end
    Supp=_sym_canon.(Supp)
    Supp=unique(Supp)
    lSupp=length(Supp)
    
    
    #Need to remove
    IndA=[UInt64[] for i in 1:p]
    for j in 1:lSupp
        ind=unique(Supp[j])
        for i in 1:p
            if issubset(ind,I[i])
                push!(IndA[i],j)
            end
        end
    end
    
    
    v,Usupp,lUsupp,block_g0,block_g,block_h,lblock_g0,lblock_g,lblock_h,lt_block_g0,lt_block_g,lt_block_h,sk,sk_g,sk_h,s2k,s2k_g,s2k_h=run_get_blocks_Cliq_mix_CGAL(n,p,t,k,lI,lJ,lW,I,J,W,Supp,IndA,lmon_g,lmon_h,supp_g,supp_h,coe_g,coe_h,dg,dh)
        
   
    
    
    
    mex,lt_mex,l_mex,mex_cliq,lmex_cliq=get_mex_mix_CGAL(n,lUsupp,Usupp,p,I,lI)
    
    Indf,lIndf=decomp_obj(supp_f,lmon_f,I,p,n)
    
    omega_cliq,a0_cliq,a_ind1_cliq,a_ind2_cliq,a_val_cliq,a_len_cliq,a_mex_ind,a_mex_val,r_cliq,u_cliq,s_cliq,zeta_cliq,d_cliq,ak=run_model_cliq_mix_CGAL(p,lI,lJ,lW,I,J,W,lIndf,Indf,supp_f,coe_f,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,dg,dh,mex_cliq,lmex_cliq,v,Usupp,lUsupp,block_g0,block_g,block_h,lblock_g0,lblock_g,lblock_h,lt_block_g0,lt_block_g,lt_block_h,sk,sk_g,sk_h,s2k,s2k_g,s2k_h,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)
    
    
    d=sum(d_cliq)
    zeta=sum(zeta_cliq)+sum(l_mex)-lt_mex+1
    
    
    println("  Number of blocks: omega=", sum(omega_cliq))
    println("  Number of equality consraints: zeta=", zeta)
    println("  Size of the largest block: s^max=", maximum([s_cliq[j][1] for j in 1:p]))
    
    a_ind1=Vector{UInt64}([])
    a_ind2=Vector{UInt64}([])
    a_val=Vector{Float64}([])
    a_len=UInt64(0)  
    
    a0=zeros(Float64,d)
    r=Vector{Int64}(undef,p)
    
    t1=0
    t2=0
    for j in 1:p
        r[j]=t1
        
        @fastmath @inbounds @simd for i in 1:a_len_cliq[j]
            push!(a_ind1,a_ind1_cliq[j][i]+t1)
            push!(a_ind2,a_ind2_cliq[j][i]+t2)
            push!(a_val,a_val_cliq[j][i])
            a_len+=1
        end
        
        a0[t1+1:t1+d_cliq[j]]=a0_cliq[j]
        t1+=d_cliq[j]
        t2+=zeta_cliq[j]
    end
    
    Ind=zeros(UInt64,2)
    Ind1=zeros(UInt64,2)
    w=UInt64(0)
    I=UInt64(0)
    V=Float64(0)

    for j in 1:lt_mex
        Ind=mex[j][1]
        w=findfirst(y->y==Ind[2],mex_cliq[Ind[1]])
        I=a_mex_ind[Ind[1]][w]
        V=a_mex_val[Ind[1]][w]
        
        @fastmath @inbounds @simd for i in 2:l_mex[j]
            
            push!(a_ind1,I+r[Ind[1]])
            push!(a_ind2,1+t2)
            push!(a_val,-V)
            a_len+=1
            
            Ind1=mex[j][i]
            w1=findfirst(y->y==Ind1[2],mex_cliq[Ind1[1]])
            
            push!(a_ind1,r[Ind1[1]]+a_mex_ind[Ind1[1]][w1])
            push!(a_ind2,1+t2)
            push!(a_val,a_mex_val[Ind1[1]][w1])
            a_len+=1
            
            t2+=1
        end
        
    end
   
    
    push!(a_ind1,1)
    push!(a_ind2,zeta)
    push!(a_val,ak)
    a_len+=1
    
    
    _,_,a_val,_,a0,norm_a0,opnorm_a,_,_=rescale_dense(a_ind1,a_ind2,a_val,a_len,a0,zeta)
    
    
    Ind=Vector{Vector{Vector{UInt64}}}(undef,p)
    a_block=Vector{Vector{SparseMatrixCSC{Float64}}}(undef,p)
    a0_block=Vector{Vector{Adjoint{Float64,Array{Float64,1}}}}(undef,p)
    Index=Vector{UInt64}(1:a_len)
    pos=zeros(UInt64,1)
 
    @fastmath @inbounds @simd for i in 1:p
        
        Ind[i]=[@fastmath @inbounds Vector{UInt64}([]) for j in 1:omega_cliq[i]]
        @fastmath @inbounds @simd for j in 1:omega_cliq[i]
            pos=findall(t->(a_ind1[t]>=(r[i]+r_cliq[i][j]+1)) && (a_ind1[t]<=(r[i]+r_cliq[i][j]+u_cliq[i][j])),Index)
            Ind[i][j]=Index[pos]
            deleteat!(Index, pos)
        end

        a_block[i]=Vector{SparseMatrixCSC{Float64}}(undef,omega_cliq[i])
       
        a0_block[i]=Vector{Adjoint{Float64,Array{Float64,1}}}(undef,omega_cliq[i])
        
        @fastmath @inbounds @simd for ind in 1:omega_cliq[i]
            a_block[i][ind]=sparse(a_ind2[Ind[i][ind]],a_ind1[Ind[i][ind]].-(r[i]+r_cliq[i][ind]),a_val[Ind[i][ind]],zeta,u_cliq[i][ind])
            
            a0_block[i][ind]=a0[r[i]+r_cliq[i][ind]+1:r[i]+r_cliq[i][ind]+u_cliq[i][ind]]'
        end
    end
    
    
    
    println("Modeling time:")
 
    return omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p
end





function POP_mix_CGAL(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Vector{Vector{UInt64}}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Vector{Vector{UInt64}}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::Vector{Vector{UInt64}},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64,t::Int64;maxit::Int64=Int64(1e6),tol::Float64=1e-3,use_eqcons_to_get_constant_trace::Bool=true,check_tol_each_iter::Bool=true)

    @time begin
    
    
    @time omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p=model_mix_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,t,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)

      
    @time opt_val=solve_CS_CGAL(omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p,maxit=maxit,tol=tol,check_tol_each_iter=check_tol_each_iter)
    
    println("Total time:")
                    end
    return opt_val
end
                
   