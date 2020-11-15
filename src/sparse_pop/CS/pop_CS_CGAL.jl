function model_CS_CGAL(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{SparseMatrixCSC{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{SparseMatrixCSC{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::SparseMatrixCSC{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64;use_eqcons_to_get_constant_trace::Bool=true)

    
    I,p,lI=clique_decomp(n,m+l,[dg;dh],[[supp_f];supp_g;supp_h],order=k,alg="MD",minimize=true)
    
    #J,lJ,~=get_indcons(m,supp_g,I,p,lI,assign="all")
    
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
    
    #=println(" I=",I)
    println(" J=",J)
    println(" W=",W)=#
    
    mex,lt_mex,l_mex,mex_cliq,lmex_cliq,v,s2k,sort_v,re_ind=get_mex_CS_CGAL(k,n,lI,p,I)
    
    Indf,lIndf=decomp_obj(supp_f,lmon_f,I,p,n)
    
     omega_cliq,a0_cliq,a_ind1_cliq,a_ind2_cliq,a_val_cliq,a_len_cliq,a_mex_ind,a_mex_val,r_cliq,u_cliq,s_cliq,zeta_cliq,d_cliq,ak=run_model_cliq_CS_CGAL(p,lI,lJ,lW,I,J,W,lIndf,Indf,supp_f,coe_f,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,dg,dh,v,s2k,sort_v,re_ind,mex_cliq,lmex_cliq,k,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)

    
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
  
    @fastmath @inbounds @simd for j in 1:p
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
    
    
    
    
    @fastmath @inbounds @simd for j in 1:lt_mex

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
    
    
    a_val,a0,norm_a0,opnorm_a=rescale_dense(a_ind1,a_ind2,a_val,a_len,a0,zeta)
    
   
    
    
   
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


function update_CS_CGAL!(y::Vector{Float64},z::Vector{Float64},val::Vector{Float64},t::Int64,a_block::Vector{Vector{SparseMatrixCSC{Float64}}},a0_block::Vector{Vector{Adjoint{Float64,Array{Float64,1}}}},s_cliq::Vector{Vector{UInt64}},omega_cliq::Vector{Int64},norm_a0::Float64,opnorm_a::Float64,ak::Float64,p::Int64,zeta::UInt64)
    @fastmath @inbounds beta=sqrt(t+1)
    @fastmath @inbounds eta=2/(t+1)
    
    @fastmath @inbounds q=y+beta*z
                
    @fastmath @inbounds axpy!(-eta,z,z)#z=(1-eta)*z
    @fastmath @inbounds z[zeta]-=eta/opnorm_a/ak
    @fastmath @inbounds axpy!(-eta,val,val)#val=(1-eta)*val
    @fastmath @inbounds dualval=-y[end]/opnorm_a/ak
    
    @fastmath @inbounds mateigvec=zeros(Float64,1)
    @fastmath @inbounds eigval=Float64(0.0)
    @fastmath @inbounds ind=0

    
    @fastmath @inbounds @simd for j in 1:p
        eigval,mateigvec,ind=SmallEigBlocks_dense(a0_block[j],a_block[j],q,s_cliq[j],omega_cliq[j])
        dualval+=eigval
        axpy!(eta,a_block[j][ind]*mateigvec,z)#z+=eta*a_block[j][ind]*eigmat
        axpy!(eta,norm_a0*a0_block[j][ind]*mateigvec*ones(Float64,1),val)#val+=eta*(a0_block[j][ind]*eigmat)[1]*norm_a0
    end
                
    @fastmath @inbounds dualval*=norm_a0
                
    @fastmath @inbounds feas=norm(z)
    @fastmath @inbounds axpy!(minimum([1;beta*eta^2/feas^2]), z, y)#y+=minimum([1;beta*eta^2/feas^2])*z
                
    @fastmath @inbounds feas*=opnorm_a*ak
    
    @fastmath @inbounds gap=abs(val[1]-dualval)/(1+maximum([abs(val[1]);abs(dualval)]))           
     

    return feas,gap
end
         
                         
function solve_CS_CGAL(omega_cliq::Vector{Int64},a0_block::Vector{Vector{Adjoint{Float64,Array{Float64,1}}}},a_block::Vector{Vector{SparseMatrixCSC{Float64}}},s_cliq::Vector{Vector{UInt64}},zeta::UInt64,norm_a0::Float64,opnorm_a::Float64,ak::Float64,p::Int64;maxit::Int64=Int64(1e6),tol::Float64=1e-3,check_tol_each_iter::Bool=true)
  
    y=zeros(Float64,zeta)
    z=zeros(Float64,zeta)
    z[end]=-1/opnorm_a/ak
                
    val=zeros(Float64,1)
    feas=Float64(0)
    gap=Float64(0)
    time=Float64(0)

    if check_tol_each_iter
        i=1
        @fastmath @inbounds for t in 1:maxit
            feas,gap=update_CS_CGAL!(y,z,val,t,a_block,a0_block,s_cliq,omega_cliq,norm_a0,opnorm_a,ak,p,zeta)

            if  gap <tol && feas<tol
                println("iter=",t,"   val=",val[1],"   gap=",gap,"   feas=",feas)
                println("tol satisfies!!!")
                break
            end
            if t==i || t==maxit
                println("iter=",t,"   val=",val[1],"   gap=",gap,"   feas=",feas)
                i*=2
            end

        end
                
    else
        i=1
        while 2^i < maxit
            time+= @elapsed begin
            @fastmath @inbounds for t in 2^(i-1):2^i-1
                feas,gap=update_CS_CGAL!(y,z,val,t,a_block,a0_block,s_cliq,omega_cliq,norm_a0,opnorm_a,ak,p,zeta)
            end
                                    end
            println("iter=",2^i-1,"   val=",val[1],"   gap=",gap,"   feas=",feas,"   time=",time)
            if feas<tol && gap<tol
                println("tol satisfies!!!")
                break
            end
            i+=1
        end
    end
 
                
    println()
    println("####################################")
    println("opt_val = ",val[1])
    println("####################################")
    println("Solving time:")
    return val[1]
end
    


function POP_CS_CGAL(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{SparseMatrixCSC{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{SparseMatrixCSC{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::SparseMatrixCSC{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64;maxit::Int64=Int64(1e6),tol::Float64=1e-3,use_eqcons_to_get_constant_trace::Bool=true,check_tol_each_iter::Bool=true)

    @time begin
    
    
    @time omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p=model_CS_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)

      
    @time opt_val=solve_CS_CGAL(omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p,maxit=maxit,tol=tol,check_tol_each_iter=check_tol_each_iter)
    
    println("Total time:")
                    end
    return opt_val
end
                
                
           