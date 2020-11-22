function get_constant_trace_dense(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Matrix{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Matrix{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::Matrix{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64;use_eqcons_to_get_constant_trace::Bool=true)
    
    
   
            
    v=get_basis(n,2*k)
    
    sk=binomial(k+n,n)
    sk_g=Vector{UInt64}(undef,m)
    s2k_g=Vector{UInt64}(undef,m)
    sk_h=Vector{UInt64}(undef,l)
    s2k_h=Vector{UInt64}(undef,l)
    
   
        
    
    ceil_g=Int64(0)
    @fastmath @inbounds @simd for i in 1:m
        ceil_g=ceil(Int64,dg[i]/2)
        sk_g[i]=binomial(k-ceil_g+n,n)
        s2k_g[i]=binomial(2*(k-ceil_g)+n,n)
    end
                    
    @fastmath @inbounds @simd for i in 1:l
        ceil_g=ceil(Int64,dh[i]/2)
        sk_h[i]=binomial(k-ceil_g+n,n)
        s2k_h[i]=binomial(2*(k-ceil_g)+n,n)
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
    V=Matrix{UInt64}(undef,n,lV)
    V[:,1:sk]=2*v[:,1:sk]
    t=sk+1
    @fastmath @inbounds @simd for i in 1:m
        @fastmath @inbounds @simd for r in 1:lmon_g[i]
            @fastmath @inbounds @simd for j in 1:sk_g[i]
                V[:,t]=2*v[:,j]+supp_g[i][:,r]
                t+=1
            end
        end
    end
    if use_eqcons_to_get_constant_trace               
        @fastmath @inbounds @simd for i in 1:l
            @fastmath @inbounds @simd for r in 1:lmon_h[i]
                @fastmath @inbounds @simd for j in 1:s2k_h[i]
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
    

    
    @fastmath @inbounds @simd for j in 1:sk
        add_to_expression!(cons[bfind(V,lV,2*v[:,j],n)],p0[j])
    end


    @fastmath @inbounds @simd for i in 1:m
        p[i]=@variable(model, [1:sk_g[i]], lower_bound=1)
        @fastmath @inbounds @simd for j in 1:sk_g[i]
            @fastmath @inbounds @simd for r in 1:lmon_g[i]
                add_to_expression!(cons[bfind(V,lV,2*v[:,j]+supp_g[i][:,r],n)],p[i][j]*coe_g[i][r])
            end
        end
    end
    
    if use_eqcons_to_get_constant_trace
        q=Vector{Vector{VariableRef}}(undef, l)
        @fastmath @inbounds @simd for i in 1:l
            q[i]=@variable(model, [1:s2k_h[i]])
            @fastmath @inbounds @simd for j in 1:s2k_h[i]
                @fastmath @inbounds @simd for r in 1:lmon_h[i]
                    add_to_expression!(cons[bfind(V,lV,v[:,j]+supp_h[i][:,r],n)],q[i][j]*coe_h[i][r])
                end
            end
        end
    end
    @constraint(model, cons[2:end].==0)
    @variable(model, lambda)
    @constraint(model, cons[1]==lambda)
    @objective(model, Min, lambda)
    optimize!(model)
                                    
    println("  Computing constant trace status: ", termination_status(model))
                                    
    ak=value(lambda)
    println("  Constant trace: ak = ",ak)
    
    P=sqrt.(value.(p0))
    Pg=[sqrt.(value.(p[i])) for i in 1:m]

    
    return v,sk,sk_g,s2k_g,sk_h,s2k_h,ak,P,Pg
end


function model_POP_dense(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Matrix{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Matrix{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::Matrix{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64;use_eqcons_to_get_constant_trace::Bool=true)
    
    v,sk,sk_g,s2k_g,sk_h,s2k_h,ak,P,Pg=get_constant_trace_dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)
    
    s2k=size(v,2)
    
    sort_v=sortslices(v,dims=2)
    re_ind=Vector{UInt64}(undef,s2k)
    @fastmath @inbounds @simd for j in 1:s2k
        re_ind[bfind(sort_v,s2k,v[:,j],n)]=j
    end
    
    Order(alpha::Vector{UInt64})=re_ind[bfind(sort_v,s2k,alpha,n)]        
    
    println("  Number of blocks: omega=",m+1)
    println("  Size of the largest block: s^max=",sk)
    
    u=Int64(0.5*sk*(sk+1))
    u_g=[@fastmath @inbounds Int64(0.5*sk_g[i]*(sk_g[i]+1)) for i in 1:m]
    d=u+sum(u_g)
    zeta=d-s2k+sum(s2k_h)+1   
    println("  Number of equality trace constraints: zeta=",zeta) 
            

    IndM=[@fastmath @inbounds Vector{Vector{UInt64}}([]) for j in 1:s2k]
    invIndeM=[Vector{UInt64}([]) for i in 1:sk]
    r=UInt64(0)
    t_iter=UInt64(1)
    

    a_ind1=Vector{UInt64}([])
    a_ind2=Vector{UInt64}([])
    a_val=Vector{Float64}([])
    a_len=UInt64(0)   
    
    
    @fastmath @inbounds for i in 1:sk, j in 1:i
        r=Order(v[:,i]+v[:,j])
        push!(IndM[r],[i;j])
        push!(invIndeM[i],t_iter)
        t_iter+=1
    end
    
    l_IndM=[length(IndM[r]) for r in 1:s2k]
 
         
    t_a=UInt64(1)
    I1=UInt64(0)
    I2=UInt64(0)
    
    @fastmath @inbounds @simd for r in 1:s2k
        if l_IndM[r]>1
            @fastmath @inbounds @simd for i in 2:l_IndM[r]
                I1,I2=IndM[r][1]
                push!(a_ind1,invIndeM[I1][I2])
                push!(a_ind2,t_a)
                push!(a_val,ak/P[I1]/P[I2]*((0.5*sqrt(2)-1)*(I2<I1)+1))
                a_len+=1
                
                I1,I2=IndM[r][i]
                push!(a_ind1,invIndeM[I1][I2])
                push!(a_ind2,t_a)
                push!(a_val,-ak/P[I1]/P[I2]*((0.5*sqrt(2)-1)*(I2<I1)+1))
                a_len+=1
                    
                    
                t_a+=1
             end
             IndM[r]=Vector{Int64}[IndM[r][1]]
         end
    end
   
    
    
    
    
    IndMg=[@fastmath @inbounds [Vector{Vector{Int64}}([]) for j in 1:s2k_g[i]] for i in 1:m]
    invIndeMg=[@fastmath @inbounds [Vector{UInt64}([]) for j in 1:sk_g[i]] for i in 1:m]
    t_Blo=u
    
    
 
    @fastmath @inbounds @simd for j in 1:m
        r=UInt64(0)
        t_iter=UInt64(1)
        @fastmath @inbounds for p in 1:sk_g[j], q in 1:p
            r=Order(v[:,p]+v[:,q])
            push!(IndMg[j][r],[p;q])
            push!(invIndeMg[j][p],t_iter)
            t_iter+=1
        end
                
        
        l_IndM=[@fastmath @inbounds length(IndMg[j][q]) for q in 1:s2k_g[j]]
       
        @fastmath @inbounds @simd for r in 1:s2k_g[j]
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
            
                    I1,I2=IndMg[j][r][1]
                    push!(a_ind1,invIndeMg[j][I1][I2]+t_Blo)
                    push!(a_ind2,t_a)
                    push!(a_val,-ak/Pg[j][I1]/Pg[j][I2]*((0.5*sqrt(2)-1)*(I2<I1)+1))
                    a_len+=1
            
              
            @fastmath @inbounds @simd  for p in 1:lmon_g[j]  
                    I1,I2=IndM[Order(supp_g[j][:,p]+v[:,r])][end]
                    push!(a_ind1,invIndeM[I1][I2])
                    push!(a_ind2,t_a)
                    push!(a_val,coe_g[j][p]*ak/P[I1]/P[I2]*((0.5*sqrt(2)-1)*(I2<I1)+1))
                    a_len+=1
                 
                   end
                   t_a+=1
            
         end
        
         
        
         t_Blo+=u_g[j]
    end
   
    
               
            
                     
            
    
    
  
            
            
  
    @fastmath @inbounds @simd for j in 1:l
        @fastmath @inbounds @simd  for r in 1:s2k_h[j]
            @fastmath @inbounds @simd  for p in 1:lmon_h[j]
                      I1,I2=IndM[Order(supp_h[j][:,p]+v[:,r])][1]
                        push!(a_ind1,invIndeM[I1][I2])
                        push!(a_ind2,t_a)
                        push!(a_val,coe_h[j][p]*ak/P[I1]/P[I2]*((0.5*sqrt(2)-1)*(I2<I1)+1))
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
    
    
    @fastmath @inbounds @simd  for p in 1:lmon_f
        I1,I2=IndM[Order(supp_f[:,p])][1]
        a0[invIndeM[I1][I2]]=coe_f[p]*ak/P[I1]/P[I2]*((0.5*sqrt(2)-1)*(I2<I1)+1)
   end
    
    a_val,a0,norm_a0,opnorm_a=rescale_dense(a_ind1,a_ind2,a_val,a_len,a0,zeta)
    
    omega=m+1
    
    
    r_blo=Vector{UInt64}(undef,omega)
    u_vec=Vector{UInt64}(undef,omega)
    s=Vector{UInt64}(undef,omega)
    
    t=1
    u_vec[t]=u
    s[t]=sk
    r_blo[t]=0
    t+=1
    
    t_Blo=u
    @fastmath @inbounds @simd  for j in 1:m
        r_blo[t]=t_Blo
        u_vec[t]=u_g[j]
        s[t]=sk_g[j]
        t+=1
        t_Blo+=u_g[j]
    end

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
     

function rescale_dense(a_ind1::Vector{UInt64},a_ind2::Vector{UInt64},a_val::Vector{Float64},a_len::UInt64,a0::Vector{Float64},zeta::UInt64)
    

    norm_a=zeros(Float64,zeta)
    
    @fastmath @inbounds @simd  for t in 1:a_len
        norm_a[a_ind2[t]]+=a_val[t]^2
    end
        
    norm_a=sqrt.(norm_a)
   
    
   
    @fastmath @inbounds @simd  for t in 1:a_len
        a_val[t]/=norm_a[a_ind2[t]]
    end
    
    
    opnorm_a=svds(sparse(a_ind1,a_ind2,a_val), nsv = 1)[1].S[1]
    
    a_val=a_val./opnorm_a

    
    norm_a0=norm(a0)
    a0=a0./norm_a0
    
    
    return a_val,a0,norm_a0,opnorm_a
end



function SmallEig_dense(mat::Matrix{Float64},s::UInt64)
    try
       @fastmath @inbounds E=eigs(mat,nev = 1,which=:SR,tol=1e-2) 
       return E[1][1],E[2][:,1]
    catch
       @fastmath @inbounds E=eigen(Symmetric(mat),1:1)
       return E.values[1],E.vectors[:,1]
    end
end

function getmat_dense(vec::Adjoint{Float64,Array{Float64,1}},sk::UInt64)
    B=zeros(Float64,sk,sk)
    r=1
    @fastmath @inbounds for i in 1:sk, j in 1:i
        B[i,j]=vec[r]/sqrt(1+(j<i))
        B[j,i]= copy(B[i,j])
        r+=1
    end
    return B
end

getmatvec_dense(vec::Array{Float64,1},s::UInt64)=[@fastmath @inbounds vec[i]*vec[j]*sqrt(1+(j<i)) for i in 1:s for j in 1:i]
                                        
function SmallEig_block_dense(vec::Adjoint{Float64,Array{Float64,1}},sk::UInt64)
    if sk==1
        return vec[1], ones(Float64,1)
    else
        return SmallEig_dense(getmat_dense(vec,sk),sk)
    end

end



function SmallEigBlocks_dense(a0_block::Vector{Adjoint{Float64,Array{Float64,1}}},a_block::Vector{SparseMatrixCSC{Float64}},y::Vector{Float64},s::Vector{UInt64},omega::Int64)
    
    
    eigval=Vector{Float64}(undef,omega)
    eigvec=Vector{Vector{Float64}}(undef,omega)


    @fastmath @inbounds @simd for p in 1:omega
        eigval[p],eigvec[p]=SmallEig_block_dense(a0_block[p]+y'*a_block[p],s[p])
    end
               
    ~,ind=findmin(eigval)

    return eigval[ind],getmatvec_dense(eigvec[ind],s[ind]),ind
end       



function update_dense_CGAL!(y::Vector{Float64},z::Vector{Float64},val::Vector{Float64},t::Int64,a_block::Vector{SparseMatrixCSC{Float64}},a0_block::Vector{Adjoint{Float64,Array{Float64,1}}},s::Vector{UInt64},omega::Int64,norm_a0::Float64,opnorm_a::Float64,ak::Float64)

    @fastmath @inbounds beta=sqrt(t+1)
    @fastmath @inbounds eta=2/(t+1)
    
    @fastmath @inbounds eigval,matveceig,ind=SmallEigBlocks_dense(a0_block,a_block,y+beta*z,s,omega)
    
    @fastmath @inbounds dualval=(eigval-y[end]/opnorm_a/ak)*norm_a0 
      
    @fastmath @inbounds axpy!(eta, a_block[ind]*matveceig-z, z) #z+=eta*(a_block[ind]*matveceig-z)            
    z[end]-=eta/opnorm_a/ak
                
    @fastmath @inbounds feas=norm(z)

    @fastmath @inbounds axpy!(minimum([1;beta*eta^2/feas^2]), z, y)#y+=minimum([1;beta*eta^2/feas^2])*z
                
    @fastmath @inbounds feas*=opnorm_a*ak
                
    @fastmath @inbounds axpy!(eta,norm_a0*(a0_block[ind]*matveceig).-val,val)#val+=eta*((a0_block[ind]*matveceig)[1]*norm_a0-val)
    
    @fastmath @inbounds gap=abs(val[1]-dualval)/(1+maximum([abs(val[1]);abs(dualval)]))

               
    return feas,gap
end

function solve_POP_dense_CGAL(omega::Int64,a0_block::Vector{Adjoint{Float64,Array{Float64,1}}},a_block::Vector{SparseMatrixCSC{Float64}},s::Vector{UInt64},zeta::UInt64,norm_a0::Float64,opnorm_a::Float64,ak::Float64;maxit::Int64=1e6,tol::Float64=1e-4,check_tol_each_iter::Bool=true)
    
    y=zeros(Float64,zeta)
    z=zeros(Float64,zeta)
    z[end]=-1/opnorm_a/ak
                
    val=zeros(Float64,1)
    feas=Float64(1)
    gap=Float64(1)
    time=Float64(0)  
                
    if check_tol_each_iter
        i=1
        @fastmath @inbounds for t in 1:maxit
            feas,gap=update_dense_CGAL!(y,z,val,t,a_block,a0_block,s,omega,norm_a0,opnorm_a,ak)

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
                feas,gap=update_dense_CGAL!(y,z,val,t,a_block,a0_block,s,omega,norm_a0,opnorm_a,ak)           
            end     

                        end
            println("iter=",2^i-1,"   val=",val[1],"   gap=",gap,"   feas=",feas,"   time=",time)
            if  gap <tol && feas<tol
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
    


function POP_dense_CGAL(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Matrix{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Matrix{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::Matrix{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64;maxit::Int64=1e6,tol::Float64=1e-4,use_eqcons_to_get_constant_trace::Bool=true,check_tol_each_iter::Bool=true)

    @time begin
    
    
    @time omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak=model_POP_dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)
    

      
    @time opt_val=solve_POP_dense_CGAL(omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak,maxit=maxit,tol=tol,check_tol_each_iter=check_tol_each_iter)
    
              
       println("Total time:")
                    end
    return opt_val
end
                    