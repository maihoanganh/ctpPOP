function get_constant_trace_DSMA(x::Vector{PolyVar{true}},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64;UseEq::Bool=true)
    
    
    n=length(x)
    m=length(g)
    l=length(h)
            
    v=get_basis(n,2*k)
    
    sk=binomial(k+n,n)
    sk_g=Vector{UInt64}(undef,m)
    s2k_g=Vector{UInt64}(undef,m)
    s2k_h=Vector{UInt64}(undef,l)
    
    lmon_g=Vector{UInt64}(undef,m)
    supp_g=Vector{Matrix{UInt64}}(undef,m)
    coe_g=Vector{Vector{Float64}}(undef,m)
    
    lmon_h=Vector{UInt64}(undef,l)
    supp_h=Vector{Matrix{UInt64}}(undef,l)
    coe_h=Vector{Vector{Float64}}(undef,l)
        
    
    ceil_g=Int64(0)
    @fastmath @inbounds @simd for i in 1:m
        ceil_g=ceil(Int64,maxdegree(g[i])/2)
        sk_g[i]=binomial(k-ceil_g+n,n)
        s2k_g[i]=binomial(2*(k-ceil_g)+n,n)
        lmon_g[i],supp_g[i],coe_g[i]=info(g[i],x,n)
    end
                    
    @fastmath @inbounds @simd for i in 1:l
        s2k_h[i]=binomial(2*(k-ceil(Int64,maxdegree(h[i])/2))+n,n)
        lmon_h[i],supp_h[i],coe_h[i]=info(h[i],x,n)
    end

    lV=sk
    
    if m!=0
        lV+=sum(sk_g[i]*lmon_g[i] for i in 1:m)
    end
    if UseEq
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
    if UseEq               
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
    
    if UseEq
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

    
    return n,m,l,v,sk,sk_g,s2k_g,s2k_h,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,ak,P,Pg
end


function model_POP_DSMA(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64;UseEq::Bool=true)
    
    n,m,l,v,sk,sk_g,s2k_g,s2k_h,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,ak,P,Pg=get_constant_trace_DSMA(x,g,h,k,UseEq=UseEq)
    
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
    
    
    lmon_f,supp_f,coe_f=info(f,x,n)
    @fastmath @inbounds @simd  for p in 1:lmon_f
        I1,I2=IndM[Order(supp_f[:,p])][1]
        a0[invIndeM[I1][I2]]=coe_f[p]*ak/P[I1]/P[I2]*((0.5*sqrt(2)-1)*(I2<I1)+1)
   end
    
    a_val,a0,norm_a0,opnorm_a=rescale_CGAL(a_ind1,a_ind2,a_val,a_len,a0,zeta)
    
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
    
    a_block=[@fastmath @inbounds sparse(a_ind2[Ind[i]],a_ind1[Ind[i]].-r_blo[i],a_val[Ind[i]],zeta,u_vec[i]) for i in 1:omega]
    
    #@time a_block=[SparseMatrixCSC{Float64}(a[r_blo[i]+1:r_blo[i]+u_vec[i],:]') for i in 1:omega]
    a0_block=[@fastmath @inbounds a0[r_blo[i]+1:r_blo[i]+u_vec[i]]' for i in 1:omega]
 
    
    println("Modeling time:")
 
    return omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak
end
     

function rescale_DSMA(a_ind1::Vector{UInt64},a_ind2::Vector{UInt64},a_val::Vector{Float64},a_len::UInt64,a0::Vector{Float64},zeta::UInt64)
    

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



function SmallEig_DSMA(mat::Matrix{Float64},s::UInt64;EigAlg::String="Arpack")
    if s==1
        return mat[1,1], ones(Float64,1)
    #=elseif norm(mat)==0
        return 0.0, [1;zeros(Float64,s-1)]=#
    elseif EigAlg=="Arpack"
       @fastmath @inbounds E=eigs(mat,nev = 1,which=:SR,tol=1e-2) 
       return E[1][1],E[2][:,1]
    elseif EigAlg=="ArnoldiMethod"      
       @fastmath @inbounds E=partialeigen(partialschur(mat, nev=1,tol=1e-2, which=SR())[1])
       return E[1][1],E[2][:,1]
    elseif EigAlg=="Normal"
       @fastmath @inbounds E=eigen(Symmetric(mat),1:1)
       return E.values[1],E.vectors[:,1]
    elseif EigAlg=="Mix"
       try 
           @fastmath @inbounds E=partialeigen(partialschur(mat, nev=1,tol=1e-2, which=SR())[1])
           return E[1][1],E[2][:,1]
       catch
           @fastmath @inbounds E=eigs(mat,nev = 1,which=:SR,tol=1e-2) 
           return E[1][1],E[2][:,1]
       end
    else
       println("No eigenvalue algorithm!!!")
    end  
    
end

                                        
function SmallEig_block_DSMA(vec::Adjoint{Float64,Array{Float64,1}},sk::UInt64;EigAlg::String="Arpack")
    
    B=zeros(Float64,sk,sk)
    r=1
    @fastmath @inbounds for i in 1:sk, j in 1:i
        B[i,j]=vec[r]/sqrt(1+(j<i))
        B[j,i]= copy(B[i,j])
        r+=1
    end
 
    return SmallEig_DSMA(B,sk,EigAlg=EigAlg)

end




function SmallEigBlocks_DSMA(a0_block::Vector{Adjoint{Float64,Array{Float64,1}}},a_block::Vector{SparseMatrixCSC{Float64,UInt64}},y::Vector{Float64},s::Vector{UInt64},omega::Int64;EigAlg::String="Arpack")
    
    
    eigval,eigvec=Float64(Inf),zeros(Float64,1)

    smalleigval=eigval
    smalleigvec=eigvec
    ind=0

 
    @fastmath @inbounds @simd for p in 1:omega
        eigval,eigvec=SmallEig_block_DSMA(a0_block[p]+y'*a_block[p],s[p],EigAlg=EigAlg)
        if smalleigval>eigval
            smalleigval=eigval
            smalleigvec=eigvec
            ind=p
        end
    end

    return smalleigval,smalleigvec,ind
end       



function update_DSMA(y_plus::Vector{Float64},y_minus::Vector{Float64},z_plus::Vector{Float64},z_minus::Vector{Float64},val::Float64,t::Float64,Gamma::Float64,a_block::Vector{SparseMatrixCSC{Float64,UInt64}},a0_block::Vector{Adjoint{Float64,Array{Float64,1}}},s::Vector{UInt64},omega::Int64,norm_a0::Float64,opnorm_a::Float64,ak::Float64,eps::Float64;EigAlg::String="Arpack")

    @fastmath @inbounds eta=t/(t+1)
    @fastmath @inbounds theta=(t+1)/(t+2)
    @fastmath @inbounds gamma=1/sqrt(t+1)
    @fastmath @inbounds Gamma=(t*Gamma+gamma)./(t+1)
    @fastmath @inbounds beta=(1-theta)*(1/Gamma)
    
    @fastmath @inbounds eigval,eigvec,ind=SmallEigBlocks_DSMA(a0_block,a_block,y_plus-y_minus,s,omega,EigAlg=EigAlg)
    
    @fastmath @inbounds dualval=(eigval-(y_plus[end]*(1+eps)+y_minus[end]*(-1+eps))/opnorm_a/ak)*norm_a0 
    

    p=[@fastmath @inbounds eigvec[i]*eigvec[j]*sqrt(1+(j<i)) for i in 1:s[ind] for j in 1:i]
                
    #obj=val+(y[end]+beta*z[end])+0.5*beta*feas^2-dot(a0[Ind[ind]],p)-dot(y+beta*z,a[Ind[ind],:]'*p)
                
    @fastmath @inbounds u=(1-eta)*(a_block[ind]*p)            
                
    @fastmath @inbounds z_plus=eta*z_plus+u            
    z_plus[end]-=(1-eta)*(1+eps)/opnorm_a/ak
                
    @fastmath @inbounds z_minus=eta*z_minus-u            
    z_minus[end]-=(1-eta)*(-1+eps)/opnorm_a/ak
                
    

    @fastmath @inbounds y_plus=theta*y_plus
    @fastmath @inbounds ind_plus=findall(i->i>0,z_plus)
    @fastmath @inbounds y_plus[ind_plus]+=beta*z_plus[ind_plus]
                
    @fastmath @inbounds y_minus=theta*y_minus
    @fastmath @inbounds ind_minus=findall(i->i>0,z_minus)
    @fastmath @inbounds y_minus[ind_minus]+=beta*z_minus[ind_minus]

    @fastmath @inbounds feas=opnorm_a*ak*norm([z_plus[ind_plus];z_minus[ind_minus]])
                
    @fastmath @inbounds val=eta*val+(1-eta)*(a0_block[ind]*p)[1]*norm_a0
    @fastmath @inbounds gap=abs(val-dualval)
               
             
    return y_plus,y_minus,z_plus,z_minus,val,feas,gap
end

function solve_POP_DSMA(omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak;EigAlg="Arpack",maxit=1e6,tol=1e-4,eps::Float64=1e-3)
  
    y_plus=zeros(Float64,zeta)
    z_plus=zeros(Float64,zeta)
    z_plus[end]=-(1+eps)/opnorm_a/ak
                
    y_minus=zeros(Float64,zeta)
    z_minus=zeros(Float64,zeta)
    z_minus[end]=-(-1+eps)/opnorm_a/ak
                
    val=Float64(0)
    feas=Float64(0)
                
    gap=Float64(0)
    Gamma=Float64(0)
                
    #stop=1.0
    i=1
    @fastmath @inbounds for t in 1:maxit
        y_plus,y_minus,z_plus,z_minus,val,feas,gap=update_DSMA(y_plus,y_minus,z_plus,z_minus,val,t,Gamma,a_block,a0_block,s,omega,norm_a0,opnorm_a,ak,eps,EigAlg=EigAlg)
                   
        if  gap <tol && feas<tol
            println("iter=",t,"   val=",val,"   gap=",gap,"   feas=",feas)
            println("tol satisfies!!!")
            break
        end
        if t==i || t==maxit
            println("iter=",t,"   val=",val,"   gap=",gap,"   feas=",feas)
            i*=2
        end
                    
    end
    
 
                
    println()
    println("####################################")
    println("opt_val = ",val)
    println("####################################")
    println("Solving time:")
    return val
end
    


function POP_DSMA(x::Vector{PolyVar{true}},f::Polynomial{true},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64;EigAlg::String="Arpack",maxit::Float64=1e5,tol::Float64=1e-4,UseEq::Bool=true,eps::Float64=1e-3)

    @time begin
    
    
    @time omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak=model_POP_DSMA(x,f,g,h,k,UseEq=UseEq)
    

      
    @time opt_val=solve_POP_DSMA(omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak,maxit=maxit,EigAlg=EigAlg,tol=tol,eps=eps)
    
    #@time opt_sol=extract_optimizer_moment_matrix2(X_sol,sk,vsk,n,m,l,opt_val,f,g,h,x)
                    
       println("Total time:")
                    end
    return opt_val
end
                    