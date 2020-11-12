function get_constant_trace_LMBM_arb_cons(x::Vector{PolyVar{true}},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64;UseEq::Bool=true)
    
    
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
    

    lV+=sk_g[m]*lmon_g[m]

    if UseEq
        if l!=0
            lV+=sum(lmon_h[i]*s2k_h[i] for i in 1:l)
        end
    end
    V=Matrix{UInt64}(undef,n,lV)
    V[:,1:sk]=2*v[:,1:sk]
    t=sk+1
    
    @fastmath @inbounds @simd for r in 1:lmon_g[m]
        @fastmath @inbounds @simd for j in 1:sk_g[m]
            V[:,t]=2*v[:,j]+supp_g[m][:,r]
            t+=1
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
    

    
    @fastmath @inbounds @simd for j in 1:sk
        add_to_expression!(cons[bfind(V,lV,2*v[:,j],n)],p0[j])
    end


    
    p=@variable(model, [1:sk_g[m]], lower_bound=1)
    @fastmath @inbounds @simd for j in 1:sk_g[m]
        @fastmath @inbounds @simd for r in 1:lmon_g[m]
            add_to_expression!(cons[bfind(V,lV,2*v[:,j]+supp_g[m][:,r],n)],p[j]*coe_g[m][r])
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
    
    G=value.(p0)
    
    Pg=Vector{Vector{Float64}}(undef,m)
    Pg[m]=sqrt.(value.(p))
    
    
    lsupp_other_cons=sum(sk_g[i]*lmon_g[i] for i in 1:m-1)
    
    
    supp_other_cons=Matrix{UInt64}(undef,n,lsupp_other_cons)

    t=1
    @fastmath @inbounds @simd for i in 1:m-1
        @fastmath @inbounds @simd for r in 1:lmon_g[i]
            @fastmath @inbounds @simd for j in 1:sk_g[i]
                supp_other_cons[:,t]=2*v[:,j]+supp_g[i][:,r]
                t+=1
            end
        end
    end
    
    supp_other_cons=unique(supp_other_cons,dims=2)
    supp_other_cons=sortslices(supp_other_cons,dims=2)
    lsupp_other_cons=size(supp_other_cons,2)
    
    
    coe_other_cons=zeros(Float64,lsupp_other_cons)
    
    @fastmath @inbounds @simd for i in 1:m-1
        @fastmath @inbounds @simd for j in 1:sk_g[i]
            @fastmath @inbounds @simd for r in 1:lmon_g[i]
                coe_other_cons[bfind(supp_other_cons,lsupp_other_cons,2*v[:,j]+supp_g[i][:,r],n)]+=coe_g[i][r]
            end
        end
    end


    
    return n,m,l,v,sk,sk_g,s2k_g,s2k_h,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,ak,G,Pg,lsupp_other_cons,supp_other_cons,coe_other_cons
end


function model_POP_LMBM_arb_cons(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64;UseEq::Bool=true)
    
    n,m,l,v,sk,sk_g,s2k_g,s2k_h,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,ak,G,Pg,lsupp_other_cons,supp_other_cons,coe_other_cons=get_constant_trace_LMBM_arb_cons(x,g,h,k,UseEq=UseEq)
    
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
    

    a=spzeros(Float64,d,zeta)
    
    
    @fastmath @inbounds for i in 1:sk, j in 1:i
        r=Order(v[:,i]+v[:,j])
        push!(IndM[r],[i;j])
        push!(invIndeM[i],t_iter)
        t_iter+=1
    end
    
    l_IndM=[length(IndM[r]) for r in 1:s2k]
    
    
 
         
    t_a=UInt64(1)
    
    
    @fastmath @inbounds @simd for r in 1:s2k
        if l_IndM[r]>1
            @fastmath @inbounds @simd for i in 2:l_IndM[r]
                I1,I2=IndM[r][1]
                a[invIndeM[I1][I2],t_a]=ak*((0.5*sqrt(2)-1)*(I2<I1)+1)
            
                I1,I2=IndM[r][i]
                a[invIndeM[I1][I2],t_a]=-ak*((0.5*sqrt(2)-1)*(I2<I1)+1)

                t_a+=1
             end
             IndM[r]=Vector{Int64}[IndM[r][1]]
         end
    end
    
    
    
    
    sqrtG=sqrt.(G)
    
    invP=zeros(Float64,sk,sk)
    W_tilde=zeros(Float64,sk,sk)
    I1,I2=UInt64(0),UInt64(0)
   @fastmath @inbounds @simd  for p in 1:lsupp_other_cons
        I1,I2=IndM[Order(supp_other_cons[:,p])][1]
        invP[I1,I2]=coe_other_cons[p]/((I2<I1)+1)
        invP[I2,I1]=invP[I1,I2]
        
        W_tilde[I1,I2]=invP[I1,I2]/sqrtG[I1]/sqrtG[I2]
        W_tilde[I2,I1]=W_tilde[I1,I2]
   end
    
    
    
    delta=1/(0.1+abs(eigs(W_tilde,nev = 1,which=:LR)[1][1]))
    
    invP=-delta*invP
    @fastmath @inbounds @simd  for j in 1:sk
        invP[j,j]+=G[j]
    end
    
 
    E=eigen(Symmetric(invP))
    
    
    
    invP=E.vectors*diagm(E.values.^(-0.5))*E.vectors'
    

    
    delta=sqrt(delta)
    @fastmath @inbounds @simd  for j in 1:m-1
        Pg[j]=delta*ones(Float64,sk_g[j])
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
                    a[invIndeMg[j][I1][I2],t_a]=ak/Pg[j][I1]/Pg[j][I2]*((0.5*sqrt(2)-1)*(I2<I1)+1)

                    I1,I2=IndMg[j][r][i]
                    a[invIndeMg[j][I1][I2],t_a]=-ak/Pg[j][I1]/Pg[j][I2]*((0.5*sqrt(2)-1)*(I2<I1)+1)
                   
                    t_a+=1
                end
                IndMg[j][r]=Vector{Int64}[IndMg[j][r][1]]
            end
            
                    I1,I2=IndMg[j][r][1]
                    a[invIndeMg[j][I1][I2],t_a]=-ak/Pg[j][I1]/Pg[j][I2]*((0.5*sqrt(2)-1)*(I2<I1)+1)
            
              
            @fastmath @inbounds @simd  for p in 1:lmon_g[j]  
                    I1,I2=IndM[Order(supp_g[j][:,p]+v[:,r])][end]
                    a[invIndeM[I1][I2],t_a]=coe_g[j][p]*ak*((0.5*sqrt(2)-1)*(I2<I1)+1)
                 
                   end
                   t_a+=1
            
         end
        
         
        
         t_Blo+=u_g[j]
    end        
            
  
    @fastmath @inbounds @simd for j in 1:l
        @fastmath @inbounds @simd  for r in 1:s2k_h[j]
            @fastmath @inbounds @simd  for p in 1:lmon_h[j]
                      I1,I2=IndM[Order(supp_h[j][:,p]+v[:,r])][1]
                       a[invIndeM[I1][I2],t_a]=coe_h[j][p]*ak*((0.5*sqrt(2)-1)*(I2<I1)+1)
                     
                   end
                    @inbounds t_a+=1
               end       
    end
    
    a[1,t_a]=ak

    
    a0=zeros(Float64,d)
    
    
    lmon_f,supp_f,coe_f=info(f,x,n)
    @fastmath @inbounds @simd  for p in 1:lmon_f
        I1,I2=IndM[Order(supp_f[:,p])][1]
        a0[invIndeM[I1][I2]]=coe_f[p]*ak*((0.5*sqrt(2)-1)*(I2<I1)+1)
   end
    
    
    for j in 1:zeta
        a[1:u,j]=get_vec_LMBM_arb_cons(invP*get_mat_LMBM_arb_cons(a[1:u,j],sk)*invP,sk,u)
    end
 
    a0[1:u]=get_vec_LMBM_arb_cons(invP*get_mat_LMBM_arb_cons(a0[1:u],sk,sparsemat=false)*invP,sk,u,sparsevec=false)
    
    a_ind1,a_ind2,a_val,a_len,a0,norm_a0,opnorm_a=rescale_LMBM_arb_cons(a,a0,zeta)
    
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
 
    return omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak,invP,v[:,1:sk],n,m,l
end

function get_mat_LMBM_arb_cons(vec,sk;sparsemat=true)
    if sparsemat
        B=spzeros(Float64,sk,sk)
    else
        B=zeros(Float64,sk,sk)
    end
    
    r=1
    @fastmath @inbounds for i in 1:sk, j in 1:i
        B[i,j]=vec[r]/sqrt(1+(j<i))
        B[j,i]= copy(B[i,j])
        r+=1
    end
    return B
end

function get_vec_LMBM_arb_cons(mat,sk,u;sparsevec=true)
    r=1
    if sparsevec
        vec=spzeros(Float64,u)
    else
        vec=zeros(Float64,u)
    end
    @fastmath @inbounds for i in 1:sk, j in 1:i
        vec[r]=mat[i,j]*sqrt(1+(j<i))
        r+=1
    end
    return vec
end
     

function rescale_LMBM_arb_cons(a::SparseMatrixCSC{Float64},a0::Vector{Float64},zeta::UInt64)
    
    a_ind1,a_ind2,a_val=findnz(a)
    a_len=length(a_ind1)

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
    
    
    return a_ind1,a_ind2,a_val,a_len,a0,norm_a0,opnorm_a
end



function SmallEig_LMBM_arb_cons(mat::Matrix{Float64},s::UInt64;EigAlg::String="Arpack")
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



                                        
function SmallEig_block_LMBM_arb_cons(vec::Adjoint{Float64,Array{Float64,1}},sk::UInt64;EigAlg::String="Arpack")
    
    B=zeros(Float64,sk,sk)
    r=1
    @fastmath @inbounds for i in 1:sk, j in 1:i
        B[i,j]=vec[r]/sqrt(1+(j<i))
        B[j,i]= copy(B[i,j])
        r+=1
    end
 
    return SmallEig_LMBM_arb_cons(B,sk,EigAlg=EigAlg)

end




function SmallEigBlocks_LMBM_arb_cons(a0_block::Vector{Adjoint{Float64,Array{Float64,1}}},a_block::Vector{SparseMatrixCSC{Float64,Int64}},y::Vector{Float64},s::Vector{UInt64},omega::Int64;EigAlg::String="Arpack")
    
    
    eigval,eigvec=Float64(Inf),zeros(Float64,1)

    smalleigval=eigval
    smalleigvec=eigvec
    ind=0

 
    @fastmath @inbounds @simd for p in 1:omega
        eigval,eigvec=SmallEig_block_LMBM_arb_cons(a0_block[p]+y'*a_block[p],s[p],EigAlg=EigAlg)
        if smalleigval>eigval
            smalleigval=eigval
            smalleigvec=eigvec
            ind=p
        end
    end

    return smalleigval,smalleigvec,ind
end       



function solve_POP_LMBM_arb_cons(omega::Int64,a0_block::Vector{Adjoint{Float64,Array{Float64,1}}},a_block::Vector{SparseMatrixCSC{Float64,Int64}},s::Vector{UInt64},zeta::UInt64,norm_a0::Float64,opnorm_a::Float64,ak::Float64;EigAlg::String="Arpack",tol::Float64=1e-4)
                
    function phi(nvar::Cint,xp::Ptr{Cdouble},gp::Ptr{Cdouble})
        zvar=unsafe_wrap(Array, xp, (convert(Int, nvar),))
        grad=unsafe_wrap(Array, gp, (convert(Int, nvar),))  
        eigval,eigvec,ind=SmallEigBlocks_LMBM_arb_cons(a0_block,a_block,zvar,s,omega,EigAlg=EigAlg)
        grad[:]=-a_block[ind]*[@fastmath @inbounds eigvec[i]*eigvec[j]*sqrt(1+(j<i)) for i in 1:s[ind] for j in 1:i]                
        grad[nvar]+=1/opnorm_a/ak
        return(convert(Cdouble,-eigval+zvar[nvar]/opnorm_a/ak))
    end

                        
                        
    opt_val,z= lmbm(phi,zeros(Float64,zeta);printinfo=true,tol=tol)
                     
                        
    opt_val*=-norm_a0
    println()
    println("####################################")
    println("opt_val = ",opt_val)
    println("####################################")
    println("Solving time:")
    return opt_val,z            
                
end
    


function POP_LMBM_arb_cons(x::Vector{PolyVar{true}},f::Polynomial{true},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64;EigAlg="Arpack",tol=1e-4,UseEq=true)

    @time begin
    
    
    @time omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak,invP,vsk,n,m,l=model_POP_LMBM_arb_cons(x,f,g,h,k,UseEq=UseEq)
    
 
    @time opt_val,z=solve_POP_LMBM_arb_cons(omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak,EigAlg=EigAlg,tol=tol)
    
    #@time opt_sol=extract_optimizer_LMBM_arb_cons(a0_block[1]+z'*a_block[1],P,s[1],vsk,n,m,l,opt_val,f,g,h,x,EigAlg=EigAlg) 
       println("Total time:")
                    end
    return opt_val#,opt_sol
end
                    
function get_Grammat_LMBM_arb_cons(vec::Adjoint{Float64,Array{Float64,1}},sk::UInt64,P::Vector{Float64};EigAlg="Arpack")
    t_iter=1
    Gr=zeros(Float64,sk,sk)
    @fastmath @inbounds for i in 1:sk, j in i:sk
        Gr[i,j]=vec[t_iter]/sqrt(1+(j<i))
        Gr[j,i]= copy(Gr[i,j])
        t_iter+=1
    end
    eigval,~=SmallEig_LMBM_arb_cons(-Gr,sk,EigAlg=EigAlg)
    @fastmath @inbounds @simd for i in 1:sk
        Gr[i,i]+=eigval
    end
    t_iter=1               
    @fastmath @inbounds for i in 1:sk, j in i:sk
        Gr[i,j]*=P[i]*P[j]
        Gr[j,i]= copy(Gr[i,j])
        t_iter+=1
    end                   
    return Gr
end
                
                
                
function extract_optimizer_LMBM_arb_cons(vec::Adjoint{Float64,Array{Float64,1}},P::Vector{Float64},lu0::UInt64,basis_sigma0::Matrix{UInt64},n::Int64,l_g::Int64,l_h::Int64,opt_val::Float64,f::Polynomial{true},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},x::Vector{PolyVar{true}};EigAlg="Arpack")
     
    Gr=get_Grammat_LMBM_arb_cons(vec,lu0,P,EigAlg=EigAlg)                
                    
    #extraction of optimizers
    V=nullspace(Gr,atol=1e-5)
    rk=size(V,2)

    println("Dimension of the null space of Gram matrix = ", rk)
    sol=Vector{Float64}[]
    if rk >0
        # extraction of Henrion and Lasserre
        V= rref_with_pivots!(Matrix(V'),1e-2)
        U=Matrix(V[1]')

        w=basis_sigma0[:,V[2]]
        N=Vector{Array{Float64}}(undef,n)
        flag=UInt8(1)
 
        
        for i in 1:n
            @inbounds kk=UInt16[]
            for j=1:size(w)[2]
                @inbounds xwj=w[:,j]
                @inbounds xwj[i]+=1
                for t=1:lu0
                    if xwj==basis_sigma0[:,t]
                        if kk==UInt16[]
                            @inbounds kk=t
                        else
                            @inbounds kk=[kk;t]
                        end
                        break
                    end
                end
            end
            
            if rk!=length(kk)
                @inbounds flag=0
                break
            else
                @inbounds N[i]=U[kk,:]
            end
        end

        if flag==1

            # Create random convex combination
            rands = rand(n);rands = rands/sum(rands)
            M = zeros(length(V[2]),rk)
            for i=1:n
                @inbounds M+=rands[i]*N[i]
            end

            F= schur(M);
            L=F.Z
            # Extract solution
            for i=1:rk
                @inbounds atom=Float64[]
                for j = 1:n
                    @inbounds coordinatej=L[:,i]'*N[j]*L[:,i]
                    @inbounds push!(atom,coordinatej[1])
                end
                println("------------------------------------")
                println("atom ",i," = ",atom)
                @inbounds flag=1
                
                @inbounds check=polynomial(f)(x => atom)-opt_val
                
                println("  check gap of lower bound  = ",check)
                if abs(check)>1e-1
                    @inbounds flag=0
                end


                for i=1:l_g
                    @inbounds check=polynomial(g[i])(x => atom)
                    println("  check inequality constraint ",i," = ",check)
                    if check<-1e-1
                        @inbounds flag=0
                    end
                end
                
                for i=1:l_h
                    @inbounds check=polynomial(h[i])(x => atom)
                    println("  check equality constraint ",i," = ",check)
                    if abs(check)>1e-1
                        @inbounds flag=0
                    end
                end

                if flag ==1
                    @inbounds sol=atom
                    println("####################################")
                    println("Optimal solution: opt_sol = ",atom)
                    println("####################################")
                end

            end
        end
    end
    
    println("Extracting solutuion time:")

    return sol

end
