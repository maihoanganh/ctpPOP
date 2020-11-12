function model_on_each_clique(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64;first_clique=false,final_clique=false)
    
    n=length(x)
    m=length(g)
    l=length(h)
    
    v=get_basis(n,2*k)
    s2k=size(v,2)
    sort_v=sortslices(v,dims=2)
    re_ind=Vector{UInt64}(undef,s2k)
    @simd for j in 1:s2k
        @inbounds re_ind[bfind(sort_v,s2k,v[:,j],n)]=j
    end
    sk=binomial(k+n,n)      
    Order(alpha::Vector{UInt64})=re_ind[bfind(sort_v,s2k,alpha,n)]
    
    #basis on intersection of variables
    v_1=get_basis(1,2*k)
    s2k_1=size(v_1,2)
            
    #Need to improve
    ceil_g=[@inbounds ceil(Int64,maxdegree(g[i])/2) for i in 1:m]
    sk_g=[@inbounds binomial(k-ceil_g[i]+n,n) for i in 1:m]
    s2k_g=[@inbounds binomial(2*(k-ceil_g[i])+n,n) for i in 1:m]
    s2k_h=[@inbounds binomial(2*(k-ceil(Int64,maxdegree(h[j])/2))+n,n) for j in 1:l]

    #println("  Number of blocks: omega=",m+1)
    #println("  Size of the largest block: s^max=",sk)
    
    u=Int64(0.5*sk*(sk+1))
    u_g=[@inbounds Int64(0.5*sk_g[i]*(sk_g[i]+1)) for i in 1:m]
    d=u+sum(u_g)
    
    if first_clique
        zeta=d-s2k+sum(s2k_h)+1+s2k_1  
    elseif final_clique
        zeta=d-s2k+sum(s2k_h)+s2k_1
    else
        zeta=d-s2k+sum(s2k_h)+2*s2k_1
    end
    #println("  Number of equality trace constraints: zeta=",zeta) 
            
    ak=2^k
    #println("  Constant trace: ak=",ak)
    

    IndM=[@inbounds Vector{Int64}[] for j in 1:s2k]
    invIndeM=spzeros(UInt64,sk,sk)
    diagM=Vector{UInt64}(undef,sk)
    r=UInt64(0)
    t_iter=UInt64(1)
    
    a=spzeros(Float64,d,zeta) 
    
            
    for i in 1:sk, j in i:sk
        @inbounds r=Order(v[:,i]+v[:,j])
        @inbounds append!(IndM[r],[[i,j]])
        @inbounds invIndeM[i,j]=t_iter
        if i==j
            @inbounds diagM[i]=r
        end
        @inbounds t_iter+=1
    end
    l_IndM=[length(IndM[r]) for r in 1:s2k]
    
    lmon,supp,coe=info((1.0+sum(x.^2))^k,x,n)
    invP0=invdiaP2([@inbounds Order(supp[:,j]) for j in 1:lmon],coe,sk,diagM)
            
    
    t_a=1
    if first_clique
        a[1,t_a]=-ak
        t_a+=1
    else
        t_a=s2k_1+1
    end
    
    I=zeros(UInt64,2)
            
    for r in 1:s2k
        if l_IndM[r]>1
            for i in 2:l_IndM[r]
                I=IndM[r][1]
                if I[1]==I[2]
                    a[invIndeM[I[1],I[2]],t_a]=invP0[I[1]]^2*ak
                else
                    a[invIndeM[I[1],I[2]],t_a]=0.5*invP0[I[1]]*invP0[I[2]]*ak
                end
                
                I=IndM[r][i]
                if I[1]==I[2]
                    a[invIndeM[I[1],I[2]],t_a]=-invP0[I[1]]^2*ak
                else
                    a[invIndeM[I[1],I[2]],t_a]=-0.5*invP0[I[1]]*invP0[I[2]]*ak
                end
                t_a+=1
             end
             IndM[r]=Vector{Int64}[IndM[r][end]]
         end
    end
   
    
    diagMg=Vector{Vector{UInt64}}(undef,m)
    invPg=Vector{Vector{Float64}}(undef,m)
    
    IndMg=[[Vector{Int64}[] for j in 1:s2k_g[i]] for i in 1:m]
    invIndeMg=[spzeros(UInt64,sk_g[i],sk_g[i]) for i in 1:m]
    t_Blo=u
 
    for j in 1:m
        r=UInt64(0)
        t_iter=UInt64(1)
        diagMg[j]=Vector{UInt64}(undef,sk_g[j])
        for p in 1:sk_g[j], q in p:sk_g[j]
            @inbounds r=Order(v[:,p]+v[:,q])
            @inbounds append!(IndMg[j][r],[[p,q]])
            if p==q
                @inbounds diagMg[j][p]=r
            end
            @inbounds invIndeMg[j][p,q]=t_iter
            t_iter+=1
        end
                
        lmon,supp,coe=info(sum(2^j*(1.0+sum(x.^2))^(k-j-1) for j in 0:k-1),x,n)
        invPg[j]=invdiaP2([@inbounds Order(supp[:,i]) for i in 1:lmon],coe,sk_g[j],diagMg[j])
    
        l_IndM=[length(IndMg[j][q]) for q in 1:s2k_g[j]]

        for r in 1:s2k_g[j]
            if l_IndM[r]>1
                for i in 2:l_IndM[r]
                    I=IndMg[j][r][1]
                    if I[1]==I[2]
                        a[invIndeMg[j][I[1],I[2]]+t_Blo,t_a]=invPg[j][I[1]]^2*ak
                    else
                        a[invIndeMg[j][I[1],I[2]]+t_Blo,t_a]=0.5*invPg[j][I[1]]*invPg[j][I[2]]*ak
                    end
                    I=IndMg[j][r][i]
                    if I[1]==I[2]
                        a[invIndeMg[j][I[1],I[2]]+t_Blo,t_a]=-invPg[j][I[1]]^2*ak
                    else
                        a[invIndeMg[j][I[1],I[2]]+t_Blo,t_a]=-0.5*invPg[j][I[1]]*invPg[j][I[2]]*ak
                    end
                    t_a+=1
                end
                IndMg[j][r]=Vector{Int64}[IndMg[j][r][end]]
            end
         end
         t_Blo+=u_g[j]
    end
  
    
    lmon,supp,coe=Int64(0),zeros(UInt64,1,1),zeros(Float64,1)
    
    t_Blo=u
    @simd for j in 1:m
        @inbounds lmon,supp,coe=info(g[j],x,n)
        @simd for r in 1:s2k_g[j]
                   @inbounds I=IndMg[j][Order(v[:,r])][1]
                  if I[1]==I[2]
                      @inbounds a[invIndeMg[j][I[1],I[2]]+t_Blo,t_a]=-invPg[j][I[1]]^2*ak
                  else
                      @inbounds a[invIndeMg[j][I[1],I[2]]+t_Blo,t_a]=-0.5*invPg[j][I[1]]*invPg[j][I[2]]*ak
                  end
            @simd for p in 1:lmon  
                      @inbounds I=IndM[Order(supp[:,p]+v[:,r])][1]
                      if I[1]==I[2]
                          @inbounds a[invIndeM[I[1],I[2]],t_a]=coe[p]*invP0[I[1]]^2*ak
                      else
                          @inbounds a[invIndeM[I[1],I[2]],t_a]=0.5*coe[p]*invP0[I[1]]*invP0[I[2]]*ak
                      end
                   end
                   t_a+=1
               end
        t_Blo+=u_g[j]            
    end             

    @simd for j in 1:l
        @inbounds lmon,supp,coe=info(h[j],x,n)
        @simd for r in 1:s2k_h[j]
            @simd for p in 1:lmon  
                      @inbounds I=IndM[Order(supp[:,p]+v[:,r])][1]
                      if I[1]==I[2]
                          @inbounds a[invIndeM[I[1],I[2]],t_a]=coe[p]*invP0[I[1]]^2*ak
                      else
                          @inbounds a[invIndeM[I[1],I[2]],t_a]=0.5*coe[p]*invP0[I[1]]*invP0[I[2]]*ak
                      end
                   end
                    @inbounds t_a+=1
               end       
    end 
    
    
    
    alpha=zeros(UInt64,n)
    if !final_clique
        for i in 1:s2k_1
            alpha[n]=v_1[1,i]
            @inbounds I=IndM[Order(alpha)][1]
            if I[1]==I[2]
                @inbounds a[invIndeM[I[1],I[2]],t_a]=invP0[I[1]]^2*ak
            else
                @inbounds a[invIndeM[I[1],I[2]],t_a]=0.5*invP0[I[1]]*invP0[I[2]]*ak
            end
            @inbounds t_a+=1
        end
    end
    
    if !first_clique
        alpha[n]=0
        t_a=1

        for i in 1:s2k_1
            alpha[1]=v_1[1,i]
            @inbounds I=IndM[Order(alpha)][1]
            if I[1]==I[2]
                @inbounds a[invIndeM[I[1],I[2]],t_a]=-invP0[I[1]]^2*ak
            else
                @inbounds a[invIndeM[I[1],I[2]],t_a]=-0.5*invP0[I[1]]*invP0[I[2]]*ak
            end
            @inbounds t_a+=1
        end
    end

    a0=zeros(Float64,d)
    
    
    lmon,supp,coe=info(f,x,n)
    @simd for p in 1:lmon
        @inbounds I=IndM[Order(supp[:,p])][1]
        if I[1]==I[2]
          @inbounds a0[invIndeM[I[1],I[2]]]=-coe[p]*invP0[I[1]]^2*ak
        else
          @inbounds a0[invIndeM[I[1],I[2]]]=-0.5*coe[p]*invP0[I[1]]*invP0[I[2]]*ak
        end
   end
    
    Ind=Vector{UnitRange{Int64}}(undef,m+2)
    t_Blo=u
    for j in 1:m
        Ind[j]=1+t_Blo:u_g[j]+t_Blo
        t_Blo+=u_g[j]
    end
    Ind[m+1]=1:u      
    
    return m,a0,a,[sk_g;sk],u,u_g,zeta,Ind
end

function model_correlative_POP_single_ball(x::Vector{PolyVar{true}},f::Vector{Polynomial{true,Float64}},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64,I::Vector{UnitRange{Int64}},J::Vector{UnitRange{Int64}},W::Vector{UnitRange{Int64}};)
    
    p=length(I)
    
    s2k_1=binomial(2*k+1,1)
    
   
    m=Vector{Int64}(undef,p)
    a0=Vector{Vector{Float64}}(undef,p)
    a=Vector{SparseMatrixCSC{Float64}}(undef,p)
    IndClique=Vector{UnitRange{Int64}}(undef,p)
    s=Vector{Vector{Int64}}(undef,p)
    u=Vector{Int64}(undef,p)
    u_g=Vector{Vector{Int64}}(undef,p)
    
    zeta=Vector{Int64}(undef,p)
    Ind=Vector{Vector{UnitRange{Int64}}}(undef,p)
    #println("Clique 1:") 
    m[1],a0[1],a[1],s[1],u[1],u_g[1],zeta[1],Ind[1]=model_on_each_clique(x[I[1]],f[1],g[J[1]],h[W[1]],k,first_clique=true,final_clique=false)
    #println("----------------------------")
    for j in 2:p-1
        #println("Clique ",j,":") 
        @inbounds m[j],a0[j],a[j],s[j],u[j],u_g[j],zeta[j],Ind[j]=model_on_each_clique(x[I[j]],f[j],g[J[j]],h[W[j]],k,first_clique=false,final_clique=false)
        #println("----------------------------")
    end
    #println("Clique ",p,":")
    m[p],a0[p],a[p],s[p],u[p],u_g[p],zeta[p],Ind[p]=model_on_each_clique(x[I[p]],f[p],g[J[p]],h[W[p]],k,first_clique=false,final_clique=true)
    
    #println("----------------------------")
    
    println("Modeling part:")
    
    println("  Number of blocks: omega=",sum(m)+p)
    println("  Size of the largest block: s^max=",maximum([s[j][m[j]+1] for j in 1:p]))
    
    
    IndClique[1]=1:zeta[1]
    t_Cliq=zeta[1]-s2k_1

    
    
    for j in 2:p
        IndClique[j]=1+t_Cliq:zeta[j]+t_Cliq
        t_Cliq+=zeta[j]-s2k_1
    end
    
    zeta0=t_Cliq+s2k_1
    println("  Number of equality trace constraints: zeta=",zeta0) 
    
    
    return m,a0,a,s,u,u_g,zeta0,Ind,IndClique,p
end


function LargEigBlock(vec::Vector{Float64},s::Vector{Int64},u::Int64,u_g::Vector{Int64},m::Int64;EigAlg="Arpack")
    
    
    eigval,eigvec=LargEig_block(vec[1:u],s[m+1],EigAlg=EigAlg)

    largeigval=eigval
    largeigvec=eigvec
    ind=m+1
    
    

    t_Blo=u
    for p in 1:m
        @inbounds eigval,eigvec=LargEig_block(vec[t_Blo+1:t_Blo+u_g[p]],s[p],EigAlg=EigAlg)
        if largeigval<eigval
            @inbounds largeigval=eigval
            @inbounds largeigvec=eigvec
            @inbounds ind=p
        end
        t_Blo+=u_g[p]
    end

    return largeigval,largeigvec,ind
end


function eigcom_cliques(y::Vector{Float64},p,m,a0,a,s,u,u_g,IndClique;EigAlg="Arpack")
    eigval=Vector{Float64}(undef,p)
    eigvec=Vector{Vector{Float64}}(undef,p)
    ind=Vector{Int64}(undef,p)

    for j in 1:p eigval[j],eigvec[j],ind[j]=LargEigBlock(a0[j]+a[j]*y[IndClique[j]],s[j],u[j],u_g[j],m[j],EigAlg=EigAlg)
    end
    return eigval,eigvec,ind
end
function solveNSP_correlative_POP_single_ball(p,m,a0,a,s,u,u_g,zeta0,Ind,IndClique;tol=1e-3,EigAlg="Arpack")
    println("**LMBM solver:")
    function phi(nvar::Cint,xp::Ptr{Cdouble},gp::Ptr{Cdouble})
        zvar=unsafe_wrap(Array, xp, (convert(Int, nvar),))
        grad=unsafe_wrap(Array, gp, (convert(Int, nvar),))  
        eigval,eigvec,ind=eigcom_cliques(zvar,p,m,a0,a,s,u,u_g,IndClique,EigAlg=EigAlg)
        grad[:]=zeros(Float64,zeta0)
        grad[1]+=1
        for j in 1:p
            grad[IndClique[j]]+=a[j][Ind[j][ind[j]],:]'*[@inbounds eigvec[j][i]*eigvec[j][r]*(2-0^(r-i)) for i in 1:s[j][ind[j]] for r in i:s[j][ind[j]]]   
        end
        return (convert(Cdouble,sum(eigval)+zvar[1]))
    end

                        
                        
    opt_val,~= lmbm(phi,zeros(Float64,zeta0);printinfo=false,tol=tol)
                     
                        
    opt_val*=-1
 
    println("####################################")
    println("opt_val = ",opt_val)
    println("####################################")
    return opt_val
end

function correlative_POP_single_ball(x::Vector{PolyVar{true}},f::Vector{Polynomial{true,Float64}},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64,I::Vector{UnitRange{Int64}},J::Vector{UnitRange{Int64}},W::Vector{UnitRange{Int64}};EigAlg="Arpack",tol=1e-3)
    
    @time begin
    
    
    @time m,a0,a,s,u,u_g,zeta0,Ind,IndClique,p=model_correlative_POP_single_ball(x,f,g,h,k,I,J,W)
                
    @time opt_val=solveNSP_correlative_POP_single_ball(p,m,a0,a,s,u,u_g,zeta0,Ind,IndClique,tol=tol,EigAlg=EigAlg)
    
                    end
    return opt_val
end
                    
                    
function invdiaP(supp_theta::Vector{Int64}, coe_theta::Vector{Int64},size_mommat::Int64,diagu::Array{Int64,1})
    
    D=zeros(Float64,size_mommat)
    
    
    @simd for j in 1:size_mommat
        @inbounds i=findfirst(x->x==diagu[j],supp_theta)
        if supp_theta[i]==diagu[j]
            @inbounds D[j]=1/sqrt(coe_theta[i])
        end
    end
    
    return D
end

function invdiaP2(supp_theta::Vector{UInt64}, coe_theta::Vector{Float64},size_mommat::Int64,diagu::Array{UInt64,1})
    
    D=zeros(Float64,size_mommat)
    
    
    @simd for j in 1:size_mommat
        @inbounds i=findfirst(x->x==diagu[j],supp_theta)
        if supp_theta[i]==diagu[j]
            @inbounds D[j]=1/sqrt(coe_theta[i])
        end
    end
    
    return D
end
                    
                    
function LargEig(mat::Matrix{Float64},s::Int64;EigAlg="Arpack",showEvaluation=false)
    if showEvaluation
        global num_eig+=1
        global max_size=maximum([max_size,s])
    end
    if s==1
        return mat[1,1],ones(Float64,1)
    elseif EigAlg=="Arpack"
       E=eigs(mat,nev = 1,which=:LR) 
       return E[1][1],E[2][:,1]
    elseif EigAlg=="ArnoldiMethod"      
       E=partialeigen(partialschur(mat, nev=1, which=LR())[1])
       return E[1][1],E[2][:,1]
    elseif EigAlg=="Normal"
       E=eigen(Symmetric(mat),s:s)
       return E.values[1],E.vectors[:,1]
    elseif EigAlg=="Mix"
       try 
           E=partialeigen(partialschur(mat, nev=1,tol=1e-2, which=LR())[1])
           return E[1][1],E[2][:,1]
       catch
           E=eigs(mat,nev = 1,which=:LR,tol=1e-2) 
           return E[1][1],E[2][:,1]
       end
    else
       println("No eigenvalue algorithm!!!")
    end  
    
end
                    
                    
function LargEig_block(vec,sk;EigAlg="Arpack")
    B=zeros(Float64,sk,sk)
    t=1
    for i in 1:sk, j in i:sk
        @inbounds B[i,j]=vec[t]
        @inbounds B[j,i]= copy(B[i,j])
        t+=1
    end
    return LargEig(B,sk,EigAlg=EigAlg)
end


function LargEigBlocks(vec::Vector{Float64},sk::Int64,sk_g::Vector{Int64},u::Int64,u_g::Vector{Int64},m::Int64;EigAlg="Arpack")
    
    
    eigval,eigvec=LargEig_block(vec[1:u],sk,EigAlg=EigAlg)

    largeigval=eigval
    largeigvec=eigvec
    ind=m+1
    size=sk

    t_Blo=u
    for p in 1:m
        @inbounds eigval,eigvec=LargEig_block(vec[t_Blo+1:t_Blo+u_g[p]],sk_g[p],EigAlg=EigAlg)
        if largeigval<eigval
            @inbounds largeigval=eigval
            @inbounds largeigvec=eigvec
            @inbounds ind=p
            @inbounds size=sk_g[p]
        end
        t_Blo+=u_g[p]
    end

    return largeigval,largeigvec,ind,size
end       