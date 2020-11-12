function get_constant_trace_pmsv(x,g,h,k)
    n=length(x)
    m=length(g)
    l=length(h)
            
    v=get_basis(n,2*k)
    
    sk=binomial(k+n,n)
    sk_g=Vector{Int64}(undef,m)
    s2k_g=Vector{Int64}(undef,m)
    s2k_h=Vector{Int64}(undef,l)
    
    lmon_g=Vector{UInt64}(undef,m)
    supp_g=Vector{Matrix{UInt64}}(undef,m)
    coe_g=Vector{Vector{Float64}}(undef,m)
    
    lmon_h=Vector{UInt64}(undef,l)
    supp_h=Vector{Matrix{UInt64}}(undef,l)
    coe_h=Vector{Vector{Float64}}(undef,l)
        
    
    ceil_g=Int64(0)
    for i in 1:m
        ceil_g=ceil(Int64,maxdegree(g[i])/2)
        sk_g[i]=binomial(k-ceil_g+n,n)
        s2k_g[i]=binomial(2*(k-ceil_g)+n,n)
        lmon_g[i],supp_g[i],coe_g[i]=info(g[i],x,n)
    end
                    
    for i in 1:l
        s2k_h[i]=binomial(2*(k-ceil(Int64,maxdegree(h[i])/2))+n,n)
        lmon_h[i],supp_h[i],coe_h[i]=info(h[i],x,n)
    end
    
    lV=sk
    if m!=0
        lV+=sum(sk_g[i]*lmon_g[i] for i in 1:m)
    end
    if l!=0
        lV+=sum(lmon_h[i]*s2k_h[i] for i in 1:l)
    end
     
    V=Matrix{UInt64}(undef,n,lV)
    V[:,1:sk]=2*v[:,1:sk]
    t=sk+1
    for i in 1:m
        for r in 1:lmon_g[i]
            for j in 1:sk_g[i]
                @inbounds V[:,t]=2*v[:,j]+supp_g[i][:,r]
                @inbounds t+=1
            end
        end
    end
                    
    for i in 1:l
        for r in 1:lmon_h[i]
            for j in 1:s2k_h[i]
                @inbounds V[:,t]=v[:,j]+supp_h[i][:,r]
                @inbounds t+=1
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
    q=Vector{Vector{VariableRef}}(undef, l)

    
    for j in 1:sk
        @inbounds add_to_expression!(cons[bfind(V,lV,2*v[:,j],n)],p0[j])
    end


    for i in 1:m
        p[i]=@variable(model, [1:sk_g[i]], lower_bound=1)
        for j in 1:sk_g[i]
            for r in 1:lmon_g[i]
                @inbounds add_to_expression!(cons[bfind(V,lV,2*v[:,j]+supp_g[i][:,r],n)],p[i][j]*coe_g[i][r])
            end
        end
    end

    for i in 1:l
        q[i]=@variable(model, [1:s2k_h[i]])
        for j in 1:s2k_h[i]
            for r in 1:lmon_h[i]
                @inbounds add_to_expression!(cons[bfind(V,lV,v[:,j]+supp_h[i][:,r],n)],q[i][j]*coe_h[i][r])
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


function model_POP_ball(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64)
    
    n,m,l,v,sk,sk_g,s2k_g,s2k_h,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,ak,P,Pg=get_constant_trace(x,g,h,k)
    
    s2k=size(v,2)
    sort_v=sortslices(v,dims=2)
    re_ind=Vector{UInt64}(undef,s2k)
    @simd for j in 1:s2k
        @inbounds re_ind[bfind(sort_v,s2k,v[:,j],n)]=j
    end
        
    Order(alpha::Vector{UInt64})=re_ind[bfind(sort_v,s2k,alpha,n)]        
    
    println("  Number of blocks: omega=",m+1)
    println("  Size of the largest block: s^max=",sk)
    
    u=Int64(0.5*sk*(sk+1))
    u_g=[@inbounds Int64(0.5*sk_g[i]*(sk_g[i]+1)) for i in 1:m]
    d=u+sum(u_g)
    zeta=d-s2k+sum(s2k_h)+1   
    println("  Number of equality trace constraints: zeta=",zeta) 
            

    IndM=[@inbounds Vector{Int64}[] for j in 1:s2k]
    invIndeM=spzeros(UInt64,sk,sk)
    r=UInt64(0)
    t_iter=UInt64(1)

    a=spzeros(Float64,d,zeta)    
            
    for i in 1:sk, j in i:sk
        @inbounds r=Order(v[:,i]+v[:,j])
        @inbounds append!(IndM[r],[[i,j]])
        @inbounds invIndeM[i,j]=t_iter
        @inbounds t_iter+=1
    end
    l_IndM=[length(IndM[r]) for r in 1:s2k]
    
         
    t_a=UInt64(1)
    I=zeros(UInt64,2)
            
    for r in 1:s2k
        if l_IndM[r]>1
            for i in 2:l_IndM[r]
                I=IndM[r][1]
                if I[1]==I[2]
                    a[invIndeM[I[1],I[2]],t_a]=ak/P[I[1]]^2
                else
                    a[invIndeM[I[1],I[2]],t_a]=0.5*ak/P[I[1]]/P[I[2]]*sqrt(2)
                end
                
                I=IndM[r][i]
                if I[1]==I[2]     
                    a[invIndeM[I[1],I[2]],t_a]=-ak/P[I[1]]^2
                else
                    a[invIndeM[I[1],I[2]],t_a]=-0.5*ak/P[I[1]]/P[I[2]]*sqrt(2)
                end
                t_a+=1
             end
             IndM[r]=Vector{Int64}[IndM[r][end]]
         end
    end
   
    
    
    invPg=Vector{Vector{Float64}}(undef,m)
    
    IndMg=[[Vector{Int64}[] for j in 1:s2k_g[i]] for i in 1:m]
    invIndeMg=[spzeros(UInt64,sk_g[i],sk_g[i]) for i in 1:m]
    t_Blo=u
 
    for j in 1:m
        r=UInt64(0)
        t_iter=UInt64(1)
        for p in 1:sk_g[j], q in p:sk_g[j]
            @inbounds r=Order(v[:,p]+v[:,q])
            @inbounds append!(IndMg[j][r],[[p,q]])
            @inbounds invIndeMg[j][p,q]=t_iter
            t_iter+=1
        end
                
        
        l_IndM=[length(IndMg[j][q]) for q in 1:s2k_g[j]]

        for r in 1:s2k_g[j]
            if l_IndM[r]>1
                for i in 2:l_IndM[r]
                    I=IndMg[j][r][1]
                    if I[1]==I[2]
                        @inbounds a[invIndeMg[j][I[1],I[2]]+t_Blo,t_a]=ak/Pg[j][I[1]]^2
                    else
                        @inbounds a[invIndeMg[j][I[1],I[2]]+t_Blo,t_a]=0.5*ak/Pg[j][I[1]]/Pg[j][I[2]]*sqrt(2)
                    end
                    I=IndMg[j][r][i]
                    if I[1]==I[2]
                        @inbounds a[invIndeMg[j][I[1],I[2]]+t_Blo,t_a]=-ak/Pg[j][I[1]]^2
                    else
                        @inbounds a[invIndeMg[j][I[1],I[2]]+t_Blo,t_a]=-0.5*ak/Pg[j][I[1]]/Pg[j][I[2]]*sqrt(2)
                    end
                    t_a+=1
                end
                IndMg[j][r]=Vector{Int64}[IndMg[j][r][end]]
            end
         end
         t_Blo+=u_g[j]
    end
   
    
            
            
            
      
    t_Blo=u
    @simd for j in 1:m
        @simd for r in 1:s2k_g[j]
                   @inbounds I=IndMg[j][Order(v[:,r])][1]
                  if I[1]==I[2]
                      @inbounds a[invIndeMg[j][I[1],I[2]]+t_Blo,t_a]=-ak/Pg[j][I[1]]^2
                  else
                      @inbounds a[invIndeMg[j][I[1],I[2]]+t_Blo,t_a]=-0.5*ak/Pg[j][I[1]]/Pg[j][I[2]]*sqrt(2)
                  end
            @simd for p in 1:lmon_g[j]  
                      @inbounds I=IndM[Order(supp_g[j][:,p]+v[:,r])][1]
                      if I[1]==I[2]
                          @inbounds a[invIndeM[I[1],I[2]],t_a]=coe_g[j][p]*ak/P[I[1]]^2
                      else
                          @inbounds a[invIndeM[I[1],I[2]],t_a]=0.5*coe_g[j][p]*ak/P[I[1]]/P[I[2]]*sqrt(2)
                      end
                   end
                   t_a+=1
               end
        t_Blo+=u_g[j]            
    end             
            
    
    
  
            
            

    @simd for j in 1:l
        @inbounds lmon,supp,coe=info(h[j],x,n)
        @simd for r in 1:s2k_h[j]
            @simd for p in 1:lmon_h[j]
                      @inbounds I=IndM[Order(supp_h[j][:,p]+v[:,r])][1]
                      if I[1]==I[2]
                          @inbounds a[invIndeM[I[1],I[2]],t_a]=coe_h[j][p]*ak/P[I[1]]^2
                      else
                          @inbounds a[invIndeM[I[1],I[2]],t_a]=0.5*coe_h[j][p]*ak/P[I[1]]/P[I[2]]*sqrt(2)
                      end
                   end
                    @inbounds t_a+=1
               end       
    end 
    
    
    a[1,t_a]=-ak
            


    a0=zeros(Float64,d)
    
    
    lmon_f,supp_f,coe_f=info(f,x,n)
    @simd for p in 1:lmon_f
        @inbounds I=IndM[Order(supp_f[:,p])][1]
        if I[1]==I[2]
          @inbounds a0[invIndeM[I[1],I[2]]]=-coe_f[p]*ak/P[I[1]]^2
        else
          @inbounds a0[invIndeM[I[1],I[2]]]=-0.5*coe_f[p]*ak/P[I[1]]/P[I[2]]*sqrt(2)
        end
   end
    
    
    Ind=Vector{Vector{Int64}}(undef,m+2)
    t_Blo=u
    for j in 1:m
        Ind[j]=1+t_Blo:u_g[j]+t_Blo
        t_Blo+=u_g[j]
    end
    Ind[m+1]=1:u
    
    for j in 1:zeta
        a[:,j]=a[:,j]./norm(a[:,j])
    end
    
    opnorm_a=opnorm(Matrix(a'), 2)
    
    a=a./opnorm_a

    
    norm_a0=norm(a0)
    a0=a0./norm(a0)
   
    println("Modeling time:")
 
    return n,m,l,a0,a,v[:,1:sk],Ind,sk,sk_g,P,zeta,u,u_g,norm_a0,opnorm_a,ak
end
                                        

                                        
function LargEig_block(vec,sk;EigAlg="Arpack")
    B=zeros(Float64,sk,sk)
    t=1
    for i in 1:sk, j in i:sk
        @inbounds B[i,j]=vec[t]/sqrt(1+(j>i))
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


function solve_POP_ball(m,a0,a,Ind,sk,sk_g,u,u_g,zeta,norm_a0,opnorm_a,ak;EigAlg="Arpack",tol=1e-3)

    function phi(nvar::Cint,xp::Ptr{Cdouble},gp::Ptr{Cdouble})
        zvar=unsafe_wrap(Array, xp, (convert(Int, nvar),))
        grad=unsafe_wrap(Array, gp, (convert(Int, nvar),))  
        eigval,eigvec,ind,s=LargEigBlocks(a0+a*zvar,sk,sk_g,u,u_g,m,EigAlg=EigAlg)
        grad[:]=a[Ind[ind],:]'*[@inbounds eigvec[i]*eigvec[j]*(2-0^(j-i)) for i in 1:s for j in i:s]                
        grad[nvar]+=1/opnorm_a/ak
        return(convert(Cdouble,eigval+zvar[nvar]/opnorm_a/ak))
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
    


function POP_ball(x::Vector{PolyVar{true}},f::Polynomial{true},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64;EigAlg="Arpack",tol=1e-3)

    @time begin
    
    
    @time n,m,l,a0,a,vsk,Ind,sk,sk_g,P,zeta,u,u_g,norm_a0,opnorm_a,ak=model_POP_ball(x,f,g,h,k)
    

      
    @time opt_val,z=solve_POP_ball(m,a0,a,Ind,sk,sk_g,u,u_g,zeta,norm_a0,opnorm_a,ak,tol=tol,EigAlg=EigAlg)
    
    @time opt_sol=extract_optimizer2(get_Grammat((a0+a*z)[1:u],sk,P,EigAlg=EigAlg),sk,vsk,n,m,l,opt_val,f,g,h,x) 
       println("Total time:")
                    end
    return opt_val,opt_sol
end
                    

                    
                    
function get_Grammat(vec::Vector{Float64},sk::Int64,P::Vector{Float64};EigAlg="Arpack")
    t_iter=1
    Gr=zeros(Float64,sk,sk)
    for i in 1:sk, j in i:sk
        @inbounds Gr[i,j]=vec[t_iter]
        @inbounds Gr[j,i]= copy(Gr[i,j])
        @inbounds t_iter+=1
    end
    eigval,~=LargEig(Gr,sk,EigAlg=EigAlg)
    for i in 1:sk
        @inbounds Gr[i,i]-=eigval
    end
    matP=diagm(P)                    
    return matP*Gr*matP
end