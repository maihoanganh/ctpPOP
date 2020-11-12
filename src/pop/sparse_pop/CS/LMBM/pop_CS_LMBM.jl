function model_CS_LMBM(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{SparseMatrixCSC{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{SparseMatrixCSC{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::SparseMatrixCSC{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64;UseEq=true)
    
    
    
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
    
    #println("  Number of cliques: p=", p)
    
    #=println(" I=",I)
    println(" J=",J)
    println(" W=",W)=#
    
    mex,lt_mex,l_mex,mex_cliq,lmex_cliq,v,s2k,sort_v,re_ind=get_mex_CS_LMBM(k,n,lI,p,I)
    
    Indf,lIndf=decomp_obj(supp_f,lmon_f,I,p,n)
    
     omega_cliq,a0_cliq,a_ind1_cliq,a_ind2_cliq,a_val_cliq,a_len_cliq,a_mex_ind,a_mex_val,r_cliq,u_cliq,s_cliq,zeta_cliq,d_cliq,ak=run_model_cliq_CS_LMBM(p,lI,lJ,lW,I,J,W,lIndf,Indf,supp_f,coe_f,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,dg,dh,v,s2k,sort_v,re_ind,mex_cliq,lmex_cliq,k,UseEq=UseEq)

    
    d=sum(d_cliq)
    zeta=sum(zeta_cliq)+sum(l_mex)-lt_mex+1
    
  
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
        
        for i in 1:a_len_cliq[j]
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
        
        for i in 2:l_mex[j]
            
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
    
    
    #@time test_Mosek_CS(sparse(a_ind1,a_ind2,a_val),a0,d,p,omega_cliq,s_cliq,zeta)
    
    a_val,a0,norm_a0,opnorm_a=rescale_CS_LMBM(a_ind1,a_ind2,a_val,a_len,a0,zeta)
    
   
    #@time test_Mosek_CS(sparse(a_ind1,a_ind2,a_val),a0,d,p,omega_cliq,s_cliq,zeta,norm_a0,opnorm_a,ak)
    
   
    Ind=Vector{Vector{Vector{UInt64}}}(undef,p)
    a_block=Vector{Vector{SparseMatrixCSC{Float64}}}(undef,p)
    a0_block=Vector{Vector{Adjoint{Float64,Array{Float64,1}}}}(undef,p)
    Index=Vector{UInt64}(1:a_len)
    pos=zeros(UInt64,1)
 
    for i in 1:p
        
        Ind[i]=[Vector{UInt64}([]) for j in 1:omega_cliq[i]]
        for j in 1:omega_cliq[i]
            pos=findall(t->(a_ind1[t]>=(r[i]+r_cliq[i][j]+1)) && (a_ind1[t]<=(r[i]+r_cliq[i][j]+u_cliq[i][j])),Index)
            Ind[i][j]=Index[pos]
            deleteat!(Index, pos)
        end

        a_block[i]=[sparse(a_ind2[Ind[i][ind]],a_ind1[Ind[i][ind]].-(r[i]+r_cliq[i][ind]),a_val[Ind[i][ind]],zeta,u_cliq[i][ind]) for ind in 1:omega_cliq[i]]
       
        a0_block[i]=[a0[r[i]+r_cliq[i][ind]+1:r[i]+r_cliq[i][ind]+u_cliq[i][ind]]' for ind in 1:omega_cliq[i]]
    end
   
    println("Modeling time:")
 
    return omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p
end

function test_Mosek_CS(a,a0,d,p,omega_cliq,s_cliq,zeta,norm_a0,opnorm_a,ak)
    b=spzeros(Float64,zeta)
    b[end]=1/opnorm_a/ak
    println("Mosek:-----------------")
    model=JuMP.Model(with_optimizer(Mosek.Optimizer, QUIET=true))

    xvec=@variable(model, [1:d])
    
    @constraint(model, a'*xvec.==b)
    
    X=Vector{Vector{Matrix{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}}}}(undef,p)
    t=1
    for j=1:p
        X[j]=Vector{Matrix{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}}}(undef,omega_cliq[j])
        for i=1:omega_cliq[j]
            X[j][i]=Matrix{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}}(undef,s_cliq[j][i],s_cliq[j][i])
            for ii=1:s_cliq[j][i]
                for jj=1:ii
                    X[j][i][ii,jj]=xvec[t]/sqrt(1+(jj<ii))
                    X[j][i][jj,ii]= X[j][i][ii,jj]
                    t+=1
                end
            end
            if s_cliq[j][i]!=1
                @constraint(model, X[j][i] in PSDCone())
            else
                @constraint(model, X[j][i].>=0)
            end
        end
    end
    
    @objective(model, Min, a0'*xvec)
    optimize!(model)
    println(termination_status(model))
    println(objective_value(model)*norm_a0)
    println("-----------------")
end

function rescale_CS_LMBM(a_ind1,a_ind2,a_val,a_len,a0,zeta)


    norm_a=zeros(Float64,zeta)
    
    for t in 1:a_len
        norm_a[a_ind2[t]]+=a_val[t]^2
    end
        
    norm_a=sqrt.(norm_a)
   
    
   
    for t in 1:a_len
        a_val[t]/=norm_a[a_ind2[t]]
    end
    
    
    opnorm_a=1#svds(sparse(a_ind1,a_ind2,a_val), nsv = 1)[1].S[1]
    
    a_val=a_val./opnorm_a

    
    norm_a0=1#norm(a0)
    a0=a0./norm_a0
    

    return a_val,a0,norm_a0,opnorm_a
end

function SmallEig_CS_LMBM(mat,s;EigAlg="Arpack")
    if s==1
        return mat[1,1], ones(Float64,1)
    elseif EigAlg=="Arpack"
       E=eigs(mat,nev = 1,which=:SR,tol=1e-2) 
       return E[1][1],E[2][:,1]
    elseif EigAlg=="ArnoldiMethod"      
       E=partialeigen(partialschur(mat, nev=1,tol=1e-2 , which=SR())[1])
       return E[1][1],E[2][:,1]
    elseif EigAlg=="Normal"
       
       E=eigen(Symmetric(mat),1:1)
       return E.values[1],E.vectors[:,1]
    elseif EigAlg=="Mix"
       try 
           E=eigs(mat,nev = 1,which=:SR) 
           return E[1][1],E[2][:,1]
       catch
           E=eigen(Symmetric(mat),1:1)
           return E.values[1],E.vectors[:,1]
       end
    else
       println("No eigenvalue algorithm!!!")
    end  
    
end

                                        
function SmallEig_block_CS_LMBM(vec,sk;EigAlg="Arpack")
    B=zeros(Float64,sk,sk)
    t=1
    for i in 1:sk, j in 1:i
        @inbounds B[i,j]=vec[t]/sqrt(1+(j<i))
        @inbounds B[j,i]= copy(B[i,j])
        t+=1
    end
    return SmallEig_CS_LMBM(B,sk,EigAlg=EigAlg)
end


function SmallEigBlocks_CS_LMBM(a0_block,a_block,y,s,omega;EigAlg="Arpack")
    
    eigval,eigvec=Float64(Inf),zeros(Float64,1)

    smalleigval=eigval
    smalleigvec=eigvec
    ind=0

 
    for t in 1:omega
        @inbounds eigval,eigvec=SmallEig_block_CS_LMBM(a0_block[t]+y'*a_block[t],s[t],EigAlg=EigAlg)
        if smalleigval>eigval
            @inbounds smalleigval=eigval
            @inbounds smalleigvec=eigvec
            @inbounds ind=t
        end
    end

    return smalleigval,smalleigvec,ind
end



function solve_CS_LMBM(omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p;EigAlg="Arpack",tol=1e-4)
                    
                    
    println("**LMBM solver:")
    function phi(nvar::Cint,xp::Ptr{Cdouble},gp::Ptr{Cdouble})
        zvar=unsafe_wrap(Array, xp, (convert(Int, nvar),))
        grad=unsafe_wrap(Array, gp, (convert(Int, nvar),))  
        
        eigval=Float64(0)
        eigvec=zeros(Float64,1)
        ind=0              
          
        grad[:]=zeros(Float64,zeta)
        grad[nvar]=1/opnorm_a/ak
        obj=zvar[nvar]/opnorm_a/ak
        for j in 1:p             
            eigval,eigvec,ind=SmallEigBlocks_CS_LMBM(a0_block[j],a_block[j],zvar,s_cliq[j],omega_cliq[j],EigAlg=EigAlg)         
            grad-=a_block[j][ind]*[@inbounds eigvec[i]*eigvec[w]*sqrt(1+(w<i)) for i in 1:s_cliq[j][ind] for w in 1:i]
            obj-=eigval
        end
        
                    
        return (convert(Cdouble,obj))
    end

                        
                        
    opt_val,~= lmbm(phi,zeros(Float64,zeta);printinfo=true,tol=tol)
                     
                        
    opt_val*=-norm_a0
    
    println()
    println("####################################")
    println("opt_val = ",opt_val)
    println("####################################")
    return opt_val
                                  
end
    


function POP_CS_LMBM(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{SparseMatrixCSC{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{SparseMatrixCSC{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::SparseMatrixCSC{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64;EigAlg="Arpack",tol=1e-3,UseEq=true)

    @time begin
    
    
    @time omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p=model_CS_LMBM(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,UseEq=UseEq)

      
    @time opt_val=solve_CS_LMBM(omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p,EigAlg=EigAlg,tol=tol)
    
    println("Total time:")
                    end
    return opt_val
end
                
  