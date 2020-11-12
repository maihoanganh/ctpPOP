function make_vec_PMSV_CGAL(mat,n,s)
    vec=Vector{Float64}(undef,s)
    t=1
    for i in 1:n, j in i:n
        vec[t]=mat[i,j]*sqrt(1+(j>i))
        t+=1
    end
    return vec
end
            
function make_mat_PMSV_CGAL(vec,n)
    mat=Matrix{Float64}(undef,n,n)
    t=1
    val=0.0
    for i in 1:n, j in i:n
        val=vec[t]/sqrt(1+(j>i))
        mat[i,j]=val
        mat[j,i]=val
        t+=1
    end
    return mat
end  

function largest_eig_PMSV_CGAL(mat)
    E=eigs(mat,nev = 1,which=:LR, tol=1e-2)
    return E[1][1],E[2][:,1]
end

function update_PMSV_CGAL(x,y,val,feas,t,Gamma,c,n;EigAlg="Arpack")
    beta=sqrt(t+1)
    eta=2/(t+1)
    
    w=x+y./beta
    
    
    eigval,eigvec=largest_eig_PMSV_CGAL(make_mat_PSMV_CGAL(c+y+beta*(x-w),n))

    h=[eigvec[i]*eigvec[j]*sqrt(1+(j>i)) for i in 1:n for j in i:n]
   
    x=eta*x+(1-eta)*h

    val=dot(c,x)
                
    ind_xmin=findall(z->z<0,x)
              
                            
    feas=norm(x[ind_xmin])
                
    
    
                
    y=theta*y
    y[ind_xmin]-=beta*x[ind_xmin]
         
    return x,y,val,feas,gap
end


function PMSV_CGAL(M::Matrix{Float64};EigAlg="Arpack",maxit=1e5, tol=1e-4)
    
    n=size(M,1)
    C=M'*M
    
                
    s=Int64(n*(n+1)/2)
    
    c=make_vec_PSMV_CGAL(C,n,s)
                
    norm_c=norm(c)
    c=c./norm_c
                
    x=zeros(Float64,s)            
    y=zeros(Float64,s)
                            

    val,feas,gap=0.0,0.0,0.0
    Gamma=0.0
                
    i=1
    for t in 0:maxit
        
        x,y,val,feas,gap=update_PMSV_CGAL(x,y,val,feas,t,Gamma,c,n,EigAlg=EigAlg)
        if feas < tol# && gap<tol
            println("iter=",t,"   val=",sqrt(val*norm_c),"   feas=",feas, "   gap=",gap)
            println("feas < tol!!!")
            break
        end
        if t==i || t==maxit
            println("iter=",t,"   val=",sqrt(val*norm_c),"   feas=",feas, "   gap=",gap)
            i*=2
        end
                    
    end            
   
                
    return sqrt(val*norm_c)
end