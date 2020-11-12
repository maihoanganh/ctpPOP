function make_vec(mat,n,s)
    vec=Vector{Float64}(undef,s)
    t=1
    for i in 1:n, j in i:n
        vec[t]=mat[i,j]
        t+=1
    end
    return vec
end
            
function make_mat(vec,n)
    mat=Matrix{Float64}(undef,n,n)
    t=1
    val=0.0
    for i in 1:n, j in i:n
        val=vec[t]
        mat[i,j]=val
        mat[j,i]=val
        t+=1
    end
    return mat
end  

function largest_eig(mat)
    E=eigs(mat,nev = 1,which=:LR)
    return E[1][1],E[2][:,1]
end


function psv(C::Matrix{Float64})
    n=size(C,1)
    s=Int64(n*(n+1)/2)
    
    c=make_vec(C,n,s)
    
    function fun(nvar::Cint,xp::Ptr{Cdouble},gp::Ptr{Cdouble})
      x=unsafe_wrap(Array, xp, (convert(Int, nvar),))
      g=unsafe_wrap(Array, gp, (convert(Int, nvar),))
        
      val,vec=largest_eig(make_mat(x+c,n))

      g[:]=[vec[i]*vec[j] for i in 1:n for j in i:n]
      return(convert(Cdouble,val))
    end

    x0=zeros(Float64,s)
                
    return lmbmb(fun,x0, zeros(Float64,s), 100*ones(Float64,s);printinfo=true)
end