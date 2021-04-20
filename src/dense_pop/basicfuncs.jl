numlet(n,d)=sum(n^j for j in 0:d)

function get_info(x::Vector{PolyVar{false}},f::Polynomial{false,Float64},g::Vector{Polynomial{false,Float64}},h::Vector{Polynomial{false,Float64}})
    n=length(x)
    m=length(g)
    l=length(h)
    
    lmon_g=Vector{UInt64}(undef,m)
    coe_g=Vector{Vector{Float64}}(undef,m)
    
    lmon_h=Vector{UInt64}(undef,l)
    coe_h=Vector{Vector{Float64}}(undef,l)
    
    supp_g=Vector{Vector{Vector{UInt64}}}(undef,m)
    supp_h=Vector{Vector{Vector{UInt64}}}(undef,l)
    
        
    dg=Vector{Int64}(undef,m)
    dh=Vector{Int64}(undef,l)
    
    lmon_f,supp_f,coe_f=info(f,x,n;sparse=sparse)
    
    for i in 1:m
        dg[i]=maxdegree(g[i])
        lmon_g[i],supp_g[i],coe_g[i]=info(g[i],x,n;sparse=sparse)
    end
                    
    for i in 1:l
        dh[i]=maxdegree(h[i])
        lmon_h[i],supp_h[i],coe_h[i]=info(h[i],x,n;sparse=sparse)
    end
    return n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh
end

function info(f::Polynomial{false},x::Vector{PolyVar{false}},n::Int64;sparse=false)
    n=length(x)
    mon=monomials(f)
    coe=coefficients(f)
    lm=length(mon)
    supp=[UInt64[] for i in 1:lm]
    for i in 1:lm
        ind=mon[i].z .>0
        vars=mon[i].vars[ind]
        exp=mon[i].z[ind]
        for j in 1:length(vars)
            k=ncbfind(x, n, vars[j], rev=true)
            append!(supp[i], k*ones(UInt64, exp[j]))
        end
    end
    return lm, supp, coe
end
       



function get_basis(n,d)
    lb=binomial(n+d,d)
    basis=zeros(UInt8,n,lb)
    i=0
    t=1
    while i<d+1
        t+=1
        if basis[n,t-1]==i
           if i<d
              basis[1,t]=i+1
           end
           i+=1
        else
            j=findfirst(x->basis[x,t-1]!=0,1:n)
            basis[:,t]=basis[:,t-1]
            if j==1
               basis[1,t]-=1
               basis[2,t]+=1
            else
               basis[1,t]=basis[j,t]-1
               basis[j,t]=0
               basis[j+1,t]+=1
            end
        end
    end
    return basis
end


function get_ncbasis(n, d; ind=UInt64[i for i=1:n])
    basis=[UInt64[]]
    for i=1:d
        append!(basis, _get_ncbasis_deg(n, i, ind=ind))
    end
    return basis
end

function _get_ncbasis_deg(n, d; ind=UInt64[i for i=1:n])
    if d>0
        basis=Vector{UInt64}[]
        for i=1:n
            temp=_get_ncbasis_deg(n, d-1, ind=ind)
            push!.(temp, ind[i])
            append!(basis, temp)
        end
        return basis
    else
        return [UInt64[]]
    end
end

#function bfind(A::Matrix{UInt64},l::Int64,a::Vector{UInt64},n::Int64)
function ncbfind(A, l, a; rev=false)
    if l==0
        return 0
    end
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        if isequal(A[mid], a)
           return mid
        elseif isless(A[mid], a)
            if rev==false
                low=mid+1
            else
                high=mid-1
            end
        else
            if rev==false
                high=mid-1
            else
                low=mid+1
            end
        end
    end
    return 0
end


function _sym_canon(a::Vector{UInt64})
    i=1
    while i<=Int(ceil((length(a)-1)/2))
        if a[i]<a[end+1-i]
            return a
        elseif a[i]>a[end+1-i]
            return reverse(a)
        else
            i+=1
        end
    end
    return a
end
