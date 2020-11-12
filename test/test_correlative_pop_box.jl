function test_correlative_POP_box(n::Int64,umin::Int64,umax::Int64,eqcons::Bool)
    k=2
  
    
    @polyvar x[1:n]
    f=Polynomial{true,Float64}(x[1]+0.0)
    g=Vector{Polynomial{true,Float64}}([x[1]+0.0])
    h=Vector{Polynomial{true,Float64}}([x[1]+0.0])
    
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=ctpPOP.get_info(x,f,g,h,sparse=true)
    
    
    for u in umin:5:umax
        @time begin
        n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=correlative_POP_box(n,u,k,eqcons=eqcons)
        println("Preparation time:")
        end
        println()
        println("--------------------------------------------------")
        println()
        POP_CS_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,EigAlg="Arpack",maxit=1e10,tol=1e-2,UseEq=false)
        println()
        println("--------------------------------------------------")
        println()
        POP_CS_LMBM(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,EigAlg="Arpack",tol=1e-3,UseEq=false)
        println()
        println("--------------------------------------------------")
        println()
        try
            @time cs_tssos_first(Vector{SparseMatrixCSC{UInt8,UInt32}}([[supp_f];supp_g;supp_h]),[[coe_f];coe_g;coe_h],n,k,[dg;dh],numeq=l,CS="MD",TS=false,CTP=false)
        catch
            println("Mosek is out of space!!!")
        end
        println()
        println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
        println()
    end
end
    
function correlative_POP_box(n::Int64,u::Int64,k::Int64;eqcons::Bool=false)
    
    println("***Problem setting***")
    println("Number of variable: n=",n)
    println("====================")

    @polyvar x[1:n]# variables

    function generate_random_poly(T::UnitRange{Int64})
        v=reverse(monomials(x[T],0:2))
        c=2*rand(Float64,length(v)).-1
        return c'*v
    end

    # unit sphere constraint
    p=floor(Int64,n/u)+1

    I=Vector{UnitRange{Int64}}(undef,p)
    I[1]=1:u
    I[2:p-1]=[u*(j-1):u*j for j in 2:p-1]
    I[p]=u*(p-1):n

    # random quadratic objective function f
    f=sum(generate_random_poly(I[j]) for j in 1:p)

    g=-(x.^2).+(1/u)

    m=length(g)
    println("Number of inequality constraints: m=",m)
    println("====================")
    
    if eqcons
        l=ceil(Int64, n/7)
    else
        l=0
    end
    r=floor(Int64,l/p)
    W=[(j-1)*r+1:j*r for j in 1:p-1]
    append!(W,[(p-1)*r+1:l])

    h=Vector{Polynomial{true,Float64}}(undef,l)

    #=randx=[2*rand(length(I[j])).-1 for j in 1:p]# create a feasible solution
    randx=[randx[j]/norm(randx[j]) for j in 1:p]=#

    randx=2*rand(n).-1
    randx=randx./sqrt(u)
   

    for j in 1:p
        for i in W[j]
            h[i]=generate_random_poly(I[j])
            h[i]-=h[i](x => randx) #make constraints feasible
        end
    end

    l=length(h)
    println("Number of equality constraints: l=",l)
    println("====================")
    println("Relaxed order: k=",k)
    println("====================")

        
    return get_info(x,f,g,h,sparse=true)
end
