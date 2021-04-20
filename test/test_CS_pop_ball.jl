function run_CS_POP(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64)
    
    println("Time to get information:")
    @time n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=ctpNCPOP.get_info(x,f,g,h);
    println()
    println("--------------------------------------------------")
    println()
    println("**CTP+CGAL**")
    POP_CS_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,
    maxit=Int64(1e10),tol=1e-3,use_eqcons_to_get_constant_trace=false,check_tol_each_iter=true)
    
    println()
    println("--------------------------------------------------")
    println()
    try
        println("**CS+Mosek**")
        @time POP_CS_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,
    maxit=Int64(1e10),tol=5e-3,use_eqcons_to_get_constant_trace=false,check_tol_each_iter=true,check_by_mosek=true)
    catch
        println("Mosek is out of space!!!")
    end
    println()
    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    println()
end




function test_CS_POP_ball(n::Int64,umin::Int64,umax::Int64,k::Int64;have_eqcons::Bool=true)
    @ncpolyvar x[1:n]
    f=Polynomial{false,Float64}(x[1]+0.0)
    g=Vector{Polynomial{false,Float64}}([x[1]+0.0])
    h=Vector{Polynomial{false,Float64}}([x[1]+0.0])
    
    for u in umin:5:umax
        x,f,g,h=generate_CS_POP_ball(n,u,have_eqcons=have_eqcons)
        println("Relaxation order: k=",k)
        println("====================")
        run_CS_POP(x,f,g,h,k)
    end
end
    
function generate_CS_POP_ball(n::Int64,u::Int64;have_eqcons::Bool=false)
    println("***Problem setting***")
    println("Number of variable: n=",n)
    println("====================")

    @ncpolyvar x[1:n]# variables

    function star_algebra(mom)
        if mom==1
            return mom
        else
            ind=mom.z .>0
            vars=mom.vars[ind]
            exp=mom.z[ind]
            return prod(vars[i]^exp[i] for i in length(exp):-1:1)
        end
    end

    #function to get a random quadratic polynomial of variables x(T)
    function generate_random_poly(T::UnitRange{Int64})
        v=reverse(monomials(x[T],0:2))
        v+=star_algebra.(v)
        v=v./2
        c=2*rand(Float64,length(v)).-1
        return c'*v
    end

    # ball constraints on subsets of variables
    #u=10# clique size
    p=floor(Int64,n/u) #number of cliques

    # indices of variables on each clique
    I=Vector{UnitRange{Int64}}(undef,p)
    I[1]=1:u 
    I[2:p-1]=[u*(j-1):u*j for j in 2:p-1]
    I[p]=u*(p-1):n

    # random quadratic objective function f
    vecf=[generate_random_poly(I[j]) for j in 1:p] #vector of separable polynomials on each clique
    f=sum(vecf)/100

    # ball constraints on each clique
    g=[1.0-sum(x[I[j]].^2) for j in 1:p]
    J=[j:j for j in 1:p] # assign inequality constraints



    m=length(g)
    println("Number of inequality constraints: m=",m)
    println("====================")

    if have_eqcons
        l=ceil(Int64, n/7)
    else
        l=0
    end

    r=floor(Int64,l/p)
    W=[(j-1)*r+1:j*r for j in 1:p-1]# assign equality constraints
    append!(W,[(p-1)*r+1:l])

    h=Vector{Polynomial{false,Float64}}(undef,l)


    # get a random point satisfies the inequality constraints
    randx=2*rand(Float64,n).-1

    for j in 1:p
        randx[I[j]]=randx[I[j]]./norm(randx[I[j]])
        randx[I[j]]=rand(Float64,1)[1]*randx[I[j]]
    end

    for j in 1:p
        for i in W[j]
            h[i]=generate_random_poly(I[j])
            h[i]-=h[i](x => randx) #make the random point satisfy the equality constraint h[i](randx) = 0
        end
    end

    l=length(h)
    println("Number of equality constraints: l=",l)
    println("====================")

    return x,f,g,h
end
