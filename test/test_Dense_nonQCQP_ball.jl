function run_dense_nonQCQP_ball(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64)
    
    println("Time to get information:")
    @time n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=ctpPOP.get_info(x,f,g,h,sparse=false);
    println()
    println("--------------------------------------------------")
    println()
    println("**CTP+CGAL**")
    POP_dense_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,
                     maxit=Int64(1e6),tol=1e-3,
                     use_eqcons_to_get_constant_trace=false,
                     check_tol_each_iter=true)
    println()
    println("--------------------------------------------------")
    println()
    #println("**Ipopt**")
    #POP_NLP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f)
    println("**SOS+Mosek**")
    try
        POP_dense_SOS(x,f,g,h,k,tol=1e-2)
    catch
        println("Mosek is out of space!!!")
    end

    println()
    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    println()
end

function test_dense_nonQCQP_ball()
    
    n=15
    d=3
    x,f,g,h=generate_dense_nonQCQP_ball(n,d,have_eqcons=true)
    for k in [2;3]
        println("Relaxed order: k=",k)
        println("====================")
        run_dense_nonQCQP_ball(x,f,g,h,k)
    end
    
    n=10
    d=4
    x,f,g,h=generate_dense_nonQCQP_ball(n,d,have_eqcons=true)
    for k in [2;3]
        println("Relaxed order: k=",k)
        println("====================")
        run_dense_nonQCQP_ball(x,f,g,h,k)
    end
    
#     n=5
#     d=5
#     x,f,g,h=generate_dense_nonQCQP_ball(n,d,have_eqcons=true)
#     for k in [3;4]
#         println("Relaxed order: k=",k)
#         println("====================")
#         run_dense_nonQCQP_ball(x,f,g,h,k)
#     end
end
    

function generate_dense_nonQCQP_ball(n::Int64,d;have_eqcons::Bool=false)
    println("***Problem setting***")
    println("Number of variable: n=",n)
    println("====================")

    @polyvar x[1:n]# variables

    function generate_random_poly(v)
        c=2*rand(Float64,length(v)).-1
        return c'*v
    end
    # random quadratic objective function f
    v=reverse(monomials(x,0:d))
    f=generate_random_poly(v)


    # unit sphere constraint
    g=[1.0-sum(x.^2)] #type of coefficients of each polynomial must be float

    m=length(g)
    println("Number of inequality constraints: m=",m)
    println("====================")

    if have_eqcons
        l=ceil(Int64, n/4)
    else
        l=0
    end
    h=Vector{Polynomial{true,Float64}}(undef,l)
    randx=2*rand(n).-1# create a feasible solution
    randx=rand(1)[1]*randx./sqrt(sum(randx.^2))


    for j in 1:l
        h[j]=generate_random_poly(v[2:end])
        h[j]-=h[j](x => randx) #make constraints feasible
    end

    l=length(h)
    println("Number of equality constraints: l=",l)
    println("====================")


    return x,f,g,h
end