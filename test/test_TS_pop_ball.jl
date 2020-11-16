function run_TS_POP(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64,t::Int64)
    
    println("Time to get information:")
    @time n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=ctpPOP.get_info(x,f,g,h,sparse=false);
    println()
    println("--------------------------------------------------")
    println()
    println("**CTP+CGAL**")
    POP_TS_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,t,
                     maxit=Int64(1e6),tol=1e-2,
                     use_eqcons_to_get_constant_trace=false,
                     check_tol_each_iter=true)
    println()
    println("--------------------------------------------------")
    println()
    println("**CTP+LMBM**")
    POP_TS_LMBM(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,t,
                     tol=1e-2,use_eqcons_to_get_constant_trace=false)
    println()
    println("--------------------------------------------------")
    println()
    try
        println("**TSSOS+Mosek**")
        @time begin
            ~,~,data=tssos_first([[f];g;h],x,k,numeq=length(h),TS="block",quotient=false);
            for j in 1:t-1
                ~,~,data=tssos_higher!(data,TS="block")
            end
        end
    catch
        println("Mosek is out of space!!!")
    end
    println()
    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    println()
end



function test_TS_POP_ball(nvarmin::Int64,nvarmax::Int64,k::Int64,t::Int64;have_eqcons::Bool=true)
    @polyvar x[1:1]
    f=Polynomial{true,Float64}(x[1]+0.0)
    g=Vector{Polynomial{true,Float64}}([])
    h=Vector{Polynomial{true,Float64}}([])
    
    for n in nvarmin:10:nvarmax
        x,f,g,h=generate_TS_POP_ball(n,have_eqcons=have_eqcons)
        println("Relaxation order: k=",k)
        println("====================")
        println("Sparse order: t=",t)
        println("====================")
        run_TS_POP(x,f,g,h,k,t)
    end
end
    
function generate_TS_POP_ball(n::Int64;have_eqcons::Bool=false)
    println("***Problem setting***")
    println("Number of variable: n=",n)
    println("====================")

    @polyvar x[1:n]# variables

    function generate_random_poly(v)
        c=2*rand(Float64,length(v)).-1
        return c'*v
    end
    # random quadratic objective function f
    v=reverse(monomials(x,2))
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
