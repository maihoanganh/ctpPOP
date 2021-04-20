function run_dense_POP(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64)
    
    println("Time to get information:")
    @time n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=ctpNCPOP.get_info(x,f,g,h);
    println()
    println("--------------------------------------------------")
    println()
    println("**CTP+CGAL**")
    opt_val1=POP_dense_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,
                                 maxit=Int64(1e6),tol=1e-4,
                                 use_eqcons_to_get_constant_trace=false,
                                 check_tol_each_iter=true)
    
    println()
    println("--------------------------------------------------")
    println()
    println("**SOS+Mosek**")
    try
        nctssos_first(Vector{Vector{Vector{UInt16}}}([[supp_f];supp_g;supp_h]),[[coe_f];coe_g;coe_h],
    n,k,[dg;dh],numeq=l,reducebasis=false,obj="eigen",merge=false,TS=false,QUIET=false);
    catch
        println("Mosek is out of space!!!")
    end
    println()
    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    println()
end

function test_dense_POP_ball(nvarmin::Int64,nvarmax::Int64,k::Int64;have_eqcons::Bool=true)
    @ncpolyvar x[1:1]
    f=Polynomial{false,Float64}(x[1]+0.0)
    g=Vector{Polynomial{false,Float64}}([])
    h=Vector{Polynomial{false,Float64}}([])
    
    for n in nvarmin:10:nvarmax
        x,f,g,h=generate_dense_POP_ball(n,have_eqcons=have_eqcons)
        println("Relaxed order: k=",k)
        println("====================")
        run_dense_POP(x,f,g,h,k)
    end
end
    
function generate_dense_POP_ball(n::Int64;have_eqcons::Bool=false)
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
    function generate_random_poly(v)
        c=2*rand(Float64,length(v)).-1
        return c'*v
    end
    # random quadratic objective function f
    v=reverse(monomials(x,0:2))
    v+=star_algebra.(v)
    v=v./2

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
    h=Vector{Polynomial{false,Float64}}(undef,l)
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
