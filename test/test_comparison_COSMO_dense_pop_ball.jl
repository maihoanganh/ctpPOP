function run_comparison_dense_POP(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64)
    
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
    println("**Ipopt**")
    POP_NLP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f)
    println()
    println("--------------------------------------------------")
    println()
    println("**SOS+COSMO**")

    POP_dense_SOS(x,f,g,h,k,tol=1e-3,solver="COSMO")

    println()
    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    println()
end

function test_comparison_dense_POP_ball(nvarmin::Int64,nvarmax::Int64,k::Int64;have_eqcons::Bool=true)
    @polyvar x[1:1]
    f=Polynomial{true,Float64}(x[1]+0.0)
    g=Vector{Polynomial{true,Float64}}([])
    h=Vector{Polynomial{true,Float64}}([])
    
    for n in nvarmin:10:nvarmax
        x,f,g,h=generate_dense_POP_ball(n,have_eqcons=have_eqcons)
        println("Relaxed order: k=",k)
        println("====================")
        run_comparison_dense_POP(x,f,g,h,k)
    end
end
    
