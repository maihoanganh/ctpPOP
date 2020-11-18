function POP_dense_SOS(x::Vector{PolyVar{true}},f::Polynomial{true},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64;tol::Float64=1e-5,solver="Mosek")
@time begin
    n=length(x)
    l_g=length(g)
    l_h=length(h)

    

    sk=binomial(k+n,n)
    sk_g=[@inbounds binomial(k-ceil(Int64,maxdegree(g[i])/2)+n,n) for i in 1:l_g]


    s2k_h=[@inbounds binomial(2*(k-ceil(Int64,maxdegree(h[i])/2))+n,n) for i in 1:l_h]

    if solver=="Mosek"
        model = SOSModel(with_optimizer(Mosek.Optimizer, QUIET=false))

        set_optimizer_attribute(model, "INTPNT_CO_TOL_REL_GAP", tol)
    elseif solver=="COSMO"
        #model = SOSModel(with_optimizer(COSMO.Optimizer, decompose = true, merge_strategy = COSMO.NoMerge))
        model = SOSModel(with_optimizer(COSMO.Optimizer, decompose = false, merge_strategy = COSMO.NoMerge))
    else
        println("No solver!!!")
    end
        
    @variable(model, lambda)
    @objective(model, Max, lambda)

    wsos=f-lambda
    psi_monos = reverse(monomials(x, 0:2*k))
    


    sigma0=@variable(model, [1:1], SOSPoly(psi_monos[1:sk]))
    
    
    wsos-=sigma0[1]
    for i in 1:l_g
        sigma= @variable(model,[1:1], SOSPoly(psi_monos[1:sk_g[i]]))
        
        wsos-=sigma[1]*g[i]
    end


    for i in 1:l_h
        psi=@variable(model, [1:1],Poly(psi_monos[1:s2k_h[i]]))
        wsos-=psi[1]*h[i]
    end

    
    @constraint(model, wsos==0)
    optimize!(model)
    println(termination_status(model))
    opt_val=objective_value(model)
    println("opt_val=",opt_val)
    end
    return opt_val
end

