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


function POP_dense_SOS(x::Vector{PolyVar{false}},f::Polynomial{false},g::Vector{Polynomial{false,Float64}},h::Vector{Polynomial{false,Float64}},k::Int64;tol::Float64=1e-5,solver="Mosek")
@time begin
    n=length(x)
    l_g=length(g)
    l_h=length(h)

    

    sk=numlet(n,k)
    sk_g=[@inbounds numlet(n,k-ceil(Int64,maxdegree(g[i])/2)) for i in 1:l_g]


    sk_h=[@inbounds numlet(n,k-ceil(Int64,maxdegree(h[i])/2)) for i in 1:l_h]

    if solver=="Mosek"
        model = Model(with_optimizer(Mosek.Optimizer, QUIET=false))

        set_optimizer_attribute(model, "INTPNT_CO_TOL_REL_GAP", tol)
    elseif solver=="COSMO"
        #model = SOSModel(with_optimizer(COSMO.Optimizer, decompose = true, merge_strategy = COSMO.NoMerge))
        model = Model(with_optimizer(COSMO.Optimizer, decompose = false, merge_strategy = COSMO.NoMerge))
    else
        println("No solver!!!")
    end
        
    @variable(model, lambda)
    @objective(model, Max, lambda)

    wsos=f-lambda
    monos = reverse(monomials(x, 0:k))
    monos_star=star_algebra.(monos)
    


    G0=@variable(model, [1:sk,1:sk], PSD)
    
    
    wsos-=dot(G0,monos_star*monos')
        
    G=Vector{LinearAlgebra.Symmetric{VariableRef,Array{VariableRef,2}}}(undef,l_g)
    for i in 1:l_g
        G[i]= @variable(model,[1:sk_g[i],1:sk_g[i]], PSD)
        wsos-=dot(G[i],monos_star[1:sk_g[i]]*g[i]*monos[1:sk_g[i]]')
    end


    H=Vector{LinearAlgebra.Symmetric{VariableRef,Array{VariableRef,2}}}(undef,l_h)
    for i in 1:l_h
        H[i]= @variable(model,[1:sk_h[i],1:sk_h[i]], Symmetric)
        wsos-=dot(H[i],monos_star[1:sk_h[i]]*h[i]*monos[1:sk_h[i]]')
    end

    
    @constraint(model, wsos==0)
    optimize!(model)
    println(termination_status(model))
    opt_val=objective_value(model)
    println("opt_val=",opt_val)
    end
    return opt_val
end

