function PMSV_Mosek(M)
    n=size(M,1)
    C=M'*M
    
    model=JuMP.Model(with_optimizer(Mosek.Optimizer, QUIET=true))
    
    X=@variable(model, [1:n, 1:n], PSD)
    
    @constraint(model, X.>=0)
    
    @constraint(model, sum(X[i,i] for i=1:n)==1)
    
    @objective(model, Max, dot(C,X))
    
    optimize!(model)
    
    println("Status: ",termination_status(model))
    
    val=sqrt(objective_value(model))
    
    println("val=",val)
    return val
end