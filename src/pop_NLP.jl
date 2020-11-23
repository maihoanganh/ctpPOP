function POP_NLP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f)#;set_initial_point_opf=false)
    @time begin
    model = Model(with_optimizer(Ipopt.Optimizer))
    @variable(model, x[1:n])

    #=if set_initial_point_opf
        set_start_value.(x[1:l], ones(Int64,l))
        end=#

    function get_func(x,lmon_f,supp_f,coe_f)
        f=0
        for j in 1:lmon_f
            ind=findall(r->r>0,supp_f[:,j])
            lind=length(ind)
            if lind==0
                f+=coe_f[j]
            elseif lind==1
                f+=coe_f[j]*x[ind[1]]^supp_f[ind[1],j]
            else
                f+=coe_f[j]*x[ind[1]]*x[ind[2]]
            end
        end
        return f
    end

    @objective(model, Min, get_func(x,lmon_f,supp_f,coe_f))
    for i in 1:m
        @constraint(model, get_func(x,lmon_g[i],supp_g[i],coe_g[i])>=0)
    end
    for i in 1:l
        @constraint(model, get_func(x,lmon_h[i],supp_h[i],coe_h[i])==0)
    end
    optimize!(model)
    println(termination_status(model))
    opt_val=objective_value(model)
    println("opt_val=",opt_val)   
                    end
    return opt_val
end
