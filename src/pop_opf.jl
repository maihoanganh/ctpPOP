function pop_opf_classic(data::Dict{String, Any}; normalize = true)

    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    model = PolyModel()   
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    @assert isempty(ref[:dcline])
    
    
    

    vr = Dict()
    vi = Dict()
    mvr = Dict()
    mvi = Dict()

    for key in keys(ref[:bus])
        # real part voltage variables
        mvr[key] = PolyPowerModels.new_polyvar("mvr"*string(key))
        
        vr[key]=mvr[key]*ref[:bus][key]["vmax"]
       

        # imaginary part voltage variables
        mvi[key] = PolyPowerModels.new_polyvar("mvi"*string(key))
        vi[key]=mvi[key]*ref[:bus][key]["vmax"]
        # voltage magnitude constraints
        
        
        
        add_constraint!(model, ref[:bus][key]["vmin"]^2/ref[:bus][key]["vmax"]^2, LT, mvr[key]^2 + mvi[key]^2; normalize = false )
        add_constraint!(model, mvr[key]^2 + mvi[key]^2, LT, 1; normalize = false )

    end

    p = Dict()
    q = Dict()

    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        vr_fr = vr[branch["f_bus"]]
        vr_to = vr[branch["t_bus"]]
        vi_fr = vi[branch["f_bus"]]
        vi_to = vi[branch["t_bus"]]

        # Line Flow
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]

        p[f_idx] = (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to)+ (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)

        q[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2)  - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)

        p[t_idx] = (g+g_to)*(vr_to^2 + vi_to^2)  + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to)  + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))

        q[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
    end



    pg = Dict()
    qg = Dict()
    mpg = Dict()
    mqg = Dict()

    for (i,bus) in ref[:bus]

        # active/reactive power

        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        
        for gen_id in ref[:bus_gens][i]

            if ref[:gen][gen_id]["pmin"]<ref[:gen][gen_id]["pmax"]   
                mpg[gen_id] = PolyPowerModels.new_polyvar("mpg"*string(gen_id))
                pg[gen_id] = mpg[gen_id]*(ref[:gen][gen_id]["pmax"]-ref[:gen][gen_id]["pmin"])/2+(ref[:gen][gen_id]["pmax"]+ref[:gen][gen_id]["pmin"])/2
                add_constraint!( model, mpg[gen_id]^2, LT, 1.0; normalize = false )
            else
                pg[gen_id] = 0.0
            end
            if ref[:gen][gen_id]["qmin"]<ref[:gen][gen_id]["qmax"]
                mqg[gen_id] = PolyPowerModels.new_polyvar("mqg"*string(gen_id))
                qg[gen_id] = mqg[gen_id]*(ref[:gen][gen_id]["qmax"]-ref[:gen][gen_id]["qmin"])/2+(ref[:gen][gen_id]["qmax"]+ref[:gen][gen_id]["qmin"])/2
                add_constraint!( model, mqg[gen_id]^2, LT, 1.0; normalize = false )
            else
                qg[gen_id] = 0.0
            end
        end
        add_constraint!( model,
                        PolyPowerModels.fl_sum(p[a] for a in ref[:bus_arcs][i]), EQ, 
                        PolyPowerModels.fl_sum(pg[g] for g in ref[:bus_gens][i]) -
                        PolyPowerModels.fl_sum(load["pd"] for load in bus_loads) -
                        PolyPowerModels.fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = false
                       )  
        add_constraint!( model,
                        PolyPowerModels.fl_sum(q[a] for a in ref[:bus_arcs][i]), EQ,
                        PolyPowerModels.fl_sum(qg[g] for g in ref[:bus_gens][i]) -
                        PolyPowerModels.fl_sum(load["qd"] for load in bus_loads) +
                        PolyPowerModels.fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = false
                       )
        
    end
    
    

    # objective
    set_objective!( model, MIN, sum(gen["cost"][1]*pg[i]^2+gen["cost"][2]*pg[i]+gen["cost"][3] for (i,gen) in ref[:gen]))

    var = Dict(:mvr => mvr, :mvi => mvi, :mpg => mpg, :mqg => mqg)

    return PolyPowerModel(model, data, ref, var, Dict())
end



function get_CTP_POP_OPF(data::Dict{String, Any})
    pm = pop_opf_classic(data,normalize = true)
    #pm = pop_opf_deg2(data, normalize = true)
    
    f=objective_function(pm)
    g=Vector{Polynomial{true,Float64}}([])
    h=Vector{Polynomial{true,Float64}}([])
    
    cons=constraints(pm)
    for j in 1:length(cons)
        if sense(cons[j])==PolyPowerModels.EQ_sense()
            push!(h,constraint_function(cons[j]))
        else
            push!(g,-constraint_function(cons[j]))
        end
    end
    x = sort!(union(variables(f),variables.([g;h])...), rev = true)
    return x,f,g,h
end


