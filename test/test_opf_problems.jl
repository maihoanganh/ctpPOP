

function get_poly_OPF(data::Dict{String,Any})
    
    x,f,g,h=get_CTP_POP_OPF(data);

    println("Number of variable: n=",length(x))
    println("Number of inequality constraints: m=",length(g))
    println("Number of equality constraints: l=",length(h))
    println("====================")
    k=2; t=2;
    println("Relaxed order: k=",k)
    println("Sparse order: t=",t)
    println("====================")

 

    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=ctpPOP.get_info(x,f,g,h,sparse=true);
    
    return k,t,n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh
end


function run_OPF(data::Dict{String,Any})


   k,t,n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_poly_OPF(data)

    #@time POP_NLP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,set_initial_point_opf=true)
    
    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    println("Shor's relaxation:")
    @time begin
    opt,sol,data=cs_tssos_first(Vector{SparseMatrixCSC{UInt8,UInt32}}([[supp_f];supp_g;supp_h]),[[coe_f];coe_g;coe_h],n,1,[dg;dh],numeq=l,CS="MD",TS="block");
    for j in 1:t-1
        opt,sol,data=cs_tssos_higher!(data,TS="block");
    end
    end

    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
   

    POP_mix_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,t,
                 maxit=Int64(1e10),tol=1e-2,use_eqcons_to_get_constant_trace=false,check_tol_each_iter=false)
    
    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")

    @time begin
    opt,sol,data=cs_tssos_first(Vector{SparseMatrixCSC{UInt8,UInt32}}([[supp_f];supp_g;supp_h]),[[coe_f];coe_g;coe_h],n,k,[dg;dh],numeq=l,CS="MD",TS="block");
    for j in 1:t-1
        opt,sol,data=cs_tssos_higher!(data,TS="block");
    end
    end

    
end
