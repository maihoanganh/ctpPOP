function test_OPF_case39_epri(data::Dict{String,Any})

    x,f,g,h=ctpPOP.get_CTP_POP_OPF(data);

    println("Number of variable: n=",length(x))
    println("Number of inequality constraints: m=",length(g))
    println("Number of equality constraints: l=",length(h))
    println("====================")
    k=1; t=2;
    println("Relaxed order: k=",k)
    println("Sparse order: t=",t)
    println("====================")

    #f,g,h=ctpPOP.rescale_POP_OPF3(x,f,g,h,k);

    

    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=ctpPOP.get_info(x,f,g,h,sparse=true);
    
    #=norm_coe_f=1#norm(coe_f,Inf)
    
    println("norm_coe_f=",norm_coe_f)
    
    f/=norm_coe_f
    g=[g[i]./abs(coe_g[i][1]) for i in 1:m];
    #h=[h[i]./norm(coe_h[i],Inf) for i in 1:l];
    
    coe_f=coe_f./norm_coe_f
    coe_g=[coe_g[i]./abs(coe_g[i][1]) for i in 1:m];
    #coe_h=[coe_h[i]./norm(coe_h[i],Inf) for i in 1:l];=#
    

    @time begin
    opt,sol,data=cs_tssos_first(Vector{SparseMatrixCSC{UInt8,UInt32}}([[supp_f];supp_g;supp_h]),[[coe_f];coe_g;coe_h],n,k,[dg;dh],numeq=l,CS="MD",TS="block",CTP=false);
    #opt,sol,data=cs_tssos_higher!(data,TS="block");
    end

    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")


    #ctpPOP.POP_CGAL(x,f,g,h,k,maxit=1.1e6,tol=1e-2,UseEq=false)
 
    #ctpPOP.POP_TS_CGAL(x,f,g,h,k,t,EigAlg="Arpack",maxit=1e10,tol=1e-2,UseEq=false)

    #ctpPOP.POP_mix_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,t,maxit=1e10,tol=1e-2,UseEq=false)
    #ctpPOP.POP_CS_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,maxit=1e10,tol=1e-2,UseEq=false)
    
    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    
    @time NLP_opf(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f)
    
end


function run_opf_case39_epri()
    data=PolyPowerModels.parse_file("/home/hoanganh/Desktop/math-topics/ctpPOP/codes/OPFproblems/pglib_opf_case39_epri.m")
    test_OPF_case39_epri(data)

end
