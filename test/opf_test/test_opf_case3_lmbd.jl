



function test_OPF_case3_lmbd(data::Dict{String,Any})

    x,f,g,h=ctpPOP.get_CTP_POP_OPF(data);

    println("f=",f)
println("g=",g)
println("h=",h)

    println("Number of variable: n=",length(x))
    println("Number of inequality constraints: m=",length(g))
    println("Number of equality constraints: l=",length(h))
    println("====================")
    k=2; t=2;
    println("Relaxed order: k=",k)
    println("Sparse order: t=",t)
    println("====================")

    #f,g,h=ctpPOP.rescale_POP_OPF3(x,f,g,h,k);



    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=ctpPOP.get_info(x,f,g,h,sparse=false);
    g=[g[i]./abs(coe_g[i][1]) for i in 1:m];
    #h=[h[i]./abs(coe_h[i][end]) for i in 1:l];
    
    coe_g=[coe_g[i]./abs(coe_g[i][1]) for i in 1:m];
    #coe_h=[coe_h[i]./abs(coe_h[i][end]) for i in 1:l];


    @time begin
    opt,sol,data=cs_tssos_first(Vector{SparseMatrixCSC{UInt8,UInt32}}([[supp_f];supp_g;supp_h]),[[coe_f];coe_g;coe_h],n,k,[dg;dh],numeq=l,CS="MD",TS="block",CTP=false);
    #opt,sol,data=cs_tssos_higher!(data,TS="block");
    end

    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")


    ctpPOP.POP_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,maxit=Int64(1e10),tol=1e-3)
 
    #ctpPOP.POP_TS_CGAL(x,f,g,h,k,t,EigAlg="Arpack",maxit=1e10,tol=1e-2,UseEq=false)

    #ctpPOP.POP_CS_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,EigAlg="Arpack",maxit=1e10,tol=1e-2,UseEq=true)
    #ctpPOP.POP_mix_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,t,EigAlg="Arpack",maxit=1e10,tol=1e-2,UseEq=false)
    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    
    @time NLP_opf(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f)
    
    
end


function run_opf_case3_lmbd()
    data=Vector{Dict{String,Any}}(undef,4)
    data[1] = PolyPowerModels.parse_file("/home/hoanganh/Desktop/math-topics/ctpPOP/codes/OPFproblems/pglib_opf_case3_lmbd.m")
    #data[2] = PolyPowerModels.parse_file("/home/hoanganh/Desktop/math-topics/ctpPOP/codes/OPFproblems/pglib_opf_case5_pjm.m")
    #data[3] = PolyPowerModels.parse_file("/home/hoanganh/Desktop/math-topics/ctpPOP/codes/OPFproblems/pglib_opf_case118_ieee__sad.m")
    #data[4] = PolyPowerModels.parse_file("/home/hoanganh/Desktop/math-topics/ctpPOP/codes/OPFproblems/pglib_opf_case2312_goc.m")

    for j in 1:1
        #try 
        test_OPF_case3_lmbd(data[j])
        #=catch
            println("Error!!!!!")
        end=#
        println("====================")
    end
end
