function solve_POP_dense_LMBM(omega::Int64,a0_block::Vector{Adjoint{Float64,Array{Float64,1}}},a_block::Vector{SparseMatrixCSC{Float64}},s::Vector{UInt64},zeta::UInt64,norm_a0::Float64,opnorm_a::Float64,ak::Float64;tol::Float64=1e-4)
                
    function phi(nvar::Cint,xp::Ptr{Cdouble},gp::Ptr{Cdouble})
        zvar=unsafe_wrap(Array, xp, (convert(Int, nvar),))
        grad=unsafe_wrap(Array, gp, (convert(Int, nvar),))  
        eigval,eigvecmat,ind=SmallEigBlocks_dense(a0_block,a_block,zvar,s,omega)
        grad[:]=-a_block[ind]*eigvecmat                
        grad[nvar]+=1/opnorm_a/ak
        return(convert(Cdouble,-eigval+zvar[nvar]/opnorm_a/ak))
    end

                        
                        
    opt_val,~= lmbm(phi,zeros(Float64,zeta);printinfo=true,tol=tol)
                     
                        
    opt_val*=-norm_a0
    println()
    println("####################################")
    println("opt_val = ",opt_val)
    println("####################################")
    println("Solving time:")
    return opt_val       
                
end
    


function POP_dense_LMBM(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Matrix{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Matrix{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::Matrix{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64;tol::Float64=1e-4,use_eqcons_to_get_constant_trace::Bool=true)

    @time begin
    
    
    @time omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak=model_POP_dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)
    
 
    @time opt_val=solve_POP_dense_LMBM(omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak,tol=tol)

       println("Total time:")
                    end
    return opt_val
end
                    
