function POP_mix_LMBM(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{SparseMatrixCSC{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{SparseMatrixCSC{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::SparseMatrixCSC{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64,t::Int64;tol::Float64=1e-3,use_eqcons_to_get_constant_trace::Bool=true)

    @time begin
    
    
    @time omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p=model_mix_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,t,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)

      
    @time opt_val=solve_CS_LMBM(omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p,tol=tol)
    
    println("Total time:")
                    end
    return opt_val
end