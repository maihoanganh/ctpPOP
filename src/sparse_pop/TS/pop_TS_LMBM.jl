



function POP_TS_LMBM(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Matrix{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Matrix{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::Matrix{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64,t::Int64;tol::Float64=1e-4,use_eqcons_to_get_constant_trace::Bool=true)

    @time begin
    
    
    @time omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak=model_POP_TS(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,t,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)
    
 
    @time opt_val=solve_POP_dense_LMBM(omega,a0_block,a_block,s,zeta,norm_a0,opnorm_a,ak,tol=tol)
    
    println("Total time:")
                    end
    return opt_val
end
                    
