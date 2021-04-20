function test_Mosek_CS(a_block,a0_block,p,omega_cliq,s_cliq,zeta)
    b=spzeros(Float64,zeta)
    b[end]=1
    println("Mosek:-----------------")
    model=JuMP.Model(with_optimizer(Mosek.Optimizer, QUIET=true))

    xvec=[@variable(model, [1:d]) for j in ]
    
    @constraint(model, a'*xvec.==b)
    
    X=Vector{Vector{Matrix{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}}}}(undef,p)
    t=1
    for j=1:p
        X[j]=Vector{Matrix{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}}}(undef,omega_cliq[j])
        for i=1:omega_cliq[j]
            X[j][i]=Matrix{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}}(undef,s_cliq[j][i],s_cliq[j][i])
            for ii=1:s_cliq[j][i]
                for jj=1:ii
                    X[j][i][ii,jj]=xvec[t]/sqrt(1+(jj<ii))
                    X[j][i][jj,ii]= X[j][i][ii,jj]
                    t+=1
                end
            end
            if s_cliq[j][i]!=1
                @constraint(model, X[j][i] in PSDCone())
            else
                @constraint(model, X[j][i].>=0)
            end
        end
    end
    
    @objective(model, Min, a0'*xvec)
    optimize!(model)
    println(termination_status(model))
    println(objective_value(model))
    println("-----------------")
end





function solve_CS_Mosek(omega_cliq::Vector{Int64},a0_block::Vector{Vector{Adjoint{Float64,Array{Float64,1}}}},a_block::Vector{Vector{SparseMatrixCSC{Float64}}},s_cliq::Vector{Vector{UInt64}},zeta::UInt64,norm_a0::Float64,opnorm_a::Float64,ak::Float64,p::Int64;tol::Float64=1e-4)
                    
                    
test_Mosek_CS(sparse(a_ind1,a_ind2,a_val),a0,d,p,omega_cliq,s_cliq,zeta)

    println("**LMBM solver:")
    function phi(nvar::Cint,xp::Ptr{Cdouble},gp::Ptr{Cdouble})
        zvar=unsafe_wrap(Array, xp, (convert(Int, nvar),))
        grad=unsafe_wrap(Array, gp, (convert(Int, nvar),))  
        
        eigval=Float64(0)
        eigvecmat=zeros(Float64,1)
        ind=0              
          
        grad[:]=zeros(Float64,zeta)
        grad[nvar]=1/opnorm_a/ak
        obj=zvar[nvar]/opnorm_a/ak
        for j in 1:p             
            eigval,eigvecmat,ind=SmallEigBlocks_dense(a0_block[j],a_block[j],zvar,s_cliq[j],omega_cliq[j])       
            grad-=a_block[j][ind]*eigvecmat
            obj-=eigval
        end
        
                    
        return convert(Cdouble,obj)
    end

                        
                        
    opt_val,~= lmbm(phi,zeros(Float64,zeta);printinfo=true,tol=tol)
                     
                        
    opt_val*=-norm_a0
    
    println()
    println("####################################")
    println("opt_val = ",opt_val)
    println("####################################")
    return opt_val
                                  
end

    


function POP_CS_LMBM(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{SparseMatrixCSC{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{SparseMatrixCSC{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::SparseMatrixCSC{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64;tol::Float64=1e-3,use_eqcons_to_get_constant_trace::Bool=true)

    @time begin
    
    
    @time omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p=model_CS_CGAL(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,use_eqcons_to_get_constant_trace=use_eqcons_to_get_constant_trace)

      
    @time opt_val=solve_CS_LMBM(omega_cliq,a0_block,a_block,s_cliq,zeta,norm_a0,opnorm_a,ak,p,tol=tol)
    
    println("Total time:")
                    end
    return opt_val
end

