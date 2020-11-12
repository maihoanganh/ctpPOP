function get_blocks(r,n,m,l,supp_U,lu0,lu,lw,lg,lh,supp_g,supp_h,coe_g,coe_h,basis_sigma0,basis_sigma,basis_psi)
    
    supp_U=sortslices(supp_U,dims=2)
    supp_U=unique(supp_U,dims=2)
    lsupp_U=size(supp_U,2)

    block_sigma0=Vector{Vector{Int64}}(undef,1)
    block_sigma=Vector{Vector{Vector{Int64}}}(undef,m)
    block_psi=Vector{Vector{Vector{Int64}}}(undef,l)

    lblock_sigma0=Vector{Int64}(undef,1)
    lblock_sigma=Vector{Int64}(undef,m)
    lblock_psi=Vector{Int64}(undef,l)

    lt_block_sigma0=Vector{Int64}(undef,1)
    lt_block_sigma=Vector{Vector{Int64}}(undef,m)
    lt_block_psi=Vector{Vector{Int64}}(undef,l)

    old_block_sigma0=Vector{Vector{Int64}}(undef,1)
    old_block_sigma=Vector{Vector{Vector{Int64}}}(undef,m)
    old_block_psi=Vector{Vector{Vector{Int64}}}(undef,l)

    iter=1
    while iter<=r
        graph_G0=SimpleGraph(lu0)
        for p=1:lu0
            for q=p:lu0
                @inbounds nota=basis_sigma0[:,p]+basis_sigma0[:,q]
                if bfind(supp_U,lsupp_U,nota,n)!=0
                    add_edge!(graph_G0,p,q)
                end
            end
        end
        block_sigma0=connected_components(graph_G0)
        lblock_sigma0=length(block_sigma0)
        for j=1:lblock_sigma0
            if j==1
                lt_block_sigma0=Int64(length(block_sigma0[j]))
            else
                lt_block_sigma0=[lt_block_sigma0;Int64(length(block_sigma0[j]))]
            end
        end

        
        for i=1:m
            graph_G=SimpleGraph(lu[i])
            for p=1:lu[i]
                for q=p:lu[i]
                    y=1
                    while y<=lg[i]
                        @inbounds nota=basis_sigma[i][:,p]+basis_sigma[i][:,q]+supp_g[i][:,y]
                        if bfind(supp_U,lsupp_U,nota,n)!=0
                            break
                        else
                            y+=1
                        end
                    end
                    if y<=lg[i]
                       add_edge!(graph_G,p,q)
                    end
                end
            end
            block_sigma[i]=connected_components(graph_G)
            lblock_sigma[i]=length(block_sigma[i])
            lt_block_sigma[i]=Vector{Int64}(undef,lblock_sigma[i])
            for j=1:lblock_sigma[i]
                lt_block_sigma[i][j]=length(block_sigma[i][j])
            end
        end


       
        for i=1:l
            graph_H=SimpleGraph(lw[i])
            for p=1:lw[i]
                for q=p:lw[i]
                    y=1
                    while y<=lh[i]
                        @inbounds nota=basis_psi[i][:,p]+basis_psi[i][:,q]+supp_h[i][:,y]
                        if bfind(supp_U,lsupp_U,nota,n)!=0
                            break
                        else
                            y+=1
                        end
                    end
                    if y<=lh[i]
                       add_edge!(graph_H,p,q)
                    end
                end
            end
            block_psi[i]=connected_components(graph_H)
            lblock_psi[i]=length(block_psi[i])
            lt_block_psi[i]=Vector{Int64}(undef,lblock_psi[i])
            for j=1:lblock_psi[i]
                lt_block_psi[i][j]=Int64(length(block_psi[i][j]))
            end
        end

        supp_U=[]

        for j=1:lblock_sigma0
            for p=1:lt_block_sigma0[j]
                for q=p:lt_block_sigma0[j]
                    @inbounds nota=basis_sigma0[:,block_sigma0[j][p]]+basis_sigma0[:,block_sigma0[j][q]]
                    if supp_U==[]
                        supp_U=nota
                    else
                        supp_U=[supp_U nota]
                    end
                end
            end
        end
        """
        for i=1:m
            for j=1:lblock_sigma[i]
                for p=1:lt_block_sigma[i][j]
                    for q=p:lt_block_sigma[i][j]
                        for z=1:lg[i]
                            nota=basis_sigma[i][:,block_sigma[i][j][p]]+basis_sigma[i][:,block_sigma[i][j][q]]+supp_g[i][:,z]
                            supp_U=[supp_U nota]
                        end
                    end
                end
            end
        end
        for i=1:l
            for j=1:lblock_psi[i]
                for p=1:lt_block_psi[i][j]
                    for q=p:lt_block_psi[i][j]
                        for z=1:lh[i]
                            nota=basis_psi[i][:,block_psi[i][j][p]]+basis_psi[i][:,block_psi[i][j][q]]+supp_h[i][:,z]
                            supp_U=[supp_U nota]
                        end
                    end
                end
            end
        end
        """
        supp_U=sortslices(supp_U,dims=2)
        supp_U=unique(supp_U,dims=2)
        lsupp_U=size(supp_U,2)


        if iter==1
            old_block_sigma0=block_sigma0
            old_block_sigma=block_sigma
            old_block_psi=block_psi
        else
            if old_block_sigma0==block_sigma0 && old_block_sigma==block_sigma && old_block_psi==block_psi
                println("The block-closure operation is stable at sparse order ",iter-1,"!")
                break
            else
                old_block_sigma0=block_sigma0
                old_block_sigma=block_sigma
                old_block_psi=block_psi
            end
        end
        iter+=1
    end

    return supp_U,lsupp_U,block_sigma0,block_sigma,block_psi,lblock_sigma0,lblock_sigma,lblock_psi,lt_block_sigma0,lt_block_sigma,lt_block_psi

end

####################################################
function extract_optimizers_term(Gr,lu0,basis_sigma0,n,n_sub,m,l,opt_val,f,g,h,x,I_sub)
    #extraction of optimizers
    V=nullspace(Gr,atol=1e-4)
    rk=size(V,2)

    println("Dimension of the null space of Gram matrix = ", rk)
    sol=Vector{Float64}[]
    if rk >0
        # extraction of Henrion and Lasserre
        V=Matrix(V')
        V= rref_with_pivots!(V,1e-3);
        U=V[1]

        U=Matrix(U')

        w=basis_sigma0[:,V[2]];
        N=Vector{Array{Float64}}(undef,n_sub)
        flag=1
        Idmat=zeros(n,n_sub)
        Idmat[I_sub,:]=Matrix{Int64}(LinearAlgebra.I,n_sub,n_sub)
        for i in 1:n_sub
            kk=Int64[]
            for j=1:size(w)[2]
                xwj=w[:,j]+Idmat[:,i]
                for t=1:lu0
                    if xwj==basis_sigma0[:,t]
                        if kk==Int64[]
                            kk=t
                        else
                            kk=[kk;t]
                        end
                        break
                    end
                end
            end
            
            if rk!=length(kk)
                flag=0
                break
            else
                N[i]=U[kk,:]
            end
        end

        if flag==1

            # Create random convex combination
            rands = rand(n_sub,1);rands = rands/sum(rands)
            M = zeros(length(V[2]),rk)
            for i=1:n_sub
                M+=rands[i]*N[i]
            end

            F= schur(M);
            L=F.Z
            # Extract solution
            for i=1:rk
                atom=Float64[]
                for j = 1:n_sub
                    coordinatej=L[:,i]'*N[j]*L[:,i]
                    atom=push!(atom,coordinatej[1])
                end
                println("------------------------------------")
                println("atom ",i," = ",atom)
                flag=1
                if f!=0
                    check=polynomial(f)(x[I_sub] => atom)-opt_val
                else
                    check=0
                end
                println("check lower bound  = ",check)
                if abs(check)>1e-2
                    flag=0
                end

                for i=1:m
                    check=polynomial(g[i])(x[I_sub] => atom)
                    println("check inequality ",i," = ",check)
                    if check<-1e-3
                        flag=0
                    end
                end


                for i=1:l
                    check=polynomial(h[i])(x[I_sub] => atom)
                    println("check equality ",i," = ",check)
                    if abs(check)>1e-3
                        flag=0
                    end
                end

                if flag ==1
                    sol=atom
                    println("####################################")
                    println("Solution = ",atom)
                    println("####################################")
                end

            end
        end
    end

    return sol

end

#######################################################

function Pu_block_hierarchy(x,f,g,h,k,r)

    n=length(x) # Number of variables
    m=length(g) # Number of constraints
    l=length(h) # Number of constraints

    u = Vector{Int64}(undef,m)
    w = Vector{Int64}(undef,l)

    lf,supp_f,coe_f=info(f,x,n)
   
    A=supp_f

    supp_g=Vector{Array{Int64,2}}(undef,m)
    coe_g=Vector{Vector{Float64}}(undef,m)
    lg=Vector{Int64}(undef,m)

    supp_h=Vector{Array{Int64,2}}(undef,l)
    coe_h=Vector{Vector{Float64}}(undef,l)
    lh=Vector{Int64}(undef,l)

    lu0=binomial(k+n,n)
    lu=Vector{Int64}(undef,m)
    lw=Vector{Int64}(undef,l)

    basis_sigma0=get_basis(n,k)
    basis_sigma=Vector{Array{Int64,2}}(undef,m)
    basis_psi=Vector{Array{Int64,2}}(undef,l)

    
    #A=[A basis_sigma0.*2]
    
    for i=1:m
        lg[i],supp_g[i],coe_g[i]=info(g[i],x,n)
        @inbounds A=[A supp_g[i]]

        @inbounds u[i]=ceil(Int64,maxdegree(g[i])/2)
        @inbounds lu[i]=binomial(k-u[i]+n,n)
        @inbounds basis_sigma[i]=basis_sigma0[:,1:lu[i]]
    end


    for i=1:l
        lh[i],supp_h[i],coe_h[i]=info(h[i],x,n)
        @inbounds A=[A supp_h[i]]

        @inbounds w[i]=ceil(Int64,maxdegree(h[i])/2)
        @inbounds lw[i]=binomial(k-w[i]+n,n)
        @inbounds basis_psi[i]=basis_sigma0[:,1:lw[i]]
    end

    supp_U=[A zeros(n,1)]

    supp_U,lsupp_U,block_sigma0,block_sigma,block_psi,lblock_sigma0,lblock_sigma,lblock_psi,lt_block_sigma0,lt_block_sigma,lt_block_psi=get_blocks(r,n,m,l,supp_U,lu0,lu,lw,lg,lh,supp_g,supp_h,coe_g,coe_h,basis_sigma0,basis_sigma,basis_psi)

    println("#block_sigma0 = ",[length(block_sigma0[i]) for i in 1:lblock_sigma0])
    println("-----------------------------")
    println("#block_sigma = ",[[length(block_sigma[i][j]) for j in 1:lblock_sigma[i]] for i in 1:m])
    println("-----------------------------")
    println("#block_psi = ",[[length(block_psi[i][j]) for j in 1:lblock_psi[i]] for i in 1:l])

    model=JuMP.Model(with_optimizer(Mosek.Optimizer, QUIET=true))
    cons=[AffExpr(0) for i=1:lsupp_U]

    G0=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, lblock_sigma0)
    G=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m)
    H=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, l)

    for j=1:lblock_sigma0
        if lt_block_sigma0[j]==1
            @inbounds G0[j]=@variable(model, lower_bound=0)
            @inbounds nota=Int64(2)*basis_sigma0[:,block_sigma0[j]]
            Locb=bfind(supp_U,lsupp_U,nota,n)
            @inbounds add_to_expression!(cons[Locb],G0[j])
        else
            @inbounds G0[j]=@variable(model, [1:lt_block_sigma0[j], 1:lt_block_sigma0[j]],PSD)
            for p=1:lt_block_sigma0[j]
                for q=p:lt_block_sigma0[j]
                    @inbounds nota=basis_sigma0[:,block_sigma0[j][p]]+basis_sigma0[:,block_sigma0[j][q]]
                    Locb=bfind(supp_U,lsupp_U,nota,n)
                    if p==q
                        @inbounds add_to_expression!(cons[Locb],G0[j][p,q])
                    else
                        @inbounds add_to_expression!(cons[Locb],2*G0[j][p,q])
                    end
                end
            end
        end
    end


    for i=1:m
        G[i]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, lblock_sigma[i])
        for j=1:lblock_sigma[i]
            if lt_block_sigma[i][j]==1
                G[i][j]=@variable(model, lower_bound=0)
                for z=1:lg[i]
                    @inbounds nota=supp_g[i][:,z]+Int64(2)*basis_sigma[i][:,block_sigma[i][j]]
                    Locb=bfind(supp_U,lsupp_U,nota,n)
                    @inbounds add_to_expression!(cons[Locb],coe_g[i][z]*G[i][j])
                end
            else
                G[i][j]=@variable(model, [1:lt_block_sigma[i][j], 1:lt_block_sigma[i][j]],PSD)
                for p=1:lt_block_sigma[i][j]
                    for q=p:lt_block_sigma[i][j]
                        for z=1:lg[i]
                            @inbounds nota=basis_sigma[i][:,block_sigma[i][j][p]]+basis_sigma[i][:,block_sigma[i][j][q]]+supp_g[i][:,z]
                            Locb=bfind(supp_U,lsupp_U,nota,n)
                            if p==q
                              @inbounds add_to_expression!(cons[Locb],coe_g[i][z]*G[i][j][p,q])
                            else
                              @inbounds add_to_expression!(cons[Locb],2*coe_g[i][z]*G[i][j][p,q])
                            end
                        end
                    end
                end
            end
        end
    end

    for i=1:l
        H[i]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, lblock_psi[i])
        for j=1:lblock_psi[i]
            if lt_block_psi[i][j]==1
                H[i][j]=@variable(model)
                for z=1:lh[i]
                    @inbounds nota=supp_h[i][:,z]+Int64(2)*basis_psi[i][:,block_psi[i][j]]
                    Locb=bfind(supp_U,lsupp_U,nota,n)
                    @inbounds add_to_expression!(cons[Locb],coe_h[i][z]*H[i][j])
                end
            else
                H[i][j]=@variable(model, [1:lt_block_psi[i][j], 1:lt_block_psi[i][j]],Symmetric)
                for p=1:lt_block_psi[i][j]
                    for q=p:lt_block_psi[i][j]
                        for z=1:lh[i]
                            @inbounds nota=basis_psi[i][:,block_psi[i][j][p]]+basis_psi[i][:,block_psi[i][j][q]]+supp_h[i][:,z]
                            Locb=bfind(supp_U,lsupp_U,nota,n)
                            if p==q
                              @inbounds add_to_expression!(cons[Locb],coe_h[i][z]*H[i][j][p,q])
                            else
                              @inbounds add_to_expression!(cons[Locb],2*coe_h[i][z]*H[i][j][p,q])
                            end
                        end
                    end
                end
            end
        end
    end

    bc=zeros(lsupp_U,1)
    for i=1:lf
        Locb=bfind(supp_U,lsupp_U,supp_f[:,i],n)
        bc[Locb]=coe_f[i]
    end
    @constraint(model, cons[2:end].==bc[2:end])
    @variable(model, lambda)
    @constraint(model, cons[1]+lambda==bc[1])
    @objective(model, Max, lambda)
    optimize!(model)

    opt_val = value(lambda)
    println("Termination status = ", termination_status(model))
    println("Optimal value = ",opt_val)

    Gr=zeros(lu0,lu0)
    for j=1:lblock_sigma0
        if lt_block_sigma0[j]==1
            Gr[block_sigma0[j][1],block_sigma0[j][1]]=value(G0[j])
        else
            Gr[block_sigma0[j],block_sigma0[j]]=value.(G0[j])
        end
    end

    sol=extract_optimizers_term(Gr,lu0,basis_sigma0,n,n,m,l,opt_val,f,g,h,x,1:n)

    return opt_val, sol

end