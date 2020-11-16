function test_mix_POP_box(n::Int64,umin::Int64,umax::Int64,k::Int64,t::Int64;have_eqcons::Bool=true)

    @polyvar x[1:n]
    f=Polynomial{true,Float64}(x[1]+0.0)
    g=Vector{Polynomial{true,Float64}}([x[1]+0.0])
    h=Vector{Polynomial{true,Float64}}([x[1]+0.0])
    
    
    for u in umin:5:umax
        x,f,g,h=generate_mix_POP_box(n,u,have_eqcons=have_eqcons)
        println("Relaxation order: k=",k)
        println("====================")
        println("Sparse order: t=",k)
        println("====================")
        run_mix_POP(x,f,g,h,k,t)
    end
end
    
function generate_mix_POP_box(n::Int64,u::Int64;have_eqcons::Bool=false)
    
    println("***Problem setting***")
    println("Number of variable: n=",n)
    println("====================")

    @polyvar x[1:n]# variables

    function generate_random_poly(T::UnitRange{Int64})
        v=reverse(monomials(x[T],2))
        c=2*rand(Float64,length(v)).-1
        return c'*v
    end

    # unit sphere constraint
    p=floor(Int64,n/u)+1

    I=Vector{UnitRange{Int64}}(undef,p)
    I[1]=1:u
    I[2:p-1]=[u*(j-1):u*j for j in 2:p-1]
    I[p]=u*(p-1):n

    # random quadratic objective function f
    f=sum(generate_random_poly(I[j]) for j in 1:p)

    g=-(x.^2).+(1/u)

    m=length(g)
    println("Number of inequality constraints: m=",m)
    println("====================")
    
    if have_eqcons
        l=ceil(Int64, n/7)
    else
        l=0
    end
    r=floor(Int64,l/p)
    W=[(j-1)*r+1:j*r for j in 1:p-1]
    append!(W,[(p-1)*r+1:l])

    h=Vector{Polynomial{true,Float64}}(undef,l)

    #=randx=[2*rand(length(I[j])).-1 for j in 1:p]# create a feasible solution
    randx=[randx[j]/norm(randx[j]) for j in 1:p]=#

    randx=2*rand(n).-1
    randx=randx./sqrt(u)
   

    for j in 1:p
        for i in W[j]
            h[i]=generate_random_poly(I[j])
            h[i]-=h[i](x => randx) #make constraints feasible
        end
    end

    l=length(h)
    println("Number of equality constraints: l=",l)
    println("====================")


        
    return x,f,g,h
end
