function test_TS_POP_box(nvarmin::Int64,nvarmax::Int64,k::Int64,t::Int64;have_eqcons::Bool=true)
    @polyvar x[1:1]
    f=Polynomial{true,Float64}(x[1]+0.0)
    g=Vector{Polynomial{true,Float64}}([])
    h=Vector{Polynomial{true,Float64}}([])
    
    for n in nvarmin:10:nvarmax
        x,f,g,h=generate_TS_POP_box(n,have_eqcons=have_eqcons)
        println("Relaxation order: k=",k)
        println("====================")
        println("Sparse order: t=",t)
        println("====================")
        run_TS_POP(x,f,g,h,k,t)
    end
end
    
function generate_TS_POP_box(n::Int64;have_eqcons::Bool=false)
 
    println("***Problem setting***")
    println("Number of variable: n=",n)
    println("====================")

    @polyvar x[1:n]# variables

    function generate_random_poly(v)
        c=2*rand(Float64,length(v)).-1
        return c'*v
    end
    # random quadratic objective function f
    v=reverse(monomials(x,2))
    f=generate_random_poly(v)


    # unit sphere constraint
    m=n
    q=floor(Int64,n/m)
    R=ones(Float64,m)./n
    T=[(j-1)*q+1:j*q for j in 1:m-1]
    append!(T,[(m-1)*q+1:n])

    g=[R[j]-sum(x[T[j]].^2) for j in 1:m]

    m=length(g)
    println("Number of inequality constraints: m=",m)
    println("====================")

    if have_eqcons
        l=ceil(Int64, n/7)
    else
        l=0
    end

    h=Vector{Polynomial{true,Float64}}(undef,l)
    randx=[2*rand(length(T[j])).-1 for j in 1:m]# create a feasible solution
    randx=[sqrt(R[j])*rand(1)[1]*randx[j]/norm(randx[j]) for j in 1:m]
    randx=vcat(randx...)

    for j in 1:l
        h[j]=generate_random_poly(v[2:end])
        h[j]-=h[j](x => randx) #make constraints feasible
    end

    l=length(h)
    println("Number of equality constraints: l=",l)
    println("====================")

        
    return x,f,g,h
end
