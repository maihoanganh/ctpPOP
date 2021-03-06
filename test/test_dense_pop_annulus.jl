function test_dense_POP_annulus(nvarmin::Int64,nvarmax::Int64,k::Int64;have_eqcons::Bool=true)
    @polyvar x[1:1]
    f=Polynomial{true,Float64}(x[1]+0.0)
    g=Vector{Polynomial{true,Float64}}([])
    h=Vector{Polynomial{true,Float64}}([])
    
    for n in nvarmin:10:nvarmax
        x,f,g,h=generate_dense_POP_annulus(n,have_eqcons=have_eqcons)
        println("Relaxed order: k=",k)
        println("====================")
        run_dense_POP(x,f,g,h,k)
    end
end
    
function generate_dense_POP_annulus(n::Int64;have_eqcons::Bool=false)
    println("***Problem setting***")


    println("Number of variable: n=",n)
    println("====================")

    @polyvar x[1:n]# variables

    function generate_random_poly(v)
        c=2*rand(Float64,length(v)).-1
        return c'*v
    end
    # random quadratic objective function f
    v=reverse(monomials(x,0:2))
    f=generate_random_poly(v)


    # unit sphere constraint
    R_small=0.5
    R_big=1.0

    g=[sum(x.^2)-R_small;R_big-sum(x.^2)] #type of coefficients of each polynomial must be float

    m=length(g)
    println("Number of inequality constraints: m=",m)
    println("====================")
    
    if have_eqcons
        l=ceil(Int64, n/4)
        h=Vector{Polynomial{true,Float64}}(undef,l)
    else
        l=0
        h=Vector{Polynomial{true,Float64}}([])
    end
    
    
    randx=2*rand(n).-1# create a feasible solution
    randx=(sqrt(R_small)+(sqrt(R_big)-sqrt(R_small))*rand(1)[1])*randx./sqrt(sum(randx.^2))


    for j in 1:l
        h[j]=generate_random_poly(v[2:end])
        h[j]-=h[j](x => randx) #make constraints feasible
    end

    l=length(h)
    println("Number of equality constraints: l=",l)
    println("====================")

        
    return x,f,g,h
end
