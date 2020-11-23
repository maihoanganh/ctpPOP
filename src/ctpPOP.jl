module ctpPOP

using Libdl, Printf, Compat

using DynamicPolynomials, LinearAlgebra, MosekTools, SparseArrays, JuMP, Arpack, SumOfSquares, LightGraphs, PolyPowerModels, TSSOS, COSMO, PowerModels#, Ipopt


#export CTP_POP, ASC_PolySys

#solvers
include("../deps/LMBM/build.jl")
include("../deps/LMBM/deps.jl")
include("../deps/LMBM/LMBM.jl")

# src
#include("./pop_NLP.jl")
include("./pop_opf.jl")

include("./dense_pop/basicfuncs.jl")
include("./dense_pop/pop_dense_SOS.jl")
include("./dense_pop/pop_dense_CGAL.jl")
include("./dense_pop/pop_dense_LMBM.jl")

include("./sparse_pop/CS/basic_func_sparse.jl")
include("./sparse_pop/CS/clique_merge.jl")
include("./sparse_pop/CS/chordal_extension.jl")

include("./sparse_pop/CS/pop_CS_CGAL.jl")
include("./sparse_pop/CS/pop_CS_LMBM.jl")
include("./sparse_pop/CS/basic_func_CS_CGAL.jl")
include("./sparse_pop/CS/model_cliq_CS_CGAL.jl")

include("./sparse_pop/TS/pop_TS_CGAL.jl")
include("./sparse_pop/TS/pop_TS_LMBM.jl")

include("./sparse_pop/mix/pop_mix_CGAL.jl")
include("./sparse_pop/mix/pop_mix_LMBM.jl")
include("./sparse_pop/mix/basic_func_mix_CGAL.jl")
include("./sparse_pop/mix/model_cliq_mix_CGAL.jl")




#test
include("../test/test_dense_pop_annulus.jl")
include("../test/test_dense_pop_ball.jl")
include("../test/test_dense_pop_box.jl")
include("../test/test_dense_pop_simplex.jl")

include("../test/test_CS_pop_ball.jl")
include("../test/test_CS_pop_box.jl")

include("../test/test_TS_pop_ball.jl")
include("../test/test_TS_pop_box.jl")

include("../test/test_mix_pop_ball.jl")
include("../test/test_mix_pop_box.jl")

include("../test/test_comparison_COSMO_dense_pop_ball.jl")

#OPF problems
include("../test/test_opf_problems.jl")

end


