module ctpPOP

using DynamicPolynomials, LinearAlgebra, MosekTools, SparseArrays, JuMP, Arpack, SumOfSquares, OSQP, LightGraphs, PolyPowerModels, TSSOS, COSMO, PowerModels, Ipopt


#export CTP_POP, ASC_PolySys



# src
include("./dense_pop/basicfuncs.jl")
include("./dense_pop/pop_dense_SOS.jl")
include("./dense_pop/pop_dense_CGAL.jl")
include("./pop_NLP.jl")
include("./pop_opf.jl")

include("./sparse_pop/CS/basic_func_sparse.jl")
include("./sparse_pop/CS/clique_merge.jl")
include("./sparse_pop/CS/chordal_extension.jl")

include("./sparse_pop/CS/pop_CS_CGAL.jl")
include("./sparse_pop/CS/basic_func_CS_CGAL.jl")
include("./sparse_pop/CS/model_cliq_CS_CGAL.jl")

include("./sparse_pop/TS/pop_TS_CGAL.jl")

include("./sparse_pop/mix/pop_mix_CGAL.jl")
include("./sparse_pop/mix/basic_func_mix_CGAL.jl")
include("./sparse_pop/mix/model_cliq_mix_CGAL.jl")




#test
include("../test/test_pop_annulus.jl")
include("../test/test_pop_ball.jl")
include("../test/test_pop_box.jl")
include("../test/test_pop_simplex.jl")

include("../test/test_correlative_pop_ball.jl")
include("../test/test_correlative_pop_box.jl")

include("../test/test_term_pop_ball.jl")
include("../test/test_term_pop_box.jl")

include("../test/test_mix_pop_ball.jl")
include("../test/test_mix_pop_box.jl")

#OPF problems
include("../test/opf_test/test_opf_case3_lmbd.jl")
include("../test/opf_test/test_opf_case5_pjm.jl")
include("../test/opf_test/test_opf_case2312_goc.jl")
include("../test/opf_test/test_opf_case1354_pegase.jl")
include("../test/opf_test/test_opf_case14_ieee.jl")
include("../test/opf_test/test_opf_case89_pegase.jl")
include("../test/opf_test/test_opf_case39_epri.jl")
end


