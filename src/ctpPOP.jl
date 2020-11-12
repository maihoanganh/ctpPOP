module ctpPOP

using DynamicPolynomials, LinearAlgebra, MosekTools, SparseArrays, JuMP, Arpack, RowEchelon, SumOfSquares, Libdl, Printf, Compat, OSQP, Distributed, ArnoldiMethod,LightGraphs, PolyPowerModels, TSSOS, COSMO, PowerModels, Ipopt


#export CTP_POP, ASC_PolySys


#solvers
include("../solvers/LMBM/build.jl")
include("../solvers/LMBM/deps.jl")
include("../solvers/LMBM/LMBM.jl")

#include("../solvers/ProximalBundleMethod/ProximalMethod.jl")

# src
include("./pop/basicfuncs.jl")

#pop
include("./pop/dense_pop/POP_SumOfSquares.jl")
include("./pop/dense_pop/pop_CGAL.jl")
include("./pop/dense_pop/pop_DSMA.jl")
include("./pop/dense_pop/pop_LMBM.jl")

include("./pop/dense_pop/pop_CGAL_arb_cons.jl")
include("./pop/dense_pop/pop_LMBM_arb_cons.jl")

include("./pop/sparse_pop/basic_func_sparse.jl")

include("./pop/sparse_pop/clique_merge.jl")
include("./pop/sparse_pop/chordal_extension.jl")

include("./pop/sparse_pop/old_funcs/correlative_pop_single_ball.jl")

include("./pop/sparse_pop/CS/CGAL/pop_CS_CGAL.jl")
include("./pop/sparse_pop/CS/CGAL/basic_func_CS_CGAL.jl")
include("./pop/sparse_pop/CS/CGAL/model_cliq_CS_CGAL.jl")

include("./pop/sparse_pop/CS/LMBM/pop_CS_LMBM.jl")
include("./pop/sparse_pop/CS/LMBM/basic_func_CS_LMBM.jl")
include("./pop/sparse_pop/CS/LMBM/model_cliq_CS_LMBM.jl")

include("./pop/sparse_pop/TS/pop_TS_LMBM.jl")
include("./pop/sparse_pop/TS/pop_TS_CGAL.jl")
include("./pop/sparse_pop/TS/pop_TS_Mosek.jl")

include("./pop/sparse_pop/mix/CGAL/pop_mix_CGAL.jl")
include("./pop/sparse_pop/mix/CGAL/basic_func_mix_CGAL.jl")
include("./pop/sparse_pop/mix/CGAL/model_cliq_mix_CGAL.jl")

include("./pop/sparse_pop/mix/LMBM/pop_mix_LMBM.jl")
include("./pop/sparse_pop/mix/LMBM/basic_func_mix_LMBM.jl")
include("./pop/sparse_pop/mix/LMBM/model_cliq_mix_LMBM.jl")

include("./pop/pop_opf.jl")
#pmsv
include("../solvers/LMBMB/build.jl")
include("../solvers/LMBMB/deps.jl")
include("../solvers/LMBMB/LMBMB.jl")

include("./pmsv/PSV.jl")
include("./pmsv/pmsv_mosek.jl")
include("./pmsv/pmsv_dsma.jl")

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

end


