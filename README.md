# ctpPOP
- ctpPOP is a Julia package of solving polynomial optimization problem (POP):

**inf_{x in R^n} { f(x) : gi(x) = 0, hj(x) = 0 },**

with some special cases of inequality constraints **gi**:

Case 1: Annulus constraints on subsets of variables: **Ui >= ||x(Ti)||^2 >= Li**.

Case 2: Simlex constraints: **xi >= 0, 1 - x1 -...- xn >= 0**.

- The main idea of ctpPOP is to solve the Moment-SOS relaxation of the form:

**v = inf_X { <C,X> : X is psd, AX = b },**

which has constant trace property (CTP):

**AX = b => trace(X) = a,**

by using Conditional gradient-based augmented Lagrangian (CGAL) and Limited memory bundle method (LMBM).

- Although possibly slower than the other method on the sparse POPs, ctpPOP is much more robust on the dense ones.

- ctpPOP combines CTP with term sparity (TS), correlative sparsity (CS) and correlative sparsity-term sparsity (CS-TS) to avoid memory issue of the large-scale SDP relaxations for POPs.


# Required softwares
ctpPOP has been implemented on a desktop compute with the following softwares:
- Ubuntu 18.04.4
- Julia 1.3.1

The following sofwares are used for comparison purposes:
- [Mosek 9.1](https://www.mosek.com)
- [COSMO](https://github.com/oxfordcontrol/COSMO.jl)


# Remark
- LMBM is only supported on Ubuntu with Fortran 2018.

# Installation
- To use ctpPOP in Julia, run
```ruby
Pkg> add https://github.com/maihoanganh/ctpPOP.git
```

# Usage
The following examples briefly guide to use ctpPOP:

## Polynomial optimization
Consider the following POP on the unit ball:
```ruby
using DynamicPolynomials

@polyvar x[1:2] # variables

f=x[1]^2+0.5*x[1]*x[2]-0.25*x[2]^2+0.75*x[1]-0.3*x[2] # the objective polynomial to minimize

g=[1.0-sum(x.^2)] # the inequality constraints
h=[(x[1]-1.0)*x[2]] # the equality constraints

k=2 # relaxation order

using ctpPOP

# get information from the input data f,gi,hj
n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=ctpPOP.get_info(x,f,g,h,sparse=false);

# get the optimal value of the Moment-SOS relaxation of order k
opt_val=ctpPOP.POP_dense_CGAL(  n, # the number of variables
                                m, # the number of the inequality constraints
                                l, # the number of the equality constraints
                                lmon_g, # the number of terms in each inequality constraint
                                supp_g, # the support of each inequality constraint
                                coe_g, # the coefficients of each inequality constraint
                                lmon_h, # the number of terms in each equality constraint
                                supp_h, # the support of each equality constraint
                                coe_h, # the coefficients of each equality constraint
                                lmon_f, # the number of terms in the objective polynomial
                                supp_f, # the support of the objective polynomial
                                coe_f, # the coefficients of the objective polynomial
                                dg, # the degree of each inequality constraint
                                dh, # the degree of each equality constraint
                                k,
                                maxit=Int64(1e6), # the maximal iteration of CGAL solver
                                tol=1e-3, # the tolerance of CGAL solver
                                use_eqcons_to_get_constant_trace=false, # use the equality constraints to get constant trace
                                check_tol_each_iter=true ) # check the tolerance at each iteration
```

See other examples in the [link](https://github.com/maihoanganh/ctpPOP/tree/main/examples).


# References
For more details, please refer to:

**N. H. A. Mai, J.-B. Lasserre, V. Magron and J. Wang. Exploiting constant trace property in large-scale polynomial optimization, 2020. Forthcoming.**

The following codes are to run the paper's benchmarks:
```ruby
using ctpPOP

ctpPOP.test_dense_POP_ball(10,100,2,have_eqcons=false) # Table 4
ctpPOP.test_dense_POP_ball(10,70,2,have_eqcons=true) # Table 5

ctpPOP.test_dense_POP_annulus(10,90,have_eqcons=false) # Table 6
ctpPOP.test_dense_POP_annulus(10,70,have_eqcons=true) # Table 7

ctpPOP.test_dense_POP_box(10,70,2,have_eqcons=false) # Table 8
ctpPOP.test_dense_POP_box(10,50,2,have_eqcons=true) # Table 9

ctpPOP.test_dense_POP_simplex(10,50,2,have_eqcons=false) # Table 10
ctpPOP.test_dense_POP_simplex(10,70,2,have_eqcons=true) # Table 11

ctpPOP.test_TS_POP_ball(10,80,2,1,have_eqcons=false) # Table 12
ctpPOP.test_TS_POP_ball(10,60,2,1,have_eqcons=true) # Table 13

ctpPOP.test_TS_POP_box(10,60,2,1,have_eqcons=false) # Table 14
ctpPOP.test_TS_POP_box(10,50,2,1,have_eqcons=true) # Table 15

ctpPOP.test_CS_POP_ball(1000,11,41,2,have_eqcons=false) # Table 16
ctpPOP.test_CS_POP_ball(1000,11,36,2,have_eqcons=true) # Table 17

ctpPOP.test_CS_POP_box(1000,11,26,2,have_eqcons=false) # Table 18
ctpPOP.test_CS_POP_box(1000,11,26,2,have_eqcons=true) # Table 19

ctpPOP.test_mix_POP_ball(1000,11,26,2,1,have_eqcons=false) # Table 20
ctpPOP.test_mix_POP_ball(1000,11,26,2,1,have_eqcons=true) # Table 21

ctpPOP.test_mix_POP_box(1000,11,26,2,1,have_eqcons=false) # Table 22
ctpPOP.test_mix_POP_box(1000,11,26,2,1,have_eqcons=true) # Table 23

ctpPOP.test_comparison_dense_POP_ball(10,60,2,have_eqcons=true)



```
