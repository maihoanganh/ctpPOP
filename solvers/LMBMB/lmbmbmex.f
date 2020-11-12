*************************************************************************
*                                                                       *
*     LMBM-B - Limited Memory Bundle Method for Large-Scale             *
*              Nonsmooth Bound Constrained Optimization (version 2.0)   *
*                                                                       *
*************************************************************************
*
*
*     Codes included:
*
*     tlmbmb.f    - testprogram for LMBM-B.
*     lmbmb.f     - limited memory bundle method.
*     lmbmbs.f    - subprograms for limited memory bundle method.
*     matca2.f    - matrix and vector calculus.
*     Makefile    - makefile.
*
*     tnsboc.f    - large-scale bound constrained nonsmooth test
*                   problems.
*
*
*     The software if free for academic teaching and research purposes
*     but I ask you to refer at least one of the references given below
*     if you use it.
*
*
*     Napsu Karmitsa (maiden name Haarala) 2002 - 2005.
*     Last modified 2009.
*
*
*
*     References:
*
*     N. Karmitsa, M. M. M�kel�, "Limited Memory Bundle Method for Large
*     Bound Constrained Nonsmooth Optimization: Convergence Analysis",
*     Optimization Methods and Software, Vol. 25, No. 6, pp. 895-916,
*     2010.
*
*     N. Karmitsa, M. M. M�kel�, "Adaptive Limited Memory Bundle Method
*     for Bound Constrained Large-Scale Nonsmooth Optimization",
*     Optimization: A Journal of Mathematical Programming and Operations
*     Research, Vol. 59, No. 6, pp. 945-962, 2010.
*
*
*
*************************************************************************
*
*     * SUBROUTINE FUNDER *
*
*
*     * Purpose *
*
*     Computation of the value and the subgradient of the objective
*     function.
*
*
*     * Calling sequence *
*
*     CALL FUNDER(N,X,F,G,ITERM)
*
*
*     * Parameters *
*
*     II  N             Number of variables.
*     RI  X(N)          Vector of variables.
*     RO  F             Value of the objective function.
*     RO  G(N)          Subgradient of the objective function.
*     IO  ITERM         Cause of termination:
*                          0  - Everything is ok.
*                         -3  - Failure in function or subgradient
*                               calculations (assigned by the user).
*
*
*     * Variables in COMMON /PROB/ *
*
*     I   NEXT          Number of the test problem.
*
*
*     * Subprograms used *
*
*     S   FUNC          Computation of the value and the subgradient for
*                       problem next.
*
*
       SUBROUTINE FUNDER(N,X,F,G,ITERM)

*     Scalar Arguments
       INTEGER N,ITERM
       DOUBLE PRECISION F

*     Array Arguments
       DOUBLE PRECISION G(*),X(*)


*
*     Function and subgradient evaluation
*

       CALL objfunc(N, X, F, G)

       RETURN
       END
