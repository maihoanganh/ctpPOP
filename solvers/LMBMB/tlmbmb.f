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
*
*************************************************************************
*
*     Remark:
*
*     At the beginning of each file, there is a list of the subroutines
*     and functions included to that file. Moreover, at the beginning
*     of each subroutine, there is a description of parameters used in
*     the routine (at least those needed in the calling sequence). The
*     types of the parameters (or arguments) are introduced with two
*     letters. The first letter is either I for integer arguments or R
*     for double precision real arguments.  The second letter specifies
*     whether the argument must have a value defined on the entry to
*     the subroutine (I), whether it is a value which will be returned
*     (O), or both (U), or whether it is an auxiliary value (A). Note
*     that the arguments of the types II and RI can be changed on output
*     under some circumstances: especially, if improper input values are
*     given or if set zero. In the latter case the default values will
*     be used (if applicable).
*

*************************************************************************
*************************************************************************
*************************************************************************
*
*
*     TLMBMB includes the following subroutines:
*
*     P   TLMBMB        Test program for globally convergent limited
*                         memory bundle subroutine LMBM-B.
*     S   FUNDER        Computation of the value and the subgradient
*                         of the objective function (supplied by user).
*
*
*
*************************************************************************
*
*     * PROGRAM TLMBMB *
*
*
*     * Purpose *
*
*     Test program for globally convergent limited memory bundle
*     subroutine for large bound constrained nonsmooth optimization.
*
*
*     * Parameters *
*
*     I   N             Number of variables.
*     I   NA            Maximum bundle dimension, NA >= 2 (NA may be
*                         set to zero if IPAR(7) = 1).
*     I   MCU           Upper limit for maximum number of stored
*                         corrections, MCU >= 3.
*     I   NW            Dimension of the work vector W:
*                         NW >= 1 + 10*N + 2*N*NA + 3*NA + 2*N*MCU
*                               + 5*MCU*(MCU+1)/2 + 11*MCU
*                               + (2*MCU+1)*MCU.
*
*
*     * Variables *
*
*     I   MC            Maximum number of stored corrections,
*                         MCU >= MC >= 3.
*     R   X(N)          Vector of variables.
*     R   XL(N)         Lower bounds for variables.
*     R   XU(N)         Upper bounds for variables.
*     I   IB(N)         Type of bound constraints:
*                         0  - X(I) is unbounded,
*                         1  - X(I) has only a lower bound,
*                         2  - X(I) has both lower and upper bounds,
*                         3  - X(I) has only an upper bound.
*     I   IACT(N)       Index set of active and free variables.
*     R   F             Value of the objective function.
*     R   TIME          Maximum CPU-time in seconds. If TIME <= 0.0
*                         the maximum time is ignored. REAL argument.
*     R   RTIM(2)       Auxiliary array. REAL array.
*                         On output RTIM(1) contains the CPU-time used.
*     R   RPAR(8)       Real parameters:
*           RPAR(1)       Tolerance for change of function values.
*           RPAR(2)       Second Tolerance for change of function values.
*           RPAR(3)       Tolerance for the function value.
*           RPAR(4)       Tolerance for the first termination criterion.
*           RPAR(5)       Tolerance for the second termination criterion.
*           RPAR(6)       Distance measure parameter, 0 <= RPAR(6).
*           RPAR(7)       Line search parameter, 0 < RPAR(7) < 0.25.
*           RPAR(8)       Maximum stepsize, 1 < RPAR(8).
*                           If RPAR(I) <= 0 for I=1,3,4,5,7, and 8 the
*                           default value of the parameter will be used.
*                           If RPAR(2) < 0 the the parameter and the
*                           corresponding termination criterion will be
*                           ignored. If RPAR(2) = 0 default value will
*                           be used. If RPAR(6) < 0 the default value
*                           will be used.
*     I   IPAR(7)       Integer paremeters:
*           IPAR(1)       Exponent for distance measure.
*           IPAR(2)       Maximum number of iterations.
*           IPAR(3)       Maximum number of function evaluations.
*           IPAR(4)       Maximum number of iterations with changes of
*                           function values smaller than RPAR(1).
*           IPAR(5)       Printout specification:
*                            -1  - No printout.
*                             0  - Only the error messages.
*                             1  - The final values of the objective
*                                  function.
*                             2  - The final values of the objective
*                                  function and the most serious warning
*                                  messages.
*                             3  - The whole final solution.
*                             4  - At each iteration values of the
*                                  objective function.
*                             5  - At each iteration the whole solution.
*           IPAR(6)       Selection of the scaling:
*                             0  - Scaling at every iteration with STU/UTU.
*                             1  - Scaling at every iteration with STS/STU.
*                             2  - Interval scaling with STU/UTU.
*                             3  - Interval scaling with STS/STU.
*                             4  - Preliminary scaling with STU/UTU.
*                             5  - Preliminary scaling with STS/STU.
*                             6  - No scaling.
*           IPAR(7)       Selection of initial stepsize for line search
*                           procedure:
*                             0  - Stepsize selection using polyhedral
*                                  approximation (NA >= 2).
*                             1  - Stepsize = min(1.0,TMAX), where TMAX is
*                                  the upper limit for step size assuring
*                                  the feasibility of produced point (no
*                                  additional bundle is needed, NA may be
*                                  set to zero).
*                         If IPAR(I) <= 0 the default value of the
*                         parameter will be used.
*     I   IOUT(4)         Integer parameters:
*           IOUT(1)         Number of used iterations.
*           IOUT(2)         Number of used function evaluations.
*           IOUT(3)         Cause of termination:
*                              1  - The problem has been solved.
*                                   with desired accuracy.
*                              2  - Changes in function values < RPAR(1)
*                                   in IPAR(4) subsequent iterations.
*                              3  - Changes in function values < RPAR(2)
*                                   *SMALL*MAX(|F_k|,|F_k+1|,1), where
*                                   SMALL is the smallest positive
*                                   number such that 1.0 + SMALL > 1.0.
*                              4  - Number of function calls > IPAR(3).
*                              5  - Number of iterations > IPAR(2).
*                              6  - Time limit exceeded.
*                              7  - F < RPAR(3).
*                             -1  - Two consecutive restarts or the value
*                                   of the function remains the same
*                                   between two sequential restarts.
*                             -2  - Number of restarts > maximum number
*                                   of restarts.
*                             -3  - Failure in function or subgradient
*                                   calculations (assigned by the user).
*                             -4  - Failure in attaining the demanded
*                                   accuracy.
*                             -5  - Invalid input parameters.
*                             -6  - Not enough working space.
*           IOUT(4)         Number of active variables at solution.
*     R   W(NW)           Work vector.
*
*
*     * Variables in COMMON /PROB/ *
*
*     I   NEXT            Number of the test problem.
*
*
*
*     * Subprograms used *
*
*     S   BOUNDS          Definition of bound constraints.
*     S   LMBMBI          Initialization of LMBM-B for nonsmooth
*                           optimization.
*     S   STARTX          Initiation of X.
*
*
*

      PROGRAM TLMBMB


*     Parameters
      INTEGER N,NA,MCU,NW
      PARAMETER(
     &     N = 1000,
     &     NA = 2,
     &     MCU = 15,
     &     NW = 1 + 13*N + 2*N*NA + 3*NA + 2*N*MCU
     &     + 5*MCU*(MCU+1)/2 + 13*MCU + (2*MCU+1)*MCU)


*     Scalar Arguments
      INTEGER MC
      DOUBLE PRECISION F


*     Array Arguments
      INTEGER IPAR(7),IOUT(4),IB(N),IACT(N)
      DOUBLE PRECISION W(NW),X(N),XL(N),XU(N),RPAR(8)


*     Local Scalars
      INTEGER MCINIT,I,J


*     CPU-time
      REAL TIME,RTIM(2)


*     Scalars depending on problem
      INTEGER NEXT
      COMMON /PROB/NEXT


*     External Subroutines
      EXTERNAL LMBMBI,STARTX,BOUNDS


*
*     Maximum time
*

      TIME=1800.0E+00


*
*     Loop for test problems
*

      DO 100 J=1,10


*
*     Initial number of stored corrections
*

      MC = 7
      MCINIT = MC


*
*     Number for the test problems
*

      NEXT=J


*
*     Initiation of X and bounds.
*

      CALL STARTX(N,X,NEXT)
      IF (NEXT .EQ. -1) GOTO 999

      CALL BOUNDS(N,IB,XL,XU,NEXT)

c     no need for this one
c      CALL FEASIX(N,X,XL,XU,IB,0)


*
*     Choice of integer parameters
*

      DO 10 I = 1,7
         IPAR(I) = 0
 10   CONTINUE


*     Maximum numbers of iterations and function evaluations

      IPAR(2) = 5000000
      IPAR(3) = 5000000


*     Printout specification

      IPAR(5) = 1


*     Selection of the scaling

      IPAR(6) = 0


*     Selection of stepsize

c     If NA /= 0, the bundle will be construated. If
c     IPAR(7) = 0, it will be used at every iteration.
c     Otherwise, it will be used only if WK = 0

c      IPAR(7) = 0
      IPAR(7) = 1


*
*     Choice of real parameters
*

      DO 20 I = 1,8
         RPAR(I) = 0.0D0
 20   CONTINUE


*     Second Tolerance for change of function values.

      RPAR(2)= 1.0D+03


*     Desired accuracy

      RPAR(4) = 1.0D-05
      RPAR(5) = 1.0D+06


*     Locality measure

      IF (NEXT .GT. 5) RPAR(6) = 0.50D+00


*     Line search parameter

      RPAR(7) = 0.0001D+00


*     Stepsize

      RPAR(8) = 1.50D+00


*
*     Solution
*

      CALL LMBMBI(N,NA,MC,MCU,NW,X,XL,XU,F,IB,IACT,IPAR,IOUT,RPAR,TIME,
     &     RTIM,W)


*
*     Result (additional printout)
*

      PRINT*
      PRINT*,'NEXT    = ',NEXT
      PRINT*,'ITERM   = ',IOUT(3)
      PRINT*
      PRINT*,'F(X)    = ',F
      DO 30 I=1,N
         IF (IB(I) .NE. 0 .AND. IB(I) .LE. 2 .AND.
     &        X(I)-XL(I) .LE. -0.10D-12) THEN
            PRINT*,'XL = ',I,X(I),XL(I),X(I)-XL(I)
         END IF
         IF (IB(I) .GE. 2 .AND. XU(I)-X(I) .LE. -0.10D-12) THEN
            PRINT*,'XU = ',I,X(I),XU(I),XU(I)-X(I)
         END IF
 30   CONTINUE
      PRINT*
      PRINT*,'N       = ',N
      PRINT*,'NA      = ',NA
      PRINT*,'NACT    = ',IOUT(4)
      PRINT*,'MCINIT  = ',MCINIT
      PRINT*,'MC      = ',MC
      PRINT*,'MCU     = ',MCU
      PRINT*,'NIT     = ',IOUT(1)
      PRINT*,'NFV     = ',IOUT(2)
      PRINT*,'XMAX    = ',RPAR(8)
      PRINT*,'GAM     = ',RPAR(6)
      PRINT*,'EPSL    = ',RPAR(7)
      PRINT*,'EPS     = ',RPAR(4)
      PRINT*,'SCALING = ',IPAR(6)
      PRINT*,'ISTEP   = ',IPAR(7)

      PRINT*
      PRINT*,'Used time = ',RTIM(1)
      PRINT*

 100  CONTINUE


 999  CONTINUE

      STOP
      END



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

*     External Subroutines
      EXTERNAL FUNC

*     Common blocks
      INTEGER NEXT
      COMMON /PROB/NEXT


      ITERM=0


*
*     Function and subgradient evaluation
*

      CALL FUNC(N,X,F,G,NEXT)

      IF (NEXT .LT.1) ITERM = -3

      RETURN
      END
