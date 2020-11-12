*************************************************************************
*
*
*     LMBMB includes the following subroutines:
*
*     S   LMBMBI          Initialization for LMBM-B.
*     S   LMBMB           LMBM-B.
*
*
*     Napsu Karmitsa (maiden name Haarala) 2002 - 2005. 
*     Last modified 2009.
*
*************************************************************************
*      
*     * SUBROUTINE LMBMBI *
*
*      
*     * Purpose *
*
*     Initialization for globally convergent LMBM-B for large bound 
*     constrained nonsmooth optimization.
*
*     
*     * Calling sequence *
*     
*     CALL LMBMBI(N,NA,MC,MCU,NW,X,XL,XU,F,IB,IACT,IPAR,IOUT,RPAR,TIME,
*    &     RTIM,W)
*     
*     
*     * Parameters *
*      
*     II  N             Number of variables.
*     II  NA            Maximum bundle dimension, NA >= 2.
*     II  MCU           Upper limit for maximum number of stored
*                         corrections, MCU >= 3.
*     IU  MC            Maximum number of stored corrections, MC <= MCU.
*     RU  X(N)          Vector of variables.
*     RI  XL(N)         Lower bounds for variables. 
*     RI  XU(N)         Upper bounds for variables. 
*     II  IB(N)         Type of bound constraints:
*                         0  - X(I) is unbounded,
*                         1  - X(I) has only a lower bound,
*                         2  - X(I) has both lower and upper bounds,
*                         3  - X(I) has only an upper bound. 
*     IA  IACT(N)       Index set of active and free variables.
*     RO  F             Value of the objective function.
*     RI  TIME          Maximum CPU-time in seconds. If TIME <= 0.0
*                         the maximum time is ignored. REAL argument.
*     RI  RTIM(2)       Auxiliary array. REAL array
*                         On output RTIM(1) contains the execution time.
*     RI  RPAR(8)       Real parameters:
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
*     II  IPAR(7)       Integer paremeters:
*           IPAR(1)       Exponent for distance measure.
*           IPAR(2)       Maximum number of iterations.
*           IPAR(3)       Maximum number of function evaluations.
*           IPAR(4)       Maximum number of iterations with changes of
*                           function values smaller than RPAR(1).
*           IPAR(5)       Printout specification:
*                          -1  - No printout.
*                           0  - Only the error messages.
*                           1  - The final values of the objective
*                                function.
*                           2  - The final values of the objective
*                                function and the most serious warning 
*                                messages.
*                           3  - The whole final solution. 
*                           4  - At each iteration values of the
*                                objective function.
*                           5  - At each iteration the whole solution.
*           IPAR(6)       Selection of the scaling:
*                           0  - Scaling at every iteration with STU/UTU.
*                           1  - Scaling at every iteration with STS/STU.
*                           2  - Interval scaling with STU/UTU.
*                           3  - Interval scaling with STS/STU.
*                           4  - Preliminary scaling with STU/UTU.
*                           5  - Preliminary scaling with STS/STU.
*                           6  - No scaling.      
*           IPAR(7)       Selection of initial stepsize before line 
*                           search procedure:
*                             0  - Stepsize selection using polyhedral 
*                                  approximation (NA >= 2).
*                             1  - Stepsize = min(1.0,TMAX), where TMAX is
*                                  the upper limit for step size assuring 
*                                  the feasibility of produced point (no 
*                                  additional bundle is needed, NA may be 
*                                  set to zero).
*                           If IPAR(I) <= 0 the default value of the
*                           parameter will be used.
*     IO  IOUT(4)       Integer parameters:
*           IOUT(1)       Number of used iterations.
*           IOUT(2)       Number of used function evaluations.
*           IOUT(3)       Cause of termination:
*                           1  - The problem has been solved with desired
*                                accuracy.
*                           2  - Changes in function values < RPAR(1)
*                                in IPAR(4) subsequent iterations.
*                           3  - Changes in function values < RPAR(2)
*                                *SMALL*MAX(|F_k|,|F_k+1|,1), where
*                                SMALL is the smallest positive number
*                                such that 1.0 + SMALL > 1.0.
*                           4  - Number of function calls > IPAR(3).
*                           5  - Number of iterations > IPAR(2).
*                           6  - Time limit exceeded. 
*                           7  - F < RPAR(3).
*                          -1  - Two consecutive restarts or the value 
*                                of the function remains the same 
*                                between two sequential restarts.
*                          -2  - Number of restarts > maximum number
*                                of restarts.
*                          -3  - Failure in function or subgradient
*                                calculations (assigned by the user).
*                          -4  - Failure in attaining the demanded
*                                accuracy.
*                          -5  - Invalid input parameters.
*                          -6  - Not enough working space.
*           IOUT(4)       Number of active variables at solution.
*     RA  W(NW)         Work vector.
*     II  NW            Dimension of the work vector W:
*                         NW >= 1 + 10*N + 2*N*NA + 3*NA + 2*N*MCU
*                               + 5*MCU*(MCU+1)/2 + 11*MCU
*                               + (2*MCU+1)*MCU.
*     
*
*
*     * Subprograms used *
*      
*     S   LMBMB         Limited memory bundle method LMBM-B for nonsmooth
*                         bound constrained optimization.
*     S   WPRINT        Printout the error and warning messages.
*     S   PROJX         Projection of the initial X to the feasible 
*                         region if necessary.
*     S   GETIME       Execution time.
*
*
*
      
      SUBROUTINE LMBMBI(N,NA,MC,MCU,NW,X,XL,XU,F,IB,IACT,IPAR,IOUT,
     &     RPAR,TIME,RTIM,W)


*     Scalar Arguments
      INTEGER N,NA,MC,MCU,NW
      DOUBLE PRECISION F

      
*     Array Arguments
      INTEGER IPAR(*),IOUT(*),IB(*),IACT(*)
      DOUBLE PRECISION X(*),XL(*),XU(*),RPAR(*),W(*)

      
*     Local Scalars
      INTEGER LXCP,LXO,LS,LG,LGM,LGA,LPGA,LPGM,LU,LD,LAX,LAG,LAF,LSM,
     &     LUM,LRM,LLM,LUMTUM,LSMTSM,LC,LSMTGM,LUMTGM,LSMTPG,LUMTPG,
     &     LVMC1,LVMC2,LVMC3,LVMC4,LVMC5,LVMC6,LVMC7,LVMC8,LPGAH,
     &     LVN1,LVN2,LTMAT,LKMAT,ITYPE,ITMP


*     CPU-time
      REAL TIME,START,FINI
      REAL RTIM(2)


*     External Subroutines
      EXTERNAL LMBMB,WPRINT,PROJX,GETIME


*     
*     CPU-time
*      

      CALL GETIME(START,RTIM)


*     
*     Initialization and error checking
*

      IOUT(3) = 0
      
      IF (N .LE. 0) THEN
         IOUT(3) = -5
         CALL WPRINT(IOUT(3),IPAR(5),1)
         RETURN
      END IF


      IF (MCU .LT. 3) THEN
         IOUT(3) = -5
         CALL WPRINT(IOUT(3),IPAR(5),2)
         RETURN
      END IF


      IF (IPAR(7) .GE. 1) THEN
         IF (NA .LT. 0) THEN
            IOUT(3) = -5
            CALL WPRINT(IOUT(3),IPAR(5),3)
            RETURN
         END IF
      ELSE IF (NA .LT. 2) THEN
         IOUT(3) = -5
         CALL WPRINT(IOUT(3),IPAR(5),3)
         RETURN
      END IF

      
      IF (RPAR(7) .GE. 0.25D+00) THEN
         IOUT(3) = -5
         CALL WPRINT(IOUT(3),IPAR(5),4)
         RETURN
      END IF


      ITMP = 1 + 13*N + 2*N*NA + 3*NA + 2*N*MCU
     &     + 5*MCU*(MCU+1)/2 + 13*MCU + (2*MCU+1)*MCU

      IF (NW .LT. ITMP) THEN
         IOUT(3) = -6
         CALL WPRINT(IOUT(3),IPAR(5),0)
         RETURN
      END IF


      IF (IPAR(6) .GT. 6 .OR. IPAR(6) .LT. 0) IPAR(6) = 0

      
      IF (MC .GT. MCU) THEN
         MC = MCU
         CALL WPRINT(IOUT(3),IPAR(5),-1)
      END IF
         
      IF (MC .LE. 0) MC = 3


      CALL PROJX(N,X,XL,XU,IB,IOUT(4),ITYPE)
      IF (IOUT(4) .NE. 0) THEN
         CALL WPRINT(IOUT(3),IPAR(5),IOUT(4))
      END IF
      

*     
*     Pointers for working array W
*

      LXCP   = 1
      LXO    = LXCP   + N
      LS     = LXO    + N
      LG     = LS     + N
      LGM    = LG     + N
      LGA    = LGM    + N
      LPGA   = LGA    + N
      LPGM   = LPGA   + N
      LU     = LPGM   + N
      LD     = LU     + N
      LAX    = LD     + N
      LAG    = LAX    + N*NA
      LAF    = LAG    + N*NA
      LSM    = LAF    + 3*NA
      LUM    = LSM    + N*MCU
      LRM    = LUM    + N*MCU
      LLM    = LRM    + MCU*(MCU+1)/2
      LUMTUM = LLM    + MCU*(MCU+1)/2
      LSMTSM = LUMTUM + MCU*(MCU+1)/2
      LC     = LSMTSM + MCU*(MCU+1)/2
      LSMTGM = LC     + MCU
      LUMTGM = LSMTGM + MCU
      LSMTPG = LUMTGM + MCU
      LUMTPG = LSMTPG + MCU
      LVMC1  = LUMTPG + MCU
      LVMC2  = LVMC1  + MCU
      LVMC3  = LVMC2  + MCU
      LVMC4  = LVMC3  + MCU
      LVMC5  = LVMC4  + MCU
      LVMC6  = LVMC5  + MCU
      LVMC7  = LVMC6  + MCU
      LVMC8  = LVMC7  + MCU
      LPGAH  = LVMC8  + MCU
      LVN1   = LPGAH  + N
      LVN2   = LVN1   + N
      LTMAT  = LVN2   + N
      LKMAT  = LTMAT  + MCU*(MCU+1)/2
*     size of LKMAT = (2*MCU+1)*MCU)


*     
*     Solution
*

      CALL LMBMB(N,NA,MC,MCU,IOUT(4),X,XL,XU,ITYPE,IB,IACT,W(LXCP),
     &     W(LXO),W(LS),W(LG),W(LGM),W(LGA),W(LPGA),W(LPGM),W(LU),
     &     W(LD),F,W(LAX),W(LAG),W(LAF),W(LSM),W(LUM),W(LRM),W(LLM),
     &     W(LUMTUM),W(LSMTSM),W(LC),W(LSMTGM),W(LUMTGM),W(LSMTPG),
     &     W(LUMTPG),W(LVMC1),W(LVMC2),W(LVMC3),W(LVMC4),W(LVMC5),
     &     W(LVMC6),W(LVMC7),W(LVMC8),W(LPGAH),W(LVN1),W(LVN2),
     &     W(LTMAT),W(LKMAT),RPAR(1),RPAR(2),RPAR(3),RPAR(4),RPAR(5),
     &     RPAR(6),RPAR(7),RPAR(8),IPAR(2),IPAR(3),IPAR(1),IPAR(4),
     &     IPAR(5),IPAR(6),IPAR(7),IOUT(1),IOUT(2),IOUT(3),TIME,RTIM)


*     
*     CPU-time
*      

      CALL GETIME(FINI,RTIM)
      RTIM(1) = FINI - START


      RETURN
      END


*************************************************************************
*
*     * SUBROUTINE LMBMB *
*
*      
*     * Purpose *
*      
*     Limited memory bundle subroutine for large bound constrained 
*     nonsmooth optimization.
*
*      
*     * Calling sequence *
*     
*     CALL LMBMB(N,NA,MC,MCU,NACT,X,XL,XU,ITYPE,IB,IACT,XCP,XO,S,
*    &     G,GM,GA,U,D,F,AX,AG,AF,SM,UM,RM,LM,UMTUM,SMTSM,CDIAG,SMTGM,
*    &     UMTGM,VMC1,VMC2,VMC3,VMC4,VMC5,VMC6,VMC7,VMC8,PGAH,
*    &     VN1,VN2,TMPMAT,KMAT,TOLF,TOLF2,TOLB,TOLG,TOLG2,ETA,EPSL,
*    &     XMAX,MIT,MFE,MOS,MTESF,IPRINT,ISCALE,NIT,NFE,ITERM,TIME,
*    &     RTIM)
*
*      
*     * Parameters *
*      
*     II  N               Number of variables.
*     II  NA              Maximum bundle dimension.
*     II  MCU             Upper limit for maximum number of stored
*                           corrections, MCU >= 3.
*     IU  MC              Maximum number of stored corrections, MC<=MCU.
*     RI  TIME            Maximum CPU-time in seconds. If TIME <= 0.0
*                           the maximum time is ignored. REAL argument.
*     RI  RTIM(2)         Auxiliary array. REAL array.
*                           On input RTIM(1) contains the starting time.
*     IU  NACT            Number of active variables.
*     RU  X(N)            Vector of variables.
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables. 
*     II  ITYPE           Type of problem:
*                           0  - problem is unbounded,
*                           1  - constrained problem.
*     II  IB(N)           Type of bound constraints:
*                           0  - X(I) is unbounded,
*                           1  - X(I) has only a lower bound,
*                           2  - X(I) has both lower and upper bounds,
*                           3  - X(I) has only an upper bound. 
*     IA  IACT(N)         Index set of active and free variables.
*     RA  XCP(N)          Generalized Cauchy point.
*     RA  XO(N)           Previous vector of variables.
*     RA  G(N)            Subgradient of the objective function.
*     RA  GM(N)           Previous subgradient of the objective function.
*     RA  GA(N)           Aggregate subgradient.
*     RA  PGA(N)          Projected aggregate subgradient.
*     RA  PGM(N)          Projected subgradient of the objective function.
*     RA  S(N)            Difference of current and previous variables.
*     RA  U(N)            Difference of current and previous
*                           subgradients.
*     RA  D(N)            Search direction.
*     RO  F               Value of the objective function.
*     RA  AX(N*NA)        Matrix whose columns are bundle points.
*     RA  AG(N*NA)        Matrix whose columns are bundle subgradients.
*     RA  AF(3*NA)        Vector of bundle values.
*     RA  SM(N*MC)        Matrix whose columns are stored differences of
*                           variables.
*     RA  UM(N*MC)        Matrix whose columns are stored subgradient
*                           differences.
*     RA  LM(MC*(MC+1)/2) Lower triangular matrix stored rowwise in
*                           one-dimensional array.
*     RA  RM(MC*(MC+1)/2) Upper triangular matrix stored columnwise in
*                           one-dimensional array.
*     RA  UMTUM(MC*(MC+1)/2)  Auxiliary matrix: UMTUM = UM'*UM.
*     RA  SMTSM(MC*(MC+1)/2)  Auxiliary matrix: SMTSM = SM'*SM.
*     RA  CDIAG(MC)       Diagonal matrix.
*     RA  PGAH            Auxiliary matrix: PGAH = -H * PGA.
*     RA  SMTGM(MC)       Auxiliary vector.
*     RA  UMTGM(MC)       Auxiliary vector.
*     RA  VMC#(MC)        Auxiliary arrays; # = 1,...,8.
*     RA  VN#(N)          Auxiliary arrays; # = 1,2.
*     RA  TMPMAT(MC*(MC+1)/2)  Auxiliary matrix.
*     RA  KMAT((2*MC+1)*MC)    Auxiliary 2*mc matrix.
*     RI  TOLF            Tolerance for change of function values.
*     RI  TOLF2           Second tolerance for change of function
*                           values.
*     RI  TOLB            Tolerance for the function value.
*     RI  TOLG            Tolerance for the first termination criterion.
*     RI  TOLG2           Tolerance for the second termination criterion.
*     RI  ETA             Distance measure parameter, ETA = 0 if all the
*                           functions involved are convex.
*     RI  EPSL            Line search parameter, 0 < EPSL < 0.25.
*     RI  XMAX            Maximum stepsize, 1 < XMAX.
*     II  MIT             Maximun number of iterations.
*     II  MFE             Maximun number of function evaluations.
*     II  MOS             Exponent for distance measure.
*     II  MTESF           Maximum number of iterations with changes of
*                           function values smaller than TOLF.
*     II  IPRINT          Printout specification:
*                          -1  - No printout.
*                           0  - Only the error messages.
*                           1  - The final values of the objective
*                                function.
*                           2  - The final values of the objective
*                                function and the most serious
*                                warning messages.
*                           3  - The whole final solution. 
*                           4  - At each iteration values of the
*                                objective function.
*                           5  - At each iteration the whole solution.
*     II  ISCALE          Selection of the scaling:
*                           0  - Scaling at every iteration with STU/UTU.
*                           1  - Scaling at every iteration with STS/STU.
*                           2  - Interval scaling with STU/UTU.
*                           3  - Interval scaling with STS/STU.
*                           4  - Preliminary scaling with STU/UTU.
*                           5  - Preliminary scaling with STS/STU.
*                           6  - No scaling.      
*     II  ISTEP           Selection of initial stepsize before line 
*                         search procedure:
*                           0  - Stepsize selection using polyhedral 
*                                approximation (NA >= 2).
*                           1  - Stepsize = min(1.0,TMAX), where TMAX is
*                                the upper limit for step size assuring 
*                                the feasibility of produced point (no 
*                                additional bundle is needed, NA may be 
*                                set to zero).
*     IO  NIT             Number of used iterations.
*     IO  NFE             Number of used function evaluations.
*     IO  ITERM           Cause of termination:
*                           1  - The problem has been solved with 
*                                the desired accuracy.
*                           2  - (F - FO) < TOLF in MTESF subsequent
*                                iterations.
*                           3  - (F - FO) < TOLF2*SMALL*MAX(|F|,|FO|,1).
*                           4  - Number of function calls > MFE.
*                           5  - Number of iterations > MIT.
*                           6  - Time limit exceeded. 
*                           7  - F < TOLB.
*                          -1  - Two consecutive restarts or the value 
*                                of the function remains the same 
*                                between two sequential restarts.
*                          -2  - Number of restarts > maximum number
*                                of restarts.
*                          -3  - Failure in function or subgradient
*                                calculations (assigned by the user).
*                          -4  - Failure in attaining the demanded
*                                accuracy.
*                          -5  - Invalid input parameters.
*                          -6  - Not enough working space.
*
*
*
*     * Local parameters *
*
*     I   MAXEPS          Maximum number of consecutive equal stopping
*                           criterions.
*     I   MAXNRS          Maximum number of restarts.
*     R   ETA9            Maximum for real numbers.
*     R   FMIN            Smallest acceptable value of the function.
*     R   TMIN            Minimum stepsize.
*     R   LENGTHD         Direction vector length.
*     R   RHO             Correction parameter.
*     
*      
*     * Local variables *
*
*     I   MN              Current number of stored corrections.
*     I   INEW            Index for the circular arrays.
*     I   IOLD            Index for the oldest element in circular 
*                           arrays.
*     I   IBUN            Index for the circular arrays in bundle.
*     I   IBFGS           Index of the type of BFGS update.
*                           0  - SR1 update or not updated yet.
*                           1  - BFGS update: the corrections are
*                                  stored.
*                           3  - BFGS update is skipped.
*     I   ISR1            Index of the type of SR1 update.
*     I   ITERS           Null step indicator.
*                           0  - Null step.
*                           1  - Serious step.
*     I   ICN             Correction indicator for null steps.
*     I   ISGNGA          ISGNGA = 0, if GA(I)<=0 for all I such that 
*                           X(I)=XU(I) and GA(I)>=0 for all I such that
*                           X(I)=XL(I). Otherwise ISGNGA = 1.
*     I   IFLAG           Index for adaptive version:
*                           0  - Maximum number of stored corrections
*                                has not been changed.
*                           1  - Maximum number of stored corrections
*                                has been changed.
*     I   IERR            Error indicador: 
*                           0  - Everything is ok.
*                          -1  - Error in CHOFA; Restart.
*                          -2  - Error in LINEQ; Restart.
*                          -3  - Error in TRLIEQ; Restart.
*                          -4  - Error in FRMLEL; Restart.
*                          -5  - Error in PRDKMV; Restart.
*                         -10  - Warning: MN=0 and IBFGS=3; Restart.
*                         -11  - Warning: MN=0 and ISR1=3; Restart.
*     I   NAC             Current size of the bundle.
*     I   NEPS            Number of consecutive equal stopping
*                           criterions.
*     I   NNK             Consecutive null steps counter.
*     I   NCRES1          Number of consecutive restarts.
*     I   NCRES2          Number of consecutive restarts in case of
*                           TMAX < TMIN.
*     I   NCRES3          Number of consecutive restarts in case of
*                           DNORM < SMALL.
*     I   NRES            Number of restars.      
*     I   NTESF           Number of tests on function decrease.
*     I   NTESF2          Number of tests on function decrease after 
*                           restart.
*     R   BETA            Locality measure.
*     R   AGBETA          Aggregate locality measure.
*     R   EPSR            Line search parameter.
*     R   GAMMA           Scaling parameter.
*     R   AMUGAD          AMUGAD = -(A*MU + GA)'* D, where A*MU denotes
*                           the Lagrange multipliers for problem.
*     R   XBX             XBX = (XCP-X)'B(XCP-X).
*     R   ALPHA           Backtrack multiplier.
*     R   P               Directional derivative.
*     R   FO              Previous value of objective function.
*     R   DNORM           Euclidean norm of the direction vector.
*     R   DSTRN           Euclidean norm of the direction vector without 
*                           backtracking.
*     R   PGANRM          Euclidean norm of the projected aggregate 
*                           subgradient vector.
*     R   WK              Stopping criterion.
*     R   PWK             Previous stopping criterion.
*     R   QK              Second stopping criterion.
*     R   T               Stepsize.
*     R   TMAX            Maximum stepsize.
*     R   THETA           Correction parameter for stepsize.
*     R   TTHETA          TTHETA = T * THETA.
*     R   SMALL           The smallest positive number such that
*                           1.0 + SMALL > 1.0.
*
*     
*      
*     * Subprograms used *
*      
*     S   ACTVAR          Finding the index set of free and active
*                           variables at the generalized Cauchy point.
*     S   AGBFGS          Simplified subgradient aggregation.
*     S   AGRSR1          Subgradient aggregation.
*     S   BFGSXV          Computation of the product H*V or -H*V, 
*                           where H is the inverse approximation of 
*                           the Hessian calculated by the L-BFGS 
*                           formula and V is an arbitrary vector.
*     S   COPY            Copying of a vector.
*     S   COPY2           Copying of two vectors.
*     S   DLBFGS          Computing the search direction by limited
*                           memory BFGS update.
*     S   DLSR1           Computing the search direction by limited
*                           memory SR1 update.    
*     S   DOBUN           Bundle selection.
*     S   LLS3            Line search using function values and
*                           derivatives.
*     S   XDIFFY          Difference of two vectors.
*     S   VNEG            Copying of a vector with change of the sign.
*     S   PROJGR          Simple projection of the subgradient and 
*                           calculation of Eucleidean norm.
*     S   RESTAB          Initialization of LMBM-B.
*     S   GETIME          Execution time.
*     S   RPRINT          Printout the results.
*     S   TERMIB          Calculation of stopping criterion and test 
*                           for termination with desired accuracy.
*     S   TINIT           Calculation of initial step size.
*     S   WPRINT          Printout the error and warning messages.
*     RF  EPS0            The smallest positive number such that
*                           1.0 + EPS0 > 1.0. 
*     RF  VDOT            Dot product of two vectors.
*     
*
*
*     * External subroutines *
*      
*     SE  FUNDER          Computation of the value and the subgradient 
*                         of the objective function. Calling sequence:
*                         CALL FUNDER(N,X,F,G,ITERM), where N is a 
*                         number of variables, X(N) is a vector of 
*                         variables, F is the value of the objective 
*                         function, G(N) is the subgradient of the 
*                         objective function, and ITERM is the error
*                         indicator.
*
*      

      
      SUBROUTINE LMBMB(N,NA,MC,MCU,NACT,X,XL,XU,ITYPE,IB,IACT,XCP,XO,S,
     &     G,GM,GA,PGA,PGM,U,D,F,AX,AG,AF,SM,UM,RM,LM,UMTUM,SMTSM,CDIAG,
     &     SMTGM,UMTGM,SMTPGM,UMTPGM,VMC1,VMC2,VMC3,VMC4,VMC5,VMC6,VMC7,
     &     VMC8,PGAH,VN1,VN2,TMPMAT,KMAT,TOLF,TOLF2,TOLB,TOLG,
     &     TOLG2,ETA,EPSL,XMAX,MIT,MFE,MOS,MTESF,IPRINT,ISCALE,ISTEP,
     &     NIT,NFE,ITERM,TIME,RTIM)

*     Scalar Arguments
      INTEGER N,NA,MC,MCU,MIT,MFE,MOS,MTESF,IPRINT,NIT,NFE,
     &     ITERM,ISCALE,ISTEP,ITYPE,NACT
      DOUBLE PRECISION F,ETA,EPSL,TOLF,TOLF2,TOLB,TOLG,TOLG2,XMAX

*     Array Arguments
      INTEGER IB(*),IACT(*)
      DOUBLE PRECISION X(*),XCP(*),XL(*),XU(*),XO(*),S(*),G(*),GM(*),
     &     GA(*),PGA(*),PGM(*),U(*),D(*),AX(*),AG(*),AF(*),SM(*),UM(*),
     &     RM(*),LM(*),UMTUM(*),SMTSM(*),CDIAG(*),SMTGM(*),UMTGM(*),
     &     SMTPGM(*),UMTPGM(*),
     &     PGAH(*),VMC1(*),VMC2(*),VMC3(*),VMC4(*),VMC5(*),VMC6(*),
     &     VMC7(*),VMC8(*),VN1(*),VN2(*),TMPMAT(*),KMAT(*)

*     Local Scalars
      INTEGER I,INEW,IOLD,IBFGS,ISR1,ITERS,MAL,MN,NNK,NTESF,NTESF2,
     &     NCRES1,NCRES2,NCRES3,NRES,ICN,NEPS,IFLAG,IBUN,NOUT,IERR,
     &     ITMAX,IPD,ISGNGA,NEPS2,IBADD
      DOUBLE PRECISION BETA,AGBETA,AMUGAD,XBX,EPSR,DNORM,PGANRM,WK,QK,
     &     P,TMAX,T,FO,GAMMA,PWK,THETA,TTHETA,SMALL,TMP,PGAHGA,PGAHG,
     &     PGNRM,PGMNRM,PREVF,ALPHA,DSTRN

*     External Functions
      DOUBLE PRECISION VDOT,EPS0
      EXTERNAL VDOT,EPS0

*     External Subroutines
      EXTERNAL FUNDER,COPY,XDIFFY,VNEG,DLBFGS,DLSR1,
     &     LLS3,AGBFGS,AGRSR1,DOBUN,TINIT,PROJGR,RESTAB,COPY2,
     &     RPRINT,WPRINT,ACTVAR,TERMIB,GETIME,BFGSXV,RWAXV2

*     Intrinsic Functions
      INTRINSIC ABS,MAX,SQRT

*     Computational Time
      REAL TIME,STRTIM,CTIM,RTIM(2)
       
*     Parameters
      INTEGER MAXEPS,MXEPS2,MAXNRS
      DOUBLE PRECISION ETA9,FMIN,TMIN,LENGTHD,RHO
      PARAMETER(
     &     MAXEPS = 20,
     &     MXEPS2 = 1,
     &     MAXNRS = 1000,
     &     ETA9 = 1.0D+60,
     &     FMIN = -1.0D+60,
     &     TMIN = 1.0D-12,
     &     LENGTHD = 1.0D+20,
     &     RHO = 1.0D-12)

      
      IF (IPRINT .GT. 3) WRITE (6,FMT='(1X,''Entry to LMBM-B:'')')


*     
*     Initialization
*

      NOUT   = 0
      NIT    = 1
      NFE    = 0
      NTESF  = 0
      NTESF2 = 0
      NCRES1 = 1
      NCRES2 = 0
      NCRES3 = 0
      NRES   = -1
      NEPS   = 0
      NEPS2  = 0
      IBADD  = 0
      ITERM  = 0
      ITERS  = 1
      NNK    = 0
      ISR1   = 0
      IPD    = 0
      BETA   = 0.0D+00
      AGBETA = 0.0D+00

      SMALL = EPS0()
      STRTIM = RTIM(1)

      IF (TOLF  .LE. 0.0D+00) TOLF = 1.0D-8
      IF (TOLF2 .EQ. 0.0D+00) TOLF2 = 1.0D+04
      IF (TOLB  .EQ. 0.0D+00) TOLB = FMIN + SMALL
      IF (TOLG  .LE. 0.0D+00) TOLG = 1.0D-05
      IF (TOLG2 .LE. 0.0D+00) TOLG2 = TOLG
      IF (XMAX  .LE. 0.0D+00) XMAX = 1.5D+00
      IF (ETA   .LT. 0.0D+00) ETA = 0.50D+00
      IF (EPSL  .LE. 0.0D+00) EPSL = 1.0D-04
      IF (MOS   .LE. 0) MOS = 2
      IF (MTESF .LE. 0) MTESF = 10
      IF (MIT   .LE. 0) MIT = 10000
      IF (MFE   .LE. 0) MFE = 20000
      
      TMAX = XMAX
      WK = ETA9
      PREVF = ETA9


      EPSR = 0.25D+00+SMALL
      IF (2.0D+00*EPSL .GE. EPSR) THEN
         EPSR = 2.0D+00*EPSL + SMALL
         IF (EPSR .GE. 0.5D+00) THEN
            CALL WPRINT(ITERM,IPRINT,-2)
         END IF
      END IF
            

*     
*     Computation of the value and the subgradient of the objective
*     function and the search direction for the first iteration
*
      
      CALL FUNDER(N,X,F,G,ITERM)
      NFE = NFE + 1
      
      IF (ITERM .NE. 0) GOTO 900

      GOTO 800
      

*     
*     Start of the iteration
*

 100  CONTINUE

      
*     
*     Direction finding
*

      IF (ITERS.GT.0) THEN
         

*     
*     Serious step initialization
*

         ICN = 0
         BETA = 0.0D+00
         AGBETA = 0.0D+00


*     
*     L-BFGS update
*

         CALL UPSERI(N,MC,MN,INEW,IOLD,IPD,IFLAG,IBFGS,G,GM,S,U,
     &        SM,UM,LM,RM,CDIAG,SMTSM,UMTUM,SMTGM,UMTGM,VMC1,VMC2,
     &        VMC3,VMC4,VMC5,VMC6,AMUGAD,XBX,TTHETA,GAMMA,
     &        SMALL,ISCALE,IERR)

         IF (IERR .NE. 0) THEN
            IF (ITYPE .EQ. 0) THEN
               IFLAG = 0

               CALL VNEG(N,G,D)
               CALL VNEG(N,G,PGAH)

               AMUGAD = VDOT(N,G,G)
               XBX = 0.0D+00
               P = -AMUGAD
               DNORM  = AMUGAD
               PGANRM = AMUGAD
               PGAHGA = AMUGAD
               PGMNRM = AMUGAD

               CALL TERMIB(MC,MCU,PGAHGA,PGANRM,WK,PWK,QK,AGBETA,TOLG,
     &              TOLG2,ISGNGA,ITERM,IFLAG)
               IF (ITERM .NE. 0) GOTO 900

               GOTO 200

            ELSE
               GOTO 800

            END IF

         END IF


*
*     Calculation of the vector PGAH = -H*PGM and the scalars
*     PGMNRM = PGANRM = PGM'*PGM, and  PGAHGA = PGM'*H*PGM, where H 
*     denotes the L-BFGS approximation of the inverse of the Hessian 
*     matrix and PGM denotes the simple projection of subgradient.
*

         CALL PROJGR(N,X,XL,XU,G,PGM,PGMNRM,ISGNGA,IB,SMALL)

         CALL RWAXV2(N,MN,SM,UM,PGM,PGM,SMTPGM,UMTPGM)

         CALL BFGSXV(N,MN,IOLD,PGAH,PGM,VN2,SM,UM,RM,UMTUM,CDIAG,
     &        SMTPGM,UMTPGM,VMC3,VMC4,VMC5,GAMMA,SMALL,IERR,1)
         IF (IERR .NE. 0) GOTO 800


         PGAHGA = - VDOT(N,PGM,PGAH)
         PGANRM = PGMNRM


*     
*     Test of positive definiteness
*

         IF (PGAHGA .LT. 0.0D+00) GOTO 800


*     
*     Tests for termination
*

         CALL TERMIB(MC,MCU,PGAHGA,PGANRM,WK,PWK,QK,AGBETA,TOLG,TOLG2,
     &        ISGNGA,ITERM,IFLAG)
         IF (ITERM .NE. 0) GOTO 900


*
*     Correction.
*      

         IF (PGAHGA .LT. RHO*PGANRM) THEN

            CALL WPRINT(ITERM,IPRINT,-7)

c     Note! If this happens we may lose the global convegence. We should
c     use correction - not restart. However, in practice the correction is
c     rarely needed.

            GOTO 800

         END IF

         
*     
*     Direction determination
*

         CALL DLBFGS(N,NACT,MN,IOLD,ITYPE,IB,IACT,X,XCP,XL,XU,G,D,
     &        SM,UM,LM,RM,CDIAG,SMTSM,UMTUM,SMTGM,UMTGM,VN1,VN2,
     &        KMAT,TMPMAT,VMC1,VMC2,VMC3,VMC4,VMC5,VMC6,VMC7,VMC8,
     &        AMUGAD,XBX,DSTRN,ALPHA,GAMMA,SMALL,IERR)

         IF (IERR .NE. 0) THEN
            IF (IERR .GT. -10) THEN
               CALL WPRINT(ITERM,IPRINT,-6)
            END IF
         
            GOTO 800

         END IF
         
         IF (DSTRN .LE. TOLG) GOTO 800


*
*     Computation of norms
*

         P = VDOT(N,G,D)
         IF (ALPHA .EQ. 1.0D+00) THEN
            DNORM = DSTRN
         ELSE
            DNORM = VDOT(N,D,D)
         END IF


      ELSE

         IF (IPD .NE. 0) THEN

*     Non-positive definite matrix: Restart.

            CALL WPRINT(ITERM,IPRINT,-6)
            GOTO 800
         END IF
            

*
*     Computation of projected subgradient PGA and the norm 
*     PGANRM = PGA'*PGA
*
            
         CALL PROJGR(N,X,XL,XU,GA,PGA,PGANRM,ISGNGA,IB,SMALL)


*     
*     L-SR1 update and calculation of the vector PGAH = -H*PGA and 
*     the scalar PGAHGA = PGA'*H*PGA, where H denotes the L-SR1 
*     approximation of the inverse of the Hessian matrix.
*

         CALL UPNULL(N,MC,MN,INEW,IOLD,IFLAG,ISR1,NNK,PGAH,PGA,PGM,
     &        GA,GM,S,U,SM,UM,LM,RM,CDIAG,SMTSM,UMTUM,VMC1,VMC2,VMC7,
     &        VMC8,SMTGM,UMTGM,SMTPGM,UMTPGM,VN1,VN2,VMC3,VMC4,
     &        VMC5,VMC6,TMPMAT,GAMMA,PGAHGA,AMUGAD,XBX,TTHETA,SMALL,
     &        IPRINT,IERR)

         IF (IERR .NE. 0) THEN
            IF (ITYPE .EQ. 0 .AND. IERR .EQ. -11) THEN
 
               IFLAG = 0

               CALL VNEG(N,GA,D)
               CALL VNEG(N,GA,PGAH)

               AMUGAD = VDOT(N,GA,GA)
               XBX = 0.0D+00
               P = -AMUGAD
               DNORM = AMUGAD
               PGANRM = AMUGAD
               PGAHGA = AMUGAD

               CALL TERMIB(MC,MCU,PGAHGA,PGANRM,WK,PWK,QK,AGBETA,TOLG,
     &              TOLG2,ISGNGA,ITERM,IFLAG)
               IF (ITERM .NE. 0) GOTO 900

               GOTO 200

            ELSE
               GOTO 800
               
            END IF
            
         END IF


*     
*     Test of positive definiteness
*

         IF (PGAHGA .LT. 0.0D+00) GOTO 800


*     
*     Tests for termination
*

         CALL TERMIB(MC,MCU,PGAHGA,PGANRM,WK,PWK,QK,AGBETA,TOLG,TOLG2,
     &        ISGNGA,ITERM,IFLAG)
         IF (ITERM .NE. 0) GOTO 900


*
*     Correction.
*      

         IF (PGAHGA .LT. RHO*PGANRM .OR. ICN .EQ. 1) THEN
               
c     Note! If this happens we may lose the global convegence. We should
c     use correction - not restart. However, in practice the correction is
c     rarely needed.

            CALL WPRINT(ITERM,IPRINT,-7)

            ICN = 1
            GOTO 800

         END IF


*     
*     Direction determination
*

         CALL DLSR1(N,MN,IOLD,X,XCP,XL,XU,GA,D,ITYPE,NACT,IB,
     &           IACT,IFLAG,SM,UM,LM,RM,CDIAG,SMTSM,UMTUM,VMC1,VMC2,
     &           VN1,VN2,S,VMC3,VMC4,VMC5,VMC6,TMPMAT,
     &           KMAT,GAMMA,AMUGAD,XBX,DSTRN,ALPHA,SMALL,IPRINT,IERR)

         IBFGS=0
            
         IF (IERR .NE. 0) THEN

            CALL WPRINT(ITERM,IPRINT,-6)
            GOTO 800

         END IF

            
         IF (DSTRN .LE. TOLG) GOTO 800


*
*     Computation of norms
*
            
         P = VDOT(N,GA,D)
         IF (ALPHA .EQ. 1.0D+00) THEN
            DNORM = DSTRN
         ELSE
            DNORM = VDOT(N,D,D)
         END IF


      END IF


      
 200  CONTINUE


*     
*     Test on descent direction
*

      IF (P+SMALL*SQRT(PGANRM*DNORM) .LE. 0.0D+00) THEN
         NCRES1 = 0

      ELSE
         NCRES1 = NCRES1 + 1

         IF (NCRES1 .EQ. 1) THEN

*     Nondescent search direction: Restart.

            CALL WPRINT(ITERM,IPRINT,-3)
            GOTO 800

         END IF

         NOUT = -1
         ITERM = -1
         GOTO 900
      END IF


*     
*     Tests for termination without desired accuracy
*

      IF (TIME .GT. 0.0E+00) THEN
         CALL GETIME(CTIM,RTIM)
         IF (CTIM-STRTIM .GT. TIME) THEN
            ITERM = 6
            GOTO 900
         END IF
      END IF


      IF (NFE .GE. MFE) THEN
         NOUT = MFE
         ITERM = 4
         GOTO 900
      END IF

      
      IF (NIT .GE. MIT) THEN
         NOUT = MIT
         ITERM = 5
         GOTO 900
      END IF

      
      IF (F .LE. TOLB) THEN
         ITERM = 7
         GOTO 900
      END IF


      IF (ITERS .EQ. 0) THEN
         IF (ABS(WK - PWK) .LE. SMALL) THEN
            NEPS = NEPS + 1

            IF (NEPS .GT. MAXEPS) THEN
               NEPS2 = NEPS2 + 1

               IF (NEPS2 .GT. MXEPS2) THEN
                  ITERM = -4
                  GOTO 900

               ELSE
                  GOTO 800

               END IF

            END IF

         ELSE
            NEPS = 0

         END IF

      ELSE
         NEPS = 0

      END IF


      IF (PWK .LT. WK .AND. NNK .GT. 2) THEN

*     Does not converge

         CALL WPRINT(ITERM,IPRINT,-4)

      END IF

      CALL RPRINT(N,NIT,NFE,NACT,X,F,WK,QK,ITERM,IPRINT)
      
      
*     
*     Preparation of line search
*

      FO = F
      
      IF (ITERS .GT. 0) THEN
         CALL COPY2(N,X,XO,G,GM)
      END IF

c      IF (DNORM .GT. 1.0D+5*SMALL) THEN
      IF (DNORM .GT. SMALL) THEN
         NCRES3 = 0

         TMAX = XMAX/SQRT(DNORM)
         IF (ALPHA .LT. 1.0D+00) THEN
            TMAX = MIN(TMAX,1.0D+00)
         END IF

      ELSE

         IF (WK .LE. TOLG) THEN
            ITERM = 1
            GOTO 900
         ELSE

            NCRES3 = NCRES3 + 1
            IF (NCRES3 .EQ. 1) GOTO 800

         END IF
         
         ITERM = -1
         GOTO 900
            
      END IF


*
*     Defining the maximum step length to the closest bound along the
*     current search direction.
*
      IF (TMAX .GT. 1.0D+00) THEN
         IF (ITYPE .EQ. 1) THEN

            DO 260 I = 1,N
               IF (IB(I) .NE. 0) THEN

*     Otherwise the lower bound is never met
                  IF (IB(I) .LE. 2 .AND. D(I) .LT. 0.0D+00) THEN
                     IF (XL(I) - X(I) .LT. -1.0D+01*SMALL .AND. 
     &                    D(I) .LE. -1.0D+01*SMALL) THEN
                        TMP = (XL(I) - X(I))/D(I)
                        IF (TMP .LT. TMAX) THEN
                           ITMAX=I
                           TMAX=TMP
                        END IF
                     ELSE
                        D(I) = 0.0D+00
                     END IF

c     Otherwise the upper bound is never met
                  ELSE IF (IB(I) .GE. 2 .AND. D(I) .GT. 0.0D+00) THEN 
                     IF (XU(I) - X(I) .GT. 1.0D+01*SMALL .AND. 
     &                    D(I) .GE. 1.0D+01*SMALL) THEN
                        TMP = (XU(I) - X(I))/D(I)
                        IF (TMP .LT. TMAX) THEN
                           ITMAX=I
                           TMAX=TMP
                        END IF
                     ELSE
                        D(I) = 0.0D+00
                     END IF
                     
                  END IF
               END IF
 260        CONTINUE

         END IF
      END IF

      
      IF (TMAX .GE. TMIN) THEN
         NCRES2 = 0

      ELSE
         NCRES2 = NCRES2 + 1
         IF (NCRES2 .EQ. 1) THEN

            CALL WPRINT(ITERM,IPRINT,-5)

            GOTO 800

         END IF

         ITERM = -1
         GOTO 900

      END IF
      
*
*     Initial step size
*

      IF (ISTEP .EQ. 0) THEN
         CALL TINIT(N,NA,MAL,X,AF,AG,AX,IBUN,D,F,P,T,TMAX,TMIN,
     &        ETA,ETA9,MOS,SMALL,ITERS)

      ELSE
         IF (WK .EQ. 0.0D+00 .AND. NA .NE. 0) THEN
            CALL TINIT(N,NA,MAL,X,AF,AG,AX,IBUN,D,F,P,T,TMAX,TMIN,
     &           ETA,ETA9,MOS,SMALL,ITERS)
            
         ELSE IF (ITERS .EQ. 1) THEN
            T = MIN(2.0D+00,TMAX)
         ELSE
            T = MIN(1.0D+00,TMAX)
         END IF

      END IF

      
*     
*     Line search with directional derivatives which allows null steps
*

      THETA=1.0D+00
      IF (SQRT(DNORM) .GT. LENGTHD) THEN
         THETA=LENGTHD/SQRT(DNORM)
      END IF


      CALL LLS3(N,X,XO,XL,XU,IB,G,VN2,PGAH,PGAHG,PGNRM,D,T,F,FO,BETA,
     &     TMIN,SQRT(DNORM),WK,THETA,TTHETA,EPSL,EPSR,ETA,MOS,
     &     ITERS,NFE,NNK,SMALL,ITERM)

      IF (ITERM .EQ. -1) THEN


*     Deficient search direction.

         IBADD = IBADD + 1
         CALL COPY(N,XO,X)
         F = FO
         
         IF (IBADD .EQ. 1) THEN
            ITERM = 0
            CALL WPRINT(ITERM,IPRINT,-8)
            GOTO 800

         ELSE

            GOTO 900
         END IF

      ELSE

         IBADD = 0

      END IF

      IF (ITERM .NE. 0) GOTO 900

      IF (TOLF2 .GE. 0) THEN
         IF (ABS(FO-F) .LE. TOLF2*SMALL*MAX(ABS(F),ABS(FO),1.0D+00)
     &        .AND. ITERS .EQ. 1) THEN
         
            ITERM = 3
            GOTO 900
         
         END IF
      END IF

      IF (ABS(FO-F) .LE. TOLF) THEN
         NTESF = NTESF + 1
         
         IF (NTESF .GE. MTESF .AND. ITERS .EQ. 1) THEN
            ITERM = 2
            GOTO 900
         END IF
         
      ELSE
         NTESF = 0
      END IF

  
*
*     Bundle updating
*      

      IF (NA .NE. 0) THEN
         CALL DOBUN(N,NA,MAL,X,G,F,AX,AG,AF,ITERS,IBUN)
      END IF
      

*
*     Computation of variables difference 
*

      CALL XDIFFY(N,X,XO,S)

     
*
*     Computation of aggregate values and subgradients difference
*

      IF (ITERS.EQ.0) THEN
         NNK = NNK + 1

         IF (NNK.EQ.1) THEN
            CALL XDIFFY(N,G,GM,U)
            CALL AGBFGS(N,MN,IOLD,PGAHG,PGAHGA,PGNRM,G,GM,GA,VN2,
     &           SM,UM,RM,CDIAG,UMTUM,GAMMA,BETA,AGBETA,VMC1,VMC2,
     &           SMALL)

         ELSE
            
            CALL AGRSR1(N,MN,IOLD,G,GM,GA,VN2,PGM,PGA,PGAH,PGAHG,PGAHGA,
     &           PGNRM,PGMNRM,BETA,AGBETA,SM,UM,TMPMAT,UMTUM,RM,
     &           SMTPGM,UMTPGM,VN1,X,VMC1,VMC2,GAMMA,SMALL,IERR)

            IF(IERR .NE. 0) THEN
               CALL WPRINT(ITERM,IPRINT,-6)
               
            END IF
            CALL XDIFFY(N,G,GM,U)
         END IF

         CALL COPY(N,XO,X)
         F = FO
         
      ELSE
         NNK = 0
         CALL XDIFFY(N,G,GM,U)
      END IF

      NIT = NIT + 1
      GOTO 100


 800  CONTINUE

*
*     Restart
*

      IF (PREVF .GT. F) THEN
         PREVF = F
         NTESF2 = 0
         
      ELSE
         NTESF2 = NTESF2 + 1
         IF(NTESF2  .GT. MTESF) THEN
            ITERM = 2
            GOTO 900
         END IF
      END IF


      CALL RESTAB(N,MC,MCU,MN,INEW,NA,MAL,IBUN,ISTEP,ITERS,NNK,
     &     ITYPE,IB,IACT,NACT,IBFGS,ICN,IFLAG,IPD,F,X,XCP,
     &     XL,XU,D,GM,G,PGM,PGMNRM,VN2,AX,AG,AF,AGBETA,BETA,
     &     GAMMA,AMUGAD,XBX,ALPHA,WK,PWK,QK,TOLG,TOLG2,SMALL,
     &     NRES,MAXNRS,NOUT,IERR,ITERM)

      IF (ITERM .NE. 0) GOTO 900

      P = VDOT(N,G,D)
      DSTRN = AMUGAD

      IF (DSTRN .LE. TOLG*TOLG) THEN
         ITERM = 1
         GOTO 900
      END IF
      
      IF (ALPHA .EQ. 1.0D+00) THEN
         DNORM = AMUGAD
      ELSE
         DNORM = VDOT(N,D,D)
      END IF
      
      CALL VNEG(N,PGM,PGAH)
      PGAHGA = PGMNRM
      PGANRM = PGMNRM
      
      GOTO 200


 900  CONTINUE
      

*
*     Printout the final results
*

      IF (IPRINT .GT. 3) WRITE (6,FMT='(1X,''Exit from LMBM-B:'')')

      CALL WPRINT(ITERM,IPRINT,NOUT)
      CALL RPRINT(N,NIT,NFE,NACT,X,F,WK,QK,ITERM,IPRINT)

      RETURN
      END

************************************************************************
      
