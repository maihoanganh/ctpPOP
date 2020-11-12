*************************************************************************
*
*     * LMBMBS includes the following subroutines *
*
*     S   ACTVAR          Finding the index set of free and active
*                           variables at the generalized Cauchy point.
*     S   AGBFGS          Simplified subgradient aggregation with L-BFGS
*                           formula.
*     S   AGRSR1          Subgradient aggregation.
*     S   BFGSXV          Computation of the product H*V or -H*V, 
*                           where H is the inverse approximation of 
*                           the Hessian calculated by the L-BFGS 
*                           formula and V is an arbitrary vector.
*     S   BMXV2           Computation of the product of the 2MC x 2MC
*                           middle matrix in the compact L-BFGS formula
*                           and a 2MC vector.
*     S   CP1ST           Computing the generalized Cauchy point and the
*                           search direction at the first iteration.
*     S   CPBFGS          Computing the generalized Cauchy point by the
*                           limited memory BFGS update formula.
*     S   CPSR1           Computing the generalized Cauchy point by the
*                           limited memory SR1 update formula.
*     S   DESTEP          Stepsize determination for descent steps.
*     S   DLBFGS          Computing the search direction by the limited 
*                           memory BFGS update.
*     S   DLSR1           Computing the search direction by the limited 
*                           memory SR1 update.
*     S   DOBUN           Bundle selection.
*     S   FRMJ2           Formation of the JJ' Cholesky factorization of
*                           matrix J = SM' SM / GAMMA + LC^(-1)L'.
*     S   FRMLEL          Formation of the LEL' factorization of the
*                           indefinite (BFGS) matrix.
*     S   GETIME          Execution time.
*     S   INDIC2          Initialization of indices.
*     S   KMXV2           Computation of the product of the inverse of
*                           2MC x 2MC matrix and and a 2MC vector.
*     S   LLS3            Special line search for LMBM-B.
*     S   NULSTP          Stepsize determination for null steps.
*     S   PROJGR          Simple projection of the subgradient and 
*                           calculation of Eucleidean norm.
*     S   PROJX           Projection of the initial X to the feasible
*                           region if necessary.
*     S   QINT            Quadratic interpolation for line search.
*     S   RESTAB          Initialization of LMBM-B.
*     RF  SCLPAR          Calculation of the scaling parameter.
*     S   SR1XV           Computation of the product H*V or -H*V, 
*                           where H is the inverse approximation of 
*                           the Hessian calculated by the L-SR1 
*                           formula and V is an arbitrary vector.
*     S   SSBFGS          Subspace minimization and direction finding
*                           using limited memory BFGS matrices.
*     S   SSSR1           Subspace minimization and direction finding
*                           using limited memory SR1 matrices.
*     S   TERMIB          Calculation of stopping criterion and test 
*                           for termination with desired accuracy.
*     S   TINIT           Calculation of the initial step size.
*     S   RPRINT          Printout the results.
*     S   UPDATE          Updating the matrices RM, LM, CDIAG, SMTSM,
*                           and UMTUM for LMBM-B.
*     S   UPNULL          Updating limited memory SR1 matrices.
*     S   UPSERI          Updating limited memory BFGS matrices.
*     S   WPRINT          Printout error and warning messages.
*      
*
*      
*     Napsu Karmitsa (maiden name Haarala) 2002 - 2007, 
*     last modified 2009. 
*
*
************************************************************************
*
*     * SUBROUTINE PROJX *
*
*      
*     * Purpose *
*      
*     Projection of the initial X to the feasible region if necessary.
*
*      
*     * Calling sequence *
*     
*     CALL PROJX(N,X,XL,XU,IB,NACT,ITYPE)
*
*      
*     * Parameters *
*     
*     II  N               Number of variables.
*     RU  X(N)            Vector of variables. 
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables XU > XL.
*     II  IB(N)           Type of bound constraints:
*                           0  - X(I) is unbounded,
*                           1  - X(I) has only a lower bound,
*                           2  - X(I) has both lower and upper bounds,
*                           3  - X(I) has only an upper bound. 
*     IO  NACT            Number of active variables.
*     IO  ITYPE           Type of problem:
*                           0  - problem is unbounded,
*                           1  - problem is bounded.
*     
*      
*     Napsu Karmitsa (2005, last modified 2009)
*
*      
      
      SUBROUTINE PROJX(N,X,XL,XU,IB,NACT,ITYPE)

      
*     Scalar Arguments
      INTEGER N,NACT,ITYPE
      
*     Array Arguments
      INTEGER IB(*)
      DOUBLE PRECISION X(*),XL(*),XU(*)

*     Local Scalars
      INTEGER I


*
*     Project the initial x to the feasible set.
*

      NACT = 0
      ITYPE = 0
      
      DO 10 I=1,N
         IF (IB(I) .GT. 0) THEN
            ITYPE = 1

            IF (IB(I) .LE. 2) THEN
               IF (X(I) .LE. XL(I)) THEN
                  X(I)=XL(I)
                  NACT = NACT+1
               END IF
            END IF

            IF (IB(I) .GE. 2) THEN
               IF (X(I) .GE. XU(I)) THEN
                  X(I)=XU(I)
                  NACT = NACT+1
               END IF
            END IF

         END IF
 10   CONTINUE


      RETURN
      END



*************************************************************************
*
*     * SUBROUTINE RESTAB *
*
*      
*     * Purpose *
*      
*     Initialization of LMBM-B.
*
*      
*     * Calling sequence *
*     
*      CALL RESTAB(N,MC,MCU,MN,INEW,NA,MAL,IBUN,ISTEP,ITERS,NNK,ITYPE,
*    &     IB,IACT,NACT,IBFGS,ICN,IFLAG,IPD,F,X,XCP,XL,XU,D,GM,G,PGA,
*    &     PGANRM,VN,AX,AG,AF,AGBETA,BETA,GAMMA,AMUGAD,XBX,WK,PWK,QK,
*    &     TOLG,TOLG2,SMALL,NRES,MAXNRS,NOUT,IERR,ITERM)
*     
*     
*     * Parameters *
*      
*     II  N               Number of variables.
*     II  MCU             Upper limit for maximum number of stored
*                           corrections, MCU >= 3.
*     IU  MC              Current maximum number of stored corrections, 
*                           MC <= MCU.
*     IO  MN              Current number of depositories used.
*     IO  INEW            Index for the circular arrays.
*     II  NA              Maximum bundle dimension.
*     IO  MAL             Current size of the bundle.
*     IO  IBUN            Index for the circular arrays in bundle
*                           updating.
*     II  ISTEP           Selection of initial stepsize before line 
*                         search procedure:
*                           0  - Stepsize selection using polyhedral 
*                                approximation.
*                           1  - Stepsize = min(1.0,TMAX), where TMAX is
*                                the upper limit for step size assuring 
*                                the feasibility of produced point (no 
*                                additional bundle or bundle updating is 
*                                needed).
*     IU  ITERS           Null step indicator.
*                           0  - Null step.
*                           1  - Serious step.
*     IO  NNK             Consecutive null steps counter.
*     II  ITYPE           Type of problem:
*                           0  - problem is unbounded,
*                           1  - constrained problem.
*     II  IB(N)           Type of bound constraints:
*                           0  - X(I) is unbounded,
*                           1  - X(I) has only a lower bound,
*                           2  - X(I) has both lower and upper bounds,
*                           3  - X(I) has only an upper bound. 
*     IO  IACT(N)         Index set of active and free variables.
*     IO  NACT            Number of active variables.
*     IO  IBFGS           Index of the type of BFGS update.
*     IO  ICN             Correction indicator for null steps.
*     IO  IFLAG           Index for adaptive version.
*     IO  IPD             Index for positive definiteness for SR1 update:
*                           0  - Positive definiteness is preserved.
*                          -1  - Positive definiteness is not preserved.
*     RI  F               Value of the objective function.
*     RI  X(N)            Vector of variables.
*     RO  XCP(N)          Generalized Cauchy point.
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables. 
*     RO  D(N)            Search direction.
*     RI  GM(N)           Basic subgradient of the objective function.
*     RO  G(N)            Current subgradient of the objective function.
*     RO  PGM(N)          Projected subgradient of the objective function.
*     RO  PGANRM          Euclidean norm of the projected subgradient 
*                           vector.
*     RA  VN(N)           Auxiliary array.
*     RO  AX(N*NA)        Matrix whose columns are bundle points.
*     RO  AG(N*NA)        Matrix whose columns are bundle subgradients.
*     RO  AF(3*NA)        Vector of bundle values.
*     RO  BETA            Locality measure.
*     RO  AGBETA          Aggregate locality measure.
*     RO  GAMMA           Scaling parameter.
*     RO  AMUGAD          AMUGAD = D'*D.
*     RO  XBX             XBX = (XCP-X)'*(XCP-X)
*     RO  ALPHA           Backtrack multiplier.
*     RU  WK              Stopping criterion.
*     RO  PWK             Previous stopping criterion.
*     RO  QK              Second stopping criterion.
*     RI  TOLG            Tolerance for the first termination criterion.
*     RI  TOLG2           Tolerance for the second termination criterion.
*     RI  SMALL           The smallest positive number such that
*                           1.0 + SMALL > 1.0.
*     IO  NRES            Number of restarts.
*     II  MAXNRS          Maximum number of restarts.
*     II  NOUT            Auxilary printout specification.
*     IO  IERR            Error indicator.
*     IO  ITERM           Cause of termination:
*                           0  - Everything is ok.
*                          -2  - Number of restarts > maximum number
*                                of restarts.
*
*
*     
*     * Local Parameters *
*
*     I   ISGNGA          ISGNGA = 0, if G(I)<=0 for all I such that 
*                           X(I)=XU(I) and G(I)>=0 for all I such that
*                           X(I)=XL(I). Otherwise ISGNGA = 1.
*
*
*      
*     * Subprograms used *
*      
*     S   ACTVAR          Finding the index set of free and active
*                           variables at the generalized Cauchy point.
*     S   COPY            Copying of a vector.
*     S   CP1ST           Computing the generalized Cauchy point and the 
*                           search direction at the first iteration and 
*                           in case of restart.
*     S   DOBUN           Bundle updating.
*     S   PROJGR          Simple projection of the subgradient and 
*                           calculation of Eucleidean norm.
*     S   TERMIB          Calculation of stopping criterion and test 
*                           for termination with desired accuracy.
*     S   VNEG            Copying of a vector with change of the sign.
*
*     RF  VDOT            Dot product of two vectors.
*     
*      

      SUBROUTINE RESTAB(N,MC,MCU,MN,INEW,NA,MAL,IBUN,ISTEP,ITERS,NNK,
     $     ITYPE,IB,IACT,NACT,IBFGS,ICN,IFLAG,IPD,F,X,XCP,XL,XU,D,GM,G,
     $     PGM,PGANRM,VN,AX,AG,AF,AGBETA,BETA,GAMMA,AMUGAD,XBX,ALPHA,
     $     WK,PWK,QK,TOLG,TOLG2,SMALL,NRES,MAXNRS,NOUT,IERR,ITERM)

*     Scalar Arguments
      INTEGER N,MC,MCU,MN,INEW,NA,MAL,IBUN,ISTEP,ITERS,NNK,ITYPE,NACT,
     &     IBFGS,ICN,IFLAG,IPD,NRES,MAXNRS,NOUT,IERR,ITERM
      DOUBLE PRECISION F,PGANRM,BETA,AGBETA,GAMMA,AMUGAD,XBX,SMALL,
     &     WK,PWK,QK,TOLG,TOLG2,ALPHA

*     Array Arguments
      INTEGER IB(*),IACT(*)
      DOUBLE PRECISION X(*),XCP(*),XL(*),XU(*),D(*),G(*),GM(*),PGM(*),
     &     AX(*),AG(*),AF(*),VN(*)

*     Local scalars
      INTEGER ISGNGA

*     External Functions
      DOUBLE PRECISION VDOT
      EXTERNAL VDOT

*     External Subroutines
      EXTERNAL ACTVAR,COPY,CP1ST,DOBUN,VNEG,PROJGR,TERMIB


      MN    = 0
      INEW  = 1
      IBUN  = 1
      IBFGS = 0
      ICN   = 0
      IFLAG = 0
      IPD   = 0
      IERR  = 0
      MAL   = 0
      NRES  = NRES + 1


      IF (NRES .GT. MAXNRS) THEN
         NOUT = MAXNRS
         ITERM = -2
         GOTO 900

      END IF

      IF (ITERS .EQ. 0) THEN

         CALL COPY(N,GM,G)
         ITERS = 1
         NNK = 0
         AGBETA=0.0D+00
         BETA=0.0D+00

      END IF

      GAMMA=1.0D+00

      CALL PROJGR(N,X,XL,XU,G,PGM,PGANRM,ISGNGA,IB,SMALL)


*
*     Termination with the desired accuracy
*

      CALL TERMIB(MC,MCU,PGANRM,PGANRM,WK,PWK,QK,AGBETA,TOLG,TOLG2,
     &     ISGNGA,ITERM,IFLAG)


      IF (ITERM .NE. 0) GOTO 900


*     No need for correction, since PGANRM > RHO*PGANRM.

      IF (ITYPE .EQ. 0) THEN
         CALL VNEG(N,G,D)
         AMUGAD = VDOT(N,G,G)
         XBX = 0.0D+00

      ELSE
         CALL CP1ST(N,X,XCP,XL,XU,IB,IACT,G,D,VN,AMUGAD,XBX,ALPHA,SMALL)
         CALL ACTVAR(N,NACT,IACT,IB,XCP,XL,XU,SMALL)
               
      END IF

      IF (NA .NE. 0) THEN
         CALL DOBUN(N,NA,MAL,X,G,F,AX,AG,AF,ITERS,IBUN)
      END IF


 900  CONTINUE

      RETURN
      END


*************************************************************************
*
*     * SUBROUTINE CP1ST *
*
*      
*     * Purpose *
*      
*     Computation of generalized Cauchy point and the search direction 
*     at the first iteration.
*
*      
*     * Calling sequence *
*     
*     CALL CP1ST(N,X,XCP,XL,XU,IB,IT,G,DC,T,AMUGAD,XBX,ALPHA,SMALL)
*
*      
*     * Parameters *
*     
*     II  N               Number of variables.
*     RI  X(N)            Vector of variables. 
*     RO  XCP(N)          Generalized Cauchy point.
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables. 
*     II  IB(N)           Type of bound constraints:
*                           0  - X(I) is unbounded,
*                           1  - X(I) has only a lower bound,
*                           2  - X(I) has both lower and upper bounds,
*                           3  - X(I) has only an upper bound. 
*     IA  IT(N)           Index set of breakpoints:
*                           for I=1,...,NLEFT, IT(I) are the indices
*                             of breakpoints that have not been
*                             encountered,
*                           for i=NLEFT+1,...,NBREAK, IT(I) are the
*                             indices of encountered breakpoints and
*                           for i=IFREE,...,N, IT(I) are the indices
*                             of variables that have no bound
*                             constraints along the search direction.
*     RI  G(N)            Current subgradient of the objective function.
*     RO  DC(N)           Search direction.
*     RA  T(N)            Auxiliary array used to store the breakpoints.
*     RO  AMUGAD          AMUGAD = DC'*DC.
*     RO  XBX             XBX = (XCP-X)'*(XCP-X)
*     RO  ALPHA           Backtrack multiplier.
*     RI  SMALL           The smallest positive number such that
*                           1.0 + SMALL > 1.0.
*     
*
*     * Subprograms used *
*      
*     S   SCSUM           Sum of a vector and a scaled vector.
*     S   HPSRT           Heapsort algorithm.
*     S   XDIFFY          Difference of two vectors.
*
*     RF  VDOT            Dot product of two vectors.
*
*
      
      SUBROUTINE CP1ST(N,X,XCP,XL,XU,IB,IT,G,DC,T,AMUGAD,XBX,ALPHA,
     $     SMALL)


*     Scalar Arguments
      INTEGER N
      DOUBLE PRECISION AMUGAD,XBX,SMALL,ALPHA
      
*     Array Arguments
      INTEGER IB(*),IT(*)
      DOUBLE PRECISION G(*),X(*),XCP(*),XL(*),XU(*),DC(*),T(*)
      
*     Local Scalars
      INTEGER I,J,NBREAK,IFREE,ITMIN,NINT,NLEFT,ISBP,IBOUND
      DOUBLE PRECISION F1,F2,TMIN,DLTMIN,TOLD,DELTAT,GISBP,GISBP2,
     &     ZISBP,DFREE,SMALL2

*     Intrinsic Functions
      INTRINSIC MAX,MIN,ABS,SQRT

*     External Functions
      DOUBLE PRECISION VDOT
      EXTERNAL VDOT

*     External Subroutines
      EXTERNAL HPSRT,SCSUM,XDIFFY


*
*     Initialization
*

      NBREAK = 0
      NLEFT  = 0
      IFREE  = N+1
      ITMIN  = 0
      IBOUND = 0   
      F1     = 0.0D+00
      F2     = 0.0D+00
      TMIN   = 0.0D+00
      XBX    = 0.0D+00

      SMALL2 = SQRT(SMALL)



*
*     Computation of the Cauchy direction DC, the first derivative
*     F1= -DC'*DC, and the breakpoints T in each coordinate
*     direction. Identification of the smallest breakpoint TMIN. 
*      

      DO 10 I=1,N

         IF (IB(I) .NE. 0) THEN
            IF (IB(I).LE.2 .AND. ABS(G(I)) .LE. SMALL2 .AND. 
     &           X(I)-XL(I) .LE. SMALL2) THEN
               DC(I) = 0.0D+00
               XCP(I) = X(I)
               GOTO 10
            END IF

            IF (IB(I).GE.2 .AND. ABS(G(I)) .LE. SMALL2 .AND. 
     &           XU(I)-X(I) .LE. SMALL2) THEN
               DC(I) = 0.0D+00
               XCP(I) = X(I)
               GOTO 10
            END IF
         END IF


         IF (IB(I).NE.0 .AND. IB(I).LE.2 .AND. G(I) .GT. 0.0D+00) THEN

              
*
*     Breakpoint in lower bound
*

            T(NBREAK+1) = (X(I) - XL(I))/G(I)
            
            IF (T(NBREAK+1) .LE. SMALL2) THEN
               DC(I) = 0.0D+00

            ELSE
               NBREAK = NBREAK + 1
               IT(NBREAK) = I
               DC(I) = -G(I)
               F1 = F1 - DC(I)*DC(I)
               

*     
*     Determination of the smallest breakpoint
*     
          
               IF (NBREAK .EQ. 1 .OR. T(NBREAK) .LT. TMIN) THEN
                  TMIN=T(NBREAK)
                  ITMIN=NBREAK
               END IF
            END IF
      
         ELSE IF (IB(I) .GE. 2 .AND. G(I) .LT. 0.0D+00) THEN


*
*     Breakpoint in upper bound
*               
 
            T(NBREAK+1) = (X(I) - XU(I))/G(I)
            
            IF (T(NBREAK+1) .LE. SMALL2) THEN
               DC(I) = 0.0D+00

            ELSE
               NBREAK = NBREAK + 1
               IT(NBREAK) = I
               DC(I) = -G(I)
               F1 = F1 - DC(I)*DC(I)


*     
*     Determination of the smallest breakpoint
*               

               IF (NBREAK .EQ. 1 .OR. T(NBREAK) .LT. TMIN) THEN
                  TMIN=T(NBREAK)
                  ITMIN=NBREAK
               END IF
            END IF
            
         ELSE
            

*
*     No breakpoint 
*

            IFREE = IFREE - 1
            IT(IFREE) = I
            DC(I) = -G(I)
            F1 = F1 - DC(I)*DC(I)
            IF(ABS(G(I)) .LE. 0.0D+00) IBOUND=-1
         END IF


*     
*     Initialization of Cauchy point
*

         XCP(I) = X(I)
           
 10   CONTINUE


*
*     The indices corresponding to free variables are located in
*     IT(1),...,IT(NBREAK) and IT(IFREE),...,IT(N). The smallest of the
*     NBREAK breakpoints is in T(ITMIN)=TMIN.
*

      
*     DC is a zero vector. Return with the initial XCP.

      NLEFT = NBREAK

      IF (NBREAK .EQ. 0 .AND. IBOUND .EQ. 0) GOTO 999


*
*     Calculation of derivative F2 = DC'*DC, 
*

      F2 =  -F1
      DLTMIN = 1.0D+00

      TOLD = 0.0D+00
      NINT = 1
      

*     There is no breakpoints, locate the GCP and return. 
      
      IF (NBREAK .EQ. 0) GOTO 888

     
*
*     Begining of the loop.
*      

 100  CONTINUE


      IF (NINT .EQ. 1) THEN


*     the smallest breakpoint is TMIN=T(ITMIN) and its index is
*     in IT(ITMIN)

         ISBP = IT(ITMIN)

      ELSE


*     the smallest breakpoint is chosen by heapsort algorithm

         IF (NINT .EQ. 2) THEN


*     Remove the already used smallest breakpoint before heapsort call.
            
            IF (ITMIN .NE. NBREAK) THEN
               T(ITMIN) = T(NBREAK)
               IT(ITMIN) = IT(NBREAK)
            END IF


*     
*     Heapsort with initialization of the heap.
*

            CALL HPSRT(NLEFT,T,IT,0)

         ELSE


*            
*     Heapsort with updating the heap.
*            

            CALL HPSRT(NLEFT,T,IT,1)

         END IF


         TMIN = T(NLEFT)
         ISBP = IT(NLEFT)

      END IF


      DELTAT = TMIN - TOLD
      
      IF (DELTAT .LT. 0.0D+00) THEN
         PRINT*,' no nyt voit jo hirttää ihtes.'
      END IF


*
*     Minimizer is within this interval, locate the GCP and return.
*

      IF (DLTMIN .LT. DELTAT) GOTO 888


*      
*     Examination of subsequent seqment
*

      NLEFT = NLEFT - 1
      NINT = NINT + 1
      GISBP = G(ISBP)
      DC(ISBP) = 0.0D+00

      IF (GISBP .LT. 0.0D+00) THEN
         XCP(ISBP) = XU(ISBP)
         ZISBP = XU(ISBP) - X(ISBP)
      ELSE 
         XCP(ISBP) = XL(ISBP)
         ZISBP = XL(ISBP) - X(ISBP)
      END IF


*         
*     All  variables are fixed, return with current XCP as GCP.
*         

      IF (NLEFT .EQ. 0 .AND. NBREAK .EQ. N) GOTO 999


*      
*     Update F1 and F2
*      

      GISBP2 = GISBP*GISBP

      F1 = F1 + DELTAT*F2 + GISBP2 + GISBP*ZISBP
      F2 = F2 - GISBP2

      XBX = XBX + TMIN*TMIN * GISBP2


*
*     DC is a zero vector. Return with the current XCP.
*

      IF (F2 .LE. SMALL2) GOTO 999

      DLTMIN = -F1/F2
      TOLD = TMIN



*
*     Repeat the loop for unsearched intervals.          
*

      IF (NLEFT .GT. 0) GOTO 100

      IF (IBOUND .EQ. 0) DLTMIN = 0.0D+00   


*     
*     End of the loop.
*

 888  CONTINUE

      DLTMIN = MAX(0.0D+00,DLTMIN)
      TOLD = TOLD + DLTMIN


*      
*     Update free variables 
* 

      DO 20 I=1,NLEFT
         XCP(IT(I))=X(IT(I))+TOLD*DC(IT(I))
 20   CONTINUE

      DO 30 I=IFREE,N
         XCP(IT(I))=X(IT(I))+TOLD*DC(IT(I))
 30   CONTINUE

      XBX = XBX + TOLD*TOLD * F2

     
*
*     Calculate the search direction for the first iteration.
*

 999  CONTINUE

  
*     The indices corresponding to free variables are located in
*     IT(1),...,IT(NLEFT) and IT(IFREE),...,IT(N). The 
*     considered breakpoints are in IT(NLEFT+1),...,IT(NBREAK).
*     Otherwise, DC(I) = 0.

      DO 40 I=NLEFT+1,NBREAK
         DC(IT(I))=XCP(IT(I))-X(IT(I))
 40   CONTINUE

      AMUGAD = VDOT(N,DC,DC)


*     
*     Backtrack to the feasible region if X + DC violates bounds
*     (among free variables).
*

      ALPHA = 1.0D+00
      DO 50 I=1,NLEFT
         J=IT(I)
         DFREE=DC(J)
         IF (IB(J) .NE. 0) THEN
            IF (IB(J) .LE. 2 .AND. X(J) + DFREE .LT. XL(J)) THEN

*     Otherwise the lower bound is never met
               IF (ABS(X(J)-XCP(J)+DFREE) .GT. SMALL) THEN
                  ALPHA=MIN(ALPHA,(XL(J)-XCP(J))/(X(J)-XCP(J)+DFREE))
               END IF
 
            ELSE IF (IB(J) .GE. 2 .AND. X(J) + DFREE  .GT. XU(J)) THEN

*     Otherwise the upper bound is never met
               IF (ABS(X(J)-XCP(J)+DFREE) .GT. SMALL) THEN
                  ALPHA=MIN(ALPHA,(XU(J)-XCP(J))/(X(J)-XCP(J)+DFREE))
               END IF

            END IF
         END IF
 50   CONTINUE

      IF (ALPHA .LT. 1.0D+00) THEN
         
         DO 70 I=1,NLEFT
            J=IT(I)
            DC(J) = XCP(J)-X(J)+ALPHA*(X(J)-XCP(J)+DC(J))
 70      CONTINUE

         DO 80 I=IFREE,N
            J=IT(I)
            DC(J) = XCP(J)-X(J)+ALPHA*(X(J)-XCP(J)+DC(J))
 80      CONTINUE

      ELSE

         XBX = 0.0D+00
            
      END IF

      RETURN
      END


************************************************************************
*
*     * SUBROUTINE ACTVAR *
*
*      
*     * Purpose *
*      
*     Finding the index set of free and active variables at the
*     generalized Cauchy point.
*
*      
*     * Calling sequence *
*     
*     CALL ACTVAR(N,NACT,IACT,IB,XCP,XL,XU,SMALL)
*
*      
*     * Parameters *
*     
*     II  N               Number of variables.
*     IO  NACT            Number of active variables.
*     IO  IACT(N)         Index set of active and free variables:
*                           for I=1,...,NACT, IACT(I) are the indices
*                             of active variables (-J, if XCP(J)=XL(J)
*                             and +J, otherwise)
*                           for i=NACT+1,...,N, IACT(I) are the indices
*                             of free variables.
*     II  IB(N)           Type of bound constraints:
*                           0  - X(I) is unbounded,
*                           1  - X(I) has only a lower bound,
*                           2  - X(I) has both lower and upper bounds,
*                           3  - X(I) has only an upper bound. 
*     RI  XCP(N)          Generalized Cauchy point.
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables. 
*     RI  SMALL           Small positive value.
*     
*

      SUBROUTINE ACTVAR(N,NACT,IACT,IB,XCP,XL,XU,SMALL)
      
*     Scalar Arguments
      INTEGER N,NACT
      DOUBLE PRECISION SMALL
      
*     Array Arguments
      INTEGER IB(*),IACT(*)
      DOUBLE PRECISION XCP(*),XL(*),XU(*)

*     Local Scalars
      INTEGER I,IFREE
      DOUBLE PRECISION SMALL2

      SMALL2 = 1.0D+01*SMALL
c      SMALL2 = 1.0D+03*SMALL
c      SMALL2 = SQRT(SMALL)

      NACT = 0
      IFREE = N+1

      DO 10 I = 1,N

         IF (IB(I) .EQ. 0) THEN
            IFREE = IFREE - 1
            IACT(IFREE) = I

         ELSE IF (IB(I) .LE. 2 .AND. XCP(I) .LE. XL(I)+SMALL2) THEN
            NACT = NACT + 1
            IACT(NACT) = -I

         ELSE IF (IB(I) .GE. 2 .AND. XCP(I) .GE. XU(I)-SMALL2) THEN
            NACT = NACT + 1
            IACT(NACT) = I

         ELSE
            IFREE = IFREE - 1
            IACT(IFREE) = I

         END IF

 10   CONTINUE
 
      RETURN
      END


************************************************************************
*
*     * SUBROUTINE PROJGR *
*
*      
*     * Purpose *
*      
*     Simple projection of the subgradient and calculation of 
*     Eucleidean norm.
*
*      
*     * Calling sequence *
*     
*     CALL PROJGR(N,X,XL,XU,GA,PGA,PGANRM,ISGNGA,IB,SMALL)
*
*      
*     * Parameters *
*
*     II  N               Number of variables.
*     RI  X(N)            Vector of variables. 
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables.
*     RI  GA(N)           Current subgradient of the objective function.
*     RO  PGA(N)          Projected subgradient.
*     RO  PGANRM          Eucleidean norm of projected subgradient.
*     IO  ISGNGA          ISGNGA = 0, if GA(I)<=0 for all I such that 
*                           X(I)=XU(I) and GA(I)>=0 for all I such that
*                           X(I)=XL(I). Otherwise ISGNGA = 1.
*     II  IB(N)           Type of bound constraints:
*                           0  - X(I) is unbounded,
*                           1  - X(I) has only a lower bound,
*                           2  - X(I) has both lower and upper bounds,
*                           3  - X(I) has only an upper bound. 
*     RI  SMALL           The smallest positive number such that
*                           1.0 + SMALL > 1.0.
*     
*      
      
      SUBROUTINE PROJGR(N,X,XL,XU,GA,PGA,PGANRM,ISGNGA,IB,SMALL)
      
*     Scalar Arguments
      INTEGER N,ISGNGA
      DOUBLE PRECISION PGANRM,SMALL
      
*     Array Arguments
      INTEGER IB(*)
      DOUBLE PRECISION X(*),XL(*),XU(*),GA(*),PGA(*)

*     Local Scalars
      INTEGER I
      DOUBLE PRECISION SMALL2

*     Intrinsic Functions
      INTRINSIC MAX,MIN


      SMALL2 = 1.0D+01*SMALL

      ISGNGA = 0
      PGANRM = 0.0D+00

      DO 10 I=1,N
      
         PGA(I) = GA(I)
         
         IF (IB(I) .GT. 0) THEN
            
            IF (IB(I) .LE. 2) THEN
               IF (X(I) .LE. XL(I) + SMALL2) THEN
                  PGA(I) = 0.0D+00
                  IF (GA(I) .LT. 0.0D+00) ISGNGA = 1
               END IF
            END IF
               
            IF (IB(I) .GE. 2) THEN
               IF (X(I) .GE. XU(I) - SMALL2) THEN
                  PGA(I) = 0.0D+00
                  IF (GA(I) .GT. 0.0D+00) ISGNGA = 1
               END IF
            END IF

         END IF

         PGANRM = PGANRM + PGA(I)*PGA(I)

 10   CONTINUE


      RETURN
      END


************************************************************************
*
*     * SUBROUTINE TERMIB *
*
*      
*     * Purpose *
*      
*     Calculation of stopping criterion and test for termination with 
*     the desired accuracy for LMBM-B.
*
*      
*     * Calling sequence *
*     
*      CALL TERMIB(MC,MCU,PGAHGA,PGANRM,WK,PWK,QK,BETALA,TOLG,TOLG2,
*     &     ISGNGA,ITERM,IFLAG)
*
*      
*     * Parameters *
*
*     IU  MC              Declared number of stored corrections.
*     II  MCU             Upper limit for maximum number of stored
*                           corrections.
*     RI  PGAHGA          PGAHGA = PGA'*H*PGA, where H denotes the 
*                           limited memory approximation of the inverse
*                           of the Hessian matrix and PGA denotes the 
*                           projected aggregate subgradient.
*     RI  PGANRM          Euclidean norm of the projected aggregate 
*                           subgradient vector.
*     RU  WK              Stopping criterion.
*     RO  PWK             Previous stopping criterion.
*     RO  QK              Second stopping criterion.
*     RI  BETALA          Aggregate locality measure.
*     II  ISGNGA          ISGNGA = 0, if GA(I)<=0 for all I such that 
*                           X(I)=XU(I) and GA(I)>=0 for all I such that
*                           X(I)=XL(I). Otherwise ISGNGA = 1.
*     RI  TOLG            Tolerance for the termination criterion.
*     RI  TOLG2           Tolerance for the second termination criterion.
*     IU  IFLAG           Index for adaptive version:
*                           0  - Maximum number of stored corrections
*                                has not been changed.
*                           1  - Maximum number of stored corrections
*                                has been changed.
*     IO  ITERM           Cause of termination:
*                           0  - The problem is not yet solved
*                                with desired accuracy.
*                           1  - The problem has been solved.
*                                with desired accuracy.
*
*
*     * Local parameters *
*
*     R   PAR             Multiplier for adaptive version. If WK <=
*                           PAR*TOLG then the adaptability may be applied.
*
*
      
      SUBROUTINE TERMIB(MC,MCU,PGAHGA,PGANRM,WK,PWK,QK,BETALA,TOLG,
     &     TOLG2,ISGNGA,ITERM,IFLAG)
      
*     Scalar Arguments
      INTEGER MC,MCU,ISGNGA,ITERM,IFLAG
      DOUBLE PRECISION PGAHGA,PGANRM,WK,PWK,QK,BETALA,TOLG,TOLG2

*     Parameters
      DOUBLE PRECISION PAR
      PARAMETER(PAR = 1.0D+04)


      PWK = WK
      WK = PGAHGA + 2.0D+00*BETALA
      QK = 0.5D+00*PGANRM + BETALA 

      IF (WK .LE. PAR*TOLG) THEN

         IF(QK .LE. TOLG2 .AND. WK .LE. TOLG) THEN
           
            IF (ISGNGA .EQ. 0) THEN 
               ITERM = 1
               GOTO 900
            END IF
         END IF
         
         IF (MC .LT. MCU .AND. IFLAG .EQ. 0) THEN
            MC = MC+1
            IFLAG = 1
         END IF
      END IF


 900  CONTINUE


      RETURN
      END


************************************************************************
*
*     * SUBROUTINE UPDATE *
*
*      
*     * Purpose *
*      
*     Updating the matrices RM, LM, CDIAG, SMTSM, and UMTUM for LMBM-B.
*
*      
*     * Calling sequence *
*     
*     CALL UPDATE(MN,INEW,IOLD,RM,LM,CDIAG,SMTSM,UMTUM,SMTU,
*    &     UMTU,SMTS,UMTS,STU,IFLAG,IFULL)
*     
*     
*     * Parameters *
*
*     II  MN              Current number of depositories used.
*     II  INEW            Index for circular arrays.
*     II  IOLD            Oldest correction.
*     II  IFLAG           Index for adaptive version:
*                           0  - Maximum number of stored corrections
*                                  has not been changed at previous
*                                  iteration.
*                           1  - Maximum number of stored corrections
*                                  has been changed at previous
*                                  iteration.
*     II  IFULL           Index for matrix updating:
*                           0  - No need for deleting previous values.
*                           1  - Old values must be deleted.
*     RU  RM(MN*(MN+1)/2) Upper triangular matrix stored columnwise
*                           in the one-dimensional array.
*     RU  LM(MN*(MN+1)/2) Lower triangular matrix stored rowwise in the
*                           one-dimensional array.
*     RU  UMTUM(MN*(MN+1)/2)  Matrix UMTUM = UM'*UM.
*     RU  SMTSM(MN*(MN+1)/2)  Matrix SMTSM = SM'*SM.
*     RU  CDIAG(MN)       Diagonal matrix.
*     RI  SMTS(MN)        Vector SMTS = SM'*S, where S is the 
*                           stored variable differences.
*     RI  UMTS(MN)        Vector UMTS = UM'*S.
*     RI  SMTU(MN)        Vector SMTS = SM'*U, where U is the 
*                           stored subgradient differences.
*     RI  UMTU(MN)        Vector UMTU = UM'*U.
*     RI  STU             Scalar STU = S'*U.
*
*      
*     * Subprograms used *
*      
*     S   COPY4           Copying of four vectors.
*     
*
*     All the MN-vectors are stored in a circular order controlled by
*     the parameter INEW.
*
*
  
      SUBROUTINE UPDATE(MN,INEW,IOLD,RM,LM,CDIAG,SMTSM,UMTUM,SMTU,
     &     UMTU,SMTS,UMTS,STU,IFLAG,IFULL)

      
*     Scalar Arguments
      INTEGER MN,INEW,IOLD,IFLAG,IFULL
      DOUBLE PRECISION STU
      
*     Array Arguments
      DOUBLE PRECISION RM(*),LM(*),CDIAG(*),SMTSM(*),UMTUM(*),SMTU(*),
     &     UMTU(*),SMTS(*),UMTS(*)

*     Local Scalars
      INTEGER I,J,K

*     External Subroutines
      EXTERNAL COPY4


*         
*     Update RM, LM, SMTSM, and UMTUM
*

      IF (IFULL.EQ.1 .OR. IFLAG.EQ.1) THEN


*     Shift the old values.

         DO 10 I=1,MN-1
            J=(I-1)*I/2+1
            K=I*(I+1)/2+2
            CALL COPY4(I,RM(K),RM(J),LM(K),LM(J),SMTSM(K),SMTSM(J),
     &           UMTUM(K),UMTUM(J))
 10      CONTINUE
      END IF


*     Add the new values.

      IF (IOLD .EQ. 1) THEN
         CALL COPY4(MN,SMTU,RM((MN-1)*MN/2+1),
     &        UMTS,LM((MN-1)*MN/2+1),
     &        SMTS,SMTSM((MN-1)*MN/2+1),
     &        UMTU,UMTUM((MN-1)*MN/2+1))

      ELSE

         CALL COPY4(MN-INEW,SMTU(IOLD),RM((MN-1)*MN/2+1),
     &        UMTS(IOLD),LM((MN-1)*MN/2+1),
     &        SMTS(IOLD),SMTSM((MN-1)*MN/2+1),
     &        UMTU(IOLD),UMTUM((MN-1)*MN/2+1))
         
         CALL COPY4(INEW,SMTU,RM((MN-1)*MN/2+MN-INEW+1),
     &        UMTS,LM((MN-1)*MN/2+MN-INEW+1),
     &        SMTS,SMTSM((MN-1)*MN/2+MN-INEW+1),
     &        UMTU,UMTUM((MN-1)*MN/2+MN-INEW+1))

      END IF

      LM((MN+1)*MN/2) = 0.0D+00


*
*     Update CDIAG
*

      CDIAG(INEW) = STU
            

      RETURN
      END
      

************************************************************************
*
*     * SUBROUTINE UPSERI *
*
*      
*     * Purpose *
*      
*     Matrix update by the limited memory BFGS update.
*
*      
*     * Calling sequence *
*     
*     CALL UPSERI(N,MC,MN,INEW,IOLD,IPD,IFLAG,IBFGS,G,GM,S,U,
*     &     SM,UM,LM,RM,CDIAG,SMTSM,UMTUM,SMTGM,UMTGM,VMC1,VMC2,
*     &     SMTS,UMTS,SMTU,UMTU,AMUGAD,XBX,TTHETA,GAMMA,SMALL,
*     &     ISCALE,IERR)
*
*      
*     * Parameters *
*
*     II  N               Number of variables.
*     II  MC              Declared number of stored corrections.
*     IU  MN              Current number of depositories used.
*     IU  INEW            Index for circular arrays.
*     IO  IOLD            Index of the oldest corrections.
*     IU  IPD             Index for positive definiteness for SR1 update:
*                           0  - Positive definiteness is preserved.
*                          -1  - Positive definiteness is not preserved.
*     IU  IFLAG           Index for adaptive version:
*                           0  - Maximum number of stored corrections
*                                has not been changed at previous
*                                iteration.
*                           1  - Maximum number of stored corrections
*                                has been changed at previous iteration.
*     IO  IBFGS           Index of the type of BFGS update:
*                           1  - BFGS update: the corrections are stored.
*                           3  - BFGS update is skipped.
*     RI  G(N)            Current subgradient of the objective function.
*     RI  GM(N)           Previous subgradient of the objective
*                           function.
*     RI  S(N)            Difference of current and previous variables.
*     RI  U(N)            Difference of current and previous
*                           subgradients of the Lagrangian.
*     RU  SM(N*MN)        Matrix whose columns are stored corrections.
*     RU  UM(N*MN)        Matrix whose columns are stored subgradient
*                           differences.
*     RU  LM(MN*(MN+1)/2) Lower triangular matrix stored rowwise
*                           in the one-dimensional array.
*     RU  RM(MN*(MN+1)/2) Upper triangular matrix stored columnwise
*                           in the one-dimensional array.
*     RU  CDIAG(MN)       Diagonal matrix.
*     RU  SMTSM(MN*(MN+1)/2)  Matrix SMTSM = SM'*SM.
*     RU  UMTUM(MN*(MN+1)/2)  Matrix UMTUM = UM'*UM.
*     RU  SMTGM(MN)       Vector SMTGM = SM'*GM.
*     RU  UMTGM(MN)       Vector UMTGM = UM'*GM.
*     RI  AMUGAD          AMUGAD = -(A*MU + G)'*D, where A*MU denotes 
*                           the Lagrange multipliers for problem.
*     RI  XBX             XBX = (XCP-X)'*BM*(XCP-X), where BM denotes 
*                           the limited memory approximation of the 
*                           Hessian matrix.
*     RI  TTHETA          Scaled stepsize: TTHETA=T*THETA.
*     RU  GAMMA           Scaling parameter.
*     RA  VMC#(MN)        Auxiliary arrays; # = 1,2.
*     RA  SMTS(MN)        Auxiliary array: Used to store vector 
*                           SMTS = SM'*S.
*     RA  UMTS(MN)        Auxiliary array: Used to store vector 
*                           UMTS = UM'*S.
*     RA  SMTU(MN)        Auxiliary array: Used to store vector 
*                           SMTS = SM'*U.
*     RA  UMTU(MN)        Auxiliary array: Used to store vector 
*                           UMTU = UM'*U.
*     RI  SMALL           Small positive value.
*     II  ISCALE          Selection of the scaling:
*                           0  - Scaling at every iteration with STU/UTU.
*                           1  - Scaling at every iteration with STS/STU.
*                           2  - Interval scaling with STU/UTU.
*                           3  - Interval scaling with STS/STU.
*                           4  - Preliminary scaling with STU/UTU.
*                           5  - Preliminary scaling with STS/STU.
*                           6  - No scaling.      
*     IO  IERR            Error indicador: 
*                           0  - Everything is ok.
*                         -10  - Warning: MN=0 and IBFGS=3.
*     
*     
*
*     * Local variables *
*
*     R   STU             STU = S'*U. 
*     R   STS             STS = S'*S. 
*     I   IFULL           Index for matrix updating.
*
*
*
*     * Subprograms used *
*
*     S   COPY2           Copying of two vectors.
*     S   INDIC2          Initialization of indices.
*     S   RWAXV2          Multiplication of two rowwise stored dense 
*                           rectangular matrices A and B by vectors X 
*                           and Y.
*     S   UPDATE          Updating limited memory matrices LM, RM, SMTSM
*                           UMTUM, and CDIAG.
*      
*     RF  VDOT            Dot product of two vectors.
*     RF  SCLPAR          Calculation of the scaling parameter.
*     
*
*     The variable and subgradient differences and all the MN-vectors are
*     stored in a circular order controlled by the parameter point INEW.
*
*      

      SUBROUTINE UPSERI(N,MC,MN,INEW,IOLD,IPD,IFLAG,IBFGS,G,GM,S,U,
     &     SM,UM,LM,RM,CDIAG,SMTSM,UMTUM,SMTGM,UMTGM,VMC1,VMC2,
     &     SMTS,UMTS,SMTU,UMTU,AMUGAD,XBX,TTHETA,GAMMA,SMALL,
     &     ISCALE,IERR)

*     Scalar Arguments
      INTEGER N,MC,INEW,IOLD,MN,IFLAG,IBFGS,ISCALE,IPD,IERR
      DOUBLE PRECISION AMUGAD,XBX,GAMMA,TTHETA,SMALL

*     Array Arguments
      DOUBLE PRECISION S(*),U(*),SM(*),UM(*),RM(*),LM(*),UMTUM(*),
     &     SMTSM(*),CDIAG(*),G(*),GM(*),SMTGM(*),UMTGM(*),VMC1(*),
     &     VMC2(*),SMTU(*),UMTU(*),SMTS(*),UMTS(*)

*     Local Scalars
      INTEGER I,IFULL
      DOUBLE PRECISION STU,STS

*     External Functions
      DOUBLE PRECISION VDOT,SCLPAR
      EXTERNAL VDOT,SCLPAR

*     External Subroutines
      EXTERNAL COPY2,INDIC2,RWAXV2,UPDATE

*     Intrinsic Functions
      INTRINSIC MAX,SQRT


      IERR = 0
      IFULL = 0

      STU = VDOT(N,S,U)
      STS = VDOT(N,S,S)


*
*     Positive definiteness
*

      IF (STU .GT. SQRT(SMALL)) THEN
     
         IF (-STU + TTHETA*TTHETA*MAX(XBX,AMUGAD) .GE. -SQRT(SMALL))
     &        THEN
            IPD = -1
         END IF

            
*     
*     Update matrices
*         

         IBFGS = 1


*
*     Initialization of indices.
*            

         IF (MN .EQ. MC) IFULL = 1
         CALL INDIC2(MC,MN,INEW,IOLD,IFLAG,IBFGS)


*     
*     Update SM and UM
*

         CALL COPY2(N,S,SM((INEW-1)*N+1),U,UM((INEW-1)*N+1))


*     
*     Computation of SM'*G, UM'*G, SM'*S, and UM'*S
*

         CALL RWAXV2(N,MN,SM,UM,G,G,VMC1,VMC2)
         CALL RWAXV2(N,MN,SM,UM,S,S,SMTS,UMTS)
            
            
*         
*     Computation of SM'*U and UM'*U
*

         IF (IOLD .EQ. 1) THEN
            DO 30 I=1,INEW-1
               SMTU(I)=VMC1(I) - SMTGM(I)
               SMTGM(I)=VMC1(I)
               UMTU(I)=VMC2(I) - UMTGM(I)
               UMTGM(I)=VMC2(I)
 30         CONTINUE

         ELSE
            DO 10 I=1,INEW-1
               SMTU(I)=VMC1(I) - SMTGM(I)
               SMTGM(I)=VMC1(I)
               UMTU(I)=VMC2(I) - UMTGM(I)
               UMTGM(I)=VMC2(I)
 10         CONTINUE

            DO 20 I=IOLD,MN
               SMTU(I)=VMC1(I) - SMTGM(I)
               SMTGM(I)=VMC1(I)
               UMTU(I)=VMC2(I) - UMTGM(I)
               UMTGM(I)=VMC2(I)
 20         CONTINUE
         END IF
            
         SMTU(INEW)=STU
         SMTGM(INEW)=VMC1(INEW)
         UMTU(INEW)=VDOT(N,U,U)
         UMTGM(INEW)=VMC2(INEW)


*         
*     Update RM, LM, SMTSM, UMTUM, and CDIAG.
*

         CALL UPDATE(MN,INEW,IOLD,RM,LM,CDIAG,SMTSM,UMTUM,SMTU,
     &        UMTU,SMTS,UMTS,STU,IFLAG,IFULL)

            
*         
*     Computation of GAMMA
*

         GAMMA = SCLPAR(MN,ISCALE,STS,STU,UMTU(INEW),SMALL)

         
         INEW = INEW + 1
         IF (INEW .GT. MC) INEW = 1


      ELSE

         
*     
*     BFGS update is skipped
*     

         IBFGS = 3

         IF (MN .EQ. 0) THEN

            IERR = -10
            RETURN
         END IF
         

*
*     Initialization of indices.
*            

         CALL INDIC2(MC,MN,INEW,IOLD,IFLAG,IBFGS)


*         
*     Computation of GAMMA
*
         
         IF (ISCALE .GE. 4) GAMMA = 1.0D+00

         
*         
*     Computation of SM'*G and UM'*G
*
         
         CALL RWAXV2(N,MN,SM,UM,G,G,SMTGM,UMTGM)

      END IF

      RETURN
      END


************************************************************************
*
*     * SUBROUTINE BFGSXV *
*
*      
*     * Purpose *
*      
*     Computation of the product H*V or -H*V, where H is the inverse
*     approximation of the Hessian calculated by the L-BFGS formula and
*     V is an arbitrary vector.
*
*      
*     * Calling sequence *
*     
*      CALL BFGSXV(N,MN,IOLD,Z,V,VN1,SM,UM,RM,UMTUM,CDIAG,
*     &     SMTV,UMTV,VMC1,VMC2,VMC3,GAMMA,SMALL,IERR,JOB)
*
*      
*     * Parameters *
*
*     II  N               Number of variables.
*     II  MN              Current size of vectors.
*     II  IOLD            Index of the oldest corrections.
*     RI  V(N)            Input vector.
*     RO  Z(N)            Output vector Z = H*V or Z = -H*V.
*     RI  SM(N*MN)        Matrix whose columns are stored corrections.
*     RI  UM(N*MN)        Matrix whose columns are stored subgradient
*                           differences.
*     RI  RM((MN+1)*MN/2) Upper triangular matrix stored columnwise in
*                           the one-dimensional array.
*     RI  UMTUM((MN+1)*MN/2)  Matrix UMTUM = UM'*UM.
*     RI  CDIAG(MN)       Diagonal matrix.
*     RI  SMTV(MN)        Vector SMTV = SM'*V.
*     RI  UMTV(MN)        Vector UMTV = UM'*V.
*     RI  GAMMA           Scaling parameter.
*     RA  VMC#(MN)        Auxiliary arrays; # = 1,...,3.
*     RA  VN1(N)          Auxiliary array.
*     II  JOB             Selection of the sign:
*                           0  - Z = H*V.
*                           1  - Z = -H*V.
*     IO  IERR            Error indicador: 
*                           0  - Everything is ok.
*                          -3  - Error in TRLIEQ.
*     
*     
*     * Subprograms used *
*      
*     S   RECMAX          Multiplication of a columnwise stored dense 
*                           rectangular matrix by a vector.
*     S   RWAXV2          Multiplication of two rowwise stored dense 
*                           rectangular matrices A and B by vectors X 
*                           and Y.
*     S   SCDIFF          Difference of the scaled vector and a vector.
*     S   SCSUM           Sum of a vector and the scaled vector.
*     S   SYMAX           Multiplication of a dense symmetric matrix
*                           by a vector.
*     S   TRLIEQ          Solving X from linear equation L*X=Y or
*                           L'*X=Y.
*     S   VXDIAG          Multiplication of a vector and a diagonal
*                           matrix.
*     S   XDIFFY          Difference of two vectors.
*     S   XSUMY           Sum of two vectors.
*     
*
*     The variable and subgradient differences and all the MN-vectors are
*     stored in a circular order controlled by the parameter point IOLD.
*
*
      
      SUBROUTINE BFGSXV(N,MN,IOLD,Z,V,VN1,SM,UM,RM,UMTUM,
     &     CDIAG,SMTV,UMTV,VMC1,VMC2,VMC3,GAMMA,SMALL,IERR,JOB)

*     Scalar Arguments
      INTEGER N,MN,IOLD,IERR,JOB
c     huom!!! smallia käytetään implisiittisesti.!!!
      DOUBLE PRECISION GAMMA,SMALL
      
*     Array Arguments
      DOUBLE PRECISION Z(*),V(*),VN1(*),SM(*),UM(*),RM(*),UMTUM(*),
     &     CDIAG(*),SMTV(*),UMTV(*),VMC1(*),VMC2(*),VMC3(*)

*     External Subroutines
      EXTERNAL RECMAX,RWAXV2,SCDIFF,SCSUM,SYMAX,TRLIEQ,VXDIAG,XDIFFY,
     &     XSUMY

      
*
*     Computation of two intermediate values VMC1 and VMC2
*

      CALL TRLIEQ(MN,MN,IOLD,RM,VMC1,SMTV,1,SMALL,IERR)
      IF (IERR .NE. 0) RETURN
      
      CALL SYMAX(MN,MN,IOLD,UMTUM,VMC1,VMC3)
      CALL VXDIAG(MN,CDIAG,VMC1,VMC2)
      CALL SCSUM(MN,GAMMA,VMC3,VMC2,VMC2)
      CALL SCSUM(MN,-GAMMA,UMTV,VMC2,VMC3)
      CALL TRLIEQ(MN,MN,IOLD,RM,VMC2,VMC3,0,SMALL,IERR)
      IF (IERR .NE. 0) RETURN

      
*
*     Computation of the product Z = H*V or Z = -H*V
*

*     Z = UM*VMC1 
      CALL RECMAX(N,MN,UM,VMC1,Z)
        

      IF (JOB .EQ. 1) THEN

*     Z = UM*VMC1 - V
         CALL XDIFFY(N,Z,V,Z)

*     VN1 = SM*VMC2
         CALL RECMAX(N,MN,SM,VMC2,VN1)

*     Z = GAMMA(UM*VMC1  - V) - SM*VMC2
         CALL SCDIFF(N,GAMMA,Z,VN1,Z)

      ELSE

*     Z = V - UM*VMC1
         CALL XDIFFY(N,V,Z,Z)
         
*     VN1 = SM*VMC2
         CALL RECMAX(N,MN,SM,VMC2,VN1)

*     Z = SM*VMC2 + GAMMA(V - UM*VMC1)
         CALL SCSUM(N,GAMMA,Z,VN1,Z)

      END IF

      RETURN
      END


*************************************************************************
*
*     * SUBROUTINE CPBFGS *
*
*      
*     * Purpose *
*      
*     Computation of generalized Cauchy point by the limited memory
*     BFGS update.
*
*      
*     * Calling sequence *
*     
*     CALL CPBFGS(N,MN,IOLD,X,XCP,XL,XU,IB,IT,G,DC,XBX,T,SMTGM,UMTGM,
*    &     SMTSM,SM,UM,LM,CDIAG,JT,PAUX1,PAUX2,CAUX1,CAUX2,V1,V2,W1,W2,
*    &     GAMMA,SMALL,IERR)
*
*      
*     * Parameters *
*
*     II  N               Number of variables.
*     II  MN              Current number of stored corrections.
*     II  IOLD            Index for circular arrays.
*     RI  X(N)            Vector of variables. 
*     RO  XCP(N)          Generalized Cauchy point.
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables. 
*     II  IB(N)           Type of bound constraints:
*                             0  - X(I) is unbounded,
*                             1  - X(I) has only a lower bound,
*                             2  - X(I) has both lower and upper bounds,
*                             3  - X(I) has only an upper bound. 
*     RI  G(N)            Current subgradient of the objective
*                           function.
*     RO  XBX             XBX = (XCP-X)'*BM*(XCP-X), where BM denotes the 
*                           limited memory approximation of the Hessian
*                           matrix.
*     RI  SM(N*MN)        Matrix whose columns are stored corrections.
*     RI  UM(N*MN)        Matrix whose columns are stored subgradient
*                           differences.
*     RI  SMTGM(MN)       Vector SMTGM = SM'*GM.
*     RI  UMTGM(MN)       Vector UMTGM = UM'*GM.
*     RI  SMTSM((MN+1)*MN/2)  Matrix SMTSM = SM'*SM.
*     RI  GAMMA           Scaling parameter.
*     RI  LM((MN+1)*MN/2) Lower triangular matrix stored rowwise
*                           in the one-dimensional array.
*     RI  CDIAG(MN)       Diagonal matrix.
*     RA  JT((MN+1)*MN/2) The upper triangular matrix J' that is
*                           the Cholesky factorization of 
*                           1/GAMMA*SM'*SM + L*C^(-1)*L'.  
*     RA  DC(N)           Auxiliary array used to store Cauchy 
*                           direction.
*     RA  T(N)            Auxiliary array used to store the breakpoints.
*     IA  IT(N)           Index set of breakpoints:
*                             for I=1,...,NLEFT, IT(I) are the indices
*                                 of breakpoints that have not been
*                                 encountered,
*                             for i=NLEFT+1,...,NBREAK, IT(I) are the
*                                 indices of encountered breakpoints and
*                             for i=IFREE,...,N, IT(I) are the indices
*                                 of variables that have no bound
*                                 constraints along the search 
*                                 direction.
*     RA  PAUX1(MN)       Auxiliary array used to store the vector
*                           UM'*DC.
*     RA  PAUX2(MN)       Auxiliary array used to store the vector
*                           1/GAMMA*SM'*DC.
*     RA  CAUX1(MN)       Auxiliary array used to store the vector
*                           UM'*(XCP-X).
*     RA  CAUX2(MN)       Auxiliary array used to store the vector
*                           1/GAMMA*SM'*(XCP-X).
*     RA  V1(MN)          Auxiliary array.
*     RA  V2(MN)          Auxiliary array.
*     RA  W1(MN)          Auxiliary array.
*     RA  W2(MN)          Auxiliary array.
*     RI  SMALL           Small positive value.
*     IO  IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -1  - Error in CHOFA (in FRMJ2).
*                            -3  - Error in TRLIEQ (in BMXV2).
*                            -8  - Error in CPBFGS, a nonpositive
*                                  definite matrix detected.
*     
*
*     * Subprograms used *
*      
*     S   SCSUM           Sum of a vector and the scaled vector.
*     S   HPSRT           Heapsort algorithm.
*     S   FRMJ2           Formation of the JJ' Cholesky factorization of
*                           matrix JT = SM'*SM/GAMMA+LM*CDIAG^(-1)*LM'.
*     S   BMXV2           The product of the 2MN x 2MN middle matrix in
*                           the compact L-BFGS formula and a 2MN vector.
*
*     RF  VDOT            Dot product of two vectors.
*     
*
*     The variable and subgradient differences and all the MN-vectors are
*     stored in a circular order controlled by the parameter point IOLD.
*
*
      
      SUBROUTINE CPBFGS(N,MN,IOLD,X,XCP,XL,XU,IB,IT,G,DC,XBX,T,SMTGM,
     &     UMTGM,SMTSM,SM,UM,LM,CDIAG,JT,PAUX1,PAUX2,CAUX1,CAUX2,V1,V2,
     &     W1,W2,GAMMA,SMALL,IERR)
      
*     Scalar Arguments
      INTEGER N,MN,IOLD,IERR
      DOUBLE PRECISION XBX,GAMMA,SMALL

*     Array Arguments
      INTEGER IB(*),IT(*)
      DOUBLE PRECISION G(*),X(*),XCP(*),XL(*),XU(*),DC(*),T(*),SMTGM(*),
     &     UMTGM(*),SM(*),UM(*),LM(*),CDIAG(*),SMTSM(*),JT(*),
     &     PAUX1(*),PAUX2(*),CAUX1(*),CAUX2(*),V1(*),V2(*),W1(*),W2(*)
      
*     Local Scalars
      INTEGER I,J,L,NBREAK,IFREE,ITMIN,NINT,ISBP,NLEFT,IBOUND
      DOUBLE PRECISION F1,F2,TMIN,DLTMIN,TOLD,DELTAT,GISBP,GISBP2,ZISBP,
     &     WMC,WMP,WMW,SMALL2,TMP
      
*     External Functions
      DOUBLE PRECISION VDOT
      EXTERNAL VDOT

*     Intrinsic Functions
      INTRINSIC MAX,ABS,SQRT

*     External Subroutines
      EXTERNAL FRMJ2,HPSRT,SCSUM,BMXV2

      
*
*     Initialization
*      

      NBREAK = 0
      NLEFT  = 0
      IFREE  = N+1
      ITMIN  = 0
      IBOUND = 0
      F1     = 0.0D+00
      F2     = 0.0D+00
      XBX    = 0.0D+00
      TMIN   =  0.0D+00
      IERR = 0
c      SMALL2=1.0D+01*SMALL
      SMALL2=SQRT(SMALL)
      

*     
*     Temporarily set PAUX = [PAUX1 PAUX2]' = -[UM  SM]'*G
*     and CAUX = [CAUX1 CAUX2]' = [UM  1/GAMMA*SM]'*(XCP-X) = 0
*      
      
      DO 10 I=1,MN
         L = I+IOLD-1
         IF (L .GT. MN) L=L-MN
         PAUX1(I) = -UMTGM(L)
         PAUX2(I) = -SMTGM(L)

         CAUX1(I)= 0.0D+00
         CAUX2(I)= 0.0D+00
 10   CONTINUE


*
*     Computation of the Cauchy direction DC, the first derivative
*     F1= -(DC)'*DC, and the breakpoints T in each coordinate
*     direction. Identification of the smallest breakpoint TMIN. 
*      

      DO 20 I=1,N

         IF (IB(I) .NE. 0) THEN
            IF (IB(I).LE.2 .AND. ABS(G(I)) .LE. SMALL2 .AND. 
     &           X(I)-XL(I) .LE. SMALL2) THEN
               DC(I) = 0.0D+00
               XCP(I) = X(I)
               GOTO 20
            END IF

            IF (IB(I).GE.2 .AND. ABS(G(I)) .LE. SMALL2 .AND. 
     &           XU(I)-X(I) .LE. SMALL2) THEN
               DC(I) = 0.0D+00
               XCP(I) = X(I)
               GOTO 20
            END IF
         END IF

         IF (IB(I) .NE. 0 .AND. IB(I) .LE. 2 .AND. G(I) .GT. 
     &        0.0D+00) THEN


*
*     Breakpoint in lower bound
*

            T(NBREAK+1) = (X(I) - XL(I))/G(I)

            IF (T(NBREAK+1) .LE. SMALL2) THEN
               DC(I) = 0.0D+00


*     
*     Correction to PAUX = [PAUX1 PAUX2]' = [UM  SM]'*DC 
*

               DO 30 J=1,MN
                  L = J + IOLD - 1
                  IF (L .GT. MN) L=L-MN
                  PAUX1(J) = PAUX1(J)+UM(N*(L-1)+I)*G(I)
                  PAUX2(J) = PAUX2(J)+SM(N*(L-1)+I)*G(I)
 30            CONTINUE
               
            ELSE
               NBREAK = NBREAK + 1
               IT(NBREAK) = I
               DC(I) = -G(I)
               F1 = F1 - DC(I)*DC(I)


*     
*     Determination of the smallest breakpoint
*               

               IF (NBREAK .EQ. 1 .OR. T(NBREAK) .LT. TMIN) THEN
                  TMIN=T(NBREAK)
                  ITMIN=NBREAK
               END IF
            END IF

         ELSE IF (IB(I) .GE. 2 .AND. G(I) .LT. 0.0D+00) THEN


*
*     Breakpoint in upper bound
*               

            T(NBREAK+1) = (X(I) - XU(I))/G(I)
            
            IF (T(NBREAK+1) .LE. SMALL2) THEN
               DC(I) = 0.0D+00


*     
*     Correction to PAUX = [PAUX1 PAUX2]' = [UM  SM]'*DC 
*

               DO 40 J=1,MN
                  L = J + IOLD - 1
                  IF (L .GT. MN) L=L-MN
                  PAUX1(J) = PAUX1(J)+UM(N*(L-1)+I)*G(I)
                  PAUX2(J) = PAUX2(J)+SM(N*(L-1)+I)*G(I)
 40            CONTINUE
               
            ELSE
               NBREAK = NBREAK + 1
               IT(NBREAK) = I
               DC(I) = -G(I)
               F1 = F1 - DC(I)*DC(I)


*
*     Determination of the smallest breakpoint
*               

               IF (NBREAK .EQ. 1 .OR. T(NBREAK) .LT. TMIN) THEN
                  TMIN=T(NBREAK)
                  ITMIN=NBREAK
               END IF
            END IF
            
         ELSE


*
*     No breakpoint 
*

            IFREE = IFREE - 1
            IT(IFREE) = I
            DC(I) = -G(I)
            F1 = F1 - DC(I)*DC(I)
            
            IF (ABS(G(I)) .GT. 0.0D+00) IBOUND=-1

         END IF
         

*     
*     Initialization of Cauchy point
*

         XCP(I) = X(I)
           
 20   CONTINUE


*     
*     Correction to PAUX2 = 1/GAMMA*SM'*DC 
*

      IF (GAMMA .NE. 1.0D+00) THEN
         DO 70 J=1,MN
            PAUX2(J) = PAUX2(J)/GAMMA
 70      CONTINUE
      END IF


*     The indices corresponding to free variables are located in
*     IT(1),...,IT(NBREAK) and IT(IFREE),...,IT(N). The smallest of the
*     NBREAK breakpoints is in T(ITMIN)=TMIN.

      
*
*     DC is a zero vector. Return with the initial XCP.
*

      IF (NBREAK .EQ. 0 .AND. IBOUND .EQ. 0) GOTO 999


*
*     Initialization of derivative F2 = (DC)' BM DC, where BM is 
*     representated by the compact limited memory BFGS formula.
*

      F2 =  -F1/GAMMA


*
*     The product of the 2MN x 2MN middle matrix in the compact
*     L-BFGS formula and the 2MN vector P; The product is returned
*     in V.
*

      CALL FRMJ2(MN,IOLD,SMTSM,JT,LM,CDIAG,GAMMA,SMALL,IERR)
      IF (IERR .NE. 0) RETURN
      
      CALL BMXV2(MN,IOLD,PAUX1,PAUX2,V1,V2,CDIAG,LM,JT,SMALL,IERR)
      IF (IERR .NE. 0) RETURN

      F2 = F2 - VDOT(MN,PAUX1,V1) - VDOT(MN,PAUX2,V2)
      

      IF (F2 .LE. SMALL2) THEN
         IERR = -8
         RETURN
      END IF

      DLTMIN = -F1/F2

      TOLD = 0.0D+00
      NINT = 1

      NLEFT = NBREAK
     
      
*     There is no breakpoints, locate the GCP and return. 
 
 
      IF (NBREAK .EQ. 0) GOTO 888

      
*
*     Begining of the loop.
*      

 100  CONTINUE

      IF (NINT .EQ. 1) THEN


*         
*     the smallest breakpoint is TMIN=T(ITMIN) and its index is
*     in IT(ITMIN)
*         

         ISBP = IT(ITMIN)

      ELSE


*         
*     the smallest breakpoint is chosen by heapsort algorithm
*

         IF (NINT .EQ. 2) THEN


*            
*     Remove the already used smallest breakpoint before heapsort call.
*            

            IF (ITMIN .NE. NBREAK) THEN
               T(ITMIN) = T(NBREAK)
               IT(ITMIN) = IT(NBREAK)
            END IF


*     
*     Heapsort with initialization of the heap.
*

            CALL HPSRT(NLEFT,T,IT,0)

         ELSE


*            
*     Heapsort with updating the heap.
*            

            CALL HPSRT(NLEFT,T,IT,1)

         END IF

         TMIN = T(NLEFT)
         ISBP = IT(NLEFT)

      END IF


      DELTAT = TMIN - TOLD

      IF (DELTAT .LT. 0.0D+00) THEN
         PRINT*,' no nyt voit jo hirttää ihtes.'
      END IF
      

*
*     Minimizer is within this interval, locate the GCP and return.
*

      IF (DLTMIN .LT. DELTAT) GOTO 888


*      
*     Examination of subsequent seqment
*
      
      NLEFT = NLEFT - 1
      NINT = NINT + 1
      GISBP = G(ISBP)
      DC(ISBP) = 0.0D+00

      IF (GISBP .LT. 0.0D+00) THEN
         XCP(ISBP) = XU(ISBP)
         ZISBP = XU(ISBP) - X(ISBP)
      ELSE 
         XCP(ISBP) = XL(ISBP)
         ZISBP = XL(ISBP) - X(ISBP)
      END IF


*         
*     All  variables are fixed, return with current XCP as GCP.
*         

      IF (NLEFT .EQ. 0 .AND. NBREAK .EQ. N) GOTO 999


*      
*     Update F1 and F2
*      

      GISBP2=GISBP*GISBP


*
*     Initial updating of F1 and F2
*      

      F1 = F1 + DELTAT*F2 + GISBP2 + GISBP*ZISBP/GAMMA
      TMP = GISBP2/GAMMA
      F2 = F2 - TMP


*         
*     Update CAUX1= CAUX1 + DELTAT*PAUX1 and CAUX2= CAUX2 + DELTAT*PAUX2.
*         

      CALL SCSUM(MN,DELTAT,PAUX1,CAUX1,CAUX1)
      CALL SCSUM(MN,DELTAT,PAUX2,CAUX2,CAUX2)


*     
*     Selection of W, the ISBP:th row of [UM  1/GAMMA*SM]. 
*
         
      DO 50 J=1,MN
         L = J + IOLD - 1
         IF (L .GT. MN) L=L-MN
         W1(J)= UM(N*(L-1)+ISBP)
         W2(J)= SM(N*(L-1)+ISBP)/GAMMA
 50   CONTINUE


*
*     Computation of products WMC, WMP and WMW
*

      CALL BMXV2(MN,IOLD,W1,W2,V1,V2,CDIAG,LM,JT,SMALL,IERR)
      IF (IERR .NE. 0) RETURN


      WMC = VDOT(MN,CAUX1,V1) + VDOT(MN,CAUX2,V2) 
      WMP = VDOT(MN,PAUX1,V1) + VDOT(MN,PAUX2,V2)
      WMW = VDOT(MN,W1,V1) + VDOT(MN,W2,V2)


*         
*     Update PAUX1 = PAUX1 + GISBP*W1 and PAUX2 = PAUX2 + GISBP*W2. 
*

      CALL SCSUM(MN,GISBP,W1,PAUX1,PAUX1)
      CALL SCSUM(MN,GISBP,W2,PAUX2,PAUX2)
         

*
*     Complete updating of F1 and F2
*      

      F1 = F1 - GISBP*WMC
      F2 = F2 - 2.0D+00*GISBP*WMP - GISBP2*WMW
      TMP = TMP + 2.0D+00*GISBP*WMP + GISBP2*WMW
      XBX = XBX + TMIN*TMIN * TMP


      IF (F2 .LE. SMALL2) THEN
c     DC is a zero vector. Return with current XCP.

         GOTO 999
      END IF

      DLTMIN = -F1/F2
      TOLD=TMIN


*     
*     End of the loop.
*

      IF (NLEFT .GT. 0) GOTO 100

      IF (IBOUND .EQ.0) DLTMIN = 0.0D+00   


 888  CONTINUE

      DLTMIN = MAX(0.0D+00,DLTMIN)
      TOLD = TOLD + DLTMIN

      XBX = XBX + TOLD*TOLD * F2


*      
*     Update free variables 
* 

      DO 60 I=1,NLEFT
         XCP(IT(I))=X(IT(I))+TOLD*DC(IT(I))
 60   CONTINUE
      DO 80 I=IFREE,N
         XCP(IT(I))=X(IT(I))+TOLD*DC(IT(I))
 80   CONTINUE


 999  CONTINUE
      
     
      RETURN
      END


*************************************************************************
*
*     * SUBROUTINE SSBFGS *
*
*      
*     * Purpose *
*      
*     Subspace minimization and computation of the search direction D
*     by the limited memory BFGS update.
*
*      
*     * Calling sequence *
*     
*     CALL SSBFGS(N,NACT,MN,IOLD,ITYPE,IACT,IB,X,XCP,XL,XU,
*    &     G,D,AMUGAD,XBX,ALPHA,SM,UM,RM,UMTUM,CDIAG,SMTGM,UMTGM,
*    &     VN1,VN2,KM,VMC1,VMC2,VMC3,VMC4,VMC5,GAMMA,SMALL,IERR)
*     
*      
*     * Parameters *
*
*     II  N               Number of variables.
*     II  NACT            Number of active variables.
*     II  MN              Current number of stored corrections.
*     II  IOLD            Index of the oldest correction.
*     II  ITYPE           Type of problem:
*                             0  - problem is unbounded,
*                             1  - constrained problem.
*     II  IACT(N)         Index set of active and free variables:
*                             for I=1,...,NACT, IACT(I) are the indices
*                                 of active variables (-J, if
*                                 XCP(J)=XL(J) and +J, otherwise) and
*                             for i=NACT+1,...,N, IACT(I) are the indices
*                                 of free variables.
*     II  IB(N)           Type of bound constraints:
*                             0  - X(I) is unbounded,
*                             1  - X(I) has only a lower bound,
*                             2  - X(I) has both lower and upper bounds,
*                             3  - X(I) has only an upper bound. 
*     RI  X(N)            Vector of variables. 
*     RI  XCP(N)          Generalized Chauchy point.
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables. 
*     RI  G(N)            Subgradient of the objective function.
*     RO  D(N)            Search direction.
*     RO  AMUGAD          AMUGAD = -(A*MU + G)'*D, where A*MU denotes
*                           the Lagrange multipliers for problem.
*     RU  XBX             XBX = (XCP-X)'*BM*(XCP-X), where BM denotes the 
*                           limited memory approximation of the Hessian
*                           matrix.
*     RO  ALPHA           Backtrack multiplier.
*     RO  DSTRN           Euclidean norm of the direction vector without 
*                           backtracking.
*     RI  SM(N*MN)        Matrix whose columns are stored corrections.
*     RI  UM(N*MN)        Matrix whose columns are stored subgradient
*                           differences.
*     RI  RM(MN*(MN+1)/2) Upper triangular matrix stored columnwise
*                           in the one-dimensional array.
*     RI  UMTUM(MN*(MN+1)/2)  Matrix UMTUM = UM'*UM.
*     RI  CDIAG(MN)       Diagonal matrix.
*     RI  SMTGM(MN)       Vector SMTGM = SM'*GM.
*     RI  UMTGM(MN)       Vector UMTGM = UM'*GM.
*     RI  GAMMA           Scaling parameter.
*     RI  SMALL           Small positive value.
*     RA  VMC#(MN)      Auxiliary arrays; # = 1,...,5.
*     RA  VN#(N)        Auxiliary array; # = 1,2.
*     RA  KM((2*MN+1)*MN) Auxiliary 2*MN matrix.
*     RI  SMALL           Small positive value.
*     IO  IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -3  - Error in TRLIEQ.
*                            -4  - Error in FRMLEL.
*                            -5  - Error in KMXV2.
*     
*
*
*     * Subprograms used *
*      
*     S   FRMLEL          Computation of Cholesky factorization of
*                           indefinite matrix KM.
*     S   KMXV2           Computation of the product of the inverse of
*                           2MN x 2MN matrix KM and and a 2MN vector P.
*     S   RECMAX          Multiplication of a vector by a dense
*                           rectangular matrix.
*     S   SCALEX          Scaling a vector.
*     S   SCDIFF          Difference of the scaled vector and a vector.
*     S   SCSUM           Sum of a vector and the scaled vector.
*     S   SYMAX           Multiplication of a dense symmetric matrix
*                           by a vector.
*     S   TRLIEQ          Solving X from linear equation L*X=Y or
*                           L'*X=Y.
*     S   VXDIAG          Multiplication of a vector and a diagonal
*                           matrix.
*     S   XDIFFY          Difference of two vectors.
*     S   XSUMY           Sum of two vectors.
*
*     RF  VDOT            Dot product of two vectors.
*
*
*     All the MN-vectors are stored in a circular order controlled by
*     the parameter point IOLD.
*
*
      
      SUBROUTINE SSBFGS(N,NACT,MN,IOLD,ITYPE,IACT,IB,X,XCP,XL,XU,
     &     G,D,AMUGAD,XBX,DSTRN,ALPHA,SM,UM,RM,UMTUM,CDIAG,SMTGM,UMTGM,
     &     VN1,VN2,KM,VMC1,VMC2,VMC3,VMC4,VMC5,GAMMA,SMALL,IERR)

*     Scalar Arguments
      INTEGER N,NACT,MN,IOLD,ITYPE,IERR
      DOUBLE PRECISION AMUGAD,XBX,DSTRN,ALPHA,GAMMA,SMALL
      
*     Array Arguments
      INTEGER IACT(*),IB(*)
      DOUBLE PRECISION X(*),XCP(*),XL(*),XU(*),D(*),G(*),
     &     SM(*),UM(*),RM(*),UMTUM(*),CDIAG(*),SMTGM(*),UMTGM(*),
     &     VN1(*),VN2(*),KM(*),VMC1(*),VMC2(*),VMC3(*),
     &     VMC4(*),VMC5(*)

*     Local Scalars
      INTEGER I,J,IFREE,IND
      DOUBLE PRECISION DFREE,SIGNUM

*     Intrinsic Functions
      INTRINSIC ABS,SIGN,MIN,DBLE

*     External Functions
      DOUBLE PRECISION VDOT
      EXTERNAL VDOT

*     External Subroutines
      EXTERNAL FRMLEL,KMXV2,RECMAX,SCALEX,SCDIFF,SCSUM,SYMAX,
     &     TRLIEQ,VXDIAG,XDIFFY,XSUMY,BFGSXV


      IERR = 0

      
      IF (NACT .EQ. 0) THEN


*
*     Computation of the search direction D = - HM*G
*

         CALL BFGSXV(N,MN,IOLD,D,G,VN1,SM,UM,RM,UMTUM,CDIAG,SMTGM,UMTGM,
     &        VMC1,VMC2,VMC3,GAMMA,SMALL,IERR,1)

         AMUGAD = -VDOT(N,G,D)

      ELSE


*
*     Computation of two intermediate vectors VMC1 and VMC2
*

         CALL TRLIEQ(MN,MN,IOLD,RM,VMC1,SMTGM,1,SMALL,IERR)
         IF (IERR .NE. 0) RETURN
         
         CALL SYMAX(MN,MN,IOLD,UMTUM,VMC1,VMC3)
         CALL VXDIAG(MN,CDIAG,VMC1,VMC2)
         CALL SCSUM(MN,GAMMA,VMC3,VMC2,VMC2)
         CALL SCSUM(MN,-GAMMA,UMTGM,VMC2,VMC3)
         CALL TRLIEQ(MN,MN,IOLD,RM,VMC2,VMC3,0,SMALL,IERR)
         IF (IERR .NE. 0) RETURN


*
*     Computation of the intermediate vector VN1= -A'DM G - A'(XCP-X)
*     

         DO 10 I=1,NACT
            VN1(I)=0.0D+00
            IND=ABS(IACT(I))
            SIGNUM = SIGN(1.0D+00,DBLE(IACT(I)))
            D(I) = SIGNUM*(XCP(IND)-X(IND))
            
            DO 11 J = 1,MN
               VN1(I)=VN1(I) + SIGNUM*(SM((J-1)*N+IND)*VMC2(J)
     &              - GAMMA*UM((J-1)*N+IND)*VMC1(J))
 11         CONTINUE

            VN1(I) = -SIGNUM*GAMMA*G(IND) - VN1(I) - D(I) 
 10      CONTINUE
         


*
*     Computation of Lagrange multipiliers VN1 = (A'DM A)^{-1} VN1.
*         
     
*     [GAMMA*UM SM]' A VN1 (Note the reorder of blocks).
            
         DO 20 J = 1, MN
            VMC3(J)= 0.0D+00
            VMC4(J)= 0.0D+00
            
            DO 21 I=1,NACT
               IND=ABS(IACT(I))
               SIGNUM = SIGN(1.0D+00,DBLE(IACT(I)))
               VMC3(J)=VMC3(J)+SIGNUM*GAMMA*UM((J-1)*N+IND)*VN1(I)
               VMC4(J)=VMC4(J)+SIGNUM*SM((J-1)*N+IND)*VN1(I)
 21         CONTINUE

 20      CONTINUE
            

*     Computation of [ VMC3 VMC4 ] = KM^{-1} [ VMC3 VMC4 ]

         CALL FRMLEL(N,MN,NACT,IOLD,KM,IACT,UM,SM,CDIAG,GAMMA,SMALL,
     &        IERR)
         IF (IERR .NE. 0) RETURN
         
         CALL KMXV2(MN,IOLD,VMC3,VMC4,VMC3,VMC4,KM,
     &        SMALL,IERR)
         IF (IERR .NE. 0) RETURN
         
         DO 50 I=1,NACT
            D(I)=0.0D+00
            IND=ABS(IACT(I))
            SIGNUM = SIGN(1.0D+00,DBLE(IACT(I)))
            DO 51 J = 1, MN
               D(I)=D(I) + SIGNUM*(SM((J-1)*N+IND)*VMC4(J)+
     &              GAMMA*UM((J-1)*N+IND)*VMC3(J))
 51         CONTINUE

            VN1(I) = (VN1(I) - D(I)/GAMMA) /GAMMA
 50      CONTINUE
      
      
*     
*     Computation of the search direction D = -DM(A VN1 + G)
*         

*     [SM UM]' A VN1
            
         DO 60 J = 1, MN
            VMC3(J)= 0.0D+00
            VMC4(J)= 0.0D+00
            DO 61 I=1,NACT
               IND=ABS(IACT(I))
               SIGNUM = SIGN(1.0D+00,DBLE(IACT(I)))
               VMC3(J)=VMC3(J)+SIGNUM*SM((J-1)*N+IND)*VN1(I)
               VMC4(J)=VMC4(J)+SIGNUM*UM((J-1)*N+IND)*VN1(I)
 61         CONTINUE
 60      CONTINUE

         DO J = 1, MN
            VMC3(J)= VMC3(J) + SMTGM(J)
            VMC4(J)= VMC4(J) + UMTGM(J)
         END DO

         CALL COPY(N,G,VN2)
         DO I=1,NACT
            IND=ABS(IACT(I))
            SIGNUM = SIGN(1.0D+00,DBLE(IACT(I)))
            
            VN2(IND) = VN2(IND) + SIGNUM*VN1(I)
         END DO

         CALL BFGSXV(N,MN,IOLD,D,VN2,VN1,SM,UM,RM,UMTUM,CDIAG,VMC3,
     &        VMC4,VMC1,VMC2,VMC5,GAMMA,SMALL,IERR,1) 

         AMUGAD = - VDOT(N,VN2,D)
         
      END IF


*     
*     Euclidean norm of the direction vector without backtracking.
*

      IF (ITYPE .EQ. 1) THEN
         J=0
         DO 220 I=1,NACT
            IND=ABS(IACT(I))
            IF(abs(D(IND)-(XCP(IND)-X(IND))) .GT. 1.0D-5) THEN
               J=J+1
               IF(J .EQ.1) then
                  PRINT*,'BFGS1',XCP(IND)-X(IND)-D(IND),XCP(IND)-X(IND),
     &                 D(IND),G(IND),MN,IOLD,NACT
               END IF
            END IF
            D(IND)=XCP(IND)-X(IND)
 220     CONTINUE
      END IF

      DSTRN = VDOT(N,D,D)

      
*
*     Backtrack to the feasible region if X + D violates bounds
*     (among free variables).
*
      
      IF (ITYPE .EQ. 1) THEN

         ALPHA = 1.0D+00
         DO 200 I = NACT+1,N
            IFREE=IACT(I)
            DFREE=D(IFREE)
            IF (IB(IFREE) .NE. 0) THEN
               IF (IB(IFREE) .LE. 2 .AND. 
     &              X(IFREE) + DFREE .LT. XL(IFREE)) THEN

*     Otherwise the lower bound is never met
 
                  IF (ABS(X(IFREE)-XCP(IFREE)+DFREE) .GT. SMALL) THEN
                     ALPHA=MIN(ALPHA,(XL(IFREE)-XCP(IFREE))/
     &                    (X(IFREE)-XCP(IFREE)+DFREE))
                  END IF

               ELSE IF (IB(IFREE) .GE. 2 .AND. 
     &                 X(IFREE) + DFREE  .GT. XU(IFREE)) THEN

*     Otherwise the upper bound is never met
 
                  IF (ABS(X(IFREE)-XCP(IFREE)+DFREE) .GT. SMALL) THEN
                     ALPHA=MIN(ALPHA,(XU(IFREE)-XCP(IFREE))/
     &                    (X(IFREE)-XCP(IFREE)+DFREE))
                  END IF
                  
               END IF
               
            END IF
 200     CONTINUE

         IF (ALPHA .LT. 1.0D+00) THEN
            DO 210 I=NACT+1,N
               IFREE=IACT(I)
               D(IFREE) = XCP(IFREE) - X(IFREE) +
     &              ALPHA*(X(IFREE) - XCP(IFREE) + D(IFREE))
 210        CONTINUE
         ELSE
            
            XBX = 0.0D+00

         END IF

            
      END IF
 
      RETURN
      END


*************************************************************************
*
*     * SUBROUTINE DLBFGS *
*
*      
*     * Purpose *
*      
*     Computation of the search direction by the limited memory BFGS 
*     update.
*
*      
*     * Calling sequence *
*     
*     CALL DLBFGS(N,NACT,MN,IOLD,ITYPE,IB,IACT,X,XCP,XL,XU,G,D,
*    &     SM,UM,LM,RM,CDIAG,SMTSM,UMTUM,SMTGM,UMTGM,VN1,VN2,KM,TMPMAT,
*    &     VMC1,VMC2,VMC3,VMC4,VMC5,VMC6,VMC7,VMC8,AMUGAD,XBX,ALPHA,
*    &     GAMMA,SMALL,IERR)
*
*      
*     * Parameters *
*
*     II  N               Number of variables.
*     IO  NACT            Number of active variables.
*     II  MN              Current number of depositories used.
*     II  IOLD            Index of the oldest corrections.
*     II  ITYPE           Type of problem:
*                           0  - problem is unbounded,
*                           1  - constrained problem.
*     II  IB(N)           Type of bound constraints:
*                           0  - X(I) is unbounded,
*                           1  - X(I) has only a lower bound,
*                           2  - X(I) has both lower and upper bounds,
*                           3  - X(I) has only an upper bound. 
*     IO  IACT(N)         Index set of active and free variables:
*                           for I=1,...,NACT, IACT(I) are the indices
*                             of active variables (-J, if XCP(J)=XL(J)
*                             and +J, otherwise)
*                           for i=NACT+1,...,N, IACT(I) are the indices
*                             of free variables.
*     RI  X(N)            Vector of variables.
*     RO  XCP(N)          Generalized Cauchy point.
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables. 
*     RU  D(N)            Search direction.
*     RI  G(N)            Current subgradient of the objective
*                           function.
*     RO  AMUGAD          AMUGAD = -(A*MU + G)'*D, where A*MU denotes
*                           the Lagrange multipliers for problem.
*     RO  XBX             XBX = (XCP-X)'*BM*(XCP-X), where BM denotes 
*                           the limited memory approximation of the 
*                           Hessian matrix.
*     RO  ALPHA           Backtrack multiplier.
*     RO  DSTRN           Euclidean norm of the direction vector without 
*                           backtracking.
*     RI  GAMMA           Scaling parameter.
*     RI  SM(N*MN)        Matrix whose columns are stored corrections.
*     RI  UM(N*MN)        Matrix whose columns are stored subgradient
*                           differences.
*     RI  LM(MN*(MN+1)/2) Lower triangular matrix stored rowwise
*                           in the one-dimensional array.
*     RI  RM(MN*(MN+1)/2) Upper triangular matrix stored columnwise
*                           in the one-dimensional array.
*     RI  CDIAG(MN)       Diagonal matrix.
*     RI  SMTSM(MN*(MN+1)/2)  Matrix SMTSM = SM' * SM.
*     RI  UMTUM(MN*(MN+1)/2)  Matrix UMTUM = UM' * UM.
*     RI  SMTGM(MN)       Vector SMTGM = SM' * G.
*     RI  UMTGM(MN)       Vector UMTGM = UM' * G.
*     RA  KM((2MN+1)*MN)  Auxiliary 2MN matrix.
*     RA  TMPMAT(MN*(MN+1)/2) Auxiliary MN+1 matrix.
*     RA  VMC#(MN)      Auxiliary arrays; # = 1,...,8.
*     RA  VN#(N)        Auxiliary arrays; # = 1,2.
*     RI  SMALL           Small positive value.
*     IO  IERR            Error indicador: 
*                           0  - Everything is ok.
*                          -1  - Error in CHOFA (in CPBFGS).
*                          -3  - Error in TRLIEQ.
*                          -4  - Error in FRMLEL (in SSBFGS).
*                          -5  - Error in KMXV2 (in SSBFGS).
*                          -8  - Error in CPBFGS.
*     
*
*
*     * Subprograms used *
*      
*     S   ACTVAR          Finding the index set of free and active
*                           variables  at the generalized Cauchy point.
*     S   CPBFGS          Computing the generalized Cauchy point by the
*                           limited memory BFGS update formula.
*     S   SSBFGS          Subspace minimization and direction finding
*                           using limited memory BFGS matrices.
*     S   XDIFFY          Difference of two vectors.
*     
*
*     The variable and subgradient differences and all the MN-vectors 
*     are stored in a circular order controlled by the parameter point
*     IOLD.
*
*
      
      SUBROUTINE DLBFGS(N,NACT,MN,IOLD,ITYPE,IB,IACT,X,XCP,XL,XU,G,D,SM,
     $     UM,LM,RM,CDIAG,SMTSM,UMTUM,SMTGM,UMTGM,VN1,VN2,KM,TMPMAT,
     $     VMC1,VMC2,VMC3,VMC4,VMC5,VMC6,VMC7,VMC8,AMUGAD,XBX,DSTRN,
     $     ALPHA,GAMMA,SMALL,IERR)
   
*     Scalar Arguments
      INTEGER N,NACT,MN,IOLD,ITYPE,IERR

      DOUBLE PRECISION AMUGAD,XBX,ALPHA,DSTRN,GAMMA,SMALL
      
*     Array Arguments
      INTEGER IACT(*),IB(*)
      DOUBLE PRECISION X(*),XCP(*),XL(*),XU(*),D(*),G(*),SM(*),UM(*),
     &     LM(*),RM(*),CDIAG(*),SMTSM(*),UMTUM(*),SMTGM(*),UMTGM(*),
     &     VN1(*),VN2(*),KM(*),TMPMAT(*),VMC1(*),VMC2(*),VMC3(*),
     &     VMC4(*),VMC5(*),VMC6(*),VMC7(*),VMC8(*)

*     External Subroutines
      EXTERNAL ACTVAR,CPBFGS,SSBFGS,XDIFFY


      IERR = 0


      IF (ITYPE .EQ. 1) THEN

         
*
*     Computation of Cauchy point.
*

         CALL CPBFGS(N,MN,IOLD,X,XCP,XL,XU,IB,IACT,G,D,XBX,VN1,SMTGM,
     &        UMTGM,SMTSM,SM,UM,LM,CDIAG,TMPMAT,VMC1,VMC2,VMC3,VMC4,
     &        VMC5,VMC6,VMC7,VMC8,GAMMA,SMALL,IERR)

         IF (IERR .NE. 0) RETURN


*
*     Determination of the active variables.
*

         CALL ACTVAR(N,NACT,IACT,IB,XCP,XL,XU,SMALL)

      ELSE
         
         XBX = 0.0D+00
         
      END IF


*
*     Subspace minimization and direction finding.
*

      IF (NACT .NE. N) THEN
     
         CALL SSBFGS(N,NACT,MN,IOLD,ITYPE,IACT,IB,X,XCP,XL,XU,G,D,
     &        AMUGAD,XBX,DSTRN,ALPHA,SM,UM,RM,UMTUM,CDIAG,SMTGM,UMTGM,
     &        VN1,VN2,KM,VMC1,VMC2,VMC3,VMC4,VMC5,GAMMA,SMALL,IERR)


      ELSE


*     If there are no free variables, then skip the subspace
*     minimization D=XCP-X.

         CALL XDIFFY(N,XCP,X,D)
         AMUGAD = XBX

      END IF

      RETURN
      END


************************************************************************
*
*     * SUBROUTINE FRMJ2 *
*
*      
*     * Purpose *
*      
*     Formation of the JJ' Cholesky factorization of the matrix
*      
*          JT = SM' SM / GAMMA + LM CDIAG^(-1) LM'.
*
*      
*     * Calling sequence *
*     
*     CALL FRMJ2(MN,IOLD,SMTSM,JT,LM,CDIAG,GAMMA,IERR)
*
*      
*     * Parameters *
*
*     II  MN              Current size of vectors. 
*     II  IOLD            Index of the oldest correction.
*     RO  JT((MN+1)*MN/2) The upper triangular matrix J' that is
*                             the Cholesky factorization of
*                             (SM'SM/GAMMA+LM C^(-1) LM').  
*     RI  LM((MN+1)*MN/2) Lower triangular matrix stored rowwise
*                             in the one-dimensional array.
*     RI  SMTSM((MN+1)*MN/2)  Matrix SMTSM = SM' * SM.
*     RI  CDIAG(MN)       Diagonal matrix stored in circular order.
*     RI  GAMMA           Scaling parameter.
*     RI  SMALL           Small positive value.
*     IO  IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -1  - Error: matrix is not positive
*                                  definite.
*
*
*     * Subprograms used *
*      
*     S   CHOFA           Cholesfy factorization.
*     
*     
      
      SUBROUTINE FRMJ2(MN,IOLD,SMTSM,JT,LM,CDIAG,GAMMA,SMALL,IERR)
      
*     Scalar Arguments
      INTEGER MN,IOLD,IERR
      DOUBLE PRECISION GAMMA,SMALL
      
*     Array Arguments
      DOUBLE PRECISION SMTSM(*),JT(*),LM(*),CDIAG(*)

*     Local Scalars
      INTEGER I,IJ,JK,IK,K,J,JC
      
*     External Subroutines
      EXTERNAL CHOFA      
    
      IERR = 0
      IF (MN .EQ. 0) RETURN
      
*     Form JT = SM' SM / GAMMA + LM C^(-1)LM'.

      IJ = 0
      DO 10 I = 1, MN
         IK = 0
         DO 20 K = 1,I
            IJ=IJ+1
            JT(IJ) = SMTSM(IJ)/GAMMA
            DO 30 J = 1, K
               JC=J+IOLD-1
               IF (JC .GT. MN) JC=JC-MN
               IK = IK+1
               JK = IJ+J-K
               JT(IJ) = JT(IJ) + LM(IK)*LM(JK)/CDIAG(JC)
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE


*     Cholesky factorizion of JT = SM'SM/GAMMA + LM C^(-1)LM' to JJ'.
      
      CALL CHOFA(MN,JT,SMALL,IERR)

      RETURN
      END
      

************************************************************************
*
*     * SUBROUTINE BMXV2 *
*
*      
*     * Purpose *
*      
*     Computation of the product of the 2MN x 2MN middle matrix
*
*     [-CDIAG^(1/2) CDIAG^(-1/2)*LM']^(-1) [ CDIAG^(1/2)      0 ]^(-1)
*     [ 0           JT'             ]      [-LM*CDIAG^(-1/2)  JT] 
*      
*     in the compact L-BFGS formula and a 2MN vector P = [ P1 P2 ]'.
*     The product is returned on V = [ V1 V2 ]' .
*
*      
*     * Calling sequence *
*     
*     CALL BMXV2(MN,IOLD,P1,P2,V1,V2,CDIAG,LM,JT,SMALL,IERR)
*
*      
*     * Parameters *
*     
*
*     II  MN              Current number of stored corrections.
*     II  IOLD            Index for circular arrays.
*     RI  P1(MN)          Input vector.
*     RI  P2(MN)          Input vector.
*     RO  V1(MN)          Output vector.
*     RO  V2(MN)          Output vector.
*     RI  LM((MN+1)*MN/2) Lower triangular matrix stored rowwise
*                           in the one-dimensional array.
*     RI  CDIAG(MN)       Diagonal matrix.
*     RI  JT((MN+1)*MN/2) Upper triangular matrix J' that is the
*                         Cholesky factor of (1/GAMMA * SM'SM+LM 
*                           * CDIAG^(-1)LM').
*     RI  SMALL           Small positive value.
*     IU  IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -3  - Error; 0 at diagonal.
*     
*
*     The input matrix CDIAG is stored in a circular order controlled 
*     by the parameter point IOLD.
*
*
      
      SUBROUTINE BMXV2(MN,IOLD,P1,P2,V1,V2,CDIAG,LM,JT,SMALL,IERR)
      
*     Scalar Arguments
      INTEGER MN,IOLD,IERR
      DOUBLE PRECISION SMALL
      
*     Array Arguments
      DOUBLE PRECISION P1(*),P2(*),CDIAG(*),LM(*),V1(*),V2(*),JT(*)

*     Local Scalars
      INTEGER I,II,K,IK,J
      DOUBLE PRECISION TMP
      
*     Intrinsic Functions
      INTRINSIC SQRT

     
*
*     Initialization
*

      IERR = 0
      IF (MN .EQ. 0) RETURN


*
*     Solving [  CDIAG^(1/2)      0  ] [ v1 ] = [ p1 ]
*             [ -LM*CDIAG^(-1/2)  JT ] [ v2 ]   [ p2 ].

*     Solving CDIAG^(1/2) V1 = P1.
      
      DO 30 I = 1, MN
         J=I+IOLD-1
         IF (J .GT. MN) J=J-MN
         V1(I) = P1(I) / SQRT(CDIAG(J))
 30   CONTINUE 


*     
*     Solving JT V2 = P2 + LM CDIAG^(-1) P1.
*      
      IK = 0
      V2(1) = P2(1)
      DO 10 I = 2, MN
         IK=IK+1
         TMP = 0.0D+00
         DO 20 K = 1, I - 1
            J=K+IOLD-1
            IF (J .GT. MN) J=J-MN
            IK = IK + 1
            TMP = TMP + LM(IK)*P1(K) / CDIAG(J)
 20      CONTINUE
      
         V2(I) = P2(I) + TMP
 10   CONTINUE  


*
*     Solving the triangular system.
*

*     Vector V2 is not in circular order.

      CALL TRLIEQ(MN,MN,1,JT,V2,V2,0,SMALL,IERR)
      IF (IERR .NE. 0) RETURN

 
*
*     Solving [ -CDIAG^(1/2)   CDIAG^(-1/2)*LM' ] [ V1 ] = [ V1 ]
*             [  0             JT'              ] [ V2 ]   [ V2 ].
*      
*
*     Solving JT'V2=V2.
*      
 
      CALL TRLIEQ(MN,MN,1,JT,V2,V2,1,SMALL,IERR)
      IF (IERR .NE. 0) RETURN


*      
*     Computation of V1=-CDIAG^(-1/2)(V1-CDIAG^(-1/2)LM'V2)
*                      =-CDIAG^(-1/2)V1+CDIAG^(-1)LM'V2.
*

      DO 40 I = 1, MN
         J=I+IOLD-1
         IF (J .GT. MN) J=J-MN
         V1(I) = -V1(I) / SQRT(CDIAG(J))
 40   CONTINUE

      II = 1
      DO 50 I = 1, MN
         TMP = 0.0D+00
         J=I+IOLD-1
         IF (J .GT. MN) J=J-MN
         II = II + I + 1
         IK = II - 1
         DO 60 K = I+1, MN
            TMP = TMP + LM(IK)*V2(K)
            IK = IK + K
  60     CONTINUE
         TMP = TMP/CDIAG(J)
         V1(I) = V1(I) + TMP
  50  CONTINUE
      

      RETURN
      END

      
*************************************************************************
*
*     * SUBROUTINE FRMLEL *
*
*      
*     * Purpose *
*      
*     Formation of the LEL' factorization of the indefinite matrix
*      
*          KM = [-CDIAG - GAMMA*U'ZZ'U    L_a' - R_z'  ]
*               [L_a - R_z                S'AA'S/GAMMA ]
*     
*     where     E = [-I  0]
*                   [ 0  I]
*
*     and L_a is the strictly lower triangular part of S'AA'U and R_z is
*     the upper triangular part of S'ZZ'U.
*      
*
*     * Calling sequence *
*     
*     CALL FRMLEL(N,MN,NACT,IOLD,KM,IACT,UM,SM,CDIAG,GAMMA,SMALL,IERR)
*
*      
*     * Parameters *
*
*     II  N               Number of variables.
*     II  MN              Current number of stored corrections.
*     II  IOLD            Index for circular arrays.
*     II  NACT            Number of active variables.
*     II  IACT(N)         Index set of active and free variables:
*                             for I=1,...,NACT, IACT(I) are the indices
*                                 of active variables (-J, if
*                                 XCP(J)=XL(J) and +J, otherwise)
*                             for i=NACT+1,...,N, IACT(I) are the indices
*                                 of free variables.
*     RI  SM(N*MN)        Matrix whose columns are stored corrections.
*     RI  UM(N*MN)        Matrix whose columns are stored subgradient
*                           differences.
*     RI  GAMMA           Scaling parameter.
*     RI  CDIAG(MN)       Diagonal matrix.
*     RO  KM((2MN+1)*MN)  LEL' factorization of indefinite matrix.
*     RI  SMALL           Small positive value.
*     IO  IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -1  - Error in CHOFA.
*                            -3  - Error in TRLIEQ.
*                            -4  - Error: nonpositive value at diagonal.
*     
*     
*     * Subprograms used *
*
*     S   CHOFA           Cholesfy factorization.
*     S   TRLIEQ          Solving X from linear equation L*X=Y or
*                           L'*X=Y.
*
*     RF  VDOT            Dot product of two vectors.
*
*     
*     The input matrices SM, UM and CDIAG are stored in a circular order
*     controlled by the parameter point IOLD.
*
*
      
      SUBROUTINE FRMLEL(N,MN,NACT,IOLD,KM,IACT,UM,SM,CDIAG,GAMMA,SMALL,
     &     IERR)

*     Scalar Arguments
      INTEGER N,MN,NACT,IOLD,IERR
      DOUBLE PRECISION GAMMA,SMALL
      
*     Array Arguments
      INTEGER IACT(*)
      DOUBLE PRECISION KM(*),UM(*),SM(*),CDIAG(*)

*     Local Scalars
      INTEGER I,II,I1,IJ,IJ2,IJ3,J,JI,JJ,J1,K,LI,LJ,IND
      DOUBLE PRECISION SUM,TMP

*     External Functions
      DOUBLE PRECISION VDOT
      EXTERNAL VDOT
      
*     Intrinsic Functions
      INTRINSIC ABS,SQRT

*     External Subroutines
      EXTERNAL CHOFA,TRLIEQ
           

*
*     Initialization
*
      IERR = 0


*
*     Form the upper triangle
*     
*     KM = [CDIAG+GAMMA*U'ZZ'U    -L_a'+R_z'   ] 
*          [-L_a+R_z              S'AA'S/GAMMA ]
*
c     Note, this could be done more efficiently.

      IJ = 0
      IJ2 = MN*(MN+1)/2
      IJ3 = IJ2
      
      DO 10 I = 1,MN
         LI=I+IOLD-1
         IF (LI .GT. MN) LI=LI-MN

*     computation of -L_a'

         DO 11 J = 1,I-1
            LJ=J+IOLD-1
            IF (LJ .GT. MN) LJ=LJ-MN
            IJ3 = IJ3+1
            KM(IJ3)=0.0D+00
            DO 12 K = 1,NACT
               IND = ABS(IACT(K))
               KM(IJ3) = KM(IJ3) - SM((LI-1)*N+IND)*UM((LJ-1)*N+IND)
 12         CONTINUE
 11      CONTINUE

*     computation of R_z'

         DO 13 J = I,MN
            LJ=J+IOLD-1
            IF (LJ .GT. MN) LJ=LJ-MN
            IJ3 = IJ3+1
            KM(IJ3) = 0.0D+00
            TMP = 0.0D+00
            
            DO 14 K = NACT+1,N
               IND = IACT(K)
               KM(IJ3) = KM(IJ3) + SM((LI-1)*N+IND)*UM((LJ-1)*N+IND)
               TMP =  TMP + SM((LJ-1)*N+IND)*UM((LI-1)*N+IND)
 14         CONTINUE
 13      CONTINUE

         IJ2=IJ2+MN
         DO 15 J=1,I
            LJ=J+IOLD-1
            IF (LJ .GT. MN) LJ=LJ-MN
            IJ=IJ+1
            IJ2=IJ2+1

*     computation of GAMMA*U'ZZ'U

            KM(IJ)=0.0D+00
            DO 16 K = NACT + 1, N
               IND = IACT(K)
               KM(IJ) = KM(IJ) + UM((LJ-1)*N+IND)*UM((LI-1)*N+IND)
 16         CONTINUE
            KM(IJ)= GAMMA*KM(IJ)

*     computation of S'AA'S/GAMMA

            KM(IJ2)=0.0D+00
            DO 17 K = 1,NACT
               IND = ABS(IACT(K))
               KM(IJ2) = KM(IJ2) + SM((LJ-1)*N+IND)*SM((LI-1)*N+IND)
 17         CONTINUE
            KM(IJ2)= KM(IJ2)/GAMMA
 15      CONTINUE

         IJ3=IJ3+I
 10   CONTINUE
      
*     computation of CDIAG+GAMMA*U'ZZ'U

      II=0
      DO 20 I=1,MN
         LI=I+IOLD-1
         IF (LI .GT. MN) LI=LI-MN
         II=II+I
         KM(II)=CDIAG(LI)+KM(II)
 20   CONTINUE


*
*     Form the upper triangle of
*     
*     KM=[ LL'              L^-1(-La'+Rz')                              ]
*        [(-L_a+R_z)L'^-1   S'AA'S/GAMMA+(L^-1(-La'+Rz'))'L^-1(-La'+Rz')]
*      

*     First Cholesky factor (1,1) block of KM to get LL'
*     with L' stored in the (1,1) upper triangle of KM.

      CALL CHOFA(MN,KM,SMALL,IERR)
      IF (IERR .NE. 0) RETURN

*     then form L^-1(-L_a'+R_z') in the (1,2) block.

      K = MN*(MN+1)/2+1
      DO 30 I = 1, MN
         CALL TRLIEQ(MN,MN,1,KM,KM(K),KM(K),0,SMALL,IERR)
         IF (IERR .NE. 0) RETURN
         K=K+MN+I
 30   CONTINUE


*     
*     Form S'AA'S/GAMMA + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
*     upper triangle of (2,2) block of KM.
*

      IJ2 = MN*(MN+1)/2
      JI = IJ2+1

      DO 40 I = 1,MN
         IJ2 = IJ2 + MN
         IJ = MN*(MN+1)/2+1
         DO 41 J = 1,I
            IJ2 = IJ2 + 1
            KM(IJ2) = KM(IJ2)+VDOT(MN,KM(IJ),KM(JI))
            IJ=IJ+J+MN
 41      CONTINUE
         JI=JI+I+MN
 40   CONTINUE


*      
*     Cholesky factorization of (2,2) block of KM.
*

      II= (MN+1)*(MN+2)/2
      I1= II

      DO 50 I=1,MN
         SUM=0.0D+00
         JI = I1
         JJ = (MN+1)*(MN+2)/2
         J1 = JJ
         DO 51 J=1, I-1
            TMP = KM(JI) - VDOT(J-1,KM(J1),KM(I1))
            TMP = TMP/KM(JJ)
            KM(JI)=TMP
            SUM=SUM+TMP*TMP
            JI = JI + 1
            J1 = JJ + MN + 1
            JJ = JJ + MN + 1 + J
 51      CONTINUE
         SUM = KM(II) - SUM

         IF (SUM .LE. SMALL) THEN
            IERR = -4
            RETURN
         END IF

         KM(II) = SQRT(SUM)
         I1 = II + MN + 1
         II = II + MN + 1 + I
 50   CONTINUE
     
      RETURN
      END
      

************************************************************************
*     
*     * SUBROUTINE KMXV2 *
*
*      
*     * Purpose *
*      
*     Computation of the product of the inverse of 2MN x 2MN matrix
*     
*         KM = [-CDIAG - GAMMA*U'ZZ'U    L_a' - R_z'  ]
*              [L_a - R_z                S'AA'S/GAMMA ]
*      
*     and a 2MN vector [ P1 P2 ]. The indefinite matrix KM is given
*     using Cholesky factorization KM=L*E*L' obtained by the subroutine
*     FRMLEL.
*
*      
*     * Calling sequence *
*     
*     CALL KMXV2(MN,IOLD,P1,P2,V1,V2,KM,SMALL,IERR)
*
*      
*     * Parameters *
*     
*
*     II  MN              Current number of corrections.
*     II  IOLD            Index of the oldest corrections.
*     RI  KM(MN*(2*MN+1)) Factorization KM=L*E*L' obtained by the
*                           subroutine FRMLEL.
*     RI  P1(MN)          Input vector.
*     RI  P2(MN)          Input vector
*     RO  V1(MN)          Output vector.
*     RO  V2(MN)          Output vector.
*     RI  SMALL           Small positive value.
*     IO  IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -3  - Error in TRLIEQ.
*                            -5  - Error: 0 at diagonal.
*     
*     
*     * Subprograms used *
*
*     S   TRLIEQ          Solving X from linear equation L*X=Y or
*                           L'*X=Y.
*
*     
*     All the input vectors are stored in a circular order
*     controlled by the parameter point IOLD.
*
*
      
      SUBROUTINE KMXV2(MN,IOLD,P1,P2,V1,V2,KM,SMALL,IERR)
      
*     Scalar Arguments
      INTEGER MN,IOLD,IERR
      DOUBLE PRECISION SMALL
      
*     Array Arguments
      DOUBLE PRECISION P1(*),V1(*),P2(*),V2(*),KM(*)

*     Local Scalars
      INTEGER I,II,IJ,L,J,K

*     External Subroutines
      EXTERNAL TRLIEQ


      IERR = 0


      IF (MN .EQ. 0) RETURN


*
*     Compute V=KM**(-1)P , where V=[V1' V2']', P=[P1' P2']' and KM=LEL'.
*

*     Solving V1 (=E1*L1'V1) from L1 V1 = P1, where L1 is the upper part 
*     of L

      CALL TRLIEQ(MN,MN,IOLD,KM,V1,P1,0,SMALL,IERR)
      IF (IERR .NE. 0) RETURN
      
*     Solving the rest of the triangular system: V=L**(-1)*P.

      II = MN*(MN+1)/2
      IJ = II
      DO 10 I = 1,MN
         II = II + I + MN
         L=I+IOLD-1
         IF (L .GT. MN) L=L-MN

         IF (KM(II) .LE. SQRT(SMALL)) THEN
            IERR = -5
            RETURN
         END IF

         V2(L) = P2(L)
         
         DO 20 J = 1,MN
            IJ = IJ + 1
            K=J+IOLD-1
            IF (K .GT. MN) K=K-MN
            V2(L) = V2(L) - KM(IJ)*V1(K)
 20      CONTINUE

         DO 30 J = 1,I - 1
            IJ = IJ + 1
            K=J+IOLD-1
            IF (K .GT. MN) K=K-MN
            V2(L) = V2(L) - KM(IJ)*V2(K)
 30      CONTINUE

         IJ = IJ + 1
         V2(L) = V2(L)/KM(II)
 10   CONTINUE

      DO 40 I=1,MN
         V1(I) = -V1(I)
 40   CONTINUE


*     Solving V2 from L2' V2 = V2, where L2 is the lower MN part of L 

      II = MN * (2*MN+1)
      DO 50 I = MN,1,-1
         L=I+IOLD-1
         IF (L .GT. MN) L=L-MN
         IF (KM(II) .LE. SQRT(SMALL)) THEN
            IERR = -5
            RETURN
         END IF
         IJ = II 

         DO 60 J = I + 1,MN
            K=J+IOLD-1
            IF (K .GT. MN) K=K-MN
            IJ = IJ + MN + J - 1
            V2(L) = V2(L) - KM(IJ)*V2(K)
 60      CONTINUE

         V2(L)=V2(L)/KM(II)
         II = II - I - MN
 50   CONTINUE

      DO 70 I = MN,1,-1
         L=I+IOLD-1
         IF (L .GT. MN) L=L-MN
         IF (KM(II) .LE. SQRT(SMALL)) THEN
            IERR = -5
            RETURN
         END IF

         IJ = II 
         DO 80 J = I + 1,MN
            K=J+IOLD-1
            IF (K .GT. MN) K=K-MN
            IJ = IJ + J - 1
            V1(L) = V1(L) - KM(IJ)*V1(K)
 80      CONTINUE

         DO 90 J = 1,MN
            K=J+IOLD-1
            IF (K .GT. MN) K=K-MN
            IJ = IJ + MN + J - 1
            V1(L) = V1(L) - KM(IJ)*V2(K)
 90      CONTINUE

         V1(L)=V1(L)/KM(II)
         II = II - I
 70   CONTINUE

      RETURN
      END


************************************************************************
*
*     * SUBROUTINE AGBFGS *
*
*      
*     * Purpose *
*      
*     Computation of aggregate values by the limited memory BFGS update.
*
*      
*     * Calling sequence *
*     
*     CALL AGBFGS(N,MN,IOLD,PGMHG,PGMHGM,PGNRM,G,GM,GA,PG,SM,UM,
*    &     RM,CDIAG,UMTUM,GAMMA,BETA,AGBETA,VMC1,VMC2,SMALL)
*     
*     
*     * Parameters *
*
*     II  N               Number of variables.
*     II  MN              Current number of stored corrections.
*     II  IOLD            Index for circular arrays.
*     RI  G(N)            Current (auxiliary) subgradient of the
*                           objective function.
*     RI  GM(N)           Previous subgradient and also aggregate
*                           subradient of the objective function
*     RO  GA(N)           Next aggregate subgradient of the objective
*                           function.
*     RI  PG(N)           Projected (auxiliary) subgradient of the 
*                           objective function.
*     RI  PGMHG           PGMHG = - trans(PGM)*H*PG, where PGM is a 
*                           projection of GM and H is the inverse
*                           approximation of the Hessian calculated
*                           by the L-BFGS formula.
*     RI  PGMHGM          PGMHGM = trans(PGM)*H*PGM.
*     RI  PGNRM           PGNRM = trans(PG)*PG.
*     RI  SM(N*MN)        Matrix whose columns are stored corrections.
*     RI  UM(N*MN)        Matrix whose columns are stored subgradient
*                           differences.
*     RI  RM(MN*(MN+1)/2) Upper triangular matrix.
*     RI  CDIAG(MN)       Diagonal matrix.
*     RI  UMTUM(MN*(MN+1)/2)  Matrix UMTUM = trans(UM)*UM.
*     RI  GAMMA           Scaling parameter.
*     RI  BETA            Locality measure.
*     RO  AGBETA          Aggregate locality measure.
*     RA  VMC#(MN)        Auxiliary arrays; # = 1,2.
*     RI  SMALL           Small positive value.
*     
*     
*     * Local variables *
*      
*     R   P               P = trans(PGM)*H*(PG-PGM) - BETA.
*     R   Q               Q = trans(PG-PGM)*H*(PG-PGM).
*     R   LAM             Multiplier used to calculate aggregate
*                            values.
*     I   IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -3  - Error in TRLIEQ.
*      
*     
*     * Subprograms used *
*      
*     S   SYMAX           Multiplication of a dense symmetric matrix
*                           by a vector.
*     S   RWAXV2          Multiplication of two rowwise stored dense 
*                           rectangular matrices A and B by vectors X 
*                           and Y.
*     S   TRLIEQ          Solving X from linear equation L*X=Y or
*                           trans(L)*X=Y.
*
*     RF  VDOT            Dot product of two vectors.
*     
*
*     The variable and subgradient differences and all the MN-vectors are
*     stored in a circular order controlled by the parameter point IOLD.
*
*      
      
      SUBROUTINE AGBFGS(N,MN,IOLD,PGMHG,PGMHGM,PGNRM,G,GM,GA,PG,
     &     SM,UM,RM,CDIAG,UMTUM,GAMMA,BETA,AGBETA,VMC1,VMC2,SMALL)

*     Scalar Arguments
      INTEGER N,MN,IOLD
      DOUBLE PRECISION BETA,AGBETA,PGMHG,PGMHGM,PGNRM,GAMMA,SMALL

*     Array Arguments
      DOUBLE PRECISION G(*),GM(*),GA(*),PG(*),SM(*),UM(*),RM(*),
     &     CDIAG(*),UMTUM(*),VMC1(*),VMC2(*)

*     Local Scalars
      INTEGER I,IERR
      DOUBLE PRECISION P,Q,LAM,SMALL2

*     External Functions
      DOUBLE PRECISION VDOT
      EXTERNAL VDOT

*     Intrinsic Functions
      INTRINSIC MAX,MIN,SIGN,SQRT

*     External Subroutines
      EXTERNAL RWAXV2,SYMAX,TRLIEQ
      

      SMALL2=SQRT(SMALL)
      IERR = 0


*     
*     Computation of P = trans(PGM)*DM*(PG-PGM)
*

      P = PGMHG - PGMHGM
      Q = PGNRM


*     
*     Computation of the product trans(PG)*DM*PG
*

      IF (MN .GT. 0) THEN

         CALL RWAXV2(N,MN,SM,UM,PG,PG,VMC1,VMC2)

         CALL TRLIEQ(MN,MN,IOLD,RM,VMC1,VMC1,1,SMALL,IERR)
            
         Q = Q - 2*VDOT(MN,VMC2,VMC1)
         Q = GAMMA*Q

         DO 10 I=1,MN
            VMC2(I) = CDIAG(I)*VMC1(I)
 10      CONTINUE

         Q = Q + VDOT(MN,VMC1,VMC2)

         CALL SYMAX(MN,MN,IOLD,UMTUM,VMC1,VMC2)

         Q = Q + GAMMA*VDOT(MN,VMC1,VMC2)

      END IF


*     
*     Computation of the product trans(PG-PGM)*DM*(PG-PGM)
*
      
      Q = Q + PGMHG - P
      P = P - BETA

      LAM = 0.5D+00 + SIGN(0.5D+00,P)
      IF (Q .GE. SMALL2) LAM = MIN(1.0D+00,MAX(0.0D+00,P/Q))
      

*
*     Computation of the aggregate values
*      

      DO 30 I=1,N
         GA(I)=LAM*G(I) + (1.0D+00 - LAM)*GM(I)
 30   CONTINUE
      
      AGBETA = LAM*BETA
      
      RETURN
      END


************************************************************************
*
*     * SUBROUTINE UPNULL *
*
*      
*     * Purpose *
*      
*     Matrix update by the limited memory SR1 update and calculation 
*     of the vector HPG = -H*PGA and the scalar PGTHPG = PGA'*H*PGA, 
*     where H denotes the L-SR1 approximation of the inverse of the 
*     Hessian matrix and PGA denotes the simple projection of 
*     subgradient.
*
*      
*     * Calling sequence *
*     
*      CALL UPNULL(N,MC,MN,INEW,IOLD,IFLAG,ISR1,NNK,HPG,PGA,PGM,
*     &     GA,GM,S,U,SM,UM,LM,RM,CDIAG,SMTSM,UMTUM,SMTGA,UMTGA,SMTPGA,
*     &     UMTPGA,SMTGM,UMTGM,SMTPGM,UMTPGM,VN1,VN2,VMC3,VMC4,VMC5,
*     &     VMC6,TMPMAT,GAMMA,PGTHPG,AMUGAD,XBX,TTHETA,SMALL,IPRINT,IERR)
*     
*     
*     * Parameters *
*
*     II  N               Number of variables.
*     II  MC              Declared number of stored corrections.
*     IU  MN              Current size of vectors.
*     IU  INEW            Index for circular arrays.
*     IO  IOLD            Index of the oldest corrections.
*     IU  IFLAG           Index for adaptive version:
*                             0  - Maximum number of stored corrections
*                                  has not been changed at previous 
*                                  iteration.
*                             1  - Maximum number of stored corrections
*                                  has been changed at previous iteration.
*     IO  ISR1            Index of the type of SR1 update.
*                             1  - SR1 update: the corrections are 
*                                  stored.
*                             3  - SR1 update is skipped.
*     II  NNK             Consecutive null steps counter.
*     RO  HPG(N)          Vector HPG = -H*PGA.
*     RI  GA(N)           Current aggregate subgradient.
*     RI  GM(N)           Basic subgradient of the objective function.
*     RI  PGA(N)          Projected aggregate subgradient.
*     RI  PGM(N)          Projected subgradient.
*     RI  S(N)            Difference of auxiliary and current variables.
*     RI  U(N)            Difference of auxiliary and current 
*                           subgradients.
*     RU  SM(N*MN)        Matrix whose columns are stored corrections.
*     RU  UM(N*MN)        Matrix whose columns are stored subgradient
*                           differences.
*     RU  LM(MN*(MN+1)/2) Lower triangular matrix stored rowwise in the
*                           one-dimensional array.
*     RU  RM(MN*(MN+1)/2) Upper triangular matrix stored columnwise
*                           in the one-dimensional array.
*     RU  CDIAG(MN)       Diagonal matrix.
*     RU  SMTSM(MN*(MN+1)/2)  Matrix SMTSM = SM'*SM.
*     RU  UMTUM(MN*(MN+1)/2)  Matrix UMTUM = UM'*UM.
*     RO  SMTGA(MN)       Vector SMTGA = SM'*GA.
*     RO  UMTGA(MN)       Vector UMTGA = UM'*GA.
*     RA  SMTPGA(MN)      Vector SMTPGA = SM'*PGA.
*     RA  UMTPGA(MN)      Vector UMTPGA = UM'*PGA.
*     RU  SMTGM(MN)       Vector SMTGM = SM'*GM.
*     RU  UMTGM(MN)       Vector UMTGM = UM'*GM.
*     RU  SMTPGM(MN)      Vector SMTPGA = SM'*PGM.
*     RU  UMTPGM(MN)      Vector UMTPGA = UM'*PGM.
*     RA  VN#(N)          Auxiliary arrays: # = 1,2.
*     RA  VMC#(MN)        Auxiliary arrays: # = 3,...,6.
*     RO  TMPMAT((MN+1)*(MN)/2)  Auxiliary matrix.
*                           On output: Factorization of matrix 
*                           GAMMA*UMTUM-RM-RM'+CDIAG
*     RU  GAMMA           Scaling parameter.
*     RO  PGTHPG          PGTHPG = PGA'*H*PGA.
*     RI  AMUGAD          AMUGAD = -(A*MU + G)'* D, where A*MU denotes
*                           the Lagrange multipliers for the problem.
*     RI  XBX             XBX = (XCP-X)'*BM*(XCP-X), where BM denotes the 
*                           limited memory approximation of the Hessian 
*                           matrix.
*     RI  TTHETA          Scaled stepsize: TTHETA = T * THETA.
*     RI  SMALL           Small positive value.
*     II  IPRINT          Printout specification.
*     IO  IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -2  - Error in CALQ.
*                           -11  - Warning: MN=0 and ISR1=3.
*     
*     
*     * Local variables *
*
*     R   STU             STU = S'*U. 
*
*      
*     * Subprograms used *
*      
*     S   CALQ            Solving X from linear equation A*X=Y.
*     S   COPY            Copying of a vector.
*     S   COPY2           Copying of two vectors.
*     S   INDIC2          Initialization of indices.
*     S   MXDPGF          Gill-Murray decomposition of a dense symmetric
*                           matrix.
*     S   RECMAX          Multiplication of a vector by a dense
*                           rectangular matrix.
*     S   RWAXV2          Multiplication of two rowwise stored dense 
*                           rectangular matrices A and B by vectors X 
*                           and Y.
*     S   SCALEX          Scaling a vector.
*     S   SCDIFF          Difference of the scaled vector and a vector.
*     S   SCSUM           Sum of the scaled vector and a vector.
*     S   SR1XV           Computation of the product H*V or -H*V, 
*                           where H is the inverse approximation of 
*                           the Hessian calculated by the L-SR1 
*                           formula and V is an arbitrary vector.
*     S   UPDATE          Updating the matrices RM, LM, CDIAG, SMTSM,
*                            and UMTUM for LMBM-B.
*     S   XDIFFY          Difference of two vectors.
*     S   XSUMY           Sum of two vectors.
*
*     RF  VDOT            Dot product of two vectors.
*     
*
*     The variable and subgradient differences and all the MN-vectors are
*     stored in a circular order controlled by the parameter point INEW.
*
*

      SUBROUTINE UPNULL(N,MC,MN,INEW,IOLD,IFLAG,ISR1,NNK,HPG,PGA,PGM,
     &     GA,GM,S,U,SM,UM,LM,RM,CDIAG,SMTSM,UMTUM,SMTGA,UMTGA,SMTPGA,
     &     UMTPGA,SMTGM,UMTGM,SMTPGM,UMTPGM,VN1,VN2,VMC3,VMC4,VMC5,
     &     VMC6,TMPMAT,GAMMA,PGTHPG,AMUGAD,XBX,TTHETA,SMALL,IPRINT,IERR)

*     Scalar Arguments
      INTEGER N,MC,MN,INEW,IOLD,IFLAG,ISR1,NNK,IPRINT,IERR
      DOUBLE PRECISION GAMMA,PGTHPG,TTHETA,AMUGAD,XBX,SMALL

*     Array Arguments
      DOUBLE PRECISION HPG(*),GA(*),PGA(*),PGM(*),GM(*),S(*),U(*),
     &     SM(*),UM(*),UMTUM(*),SMTSM(*),LM(*),RM(*),CDIAG(*),SMTGA(*),
     &     UMTGA(*),SMTPGA(*),UMTPGA(*),SMTGM(*),UMTGM(*),SMTPGM(*),
     &     UMTPGM(*),VN1(*),VN2(*),VMC3(*),VMC4(*),VMC5(*),VMC6(*),
     &     TMPMAT(*)

*     Local Scalars
      INTEGER I,J,K,IFULL,INF
      DOUBLE PRECISION STU,STGA,UTGA,STPGA,UTPGA,GHOLDG,ETA,BET

*     Intrinsic Functions
      INTRINSIC MAX,SQRT
      
*     External Functions
      DOUBLE PRECISION VDOT
      EXTERNAL VDOT

*     External Subroutines
      EXTERNAL CALQ,COPY,COPY2,INDIC2,MXDPGF,RECMAX,RWAXV2,SCALEX,
     &     SCDIFF,SCSUM,SR1XV,UPDATE,XDIFFY,XSUMY


      IFULL = 0
      ISR1 = 0
      IERR = 0
      INF = 0
      ETA = SQRT(SMALL) + SMALL


*
*     Initialization of indices
*      

      IF (MN .EQ. MC) IFULL = 1
      CALL INDIC2(MC,MN,INEW,IOLD,IFLAG,3)


*         
*     Computation of GAMMA 
*

      GAMMA = 1.0D+00

      
*     
*     Computation of SM'*GA, UM'*GA, SM'*PGA, and UM'*PGA
*

      CALL RWAXV2(N,MN,SM,UM,GA,GA,SMTGA,UMTGA)
      CALL RWAXV2(N,MN,SM,UM,PGA,PGA,SMTPGA,UMTPGA)

      
*
*     Positive definiteness
*

      STU = VDOT(N,S,U)

      IF (-STU + TTHETA*TTHETA*MAX(XBX,AMUGAD) .GE. -SQRT(SMALL) .OR.
     $     STU .LT. SQRT(SMALL)) THEN


*     
*     SR1 update is skipped
*     

         ISR1 = 3

         IF (MN .EQ. 0) THEN

            IERR = -11
            RETURN
         END IF

         GO TO 200

      END IF

      STGA = VDOT(N,S,GA)
      UTGA = VDOT(N,U,GA)
      STPGA = VDOT(N,S,PGA)
      UTPGA = VDOT(N,U,PGA)

     
*
*     Convergence conditions
*
      IF (NNK.EQ.1 .OR. ((MN.LT.MC) .AND. (IFLAG.EQ.0))) GO TO 100


*
*     Calculation of matrix TMPMAT = GAMMA*UMTUM-RM-RM'+CDIAG) from 
*     previous iteration
*      

      DO 10 I=1,MN*(MN+1)/2
         TMPMAT(I) = GAMMA * UMTUM(I) - RM(I)
 10   CONTINUE


*
*     Computation of vector VN2 = -H_(k-1)*PGA and 
*     GHOLDG = -PGA'*H_(k-1)*PGA
*

      CALL SR1XV(N,IOLD,MN,VN2,PGA,SM,UM,TMPMAT,VN1,SMTPGA,
     &     UMTPGA,VMC3,VMC4,GAMMA,SMALL,1,IERR)
      IF (IERR .NE. 0) RETURN

      GHOLDG = VDOT(N,PGA,VN2)


*
*     Calculation of the new canditate for updating matrix.
*     Updates are not necessarily saved
*

      IF (IFLAG .EQ. 0) THEN
         IOLD = INEW + 1
         IF (IOLD .GT. MC) IOLD = 1
         
      ELSE
         IOLD = INEW + 1
         IF (IOLD .GT. MC-1) IOLD = 1
         
      END IF

      
*     
*     Computation of SM'*U and UM'*U
*

      IF (IOLD .EQ. 1) THEN
         CALL RWAXV2(N,MN-1,SM,UM,U,U,VMC3,VMC4)

      ELSE
         CALL RWAXV2(N,INEW-1,SM,UM,U,U,VMC3,VMC4)
         CALL RWAXV2(N,MN-INEW,SM((IOLD-1)*N+1),UM((IOLD-1)*N+1),U,U,
     &        VMC3(IOLD),VMC4(IOLD))
         
      END IF

      VMC3(INEW) = STU
      VMC4(INEW) = VDOT(N,U,U)


*
*     Calculation of matrix TMPMAT=GAMMA*UMTUM-RM-RM'+CDIAG without 
*     updating matrices RM, UMTUM and CDIAG
*

      DO 20 I=1,MN*(MN+1)/2
         TMPMAT(I)= GAMMA * UMTUM(I) - RM(I)
 20   CONTINUE

      DO 30 I=1,MN-1
         J=(I-1)*I/2+1
         K=I*(I+1)/2+2
         CALL COPY(I,TMPMAT(K),TMPMAT(J))
 30   CONTINUE


      CALL SCDIFF(MN,GAMMA,VMC4,VMC3,VMC5)

      IF (IOLD .EQ. 1) THEN
         CALL COPY(MN,VMC5,TMPMAT((MN-1)*MN/2+1))

      ELSE
         CALL COPY(MN-INEW,VMC5(IOLD),
     &        TMPMAT((MN-1)*MN/2+1))
         CALL COPY(INEW,VMC5,
     &        TMPMAT((MN-1)*MN/2+MN-INEW+1))

      END IF

            
      CALL SCDIFF(MN,GAMMA,UMTPGA,SMTPGA,VMC5)
      VMC5(INEW) = GAMMA*UTPGA-STPGA


      CALL CALQ(MN,MN,IOLD,TMPMAT,VMC5,VMC5,SMALL,IPRINT,
     &     IERR)
      IF (IERR .NE. 0) RETURN

     
*
*     Calculation of the new canditate for vector HPG = -H*PGA
*     and the scalar PGTHPG = PGA'*H*PGA
*

      IF (IOLD .EQ. 1) THEN
         CALL SCALEX(MN,GAMMA,VMC5,VMC6)
         CALL RECMAX(N,MN-1,SM,VMC5,VN1)
         CALL SCSUM(N,VMC5(INEW),S,VN1,VN1)
         CALL SCDIFF(N,-GAMMA,PGA,VN1,HPG)
         CALL RECMAX(N,MN-1,UM,VMC6,VN1)
         CALL SCSUM(N,VMC6(INEW),U,VN1,VN1)
         CALL XSUMY(N,HPG,VN1,HPG)
         
      ELSE
         CALL SCALEX(MN,GAMMA,VMC5,VMC6)
         CALL RECMAX(N,INEW-1,SM,VMC5,VN1)
         CALL SCSUM(N,VMC5(INEW),S,VN1,VN1)
         CALL SCDIFF(N,-GAMMA,PGA,VN1,HPG)
         CALL RECMAX(N,MN-INEW,SM((IOLD-1)*N+1),VMC5(IOLD),VN1)
         CALL XDIFFY(N,HPG,VN1,HPG)
         CALL RECMAX(N,INEW-1,UM,VMC6,VN1)
         CALL SCSUM(N,VMC6(INEW),U,VN1,VN1)
         CALL XSUMY(N,HPG,VN1,HPG)
         CALL RECMAX(N,MN-INEW,UM((IOLD-1)*N+1),VMC6(IOLD),VN1)
         CALL XSUMY(N,HPG,VN1,HPG)
      END IF

      PGTHPG = -VDOT(N,PGA,HPG)

*
*     Checking the convergence conditions
*      

      IF (-PGTHPG - GHOLDG .LT. 0.0D+00) THEN

         ISR1 = 3

         CALL COPY(N,VN2,HPG)
         PGTHPG = - GHOLDG


*
*     Calculation and factorization of matrix TMPMAT =
*     GAMMA*UMTUM-RM-RM'+CDIAG from  previous iteration
*      
         
         DO 40 I=1,MN*(MN+1)/2
            TMPMAT(I) = GAMMA * UMTUM(I) - RM(I)
 40      CONTINUE

         CALL MXDPGF(MN,TMPMAT,INF,ETA,BET)
         IOLD = INEW

      ELSE

         ISR1 = 1
         

*     
*     Update SM and UM
*
         CALL COPY2(N,S,SM((INEW-1)*N+1),U,UM((INEW-1)*N+1))
         

*     
*     Update SM'*GM, UM'*GM, SM'*PGM, and UM'*PGM
*  

         SMTGM(INEW) = VDOT(N,S,GM)
         UMTGM(INEW) = VDOT(N,U,GM)
         SMTPGM(INEW) = VDOT(N,S,PGM)
         UMTPGM(INEW) = VDOT(N,U,PGM)


*     
*     Update SM'*GA, UM'*GA, SM'*PGA, and UM'*PGA
*     

         SMTGA(INEW) = STGA
         UMTGA(INEW) = UTGA
         SMTPGA(INEW) = STPGA
         UMTPGA(INEW) = UTPGA

            
*     
*     Computation of SM'*S, and UM'*S
*

         CALL RWAXV2(N,MN,SM,UM,S,S,VMC5,VMC6)


*         
*     Update RM, LM, SMTSM, UMTUM, and CDIAG
*

         CALL UPDATE(MN,INEW,IOLD,RM,LM,CDIAG,SMTSM,UMTUM,VMC3,
     &        VMC4,VMC5,VMC6,STU,IFLAG,IFULL)

         INEW = INEW + 1
         IF (INEW .GT. MC) INEW = 1

      END IF


      GO TO 300

      
 100  CONTINUE

      ISR1 = 1

      CALL INDIC2(MC,MN,INEW,IOLD,IFLAG,1)        

      SMTGA(INEW) = STGA
      UMTGA(INEW) = UTGA
      SMTPGA(INEW) = STPGA
      UMTPGA(INEW) = UTPGA

            
*     
*     Update SM and UM
*

      CALL COPY2(N,S,SM((INEW-1)*N+1),U,UM((INEW-1)*N+1))

            
*     
*     Update SM'*GM, UM'*GM, SM'*PGM, and UM'*PGM
*  

      SMTGM(INEW) = VDOT(N,S,GM)
      UMTGM(INEW) = VDOT(N,U,GM)
      SMTPGM(INEW) = VDOT(N,S,PGM)
      UMTPGM(INEW) = VDOT(N,U,PGM)

      
*     
*     Computation of SM'*U, UM'*U, SM'*S, and UM'*S.
*

      IF (IOLD .EQ. 1) THEN
         CALL RWAXV2(N,MN-1,SM,UM,U,U,VMC3,VMC4)
         CALL RWAXV2(N,MN,SM,UM,S,S,VMC5,VMC6)

      ELSE
         CALL RWAXV2(N,INEW-1,SM,UM,U,U,VMC3,VMC4)
         CALL RWAXV2(N,MN-INEW,SM((IOLD-1)*N+1),UM((IOLD-1)*N+1),
     &        U,U,VMC3(IOLD),VMC4(IOLD))
         CALL RWAXV2(N,MN,SM,UM,S,S,VMC5,VMC6)

      END IF

      VMC3(INEW) = STU
      VMC4(INEW) = VDOT(N,U,U)

      
*         
*     Update RM, LM, SMTSM, UMTUM, and CDIAG
*

      CALL UPDATE(MN,INEW,IOLD,RM,LM,CDIAG,SMTSM,UMTUM,VMC3,
     &     VMC4,VMC5,VMC6,STU,IFLAG,IFULL)

         
      INEW = INEW + 1
      IF (INEW .GT. MC) INEW = 1
      
      
 200  CONTINUE


*
*     Calculation of matrix TMPMAT = GAMMA*UMTUM - RM - RM' + CDIAG
*      

      DO 210 I=1,MN*(MN+1)/2
         TMPMAT(I)= GAMMA * UMTUM(I) - RM(I)
 210  CONTINUE


*
*     Calculation of the vector HPG = -H*PGA and the scalar 
*     PGTHPG = PGA'*H*PGA
*

      CALL SR1XV(N,IOLD,MN,HPG,PGA,SM,UM,TMPMAT,VN1,SMTPGA,
     &     UMTPGA,VMC3,VMC4,GAMMA,SMALL,1,IERR)
      IF (IERR .NE. 0) GOTO 300


      PGTHPG = -VDOT(N,PGA,HPG)

 300  CONTINUE

      
      RETURN
      END
      

************************************************************************
*
*     * SUBROUTINE SR1XV *
*
*      
*     * Purpose *
*      
*     Computation of the product H*V or -H*V, where H is the inverse
*     approximation of the Hessian calculated by the L-SR1 formula and
*     V is an arbitrary vector.
*
*      
*     * Calling sequence *
*     
*     CALL SR1XV(N,IOLD,MN,Z,V,SM,UM,TMPMAT,VN1,SMTV,UMTV,VMC3,VMC4,
*    &      GAMMA,SMALL,JOB,IERR)
*     
*     
*     * Parameters *
*
*     II  N               Number of variables.
*     II  IOLD            Index of the oldest corrections.
*     II  MN              Current size of vectors.
*     RI  V(N)            Input vector.
*     RO  Z(N)            Output vector Z = H*V or Z = -H*V.
*     RI  SM(N*MN)        Matrix whose columns are stored corrections.
*     RI  UM(N*MN)        Matrix whose columns are stored subgradient
*                           differences.
*     RI  SMTV(MN)        Vector SMTV = SM'*V.
*     RI  UMTV(MN)        Vector UMTV = UM'*V.
*     RI  GAMMA           Scaling parameter.
*     RI  TMPMAT((MN+1)*(MN)/2)  Factorization of matrix 
*                           UMTUM-RM-RM'+CDIAG.
*     RA  VN1(N)          Auxiliary array.
*     RA  VMC#(MN)        Auxiliary arrays: # = 3,4.
*     RI  SMALL           Small positive value.
*     II  JOB             Selection of the sign:
*                             0  - Z = H*V.
*                             1  - Z = -H*V.
*     IO  IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -2  - Error in LINEQ.
*     
*     
*     * Subprograms used *
*      
*     S   LINEQ           Solving X from linear equation A*X=Y.
*     S   RECMAX          Multiplication of a columnwise stored dense 
*                           rectangular matrix by a vector.
*     S   SCALEX          Scaling a vector.
*     S   SCDIFF          Difference of the scaled vector and a vector.
*     S   XDIFFY          Difference of two vectors.
*     S   XSUMY           Sum of two vectors.
*     
*
*     The variable and subgradient differences and all the MN-vectors
*     are stored in a circular order controlled by the parameter point
*     IOLD.
*
*     
  
      SUBROUTINE SR1XV(N,IOLD,MN,Z,V,SM,UM,TMPMAT,VN1,SMTV,
     &     UMTV,VMC3,VMC4,GAMMA,SMALL,JOB,IERR)

*     Scalar Arguments
      INTEGER N,IOLD,MN,IERR,JOB
      DOUBLE PRECISION GAMMA,SMALL
      
*     Array Arguments
      DOUBLE PRECISION Z(*),V(*),SM(*),UM(*),TMPMAT(*),VN1(*),SMTV(*),
     &     UMTV(*),VMC3(*),VMC4(*)
      
*     External Subroutines
      EXTERNAL RECMAX,LINEQ,SCALEX,SCDIFF,XDIFFY,XSUMY


      IERR = 0

      CALL SCDIFF(MN,GAMMA,UMTV,SMTV,VMC4)
      CALL CALQ(MN,MN,IOLD,TMPMAT,VMC4,VMC4,SMALL,0,IERR)
      IF (IERR .NE. 0) RETURN

      IF (JOB .EQ. 1) THEN

*
*     Computation of the product Z = -H*V
*

         CALL SCALEX(MN,GAMMA,VMC4,VMC3)

*     SM*(TMPMAT)^-1 * (GAMMA*UM'*V - SM'*V)

         CALL RECMAX(N,MN,SM,VMC4,VN1)

*     -GAMMA*V - SM*(TMPMAT)^-1 * (GAMMA*UM'*V - SM'*V)

         CALL SCDIFF(N,-GAMMA,V,VN1,Z)

*     GAMMA*UM*(TMPMAT)^-1 * (GAMMA*trans(UM)*V - trans(SM)*V)

         CALL RECMAX(N,MN,UM,VMC3,VN1)
         CALL XSUMY(N,Z,VN1,Z)

         
      ELSE


*
*     Computation of the product Z = H*V
*

         CALL SCALEX(MN,GAMMA,VMC4,VMC3)

*     SM*(TMPMAT)^-1 * (GAMMA*trans(UM)*V - trans(SM)*V)

         CALL RECMAX(N,MN,SM,VMC4,VN1)

*     GAMMA*V + SM*(TMPMAT)^-1 * (GAMMA*trans(UM)*V - trans(SM)*V)

         CALL SCSUM(N,GAMMA,V,VN1,Z)

*     GAMMA*UM*(TMPMAT)^-1 * (GAMMA*trans(UM)*V - trans(SM)*V)

         CALL RECMAX(N,MN,UM,VMC3,VN1)
         CALL XDIFFY(N,Z,VN1,Z)
         
      END IF

      
 300  CONTINUE
      
      RETURN
      END
      

************************************************************************
*
*     * SUBROUTINE CPSR1 *
*
*      
*     * Purpose *
*      
*     Computation of generalized Cauchy point by the limited memory
*     SR1 update.
*
*      
*     * Calling sequence *
*     
*     CALL CPSR1(N,MC,MN,IOLD,IFLAG,IPRINT,X,XCP,XL,XU,IB,IT,GA,XBX,
*    &     DC,T,SM,UM,SMTGA,UMTGA,SMTSM,LM,CDIAG,TMPM2,PAUX,CAUX,
*    &     VAUX,WAUX,GAMMA,SMALL,IERR)
*
*      
*     * Parameters *
*     
*     II  N               Number of variables.
*     II  MN              Current number of stored corrections.
*     II  IOLD            Index for circular arrays.
*     II  IFLAG           Index for adaptive version:
*                             0  - Maximum number of stored corrections
*                                    has not been changed at previous
*                                    iteration.
*                             1  - Maximum number of stored corrections
*                                    has been changed at previous
*                                    iteration.
*     RI  X(N)            Vector of variables. 
*     RO  XCP(N)          Generalized Cauchy point.
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables. 
*     II  IB(N)           Type of bound constraints:
*                             0  - X(I) is unbounded,
*                             1  - X(I) has only a lower bound,
*                             2  - X(I) has both lower and upper bounds,
*                             3  - X(I) has only an upper bound. 
*     RI  GA(N)           Current aggregate subgradient of the objective
*                           function.
*     RO  XBX             XBX = (XCP-X)'B(XCP-X), where B denotes the 
*                           limited memory approximation of the Hessian 
*                           matrix.
*     RI  SM(N*MN)        Matrix whose columns are stored corrections.
*     RI  UM(N*MN)        Matrix whose columns are stored subgradient
*                           differences.
*     RI  SMTGA(MN)       Vector SMTGA = SM'*GA.
*     RI  UMTGA(MN)       Vector UMTGA = UM'*GA.
*     RI  SMTSM(MN*(MN+1)/2)  Matrix SMTSM = SM'*SM.
*     RI  LM(MN*(MN+1)/2) Lower triangular matrix stored rowwise
*                           in the one-dimensional array.
*     RI  GAMMA           Scaling parameter.
*     RI  CDIAG(MN)       Diagonal matrix.
*     RA  TMPM2(MN*(MN+1)/2)  Matrix TMPM2= LM+LM'+CDIAG-1/GAMMA*SMTSM 
*                           or its factorization.
*     RA  DC(N)           Auxiliary array used to store Cauchy
*                           direction P(X-T*GA)-X.
*     RA  T(N)            Auxiliary array used to store the
*                           breakpoints.
*     IA  IT(N)           Index set of breakpoints:
*                             for I=1,...,NLEFT, IT(I) are the indices
*                                 of breakpoints that have not been
*                                 encountered,
*                             for i=NLEFT+1,...,NBREAK, IT(I) are the
*                                 indices of encountered breakpoints and
*                             for i=IFREE,...,N, IT(I) are the indices
*                                 of variables that have no bound
*                                 constraints along the search direction.
*     RA  PAUX(MN)        Auxiliary array used to store the vector
*                           (UM-1/GAMMA*(SM))'*DC.
*     RA  CAUX(MN)        Auxiliary array used to store the vector
*                           (UM-1/GAMMA*(SM))'*(XCP-X).
*     RA  VAUX(MN)        Auxiliary array.
*     RA  WAUX(MN)        Auxiliary array.
*     RI  SMALL           Small positive value.
*     II  IPRINT          Printout specification.
*     IO  IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -2  - Error in LINEQ or in CALQ.
*                            -7  - Error in CPSR1, a non positive
*                                  definite matrix detected.
*     
*
*     * Subprograms used *
*      
*     RF  VDOT            Dot product of two vectors.
*     S   SCSUM           Sum of a vector and the scaled vector.
*     S   HPSRT           Heapsort algorithm.
*     S   CALQ            Solving X from linear equation A*X=Y.
*     S   LINEQ           Solving X from linear equation A*X=Y, where
*                           A is given in factorized form.
*
*     
*     The variable and subgradient differences and all the input
*     MN-vectors are stored in a circular order controlled by the
*     parameter point IOLD.
*
*
      
      SUBROUTINE CPSR1(N,MN,IOLD,IFLAG,IPRINT,X,XCP,XL,XU,IB,IT,GA,XBX,
     &     DC,T,SM,UM,SMTGA,UMTGA,SMTSM,LM,CDIAG,TMPM2,PAUX,CAUX,VAUX,
     &     WAUX,GAMMA,SMALL,IERR)
      
      
*     Scalar Arguments
      INTEGER N,MN,IOLD,IFLAG,IPRINT,IERR
      DOUBLE PRECISION XBX,GAMMA,SMALL
      
*     Array Arguments
      INTEGER IB(*),IT(*)
      DOUBLE PRECISION GA(*),X(*),XCP(*),XL(*),XU(*),DC(*),T(*),
     &     SM(*),UM(*),LM(*),CDIAG(*),SMTSM(*),TMPM2(*),
     &     SMTGA(*),UMTGA(*),PAUX(*),CAUX(*),VAUX(*),WAUX(*)
      
*     Local Scalars
      INTEGER I,II,J,L,NBREAK,IFREE,ITMIN,NINT,NLEFT,ISBP,IBOUND
      DOUBLE PRECISION F1,F2,TMIN,DLTMIN,TOLD,DELTAT,GISBP,GISBP2,
     &     ZISBP,WMC,WMP,WMW,SMALL2,TMP
      
*     External Functions
      DOUBLE PRECISION VDOT
      EXTERNAL VDOT

*     Intrinsic Functions
      INTRINSIC MAX,ABS,SQRT

*     External Subroutines
      EXTERNAL CALQ,LINEQ,HPSRT,SCSUM
      

*
*     Initialization
*      

      NBREAK = 0
      NLEFT  = 0
      IFREE  = N+1
      ITMIN  = 0
      IBOUND = 0
      F1     = 0.0D+00
      F2     = 0.0D+00
      XBX    = 0.0D+00
      TMIN   =  0.0D+00
      IERR = 0

c      SMALL2=1.0D+01*SMALL
      SMALL2=SQRT(SMALL)

      
*     
*     Temporarily set PAUX = -(UM - 1/GAMMA*SM)'*GA
*     and CAUX = (UM - 1/GAMMA*SM)'*(XCP-X) = 0
*      

      DO 10 I=1,MN
         L = I+IOLD-1
         IF (L .GT. MN) L=L-MN
         PAUX(I) = -UMTGA(L)+SMTGA(L)/GAMMA
         CAUX(I)= 0.0D+00
 10   CONTINUE


*
*     Computation of the Cauchy direction DC, the first derivative
*     F1= -(DC)'*DC, and the breakpoints T in each coordinate
*     direction. Identification of the smallest breakpoint TMIN. 
*      

      DO 20 I=1,N

         IF (IB(I) .NE. 0) THEN
            IF (IB(I).LE.2 .AND. ABS(GA(I)) .LE. SMALL2 .AND. 
     &           X(I)-XL(I) .LE. SMALL2) THEN
               DC(I) = 0.0D+00
               XCP(I) = X(I)
               GOTO 20
            END IF

            IF (IB(I).GE.2 .AND. ABS(GA(I)) .LE. SMALL2 .AND. 
     &           XU(I)-X(I) .LE. SMALL2) THEN
               DC(I) = 0.0D+00
               XCP(I) = X(I)
               GOTO 20
            END IF
         END IF

         IF (IB(I) .NE. 0 .AND. IB(I) .LE. 2 .AND. GA(I) .GT. 
     &        0.0D+00) THEN

*     Breakpoint in lower bound

            T(NBREAK+1) = (X(I) - XL(I))/GA(I)
            
            IF (T(NBREAK+1) .LE. SMALL2) THEN
               DC(I) = 0.0D+00
               
*     Correction to PAUX = (UM - 1/GAMMA*SM)'*DC 

               DO 30 J=1,MN
                  L = J + IOLD - 1
                  IF (L .GT. MN) L=L-MN
                  PAUX(J) = PAUX(J) +
     &                 (UM(N*(L-1)+I)-SM(N*(L-1)+I)/GAMMA)*GA(I)
 30            CONTINUE

            ELSE
               NBREAK = NBREAK + 1
               IT(NBREAK) = I
               DC(I) = -GA(I)
               F1 = F1 - DC(I)*DC(I)
     
*     Determination of the smallest breakpoint
               
               IF (NBREAK .EQ. 1 .OR. T(NBREAK) .LT. TMIN) THEN
                  TMIN=T(NBREAK)
                  ITMIN=NBREAK
               END IF
            END IF

         ELSE IF (IB(I) .GE. 2 .AND. GA(I) .LT. 0.0D+00) THEN

*     Breakpoint in upper bound
               
            T(NBREAK+1) = (X(I) - XU(I))/GA(I)
            
            IF (T(NBREAK+1) .LE. SMALL2) THEN
               DC(I) = 0.0D+00
     
*     Correction to P = (UM - 1/GAMMA*SM)'*DC 

               DO 40 J=1,MN
                  L = J+IOLD-1
                  IF (L .GT. MN) L=L-MN
                  PAUX(J) = PAUX(J) +
     &                 (UM(N*(L-1)+I)-SM(N*(L-1)+I)/GAMMA)*GA(I)
 40            CONTINUE

            ELSE
               NBREAK = NBREAK + 1
               IT(NBREAK) = I
               DC(I) = -GA(I)
               F1 = F1 - DC(I)*DC(I)

*     Determination of the smallest breakpoint
               
               IF (NBREAK .EQ. 1 .OR. T(NBREAK) .LT. TMIN) THEN
                  TMIN=T(NBREAK)
                  ITMIN=NBREAK
               END IF
            END IF
            
         ELSE

*     No breakpoint 

            IFREE = IFREE - 1
            IT(IFREE) = I
            DC(I) = -GA(I)
            F1 = F1 - DC(I)*DC(I) 
            IF (ABS(GA(I)) .GT. 0.0D+00) IBOUND=-1

         END IF


*     
*     Initialization of Cauchy point
*

         XCP(I) = X(I)
           

 20   CONTINUE


*
*     The indices corresponding to free variables are located in
*     IT(1),...,IT(NBREAK) and IT(IFREE),...,IT(N). The smallest of the
*     NBREAK breakpoints is in T(ITMIN)=TMIN.
*
      
*
*     DC is a zero vector. Return with the initial XCP.
*

      IF (NBREAK .EQ. 0 .AND. IBOUND .EQ. 0) GOTO 999


*
*     Initialization of derivative F2 = (DC)' BM DC, where BM is 
*     representated by compact limited memory SR1 formula.
*

      F2 =  -F1/GAMMA


*
*     Calculation of matrix (LM+LM'+CDIAG-SMTSM/GAMMA)
*      

      DO 50 I=1,MN*(MN+1)/2
         TMPM2(I)= LM(I) - SMTSM(I)/GAMMA
 50   CONTINUE

      II = 0
      DO 51 I = 1, MN
         II=II+I
         J=I+IOLD-1
         IF (J .GT. MN) J=J-MN
      	 TMPM2(II) = TMPM2(II) + CDIAG(J)
 51   CONTINUE

      
*     The product of the inverse of MN x MN middle matrix in the compact
*     L-SR1 formula of BM and a MN vector PAUX (PAUX is not in circular 
*     order); The product is returned in VAUX. Note that TMPM2 will be 
*     factorized during the process.

      CALL CALQ(MN,MN,1,TMPM2,VAUX,PAUX,SMALL,IPRINT,IERR)
      IF (IERR .NE. 0) RETURN

      F2 = F2 + VDOT(MN,PAUX,VAUX)

      IF (F2 .LE. SMALL2) THEN
         IERR = -7
         RETURN
      END IF

      DLTMIN = -F1/F2

      TOLD = 0.0D+00
      NINT = 1

      NLEFT = NBREAK
 
     
*     There is no breakpoints, locate the GCP and return. 
 
      IF (NBREAK .EQ. 0) GOTO 888

      
*
*     Begining of the loop.
*      

 100  CONTINUE


      IF (NINT .EQ. 1) THEN


*         
*     the smallest breakpoint is TMIN=T(ITMIN) and its index is
*     in IT(ITMIN)
*         

         ISBP = IT(ITMIN)

      ELSE


*         
*     the smallest breakpoint is chosen by heapsort algorithm
*

         IF (NINT .EQ. 2) THEN


*            
*     Remove the already used smallest breakpoint before heapsort call.
*            

            IF (ITMIN .NE. NBREAK) THEN
               T(ITMIN) = T(NBREAK)
               IT(ITMIN) = IT(NBREAK)
            END IF


*     
*     Heapsort with initialization of the heap.
*

            CALL HPSRT(NLEFT,T,IT,0)

         ELSE


*            
*     Heapsort with updating the heap.
*            

            CALL HPSRT(NLEFT,T,IT,1)

         END IF

         TMIN = T(NLEFT)
         ISBP = IT(NLEFT)

      END IF

      DELTAT = TMIN - TOLD
      
      IF (DELTAT .LT. 0.0D+00) THEN
         PRINT*,' no nyt voit jo hirttää ihtes.'
      END IF


*
*     Minimizer is within this interval, locate the GCP and return.
*

      IF (DLTMIN .LT. DELTAT) GOTO 888


*      
*     Examination of subsequent seqment
*

      NLEFT = NLEFT - 1
      NINT = NINT + 1
      GISBP = GA(ISBP)
      DC(ISBP) = 0.0D+00


      IF (GISBP .LT. 0.0D+00) THEN
         XCP(ISBP) = XU(ISBP)
         ZISBP = XU(ISBP) - X(ISBP)
      ELSE 
         XCP(ISBP) = XL(ISBP)
         ZISBP = XL(ISBP) - X(ISBP)
      END IF


*         
*     All  variables are fixed, return with current XCP as GCP.
*         

      IF (NLEFT .EQ. 0 .AND. NBREAK .EQ. N) GOTO 999


*      
*     Update F1 and F2
*      

      GISBP2=GISBP*GISBP


*
*     Initial updating of F1 and F2
*      

      F1 = F1 + DELTAT*F2 + GISBP2 + GISBP*ZISBP/GAMMA
      TMP = GISBP2/GAMMA
      F2 = F2 - TMP


*         
*     Update CAUX = CAUX + DELTAT*PAUX.
*         

      CALL SCSUM(MN,DELTAT,PAUX,CAUX,CAUX)


*     
*     Selection of WAUX, the ISBP:th row of (UM - 1/GAMMA*SM). 
*
         
      DO 110 J=1,MN
         L = J+IOLD-1
         IF (L .GT. MN) L=L-MN
         WAUX(J)= UM(N*(L-1)+ISBP)-SM(N*(L-1)+ISBP)/GAMMA
 110  CONTINUE


*
*     Computation of products WMC, WMP and WMW
*

*     The product of the inverse of MN x MN middle matrix in the compact
*     L-SR1 formula of BM and an MN vector WAUX; The product is returned
*     in VAUX. Note that TMPM2 is already in the factorized form and WAUX
*     is not in circular order.

      CALL LINEQ(MN,MN,1,TMPM2,VAUX,WAUX,SMALL,IERR)
      IF (IERR .NE. 0) RETURN
         
      WMC = VDOT(MN,CAUX,VAUX)
      WMP = VDOT(MN,PAUX,VAUX)
      WMW = VDOT(MN,WAUX,VAUX)


*         
*     Update PAUX = PAUX + GISBP*WAUX. 
*

      CALL SCSUM(MN,GISBP,WAUX,PAUX,PAUX)
      

*
*     Complete updating of F1 and F2
*      

      F1 = F1 + GISBP*WMC
      F2 = F2 + 2.0D+00*GISBP*WMP + GISBP2*WMW
      TMP = TMP - 2.0D+00*GISBP*WMP - GISBP2*WMW
      XBX = XBX + TMIN*TMIN * TMP

      IF (F2 .LE. SMALL2) THEN
c     DC is a zero vector. Return with current XCP.

         GOTO 999
      END IF

      DLTMIN = -F1/F2
      TOLD=TMIN


*     
*     End of the loop.
*

      IF (NLEFT .GT. 0)  GOTO 100
        
      IF (IBOUND .EQ. 0) DLTMIN = 0.0D+00   
      

 888  CONTINUE


      DLTMIN = MAX(0.0D+00,DLTMIN)
      TOLD = TOLD + DLTMIN

      XBX = XBX + TOLD*TOLD * F2


*      
*     Update free variables 
* 

      DO 60 I=1,NLEFT
         XCP(IT(I))=X(IT(I))+TOLD*DC(IT(I))
 60   CONTINUE
      DO 70 I=IFREE,N
         XCP(IT(I))=X(IT(I))+TOLD*DC(IT(I))
 70   CONTINUE

 999  CONTINUE
      
      RETURN
      END
      

*************************************************************************
*
*     * SUBROUTINE SSSR1 *
*
*      
*     * Purpose *
*      
*     Subspace minimization and computation of the search direction by
*     the limited memory SR1 update.
*
*      
*     * Calling sequence *
*     
*     CALL SSSR1(N,NACT,MN,IOLD,ITYPE,X,XL,XU,XCP,GA,D,AMUGAD,XBX,DSTRN,
*    &     ALPHA,SM,UM,RM,UMTUM,SMTGA,UMTGA,MWTGA,TMPMAT,TMPMA2,VMC3,VMC4,
*    &     VN1,VN2,VN3,GAMMA,SMALL,IACT,IB,IPRINT,IERR)
*     
*      
*     * Parameters *
*
*     II  N               Number of variables.
*     II  NACT            Number of active variables.
*     II  MN              Current number of stored corrections.
*     II  IOLD            Index of the oldest corrections.
*     II  ITYPE           Type of problem:
*                           0  - problem is unbounded,
*                           1  - constrained problem.
*     II  IACT(N)         Index set of active and free variables:
*                           for I=1,...,NACT, IACT(I) are the indices
*                               of active variables (-J, if
*                               XCP(J)=XL(J) and +J, otherwise) and
*                           for I=NACT+1,...,N, IACT(I) are the indices
*                               of free variables.
*     II  IB(N)           Type of bound constraints:
*                           0  - X(I) is unbounded,
*                           1  - X(I) has only a lower bound,
*                           2  - X(I) has both lower and upper bounds,
*                           3  - X(I) has only an upper bound. 
*     RI  X(N)            Vector of variables. 
*     RI  XCP(N)          Generalized Chauchy point.
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables. 
*     RI  GA(N)           Current aggregate subgradient of the objective
*                           function.
*     RU  D(N)            Search direction.
*     RO  AMUGAD          AMUGAD = -(A*MU + G)'*D, where A*MU denotes
*                           the Lagrange multipliers for problem.
*     RO  XBX             XBX = (XCP-X)'*BM*(XCP-X), where BM denotes the 
*                           limited memory approximation of the Hessian 
*                           matrix.
*     RO  ALPHA           Backtrack multiplier.
*     RO  DSTRN           Euclidean norm of the direction vector without 
*                           backtracking.
*     RI  SMTGA(MN)       Vector SMTGA = SM'*GA.
*     RI  UMTGA(MN)       Vector UMTGA = UM'*GA.
*     RA  MWTGA(MN)       Auxiliary vector.
*     RI  SM(N*MN)        Matrix whose columns are stored corrections.
*     RI  UM(N*MN)        Matrix whose columns are stored subgradient
*                           differences.
*     RI  RM(MN*(MN+1)/2) Upper triangular matrix stored
*                           columnwise in the one-dimensional array.
*     RI  UMTUM(MN*(MN+1)/2)  Matrix UMTUM = UM'*UM.
*     RI  TMPMAT((MN+1)*(MN)/2)  Factorization of matrix 
*                           UMTUM-RM-RM'+CDIAG.
*     RA  TMPMA2(MN*(MN+1)/2) Auxiliary matrix.
*     RA  VMC#(MN)        Auxiliary arrays; # = 3,4.
*     RA  VN#(N)          Auxiliary array; # = 1,2,3.
*     RI  GAMMA           Scaling parameter.
*     RI  SMALL           Small positive value.
*     II  IPRINT          Printout specification.
*     IO  IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -2  - Error in LINEQ.
*
*
*     * Subprograms used *
*      
*     S   CALQ            Solving X from linear equation A*X=Y.
*     S   LINEQ           Solving X from linear equation A*X=Y, where
*                           A is given in factorized form.
*     S   RECMAX          Multiplication of a vector by a dense
*                           rectangular matrix.
*     S   SCALEX          Scaling a vector.
*     S   XSUMY           Sum of two vectors.
*     S   SCDIFF          Difference of the scaled vector and a vector.
*     S   SCSUM           Sum of a vector and the scaled vector.
*     S   VNEG            Copying of a vector with change of the sign
*     RF  VDOT            Dot product of two vectors.
*
*
*     All the MN-vectors are stored in a circular order controlled by
*     the parameter point IOLD.
*
*      

      SUBROUTINE SSSR1(N,NACT,MN,IOLD,ITYPE,X,XL,XU,XCP,GA,D,AMUGAD,XBX,
     &     DSTRN,ALPHA,SM,UM,RM,UMTUM,SMTGA,UMTGA,MWTGA,TMPMAT,TMPMA2,
     &     VMC3,VMC4,VN1,VN2,VN3,GAMMA,SMALL,IACT,IB,IPRINT,IERR)

      
*     Scalar Arguments
      INTEGER N,NACT,MN,IOLD,ITYPE,IPRINT,IERR
      DOUBLE PRECISION AMUGAD,XBX,DSTRN,ALPHA,GAMMA,SMALL
      
*     Array Arguments
      INTEGER IACT(*),IB(*)
      DOUBLE PRECISION X(*),XL(*),XU(*),XCP(*),GA(*),D(*),SM(*),UM(*),
     &     SMTGA(*),UMTGA(*),MWTGA(*),RM(*),UMTUM(*),TMPMAT(*),
     &     TMPMA2(*),VMC3(*),VMC4(*),VN1(*),VN2(*),VN3(*)
      
*     Local Scalars
      INTEGER I,IJ,J,K,LI,LJ,IFREE,IND
      DOUBLE PRECISION DFREE,SIGNUM

*     Intrinsic Functions
      INTRINSIC ABS,SIGN,MIN,DBLE

*     External Functions
      DOUBLE PRECISION VDOT
      EXTERNAL VDOT

*     External Subroutines
      EXTERNAL SCDIFF,CALQ,SCALEX,RECMAX,XSUMY,LINEQ,VNEG

      
      IF (NACT .EQ. 0) THEN


*     
*     Solving MWTGA from TMPMAT*MWTGA= GAMMA*UM'*GA - SM'*GA
*

         CALL SCDIFF(MN,GAMMA,UMTGA,SMTGA,MWTGA)
         CALL LINEQ(MN,MN,IOLD,TMPMAT,MWTGA,MWTGA,SMALL,IERR)
         IF (IERR .NE. 0) RETURN

      
*
*     Computation of the search direction D
*

         CALL SCALEX(MN,GAMMA,MWTGA,VMC3)
         CALL RECMAX(N,MN,SM,MWTGA,VN1)
         CALL SCDIFF(N,-GAMMA,GA,VN1,D)
         CALL RECMAX(N,MN,UM,VMC3,VN1)
         CALL XSUMY(N,D,VN1,D)
         
         AMUGAD = -VDOT(N,GA,D)
         
      ELSE


*     
*     Solving MWTGA from TMPMAT*MWTGA = GAMMA*UM'*GA - SM'*GA
*

         CALL SCDIFF(MN,GAMMA,UMTGA,SMTGA,MWTGA)
         CALL LINEQ(MN,MN,IOLD,TMPMAT,MWTGA,MWTGA,SMALL,IERR)
         IF (IERR .NE. 0) RETURN


*
*     Computation of an intermediate vector VN1 = -A'*DM*GA - A'*(XCP-X)
*         

         DO 210 I=1,NACT
            VN1(I)=0.0D+00
            IND=ABS(IACT(I))
            SIGNUM = SIGN(1.0D+00,DBLE(IACT(I)))
            VN2(I) = SIGNUM*(XCP(IND)-X(IND))
            
            DO 211 J=1,MN
               VN1(I) = VN1(I) + SIGNUM * MWTGA(J) * 
     &              (GAMMA*UM((J-1)*N+IND) - SM((J-1)*N+IND))
 211        CONTINUE
            
            VN1(I) = -SIGNUM*GAMMA*GA(IND) + VN1(I) - VN2(I)
 210     CONTINUE
         

*
*     Computation of Lagrange multipiliers VN1 = (A'DM A)^{-1}*VN1.
*
         
*     VMC3 = (GAMMA*UM - SM)'*A*VN1 

         DO 220 J = 1, MN
            VMC3(J)= 0.0D+00
            DO 221 I=1,NACT
               IND=ABS(IACT(I))
               SIGNUM = SIGN(1.0D+00,DBLE(IACT(I)))
               VMC3(J)=VMC3(J)+SIGNUM*(GAMMA*UM((J-1)*N+IND)-
     &              SM((J-1)*N+IND))*VN1(I)
 221        CONTINUE
 220     CONTINUE


*     Formulation of  TMPMA2 = TMPMAT - 1/GAMMA * W'*A*A'*W, where
*     W = GAMMA*UM - SM and TMPMAT is the middle matrix in compact SR1
*     formula.

         IJ=0
         DO 250 I = 1,MN
            LI=I+IOLD-1
            IF (LI .GT. MN) LI=LI-MN
            
            DO 251 J=1,I
               LJ=J+IOLD-1
               IF (LJ .GT. MN) LJ=LJ-MN
               IJ=IJ+1
               
               TMPMA2(IJ)=0.0D+00
               
               DO 252 K = 1,NACT
                  IND = ABS(IACT(K))
                  TMPMA2(IJ) = TMPMA2(IJ) +
     &                 (GAMMA*UM((LJ-1)*N+IND) - SM((LJ-1)*N+IND)) *
     &                 (GAMMA*UM((LI-1)*N+IND) - SM((LI-1)*N+IND))
 252           CONTINUE
               
               TMPMA2(IJ)= - TMPMA2(IJ)/GAMMA + GAMMA*UMTUM(IJ) - RM(IJ)

 251        CONTINUE
            
 250     CONTINUE

     
*     Solving VMC3 from TMPMA2*VMC3 = (GAMMA*UM - SM)' A VN1 

         CALL CALQ(MN,MN,IOLD,TMPMA2,VMC3,VMC3,SMALL,IPRINT,IERR)
         IF (IERR .NE. 0) RETURN
         
         CALL VNEG(N,GA,VN3)
         
         DO 260 I=1,NACT
            VN2(I)=0.0D+00
            IND=ABS(IACT(I))
            SIGNUM = SIGN(1.0D+00,DBLE(IACT(I)))
            
            DO 261 J=1,MN
               VN2(I) = VN2(I) + SIGNUM * VMC3(J) *
     &              (GAMMA*UM((J-1)*N+IND) - SM((J-1)*N+IND))
 261        CONTINUE
            
            VN1(I) = (VN1(I) + VN2(I)/GAMMA) /GAMMA
            VN3(IND) = VN3(IND) - SIGNUM*VN1(I)
 260     CONTINUE
            
     
*     Solving VMC3 from TMPMAT*VMC3= (GAMMA*UM-SM)'*A VN1
*     Note that TMPMAT is a factorization of the middle matrix in SR1
*     updating formula.

         DO 400 J = 1, MN
            VMC3(J)= 0.0D+00
            DO 401 I=1,NACT
               IND=ABS(IACT(I))
               SIGNUM = SIGN(1.0D+00,DBLE(IACT(I)))
               VMC3(J)=VMC3(J)+SIGNUM*(GAMMA*UM((J-1)*N+IND)-
     &              SM((J-1)*N+IND))*VN1(I)
 401        CONTINUE

 400     CONTINUE
         
         CALL LINEQ(MN,MN,IOLD,TMPMAT,VMC3,VMC3,SMALL,IERR)
         IF (IERR .NE. 0) RETURN
         
         
         DO 410 J = 1, MN
            VMC3(J) = VMC3(J) + MWTGA(J)
 410     CONTINUE
         
         
*     
*     Computation of the search direction D
*
            
         CALL SCALEX(MN,GAMMA,VMC3,VMC4)
         CALL RECMAX(N,MN,SM,VMC3,D)
         CALL RECMAX(N,MN,UM,VMC4,VN2)
         
         DO 460 I=1,N
            D(I) = -GAMMA*GA(I) + VN2(I) - D(I)
 460     CONTINUE
               
         DO 470 I=1,NACT
            IND=ABS(IACT(I))
            SIGNUM = SIGN(1.0D+00,DBLE(IACT(I)))
            
            D(IND) = D(IND) - GAMMA*SIGNUM*VN1(I)
 470     CONTINUE

         AMUGAD = VDOT(N,D,VN3)

      END IF


 600  CONTINUE


*     
*     Euclidean norm of the direction vector without backtracking.
*

      IF (ITYPE .EQ. 1) THEN
         J=0
         DO 630 I=1,NACT
            IND=ABS(IACT(I))

            IF(abs(D(IND)-(XCP(IND)-X(IND))) .GT. 1.0D-6) THEN
               J=J+1
            END IF

            D(IND)=XCP(IND)-X(IND)
 630     CONTINUE
      END IF

      DSTRN = VDOT(N,D,D)


*     
*     Backtrack to the feasible region if X + D violates bounds
*     (among free variables).
*
      
      IF (ITYPE .EQ. 1) THEN

         ALPHA = 1.0D+00
         DO 610 I = NACT+1,N
            IFREE=IACT(I)
            DFREE=D(IFREE)
            IF (IB(IFREE) .NE. 0) THEN
               IF (IB(IFREE) .LE. 2 .AND. 
     &              X(IFREE) + DFREE .LT. XL(IFREE)) THEN

*     Otherwise the lower bound is never met
 
                  IF (ABS(X(IFREE)-XCP(IFREE)+DFREE) .GT. SMALL) THEN
                     ALPHA=MIN(ALPHA,(XL(IFREE)-XCP(IFREE))/
     &                    (X(IFREE)-XCP(IFREE)+DFREE))
                  END IF

               ELSE IF (IB(IFREE) .GE. 2 .AND. 
     &                 X(IFREE) + DFREE  .GT. XU(IFREE)) THEN

*     Otherwise the upper bound is never met
 
                  IF (ABS(X(IFREE)-XCP(IFREE)+DFREE) .GT. SMALL) THEN
                     ALPHA=MIN(ALPHA,(XU(IFREE)-XCP(IFREE))/
     &                    (X(IFREE)-XCP(IFREE)+DFREE))
                  END IF
                  
               END IF
               
            END IF
 610     CONTINUE


         IF (ALPHA .LT. 1.0D+00) THEN
            
            DO 620 I=NACT+1,N
               IFREE=IACT(I)
               D(IFREE) = XCP(IFREE) - X(IFREE) +
     &              ALPHA*(X(IFREE) - XCP(IFREE) + D(IFREE))
 620        CONTINUE

         ELSE
      
            XBX = 0.0D+00
         END IF

      END IF
      
      RETURN
      END
      
      
************************************************************************
*
*     * SUBROUTINE DLSR1 *
*
*      
*     * Purpose *
*      
*     Computation of the search direction by the limited memory SR1 
*     update.
*
*      
*     * Calling sequence *
*     
*     CALL DLSR1(N,MN,IOLD,X,XCP,XL,XU,GA,D,ITYPE,NACT,IB,
*    &     IACT,IFLAG,SM,UM,LM,RM,CDIAG,SMTSM,UMTUM,SMTGA,UMTGA,
*    &     VN1,VN2,VN3,VMC3,VMC4,VMC5,VMC6,TMPMAT,TMPMA2,
*    &     GAMMA,AMUGAD,XBX,DSTRN,ALPHA,SMALL,IPRINT,IERR)
*     
*     
*     * Parameters *
*
*     II  N               Number of variables.
*     II  MN              Current number of depositories used.
*     II  IOLD            Index for circular arrays.
*     RI  X(N)            Vector of variables.
*     RO  XCP(N)          Generalized Cauchy point.
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables. 
*     RI  GA(N)           Current aggregate subgradient.
*     RO  D(N)            Search direction.
*     II  ITYPE           Type of problem:
*                           0  - problem is unbounded,
*                           1  - constrained problem.
*     IO  NACT            Number of active variables.
*     II  IB(N)           Type of bound constraints:
*                           0  - X(I) is unbounded,
*                           1  - X(I) has only a lower bound,
*                           2  - X(I) has both lower and upper bounds,
*                           3  - X(I) has only an upper bound. 
*     IO  IACT(N)         Index set of active and free variables:
*                           for I=1,...,NACT, IACT(I) are the indices
*                               of active variables (-J, if
*                               XCP(J)=XL(J) and +J, otherwise)
*                           for i=NACT+1,...,N, IACT(I) are the indices
*                               of free variables.
*     II  IFLAG           Index for adaptive version:
*                           0  - Maximum number of stored corrections
*                                  has not been changed at previous
*                                  iteration.
*                           1  - Maximum number of stored corrections
*                                  has been changed at previous
*                                  iteration.
*     RI  SM(N*MN)        Matrix whose columns are stored corrections.
*     RI  UM(N*MN)        Matrix whose columns are stored subgradient
*                           differences.
*     RI  LM(MN*(MN+1)/2) Lower triangular matrix stored rowwise in the
*                           one-dimensional array.
*     RI  RM(MN*(MN+1)/2) Upper triangular matrix stored columnwise
*                           in the one-dimensional array.
*     RI  CDIAG(MN)       Diagonal matrix.
*     RI  SMTSM(MN*(MN+1)/2)  Matrix SMTSM = SM'*SM.
*     RI  UMTUM(MN*(MN+1)/2)  Matrix UMTUM = UM'*UM.
*     RI  SMTGA(MN)       Vector SMTGA = SM'*GA.
*     RI  UMTGA(MN)       Vector UMTGA = UM'*GA.
*     RI  TMPMAT((MN+1)*(MN)/2)  Factorization of matrix 
*                           UMTUM-RM-RM'+CDIAG.
*     RA  TMPMA2(MN*(MN+1)/2)  Auxiliary matrix.
*     RA  VN#(N)          Auxiliary arrays: # = 1,3.
*     RA  VMC#(N)         Auxiliary arrays: # = 3,6.
*     RI  GAMMA           Scaling parameter.
*     RU  AMUGAD          AMUGAD = -(A*MU + G)'*D, where A*MU denotes
*                           the Lagrange multipliers for problem.
*     RU  XBX             XBX = (XCP-X)'*BM*(XCP-X), where BM denotes the 
*                           limited memory approximation of the Hessian 
*                           matrix.
*     RO  ALPHA           Backtrack multiplier.
*     RO  DSTRN           Euclidean norm of the direction vector without 
*                           backtracking.
*     RI  SMALL           Small positive value.
*     II  IPRINT          Printout specification.
*     IO  IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -2  - Error in LINEQ.
*                            -3  - Error in TRLIEQ.
*     
*     
*      
*     * Subprograms used *
*      
*     S   ACTVAR          Finding the index set of free and active
*                           variables  at the generalized Cauchy point.
*     S   CPSR1           Computing the generalized Cauchy point by
*                           limited memory SR1 update formula.
*     S   SSSR1           Subspace minimization and computation of the
*                           search direction by the L-SR1 update.
*     S   XDIFFY          Difference of two vectors.
*     
*
*     The variable and subgradient differences and all the MN-vectors are
*     stored in a circular order controlled by the parameter point IOLD.
*
*     
  
      SUBROUTINE DLSR1(N,MN,IOLD,X,XCP,XL,XU,GA,D,ITYPE,NACT,IB,
     &     IACT,IFLAG,SM,UM,LM,RM,CDIAG,SMTSM,UMTUM,SMTGA,UMTGA,
     &     VN1,VN2,VN3,VMC3,VMC4,VMC5,VMC6,TMPMAT,TMPMA2,
     &     GAMMA,AMUGAD,XBX,DSTRN,ALPHA,SMALL,IPRINT,IERR)


*     Scalar Arguments
      INTEGER N,MN,IOLD,NACT,ITYPE,IFLAG,IPRINT,IERR
      DOUBLE PRECISION AMUGAD,XBX,DSTRN,ALPHA,GAMMA,SMALL
      
*     Array Arguments
      INTEGER IB(*),IACT(*)
      DOUBLE PRECISION X(*),XCP(*),XL(*),XU(*),GA(*),D(*),SM(*),UM(*),
     &     LM(*),RM(*),CDIAG(*),SMTSM(*),UMTUM(*),SMTGA(*),UMTGA(*),
     &     VN1(*),VN2(*),VN3(*),VMC3(*),VMC4(*),VMC5(*),
     &     VMC6(*),TMPMAT(*),TMPMA2(*)

*     External Subroutines
      EXTERNAL ACTVAR,CPSR1,SSSR1,XDIFFY


      IF (ITYPE .EQ. 1) THEN


*     
*     Generalized Cauchy point.
*

         CALL CPSR1(N,MN,IOLD,IFLAG,IPRINT,X,XCP,XL,XU,IB,IACT,GA,XBX,
     &        VN2,VN1,SM,UM,SMTGA,UMTGA,SMTSM,LM,CDIAG,TMPMA2,VMC3,VMC4,
     &        VMC5,VMC6,GAMMA,SMALL,IERR)

         IF (IERR .NE. 0) RETURN


*
*     Determination of the active variables.
*

         CALL ACTVAR(N,NACT,IACT,IB,XCP,XL,XU,SMALL)

      ELSE

         XBX = 0.0D+00
      END IF


      IF (NACT .NE. N) THEN


*
*     Subspace minimization and direction finding
*      

         CALL SSSR1(N,NACT,MN,IOLD,ITYPE,
     &        X,XL,XU,XCP,GA,D,AMUGAD,XBX,DSTRN,ALPHA,
     &        SM,UM,RM,UMTUM,SMTGA,UMTGA,VMC5,TMPMAT,TMPMA2,VMC3,
     &        VMC4,VN1,VN2,VN3,GAMMA,SMALL,IACT,IB,IPRINT,IERR)
         IF (IERR .NE. 0) RETURN

      ELSE
         CALL XDIFFY(N,XCP,X,D)
         AMUGAD = XBX
         
      END IF

      RETURN
      END
      

************************************************************************
*
*
*     * SUBROUTINE AGRSR1 *
*
*      
*     * Purpose *
*      
*     Computation of aggregate values by the limited memory SR1 update.
*
*      
*     * Calling sequence *
*     
*     CALL AGRSR1(N,MN,IOLD,G,GM,GA,PG,PGM,PGA,PGAH,PGAHG,PGAHGA,
*    &     PGNRM,PGMNRM,BETA,AGBETA,SM,UM,TMPMAT,UMTUM,RM,SMTPGM,
*    &     UMTPGM,VN1,VN2,VMC1,VMC2,GAMMA,SMALL,IERR)
*
*      
*     * Parameters *
*
*     II  N               Number of variables.
*     II  MN              Current number of stored corrections.
*     II  IOLD            Index for circular arrays.
*     RI  G(N)            Current (auxiliary) subgradient of the
*                           objective function.
*     RI  GM(N)           Basic subgradient of the objective function.
*     RU  GA(N)           Current aggregate subgradient.
*     RI  PG(N)           Current projected (auxiliary) subgradient of 
*                           the objective function.
*     RI  PGM(N)          Projected basic subgradient of the objective 
*                           function.
*     RI  PGA(N)          Projected aggregate subgradient.
*     RI  PGAH(N)         PGAH = - H*PGA, where H presents the L-SR1- 
*                           approximation of Hessian.
*     RI  PGAHG           PGAHG = - PGA'*H*PG.
*     RI  PGAHGA          PGAHGA = PGA'*H*PGA.
*     RI  PGNRM           PGNRM = PG'*PG.
*     RI  PGMNRM          PGNRM = PGM'*PGM.
*     RI  BETA            Locality measure.
*     RU  AGBETA          Aggregate locality measure.
*     RI  SM(N*MN)        Matrix whose columns are stored corrections.
*     RI  UM(N*MN)        Matrix whose columns are stored subgradient.
*     RI  TMPMAT((MN+1)*(MN)/2)  Factorization of matrix 
*                           UMTUM-RM-RM'+CDIAG.
*     RI  UMTUM(MN*(MN+1)/2)  Matrix UMTUM = UM'*UM.
*     RI  RM(MN*(MN+1)/2) Upper triangular matrix.
*     RI  SMTPGM(MN)      Vector SMTPGM = SM'*PGM.
*     RI  UMTPGM(MN)      Vector UMTPGM = UM'*PGM.
*     RA  VN#(N)          Auxiliary arrays; # = 1,2.
*     RA  VMC#(MN)        Auxiliary arrays; # = 1,2.
*     RI  GAMMA           Scaling parameter.
*     RI  SMALL           Small positive value.
*     IO  IERR            Error indicador: 
*                             0  - Everything is ok.
*                            -2  - Error In LINEQ.
*
*     
*     * Local variables *
*
*     R   PR              PR = (PGM-PGA)'*H*(PGM-PGA), where H presents 
*                           the L-SR1- approximation of Hessian.
*     R   PRQR            PRQR = (PGM-PGA)'*H*(PG-PGA).
*     R   QR              QR = (PG-PGA)'*H*(PG-PGA).
*     R   RRP             RRP = (PGM-PGA)'*H*PGA - AGBETA.
*     R   RRQ             RRQ = (PG-PGA)'*H*PGA - AGBETA + BETA.
*     R   PQ              PQ = (PG-PGM)'*H*(PG-PGM).
*     R   QQP             QQP = (PG-PGM)'*H*PG + BETA.
*     R   LAM#            Multipliers: # = 1,2.
*     R   TMP#            Auxiliary scalars: # = 1,2,3.
*
*     
*     * Subprograms used *
*      
*     S   CALQ            Solving X from linear equation A*X=Y.
*     S   LINEQ           Solving X from linear equation A*X=Y, where
*                           A is in factorized form.
*     S   RECMAX          Multiplication of a columnwise stored dense 
*                           rectangular matrix by a vector.
*     S   RWAXV2          Multiplication of two rowwise stored dense 
*                           rectangular matrices by two vectors.
*     S   SCALEX          Scaling a vector.
*     S   SCDIFF          Difference of the scaled vector and a vector.
*     S   SCSUM           Sum of a vector and the scaled vector.
*     S   XDIFFY          Difference of two vectors.
*     RF  VDOT            Dot product of two vectors.
*     
*
*     The variable and subgradient differences and all the MN-vectors are
*     stored in a circular order controlled by the parameter point IOLD.
*

     
      SUBROUTINE AGRSR1(N,MN,IOLD,G,GM,GA,PG,PGM,PGA,PGAH,PGAHG,PGAHGA,
     &     PGNRM,PGMNRM,BETA,AGBETA,SM,UM,TMPMAT,UMTUM,RM,SMTPGM,UMTPGM,
     &     VN1,VN2,VMC1,VMC2,GAMMA,SMALL,IERR)

*     Scalar Arguments
      INTEGER N,MN,IOLD,IERR
      DOUBLE PRECISION PGAHG,PGAHGA,PGNRM,PGMNRM,BETA,AGBETA,GAMMA,
     &     SMALL
      
*     Array Arguments
      DOUBLE PRECISION G(*),GM(*),GA(*),PG(*),PGM(*),PGA(*),PGAH(*),
     &     TMPMAT(*),SM(*),UM(*),UMTUM(*),RM(*),SMTPGM(*),UMTPGM(*),
     &     VN1(*),VN2(*),VMC1(*),VMC2(*)
      
*     Local Scalars
      INTEGER I
      DOUBLE PRECISION PR,PRQR,PQ,QQP,QR,RRP,RRQ,LAM1,LAM2,TMP1,TMP2,
     &     TMP3,SMALL2
     
*     External Functions
      DOUBLE PRECISION VDOT
      EXTERNAL VDOT

*     External Subroutines
      EXTERNAL CALQ,LINEQ,RECMAX,RWAXV2,SCALEX,SCDIFF,SCSUM,XDIFFY

*     Intrinsic Functions
      INTRINSIC SQRT
      INTRINSIC ABS,MIN,MAX


      SMALL2 = SQRT(SMALL)
      IERR = 0
      

      RRQ = -PGAHG - PGAHGA
      
      CALL XDIFFY(N,PGM,PGA,VN1)
      RRP = - VDOT(N,VN1,PGAH)
      
      
*      
*     Calculation of VN1 = trans(GM)*H
*

      IF (MN .GT. 0) THEN

         CALL SCDIFF(MN,GAMMA,UMTPGM,SMTPGM,VMC1)
         CALL LINEQ(MN,MN,IOLD,TMPMAT,VMC2,VMC1,SMALL,IERR)
         IF (IERR .NE. 0) RETURN

         CALL SCALEX(MN,GAMMA,VMC2,VMC1)
         CALL RECMAX(N,MN,SM,VMC2,VN2)
         CALL SCSUM(N,GAMMA,PGM,VN2,VN1)
         CALL RECMAX(N,MN,UM,VMC1,VN2)
         CALL XDIFFY(N,VN1,VN2,VN1)

         PR = VDOT(N,VN1,PGM) - 2.0D+00*RRP - PGAHGA
         PRQR = VDOT(N,VN1,PG) - RRP - PGAHGA - RRQ

      ELSE

         PR = PGMNRM - 2.0D+00*RRP - PGAHGA
         PRQR = VDOT(N,PGM,PG) - RRP - PGAHGA - RRQ

      END IF

*
*     Calculation of QR = trans(G-GA)*H*(G-GA)
*
      QR = PGNRM

      IF (MN .GT. 0) THEN
         QR = GAMMA*QR

         CALL RWAXV2(N,MN,SM,UM,PG,PG,VMC1,VMC2)
         CALL SCSUM(MN,-GAMMA,VMC2,VMC1,VMC1)
         CALL LINEQ(MN,MN,IOLD,TMPMAT,VMC2,VMC1,SMALL,IERR)
            
         QR = QR - VDOT(MN,VMC1,VMC2) - 2.0D+00*RRQ - PGAHGA

      ELSE
         QR = QR - 2.0D+00*RRQ - PGAHGA
         
      END IF

      PQ = QR - 2.0D+00*PRQR + PR
      IF (PQ .LT. 0.0D+00) THEN
         PQ = 0.0D+00
         PRQR = (QR + PR) / 2.0D+00
      END IF

      QQP = PQ + PRQR + RRQ - PR - RRP + BETA
      RRP = RRP - AGBETA
      RRQ = RRQ + BETA - AGBETA

      
*     
*     Computation of multipliers LAM1 and LAM2
*

      IF (PR .LE. SMALL2 .OR. QR .LE. SMALL2) GOTO 100
      TMP1 = RRQ/QR
      TMP2 = PRQR/QR
      TMP3 = PR - PRQR*TMP2

      IF (ABS(TMP3) .LE. SMALL2) GOTO 100
      LAM1 = (TMP1*PRQR - RRP)/TMP3
      LAM2 = -TMP1 - LAM1*TMP2

      IF (LAM1*(LAM1 - 1.0D+00) .LT. 0.0D+00 .AND.
     &     LAM2*(LAM1 + LAM2 - 1.0D+00) .LT. 0.0D+00) GOTO 200

      
*
*     Minimum on the boundary
*

 100  CONTINUE

      LAM1 = 0.0D+00
      LAM2 = 0.0D+00
      IF (BETA .LE. AGBETA) LAM2 = 1.0D+00
      IF (QR .GT. SMALL2) LAM2 = MIN(1.0D+00,MAX(0.0D+00,-RRQ/QR))
      TMP3 = (LAM2*QR + 2.0D+00*RRQ)*LAM2
      TMP1 = 0.0D+00
      IF (AGBETA .GE. 0.0D+00) TMP1 = 1.0D+00
      IF (PR .GT. SMALL2) TMP1 = MIN(1.0D+00,MAX(0.0D+00,-RRP/PR))
      TMP2 = (TMP1*PR + 2.0D+00*RRP)*TMP1
      IF (TMP2 .LT. TMP3) THEN
         TMP3 = TMP2
         LAM1 = TMP1
         LAM2 = 0.0D+00
      END IF
      
      IF (QQP*(QQP - PQ) .GE. 0.0D+00) GOTO 200
      IF(PQ .LE. SMALL2) GOTO 200
      IF (QR + 2.0D+00*RRQ - QQP*QQP/PQ .GE. TMP3) GOTO 200

      LAM1 = QQP/PQ
      LAM2 = 1.0D+00 - LAM1


 200  CONTINUE

      IF (LAM1 .EQ. 0.0D+00 .AND. LAM2*(LAM2 - 1.0D+00) .LT. 0.0D+00
     &     .AND. -RRP - LAM2*PRQR .GT. 0.0D+00 .AND. PR .GT. SMALL2)
     &     LAM1 = MIN(1.0D+00 - LAM2, (-RRP-LAM2*PRQR)/PR)

      
*
*     Computation of the aggregate values
*      

      TMP1 = 1.0D+00 - LAM1 - LAM2
      DO 30 I=1,N
         GA(I)=LAM1*GM(I)+LAM2*G(I)+TMP1*GA(I)
 30   CONTINUE
      
      AGBETA = LAM2*BETA + TMP1*AGBETA
      
      RETURN
      END

      
************************************************************************
*
*     * SUBROUTINE TINIT *
*
*      
*     * Purpose *
*      
*     Initial stepsize selection for LMBM-B.
*
*
*     * Calling sequence *
*     
*     CALL TINIT(N,NA,MAL,X,AF,AG,AY,IBUN,D,F,P,T,TMAX,TMIN,ETA,ETA9,
*    &     MOS,SMALL,ITERS)
*     
*     
*     * Parameters *
*      
*     II  N             Number of variables.
*     II  NA            Maximum size of the bundle.
*     II  MAL           Current size of the bundle.
*     RI  X(N)          Vector of variables.
*     RU  AF(4*NA)      Vector of values of bundle functions.
*     RI  AG(N*NA)      Matrix whose columns are bundle subgradients.
*     RI  AY(N*NA)      Matrix whose columns are bundle points.
*     II  IBUN          Index for the circular arrays in bundle.
*     RI  D(N)          Direction vector.
*     RI  F             Value of the objective function.
*     RI  P             Directional derivative.
*     RO  T             Value of the stepsize parameter.
*     RI  TMIN          Lower limit for stepsize parameter.
*     RI  TMAX          Upper limit for stepsize parameter.
*     RI  ETA           Distance measure parameter.
*     RI  ETA9          Maximum for real numbers.
*     RI  MOS           Locality measure parameter.
*     RI  SMALL         A small positive value.
*     II  ITERS         Null step indicator.
*                            0  - Null step.
*                            1  - Serious step.
*     
*
*     * Subprograms used *
*      
*     S   DESTEP        Stepsize determination for descent steps.
*     S   NULSTP        Stepsize determination for null steps.
*
*      
*     Napsu Karmitsa (2003, last modified 2009)
*
      
      
      SUBROUTINE TINIT(N,NA,MAL,X,AF,AG,AY,IBUN,D,F,P,T,TMAX,TMIN,ETA,
     &     ETA9,MOS,SMALL,ITERS)

*     Scalar Arguments
      INTEGER N,NA,MAL,IBUN,MOS,ITERS
      DOUBLE PRECISION P,ETA,ETA9,F,T,TMAX,TMIN,SMALL

*     Array Arguments
      DOUBLE PRECISION AF(*),AG(*),AY(*),D(*),X(*)

*     Local Scalars
      DOUBLE PRECISION SMALL2

*     External Subroutines
      EXTERNAL DESTEP,NULSTP

*     Intrinsic Functions
      INTRINSIC MAX,MIN,ABS,SQRT


      SMALL2 = SQRT(SMALL)
      T = MIN(1.0D+00,TMAX)

      IF (ABS(P) .LE. SMALL2) GO TO 10

      IF (ITERS.EQ.1) THEN
         CALL DESTEP(N,NA,MAL,X,AF,AG,AY,IBUN,D,F,P,T,ETA,ETA9,MOS)

      ELSE
         CALL NULSTP(N,NA,MAL,X,AF,AG,AY,IBUN,D,F,P,T,ETA,ETA9,MOS)
      END IF

 10   CONTINUE

      T = MIN(MAX(T,TMIN),TMAX)
      
      RETURN
      END


************************************************************************
*
*     * SUBROUTINE DESTEP *
*
*      
*     * Purpose *
*      
*     Stepsize selection using polyhedral approximation for descent step
*     in limited memory bundle method.
*
*
*     * Calling sequence *
*     
*     CALL DESTEP(N,NA,MAL,X,AF,AG,AY,IBUN,D,F,DF,T,ETA,ETA9,MOS)
*     
*     
*     * Parameters *
*      
*     II  N             Number of variables.
*     II  NA            Maximum size of the bundle.
*     II  MAL           Current size of the bundle.
*     RI  X(N)          Vector of variables.
*     RU  AF(4*NA)      Vector of values of bundle functions.
*     RI  AG(N*NA)      Matrix whose columns are bundle subgradients.
*     RI  AY(N*NA)      Matrix whose columns are bundle points.
*     II  IBUN          Index for the circular arrays in bundle.
*     RI  D(N)          Direction vector.
*     RI  F             Value of the objective function.
*     RI  DF            Directional derivative.
*     RO  T             Value of the stepsize parameter.
*     RI  ETA           Distance measure parameter.
*     RI  ETA9          Maximum for real numbers.
*     RI  MOS           Locality measure parameter.
*
*      
*     * Original code *
*
*     Part of PVAR-software by L.Luksan and J. Vlcek
*
*

      SUBROUTINE DESTEP(N,MA,MAL,X,AF,AG,AY,IBUN,D,F,DF,T,ETA,ETA9,MOS)

*     Scalar Arguments
      INTEGER N,MA,MAL,MOS,IBUN
      DOUBLE PRECISION DF,ETA,ETA9,F,T

*     Array Arguments
      DOUBLE PRECISION AF(*),AG(*),AY(*),D(*),X(*)

*     Local Scalars
      INTEGER I,J,JN,K,L,LQ,IB
      DOUBLE PRECISION ALF,ALFL,ALFR,BET,BETL,BETR,DX,Q,R,W

*     Intrinsic Functions
      INTRINSIC ABS,DBLE,MAX,MIN,SQRT


      ALFL = 0.0D+00
      BETL = 0.0D+00
      
      W = DF * T * (1.0D+00-T*0.5D+00)
      
      
*     
*     Initial choice of possibly active lines
*      

      K = 0
      L = -1
      JN = (IBUN-1)*N
      BETR = -ETA9

      DO 20 J = 1,MAL - 1
         IB = IBUN - 1 + J
         IF (IB .GT. MAL) IB = IB - MAL
         IF (JN .GE. MAL*N) JN = JN - MAL*N
         R = 0.0D+00
         BET = 0.0D+00
         ALFL = AF(IB) - F

         DO 10 I = 1,N
            DX = X(I) - AY(JN+I)
            Q = AG(JN+I)
            R = R + DX*DX
            ALFL = ALFL + DX*Q
            BET = BET + D(I)*Q
 10      CONTINUE

         IF (MOS.NE.2) R = R** (DBLE(MOS)*0.5D+00)
         ALF = MAX(ABS(ALFL),ETA*R)
         R = 1.0D+00 - BET/DF
         IF (R*R+ (ALF+ALF)/DF .GT. 1.0D-6) THEN
            K = K + 1
            AF(MA+K) = ALF
            AF(MA+MA+K) = BET
            R = T*BET - ALF
            IF (R.GT.W) THEN
               W = R
               L = K
            END IF

         END IF

         BETR = MAX(BETR,BET-ALF)
         JN = JN + N
         
 20   CONTINUE

      LQ = -1
      IF (BETR .LE. DF*0.5D+00) RETURN
      LQ = 1
      IF (L .LT. 0) RETURN
      BETR = AF(MA+MA+L)
      IF (BETR .LE. 0.0D+00) THEN
         IF (T .LT. 1.0D+00 .OR. BETR .EQ. 0.0D+00) RETURN
         LQ = 2
      END IF

      ALFR = AF(MA+L)

      
*
*     Iteration loop
*

 30   CONTINUE

      IF (LQ .GE. 1) THEN
         Q = 1.0D+00 - BETR/DF
         R = Q + SQRT(Q*Q + (ALFR+ALFR)/DF)
         IF (BETR .GE. 0.0D+00) R = - (ALFR+ALFR)/ (DF*R)
         R = MIN(1.95D+00,MAX(0.0D+00,R))
         
      ELSE
         IF (ABS(BETR-BETL)+ABS(ALFR-ALFL) .LT. -1.0D-4*DF) RETURN
         R = (ALFR-ALFL)/ (BETR-BETL)
      END IF

      IF (ABS(T-R).LT.1.0D-4) RETURN
      T = R
      AF(MA+L) = -1.0D+00
      W = T*BETR - ALFR
      L = -1

      DO 40 J = 1,K
         ALF = AF(MA+J)
         IF (ALF.LT.0.0D+00) GO TO 40
         BET = AF(MA+MA+J)
         R = T*BET - ALF
         IF (R.GT.W) THEN
            W = R
            L = J
         END IF
 40   CONTINUE

      IF (L.LT.0) RETURN
      BET = AF(MA+MA+L)
      IF (BET.EQ.0.0D+00) RETURN

      
*
*     New interval selection
*

      ALF = AF(MA+L)
      IF (BET.LT.0.0D0) THEN
         IF (LQ.EQ.2) THEN
            ALFR = ALF
            BETR = BET
            
         ELSE
            ALFL = ALF
            BETL = BET
            LQ = 0
         END IF

      ELSE
         IF (LQ.EQ.2) THEN
            ALFL = ALFR
            BETL = BETR
            LQ = 0
         END IF

         ALFR = ALF
         BETR = BET
      END IF

      GO TO 30

      END


************************************************************************
*
*     * SUBROUTINE NULSTP *
*
*
*     * Purpose *
*      
*     Stepsize selection using polyhedral approximation for null step
*     in limited memory bundle method.
*
*      
*     * Calling sequence *
*     
*     CALL NULSTP(N,NA,MAL,X,AF,AG,AY,IBUN,D,F,DF,T,ETA,ETA9,MOS)
*     
*     
*     * Parameters *
*      
*     II  N             Number of variables.
*     II  NA            Maximum size of the bundle.
*     II  MAL           Current size of the bundle.
*     RI  X(N)          Vector of variables.
*     RU  AF(4*NA)      Vector of values of bundle functions.
*     RI  AG(N*NA)      Matrix whose columns are bundle subgradients.
*     RI  AY(N*NA)      Matrix whose columns are bundle points.
*     II  IBUN          Index for the circular arrays in bundle.
*     RI  D(N)          Direction vector.
*     RI  F             Value of the objective function.
*     RI  DF            Directional derivative.
*     RO  T             Value of the stepsize parameter.
*     RI  ETA           Distance measure parameter.
*     RI  ETA9          Maximum for real numbers.
*     RI  MOS           Locality measure parameter.
*
*      
*     * Original code *
*
*     Part of PVAR-software by L.Luksan and J. Vlcek
*
*
*     * Limited memory version *
*
*     Marjo Haarala (2002,2003)
*      

      
      SUBROUTINE NULSTP(N,MA,MAL,X,AF,AG,AY,IBUN,D,F,DF,T,ETA,ETA9,MOS)

*     Scalar Arguments
      INTEGER MA,MAL,MOS,N,IBUN
      DOUBLE PRECISION DF,ETA,ETA9,F,T

*     Array Arguments
      DOUBLE PRECISION AF(*),AG(*),AY(*),D(*),X(*)

*     Local Scalars
      INTEGER I,J,JN,K,L,IB
      DOUBLE PRECISION ALF,ALFL,ALFR,BET,BETL,BETR,DX,Q,R,W

*     Intrinsic Functions
      INTRINSIC ABS,DBLE,MAX,MIN,SQRT


      W = DF*T


*     
*     Initial choice of possibly active parabolas
*

      K = 0
      L = -1
      JN = (IBUN-1)*N
      BETR = -ETA9

      DO 20 J = 1,MAL - 1
         IB = IBUN - 1 + J
         IF (IB .GT. MAL) IB = IB - MAL
         IF (JN .GE. MAL*N) JN = JN - MAL*N
         BET = 0.0D+00
         R = 0.0D+00
         ALFL = AF(IB) - F

         DO 10 I = 1,N
            DX = X(I) - AY(JN+I)
            R = R + DX*DX
            Q = AG(JN+I)
            ALFL = ALFL + DX*Q
            BET = BET + D(I)*Q
 10      CONTINUE

         IF (MOS.NE.2) R = R**(DBLE(MOS)*0.5D+00)
         ALF = MAX(ABS(ALFL),ETA*R)
         BETR = MAX(BETR,BET-ALF)

         IF (ALF.LT.BET-DF) THEN
            K = K + 1
            R = T*BET - ALF
            AF(MA+K) = ALF
            AF(MA+MA+K) = BET

            IF (R.GT.W) THEN
               W = R
               L = K
            END IF

         END IF

         JN = JN + N
 20   CONTINUE

      IF (L.LT.0) RETURN
      BETR = AF(MA+MA+L)
      ALFR = AF(MA+L)
      ALF = ALFR
      BET = BETR
      ALFL = 0.0D+00
      BETL = DF


*     
*     Iteration loop
*

 30   CONTINUE
      
      W = BET/DF
      IF (ABS(BETR-BETL)+ABS(ALFR-ALFL).LT.-1.0D-4*DF) RETURN
      IF (BETR-BETL .EQ. 0.0D+00) STOP 11
      R = (ALFR-ALFL)/ (BETR-BETL)
      IF (ABS(T-W).LT.ABS(T-R)) R = W
      Q = T
      T = R
      IF (ABS(T-Q).LT.1.0D-3) RETURN
      AF(MA+L) = -1.0D+00
      W = T*BET - ALF
      L = -1

      DO 40 J = 1,K
         ALF = AF(MA+J)
         IF (ALF.LT.0.0D+00) GO TO 40
         BET = AF(MA+MA+J)
         R = T*BET - ALF
         IF (R.GT.W) THEN
            W = R
            L = J
         END IF
 40   CONTINUE

      IF (L.LT.0) RETURN
      BET = AF(MA+MA+L)
      Q = BET - T*DF
      IF (Q.EQ.0.0D+00) RETURN

      
*     
*     New interval selection
*

      ALF = AF(MA+L)
      IF (Q.LT.0.0D+00) THEN
         ALFL = ALF
         BETL = BET

      ELSE
         ALFR = ALF
         BETR = BET
      END IF

      GO TO 30

      END


************************************************************************
*
*     * SUBROUTINE LLS3 *
*
*      
*     * Purpose *
*      
*     Special line search for LMBM-B.
*
*      
*     * Calling sequence *
*     
*     CALL LLS3(N,X,XO,XL,XU,IB,G,PG,PGAH,PGAHG,PGNRM,D,T,F,FO,
*    &     BETA,TMIN,DNORM,WK,THETA,TTHETA,EPSL,EPSR,ETA,MOS,ITERS,NFE,
*    &     NNK,SMALL,ITERM)
*     
*     
*     * Parameters *
*      
*     II  N             Number of variables.
*     RU  X(N)          Vector of variables.
*     RI  XO(N)         Previous vector of variables.
*     RI  XL(N)         Lower bounds for variables. 
*     RI  XU(N)         Upper bounds for variables. 
*     II  IB(N)         Type of bound constraints:
*                         0  - X(I) is unbounded,
*                         1  - X(I) has only a lower bound,
*                         2  - X(I) has both lower and upper bounds,
*                         3  - X(I) has only an upper bound. 
*     RO  G(N)          Subgradient of the objective function.
*     RO  PG(N)         Projected subgradient of the objective function.
*     RI  PGAH(N)       PGAH = -H * PGA, where H is an inverse 
*                         approximation ot the Hessia matrix and PGA is
*                         the projected aggregate subgradient.
*     RO  PGAHG         PGAHG = - PGA'*H*PG (only on null step).
*     RO  PGNRM         PGNRM = PG'*PG (only on null step).
*     RI  D(N)          Direction vector.
*     RU  T             Stepsize.
*     RO  F             Value of the objective function.
*     RI  FO            Previous value of the objective function.
*     RO  BETA          Locality measure.
*     RI  TMIN          Minimum value of the stepsize.
*     RI  DNORM         Euclidean norm of the direction vector.
*     RI  WK            Stopping parameter.
*     RI  THETA         Scaling parameter.
*     RO  TTHETA        TTHETA = T*THETA.
*     RI  EPSL          Termination tolerance for line search (in test
*                         on the change of the function value).
*     RI  EPSR          Termination tolerance for line search (in test
*                         on the directional derivative).
*     RI  ETA           Distance measure parameter.
*     II  MOS           Locality measure parameter.
*     IO  ITERS         Null step indicator.
*                         0  - Null step.
*                         1  - Serious step.
*     IU  NFE           Number of function evaluations.
*     II  NNK           Number of consequtive null steps.
*     RI  SMALL         A small positive value.
*     IO  ITERM         Cause of termination:
*                         0  - Everything is ok.
*                        -3  - Failure in function or subgradient
*                              calculations.
*     
*      
*     * Local parameters *
*
*     I   MAXIN1        Maximum number of additional interpolations.
*     I   MAXIN2        Maximum number of interpolations.
*
*      
*     * Local variables *
*
*     I   NIN           Number of interpolations.
*     I   ISGNGA        ISGNGA = 0, if G(I)<=0 for all I such that 
*                         X(I)=XU(I) and G(I)>=0 for all I such that
*                         X(I)=XL(I). Otherwise ISGNGA = 1. Part of the 
*                         subroutine call. Not really needed here.
*     R   GAD           Directional derivative.
*     R   TL,TU         Lower and upper limits for T used in
*                         interpolation.
*     R   FL            Value of the objective function for T=TL.
*     R   FU            Value of the objective function for T=TU.
*     R   EPSA          Line search parameter.
*     R   EPST          Line search parameter.
*     R   EPSLK         Line search parameter.
*     R   EPSRK         Line search parameter.
*     R   KAPPA         Interpolation parameter.
*      
*
*     * Subprograms used *
*      
*     S   PROJGR        Simple projection of the subgradient and 
*                         calculation of Eucleidean norm.
*     S   QINT          Quadratic interpolation for line search
*                         with one directional derivative.
*     S   SCSUM         Sum of a vector and the scaled vector.
*
*     RF  VDOT          Dot product of two vectors.
*
*      
*     * External subroutines *
*      
*     SE  FUNDER        Computation of the value and the subgradient of
*                       the objective function. Calling sequence:
*                       CALL FUNDER(N,X,F,G,ITERM), where N is a number
*                       of variables, X(N) is a vector of variables, F
*                       is the value of the objective function, G(N)
*                       is the subgradient of the objective function,
*                       and ITERM is the error indicator.
*      
*

      SUBROUTINE LLS3(N,X,XO,XL,XU,IB,G,PG,PGAH,PGAHG,PGNRM,D,T,F,FO,
     &     BETA,TMIN,DNORM,WK,THETA,TTHETA,EPSL,EPSR,ETA,MOS,ITERS,NFE,
     &     NNK,SMALL,ITERM)

*     Scalar Arguments
      INTEGER N,MOS,ITERS,NFE,NNK,ITERM
      DOUBLE PRECISION PGAHG,PGNRM,T,F,FO,BETA,TMIN,DNORM,WK,THETA,
     &     TTHETA,EPSL,EPSR,ETA,SMALL

*     Array Arguments
      INTEGER IB(*)
      DOUBLE PRECISION X(*),XO(*),XL(*),XU(*),G(*),PG(*),PGAH(*),D(*)

*     Local Scalars
      INTEGER NIN,ISGNGA
      DOUBLE PRECISION GAD,FL,FU,TL,TU,EPSA,EPST,EPSLK,EPSRK,KAPPA

*     Parameters
      INTEGER MAXIN1,MAXIN2
      PARAMETER(
     &     MAXIN1 = 20,
     &     MAXIN2 = 20)


*     External Functions
      DOUBLE PRECISION VDOT
      EXTERNAL VDOT

*     Intrinsic Functions
      INTRINSIC ABS,MAX

*     External Subroutines
      EXTERNAL FUNDER,PROJGR,SCSUM,QINT

      
*
*     Initialization
*      

      NIN = 0

      EPST = 2.0D+00 * EPSL
      EPSA = (EPSR - EPST)/2.0D+00
      
      TL = 0.0D+00
      TU = T
      FL = FO

      IF (THETA .LT. 1.0D+00) THEN
         EPST=THETA*EPST
         EPSA=THETA*EPSA
         EPSLK=EPSL
         EPSL=THETA*EPSL
         EPSRK=EPSR
         EPSR=THETA*EPSR
      END IF

      KAPPA = 1.0D+00 - 0.5D+00/(1.0D+00-EPST)


*     
*     Function and gradient evaluation at a new point
*

 10   CONTINUE

      TTHETA = THETA * T
      CALL SCSUM(N,TTHETA,D,XO,X)
      CALL FUNDER(N,X,F,G,ITERM)
      NFE = NFE + 1

      IF (ITERM .NE. 0) RETURN
      

      GAD = THETA*VDOT(N,G,D)
      BETA = MAX(ABS(FO-F+GAD*T),ETA*(TTHETA*DNORM)**MOS)


*     
*     Null/descent step test (ITERS=0/1)
*

      ITERS = 1
      IF (F .LE. FO - T*EPST*WK) THEN
         TL = T
         FL = F

      ELSE
         TU = T
         FU = F
      END IF

      
*
*     Serious step
*      

      IF (F .LE. FO - T*EPSL*WK .AND. (T .GE. TMIN .OR.
     &     BETA .GT. EPSA*WK)) GO TO 40

      
*
*     Additional interpolation
*      

      IF (F .GT. FO .AND. TU-TL .GE. TMIN
     &     .AND. NNK .GE. 1 .AND. NIN .LT. MAXIN1) THEN
         GO TO 20
      ENDIF

      
*
*     Null step
*      

*     Projection (with respect to XO)
      CALL PROJGR(N,XO,XL,XU,G,PG,PGNRM,ISGNGA,IB,SMALL)
 
      PGAHG = VDOT(N,PGAH,PG)

      IF (PGAHG - BETA .GE. -EPSR*WK) GO TO 30


*
*     Emergency exit (null step)
* 

      IF (TU-TL .LT. TMIN*1.0D-01 .OR. NIN .GE. MAXIN2) THEN
         ITERM = -1
         GO TO 30
      END IF


*
*     Interpolation
*      

 20   CONTINUE
      
      NIN=NIN+1

      IF (TL .EQ. 0.0D+00 .AND. WK .GT. 0.0D+00) THEN
         CALL QINT(TU,FL,FU,WK,T,KAPPA,SMALL)

      ELSE
         T = 0.50D+00 * (TU+TL)

      END IF

      GO TO 10


 30   CONTINUE
      ITERS = 0
      

 40   CONTINUE
      IF (THETA .LT. 1.0D+00) THEN
         EPSL=EPSLK
         EPSR=EPSRK
      END IF

      RETURN
      END


************************************************************************
*
*     * SUBROUTINE QINT *
*
*      
*     * Purpose *
*      
*     Quadratic interpolation for line search with one directional
*     derivative.
*
*
*     * Calling sequence *
*     
*     CALL QINT(TU,FL,FU,WK,T,KAPPA)
*     
*     
*     * Parameters *
*
*     RI  TU            Upper value of the stepsize parameter.
*     RI  FL            Value of the objective function.
*     RI  FU            Value of the objective function for T=TU.
*     RI  WK            Directional derivative.
*     RO  T             Stepsize parameter.
*     RI  KAPPA         Interpolation parameter.
*     RI  SMALL         A small positive value.
*
*      
*
 
     
      SUBROUTINE QINT(TU,FL,FU,WK,T,KAPPA,SMALL)
      

*     Scalar Arguments
      DOUBLE PRECISION FL,FU,WK,T,TU,KAPPA,SMALL

*     Local Scalars
      DOUBLE PRECISION TMP1,TMP2,SMALL2

*     Intrinsic Functions
      INTRINSIC MAX,SQRT


      SMALL2 = SQRT(SMALL)

      IF (WK*TU .LE. SMALL2) THEN
         T = 0.50D+00 * TU
         RETURN
      END IF
 
      TMP1 = (FU-FL) / (-WK*TU)

      
*     
*     Quadratic interpolation with one directional derivative
*

      TMP2 = 2.0D+00 * (1.0D+00 - TMP1)


      IF (TMP2 .GT. 1.0D+00) THEN

         
*     
*     Interpolation accepted
*

         T = MAX(KAPPA*TU,TU/TMP2)

         RETURN

      END IF

      
*     
*     Bisection
*

      T = 0.50D+00 * TU
      
      RETURN
      END


************************************************************************
*
*     * SUBROUTINE INDIC2 *
*
*      
*     * Purpose *
*      
*     Initialization of indices.
*     
*      
*     * Calling sequence *
*     
*     CALL INDIC2(MC,MN,INEW,IOLD,IFLAG,ITYPE)
*
*      
*     * Parameters *
*
*     II  MC              Declared number of stored corrections.
*     IU  MN              Current number of depositories used.
*     II  ITYPE           Type of Initialization:
*                             1  - corrections are stored,
*                             3  - update is skipped.
*     IU  INEW            Index for circular arrays.
*     IO  IOLD            Index of the oldest corrections.
*     IU  IFLAG           Index for adaptive version:
*                             0  - Maximum number of stored corrections
*                                    has not been changed at previous
*                                    iteration.
*                             1  - Maximum number of stored corrections
*                                    has been changed at previous
*                                    iteration.
*     
*      

      SUBROUTINE INDIC2(MC,MN,INEW,IOLD,IFLAG,ITYPE)
      
*     Scalar Arguments
      INTEGER MC,INEW,IFLAG,MN,IOLD,ITYPE
      

      IF (ITYPE .EQ. 1) THEN

         IF (IFLAG .EQ. 0) THEN

            IF (MN .LT. MC) THEN
               MN = MN + 1
               IOLD = 1
            ELSE
               IOLD = INEW + 1
               IF (IOLD .GT. MC) IOLD = 1
            END IF
               
         ELSE

            IF (MN .EQ. MC-1) THEN
               IF (INEW .EQ. 1) THEN
                  INEW = MC
                  MN = MC
                  IOLD = 1
                  IFLAG = 0

               ELSE IF (INEW .EQ. MC) THEN
                  MN = MC
                  IOLD = 1
                  IFLAG = 0

               ELSE
                  IOLD = INEW + 1
                  IF (IOLD .GT. MC-1) IOLD = 1
               
               END IF

            ELSE
               IOLD = 1
               MN = MN + 1
               IFLAG = 0
            END IF

         END IF
         
     
      ELSE 
         
         IF (IFLAG .EQ. 0) THEN
            
            IF (MN .LT. MC) THEN
               IOLD = 1
            ELSE
               IOLD = INEW
            END IF

         ELSE

            IF (MN .EQ. MC-1) THEN
               IF (INEW .EQ. 1) THEN
                  IOLD = 1
                  INEW = MC
                  IFLAG = 0
            
               ELSE IF (INEW .EQ. MC) THEN
                  IOLD = 1
                  IFLAG = 0
               
               ELSE
                  IOLD = INEW
               END IF

            ELSE
               IOLD = 1
               IFLAG = 0
            END IF
            
         END IF

      END IF

      RETURN
      END


************************************************************************
*
*     * DOUBLE PRECISION FUNCTION SCLPAR *
*
*      
*     * Purpose *
*      
*     Calculation of the scaling parameter appointed by parameter 
*     ISCALE.
*
*      
*     * Calling sequence *
*      
*      GAMMA = SCLPAR(MN,ISCALE,STS,STU,UTU,SMALL)
*
*      
*     * Parameters *
*
*     II  MN              Current number of depositories used.
*     RI  STS             STS = S'*S. 
*     RI  STU             STU = S'*U. 
*     RI  UTU             UTU = U'*U. 
*     RI  SMALL           Small positive value.
*     II  ISCALE          Selection of the scaling:
*                           0  - Scaling at every iteration
*                                with STU/UTU.
*                           1  - Scaling at every iteration
*                                with STS/STU.
*                           2  - Interval scaling with STU/UTU.
*                           3  - Interval scaling with STS/STU.
*                           4  - Preliminary scaling with STU/UTU.
*                           5  - Preliminary scaling with STS/STU.
*                           6  - No scaling.      
*
*      

      DOUBLE PRECISION FUNCTION SCLPAR(MN,ISCALE,STS,STU,UTU,
     &     SMALL)
      
*     Scalar Arguments
      INTEGER MN,ISCALE
      DOUBLE PRECISION STS,STU,UTU,SMALL
      
*     Intrinsic Functions
      INTRINSIC SQRT


*        
*     Computation of scaling parameter.
*

      SCLPAR = 1.0D+00


*
*     Scaling parameter = STU/UTU
*            

      IF (ISCALE .EQ. 0 .OR. ISCALE .EQ. 2 .OR. ISCALE .EQ. 4)
     &     THEN
         IF (UTU .LT. SQRT(SMALL)) THEN
            SCLPAR = 1.0D+00
            GO TO 80
         ELSE
            SCLPAR = STU/UTU
         END IF

         
*     
*     Scaling parameter = STS/STU
*               

      ELSE IF (ISCALE .EQ. 1 .OR. ISCALE .EQ. 3 .OR.
     &        ISCALE .EQ. 5) THEN
         IF (STU .LT. SQRT(SMALL)) THEN
            SCLPAR = 1.0D+00
            GO TO 80
         ELSE
            SCLPAR = STS/STU
         END IF
      ELSE


*     
*     No scaling
*               

         SCLPAR = 1.0D+00
         GO TO 80
      END IF

               
*     
*     Scaling
*               
            
      IF (MN .EQ. 1) THEN
         IF (SCLPAR .LT. 0.01D+00) SCLPAR=0.01D+00
         IF (SCLPAR .GT. 100.0D+00) SCLPAR=100.0D+00


*               
*     Interval scaling
*               

      ELSE IF (ISCALE .EQ. 2) THEN
         IF (SCLPAR .LT. 0.6D+00 .OR. SCLPAR .GT. 6.0D+00) THEN
            SCLPAR = 1.0D+00
         END IF

      ELSE IF (ISCALE .EQ. 3) THEN
         IF (SCLPAR .LT. 0.5D+00 .OR. SCLPAR .GT. 5.0D+00) THEN
            SCLPAR = 1.0D+00
         END IF
         
               
*     
*     Preliminary scaling
*     

      ELSE IF (ISCALE .EQ. 4 .OR. ISCALE .EQ. 5) THEN
         SCLPAR = 1.0D+00
               

*     
*     Scaling at every iteration
*               

      ELSE

         CONTINUE

      END IF

 80   CONTINUE

     

*
*     Lower bound
*

      IF (SCLPAR .LE. SQRT(SMALL)) SCLPAR = SQRT(SMALL)

      RETURN
      END


************************************************************************
*
*     * SUBROUTINE DOBUN *
*
*
*     * Purpose *
*      
*     Bundle construction for limited memory bundle method
*
*      
*     * Calling sequence *
*     
*     CALL DOBUN(N,NA,MAL,X,G,F,AY,AG,AF,ITERS,IBUN)
*     
*     
*     * Parameters *
*      
*     II  N             Number of variables.
*     II  NA            Maximum size of the bundle.
*     II  MAL           Current size of the bundle.
*     RI  X(N)          Vector of variables.
*     RI  G(N)          Subgradient of the objective function.
*     RI  F             Value of the objective function.
*     RU  AY(N*NA)      Matrix whose columns are bundle points.
*     RU  AG(N*NA)      Matrix whose columns are bundle subgradients.
*     RU  AF(4*NA)      Vector of values of bundle functions.
*     IU  IBUN          Index for the circular arrays.
*     II  ITERS         Null step indicator.
*                         0  - Null step.
*                         1  - Serious step.
*
*      
*     * Subprograms used *
*
*     S   COPY2         Copying of two vectors.
*
*      
      
      SUBROUTINE DOBUN(N,MA,MAL,X,G,F,AY,AG,AF,ITERS,IBUN)

*     Scalar Arguments
      INTEGER ITERS,MA,MAL,N,IBUN
      DOUBLE PRECISION F
      
*     Array Arguments
      DOUBLE PRECISION AF(*),AG(*),AY(*),G(*),X(*)
      
*     Local Scalars
      INTEGER I,J
      
*     External Subroutines
      EXTERNAL COPY2


      IF (ITERS .EQ. 1) THEN


*     
*     Serious step
*      

         AF(IBUN) = F
         I = (IBUN-1)*N + 1
         CALL COPY2(N,G,AG(I),X,AY(I))

      ELSE


*
*     Null step
*      

         IF (MAL .LT. MA) THEN

            AF(IBUN) = AF(MAL)
            AF(MAL) = F

            I = MAL*N + 1
            CALL COPY2(N,AG(I-N),AG(I),AY(I-N),AY(I))
            CALL COPY2(N,G,AG(I-N),X,AY(I-N))
            
         ELSE
            I = IBUN-1
            IF (I .LT. 1) I = MAL
            AF(IBUN) = AF(I)
            AF(I) = F

            I = (IBUN-2)*N + 1
            IF (I .LT. 1) I = (MAL-1)*N + 1
            J = (IBUN-1)*N + 1
            CALL COPY2(N,AG(I),AG(J),AY(I),AY(J))
            CALL COPY2(N,G,AG(I),X,AY(I))
         END IF

      END IF

      
      MAL = MAL + 1
      IF (MAL .GT. MA) MAL = MA

      IBUN = IBUN + 1
      IF (IBUN .GT. MA) IBUN = 1

      RETURN
      END


************************************************************************
*
*     * SUBROUTINE RPRINT *
*
*      
*     * Purpose *
*      
*     Printout the (final) results.
*
*      
*     * Calling sequence *
*     
*     SUBROUTINE RPRINT(N,NIT,NFE,NACT,X,F,WK,QK,ITERM,IPRINT)
*     
*     
*     * Parameters *
*      
*     II  N               Number of variables.
*     II  NIT             Number of used iterations.
*     II  NFE             Number of used function evaluations.
*     II  NACT            Number of active variables at solution.
*     RI  X(N)            Vector of variables.
*     RI  F               Value of the objective function.
*     RI  WK              Value of the first stopping criterion.
*     RI  QK              Value of the second stopping criterion.
*     II  ITERM           Cause of termination:
*                           1  - The problem has been solved.
*                                with desired accuracy.
*                           2  - (F - FO) < TOLF in MTESF
*                                subsequent iterations.
*                           3  - (F - FO) < TOLF*SMALL*MAX(|F|,|FO|,1).
*                           4  - Number of function calls > MFV.
*                           5  - Number of iterations > MIT.
*                           6  - Time limit exceeded. 
*                           7  - F < TOLB.
*                          -1  - Internal error: Two consecutive 
*                                restarts.
*                          -2  - Internal error: Number of restarts 
*                                > maximum number of restarts.
*                          -3  - Failure in function or subgradient 
*                                calculations (assigned by the user).
*                          -4  - Internal error: Failure in attaining 
*                                the demanded accuracy.
*                          -5  - Invalid input parameters.
*                          -6  - Not enough working space.
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
*                           5  - At each iteration the whole
*                                solution
*     
*      
      
      SUBROUTINE RPRINT(N,NIT,NFE,NACT,X,F,WK,QK,ITERM,IPRINT)

*     Scalar Arguments
      INTEGER N,NIT,NFE,NACT,ITERM,IPRINT,I
      DOUBLE PRECISION F,WK,QK
      
*     Array Arguments
      DOUBLE PRECISION X(*)

         
*
*     Intermediate results
*

      IF (ITERM .EQ. 0) THEN
         IF (IPRINT .GT. 3) WRITE (6,FMT='(1X,''NIT='',I5,2X,
     &        ''NFE='',I5,2X,''NACT='',I5,2X,''F='',D15.8,2X,''WK=''
     &        ,D11.2,2X,''QK='',D11.2,2X)')
     &        NIT,NFE,NACT,F,WK,QK
         IF (IPRINT .EQ. 5) WRITE (6,FMT='(1X,''X ='',
     &        5D15.7:/(4X,5D15.7))')(X(I),I=1,N)
         RETURN
      END IF
         

*
*     Printout the final results
*

      IF (IPRINT .GT. 0) WRITE (6,FMT='(1X,''NIT='',I5,2X,
     &     ''NFE='',I5,2X,''NACT='',I5,2X,''F='',D15.8,2X,''WK=''
     &     ,D11.2,2X,''QK='',D11.2,2X,''ITERM='',I3)')
     &     NIT,NFE,NACT,F,WK,QK,ITERM
      IF (IPRINT .EQ. 3 .OR. IPRINT .EQ. 5)
     &     WRITE (6,FMT='(1X,''X ='',5D15.7:/(4X,5D15.7))')(X(I),I=1
     $     ,N)
      
      RETURN
      END

      
************************************************************************
*
*     * SUBROUTINE WPRINT *
*
*      
*     * Purpose *
*      
*     Printout warning and error messages.
*
*      
*     * Calling sequence *
*     
*     SUBROUTINE WPRINT(ITERM,IPRINT,NOUT)
*     
*     
*     * Parameters *
*      
*     II  ITERM           Cause of termination:
*                           1  - The problem has been solved.
*                                with desired accuracy.
*                           2  - (F - FO) < TOLF in MTESF
*                                subsequent iterations.
*                           3  - (F - FO) < TOLF*SMALL*MAX(|F|,|FO|,1).
*                           4  - Number of function calls > MFV.
*                           5  - Number of iterations > MIT.
*                           6  - Time limit exceeded. 
*                           7  - F < TOLB.
*                          -1  - Internal error: Two consecutive 
*                                restarts.
*                          -2  - Internal error: Number of restarts 
*                                > maximum number of restarts.
*                          -3  - Failure in function or subgradient 
*                                calculations (assigned by the user).
*                          -4  - Internal error: Failure in attaining 
*                                the demanded accuracy.
*                          -5  - Invalid input parameters.
*                          -6  - Not enough working space.
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
*                           5  - At each iteration the whole
*                                solution
*     II  NOUT            Auxilary printout specification.
*     
*      
      
      SUBROUTINE WPRINT(ITERM,IPRINT,NOUT)

*     Scalar Arguments
      INTEGER ITERM,IPRINT,NOUT


      IF (IPRINT .GE. 0) THEN


*     
*     Initial error messages
*

         IF (ITERM .LE. -5) THEN

            IF (ITERM .EQ. -5) THEN
               IF (NOUT .EQ. 1) WRITE (6,FMT='(1X,''Error: ''
     &              ''Number of variables (N) is too small. ITERM='',I3)
     &              ')ITERM
               IF (NOUT .EQ. 2) WRITE (6,FMT='(1X,''Error: ''
     &              ''The maximum number of stored corrections (MCU) ''
     &              ''is too small. ITERM='',I3)')ITERM
               IF (NOUT .EQ. 3) WRITE (6,FMT='(1X,''Error: ''
     &              ''The size of the bundle (NA) is too small. ITERM=''
     &              ,I3)')ITERM
               IF (NOUT .EQ. 4) WRITE (6,FMT='(1X,''Error: ''
     &              ''Line search parameter RPAR(5) >= 0.25. ITERM=''
     &              ,I3)')ITERM
            ELSE IF (ITERM .EQ. -6) THEN
               WRITE (6,FMT='(1X,''Error: ''
     &              ''Not enough working space. ITERM='',I3)')ITERM
               
            END IF
            RETURN

         END IF

   
*
*     Warning messages
*

         IF (IPRINT .EQ. 2) THEN

            IF (ITERM .EQ. 0) THEN
               IF (NOUT .GT. 0) WRITE (6,FMT='(1X,''Warning: ''
     &              ''The initial X is infeasible. ''
     &              ''Restart with its projection.'')')
               IF (NOUT .EQ. -1) WRITE (6,FMT='(1X,''Warning: ''
     &              ''MC > MCU. Assigned MC = MCU.'')')
               IF (NOUT .EQ. -2) WRITE (6,FMT='(1X,''Warning: ''
     &              ''A line search parameter EPSR >= 0.5.'')')
               IF (NOUT .EQ. -3) WRITE (6,FMT='(1X,''Warning: ''
     &              ''A nondescent search direction occured. Restart.'')
     &              ')
               IF (NOUT .EQ. -4) WRITE (6,FMT='(1X,''Warning: ''
     &              ''Does not converge.'')')
               IF (NOUT .EQ. -5) WRITE (6,FMT='(1X,''Warning: ''
     &              ''TMAX < TMIN. Restart.'')')
               IF (NOUT .EQ. -6) WRITE (6,FMT='(1X,''Warning: ''
     &              ''An indefinite matrix detected. Restart.'')')
               IF (NOUT .EQ. -7) WRITE (6,FMT='(1X,''Warning: ''
     &              ''Insufficiently positive definite matrix ''
     &              ''detected. The convergence properties may be ''
     &              ''failed. Restart.'')')
               IF (NOUT .EQ. -8) WRITE (6,FMT='(1X,''Warning: ''
     &              ''Deficient search direction. Restart.'')')
               RETURN

            END IF

         
*
*     Printout the final results
*
        
            IF (ITERM .EQ. 6) WRITE (6,FMT='(1X,''Abnormal exit: ''
     &           ''Time is up.'')')
            IF (ITERM .EQ. 7) WRITE (6,FMT='(1X,''Abnormal exit: ''
     &           ''F < TOLB.'')')
            IF (ITERM .EQ. 2) WRITE (6,FMT='(1X,''Abnormal exit: ''
     &           ''Too many steps without significant progress.'')
     &           ')
            IF (ITERM .EQ. 3) WRITE (6,FMT='(1X,''Abnormal exit: ''
     &           ''The value of the function does not change.'')')
            IF (ITERM .EQ. 5) WRITE (6,FMT='(1X,''Abnormal exit: ''
     &           ''Number of iterations > '',I5)') NOUT
            IF (ITERM .EQ. 4) WRITE (6,FMT='(1X,''Abnormal exit: ''
     &           ''Number of function evaluations > '',I5)') NOUT
            IF (ITERM .EQ. -1) THEN
               IF (NOUT .EQ. -1) THEN
                  WRITE (6,FMT='(1X,''Abnormal exit: ''
     &                 ''Two consecutive restarts.'')')
               ELSE
                  WRITE (6,FMT='(1X,''Abnormal exit: ''
     &           ''TMAX < TMIN in two subsequent iterations.'')')
               END IF
            END IF
            IF (ITERM .EQ. -2) WRITE (6,FMT='(1X,''Abnormal exit: ''
     &                 ''Number of restarts > '',I5''.'')') NOUT
            IF (ITERM .EQ. -3) WRITE (6,FMT='(1X,''Abnormal exit: ''
     &           ''Failure in function or subgradient calculations.'')')
            IF (ITERM .EQ. -4) WRITE (6,FMT='(1X,''Abnormal exit: ''
     &           ''Failure in attaining the demanded accuracy.'')')
         END IF

      END IF

      
      RETURN
      END


************************************************************************
*
*     * SUBROUTINE GETIME *
*
*      
*     * Purpose *
*      
*     Execution time.
*     
*      
*     * Calling sequence *
*     
*     CALL GETIME(CTIM,RTIM)
*
*      
*     * Parameters *
*     
*     RO  CTIM          Current time. REAL argument
*     RA  RTIM(2)       Auxiliary array. REAL array.
*     
*     
*     * Subprograms used *
*      
*     RF  ETIME         Execution time.
*      
*

      SUBROUTINE GETIME(CTIM,RTIM)
      
*     Scalar Arguments
      REAL CTIM

*     Array arguments
      REAL RTIM(2)

*     Intrinsic Functions
      REAL ETIME
      INTRINSIC ETIME


      CTIM = ETIME(RTIM)
      CTIM = RTIM(1)
      

      RETURN
      END

      
************************************************************************
