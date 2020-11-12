************************************************************************
*
*
*     Test problems for NonSmooth BOund Constrained minimization
*
*
*     TNSBOC includes the following subroutines
*
*     S   STARTX          Initiation of variables (not necessarily 
*                           feasible).
*     S   BOUNDS          Bound constraints.
*     S   FEASIX          Projection of starting points to feasible
*                           region.
*     S   FUNC            Computation of the value and the subgradient 
*                           of the objective function.
*
*
*     Napsu Karmitsa (2003, bound constrained version 2006-2007)      
*
*     Haarala M., Miettinen K. and Mäkelä M.M.: New Limited Memory
*     Bundle Method for Large-Scale Nonsmooth Optimization, Optimization
*     Methods and Software, Vol. 19, No. 6, 2004, 673-692.
*
*     Karmitsa N.: Test Problems for Large-Scale Nonsmooth Minimization,
*     Reports of the Department of Mathematical Information Technology,
*     Series B, Scientific Computing, B 4/2007, University of Jyväskylä, 
*     Jyväskylä, 2007.
*
*
************************************************************************
*      
*     * SUBROUTINE STARTX *
*
*      
*     * Purpose *
*
*
*     Initiation of X.
*
*      
*     * Calling sequence *
*
*     CALL STARTX(N,X,NEXT)
*
*
*     * Parameters *
*      
*     II  N          Number of variables.
*     RO  X(N)       Vector of variables.
*     RI  NEXT       Problem number.
*
*
*     * Problems *
*      
*     1.  Generalization of MAXQ (convex).
*     2.  Generalization of MXHILB (convex).
*     3.  Chained LQ (convex).
*     4.  Chained CB3 (convex).
*     5.  Chained CB3 2 (convex).
*     6.  Number of active faces (nonconvex).
*     7.  Nonsmooth generalization of Brown function 2 (nonconvex).
*     8.  Chained Mifflin 2 (nonconvex). 
*     9.  Chained crescent (nonconvex). 
*     10. Chained crescent 2 (nonconvex).
*     
*
*     Napsu Haarala (2003)      
*
*     Haarala M., Miettinen K. and Mäkelä M.M.: New Limited Memory
*     Bundle Method for Large-Scale Nonsmooth Optimization, Optimization
*     Methods and Software, Vol. 19, No. 6, 2004, 673-692.. 
*      
      
      SUBROUTINE STARTX(N,X,NEXT)

*     Scalar Arguments
      INTEGER N,NEXT

*     Array Arguments
      DOUBLE PRECISION X(*)

*     Local Arguments
      INTEGER I

      GOTO(10,20,30,40,40,60,70,80,90,90) NEXT


      PRINT*,'Error: Not such a problem.'
      NEXT=-1
      RETURN

      
*    
*     Generalization of MAXQ (convex)
*
 10   CONTINUE
      DO 11 I=1,N/2
         X(I) = DBLE(I)
 11   CONTINUE
      DO 12 I=N/2+1,N
         X(I) = -DBLE(I)
 12   CONTINUE

      RETURN

      
*
*     Generalization of MXHILB (convex)
*      
 20   CONTINUE
      DO 21 I=1,N
         X(I) = 1.0D+00
 21   CONTINUE

      RETURN

      
*     
*     Chained LQ (convex)
*     
 30   CONTINUE
      DO 31 I=1,N
         X(I)=-0.50D+00
 31   CONTINUE

      RETURN
      

*     
*     Chained CB3 1 and 2 (convex)
*
 40   CONTINUE
      DO 41 I=1,N
         X(I)=2.0D+00
 41   CONTINUE

      RETURN

      
*
*     Number of active faces (nonconvex)
*      
 60   CONTINUE
      DO 61 I=1,N
         X(I) = 1.0D+00
 61   CONTINUE
      
      RETURN


*      
*     Nonsmooth generalization of Brown function 2 (nonconvex)
*      
 70   CONTINUE
      DO 71 I=1,N
         IF(MOD(I,2) .EQ. 1) THEN
            X(I) = -1.0D+00
         ELSE
            X(I) = 1.0D+00
         ENDIF
 71   CONTINUE

      RETURN

      
*      
*     Chained Mifflin 2 (nonconvex)
*      
 80   CONTINUE
      DO 81 I=1,N
         X(I) = -1.0D+00
 81   CONTINUE

      RETURN


*      
*     Chained crescent (nonconvex)
*      
 90   CONTINUE
      DO 91 I=1,N
         IF(MOD(I,2) .EQ. 1) THEN
            X(I) = -1.50D+00
         ELSE
            X(I) = 2.0D+00
         ENDIF
 91   CONTINUE

      RETURN

      END


************************************************************************
*      
*     * SUBROUTINE BOUNDS *
*
*      
*     * Purpose *
*
*
*     Defination of bound constraint.
*
*      
*     * Calling sequence *
*
*     CALL BOUNDS(N,IB,XL,XU,NEXT)
*
*
*     * Parameters *
*      
*     II  N          Number of variables.
*     IO  IB(N)      Type of bound constraints:
*                      0  - X(I) is unbounded,
*                      1  - X(I) has only a lower bound,
*                      2  - X(I) has both lower and upper bounds,
*                      3  - X(I) has only an upper bound. 
*     RO  XL(N)      Lower bounds for variables. 
*     RO  XU(N)      Upper bounds for variables. 
*     RI  NEXT       Problem number.
*
*
*     * Problems *
*      
*     1.  Generalization of MAXQ (convex).
*     2.  Generalization of MXHILB (convex).
*     3.  Chained LQ (convex).
*     4.  Chained CB3 (convex).
*     5.  Chained CB3 2 (convex).
*     6.  Number of active faces (nonconvex).
*     7.  Nonsmooth generalization of Brown function 2 (nonconvex).
*     8.  Chained Mifflin 2 (nonconvex). 
*     9.  Chained crescent (nonconvex). 
*     10. Chained crescent 2 (nonconvex).
*     
*
*     Napsu Haarala (2007)      
*
      
      SUBROUTINE BOUNDS(N,IB,XL,XU,NEXT)

*     Scalar Arguments
      INTEGER N,NEXT

*     Array Arguments
      INTEGER IB(*)
      DOUBLE PRECISION XL(*),XU(*)

*     Local Arguments
      INTEGER I


      DO I=1,N
         IB(I)=0
         XL(I)=0.0D+00
         XU(I)=0.0D+00
      END DO

      GOTO(10,10,30,20,20,10,10,30,10,10) NEXT


      PRINT*,'Error: Not such a problem.'
      NEXT=-1
      RETURN


 10   CONTINUE
      
      DO 11 I=2,N,2
         IB(I)=2
         XL(I)=0.1D+00
         XU(I)=1.1D+00
 11   CONTINUE

      RETURN


 20   CONTINUE

      DO 21 I=2,N,2
         IB(I)=2
         XL(I)=1.1D+00
         XU(I)=2.1D+00
 21   CONTINUE

      RETURN


 30   CONTINUE

      DO 31 I=2,N,2
         IB(I)=2
         XL(I)=1.0D+00/SQRT(2.0D+00)+0.1D+00
         XU(I)=1.0D+00/SQRT(2.0D+00)+1.1D+00
 31   CONTINUE

      IF (NEXT .EQ. 8) THEN
         XL(2)=0.68D+00
         XU(2)=1.68D+00
         XL(N)=0.1D+00
         XU(N)=1.1D+00
      END IF
      

      RETURN

      END


************************************************************************
*      
*     * SUBROUTINE FEASIX *
*
*      
*     * Purpose *
*
*     Projection of the initial X to the feasible region.
*
*      
*     * Calling sequence *
*     
*     CALL FEASIX(N,X,XL,XU,IB,ITYPE)
*
*      
*     * Parameters *
*     
*
*     II  N               Number of variables.
*     RU  X(N)            Vector of variables. 
*     RI  XL(N)           Lower bounds for variables. 
*     RI  XU(N)           Upper bounds for variables.
*     II  IB(N)           Type of bound constraints:
*                             0  - X(I) is unbounded,
*                             1  - X(I) has only a lower bound,
*                             2  - X(I) has both lower and upper bounds,
*                             3  - X(I) has only an upper bound. 
*     II  ITYPE           Type of starting point needed:
*                             0  - feasible,
*                             1  - strictly feasible.
*     
*      
*     Napsu Karmitsa (2007)
*
     
      
      SUBROUTINE FEASIX(N,X,XL,XU,IB,ITYPE)
      
*     Scalar Arguments
      INTEGER N,ITYPE
      
*     Array Arguments
      INTEGER IB(*)
      DOUBLE PRECISION X(*),XL(*),XU(*)

*     Local Scalars
      INTEGER I
      DOUBLE PRECISION FEAS

      FEAS = 1.0D-04


c     Project the initial X to the feasible set.

      IF (ITYPE .EQ. 0) THEN
         DO 10 I=1,N
            IF (IB(I) .GT. 0) THEN
               IF (IB(I) .LE. 2) THEN
                  IF (X(I) .LT. XL(I)) THEN
                     X(I)=XL(I)
                  END IF
               END IF
               IF (IB(I) .GE. 2) THEN
                  IF (X(I) .GT. XU(I)) THEN
                     X(I)=XU(I)
                  END IF
               END IF
            END IF
 10      CONTINUE
      ELSE 

         DO 20 I=1,N
            IF (IB(I) .GT. 0) THEN
               IF (IB(I) .LE. 2) THEN
                  IF (X(I) .LE. XL(I)) THEN
                     X(I)=XL(I)+FEAS
                  END IF
               END IF
               IF (IB(I) .GE. 2) THEN
                  IF (X(I) .GE. XU(I)) THEN
                     X(I)=XU(I)-FEAS
                  END IF
               END IF
            END IF
 20      CONTINUE
      END IF
         
      RETURN
      END
      
      
************************************************************************
*      
*     * SUBROUTINE FUNC *
*
*      
*     * Purpose *
*
*
*     Computation of the value and the subgradient of the objective
*     function.
*
*      
*     * Calling sequence *
*
*     CALL FUNC(N,X,F,G,NEXT)
*
*
*     * Parameters *
*      
*     II  N          Number of variables.
*     RI  X(N)       Vector of variables.
*     RI  NEXT       Problem number.
*     RO  F          Value of the objective function.
*     RO  G(N)       Subgradient of the objective function.
*
*
*     * Problems *
*      
*     1.  Generalization of MAXQ (convex).
*     2.  Generalization of MXHILB (convex).
*     3.  Chained LQ (convex).
*     4.  Chained CB3 (convex).
*     5.  Chained CB3 2 (convex).
*     6.  Number of active faces (nonconvex).
*     7.  Nonsmooth generalization of Brown function 2 (nonconvex).
*     8.  Chained Mifflin 2 (nonconvex). 
*     9.  Chained crescent (nonconvex). 
*     10. Chained crescent 2 (nonconvex).
*     
*
*     Napsu Haarala (2003)      
*
*     Haarala M., Miettinen K. and Mäkelä M.M.: New Limited Memory
*     Bundle Method for Large-Scale Nonsmooth Optimization, Optimization
*     Methods and Software, Vol. 19, No. 6, 2004, 673-692.. 
*      
*     
      
      SUBROUTINE FUNC(N,X,F,G,NEXT)

*     Scalar Arguments
      INTEGER N,NEXT
      DOUBLE PRECISION F

*     Array Arguments
      DOUBLE PRECISION G(*),X(*)

*     Local Arguments
      INTEGER I,J,HIT
      DOUBLE PRECISION Y,TEMP2,TEMP3,A,B,C,D,P,Q

*     Intrinsic Functions
      INTRINSIC DABS,DMAX1,SIGN,DLOG,DEXP,DCOS,DSIN

      
      GOTO(10,20,30,40,50,60,70,80,90,100) NEXT

      
      PRINT*,'Error: Not such a problem.'
      NEXT=-1
      RETURN

      
*      
*     Generalization of MAXQ (convex)
*      
 10   CONTINUE
      F=X(1)*X(1)
      G(1)=0.0D+00
      HIT=1
      DO 11 I=2,N
         Y=X(I)*X(I)
         IF (Y .GT. F) THEN
            F=Y
            HIT=I
         END IF
         G(I)=0.0D+00
 11   CONTINUE

      G(HIT)=2*X(HIT)

      RETURN
      
      
*      
*     Generalization of MXHILB (convex)
*      
 20   CONTINUE
      F = 0.0D+00
      HIT=1
      DO 21 J = 1,N
         F = F + X(J)/DBLE(J)
 21   CONTINUE
      G(1)=SIGN(1.0D+00,F)
      F = DABS(F)
      DO 22 I = 2,N
         TEMP2 = 0.0D0
         DO 23 J = 1,N
            TEMP2 = TEMP2 + X(J)/DBLE(I+J-1)
 23      CONTINUE
         G(I)=SIGN(1.0D+00,TEMP2)
         TEMP2 = DABS(TEMP2)
         IF (TEMP2 .GT. F) THEN
            F=TEMP2
            HIT=I
         END IF
 22   CONTINUE
      TEMP3=G(HIT)
      DO 24 J = 1,N
         G(J) = TEMP3/DBLE(HIT+J-1)
 24   CONTINUE

      RETURN
      
      
*     
*     Chained LQ (convex)
*     
 30   CONTINUE
      F=0.0D+00
      G(1)=0.0D+00

      DO 31 I=1,N-1
         G(I+1)=0.0D+00
         A = -X(I)-X(I+1)
         B = -X(I)-X(I+1)+(X(I)*X(I)+X(I+1)*X(I+1)-1.0D+00)
         IF (A .GE. B) THEN
            F=F+A
            G(I)=G(I)-1.0D+00
            G(I+1)=-1.0D+00
         ELSE
            F=F+B
            G(I)=G(I)-1.0D+00+2.0D+00*X(I)
            G(I+1)=-1.0D+00+2.0D+00*X(I+1)
         ENDIF
 31   CONTINUE
      
      RETURN


*     
*     Chained CB3 (convex)
*     
 40   CONTINUE
      F=0.0D+00
      G(1)=0.0D+00

      DO 41 I=1,N-1
         G(I+1)=0.0D+00
         A=X(I)*X(I)*X(I)*X(I)+X(I+1)*X(I+1)
         B=(2.0D+00-X(I))*(2.0D+00-X(I))+
     &        (2.0D+00-X(I+1))*(2.0D+00-X(I+1))
         C= 2.0D+00*DEXP(-X(I)+X(I+1))
         Y=DMAX1(A,B)
         Y=DMAX1(Y,C)
         IF (Y .EQ. A) THEN
            G(I)=G(I)+4.0D+00*X(I)*X(I)*X(I)
            G(I+1)=2.0D+00*X(I+1)
         ELSE IF (Y .EQ. B) THEN
            G(I)=G(I)+2.0D+00*X(I)-4.0D+00
            G(I+1)=2.0D+00*X(I+1)-4.0D+00
         ELSE
            G(I)= G(I) - C
            G(I+1)= C
         END IF
         F=F+Y
 41   CONTINUE

      RETURN
      
      
*     
*     Chained CB3 2 (convex)
*
 50   CONTINUE
      F=0.0D+00
      G(1)=0.0D+00
      A=0.0D+00
      B=0.0D+00
      C=0.0D+00

      DO 51 I=1,N-1
         G(I+1)=0.0D+00
         A=A+X(I)*X(I)*X(I)*X(I)+X(I+1)*X(I+1)
         B=B+(2.0D+00-X(I))*(2.0D+00-X(I))+
     &        (2.0D+00-X(I+1))*(2.0D+00-X(I+1))
         C=C+2.0D+00*DEXP(-X(I)+X(I+1))
 51   CONTINUE
      F=DMAX1(A,B)
      F=DMAX1(F,C)
      IF (F .EQ. A) THEN
         DO 53 I=1,N-1
            G(I)=G(I)+4.0D+00*X(I)*X(I)*X(I)
            G(I+1)=2.0D+00*X(I+1)
 53      CONTINUE
      ELSE IF (F .EQ. B) THEN
         DO 54 I=1,N-1
            G(I)=G(I)+2.0D+00*X(I)-4.0D+00
            G(I+1)=2.0D+00*X(I+1)-4.0D+00
 54      CONTINUE
      ELSE
         DO 55 I=1,N-1
            G(I)= G(I) - 2.0D+00*DEXP(-X(I)+X(I+1))
            G(I+1)= 2.0D+00*DEXP(-X(I)+X(I+1))
 55      CONTINUE
      END IF

      RETURN
      
      
*      
*     Number of active faces (nonconvex)
*      
 60   CONTINUE

      TEMP3=1.0D+00
      Y=-X(1)
      G(1)= 0.0D+00
      F=DLOG(DABS(X(1))+1.0D+00)
      HIT=1
      TEMP2=F
      DO 62 I=2,N
         Y=Y - X(I)
         G(I)= 0.0D+00
         F=DMAX1(F,DLOG(DABS(X(I))+1.0D+00))
         IF(F .GT. TEMP2) THEN
            HIT=I
            TEMP2=F
         END IF
 62   CONTINUE
      F=DMAX1(F,DLOG(DABS(Y)+1.0D+00))
      IF(F .GT. TEMP2) THEN
         IF (Y.GE.0.0D+00) TEMP3=-1.0D+00
         DO 63 I=1,N
            G(I)= TEMP3*(1.0D+00/(DABS(Y)+1.0D+00))
 63      CONTINUE
      ELSE
         IF (X(HIT).LT.0.0D+00) TEMP3=-1.0D+00
         G(HIT)=TEMP3*(1.0D+00/(DABS(X(HIT))+1.0D+00))
      END IF
      
      RETURN
      
      
*
*     Nonsmooth generalization of Brown function 2 (nonconvex)
*      
 70   CONTINUE
      F=0.0D+00
      G(1)=0.0D+00
      DO 71 I=1,N-1
         A=DABS(X(I))
         B=DABS(X(I+1))
         C=X(I)*X(I)+1.0D+00
         D=X(I+1)*X(I+1)+1.0D+00
         F=F+B**C+A**D

         P=0.0D+00
         Q=0.0D+00
         IF (X(I).LT.0.0D+00) THEN
            IF (B .GT. P) P=DLOG(B)
            G(I)=G(I)-D*A**(D-1.0D+00)+2.0D+00*X(I)*P*B**C
         ELSE
            IF (B .GT. P) P=DLOG(B)
            G(I)=G(I)+D*A**(D-1.0D+00)+2.0D+00*X(I)*P*B**C
         ENDIF

         IF (X(I+1).EQ.0.0D+00) THEN
            G(I+1)=0.0D+00
         ELSE IF (X(I+1).LT.0.0D+00) THEN
            IF (A .GT. Q) Q=DLOG(A)
            G(I+1)=-C*B**(C-1.0D+00)+2.0D+00*X(I+1)*Q*A**D
         ELSE
            IF (A .GT. Q) Q=DLOG(A)
            G(I+1)=C*B**(C-1.0D+00)+2.0D+00*X(I+1)*Q*A**D
         ENDIF
 71   CONTINUE

      RETURN

      
*
*     Chained mifflin 2 (nonconvex)
*      
 80   CONTINUE
      F=0.0D+00
      G(1)=0.0D+00
      DO 81 I=1,N-1
         Y = X(I)*X(I) + X(I+1)*X(I+1) - 1.0D0
         F = F -X(I) + 2.0D+00*Y + 1.75D+00*DABS(Y)
         Y = SIGN(3.5D+00,Y) + 4.0D+00
         G(I) = G(I) + Y*X(I) - 1.0D+00
         G(I+1) = Y*X(I+1)
 81   CONTINUE
      
      RETURN

      
*
*     Chained crescent (nonconvex)
*
 90   CONTINUE
      TEMP2=0.0D+00
      TEMP3=0.0D+00
      DO 91 I=1,N-1
         TEMP2 = TEMP2 + X(I)*X(I) + (X(I+1)-1.0D+00)*(X(I+1)-1.0D+00)
     &        + X(I+1) - 1.0D+00
         TEMP3 = TEMP3 - X(I)*X(I) - (X(I+1)-1.0D+00)*(X(I+1)-1.0D+00)
     &        + X(I+1) + 1.0D+00
 91   CONTINUE
      F = DMAX1(TEMP2,TEMP3)

      G(1)=0.0D+00
      IF (TEMP2 .GE. TEMP3) THEN
         DO 92 I=1,N-1
            G(I)=G(I)+2.0D+00*X(I)
            G(I+1)=2.0D+00*(X(I+1)-1.0D+00) + 1.0D+00
 92      CONTINUE
      ELSE
         DO 93 I=1,N-1
            G(I)=G(I)-2.0D+00*X(I)
            G(I+1)=-2.0D+00*(X(I+1)-1.0D+00) + 1.0D+00
 93      CONTINUE
      END IF
      
      RETURN


*     
*     Chained crescent 2 (nonconvex)
*
 100  CONTINUE
      F=0.0D+00
      G(1)=0.0D+00
      
      DO 101 I=1,N-1
         TEMP2 =  X(I)*X(I) + (X(I+1)-1.0D+00)*(X(I+1)-1.0D+00)
     &        + X(I+1) - 1.0D+00
         TEMP3 =  - X(I)*X(I) - (X(I+1)-1.0D+00)*(X(I+1)-1.0D+00)
     &        + X(I+1) + 1.0D+00
         IF (TEMP2 .GE. TEMP3) THEN
            F=F+TEMP2
            G(I)=G(I)+2.0D+00*X(I)
            G(I+1)=2.0D+00*(X(I+1)-1.0D+00) + 1.0D+00
         ELSE
            F=F+TEMP3
            G(I)=G(I)-2.0D+00*X(I)
            G(I+1)=-2.0D+00*(X(I+1)-1.0D+00) + 1.0D+00
         END IF
 101  CONTINUE

      RETURN
      
      END
      
