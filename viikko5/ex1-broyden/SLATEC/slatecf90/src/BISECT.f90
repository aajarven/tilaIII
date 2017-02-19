subroutine BISECT (N, EPS1, D, E, E2, LB, UB, MM, M, W, IND, IERR, &
     RV4, RV5)
!
!! BISECT computes the eigenvalues of a symmetric tridiagonal matrix ...
!  in a given interval using Sturm sequencing.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5, D4C2A
!***TYPE      SINGLE PRECISION (BISECT-S)
!***KEYWORDS  EIGENVALUES, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the bisection technique
!     in the ALGOL procedure TRISTURM by Peters and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
!
!     This subroutine finds those eigenvalues of a TRIDIAGONAL
!     SYMMETRIC matrix which lie in a specified interval,
!     using bisection.
!
!     On INPUT
!
!        N is the order of the matrix.  N is an INTEGER variable.
!
!        EPS1 is an absolute error tolerance for the computed
!          eigenvalues.  If the input EPS1 is non-positive,
!          it is reset for each submatrix to a default value,
!          namely, minus the product of the relative machine
!          precision and the 1-norm of the submatrix.
!          EPS1 is a REAL variable.
!
!        D contains the diagonal elements of the input matrix.
!          D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the input matrix
!          in its last N-1 positions.  E(1) is arbitrary.
!          E is a one-dimensional REAL array, dimensioned E(N).
!
!        E2 contains the squares of the corresponding elements of E.
!          E2(1) is arbitrary.  E2 is a one-dimensional REAL array,
!          dimensioned E2(N).
!
!        LB and UB define the interval to be searched for eigenvalues.
!          If LB is not less than UB, no eigenvalues will be found.
!          LB and UB are REAL variables.
!
!        MM should be set to an upper bound for the number of
!          eigenvalues in the interval.  WARNING - If more than
!          MM eigenvalues are determined to lie in the interval,
!          an error return is made with no eigenvalues found.
!          MM is an INTEGER variable.
!
!     On OUTPUT
!
!        EPS1 is unaltered unless it has been reset to its
!          (last) default value.
!
!        D and E are unaltered.
!
!        Elements of E2, corresponding to elements of E regarded
!          as negligible, have been replaced by zero causing the
!          matrix to split into a direct sum of submatrices.
!          E2(1) is also set to zero.
!
!        M is the number of eigenvalues determined to lie in (LB,UB).
!          M is an INTEGER variable.
!
!        W contains the M eigenvalues in ascending order.
!          W is a one-dimensional REAL array, dimensioned W(MM).
!
!        IND contains in its first M positions the submatrix indices
!          associated with the corresponding eigenvalues in W --
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc.
!          IND is an one-dimensional INTEGER array, dimensioned IND(MM).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          3*N+1      if M exceeds MM.  In this case, M contains the
!                     number of eigenvalues determined to lie in
!                     (LB,UB).
!
!        RV4 and RV5 are one-dimensional REAL arrays used for temporary
!          storage, dimensioned RV4(N) and RV5(N).
!
!     The ALGOL procedure STURMCNT contained in TRISTURM
!     appears in BISECT in-line.
!
!     Note that subroutine TQL1 or IMTQL1 is generally faster than
!     BISECT, if more than N/4 eigenvalues are to be found.
!
!     Questions and comments should be directed to B. S. Garbow,
!     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BISECT
!
  INTEGER I,J,K,L,M,N,P,Q,R,S,II,MM,M1,M2,TAG,IERR,ISTURM
  REAL D(*),E(*),E2(*),W(*),RV4(*),RV5(*)
  REAL U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP,S1,S2
  INTEGER IND(*)
  LOGICAL FIRST
!
  SAVE FIRST, MACHEP
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BISECT
  if (FIRST) THEN
     MACHEP = R1MACH(4)
  end if
  FIRST = .FALSE.
!
  IERR = 0
  TAG = 0
  T1 = LB
  T2 = UB
!     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES ..........
  DO 40 I = 1, N
     if (I  ==  1) go to 20
     S1 = ABS(D(I)) + ABS(D(I-1))
     S2 = S1 + ABS(E(I))
     if (S2  >  S1) go to 40
   20    E2(I) = 0.0E0
   40 CONTINUE
!     .......... DETERMINE THE NUMBER OF EIGENVALUES
!                IN THE INTERVAL ..........
  P = 1
  Q = N
  X1 = UB
  ISTURM = 1
  go to 320
   60 M = S
  X1 = LB
  ISTURM = 2
  go to 320
   80 M = M - S
  if (M  >  MM) go to 980
  Q = 0
  R = 0
!     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
!                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
  100 if (R  ==  M) go to 1001
  TAG = TAG + 1
  P = Q + 1
  XU = D(P)
  X0 = D(P)
  U = 0.0E0
!
  DO 120 Q = P, N
     X1 = U
     U = 0.0E0
     V = 0.0E0
     if (Q  ==  N) go to 110
     U = ABS(E(Q+1))
     V = E2(Q+1)
  110    XU = MIN(D(Q)-(X1+U),XU)
     X0 = MAX(D(Q)+(X1+U),X0)
     if (V  ==  0.0E0) go to 140
  120 CONTINUE
!
  140 X1 = MAX(ABS(XU),ABS(X0)) * MACHEP
  if (EPS1  <=  0.0E0) EPS1 = -X1
  if (P  /=  Q) go to 180
!     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
  if (T1  >  D(P) .OR. D(P)  >=  T2) go to 940
  M1 = P
  M2 = P
  RV5(P) = D(P)
  go to 900
  180 X1 = X1 * (Q-P+1)
  LB = MAX(T1,XU-X1)
  UB = MIN(T2,X0+X1)
  X1 = LB
  ISTURM = 3
  go to 320
  200 M1 = S + 1
  X1 = UB
  ISTURM = 4
  go to 320
  220 M2 = S
  if (M1  >  M2) go to 940
!     .......... FIND ROOTS BY BISECTION ..........
  X0 = UB
  ISTURM = 5
!
  DO 240 I = M1, M2
     RV5(I) = UB
     RV4(I) = LB
  240 CONTINUE
!     .......... LOOP FOR K-TH EIGENVALUE
!                FOR K=M2 STEP -1 UNTIL M1 DO --
!                (-DO- NOT USED TO LEGALIZE -COMPUTED go to-) ..........
  K = M2
  250    XU = LB
!     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
     DO 260 II = M1, K
        I = M1 + K - II
        if (XU  >=  RV4(I)) go to 260
        XU = RV4(I)
        go to 280
  260    CONTINUE
!
  280    if (X0  >  RV5(K)) X0 = RV5(K)
!     .......... NEXT BISECTION STEP ..........
  300    X1 = (XU + X0) * 0.5E0
     S1 = 2.0E0*(ABS(XU) + ABS(X0) + ABS(EPS1))
     S2 = S1 + ABS(X0 - XU)
     if (S2  ==  S1) go to 420
!     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  320    S = P - 1
     U = 1.0E0
!
     DO 340 I = P, Q
        if (U  /=  0.0E0) go to 325
        V = ABS(E(I)) / MACHEP
        if (E2(I)  ==  0.0E0) V = 0.0E0
        go to 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
        if (U  <  0.0E0) S = S + 1
  340    CONTINUE
!
     go to (60,80,200,220,360), ISTURM
!     .......... REFINE INTERVALS ..........
  360    if (S  >=  K) go to 400
     XU = X1
     if (S  >=  M1) go to 380
     RV4(M1) = X1
     go to 300
  380    RV4(S+1) = X1
     if (RV5(S)  >  X1) RV5(S) = X1
     go to 300
  400    X0 = X1
     go to 300
!     .......... K-TH EIGENVALUE FOUND ..........
  420    RV5(K) = X1
  K = K - 1
  if (K  >=  M1) go to 250
!     .......... ORDER EIGENVALUES TAGGED WITH THEIR
!                SUBMATRIX ASSOCIATIONS ..........
  900 S = R
  R = R + M2 - M1 + 1
  J = 1
  K = M1
!
  DO 920 L = 1, R
     if (J  >  S) go to 910
     if (K  >  M2) go to 940
     if (RV5(K)  >=  W(L)) go to 915
!
     DO 905 II = J, S
        I = L + S - II
        W(I+1) = W(I)
        IND(I+1) = IND(I)
  905    CONTINUE
!
  910    W(L) = RV5(K)
     IND(L) = TAG
     K = K + 1
     go to 920
  915    J = J + 1
  920 CONTINUE
!
  940 if (Q  <  N) go to 100
  go to 1001
!     .......... SET ERROR -- UNDERESTIMATE OF NUMBER OF
!                EIGENVALUES IN INTERVAL ..........
  980 IERR = 3 * N + 1
 1001 LB = T1
  UB = T2
  return
end
