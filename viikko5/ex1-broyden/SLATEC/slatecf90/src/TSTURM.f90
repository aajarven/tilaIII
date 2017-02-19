subroutine TSTURM (NM, N, EPS1, D, E, E2, LB, UB, MM, M, W, Z, &
     IERR, RV1, RV2, RV3, RV4, RV5, RV6)
!
!! TSTURM finds those eigenvalues of a symmetric tridiagonal matrix ...
!            in a given interval and their associated eigenvectors by
!            Sturm sequencing.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5, D4C2A
!***TYPE      SINGLE PRECISION (TSTURM-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine finds those eigenvalues of a TRIDIAGONAL
!     SYMMETRIC matrix which lie in a specified interval and their
!     associated eigenvectors, using bisection and inverse iteration.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, Z, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        EPS1 is an absolute error tolerance for the computed eigen-
!          values.  It should be chosen so that the accuracy of these
!          eigenvalues is commensurate with relative perturbations of
!          the order of the relative machine precision in the matrix
!          elements.  If the input EPS1 is non-positive, it is reset
!          for each submatrix to a default value, namely, minus the
!          product of the relative machine precision and the 1-norm of
!          the submatrix.  EPS1 is a REAL variable.
!
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the symmetric
!          tridiagonal matrix in its last N-1 positions.  E(1) is
!          arbitrary.  E is a one-dimensional REAL array, dimensioned
!          E(N).
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
!          eigenvalues in the interval.  MM is an INTEGER variable.
!          WARNING -  If more than MM eigenvalues are determined to lie
!          in the interval, an error return is made with no values or
!          vectors found.
!
!     On Output
!
!        EPS1 is unaltered unless it has been reset to its
!          (last) default value.
!
!        D and E are unaltered.
!
!        Elements of E2, corresponding to elements of E regarded as
!          negligible, have been replaced by zero causing the matrix to
!          split into a direct sum of submatrices.  E2(1) is also set
!          to zero.
!
!        M is the number of eigenvalues determined to lie in (LB,UB).
!          M is an INTEGER variable.
!
!        W contains the M eigenvalues in ascending order if the matrix
!          does not split.  If the matrix splits, the eigenvalues are
!          in ascending order for each submatrix.  If a vector error
!          exit is made, W contains those values already found.  W is a
!          one-dimensional REAL array, dimensioned W(MM).
!
!        Z contains the associated set of orthonormal eigenvectors.
!          If an error exit is made, Z contains those vectors already
!          found.  Z is a one-dimensional REAL array, dimensioned
!          Z(NM,MM).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          3*N+1      if M exceeds MM no eigenvalues or eigenvectors
!                     are computed,
!          4*N+J      if the eigenvector corresponding to the J-th
!                     eigenvalue fails to converge in 5 iterations, then
!                     the eigenvalues and eigenvectors in W and Z should
!                     be correct for indices 1, 2, ..., J-1.
!
!        RV1, RV2, RV3, RV4, RV5, and RV6 are temporary storage arrays,
!          dimensioned RV1(N), RV2(N), RV3(N), RV4(N), RV5(N), and
!          RV6(N).
!
!     The ALGOL procedure STURMCNT contained in TRISTURM
!     appears in TSTURM in-line.
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
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
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  TSTURM
!
  INTEGER I,J,K,M,N,P,Q,R,S,II,IP,JJ,MM,M1,M2,NM,ITS
  INTEGER IERR,GROUP,ISTURM
  REAL D(*),E(*),E2(*),W(*),Z(NM,*)
  REAL RV1(*),RV2(*),RV3(*),RV4(*),RV5(*),RV6(*)
  REAL U,V,LB,T1,T2,UB,UK,XU,X0,X1,EPS1,EPS2,EPS3,EPS4
  REAL NORM,MACHEP,S1,S2
  LOGICAL FIRST
!
  SAVE FIRST, MACHEP
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  TSTURM
  if (FIRST) THEN
     MACHEP = R1MACH(4)
  end if
  FIRST = .FALSE.
!
  IERR = 0
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
  R = R + 1
!
  DO 160 I = 1, N
  160 Z(I,R) = 0.0E0
!
  W(R) = D(P)
  Z(P,R) = 1.0E0
  go to 940
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
!     .......... FIND VECTORS BY INVERSE ITERATION ..........
  NORM = ABS(D(P))
  IP = P + 1
!
  DO 500 I = IP, Q
  500 NORM = MAX(NORM, ABS(D(I)) + ABS(E(I)))
!     .......... EPS2 IS THE CRITERION FOR GROUPING,
!                EPS3 REPLACES ZERO PIVOTS AND EQUAL
!                ROOTS ARE MODIFIED BY EPS3,
!                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
  EPS2 = 1.0E-3 * NORM
  UK = SQRT(REAL(Q-P+5))
  EPS3 = UK * MACHEP * NORM
  EPS4 = UK * EPS3
  UK = EPS4 / SQRT(UK)
  GROUP = 0
  S = P
!
  DO 920 K = M1, M2
     R = R + 1
     ITS = 1
     W(R) = RV5(K)
     X1 = RV5(K)
!     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
     if (K  ==  M1) go to 520
     if (X1 - X0  >=  EPS2) GROUP = -1
     GROUP = GROUP + 1
     if (X1  <=  X0) X1 = X0 + EPS3
!     .......... ELIMINATION WITH INTERCHANGES AND
!                INITIALIZATION OF VECTOR ..........
  520    V = 0.0E0
!
     DO 580 I = P, Q
        RV6(I) = UK
        if (I  ==  P) go to 560
        if (ABS(E(I))  <  ABS(U)) go to 540
        XU = U / E(I)
        RV4(I) = XU
        RV1(I-1) = E(I)
        RV2(I-1) = D(I) - X1
        RV3(I-1) = 0.0E0
        if (I  /=  Q) RV3(I-1) = E(I+1)
        U = V - XU * RV2(I-1)
        V = -XU * RV3(I-1)
        go to 580
  540       XU = E(I) / U
        RV4(I) = XU
        RV1(I-1) = U
        RV2(I-1) = V
        RV3(I-1) = 0.0E0
  560       U = D(I) - X1 - XU * V
        if (I  /=  Q) V = E(I+1)
  580    CONTINUE
!
     if (U  ==  0.0E0) U = EPS3
     RV1(Q) = U
     RV2(Q) = 0.0E0
     RV3(Q) = 0.0E0
!     .......... BACK SUBSTITUTION
!                FOR I=Q STEP -1 UNTIL P DO -- ..........
  600    continue

  DO II = P, Q
        I = P + Q - II
        RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
        V = U
        U = RV6(I)
  end do
!
!  ORTHOGONALIZE WITH RESPECT TO PREVIOUS MEMBERS OF GROUP.
!
     DO JJ = 1, GROUP

        J = R - GROUP - 1 + JJ
        XU = 0.0E0

        DO I = P, Q
          XU = XU + RV6(I) * Z(I,J)
        end do

        DO I = P, Q
          RV6(I) = RV6(I) - XU * Z(I,J)
        end do

     end do

  700    NORM = 0.0E0
!
     DO 720 I = P, Q
  720    NORM = NORM + ABS(RV6(I))
!
     if (NORM  >=  1.0E0) go to 840
!     .......... FORWARD SUBSTITUTION ..........
     if (ITS  ==  5) go to 960
     if (NORM  /=  0.0E0) go to 740
     RV6(S) = EPS4
     S = S + 1
     if (S  >  Q) S = P
     go to 780
  740    XU = EPS4 / NORM
!
     DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
!     .......... ELIMINATION OPERATIONS ON NEXT VECTOR ITERATE.
!
  780    DO 820 I = IP, Q
        U = RV6(I)
!     .......... if RV1(I-1)  ==  E(I), A ROW INTERCHANGE
!                WAS PERFORMED EARLIER IN THE
!                TRIANGULARIZATION PROCESS ..........
        if (RV1(I-1)  /=  E(I)) go to 800
        U = RV6(I-1)
        RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
!
     ITS = ITS + 1
     go to 600
!     .......... NORMALIZE SO THAT SUM OF SQUARES IS
!                1 AND EXPAND TO FULL ORDER ..........
  840    U = 0.0E0
!
     DO 860 I = P, Q
  860    U = U + RV6(I)**2
!
     XU = 1.0E0 / SQRT(U)
!
     DO 880 I = 1, N
  880    Z(I,R) = 0.0E0
!
     DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
!
     X0 = X1
  920 CONTINUE
!
  940 if (Q  <  N) go to 100
  go to 1001
!     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
  960 IERR = 4 * N + R
  go to 1001
!     .......... SET ERROR -- UNDERESTIMATE OF NUMBER OF
!                EIGENVALUES IN INTERVAL ..........
  980 IERR = 3 * N + 1
 1001 LB = T1
  UB = T2
  return
end
