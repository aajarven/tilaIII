subroutine RATQR (N, EPS1, D, E, E2, M, W, IND, BD, TYPE, IDEF, IERR )
!
!! RATQR computes the largest or smallest eigenvalues of a symmetric...
!            tridiagonal matrix using the rational QR method with Newton
!            correction.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5, D4C2A
!***TYPE      SINGLE PRECISION (RATQR-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure RATQR,
!     NUM. MATH. 11, 264-272(1968) by REINSCH and BAUER.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 257-265(1971).
!
!     This subroutine finds the algebraically smallest or largest
!     eigenvalues of a SYMMETRIC TRIDIAGONAL matrix by the
!     rational QR method with Newton corrections.
!
!     On Input
!
!        N is the order of the matrix.  N is an INTEGER variable.
!
!        EPS1 is a theoretical absolute error tolerance for the
!          computed eigenvalues.  If the input EPS1 is non-positive, or
!          indeed smaller than its default value, it is reset at each
!          iteration to the respective default value, namely, the
!          product of the relative machine precision and the magnitude
!          of the current eigenvalue iterate.  The theoretical absolute
!          error in the K-th eigenvalue is usually not greater than
!          K times EPS1.  EPS1 is a REAL variable.
!
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the symmetric
!          tridiagonal matrix in its last N-1 positions.  E(1) is
!          arbitrary.  E is a one-dimensional REAL array, dimensioned
!          E(N).
!
!        E2 contains the squares of the corresponding elements of E in
!          its last N-1 positions.  E2(1) is arbitrary.  E2 is a one-
!          dimensional REAL array, dimensioned E2(N).
!
!        M is the number of eigenvalues to be found.  M is an INTEGER
!          variable.
!
!        IDEF should be set to 1 if the input matrix is known to be
!          positive definite, to -1 if the input matrix is known to
!          be negative definite, and to 0 otherwise.  IDEF is an
!          INTEGER variable.
!
!        TYPE should be set to .TRUE. if the smallest eigenvalues are
!          to be found, and to .FALSE. if the largest eigenvalues are
!          to be found.  TYPE is a LOGICAL variable.
!
!     On Output
!
!        EPS1 is unaltered unless it has been reset to its
!          (last) default value.
!
!        D and E are unaltered (unless W overwrites D).
!
!        Elements of E2, corresponding to elements of E regarded as
!          negligible, have been replaced by zero causing the matrix
!          to split into a direct sum of submatrices.  E2(1) is set
!          to 0.0e0 if the smallest eigenvalues have been found, and
!          to 2.0e0 if the largest eigenvalues have been found.  E2
!          is otherwise unaltered (unless overwritten by BD).
!
!        W contains the M algebraically smallest eigenvalues in
!          ascending order, or the M largest eigenvalues in descending
!          order.  If an error exit is made because of an incorrect
!          specification of IDEF, no eigenvalues are found.  If the
!          Newton iterates for a particular eigenvalue are not monotone,
!          the best estimate obtained is returned and IERR is set.
!          W is a one-dimensional REAL array, dimensioned W(N).  W need
!          not be distinct from D.
!
!        IND contains in its first M positions the submatrix indices
!          associated with the corresponding eigenvalues in W --
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc.
!          IND is an one-dimensional INTEGER array, dimensioned IND(N).
!
!        BD contains refined bounds for the theoretical errors of the
!          corresponding eigenvalues in W.  These bounds are usually
!          within the tolerance specified by EPS1.  BD is a one-
!          dimensional REAL array, dimensioned BD(N).  BD need not be
!          distinct from E2.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          6*N+1      if  IDEF  is set to 1 and  TYPE  to .TRUE.
!                     when the matrix is NOT positive definite, or
!                     if  IDEF  is set to -1 and  TYPE  to .FALSE.
!                     when the matrix is NOT negative definite,
!                     no eigenvalues are computed, or
!                     M is greater than N,
!          5*N+K      if successive iterates to the K-th eigenvalue
!                     are NOT monotone increasing, where K refers
!                     to the last such occurrence.
!
!     Note that subroutine TRIDIB is generally faster and more
!     accurate than RATQR if the eigenvalues are clustered.
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
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RATQR
!
  INTEGER I,J,K,M,N,II,JJ,K1,IDEF,IERR,JDEF
  REAL D(*),E(*),E2(*),W(*),BD(*)
  REAL F,P,Q,R,S,EP,QP,ERR,TOT,EPS1,DELTA,MACHEP
  INTEGER IND(*)
  LOGICAL FIRST, TYPE
!
  SAVE FIRST, MACHEP
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  RATQR
  if (FIRST) THEN
     MACHEP = R1MACH(4)
  end if
  FIRST = .FALSE.
!
  IERR = 0
  JDEF = IDEF
!     .......... COPY D ARRAY INTO W ..........
  DO 20 I = 1, N
   20 W(I) = D(I)
!
  if (TYPE) go to 40
  J = 1
  go to 400
   40 ERR = 0.0E0
  S = 0.0E0
!     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DEFINE
!                INITIAL SHIFT FROM LOWER GERSCHGORIN BOUND.
!                COPY E2 ARRAY INTO BD ..........
  TOT = W(1)
  Q = 0.0E0
  J = 0
!
  DO 100 I = 1, N
     P = Q
     if (I  ==  1) go to 60
     if (P  >  MACHEP * (ABS(D(I)) + ABS(D(I-1)))) go to 80
   60    E2(I) = 0.0E0
   80    BD(I) = E2(I)
!     .......... COUNT ALSO if ELEMENT OF E2 HAS UNDERFLOWED ..........
     if (E2(I)  ==  0.0E0) J = J + 1
     IND(I) = J
     Q = 0.0E0
     if (I  /=  N) Q = ABS(E(I+1))
     TOT = MIN(W(I)-P-Q,TOT)
  100 CONTINUE
!
  if (JDEF  ==  1 .AND. TOT  <  0.0E0) go to 140
!
  DO 110 I = 1, N
  110 W(I) = W(I) - TOT
!
  go to 160
  140 TOT = 0.0E0
!
  160 DO 360 K = 1, M
!     .......... NEXT QR TRANSFORMATION ..........
  180    TOT = TOT + S
     DELTA = W(N) - S
     I = N
     F = ABS(MACHEP*TOT)
     if (EPS1  <  F) EPS1 = F
     if (DELTA  >  EPS1) go to 190
     if (DELTA  <  (-EPS1)) go to 1000
     go to 300
!     .......... REPLACE SMALL SUB-DIAGONAL SQUARES BY ZERO
!                TO REDUCE THE INCIDENCE OF UNDERFLOWS ..........
  190    if (K  ==  N) go to 210
     K1 = K + 1
     DO 200 J = K1, N
        if (BD(J)  <=  (MACHEP*(W(J)+W(J-1))) ** 2) BD(J) = 0.0E0
  200    CONTINUE
!
  210    F = BD(N) / DELTA
     QP = DELTA + F
     P = 1.0E0
     if (K  ==  N) go to 260
     K1 = N - K
!     .......... FOR I=N-1 STEP -1 UNTIL K DO -- ..........
     DO 240 II = 1, K1
        I = N - II
        Q = W(I) - S - F
        R = Q / QP
        P = P * R + 1.0E0
        EP = F * R
        W(I+1) = QP + EP
        DELTA = Q - EP
        if (DELTA  >  EPS1) go to 220
        if (DELTA  <  (-EPS1)) go to 1000
        go to 300
  220       F = BD(I) / Q
        QP = DELTA + F
        BD(I+1) = QP * EP
  240    CONTINUE
!
  260    W(K) = QP
     S = QP / P
     if (TOT + S  >  TOT) go to 180
!     .......... SET ERROR -- IRREGULAR END OF ITERATION.
!                DEFLATE MINIMUM DIAGONAL ELEMENT ..........
     IERR = 5 * N + K
     S = 0.0E0
     DELTA = QP
!
     DO 280 J = K, N
        if (W(J)  >  DELTA) go to 280
        I = J
        DELTA = W(J)
  280    CONTINUE
!     .......... CONVERGENCE ..........
  300    if (I  <  N) BD(I+1) = BD(I) * F / QP
     II = IND(I)
     if (I  ==  K) go to 340
     K1 = I - K
!     .......... FOR J=I-1 STEP -1 UNTIL K DO -- ..........
     DO 320 JJ = 1, K1
        J = I - JJ
        W(J+1) = W(J) - S
        BD(J+1) = BD(J)
        IND(J+1) = IND(J)
  320    CONTINUE
!
  340    W(K) = TOT
     ERR = ERR + ABS(DELTA)
     BD(K) = ERR
     IND(K) = II
  360 CONTINUE
!
  if (TYPE) go to 1001
  F = BD(1)
  E2(1) = 2.0E0
  BD(1) = F
  J = 2
!     .......... NEGATE ELEMENTS OF W FOR LARGEST VALUES ..........
  400 DO 500 I = 1, N
  500 W(I) = -W(I)
!
  JDEF = -JDEF
  go to (40,1001), J
!     .......... SET ERROR -- IDEF SPECIFIED INCORRECTLY ..........
 1000 IERR = 6 * N + 1
 1001 RETURN
end
