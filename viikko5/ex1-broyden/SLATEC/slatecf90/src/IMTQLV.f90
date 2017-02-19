subroutine IMTQLV (N, D, E, E2, W, IND, IERR, RV1)
!
!! IMTQLV computes the eigenvalues of a symmetric tridiagonal matrix ...
!            using the implicit QL method.  Eigenvectors may be computed ...
!            later.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5, D4C2A
!***TYPE      SINGLE PRECISION (IMTQLV-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a variant of  IMTQL1  which is a translation of
!     ALGOL procedure IMTQL1, NUM. MATH. 12, 377-383(1968) by Martin and
!     Wilkinson, as modified in NUM. MATH. 15, 450(1970) by Dubrulle.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
!
!     This subroutine finds the eigenvalues of a SYMMETRIC TRIDIAGONAL
!     matrix by the implicit QL method and associates with them
!     their corresponding submatrix indices.
!
!     On INPUT
!
!        N is the order of the matrix.  N is an INTEGER variable.
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
!     On OUTPUT
!
!        D and E are unaltered.
!
!        Elements of E2, corresponding to elements of E regarded as
!          negligible, have been replaced by zero causing the matrix to
!          split into a direct sum of submatrices.  E2(1) is also set
!          to zero.
!
!        W contains the eigenvalues in ascending order.  If an error
!          exit is made, the eigenvalues are correct and ordered for
!          indices 1, 2, ..., IERR-1, but may not be the smallest
!          eigenvalues.  W is a one-dimensional REAL array, dimensioned
!          W(N).
!
!        IND contains the submatrix indices associated with the
!          corresponding eigenvalues in W -- 1 for eigenvalues belonging
!          to the first submatrix from the top, 2 for those belonging to
!          the second submatrix, etc.  IND is a one-dimensional REAL
!          array, dimensioned IND(N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!                     The eigenvalues should be correct for indices
!                     1, 2, ..., IERR-1.  These eigenvalues are
!                     ordered, but are not necessarily the smallest.
!
!        RV1 is a one-dimensional REAL array used for temporary storage,
!          dimensioned RV1(N).
!
!     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  IMTQLV
!
  INTEGER I,J,K,L,M,N,II,MML,TAG,IERR
  REAL D(*),E(*),E2(*),W(*),RV1(*)
  REAL B,C,F,G,P,R,S,S1,S2
  REAL PYTHAG
  INTEGER IND(*)
!
!***FIRST EXECUTABLE STATEMENT  IMTQLV
  IERR = 0
  K = 0
  TAG = 0
!
  DO 100 I = 1, N
     W(I) = D(I)
     if (I  /=  1) RV1(I-1) = E(I)
  100 CONTINUE
!
  E2(1) = 0.0E0
  RV1(N) = 0.0E0
!
  DO 290 L = 1, N
     J = 0
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
        if (M  ==  N) go to 120
        S1 = ABS(W(M)) + ABS(W(M+1))
        S2 = S1 + ABS(RV1(M))
        if (S2  ==  S1) go to 120
!     .......... GUARD AGAINST UNDERFLOWED ELEMENT OF E2 ..........
        if (E2(M+1)  ==  0.0E0) go to 125
  110    CONTINUE
!
  120    if (M  <=  K) go to 130
     if (M  /=  N) E2(M+1) = 0.0E0
  125    K = M
     TAG = TAG + 1
  130    P = W(L)
     if (M  ==  L) go to 215
     if (J  ==  30) go to 1000
     J = J + 1
!     .......... FORM SHIFT ..........
     G = (W(L+1) - P) / (2.0E0 * RV1(L))
     R = PYTHAG(G,1.0E0)
     G = W(M) - P + RV1(L) / (G + SIGN(R,G))
     S = 1.0E0
     C = 1.0E0
     P = 0.0E0
     MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
     DO 200 II = 1, MML
        I = M - II
        F = S * RV1(I)
        B = C * RV1(I)
        if (ABS(F)  <  ABS(G)) go to 150
        C = G / F
        R = SQRT(C*C+1.0E0)
        RV1(I+1) = F * R
        S = 1.0E0 / R
        C = C * S
        go to 160
  150       S = F / G
        R = SQRT(S*S+1.0E0)
        RV1(I+1) = G * R
        C = 1.0E0 / R
        S = S * C
  160       G = W(I+1) - P
        R = (W(I) - G) * S + 2.0E0 * C * B
        P = S * R
        W(I+1) = G + P
        G = C * R - B
  200    CONTINUE
!
     W(L) = W(L) - P
     RV1(L) = G
     RV1(M) = 0.0E0
     go to 105
!     .......... ORDER EIGENVALUES ..........
  215    if (L  ==  1) go to 250
!     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
     DO 230 II = 2, L
        I = L + 2 - II
        if (P  >=  W(I-1)) go to 270
        W(I) = W(I-1)
        IND(I) = IND(I-1)
  230    CONTINUE
!
  250    I = 1
  270    W(I) = P
     IND(I) = TAG
  290 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
end
