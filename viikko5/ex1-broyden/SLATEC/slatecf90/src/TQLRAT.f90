subroutine TQLRAT (N, D, E2, IERR)
!
!! TQLRAT computes the eigenvalues of symmetric tridiagonal matrix ...
!            using a rational variant of the QL method.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5, D4C2A
!***TYPE      SINGLE PRECISION (TQLRAT-S)
!***KEYWORDS  EIGENVALUES OF A SYMMETRIC TRIDIAGONAL MATRIX, EISPACK,
!             QL METHOD
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TQLRAT.
!
!     This subroutine finds the eigenvalues of a SYMMETRIC
!     TRIDIAGONAL matrix by the rational QL method.
!
!     On Input
!
!        N is the order of the matrix.  N is an INTEGER variable.
!
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E2 contains the squares of the subdiagonal elements of the
!          symmetric tridiagonal matrix in its last N-1 positions.
!          E2(1) is arbitrary.  E2 is a one-dimensional REAL array,
!          dimensioned E2(N).
!
!      On Output
!
!        D contains the eigenvalues in ascending order.  If an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1, 2, ..., IERR-1, but may not be
!          the smallest eigenvalues.
!
!        E2 has been destroyed.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
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
!               C. H. Reinsch, Eigenvalues of a real, symmetric, tri-
!                 diagonal matrix, Algorithm 464, Communications of the
!                 ACM 16, 11 (November 1973), pp. 689.
!***ROUTINES CALLED  PYTHAG, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  TQLRAT
!
  INTEGER I,J,L,M,N,II,L1,MML,IERR
  REAL D(*),E2(*)
  REAL B,C,F,G,H,P,R,S,MACHEP
  REAL PYTHAG
  LOGICAL FIRST
!
  SAVE FIRST, MACHEP
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  TQLRAT
  if (FIRST) THEN
     MACHEP = R1MACH(4)
  end if
  FIRST = .FALSE.
!
  IERR = 0
  if (N  ==  1) go to 1001
!
  DO 100 I = 2, N
  100 E2(I-1) = E2(I)
!
  F = 0.0E0
  B = 0.0E0
  E2(N) = 0.0E0
!
  DO 290 L = 1, N
     J = 0
     H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))
     if (B  >  H) go to 105
     B = H
     C = B * B
!     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
        if (E2(M)  <=  C) go to 120
!     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!
  120    if (M  ==  L) go to 210
  130    if (J  ==  30) go to 1000
     J = J + 1
!     .......... FORM SHIFT ..........
     L1 = L + 1
     S = SQRT(E2(L))
     G = D(L)
     P = (D(L1) - G) / (2.0E0 * S)
     R = PYTHAG(P,1.0E0)
     D(L) = S / (P + SIGN(R,P))
     H = G - D(L)
!
     DO 140 I = L1, N
  140    D(I) = D(I) - H
!
     F = F + H
!     .......... RATIONAL QL TRANSFORMATION ..........
     G = D(M)
     if (G  ==  0.0E0) G = B
     H = G
     S = 0.0E0
     MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
     DO 200 II = 1, MML
        I = M - II
        P = G * H
        R = P + E2(I)
        E2(I+1) = S * R
        S = E2(I) / R
        D(I+1) = H + S * (H + D(I))
        G = D(I) - E2(I) / G
        if (G  ==  0.0E0) G = B
        H = G * P / R
  200    CONTINUE
!
     E2(L) = S * G
     D(L) = H
!     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
     if (H  ==  0.0E0) go to 210
     if (ABS(E2(L))  <=  ABS(C/H)) go to 210
     E2(L) = H * E2(L)
     if (E2(L)  /=  0.0E0) go to 130
  210    P = D(L) + F
!     .......... ORDER EIGENVALUES ..........
     if (L  ==  1) go to 250
!     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
     DO 230 II = 2, L
        I = L + 2 - II
        if (P  >=  D(I-1)) go to 270
        D(I) = D(I-1)
  230    CONTINUE
!
  250    I = 1
  270    D(I) = P
  290 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
end
