subroutine TRED3 (N, NV, A, D, E, E2)
!
!! TRED3 reduces a real symmetric matrix stored in packed form to
!  symmetric tridiagonal matrix using orthogonal transformations.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B1
!***TYPE      SINGLE PRECISION (TRED3-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TRED3,
!     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     This subroutine reduces a REAL SYMMETRIC matrix, stored as
!     a one-dimensional array, to a symmetric tridiagonal matrix
!     using orthogonal similarity transformations.
!
!     On Input
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!
!        NV is an INTEGER variable set equal to the dimension of the
!          array A as specified in the calling program.  NV must not
!          be less than  N*(N+1)/2.
!
!        A contains the lower triangle, stored row-wise, of the real
!          symmetric packed matrix.  A is a one-dimensional REAL
!          array, dimensioned A(NV).
!
!     On Output
!
!        A contains information about the orthogonal transformations
!          used in the reduction in its first N*(N+1)/2 positions.
!
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the symmetric
!          tridiagonal matrix in its last N-1 positions.  E(1) is set
!          to zero.  E is a one-dimensional REAL array, dimensioned
!          E(N).
!
!        E2 contains the squares of the corresponding elements of E.
!          E2 may coincide with E if the squares are not needed.
!          E2 is a one-dimensional REAL array, dimensioned E2(N).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  TRED3
!
  INTEGER I,J,K,L,N,II,IZ,JK,NV
  REAL A(*),D(*),E(*),E2(*)
  REAL F,G,H,HH,SCALE
!
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
!***FIRST EXECUTABLE STATEMENT  TRED3
  DO  300 II = 1, N
     I = N + 1 - II
     L = I - 1
     IZ = (I * L) / 2
     H = 0.0E0
     SCALE = 0.0E0
     if (L  <  1) go to 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
     DO 120 K = 1, L
        IZ = IZ + 1
        D(K) = A(IZ)
        SCALE = SCALE + ABS(D(K))
  120    CONTINUE
!
     if (SCALE  /=  0.0E0) go to 140
  130    E(I) = 0.0E0
     E2(I) = 0.0E0
     go to 290
!
  140    DO 150 K = 1, L
        D(K) = D(K) / SCALE
        H = H + D(K) * D(K)
  150    CONTINUE
!
     E2(I) = SCALE * SCALE * H
     F = D(L)
     G = -SIGN(SQRT(H),F)
     E(I) = SCALE * G
     H = H - F * G
     D(L) = F - G
     A(IZ) = SCALE * D(L)
     if (L  ==  1) go to 290
     F = 0.0E0
!
     DO 240 J = 1, L
        G = 0.0E0
        JK = (J * (J-1)) / 2
!     .......... FORM ELEMENT OF A*U ..........
        DO 180 K = 1, L
           JK = JK + 1
           if (K  >  J) JK = JK + K - 2
           G = G + A(JK) * D(K)
  180       CONTINUE
!     .......... FORM ELEMENT OF P ..........
        E(J) = G / H
        F = F + E(J) * D(J)
  240    CONTINUE
!
     HH = F / (H + H)
     JK = 0
!     .......... FORM REDUCED A ..........
     DO 260 J = 1, L
        F = D(J)
        G = E(J) - HH * F
        E(J) = G
!
        DO 260 K = 1, J
           JK = JK + 1
           A(JK) = A(JK) - F * E(K) - G * D(K)
  260    CONTINUE
!
  290    D(I) = A(IZ+1)
     A(IZ+1) = SCALE * SQRT(H)
  300 CONTINUE
!
  return
end
