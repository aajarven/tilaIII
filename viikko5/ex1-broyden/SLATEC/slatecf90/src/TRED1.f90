subroutine TRED1 (NM, N, A, D, E, E2)
!
!! TRED1 reduces a real symmetric matrix to symmetric tridiagonal ...
!  matrix using orthogonal similarity transformations.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B1
!***TYPE      SINGLE PRECISION (TRED1-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TRED1,
!     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     This subroutine reduces a REAL SYMMETRIC matrix
!     to a symmetric tridiagonal matrix using
!     orthogonal similarity transformations.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, A, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        A contains the real symmetric input matrix.  Only the lower
!          triangle of the matrix need be supplied.  A is a two-
!          dimensional REAL array, dimensioned A(NM,N).
!
!     On Output
!
!        A contains information about the orthogonal transformations
!          used in the reduction in its strict lower triangle.  The
!          full upper triangle of A is unaltered.
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
!***END PROLOGUE  TRED1
!
  INTEGER I,J,K,L,N,II,NM,JP1
  REAL A(NM,*),D(*),E(*),E2(*)
  REAL F,G,H,SCALE
!
!***FIRST EXECUTABLE STATEMENT  TRED1
  DO 100 I = 1, N
  100 D(I) = A(I,I)
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
  DO 300 II = 1, N
     I = N + 1 - II
     L = I - 1
     H = 0.0E0
     SCALE = 0.0E0
     if (L  <  1) go to 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
     DO 120 K = 1, L
  120    SCALE = SCALE + ABS(A(I,K))
!
     if (SCALE  /=  0.0E0) go to 140
  130    E(I) = 0.0E0
     E2(I) = 0.0E0
     go to 290
!
  140    DO 150 K = 1, L
        A(I,K) = A(I,K) / SCALE
        H = H + A(I,K) * A(I,K)
  150    CONTINUE
!
     E2(I) = SCALE * SCALE * H
     F = A(I,L)
     G = -SIGN(SQRT(H),F)
     E(I) = SCALE * G
     H = H - F * G
     A(I,L) = F - G
     if (L  ==  1) go to 270
     F = 0.0E0
!
     DO 240 J = 1, L
        G = 0.0E0
!     .......... FORM ELEMENT OF A*U ..........
        DO 180 K = 1, J
  180       G = G + A(J,K) * A(I,K)
!
        JP1 = J + 1
        if (L  <  JP1) go to 220
!
        DO 200 K = JP1, L
  200       G = G + A(K,J) * A(I,K)
!     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
        F = F + E(J) * A(I,J)
  240    CONTINUE
!
     H = F / (H + H)
!     .......... FORM REDUCED A ..........
     DO 260 J = 1, L
        F = A(I,J)
        G = E(J) - H * F
        E(J) = G
!
        DO 260 K = 1, J
           A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
  260    CONTINUE
!
  270    DO 280 K = 1, L
  280    A(I,K) = SCALE * A(I,K)
!
  290    H = D(I)
     D(I) = A(I,I)
     A(I,I) = H
  300 CONTINUE
!
  return
end
