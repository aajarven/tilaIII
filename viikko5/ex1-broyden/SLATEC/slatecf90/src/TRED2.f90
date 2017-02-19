subroutine TRED2 (NM, N, A, D, E, Z)
!
!! TRED2 reduces a real symmetric matrix to a symmetric tridiagonal ...
!  matrix using and accumulating orthogonal transformations.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B1
!***TYPE      SINGLE PRECISION (TRED2-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TRED2,
!     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     This subroutine reduces a REAL SYMMETRIC matrix to a
!     symmetric tridiagonal matrix using and accumulating
!     orthogonal similarity transformations.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
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
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the symmetric
!          tridiagonal matrix in its last N-1 positions.  E(1) is set
!          to zero.  E is a one-dimensional REAL array, dimensioned
!          E(N).
!
!        Z contains the orthogonal transformation matrix produced in
!          the reduction.  Z is a two-dimensional REAL array,
!          dimensioned Z(NM,N).
!
!        A and Z may coincide.  If distinct, A is unaltered.
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
!***END PROLOGUE  TRED2
!
  INTEGER I,J,K,L,N,II,NM,JP1
  REAL A(NM,*),D(*),E(*),Z(NM,*)
  REAL F,G,H,HH,SCALE
!
!***FIRST EXECUTABLE STATEMENT  TRED2
  DO 100 I = 1, N
!
     DO 100 J = 1, I
        Z(I,J) = A(I,J)
  100 CONTINUE
!
  if (N  ==  1) go to 320
!     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
  DO 300 II = 2, N
     I = N + 2 - II
     L = I - 1
     H = 0.0E0
     SCALE = 0.0E0
     if (L  <  2) go to 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
     DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))
!
     if (SCALE  /=  0.0E0) go to 140
  130    E(I) = Z(I,L)
     go to 290
!
  140    DO 150 K = 1, L
        Z(I,K) = Z(I,K) / SCALE
        H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
!
     F = Z(I,L)
     G = -SIGN(SQRT(H),F)
     E(I) = SCALE * G
     H = H - F * G
     Z(I,L) = F - G
     F = 0.0E0
!
     DO 240 J = 1, L
        Z(J,I) = Z(I,J) / H
        G = 0.0E0
!     .......... FORM ELEMENT OF A*U ..........
        DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
!
        JP1 = J + 1
        if (L  <  JP1) go to 220
!
        DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
!     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
        F = F + E(J) * Z(I,J)
  240    CONTINUE
!
     HH = F / (H + H)
!     .......... FORM REDUCED A ..........
     DO 260 J = 1, L
        F = Z(I,J)
        G = E(J) - HH * F
        E(J) = G
!
        DO 260 K = 1, J
           Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
!
  290    D(I) = H
  300 CONTINUE
!
  320 D(1) = 0.0E0
  E(1) = 0.0E0
!     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
  DO 500 I = 1, N
     L = I - 1
     if (D(I)  ==  0.0E0) go to 380
!
     DO 360 J = 1, L
        G = 0.0E0
!
        DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
!
        DO 360 K = 1, L
           Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
!
  380    D(I) = Z(I,I)
     Z(I,I) = 1.0E0
     if (L  <  1) go to 500
!
     DO 400 J = 1, L
        Z(I,J) = 0.0E0
        Z(J,I) = 0.0E0
  400    CONTINUE
!
  500 CONTINUE
!
  return
end
