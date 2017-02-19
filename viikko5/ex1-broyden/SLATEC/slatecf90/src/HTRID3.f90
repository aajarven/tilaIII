subroutine HTRID3 (NM, N, A, D, E, E2, TAU)
!
!! HTRID3 reduces a complex Hermitian (packed) matrix to a real ...
!            symmetric tridiagonal matrix by unitary similarity
!            transformations.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B1
!***TYPE      SINGLE PRECISION (HTRID3-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of a complex analogue of
!     the ALGOL procedure TRED3, NUM. MATH. 11, 181-195(1968)
!     by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     This subroutine reduces a COMPLEX HERMITIAN matrix, stored as
!     a single square array, to a real symmetric tridiagonal matrix
!     using unitary similarity transformations.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, A, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        A contains the lower triangle of the complex Hermitian input
!          matrix.  The real parts of the matrix elements are stored
!          in the full lower triangle of A, and the imaginary parts
!          are stored in the transposed positions of the strict upper
!          triangle of A.  No storage is required for the zero
!          imaginary parts of the diagonal elements.  A is a two-
!          dimensional REAL array, dimensioned A(NM,N).
!
!     On OUTPUT
!
!        A contains some information about the unitary transformations
!          used in the reduction.
!
!        D contains the diagonal elements of the real symmetric
!          tridiagonal matrix.  D is a one-dimensional REAL array,
!          dimensioned D(N).
!
!        E contains the subdiagonal elements of the real tridiagonal
!          matrix in its last N-1 positions.  E(1) is set to zero.
!          E is a one-dimensional REAL array, dimensioned E(N).
!
!        E2 contains the squares of the corresponding elements of E.
!          E2(1) is set to zero.  E2 may coincide with E if the squares
!          are not needed.  E2 is a one-dimensional REAL array,
!          dimensioned E2(N).
!
!        TAU contains further information about the transformations.
!          TAU is a one-dimensional REAL array, dimensioned TAU(2,N).
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
!***END PROLOGUE  HTRID3
!
  INTEGER I,J,K,L,N,II,NM,JM1,JP1
  REAL A(NM,*),D(*),E(*),E2(*),TAU(2,*)
  REAL F,G,H,FI,GI,HH,SI,SCALE
  REAL PYTHAG
!
!***FIRST EXECUTABLE STATEMENT  HTRID3
  TAU(1,N) = 1.0E0
  TAU(2,N) = 0.0E0
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
  DO 300 II = 1, N
     I = N + 1 - II
     L = I - 1
     H = 0.0E0
     SCALE = 0.0E0
     if (L  <  1) go to 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
     DO 120 K = 1, L
  120    SCALE = SCALE + ABS(A(I,K)) + ABS(A(K,I))
!
     if (SCALE  /=  0.0E0) go to 140
     TAU(1,L) = 1.0E0
     TAU(2,L) = 0.0E0
  130    E(I) = 0.0E0
     E2(I) = 0.0E0
     go to 290
!
  140    DO 150 K = 1, L
        A(I,K) = A(I,K) / SCALE
        A(K,I) = A(K,I) / SCALE
        H = H + A(I,K) * A(I,K) + A(K,I) * A(K,I)
  150    CONTINUE
!
     E2(I) = SCALE * SCALE * H
     G = SQRT(H)
     E(I) = SCALE * G
     F = PYTHAG(A(I,L),A(L,I))
!     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T ..........
     if (F  ==  0.0E0) go to 160
     TAU(1,L) = (A(L,I) * TAU(2,I) - A(I,L) * TAU(1,I)) / F
     SI = (A(I,L) * TAU(2,I) + A(L,I) * TAU(1,I)) / F
     H = H + F * G
     G = 1.0E0 + G / F
     A(I,L) = G * A(I,L)
     A(L,I) = G * A(L,I)
     if (L  ==  1) go to 270
     go to 170
  160    TAU(1,L) = -TAU(1,I)
     SI = TAU(2,I)
     A(I,L) = G
  170    F = 0.0E0
!
     DO 240 J = 1, L
        G = 0.0E0
        GI = 0.0E0
        if (J  ==  1) go to 190
        JM1 = J - 1
!     .......... FORM ELEMENT OF A*U ..........
        DO 180 K = 1, JM1
           G = G + A(J,K) * A(I,K) + A(K,J) * A(K,I)
           GI = GI - A(J,K) * A(K,I) + A(K,J) * A(I,K)
  180       CONTINUE
!
  190       G = G + A(J,J) * A(I,J)
        GI = GI - A(J,J) * A(J,I)
        JP1 = J + 1
        if (L  <  JP1) go to 220
!
        DO 200 K = JP1, L
           G = G + A(K,J) * A(I,K) - A(J,K) * A(K,I)
           GI = GI - A(K,J) * A(K,I) - A(J,K) * A(I,K)
  200       CONTINUE
!     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
        TAU(2,J) = GI / H
        F = F + E(J) * A(I,J) - TAU(2,J) * A(J,I)
  240    CONTINUE
!
     HH = F / (H + H)
!     .......... FORM REDUCED A ..........
     DO 260 J = 1, L
        F = A(I,J)
        G = E(J) - HH * F
        E(J) = G
        FI = -A(J,I)
        GI = TAU(2,J) - HH * FI
        TAU(2,J) = -GI
        A(J,J) = A(J,J) - 2.0E0 * (F * G + FI * GI)
        if (J  ==  1) go to 260
        JM1 = J - 1
!
        DO 250 K = 1, JM1
           A(J,K) = A(J,K) - F * E(K) - G * A(I,K) &
                           + FI * TAU(2,K) + GI * A(K,I)
           A(K,J) = A(K,J) - F * TAU(2,K) - G * A(K,I) &
                           - FI * E(K) - GI * A(I,K)
  250       CONTINUE
!
  260    CONTINUE
!
  270    DO 280 K = 1, L
        A(I,K) = SCALE * A(I,K)
        A(K,I) = SCALE * A(K,I)
  280    CONTINUE
!
     TAU(2,L) = -SI
  290    D(I) = A(I,I)
     A(I,I) = SCALE * SQRT(H)
  300 CONTINUE
!
  return
end
