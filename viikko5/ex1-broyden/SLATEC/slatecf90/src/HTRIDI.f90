subroutine HTRIDI (NM, N, AR, AI, D, E, E2, TAU)
!
!! HTRIDI reduces a complex Hermitian matrix to a real symmetric ...
!            tridiagonal matrix using unitary similarity
!            transformations.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B1
!***TYPE      SINGLE PRECISION (HTRIDI-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of a complex analogue of
!     the ALGOL procedure TRED1, NUM. MATH. 11, 181-195(1968)
!     by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     This subroutine reduces a COMPLEX HERMITIAN matrix
!     to a real symmetric tridiagonal matrix using
!     unitary similarity transformations.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, AR and AI, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A=(AR,AI).  N is an INTEGER
!          variable. N must be less than or equal to NM.
!
!        AR and AI contain the real and imaginary parts, respectively,
!          of the complex Hermitian input matrix.  Only the lower
!          triangle of the matrix need be supplied.  AR and AI are two-
!          dimensional REAL arrays, dimensioned AR(NM,N) and AI(NM,N).
!
!     On OUTPUT
!
!        AR and AI contain some information about the unitary trans-
!          formations used in the reduction in the strict lower triangle
!          of AR and the full lower triangle of AI.  The rest of the
!          matrices are unaltered.
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
!***END PROLOGUE  HTRIDI
!
  INTEGER I,J,K,L,N,II,NM,JP1
  REAL AR(NM,*),AI(NM,*),D(*),E(*),E2(*),TAU(2,*)
  REAL F,G,H,FI,GI,HH,SI,SCALE
  REAL PYTHAG
!
!***FIRST EXECUTABLE STATEMENT  HTRIDI
  TAU(1,N) = 1.0E0
  TAU(2,N) = 0.0E0
!
  DO 100 I = 1, N
  100 D(I) = AR(I,I)
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
  DO 300 II = 1, N
     I = N + 1 - II
     L = I - 1
     H = 0.0E0
     SCALE = 0.0E0
     if (L  <  1) go to 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
     DO 120 K = 1, L
  120    SCALE = SCALE + ABS(AR(I,K)) + ABS(AI(I,K))
!
     if (SCALE  /=  0.0E0) go to 140
     TAU(1,L) = 1.0E0
     TAU(2,L) = 0.0E0
  130    E(I) = 0.0E0
     E2(I) = 0.0E0
     go to 290
!
  140    DO 150 K = 1, L
        AR(I,K) = AR(I,K) / SCALE
        AI(I,K) = AI(I,K) / SCALE
        H = H + AR(I,K) * AR(I,K) + AI(I,K) * AI(I,K)
  150    CONTINUE
!
     E2(I) = SCALE * SCALE * H
     G = SQRT(H)
     E(I) = SCALE * G
     F = PYTHAG(AR(I,L),AI(I,L))
!     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T ..........
     if (F  ==  0.0E0) go to 160
     TAU(1,L) = (AI(I,L) * TAU(2,I) - AR(I,L) * TAU(1,I)) / F
     SI = (AR(I,L) * TAU(2,I) + AI(I,L) * TAU(1,I)) / F
     H = H + F * G
     G = 1.0E0 + G / F
     AR(I,L) = G * AR(I,L)
     AI(I,L) = G * AI(I,L)
     if (L  ==  1) go to 270
     go to 170
  160    TAU(1,L) = -TAU(1,I)
     SI = TAU(2,I)
     AR(I,L) = G
  170    F = 0.0E0
!
     DO 240 J = 1, L
        G = 0.0E0
        GI = 0.0E0
!     .......... FORM ELEMENT OF A*U ..........
        DO 180 K = 1, J
           G = G + AR(J,K) * AR(I,K) + AI(J,K) * AI(I,K)
           GI = GI - AR(J,K) * AI(I,K) + AI(J,K) * AR(I,K)
  180       CONTINUE
!
        JP1 = J + 1
        if (L  <  JP1) go to 220
!
        DO 200 K = JP1, L
           G = G + AR(K,J) * AR(I,K) - AI(K,J) * AI(I,K)
           GI = GI - AR(K,J) * AI(I,K) - AI(K,J) * AR(I,K)
  200       CONTINUE
!     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
        TAU(2,J) = GI / H
        F = F + E(J) * AR(I,J) - TAU(2,J) * AI(I,J)
  240    CONTINUE
!
     HH = F / (H + H)
!     .......... FORM REDUCED A ..........
     DO 260 J = 1, L
        F = AR(I,J)
        G = E(J) - HH * F
        E(J) = G
        FI = -AI(I,J)
        GI = TAU(2,J) - HH * FI
        TAU(2,J) = -GI
!
        DO 260 K = 1, J
           AR(J,K) = AR(J,K) - F * E(K) - G * AR(I,K) &
                             + FI * TAU(2,K) + GI * AI(I,K)
           AI(J,K) = AI(J,K) - F * TAU(2,K) - G * AI(I,K) &
                             - FI * E(K) - GI * AR(I,K)
  260    CONTINUE
!
  270    DO 280 K = 1, L
        AR(I,K) = SCALE * AR(I,K)
        AI(I,K) = SCALE * AI(I,K)
  280    CONTINUE
!
     TAU(2,L) = -SI
  290    HH = D(I)
     D(I) = AR(I,I)
     AR(I,I) = HH
     AI(I,I) = SCALE * SQRT(H)
  300 CONTINUE
!
  return
end
