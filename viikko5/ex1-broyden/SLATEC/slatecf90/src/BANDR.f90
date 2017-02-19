subroutine BANDR (NM, N, MB, A, D, E, E2, MATZ, Z)
!
!! BANDR reduces a real symmetric band matrix to symmetric tridiagonal ...
!  matrix and, optionally, accumulates orthogonal similarity transformations.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B1
!***TYPE      SINGLE PRECISION (BANDR-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure BANDRD,
!     NUM. MATH. 12, 231-241(1968) by Schwarz.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 273-283(1971).
!
!     This subroutine reduces a REAL SYMMETRIC BAND matrix
!     to a symmetric tridiagonal matrix using and optionally
!     accumulating orthogonal similarity transformations.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        MB is the (half) band width of the matrix, defined as the
!          number of adjacent diagonals, including the principal
!          diagonal, required to specify the non-zero portion of the
!          lower triangle of the matrix.  MB is less than or equal
!          to N.  MB is an INTEGER variable.
!
!        A contains the lower triangle of the real symmetric band
!          matrix.  Its lowest subdiagonal is stored in the last
!          N+1-MB  positions of the first column, its next subdiagonal
!          in the last  N+2-MB  positions of the second column, further
!          subdiagonals similarly, and finally its principal diagonal
!          in the  N  positions of the last column.  Contents of storage
!          locations not part of the matrix are arbitrary.  A is a
!          two-dimensional REAL array, dimensioned A(NM,MB).
!
!        MATZ should be set to .TRUE. if the transformation matrix is
!          to be accumulated, and to .FALSE. otherwise.  MATZ is a
!          LOGICAL variable.
!
!     On OUTPUT
!
!        A has been destroyed, except for its last two columns which
!          contain a copy of the tridiagonal matrix.
!
!        D contains the diagonal elements of the tridiagonal matrix.
!          D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the tridiagonal
!          matrix in its last N-1 positions.  E(1) is set to zero.
!          E is a one-dimensional REAL array, dimensioned E(N).
!
!        E2 contains the squares of the corresponding elements of E.
!          E2 may coincide with E if the squares are not needed.
!          E2 is a one-dimensional REAL array, dimensioned E2(N).
!
!        Z contains the orthogonal transformation matrix produced in
!          the reduction if MATZ has been set to .TRUE.  Otherwise, Z
!          is not referenced.  Z is a two-dimensional REAL array,
!          dimensioned Z(NM,N).
!
!     Questions and comments should be directed to B. S. Garbow,
!     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BANDR
!
  INTEGER J,K,L,N,R,I1,I2,J1,J2,KR,MB,MR,M1,NM,N2,R1,UGL,MAXL,MAXR
  REAL A(NM,*),D(*),E(*),E2(*),Z(NM,*)
  REAL G,U,B1,B2,C2,F1,F2,S2,DMIN,DMINRT
  LOGICAL MATZ
!
!***FIRST EXECUTABLE STATEMENT  BANDR
  DMIN = 2.0E0**(-64)
  DMINRT = 2.0E0**(-32)
!     .......... INITIALIZE DIAGONAL SCALING MATRIX ..........
  DO 30 J = 1, N
   30 D(J) = 1.0E0
!
  if (.NOT. MATZ) go to 60
!
  DO 50 J = 1, N
!
     DO 40 K = 1, N
   40    Z(J,K) = 0.0E0
!
     Z(J,J) = 1.0E0
   50 CONTINUE
!
   60 M1 = MB - 1
  if (M1 - 1) 900, 800, 70
   70 N2 = N - 2
!
  DO 700 K = 1, N2
     MAXR = MIN(M1,N-K)
!     .......... FOR R=MAXR STEP -1 UNTIL 2 DO -- ..........
     DO 600 R1 = 2, MAXR
        R = MAXR + 2 - R1
        KR = K + R
        MR = MB - R
        G = A(KR,MR)
        A(KR-1,1) = A(KR-1,MR+1)
        UGL = K
!
        DO 500 J = KR, N, M1
           J1 = J - 1
           J2 = J1 - 1
           if (G  ==  0.0E0) go to 600
           B1 = A(J1,1) / G
           B2 = B1 * D(J1) / D(J)
           S2 = 1.0E0 / (1.0E0 + B1 * B2)
           if (S2  >=  0.5E0 ) go to 450
           B1 = G / A(J1,1)
           B2 = B1 * D(J) / D(J1)
           C2 = 1.0E0 - S2
           D(J1) = C2 * D(J1)
           D(J) = C2 * D(J)
           F1 = 2.0E0 * A(J,M1)
           F2 = B1 * A(J1,MB)
           A(J,M1) = -B2 * (B1 * A(J,M1) - A(J,MB)) - F2 + A(J,M1)
           A(J1,MB) = B2 * (B2 * A(J,MB) + F1) + A(J1,MB)
           A(J,MB) = B1 * (F2 - F1) + A(J,MB)
!
           DO 200 L = UGL, J2
              I2 = MB - J + L
              U = A(J1,I2+1) + B2 * A(J,I2)
              A(J,I2) = -B1 * A(J1,I2+1) + A(J,I2)
              A(J1,I2+1) = U
  200          CONTINUE
!
           UGL = J
           A(J1,1) = A(J1,1) + B2 * G
           if (J  ==  N) go to 350
           MAXL = MIN(M1,N-J1)
!
           DO 300 L = 2, MAXL
              I1 = J1 + L
              I2 = MB - L
              U = A(I1,I2) + B2 * A(I1,I2+1)
              A(I1,I2+1) = -B1 * A(I1,I2) + A(I1,I2+1)
              A(I1,I2) = U
  300          CONTINUE
!
           I1 = J + M1
           if (I1  >  N) go to 350
           G = B2 * A(I1,1)
  350          if (.NOT. MATZ) go to 500
!
           DO 400 L = 1, N
              U = Z(L,J1) + B2 * Z(L,J)
              Z(L,J) = -B1 * Z(L,J1) + Z(L,J)
              Z(L,J1) = U
  400          CONTINUE
!
           go to 500
!
  450          U = D(J1)
           D(J1) = S2 * D(J)
           D(J) = S2 * U
           F1 = 2.0E0 * A(J,M1)
           F2 = B1 * A(J,MB)
           U = B1 * (F2 - F1) + A(J1,MB)
           A(J,M1) = B2 * (B1 * A(J,M1) - A(J1,MB)) + F2 - A(J,M1)
           A(J1,MB) = B2 * (B2 * A(J1,MB) + F1) + A(J,MB)
           A(J,MB) = U
!
           DO 460 L = UGL, J2
              I2 = MB - J + L
              U = B2 * A(J1,I2+1) + A(J,I2)
              A(J,I2) = -A(J1,I2+1) + B1 * A(J,I2)
              A(J1,I2+1) = U
  460          CONTINUE
!
           UGL = J
           A(J1,1) = B2 * A(J1,1) + G
           if (J  ==  N) go to 480
           MAXL = MIN(M1,N-J1)
!
           DO 470 L = 2, MAXL
              I1 = J1 + L
              I2 = MB - L
              U = B2 * A(I1,I2) + A(I1,I2+1)
              A(I1,I2+1) = -A(I1,I2) + B1 * A(I1,I2+1)
              A(I1,I2) = U
  470          CONTINUE
!
           I1 = J + M1
           if (I1  >  N) go to 480
           G = A(I1,1)
           A(I1,1) = B1 * A(I1,1)
  480          if (.NOT. MATZ) go to 500
!
           DO 490 L = 1, N
              U = B2 * Z(L,J1) + Z(L,J)
              Z(L,J) = -Z(L,J1) + B1 * Z(L,J)
              Z(L,J1) = U
  490          CONTINUE
!
  500       CONTINUE
!
  600    CONTINUE
!
     if (MOD(K,64)  /=  0) go to 700
!     .......... RESCALE TO AVOID UNDERFLOW OR OVERFLOW ..........
     DO 650 J = K, N
        if (D(J)  >=  DMIN) go to 650
        MAXL = MAX(1,MB+1-J)
!
        DO 610 L = MAXL, M1
  610       A(J,L) = DMINRT * A(J,L)
!
        if (J  ==  N) go to 630
        MAXL = MIN(M1,N-J)
!
        DO 620 L = 1, MAXL
           I1 = J + L
           I2 = MB - L
           A(I1,I2) = DMINRT * A(I1,I2)
  620       CONTINUE
!
  630       if (.NOT. MATZ) go to 645
!
        DO 640 L = 1, N
  640       Z(L,J) = DMINRT * Z(L,J)
!
  645       A(J,MB) = DMIN * A(J,MB)
        D(J) = D(J) / DMIN
  650    CONTINUE
!
  700 CONTINUE
!     .......... FORM SQUARE ROOT OF SCALING MATRIX ..........
  800 DO 810 J = 2, N
  810 E(J) = SQRT(D(J))
!
  if (.NOT. MATZ) go to 840
!
  DO 830 J = 1, N
!
     DO 820 K = 2, N
  820    Z(J,K) = E(K) * Z(J,K)
!
  830 CONTINUE
!
  840 U = 1.0E0
!
  DO 850 J = 2, N
     A(J,M1) = U * E(J) * A(J,M1)
     U = E(J)
     E2(J) = A(J,M1) ** 2
     A(J,MB) = D(J) * A(J,MB)
     D(J) = A(J,MB)
     E(J) = A(J,M1)
  850 CONTINUE
!
  D(1) = A(1,MB)
  E(1) = 0.0E0
  E2(1) = 0.0E0
  go to 1001
!
  900 DO 950 J = 1, N
     D(J) = A(J,MB)
     E(J) = 0.0E0
     E2(J) = 0.0E0
  950 CONTINUE
!
 1001 RETURN
end
