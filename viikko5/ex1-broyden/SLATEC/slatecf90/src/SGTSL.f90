subroutine SGTSL (N, C, D, E, B, INFO)
!
!! SGTSL solves a tridiagonal linear system.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2A2A
!***TYPE      SINGLE PRECISION (SGTSL-S, DGTSL-D, CGTSL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE, TRIDIAGONAL
!***AUTHOR  Dongarra, J., (ANL)
!***DESCRIPTION
!
!     SGTSL given a general tridiagonal matrix and a right hand
!     side will find the solution.
!
!     On Entry
!
!        N       INTEGER
!                is the order of the tridiagonal matrix.
!
!        C       REAL(N)
!                is the subdiagonal of the tridiagonal matrix.
!                C(2) through C(N) should contain the subdiagonal.
!                On output, C is destroyed.
!
!        D       REAL(N)
!                is the diagonal of the tridiagonal matrix.
!                On output, D is destroyed.
!
!        E       REAL(N)
!                is the superdiagonal of the tridiagonal matrix.
!                E(1) through E(N-1) should contain the superdiagonal.
!                On output, E is destroyed.
!
!        B       REAL(N)
!                is the right hand side vector.
!
!     On Return
!
!        B       is the solution vector.
!
!        INFO    INTEGER
!                = 0 normal value.
!                = K if the K-th element of the diagonal becomes
!                    exactly zero.  The subroutine returns when
!                    this is detected.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SGTSL
  INTEGER N,INFO
  REAL C(*),D(*),E(*),B(*)
!
  INTEGER K,KB,KP1,NM1,NM2
  REAL T
!***FIRST EXECUTABLE STATEMENT  SGTSL
     INFO = 0
     C(1) = D(1)
     NM1 = N - 1
     if (NM1  <  1) go to 40
        D(1) = E(1)
        E(1) = 0.0E0
        E(N) = 0.0E0
!
        DO 30 K = 1, NM1
           KP1 = K + 1
!
!              FIND THE LARGEST OF THE TWO ROWS
!
           if (ABS(C(KP1))  <  ABS(C(K))) go to 10
!
!                 INTERCHANGE ROW
!
              T = C(KP1)
              C(KP1) = C(K)
              C(K) = T
              T = D(KP1)
              D(KP1) = D(K)
              D(K) = T
              T = E(KP1)
              E(KP1) = E(K)
              E(K) = T
              T = B(KP1)
              B(KP1) = B(K)
              B(K) = T
   10          CONTINUE
!
!              ZERO ELEMENTS
!
           if (C(K)  /=  0.0E0) go to 20
              INFO = K
              go to 100
   20          CONTINUE
           T = -C(KP1)/C(K)
           C(KP1) = D(KP1) + T*D(K)
           D(KP1) = E(KP1) + T*E(K)
           E(KP1) = 0.0E0
           B(KP1) = B(KP1) + T*B(K)
   30       CONTINUE
   40    CONTINUE
     if (C(N)  /=  0.0E0) go to 50
        INFO = N
     go to 90
   50    CONTINUE
!
!           BACK SOLVE
!
        NM2 = N - 2
        B(N) = B(N)/C(N)
        if (N  ==  1) go to 80
           B(NM1) = (B(NM1) - D(NM1)*B(N))/C(NM1)
           if (NM2  <  1) go to 70
           DO 60 KB = 1, NM2
              K = NM2 - KB + 1
              B(K) = (B(K) - D(K)*B(K+1) - E(K)*B(K+2))/C(K)
   60          CONTINUE
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
!
  return
end
