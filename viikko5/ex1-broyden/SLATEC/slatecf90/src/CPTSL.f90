subroutine CPTSL (N, D, E, B)
!
!! CPTSL solves a positive definite tridiagonal linear system.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2D2A
!***TYPE      COMPLEX (SPTSL-S, DPTSL-D, CPTSL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, POSITIVE DEFINITE, SOLVE,
!             TRIDIAGONAL
!***AUTHOR  Dongarra, J., (ANL)
!***DESCRIPTION
!
!     CPTSL given a positive definite tridiagonal matrix and a right
!     hand side will find the solution.
!
!     On Entry
!
!        N        INTEGER
!                 is the order of the tridiagonal matrix.
!
!        D        COMPLEX(N)
!                 is the diagonal of the tridiagonal matrix.
!                 On output D is destroyed.
!
!        E        COMPLEX(N)
!                 is the offdiagonal of the tridiagonal matrix.
!                 E(1) through E(N-1) should contain the
!                 offdiagonal.
!
!        B        COMPLEX(N)
!                 is the right hand side vector.
!
!     On Return
!
!        B        contains the solution.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890505  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CPTSL
  INTEGER N
  COMPLEX D(*),E(*),B(*)
!
  INTEGER K,KBM1,KE,KF,KP1,NM1,NM1D2
  COMPLEX T1,T2
!
!     CHECK FOR 1 X 1 CASE
!
!***FIRST EXECUTABLE STATEMENT  CPTSL
  if (N  /=  1) go to 10
     B(1) = B(1)/D(1)
  go to 70
   10 CONTINUE
     NM1 = N - 1
     NM1D2 = NM1/2
     if (N  ==  2) go to 30
        KBM1 = N - 1
!
!           ZERO TOP HALF OF SUBDIAGONAL AND BOTTOM HALF OF
!           SUPERDIAGONAL
!
        DO 20 K = 1, NM1D2
           T1 = CONJG(E(K))/D(K)
           D(K+1) = D(K+1) - T1*E(K)
           B(K+1) = B(K+1) - T1*B(K)
           T2 = E(KBM1)/D(KBM1+1)
           D(KBM1) = D(KBM1) - T2*CONJG(E(KBM1))
           B(KBM1) = B(KBM1) - T2*B(KBM1+1)
           KBM1 = KBM1 - 1
   20       CONTINUE
   30    CONTINUE
     KP1 = NM1D2 + 1
!
!        CLEAN UP FOR POSSIBLE 2 X 2 BLOCK AT CENTER
!
     if (MOD(N,2)  /=  0) go to 40
        T1 = CONJG(E(KP1))/D(KP1)
        D(KP1+1) = D(KP1+1) - T1*E(KP1)
        B(KP1+1) = B(KP1+1) - T1*B(KP1)
        KP1 = KP1 + 1
   40    CONTINUE
!
!        BACK SOLVE STARTING AT THE CENTER, GOING TOWARDS THE TOP
!        AND BOTTOM
!
     B(KP1) = B(KP1)/D(KP1)
     if (N  ==  2) go to 60
        K = KP1 - 1
        KE = KP1 + NM1D2 - 1
        DO 50 KF = KP1, KE
           B(K) = (B(K) - E(K)*B(K+1))/D(K)
           B(KF+1) = (B(KF+1) - CONJG(E(KF))*B(KF))/D(KF+1)
           K = K - 1
   50       CONTINUE
   60    CONTINUE
     if (MOD(N,2)  ==  0) B(1) = (B(1) - E(1)*B(2))/D(1)
   70 CONTINUE
  return
end
