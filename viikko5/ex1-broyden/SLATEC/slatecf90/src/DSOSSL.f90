subroutine DSOSSL (K, N, L, X, C, B, M)
!
!! DSOSSL is subsidiary to DSOS.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SOSSOL-S, DSOSSL-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     DSOSSL solves an upper triangular type of linear system by back
!     substitution.
!
!     The matrix C is upper trapezoidal and stored as a linear array by
!     rows. The equations have been normalized so that the diagonal
!     entries of C are understood to be unity. The off diagonal entries
!     and the elements of the constant right hand side vector B have
!     already been stored as the negatives of the corresponding equation
!     values.
!     With each call to DSOSSL a (K-1) by (K-1) triangular system is
!     resolved. For L greater than K, column L of C is included in the
!     right hand side vector.
!
!***SEE ALSO  DSOS
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DSOSSL
!
!
  INTEGER J, JKM, K, KJ, KM, KM1, KMM1, KN, L, LK, M, N, NP1
  DOUBLE PRECISION B(*), C(*), X(*), XMAX
!
!***FIRST EXECUTABLE STATEMENT  DSOSSL
  NP1 = N + 1
  KM1 = K - 1
  LK = KM1
  if (L  ==  K) LK = K
  KN = M
!
!
  DO 40 KJ = 1, KM1
     KMM1 = K - KJ
     KM = KMM1 + 1
     XMAX = 0.0D0
     KN = KN - NP1 + KMM1
     if (KM  >  LK) go to 20
        JKM = KN
!
        DO 10 J = KM, LK
           JKM = JKM + 1
           XMAX = XMAX + C(JKM)*X(J)
   10       CONTINUE
   20    CONTINUE
!
     if (L  <=  K) go to 30
        JKM = KN + L - KMM1
        XMAX = XMAX + C(JKM)*X(L)
   30    CONTINUE
     X(KMM1) = XMAX + B(KMM1)
   40 CONTINUE
!
  return
end
