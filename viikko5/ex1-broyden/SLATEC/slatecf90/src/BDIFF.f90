subroutine BDIFF (L, V)
!
!! BDIFF is subsidiary to BSKIN
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (BDIFF-S, DBDIFF-D)
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     BDIFF computes the sum of B(L,K)*V(K)*(-1)**K where B(L,K)
!     are the binomial coefficients.  Truncated sums are computed by
!     setting last part of the V vector to zero. On return, the binomial
!     sum is in V(L).
!
!***SEE ALSO  BSKIN
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   820601  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  BDIFF
  INTEGER I, J, K, L
  REAL V
  DIMENSION V(*)
!***FIRST EXECUTABLE STATEMENT  BDIFF
  if (L == 1) RETURN
  DO 20 J=2,L
    K = L
    DO 10 I=J,L
      V(K) = V(K-1) - V(K)
      K = K - 1
   10   CONTINUE
   20 CONTINUE
  return
end
