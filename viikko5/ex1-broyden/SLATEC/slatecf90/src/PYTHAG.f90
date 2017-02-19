FUNCTION PYTHAG (A, B)
!
!! PYTHAG computes the complex square root of a complex number without ...
!            destructive overflow or underflow.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PYTHAG-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Finds sqrt(A**2+B**2) without overflow or destructive underflow
!
!***SEE ALSO  EISDOC
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   811101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PYTHAG
  REAL A,B
  REAL PYTHAG
  REAL P,Q,R,S,T
!***FIRST EXECUTABLE STATEMENT  PYTHAG
  P = MAX(ABS(A),ABS(B))
  Q = MIN(ABS(A),ABS(B))
  if (Q  ==  0.0E0) go to 20
   10 CONTINUE
     R = (Q/P)**2
     T = 4.0E0 + R
     if (T  ==  4.0E0) go to 20
     S = R/T
     P = P + 2.0E0*P*S
     Q = Q*S
  go to 10
   20 PYTHAG = P
  return
end
