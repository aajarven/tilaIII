subroutine MPMLP (U, V, W, J)
!
!! MPMLP performs the inner multiplication loop for MPMUL.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DQDOTA and DQDOTI
!***LIBRARY   SLATEC
!***TYPE      ALL (MPMLP-A)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
! Performs inner multiplication loop for MPMUL. Carries are not pro-
! pagated in inner loop, which saves time at the expense of space.
!
!***SEE ALSO  DQDOTA, DQDOTI
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  MPMLP
  INTEGER U(*), V(*), W
!***FIRST EXECUTABLE STATEMENT  MPMLP
  DO 10 I = 1, J
   10 U(I) = U(I) + W*V(I)
  return
end
