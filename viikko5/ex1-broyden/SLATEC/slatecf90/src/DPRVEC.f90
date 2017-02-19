  DOUBLE PRECISION FUNCTION DPRVEC (M, U, V)
!
!! DPRVEC is subsidiary to DBVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (PRVEC-S, DPRVEC-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!  This subroutine computes the inner product of a vector U
!  with the imaginary product or mate vector corresponding to V.
!
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  DDOT
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DPRVEC
!
  DOUBLE PRECISION DDOT
  INTEGER M, N, NP
  DOUBLE PRECISION U(*), V(*), VP
!***FIRST EXECUTABLE STATEMENT  DPRVEC
  N = M/2
  NP = N + 1
  VP = DDOT(N,U(1),1,V(NP),1)
  DPRVEC = DDOT(N,U(NP),1,V(1),1) - VP
  return
end
