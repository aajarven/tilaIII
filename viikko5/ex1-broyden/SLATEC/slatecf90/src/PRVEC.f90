function PRVEC (M, U, V)
!
!! PRVEC is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PRVEC-S, DPRVEC-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!  This subroutine computes the inner product of a vector U
!  with the imaginary product or mate vector corresponding to V
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  SDOT
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  PRVEC
!
  DIMENSION U(*),V(*)
!***FIRST EXECUTABLE STATEMENT  PRVEC
  N=M/2
  NP=N+1
  VP=SDOT(N,U(1),1,V(NP),1)
  PRVEC=SDOT(N,U(NP),1,V(1),1) - VP
  return
end
