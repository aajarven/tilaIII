function PPPSF (X, IZ, C, A, BH)
!
!! PPPSF is subsidiary to CBLKTR.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PPPSF-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PPPSF
  DIMENSION       A(*)       ,C(*)       ,BH(*)
!***FIRST EXECUTABLE STATEMENT  PPPSF
  SUM = 0.
  DO 101 J=1,IZ
     SUM = SUM+1./(X-BH(J))
  101 CONTINUE
  PPPSF = SUM
  return
end
