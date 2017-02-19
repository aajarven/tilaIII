function PPSPF (X, IZ, C, A, BH)
!
!! PPSPF is subsidiary to BLKTRI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PPSPF-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  BLKTRI
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PPSPF
  DIMENSION       A(*)       ,C(*)       ,BH(*)
!***FIRST EXECUTABLE STATEMENT  PPSPF
  SUM = 0.
  DO 101 J=1,IZ
     SUM = SUM+1./(X-BH(J))
  101 CONTINUE
  PPSPF = SUM
  return
end
