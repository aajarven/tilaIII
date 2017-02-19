function PPGSF (X, IZ, C, A, BH)
!
!! PPGSF is subsidiary to CBLKTR.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PPGSF-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PPGSF
  DIMENSION       A(*)       ,C(*)       ,BH(*)
!***FIRST EXECUTABLE STATEMENT  PPGSF
  SUM = 0.
  DO 101 J=1,IZ
     SUM = SUM-1./(X-BH(J))**2
  101 CONTINUE
  PPGSF = SUM
  return
end
