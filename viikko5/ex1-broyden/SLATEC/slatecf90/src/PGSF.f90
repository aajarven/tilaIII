function PGSF (X, IZ, C, A, BH)
!
!! PGSF is subsidiary to CBLKTR.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PGSF-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PGSF
  DIMENSION       A(*)       ,C(*)       ,BH(*)
!***FIRST EXECUTABLE STATEMENT  PGSF
  FSG = 1.
  HSG = 1.
  DO 101 J=1,IZ
     DD = 1./(X-BH(J))
     FSG = FSG*A(J)*DD
     HSG = HSG*C(J)*DD
  101 CONTINUE
  if (MOD(IZ,2)) 103,102,103
  102 PGSF = 1.-FSG-HSG
  return
  103 PGSF = 1.+FSG+HSG
  return
end
