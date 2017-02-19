function PSGF (X, IZ, C, A, BH)
!
!! PSGF is subsidiary to BLKTRI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PSGF-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  BLKTRI
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PSGF
  DIMENSION       A(*)       ,C(*)       ,BH(*)
!***FIRST EXECUTABLE STATEMENT  PSGF
  FSG = 1.
  HSG = 1.
  DO 101 J=1,IZ
     DD = 1./(X-BH(J))
     FSG = FSG*A(J)*DD
     HSG = HSG*C(J)*DD
  101 CONTINUE
  if (MOD(IZ,2)) 103,102,103
  102 PSGF = 1.-FSG-HSG
  return
  103 PSGF = 1.+FSG+HSG
  return
end
