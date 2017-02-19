subroutine ZEXP (AR, AI, BR, BI)
!
!! ZEXP is subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and ZBIRY.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (ZEXP-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A)
!
!***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZEXP
  DOUBLE PRECISION AR, AI, BR, BI, ZM, CA, CB
!***FIRST EXECUTABLE STATEMENT  ZEXP
  ZM = EXP(AR)
  CA = ZM*COS(AI)
  CB = ZM*SIN(AI)
  BR = CA
  BI = CB
  return
end
