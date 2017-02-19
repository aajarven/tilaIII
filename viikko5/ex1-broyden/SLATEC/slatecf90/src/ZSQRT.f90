subroutine ZSQRT (AR, AI, BR, BI)
!
!! ZSQRT is subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and ZBIRY.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (ZSQRT-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     DOUBLE PRECISION COMPLEX SQUARE ROOT, B=CSQRT(A)
!
!***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
!***ROUTINES CALLED  ZABS
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZSQRT
  DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DRT
  DOUBLE PRECISION ZABS
  EXTERNAL ZABS
  DATA DRT , DPI / 7.071067811865475244008443621D-1, &
                   3.141592653589793238462643383D+0/
!***FIRST EXECUTABLE STATEMENT  ZSQRT
  ZM = ZABS(AR,AI)
  ZM = SQRT(ZM)
  if (AR == 0.0D+0) go to 10
  if (AI == 0.0D+0) go to 20
  DTHETA = DATAN(AI/AR)
  if (DTHETA <= 0.0D+0) go to 40
  if (AR < 0.0D+0) DTHETA = DTHETA - DPI
  go to 50
   10 if (AI > 0.0D+0) go to 60
  if (AI < 0.0D+0) go to 70
  BR = 0.0D+0
  BI = 0.0D+0
  return
   20 if (AR > 0.0D+0) go to 30
  BR = 0.0D+0
  BI = SQRT(ABS(AR))
  return
   30 BR = SQRT(AR)
  BI = 0.0D+0
  return
   40 if (AR < 0.0D+0) DTHETA = DTHETA + DPI
   50 DTHETA = DTHETA*0.5D+0
  BR = ZM*COS(DTHETA)
  BI = ZM*SIN(DTHETA)
  return
   60 BR = ZM*DRT
  BI = ZM*DRT
  return
   70 BR = ZM*DRT
  BI = -ZM*DRT
  return
end
