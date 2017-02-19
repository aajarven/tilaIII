subroutine CSROOT (XR, XI, YR, YI)
!
!! CSROOT computes the complex square root of a complex number.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CSROOT-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     (YR,YI) = complex sqrt(XR,XI)
!
!***SEE ALSO  EISDOC
!***ROUTINES CALLED  PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   811101  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  CSROOT
  REAL XR,XI,YR,YI,S,TR,TI,PYTHAG
!
!     BRANCH CHOSEN SO THAT YR  >=  0.0 AND SIGN(YI)  ==  SIGN(XI)
!***FIRST EXECUTABLE STATEMENT  CSROOT
  TR = XR
  TI = XI
  S = SQRT(0.5E0*(PYTHAG(TR,TI) + ABS(TR)))
  if (TR  >=  0.0E0) YR = S
  if (TI  <  0.0E0) S = -S
  if (TR  <=  0.0E0) YI = S
  if (TR  <  0.0E0) YR = 0.5E0*(TI/YI)
  if (TR  >  0.0E0) YI = 0.5E0*(TI/YR)
  return
end
