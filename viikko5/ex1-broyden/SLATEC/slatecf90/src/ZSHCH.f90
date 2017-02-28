subroutine ZSHCH (ZR, ZI, CSHR, CSHI, CCHR, CCHI)
!
!! ZSHCH is subsidiary to ZBESH and ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CSHCH-A, ZSHCH-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
!     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
!
!***SEE ALSO  ZBESH, ZBESK
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZSHCH
!
  DOUBLE PRECISION CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, ZI, ZR
!***FIRST EXECUTABLE STATEMENT  ZSHCH
  SH = SINH(ZR)
  CH = COSH(ZR)
  SN = SIN(ZI)
  CN = COS(ZI)
  CSHR = SH*CN
  CSHI = CH*SN
  CCHR = CH*CN
  CCHI = SH*SN
  return
end