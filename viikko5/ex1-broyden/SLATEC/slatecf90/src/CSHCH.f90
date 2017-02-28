subroutine CSHCH (Z, CSH, CCH)
!
!! CSHCH is subsidiary to CBESH and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CSHCH-A, ZSHCH-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
!     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
!
!***SEE ALSO  CBESH, CBESK
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CSHCH
  COMPLEX CCH, CSH, Z
  REAL CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, X, Y
!***FIRST EXECUTABLE STATEMENT  CSHCH
  X = REAL(Z)
  Y = AIMAG(Z)
  SH = SINH(X)
  CH = COSH(X)
  SN = SIN(Y)
  CN = COS(Y)
  CSHR = SH*CN
  CSHI = CH*SN
  CSH = CMPLX(CSHR,CSHI)
  CCHR = CH*CN
  CCHI = SH*SN
  CCH = CMPLX(CCHR,CCHI)
  return
end