function BCRH (XLL, XRR, IZ, C, A, BH, F, SGN)
!
!! BCRH is subsidiary to CBLKTR
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (BCRH-S, BSRH-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    CCBLK
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  BCRH
  DIMENSION       A(*)       ,C(*)       ,BH(*)
  COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        , &
                  NM         ,NCMPLX     ,IK
!***FIRST EXECUTABLE STATEMENT  BCRH
  XL = XLL
  XR = XRR
  DX = .5*ABS(XR-XL)
  101 X = .5*(XL+XR)
  if (SGN*F(X,IZ,C,A,BH)) 103,105,102
  102 XR = X
  go to 104
  103 XL = X
  104 DX = .5*DX
  if (DX-CNV) 105,105,101
  105 BCRH = .5*(XL+XR)
  return
end