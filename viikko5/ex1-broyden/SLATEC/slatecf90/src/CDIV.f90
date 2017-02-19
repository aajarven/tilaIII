subroutine CDIV (AR, AI, BR, BI, CR, CI)
!
!! CDIV computes the complex quotient of two complex numbers.
!
!***LIBRARY   SLATEC
!***TYPE      COMPLEX (CDIV-C)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Complex division, (CR,CI) = (AR,AI)/(BR,BI)
!
!***SEE ALSO  EISDOC
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   811101  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  CDIV
  REAL AR,AI,BR,BI,CR,CI
!
  REAL S,ARS,AIS,BRS,BIS
!***FIRST EXECUTABLE STATEMENT  CDIV
  S = ABS(BR) + ABS(BI)
  ARS = AR/S
  AIS = AI/S
  BRS = BR/S
  BIS = BI/S
  S = BRS**2 + BIS**2
  CR = (ARS*BRS + AIS*BIS)/S
  CI = (AIS*BRS - ARS*BIS)/S
  return
end
