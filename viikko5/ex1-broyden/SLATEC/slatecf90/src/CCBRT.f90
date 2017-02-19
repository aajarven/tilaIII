FUNCTION CCBRT (Z)
!
!! CCBRT computes the cube root.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C2
!***TYPE      COMPLEX (CBRT-S, DCBRT-D, CCBRT-C)
!***KEYWORDS  CUBE ROOT, ELEMENTARY FUNCTIONS, FNLIB, ROOTS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CCBRT(Z) calculates the complex cube root of Z.  The principal root
! for which -PI  <  arg(Z)  <=  +PI is returned.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CARG, CBRT
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CCBRT
!
  COMPLEX CCBRT
  COMPLEX Z
!***FIRST EXECUTABLE STATEMENT  CCBRT
  THETA = CARG(Z) / 3.0
  R = CBRT (ABS(Z))
!
  CCBRT = CMPLX (R*COS(THETA), R*SIN(THETA))
!
  return
end
