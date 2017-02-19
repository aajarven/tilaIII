FUNCTION CSINH (Z)
!
!! CSINH computes the complex hyperbolic sine.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4C
!***TYPE      COMPLEX (CSINH-C)
!***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC SINE
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CSINH(Z) calculates the complex hyperbolic sine of complex
! argument Z.  Z is in units of radians.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CSINH
  COMPLEX CSINH
  COMPLEX Z, CI
  SAVE CI
  DATA CI /(0.,1.)/
!***FIRST EXECUTABLE STATEMENT  CSINH
  CSINH = -CI*SIN(CI*Z)
!
  return
end
