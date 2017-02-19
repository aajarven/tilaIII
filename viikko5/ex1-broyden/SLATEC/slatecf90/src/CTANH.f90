FUNCTION CTANH (Z)
!
!! CTANH computes the complex hyperbolic tangent.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4C
!***TYPE      COMPLEX (CTANH-C)
!***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC TANGENT
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CTANH(Z) calculates the complex hyperbolic tangent of complex
! argument Z.  Z is in units of radians.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CTAN
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CTANH
  COMPLEX CTANH
  COMPLEX Z, CI, CTAN
  SAVE CI
  DATA CI /(0.,1.)/
!***FIRST EXECUTABLE STATEMENT  CTANH
  CTANH = -CI*CTAN(CI*Z)
!
  return
end
