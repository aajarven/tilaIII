FUNCTION CGAMMA (Z)
!
!! CGAMMA computes the complete Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      COMPLEX (GAMMA-S, DGAMMA-D, CGAMMA-C)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CGAMMA(Z) calculates the complete gamma function for COMPLEX
! argument Z.  This is a preliminary version that is portable
! but not accurate.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CLNGAM
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CGAMMA
  COMPLEX CGAMMA
  COMPLEX Z, CLNGAM
!***FIRST EXECUTABLE STATEMENT  CGAMMA
  CGAMMA = EXP (CLNGAM(Z))
!
  return
end
