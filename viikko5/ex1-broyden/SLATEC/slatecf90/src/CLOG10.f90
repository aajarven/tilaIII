FUNCTION CLOG10 (Z)
!
!! CLOG10 computes the principal value of the complex base 10 logarithm.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4B
!***TYPE      COMPLEX (CLOG10-C)
!***KEYWORDS  BASE TEN LOGARITHM, ELEMENTARY FUNCTIONS, FNLIB
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CLOG10(Z) calculates the principal value of the complex common
! or base 10 logarithm of Z for -PI  <  arg(Z)  <=  +PI.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CLOG10
  COMPLEX CLOG10
  COMPLEX Z
  SAVE ALOGE
  DATA ALOGE / 0.43429448190325182765E0 /
!***FIRST EXECUTABLE STATEMENT  CLOG10
  CLOG10 = ALOGE * LOG(Z)
!
  return
end
