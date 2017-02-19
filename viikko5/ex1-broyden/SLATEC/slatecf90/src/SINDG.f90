function SINDG (X)
!
!! SINDG computes the sine of an argument in degrees.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4A
!***TYPE      SINGLE PRECISION (SINDG-S, DSINDG-D)
!***KEYWORDS  DEGREES, ELEMENTARY FUNCTIONS, FNLIB, SINE, TRIGONOMETRIC
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! SINDG(X) evaluates the single precision sine of X where
! X is in degrees.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  SINDG
! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
  SAVE RADDEG
  DATA RADDEG / .017453292519943296E0 /
!
!***FIRST EXECUTABLE STATEMENT  SINDG
  SINDG = SIN (RADDEG*X)
!
  if (MOD(X,90.) /= 0.) RETURN
  N = ABS(X)/90.0 + 0.5
  N = MOD (N, 2)
  if (N == 0) SINDG = 0.
  if (N == 1) SINDG = SIGN (1.0, SINDG)
!
  return
end
