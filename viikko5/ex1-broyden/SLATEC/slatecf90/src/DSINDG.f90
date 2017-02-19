  DOUBLE PRECISION FUNCTION DSINDG (X)
!
!! DSINDG computes the sine of an argument in degrees.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4A
!***TYPE      DOUBLE PRECISION (SINDG-S, DSINDG-D)
!***KEYWORDS  DEGREES, ELEMENTARY FUNCTIONS, FNLIB, SINE, TRIGONOMETRIC
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DSINDG(X) calculates the double precision sine for double
! precision argument X where X is in degrees.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DSINDG
  DOUBLE PRECISION X, RADDEG
  SAVE RADDEG
  DATA RADDEG / 0.017453292519943295769236907684886D0 /
!***FIRST EXECUTABLE STATEMENT  DSINDG
  DSINDG = SIN (RADDEG*X)
!
  if (MOD(X,90.D0) /= 0.D0) RETURN
  N = ABS(X)/90.D0 + 0.5D0
  N = MOD (N, 2)
  if (N == 0) DSINDG = 0.D0
  if (N == 1) DSINDG = SIGN (1.0D0, DSINDG)
!
  return
end
