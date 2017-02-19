function COSDG (X)
!
!! COSDG computes the cosine of an argument in degrees.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4A
!***TYPE      SINGLE PRECISION (COSDG-S, DCOSDG-D)
!***KEYWORDS  COSINE, DEGREES, ELEMENTARY FUNCTIONS, FNLIB,
!             TRIGONOMETRIC
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! COSDG(X) evaluates the cosine for real X in degrees.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  COSDG
! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
  SAVE RADDEG
  DATA RADDEG / .017453292519943296E0 /
!
!***FIRST EXECUTABLE STATEMENT  COSDG
  COSDG = COS (RADDEG*X)
!
  if (MOD(X,90.) /= 0.) RETURN
  N = ABS(X)/90.0 + 0.5
  N = MOD (N, 2)
  if (N == 0) COSDG = SIGN (1.0, COSDG)
  if (N == 1) COSDG = 0.0
!
  return
end
