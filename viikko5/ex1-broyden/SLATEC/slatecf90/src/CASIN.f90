FUNCTION CASIN (ZINP)
!
!! CASIN computes the complex arc sine.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4A
!***TYPE      COMPLEX (CASIN-C)
!***KEYWORDS  ARC SINE, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CASIN(ZINP) calculates the complex trigonometric arc sine of ZINP.
! The result is in units of radians, and the real part is in the first
! or fourth quadrant.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CASIN
  COMPLEX CASIN
  COMPLEX ZINP, Z, Z2, SQZP1, CI
  LOGICAL FIRST
  SAVE PI2, PI, CI, NTERMS, RMIN, FIRST
  DATA PI2 /1.57079632679489661923E0/
  DATA PI /3.14159265358979324E0/
  DATA CI /(0.,1.)/
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  CASIN
  if (FIRST) THEN
! NTERMS = LOG(EPS)/LOG(RMAX)  WHERE RMAX = 0.1
     NTERMS = -0.4343*LOG(R1MACH(3))
     RMIN = SQRT (6.0*R1MACH(3))
  end if
  FIRST = .FALSE.
!
  Z = ZINP
  R = ABS (Z)
  if (R > 0.1) go to 30
!
  CASIN = Z
  if (R < RMIN) RETURN
!
  CASIN = (0.0, 0.0)
  Z2 = Z*Z
  DO 20 I=1,NTERMS
    TWOI = 2*(NTERMS-I) + 1
    CASIN = 1.0/TWOI + TWOI*CASIN*Z2/(TWOI+1.0)
 20   CONTINUE
  CASIN = Z*CASIN
  return
!
 30   if (REAL(ZINP) < 0.0) Z = -ZINP
!
  SQZP1 = SQRT (Z+1.0)
  if (AIMAG(SQZP1) < 0.) SQZP1 = -SQZP1
  CASIN = PI2 - CI * LOG (Z + SQZP1*SQRT(Z-1.0))
!
  if (REAL(CASIN) > PI2) CASIN = PI - CASIN
  if (REAL(CASIN) <= (-PI2)) CASIN = -PI - CASIN
  if (REAL(ZINP) < 0.) CASIN = -CASIN

  return
end
