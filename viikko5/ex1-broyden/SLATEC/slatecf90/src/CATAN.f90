FUNCTION CATAN (Z)
!
!! CATAN computes the complex arc tangent.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4A
!***TYPE      COMPLEX (CATAN-C)
!***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CATAN(Z) calculates the complex trigonometric arc tangent of Z.
! The result is in units of radians, and the real part is in the first
! or fourth quadrant.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  CATAN
  COMPLEX CATAN
  COMPLEX Z, Z2
  LOGICAL FIRST
  SAVE PI2, NTERMS, SQEPS, RMIN, RMAX, FIRST
  DATA PI2 / 1.57079632679489661923E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  CATAN
  if (FIRST) THEN
! NTERMS = LOG(EPS)/LOG(RBND) WHERE RBND = 0.1
     NTERMS = -0.4343*LOG(R1MACH(3)) + 1.0
     SQEPS = SQRT(R1MACH(4))
     RMIN = SQRT (3.0*R1MACH(3))
     RMAX = 1.0/R1MACH(3)
  end if
  FIRST = .FALSE.
!
  R = ABS(Z)
  if (R > 0.1) go to 30
!
  CATAN = Z
  if (R < RMIN) RETURN
!
  CATAN = (0.0, 0.0)
  Z2 = Z*Z
  DO 20 I=1,NTERMS
    TWOI = 2*(NTERMS-I) + 1
    CATAN = 1.0/TWOI - Z2*CATAN
 20   CONTINUE
  CATAN = Z*CATAN
  return
!
 30   if (R > RMAX) go to 50
  X = REAL(Z)
  Y = AIMAG(Z)
  R2 = R*R
  if (R2  ==  1.0 .AND. X  ==  0.0) call XERMSG ('SLATEC', 'CATAN', &
     'Z IS +I OR -I', 2, 2)
  if (ABS(R2-1.0) > SQEPS) go to 40
  if (ABS(CMPLX(1.0, 0.0)+Z*Z)  <  SQEPS) call XERMSG ('SLATEC', &
     'CATAN', 'ANSWER LT HALF PRECISION, Z**2 CLOSE TO -1', 1, 1)
!
 40   XANS = 0.5*ATAN2(2.0*X, 1.0-R2)
  YANS = 0.25*LOG((R2+2.0*Y+1.0)/(R2-2.0*Y+1.0))
  CATAN = CMPLX (XANS, YANS)
  return
!
 50   CATAN = CMPLX (PI2, 0.)
  if (REAL(Z) < 0.0) CATAN = CMPLX(-PI2,0.0)
  return
!
end
