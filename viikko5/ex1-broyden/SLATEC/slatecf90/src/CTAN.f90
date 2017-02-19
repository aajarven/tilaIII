FUNCTION CTAN (Z)
!
!! CTAN computes the complex tangent.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4A
!***TYPE      COMPLEX (CTAN-C)
!***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, TANGENT, TRIGONOMETRIC
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CTAN(Z) calculates the complex trigonometric tangent of complex
! argument Z.  Z is in units of radians.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH, XERCLR, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  CTAN
  COMPLEX CTAN
  COMPLEX Z
  SAVE SQEPS
  DATA SQEPS /0./
!***FIRST EXECUTABLE STATEMENT  CTAN
  if (SQEPS == 0.) SQEPS = SQRT (R1MACH(4))
!
  X2 = 2.0*REAL(Z)
  Y2 = 2.0*AIMAG(Z)
!
  SN2X = SIN (X2)
  call XERCLR
!
  DEN = COS(X2) + COSH(Y2)
  if (DEN  ==  0.) call XERMSG ('SLATEC', 'CTAN', &
     'TAN IS SINGULAR FOR INPUT Z (X IS PI/2 OR 3*PI/2 AND Y IS 0)', &
     2, 2)
!
  if (ABS(DEN) > MAX(ABS(X2),1.)*SQEPS) go to 10
  call XERCLR
  call XERMSG ('SLATEC', 'CTAN', &
     'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X TOO NEAR ' // &
     'PI/2 OR 3*PI/2', 1, 1)
!
 10   CTAN = CMPLX (SN2X/DEN, SINH(Y2)/DEN)
!
  return
end
