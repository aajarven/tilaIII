FUNCTION CCOT (Z)
!
!! CCOT computes the cotangent.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4A
!***TYPE      COMPLEX (COT-S, DCOT-D, CCOT-C)
!***KEYWORDS  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CCOT(Z) calculates the complex trigonometric cotangent of Z.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH, XERCLR, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  CCOT
  COMPLEX CCOT
  COMPLEX Z
  SAVE SQEPS
  DATA SQEPS /0./
!***FIRST EXECUTABLE STATEMENT  CCOT
  if (SQEPS == 0.) SQEPS = SQRT (R1MACH(4))
!
  X2 = 2.0*REAL(Z)
  Y2 = 2.0*AIMAG(Z)
!
  SN2X = SIN (X2)
  call XERCLR
!
  DEN = COSH(Y2) - COS(X2)
  if (DEN  ==  0.) call XERMSG ('SLATEC', 'CCOT', &
     'COT IS SINGULAR FOR INPUT Z (X IS 0 OR PI AND Y IS 0)', 2, 2)
!
  if (ABS(DEN) > MAX(ABS(X2),1.)*SQEPS) go to 10
  call XERCLR
  call XERMSG ('SLATEC', 'CCOT', &
     'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X TOO NEAR ' // &
     '0 OR PI', 1, 1)
!
 10   CCOT = CMPLX (SN2X/DEN, -SINH(Y2)/DEN)
!
  return
end
