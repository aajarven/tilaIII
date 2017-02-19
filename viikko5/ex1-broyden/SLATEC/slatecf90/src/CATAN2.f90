FUNCTION CATAN2 (CSN, CCS)
!
!! CATAN2 computes the complex arc tangent in the proper quadrant.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4A
!***TYPE      COMPLEX (CATAN2-C)
!***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FNLIB, POLAR ANGEL,
!             QUADRANT, TRIGONOMETRIC
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CATAN2(CSN,CCS) calculates the complex trigonometric arc
! tangent of the ratio CSN/CCS and returns a result whose real
! part is in the correct quadrant (within a multiple of 2*PI).  The
! result is in units of radians and the real part is between -PI
! and +PI.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CATAN, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  CATAN2
  COMPLEX CATAN2
  COMPLEX CSN, CCS, CATAN
  SAVE PI
  DATA PI / 3.14159265358979323846E0 /
!***FIRST EXECUTABLE STATEMENT  CATAN2
  if (ABS(CCS) == 0.) go to 10
!
  CATAN2 = CATAN (CSN/CCS)
  if (REAL(CCS) < 0.) CATAN2 = CATAN2 + PI
  if (REAL(CATAN2) > PI) CATAN2 = CATAN2 - 2.0*PI
  return
!
 10   if (ABS(CSN)  ==  0.) call XERMSG ('SLATEC', 'CATAN2', &
     'CALLED WITH BOTH ARGUMENTS ZERO', 1, 2)
!
  CATAN2 = CMPLX (SIGN(0.5*PI,REAL(CSN)), 0.0)
!
  return
end
