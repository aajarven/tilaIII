function BETA (A, B)
!
!! BETA computes the complete Beta function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7B
!***TYPE      SINGLE PRECISION (BETA-S, DBETA-D, CBETA-C)
!***KEYWORDS  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BETA computes the complete beta function.
!
! Input Parameters:
!       A   real and positive
!       B   real and positive
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  ALBETA, GAMLIM, GAMMA, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  BETA
  EXTERNAL GAMMA
  SAVE XMAX, ALNSML
  DATA XMAX, ALNSML /0., 0./
!***FIRST EXECUTABLE STATEMENT  BETA
  if (ALNSML /= 0.0) go to 10
  call GAMLIM (XMIN, XMAX)
  ALNSML = LOG(R1MACH(1))

 10   if (A  <=  0. .OR. B  <=  0.) call XERMSG ('SLATEC', 'BETA', &
     'BOTH ARGUMENTS MUST BE GT 0', 2, 2)

  if (A+B < XMAX) BETA = GAMMA(A) * GAMMA(B) / GAMMA(A+B)
  if (A+B < XMAX) RETURN

  BETA = ALBETA (A, B)
  if (BETA  <  ALNSML) call XERMSG ('SLATEC', 'BETA', &
     'A AND/OR B SO BIG BETA UNDERFLOWS', 1, 2)

  BETA = EXP (BETA)

  return
end
