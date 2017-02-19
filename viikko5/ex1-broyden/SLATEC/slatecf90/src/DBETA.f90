  DOUBLE PRECISION FUNCTION DBETA (A, B)
!
!! DBETA computes the complete Beta function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7B
!***TYPE      DOUBLE PRECISION (BETA-S, DBETA-D, CBETA-C)
!***KEYWORDS  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBETA(A,B) calculates the double precision complete beta function
! for double precision arguments A and B.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DGAMLM, DGAMMA, DLBETA, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  DBETA
  DOUBLE PRECISION A, B, ALNSML, XMAX, XMIN, DLBETA, DGAMMA, D1MACH
  LOGICAL FIRST
  EXTERNAL DGAMMA
  SAVE XMAX, ALNSML, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBETA
  if (FIRST) THEN
     call DGAMLM (XMIN, XMAX)
     ALNSML = LOG (D1MACH(1))
  end if
  FIRST = .FALSE.
!
  if (A  <=  0.D0 .OR. B  <=  0.D0) call XERMSG ('SLATEC', 'DBETA', &
     'BOTH ARGUMENTS MUST BE GT 0', 2, 2)
!
  if (A+B < XMAX) DBETA = DGAMMA(A)*DGAMMA(B)/DGAMMA(A+B)
  if (A+B < XMAX) RETURN
!
  DBETA = DLBETA (A, B)
  if (DBETA < ALNSML) go to 20
  DBETA = EXP (DBETA)
  return
!
 20   DBETA = 0.D0
  call XERMSG ('SLATEC', 'DBETA', &
     'A AND/OR B SO BIG BETA UNDERFLOWS', 1, 1)
  return
!
end
