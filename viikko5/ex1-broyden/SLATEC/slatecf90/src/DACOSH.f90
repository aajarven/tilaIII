DOUBLE PRECISION FUNCTION DACOSH (X)
!
!! DACOSH computes the arc hyperbolic cosine.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4C
!***TYPE      DOUBLE PRECISION (ACOSH-S, DACOSH-D, CACOSH-C)
!***KEYWORDS  ACOSH, ARC HYPERBOLIC COSINE, ELEMENTARY FUNCTIONS, FNLIB,
!             INVERSE HYPERBOLIC COSINE
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DACOSH(X) calculates the double precision arc hyperbolic cosine for
! double precision argument X.  The result is returned on the
! positive branch.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DACOSH
  DOUBLE PRECISION X, DLN2, XMAX,  D1MACH
  SAVE DLN2, XMAX
  DATA DLN2 / 0.69314718055994530941723212145818D0 /
  DATA XMAX / 0.D0 /
!***FIRST EXECUTABLE STATEMENT  DACOSH
  if (XMAX == 0.D0) XMAX = 1.0D0/SQRT(D1MACH(3))
!
  if (X  <  1.D0) call XERMSG ('SLATEC', 'DACOSH', &
     'X LESS THAN 1', 1, 2)
!
  if (X < XMAX) DACOSH = LOG (X+SQRT(X*X-1.0D0))
  if (X >= XMAX) DACOSH = DLN2 + LOG(X)
!
  return
end
