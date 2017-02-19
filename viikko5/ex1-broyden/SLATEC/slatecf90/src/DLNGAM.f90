  DOUBLE PRECISION FUNCTION DLNGAM (X)
!
!! DLNGAM computes the logarithm of the absolute value of the Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      DOUBLE PRECISION (ALNGAM-S, DLNGAM-D, CLNGAM-C)
!***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DLNGAM(X) calculates the double precision logarithm of the
! absolute value of the Gamma function for double precision
! argument X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, D9LGMC, DGAMMA, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  DLNGAM
  DOUBLE PRECISION X, DXREL, PI, SINPIY, SQPI2L, SQ2PIL, XMAX, &
    Y, DGAMMA, D9LGMC, D1MACH, TEMP
  LOGICAL FIRST
  EXTERNAL DGAMMA
  SAVE SQ2PIL, SQPI2L, PI, XMAX, DXREL, FIRST
  DATA SQ2PIL / 0.91893853320467274178032973640562D0 /
  DATA SQPI2L / +.225791352644727432363097614947441D+0    /
  DATA PI / 3.14159265358979323846264338327950D0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DLNGAM
  if (FIRST) THEN
     TEMP = 1.D0/LOG(D1MACH(2))
     XMAX = TEMP*D1MACH(2)
     DXREL = SQRT(D1MACH(4))
  end if
  FIRST = .FALSE.
!
  Y = ABS (X)
  if (Y > 10.D0) go to 20
!
! LOG (ABS (DGAMMA(X)) ) FOR ABS(X)  <=  10.0
!
  DLNGAM = LOG (ABS (DGAMMA(X)) )
  return
!
! LOG ( ABS (DGAMMA(X)) ) FOR ABS(X)  >  10.0
!
 20   if (Y  >  XMAX) call XERMSG ('SLATEC', 'DLNGAM', &
     'ABS(X) SO BIG DLNGAM OVERFLOWS', 2, 2)
!
  if (X > 0.D0) DLNGAM = SQ2PIL + (X-0.5D0)*LOG(X) - X + D9LGMC(Y)
  if (X > 0.D0) RETURN
!
  SINPIY = ABS (SIN(PI*Y))
  if (SINPIY  ==  0.D0) call XERMSG ('SLATEC', 'DLNGAM', &
     'X IS A NEGATIVE INTEGER', 3, 2)
!
  if (ABS((X-AINT(X-0.5D0))/X)  <  DXREL) call XERMSG ('SLATEC', &
     'DLNGAM', &
     'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', &
     1, 1)
!
  DLNGAM = SQPI2L + (X-0.5D0)*LOG(Y) - X - LOG(SINPIY) - D9LGMC(Y)
  return
!
end
