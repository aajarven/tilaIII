function ALNGAM (X)
!
!! ALNGAM computes the logarithm of the absolute value of the Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      SINGLE PRECISION (ALNGAM-S, DLNGAM-D, CLNGAM-C)
!***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! ALNGAM(X) computes the logarithm of the absolute value of the
! gamma function at X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  GAMMA, R1MACH, R9LGMC, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  ALNGAM
  LOGICAL FIRST
  EXTERNAL GAMMA
  SAVE SQ2PIL, SQPI2L, PI, XMAX, DXREL, FIRST
  DATA SQ2PIL / 0.91893853320467274E0/
  DATA SQPI2L / 0.22579135264472743E0/
  DATA PI     / 3.14159265358979324E0/
  DATA FIRST  /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ALNGAM
  if (FIRST) THEN
     XMAX = R1MACH(2)/LOG(R1MACH(2))
     DXREL = SQRT (R1MACH(4))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 10.0) go to 20
!
! LOG (ABS (GAMMA(X))) FOR  ABS(X)  <=  10.0
!
  ALNGAM = LOG (ABS (GAMMA(X)))
  return
!
! LOG (ABS (GAMMA(X))) FOR ABS(X)  >  10.0
!
 20   if (Y  >  XMAX) call XERMSG ('SLATEC', 'ALNGAM', &
     'ABS(X) SO BIG ALNGAM OVERFLOWS', 2, 2)
!
  if (X > 0.) ALNGAM = SQ2PIL + (X-0.5)*LOG(X) - X + R9LGMC(Y)
  if (X > 0.) RETURN
!
  SINPIY = ABS (SIN(PI*Y))
  if (SINPIY  ==  0.) call XERMSG ('SLATEC', 'ALNGAM', &
     'X IS A NEGATIVE INTEGER', 3, 2)
!
  if (ABS((X-AINT(X-0.5))/X)  <  DXREL) call XERMSG ('SLATEC', &
     'ALNGAM', 'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR ' // &
     'NEGATIVE INTEGER', 1, 1)
!
  ALNGAM = SQPI2L + (X-0.5)*LOG(Y) - X - LOG(SINPIY) - R9LGMC(Y)
  return
!
end
