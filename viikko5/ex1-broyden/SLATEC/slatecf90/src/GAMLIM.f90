subroutine GAMLIM (XMIN, XMAX)
!
!! GAMLIM computes minimum and maximum argument bounds for the Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A, R2
!***TYPE      SINGLE PRECISION (GAMLIM-S, DGAMLM-D)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Calculate the minimum and maximum legal bounds for X in GAMMA(X).
! XMIN and XMAX are not the only bounds, but they are the only non-
! trivial ones to calculate.
!
!             Output Arguments --
! XMIN   minimum legal value of X in GAMMA(X).  Any smaller value of
!        X might result in underflow.
! XMAX   maximum legal value of X in GAMMA(X).  Any larger value will
!        cause overflow.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  GAMLIM
!***FIRST EXECUTABLE STATEMENT  GAMLIM
  ALNSML = LOG(R1MACH(1))
  XMIN = -ALNSML
  DO 10 I=1,10
    XOLD = XMIN
    XLN = LOG(XMIN)
    XMIN = XMIN - XMIN*((XMIN+0.5)*XLN - XMIN - 0.2258 + ALNSML) &
      / (XMIN*XLN + 0.5)
    if (ABS(XMIN-XOLD) < 0.005) go to 20
 10   CONTINUE
  call XERMSG ('SLATEC', 'GAMLIM', 'UNABLE TO FIND XMIN', 1, 2)
!
 20   XMIN = -XMIN + 0.01
!
  ALNBIG = LOG(R1MACH(2))
  XMAX = ALNBIG
  DO 30 I=1,10
    XOLD = XMAX
    XLN = LOG(XMAX)
    XMAX = XMAX - XMAX*((XMAX-0.5)*XLN - XMAX + 0.9189 - ALNBIG) &
      / (XMAX*XLN - 0.5)
    if (ABS(XMAX-XOLD) < 0.005) go to 40
 30   CONTINUE
  call XERMSG ('SLATEC', 'GAMLIM', 'UNABLE TO FIND XMAX', 2, 2)
!
 40   XMAX = XMAX - 0.01
  XMIN = MAX (XMIN, -XMAX+1.)
!
  return
end
