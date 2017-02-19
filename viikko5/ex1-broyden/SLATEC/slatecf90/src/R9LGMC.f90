function R9LGMC (X)
!
!! R9LGMC computes the log Gamma correction factor so that ...
!  LOG(GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + R9LGMC(X).
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      SINGLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
!             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the log gamma correction factor for X  >=  10.0 so that
!  LOG (GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + R9LGMC(X)
!
! Series for ALGM       on the interval  0.          to  1.00000D-02
!                                        with weighted error   3.40E-16
!                                         log weighted error  15.47
!                               significant figures required  14.39
!                                    decimal places required  15.86
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  R9LGMC
  DIMENSION ALGMCS(6)
  LOGICAL FIRST
  SAVE ALGMCS, NALGM, XBIG, XMAX, FIRST
  DATA ALGMCS( 1) /    .166638948045186E0 /
  DATA ALGMCS( 2) /   -.0000138494817606E0 /
  DATA ALGMCS( 3) /    .0000000098108256E0 /
  DATA ALGMCS( 4) /   -.0000000000180912E0 /
  DATA ALGMCS( 5) /    .0000000000000622E0 /
  DATA ALGMCS( 6) /   -.0000000000000003E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  R9LGMC
  if (FIRST) THEN
     NALGM = INITS (ALGMCS, 6, R1MACH(3))
     XBIG = 1.0/SQRT(R1MACH(3))
     XMAX = EXP (MIN(LOG(R1MACH(2)/12.0), -LOG(12.0*R1MACH(1))) )
  end if
  FIRST = .FALSE.
!
  if (X  <  10.0) call XERMSG ('SLATEC', 'R9LGMC', &
     'X MUST BE GE 10', 1, 2)
  if (X >= XMAX) go to 20
!
  R9LGMC = 1.0/(12.0*X)
  if (X < XBIG) R9LGMC = CSEVL (2.0*(10./X)**2-1., ALGMCS, NALGM)/X
  return
!
 20   R9LGMC = 0.0
  call XERMSG ('SLATEC', 'R9LGMC', 'X SO BIG R9LGMC UNDERFLOWS', 2, &
     1)
  return
!
end
