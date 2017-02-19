function R9ATN1 (X)
!
!! R9ATN1 evaluates ATAN(X) from first order relative accuracy so that ...
!            ATAN(X) = X + X**3*R9ATN1(X).
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4A
!***TYPE      SINGLE PRECISION (R9ATN1-S, D9ATN1-D)
!***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FIRST ORDER, FNLIB,
!             TRIGONOMETRIC
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate  ATAN(X)  from first order, that is, evaluate
! (ATAN(X)-X)/X**3  with relative error accuracy so that
!        ATAN(X) = X + X**3*R9ATN1(X).
!
! Series for ATN1       on the interval  0.          to  1.00000D+00
!                                        with weighted error   2.21E-17
!                                         log weighted error  16.66
!                               significant figures required  15.44
!                                    decimal places required  17.32
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   780401  DATE WRITTEN
!   890206  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  R9ATN1
  DIMENSION ATN1CS(21)
  LOGICAL FIRST
  SAVE ATN1CS, NTATN1, XSML, XBIG, XMAX, FIRST
  DATA ATN1CS( 1) /   -.03283997535355202E0 /
  DATA ATN1CS( 2) /    .05833432343172412E0 /
  DATA ATN1CS( 3) /   -.00740036969671964E0 /
  DATA ATN1CS( 4) /    .00100978419933728E0 /
  DATA ATN1CS( 5) /   -.00014397871635652E0 /
  DATA ATN1CS( 6) /    .00002114512648992E0 /
  DATA ATN1CS( 7) /   -.00000317232107425E0 /
  DATA ATN1CS( 8) /    .00000048366203654E0 /
  DATA ATN1CS( 9) /   -.00000007467746546E0 /
  DATA ATN1CS(10) /    .00000001164800896E0 /
  DATA ATN1CS(11) /   -.00000000183208837E0 /
  DATA ATN1CS(12) /    .00000000029019082E0 /
  DATA ATN1CS(13) /   -.00000000004623885E0 /
  DATA ATN1CS(14) /    .00000000000740552E0 /
  DATA ATN1CS(15) /   -.00000000000119135E0 /
  DATA ATN1CS(16) /    .00000000000019240E0 /
  DATA ATN1CS(17) /   -.00000000000003118E0 /
  DATA ATN1CS(18) /    .00000000000000506E0 /
  DATA ATN1CS(19) /   -.00000000000000082E0 /
  DATA ATN1CS(20) /    .00000000000000013E0 /
  DATA ATN1CS(21) /   -.00000000000000002E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  R9ATN1
  if (FIRST) THEN
     EPS = R1MACH(3)
     NTATN1 = INITS (ATN1CS, 21, 0.1*EPS)
!
     XSML = SQRT (0.1*EPS)
     XBIG = 1.571/SQRT(EPS)
     XMAX = 1.571/EPS
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 1.0) go to 20
!
  if (Y <= XSML) R9ATN1 = -1.0/3.0
  if (Y <= XSML) RETURN
!
  R9ATN1 = -0.25 + CSEVL (2.0*Y*Y-1., ATN1CS, NTATN1)
  return
!
 20   if (Y  >  XMAX) call XERMSG ('SLATEC', 'R9ATN1', &
     'NO PRECISION IN ANSWER BECAUSE X IS TOO BIG', 2, 2)
  if (Y  >  XBIG) call XERMSG ('SLATEC', 'R9ATN1', &
     'ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG', 1, 1)
!
  R9ATN1 = (ATAN(X) - X) / X**3
  return
!
end
