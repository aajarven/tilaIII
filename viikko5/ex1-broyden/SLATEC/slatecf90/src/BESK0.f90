function BESK0 (X)
!
!! BESK0 computes the modified (hyperbolic) Bessel function of the
!  third kind of order zero.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      SINGLE PRECISION (BESK0-S, DBESK0-D)
!***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BESK0(X) calculates the modified (hyperbolic) Bessel function
! of the third kind of order zero for real argument X  >  0.0.
!
! Series for BK0        on the interval  0.          to  4.00000D+00
!                                        with weighted error   3.57E-19
!                                         log weighted error  18.45
!                               significant figures required  17.99
!                                    decimal places required  18.97
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  BESI0, BESK0E, CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  BESK0
  DIMENSION BK0CS(11)
  LOGICAL FIRST
  SAVE BK0CS, NTK0, XSML, XMAX, FIRST
  DATA BK0CS( 1) /   -.03532739323390276872E0 /
  DATA BK0CS( 2) /    .3442898999246284869E0 /
  DATA BK0CS( 3) /    .03597993651536150163E0 /
  DATA BK0CS( 4) /    .00126461541144692592E0 /
  DATA BK0CS( 5) /    .00002286212103119451E0 /
  DATA BK0CS( 6) /    .00000025347910790261E0 /
  DATA BK0CS( 7) /    .00000000190451637722E0 /
  DATA BK0CS( 8) /    .00000000001034969525E0 /
  DATA BK0CS( 9) /    .00000000000004259816E0 /
  DATA BK0CS(10) /    .00000000000000013744E0 /
  DATA BK0CS(11) /    .00000000000000000035E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BESK0
  if (FIRST) THEN
     NTK0 = INITS (BK0CS, 11, 0.1*R1MACH(3))
     XSML = SQRT (4.0*R1MACH(3))
     XMAXT = -LOG(R1MACH(1))
     XMAX = XMAXT - 0.5*XMAXT*LOG(XMAXT)/(XMAXT+0.5) - 0.01
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.) call XERMSG ('SLATEC', 'BESK0', &
     'X IS ZERO OR NEGATIVE', 2, 2)
  if (X > 2.) go to 20
!
  Y = 0.
  if (X > XSML) Y = X*X
  BESK0 = -LOG(0.5*X)*BESI0(X) - .25 + CSEVL (.5*Y-1., BK0CS, NTK0)
  return
!
 20   BESK0 = 0.
  if (X  >  XMAX) call XERMSG ('SLATEC', 'BESK0', &
     'X SO BIG K0 UNDERFLOWS', 1, 1)
  if (X > XMAX) RETURN
!
  BESK0 = EXP(-X) * BESK0E(X)
!
  return
end
