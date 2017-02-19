function BESK0E (X)
!
!! BESK0E computes the exponentially scaled modified (hyperbolic)
!  Bessel function of the third kind of order zero.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      SINGLE PRECISION (BESK0E-S, DBSK0E-D)
!***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BESK0E(X) computes the exponentially scaled modified (hyperbolic)
! Bessel function of third kind of order zero for real argument
! X  >  0.0, i.e., EXP(X)*K0(X).
!
! Series for BK0        on the interval  0.          to  4.00000D+00
!                                        with weighted error   3.57E-19
!                                         log weighted error  18.45
!                               significant figures required  17.99
!                                    decimal places required  18.97
!
! Series for AK0        on the interval  1.25000D-01 to  5.00000D-01
!                                        with weighted error   5.34E-17
!                                         log weighted error  16.27
!                               significant figures required  14.92
!                                    decimal places required  16.89
!
! Series for AK02       on the interval  0.          to  1.25000D-01
!                                        with weighted error   2.34E-17
!                                         log weighted error  16.63
!                               significant figures required  14.67
!                                    decimal places required  17.20
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  BESI0, CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  BESK0E
  DIMENSION BK0CS(11), AK0CS(17), AK02CS(14)
  LOGICAL FIRST
  SAVE BK0CS, AK0CS, AK02CS, NTK0, NTAK0, NTAK02, XSML, FIRST
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
  DATA AK0CS( 1) /   -.07643947903327941E0 /
  DATA AK0CS( 2) /   -.02235652605699819E0 /
  DATA AK0CS( 3) /    .00077341811546938E0 /
  DATA AK0CS( 4) /   -.00004281006688886E0 /
  DATA AK0CS( 5) /    .00000308170017386E0 /
  DATA AK0CS( 6) /   -.00000026393672220E0 /
  DATA AK0CS( 7) /    .00000002563713036E0 /
  DATA AK0CS( 8) /   -.00000000274270554E0 /
  DATA AK0CS( 9) /    .00000000031694296E0 /
  DATA AK0CS(10) /   -.00000000003902353E0 /
  DATA AK0CS(11) /    .00000000000506804E0 /
  DATA AK0CS(12) /   -.00000000000068895E0 /
  DATA AK0CS(13) /    .00000000000009744E0 /
  DATA AK0CS(14) /   -.00000000000001427E0 /
  DATA AK0CS(15) /    .00000000000000215E0 /
  DATA AK0CS(16) /   -.00000000000000033E0 /
  DATA AK0CS(17) /    .00000000000000005E0 /
  DATA AK02CS( 1) /   -.01201869826307592E0 /
  DATA AK02CS( 2) /   -.00917485269102569E0 /
  DATA AK02CS( 3) /    .00014445509317750E0 /
  DATA AK02CS( 4) /   -.00000401361417543E0 /
  DATA AK02CS( 5) /    .00000015678318108E0 /
  DATA AK02CS( 6) /   -.00000000777011043E0 /
  DATA AK02CS( 7) /    .00000000046111825E0 /
  DATA AK02CS( 8) /   -.00000000003158592E0 /
  DATA AK02CS( 9) /    .00000000000243501E0 /
  DATA AK02CS(10) /   -.00000000000020743E0 /
  DATA AK02CS(11) /    .00000000000001925E0 /
  DATA AK02CS(12) /   -.00000000000000192E0 /
  DATA AK02CS(13) /    .00000000000000020E0 /
  DATA AK02CS(14) /   -.00000000000000002E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BESK0E
  if (FIRST) THEN
     NTK0 = INITS (BK0CS, 11, 0.1*R1MACH(3))
     NTAK0 = INITS (AK0CS, 17, 0.1*R1MACH(3))
     NTAK02 = INITS (AK02CS, 14, 0.1*R1MACH(3))
     XSML = SQRT (4.0*R1MACH(3))
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.) call XERMSG ('SLATEC', 'BESK0E', &
     'X IS ZERO OR NEGATIVE', 2, 2)
  if (X > 2.) go to 20
!
  Y = 0.
  if (X > XSML) Y = X*X
  BESK0E = EXP(X) * (-LOG(0.5*X)*BESI0(X) &
    - .25 + CSEVL (.5*Y-1., BK0CS, NTK0) )
  return
!
 20   if (X <= 8.) BESK0E = (1.25 + CSEVL ((16./X-5.)/3., AK0CS, NTAK0)) &
    / SQRT(X)
  if (X > 8.) BESK0E = (1.25 + CSEVL (16./X-1., AK02CS, NTAK02)) &
    / SQRT(X)
!
  return
end
