function BESK1E (X)
!
!! BESK1E computes the exponentially scaled modified (hyperbolic)
!  Bessel function of the third kind of order one.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      SINGLE PRECISION (BESK1E-S, DBSK1E-D)
!***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BESK1E(X) computes the exponentially scaled modified (hyperbolic)
! Bessel function of third kind of order one for real argument
! X  >  0.0, i.e., EXP(X)*K1(X).
!
! Series for BK1        on the interval  0.          to  4.00000D+00
!                                        with weighted error   7.02E-18
!                                         log weighted error  17.15
!                               significant figures required  16.73
!                                    decimal places required  17.67
!
! Series for AK1        on the interval  1.25000D-01 to  5.00000D-01
!                                        with weighted error   6.06E-17
!                                         log weighted error  16.22
!                               significant figures required  15.41
!                                    decimal places required  16.83
!
! Series for AK12       on the interval  0.          to  1.25000D-01
!                                        with weighted error   2.58E-17
!                                         log weighted error  16.59
!                               significant figures required  15.22
!                                    decimal places required  17.16
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  BESI1, CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  BESK1E
  DIMENSION BK1CS(11), AK1CS(17), AK12CS(14)
  LOGICAL FIRST
  SAVE BK1CS, AK1CS, AK12CS, NTK1, NTAK1, NTAK12, XMIN, XSML, &
   FIRST
  DATA BK1CS( 1) /    .0253002273389477705E0 /
  DATA BK1CS( 2) /   -.353155960776544876E0 /
  DATA BK1CS( 3) /   -.122611180822657148E0 /
  DATA BK1CS( 4) /   -.0069757238596398643E0 /
  DATA BK1CS( 5) /   -.0001730288957513052E0 /
  DATA BK1CS( 6) /   -.0000024334061415659E0 /
  DATA BK1CS( 7) /   -.0000000221338763073E0 /
  DATA BK1CS( 8) /   -.0000000001411488392E0 /
  DATA BK1CS( 9) /   -.0000000000006666901E0 /
  DATA BK1CS(10) /   -.0000000000000024274E0 /
  DATA BK1CS(11) /   -.0000000000000000070E0 /
  DATA AK1CS( 1) /    .2744313406973883E0 /
  DATA AK1CS( 2) /    .07571989953199368E0 /
  DATA AK1CS( 3) /   -.00144105155647540E0 /
  DATA AK1CS( 4) /    .00006650116955125E0 /
  DATA AK1CS( 5) /   -.00000436998470952E0 /
  DATA AK1CS( 6) /    .00000035402774997E0 /
  DATA AK1CS( 7) /   -.00000003311163779E0 /
  DATA AK1CS( 8) /    .00000000344597758E0 /
  DATA AK1CS( 9) /   -.00000000038989323E0 /
  DATA AK1CS(10) /    .00000000004720819E0 /
  DATA AK1CS(11) /   -.00000000000604783E0 /
  DATA AK1CS(12) /    .00000000000081284E0 /
  DATA AK1CS(13) /   -.00000000000011386E0 /
  DATA AK1CS(14) /    .00000000000001654E0 /
  DATA AK1CS(15) /   -.00000000000000248E0 /
  DATA AK1CS(16) /    .00000000000000038E0 /
  DATA AK1CS(17) /   -.00000000000000006E0 /
  DATA AK12CS( 1) /    .06379308343739001E0 /
  DATA AK12CS( 2) /    .02832887813049721E0 /
  DATA AK12CS( 3) /   -.00024753706739052E0 /
  DATA AK12CS( 4) /    .00000577197245160E0 /
  DATA AK12CS( 5) /   -.00000020689392195E0 /
  DATA AK12CS( 6) /    .00000000973998344E0 /
  DATA AK12CS( 7) /   -.00000000055853361E0 /
  DATA AK12CS( 8) /    .00000000003732996E0 /
  DATA AK12CS( 9) /   -.00000000000282505E0 /
  DATA AK12CS(10) /    .00000000000023720E0 /
  DATA AK12CS(11) /   -.00000000000002176E0 /
  DATA AK12CS(12) /    .00000000000000215E0 /
  DATA AK12CS(13) /   -.00000000000000022E0 /
  DATA AK12CS(14) /    .00000000000000002E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BESK1E
  if (FIRST) THEN
     NTK1 = INITS (BK1CS, 11, 0.1*R1MACH(3))
     NTAK1 = INITS (AK1CS, 17, 0.1*R1MACH(3))
     NTAK12 = INITS (AK12CS, 14, 0.1*R1MACH(3))
!
     XMIN = EXP (MAX(LOG(R1MACH(1)), -LOG(R1MACH(2))) + .01)
     XSML = SQRT (4.0*R1MACH(3))
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.) call XERMSG ('SLATEC', 'BESK1E', &
     'X IS ZERO OR NEGATIVE', 2, 2)
  if (X > 2.0) go to 20
!
  if (X  <  XMIN) call XERMSG ('SLATEC', 'BESK1E', &
     'X SO SMALL K1 OVERFLOWS', 3, 2)
  Y = 0.
  if (X > XSML) Y = X*X
  BESK1E = EXP(X) * (LOG(0.5*X)*BESI1(X) + &
    (0.75 + CSEVL (.5*Y-1., BK1CS, NTK1))/X )
  return
!
 20   if (X <= 8.) BESK1E = (1.25 + CSEVL ((16./X-5.)/3., AK1CS, NTAK1)) &
    / SQRT(X)
  if (X > 8.) BESK1E = (1.25 + CSEVL (16./X-1., AK12CS, NTAK12)) &
    / SQRT(X)
!
  return
end
