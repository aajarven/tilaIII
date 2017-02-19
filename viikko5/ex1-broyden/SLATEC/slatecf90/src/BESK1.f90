function BESK1 (X)
!
!! BESK1 computes the modified (hyperbolic) Bessel function of the
!  third kind of order one.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      SINGLE PRECISION (BESK1-S, DBESK1-D)
!***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BESK1(X) computes the modified (hyperbolic) Bessel function of third
! kind of order one for real argument X, where X  >  0.
!
! Series for BK1        on the interval  0.          to  4.00000D+00
!                                        with weighted error   7.02E-18
!                                         log weighted error  17.15
!                               significant figures required  16.73
!                                    decimal places required  17.67
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  BESI1, BESK1E, CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  BESK1
  DIMENSION BK1CS(11)
  LOGICAL FIRST
  SAVE BK1CS, NTK1, XMIN, XSML, XMAX, FIRST
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
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BESK1
  if (FIRST) THEN
     NTK1 = INITS (BK1CS, 11, 0.1*R1MACH(3))
     XMIN = EXP (MAX(LOG(R1MACH(1)), -LOG(R1MACH(2))) + .01)
     XSML = SQRT (4.0*R1MACH(3))
     XMAXT = -LOG(R1MACH(1))
     XMAX = XMAXT - 0.5*XMAXT*LOG(XMAXT)/(XMAXT+0.5)
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.) call XERMSG ('SLATEC', 'BESK1', &
     'X IS ZERO OR NEGATIVE', 2, 2)
  if (X > 2.0) go to 20
!
  if (X  <  XMIN) call XERMSG ('SLATEC', 'BESK1', &
     'X SO SMALL K1 OVERFLOWS', 3, 2)
  Y = 0.
  if (X > XSML) Y = X*X
  BESK1 = LOG(0.5*X)*BESI1(X) + &
    (0.75 + CSEVL (.5*Y-1., BK1CS, NTK1))/X
  return
!
 20   BESK1 = 0.
  if (X  >  XMAX) call XERMSG ('SLATEC', 'BESK1', &
     'X SO BIG K1 UNDERFLOWS', 1, 1)
  if (X > XMAX) RETURN
!
  BESK1 = EXP(-X) * BESK1E(X)
!
  return
end
