function AI (X)
!
!! AI evaluates the Airy function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10D
!***TYPE      SINGLE PRECISION (AI-S, DAI-D)
!***KEYWORDS  AIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! AI(X) computes the Airy function Ai(X)
! Series for AIF        on the interval -1.00000D+00 to  1.00000D+00
!                                        with weighted error   1.09E-19
!                                         log weighted error  18.96
!                               significant figures required  17.76
!                                    decimal places required  19.44
!
! Series for AIG        on the interval -1.00000D+00 to  1.00000D+00
!                                        with weighted error   1.51E-17
!                                         log weighted error  16.82
!                               significant figures required  15.19
!                                    decimal places required  17.27
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  AIE, CSEVL, INITS, R1MACH, R9AIMP, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920618  Removed space from variable names.  (RWC, WRB)
!***END PROLOGUE  AI
  DIMENSION AIFCS(9), AIGCS(8)
  LOGICAL FIRST
  SAVE AIFCS, AIGCS, NAIF, NAIG, X3SML, XMAX, FIRST
  DATA AIFCS( 1) /   -.03797135849666999750E0 /
  DATA AIFCS( 2) /    .05919188853726363857E0 /
  DATA AIFCS( 3) /    .00098629280577279975E0 /
  DATA AIFCS( 4) /    .00000684884381907656E0 /
  DATA AIFCS( 5) /    .00000002594202596219E0 /
  DATA AIFCS( 6) /    .00000000006176612774E0 /
  DATA AIFCS( 7) /    .00000000000010092454E0 /
  DATA AIFCS( 8) /    .00000000000000012014E0 /
  DATA AIFCS( 9) /    .00000000000000000010E0 /
  DATA AIGCS( 1) /    .01815236558116127E0 /
  DATA AIGCS( 2) /    .02157256316601076E0 /
  DATA AIGCS( 3) /    .00025678356987483E0 /
  DATA AIGCS( 4) /    .00000142652141197E0 /
  DATA AIGCS( 5) /    .00000000457211492E0 /
  DATA AIGCS( 6) /    .00000000000952517E0 /
  DATA AIGCS( 7) /    .00000000000001392E0 /
  DATA AIGCS( 8) /    .00000000000000001E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  AI
  if (FIRST) THEN
     NAIF = INITS (AIFCS, 9, 0.1*R1MACH(3))
     NAIG = INITS (AIGCS, 8, 0.1*R1MACH(3))
!
     X3SML = R1MACH(3)**0.3334
     XMAXT = (-1.5*LOG(R1MACH(1)))**0.6667
     XMAX = XMAXT - XMAXT*LOG(XMAXT)/ &
                     (4.0*SQRT(XMAXT)+1.0) - 0.01
  end if
  FIRST = .FALSE.
!
  if (X >= (-1.0)) go to 20
  call R9AIMP (X, XM, THETA)
  AI = XM * COS(THETA)
  return
!
 20   if (X > 1.0) go to 30
  Z = 0.0
  if (ABS(X) > X3SML) Z = X**3
  AI = 0.375 + (CSEVL (Z, AIFCS, NAIF) - X*(0.25 + &
    CSEVL (Z, AIGCS, NAIG)) )
  return
!
 30   if (X > XMAX) go to 40
  AI = AIE(X) * EXP(-2.0*X*SQRT(X)/3.0)
  return
!
 40   AI = 0.0
  call XERMSG ('SLATEC', 'AI', 'X SO BIG AI UNDERFLOWS', 1, 1)
  return
!
end
