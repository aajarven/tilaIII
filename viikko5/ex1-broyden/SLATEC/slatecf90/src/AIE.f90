function AIE (X)
!
!! AIE calculates the Airy function for a negative argument...
!  and an exponentially scaled Airy function for a non-negative argument.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10D
!***TYPE      SINGLE PRECISION (AIE-S, DAIE-D)
!***KEYWORDS  EXPONENTIALLY SCALED AIRY FUNCTION, FNLIB,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! AIE(X) computes the exponentially scaled Airy function for
! non-negative X.  It evaluates AI(X) for X  <=  0.0 and
! EXP(ZETA)*AI(X) for X  >=  0.0 where ZETA = (2.0/3.0)*(X**1.5).
!
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
! Series for AIP        on the interval  0.          to  1.00000D+00
!                                        with weighted error   5.10E-17
!                                         log weighted error  16.29
!                               significant figures required  14.41
!                                    decimal places required  17.06
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, INITS, R1MACH, R9AIMP
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890206  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920618  Removed space from variable names.  (RWC, WRB)
!***END PROLOGUE  AIE
  DIMENSION AIFCS(9), AIGCS(8), AIPCS(34)
  LOGICAL FIRST
  SAVE AIFCS, AIGCS, AIPCS, NAIF, NAIG, &
   NAIP, X3SML, X32SML, XBIG, FIRST
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
  DATA AIPCS( 1) /   -.0187519297793868E0 /
  DATA AIPCS( 2) /   -.0091443848250055E0 /
  DATA AIPCS( 3) /    .0009010457337825E0 /
  DATA AIPCS( 4) /   -.0001394184127221E0 /
  DATA AIPCS( 5) /    .0000273815815785E0 /
  DATA AIPCS( 6) /   -.0000062750421119E0 /
  DATA AIPCS( 7) /    .0000016064844184E0 /
  DATA AIPCS( 8) /   -.0000004476392158E0 /
  DATA AIPCS( 9) /    .0000001334635874E0 /
  DATA AIPCS(10) /   -.0000000420735334E0 /
  DATA AIPCS(11) /    .0000000139021990E0 /
  DATA AIPCS(12) /   -.0000000047831848E0 /
  DATA AIPCS(13) /    .0000000017047897E0 /
  DATA AIPCS(14) /   -.0000000006268389E0 /
  DATA AIPCS(15) /    .0000000002369824E0 /
  DATA AIPCS(16) /   -.0000000000918641E0 /
  DATA AIPCS(17) /    .0000000000364278E0 /
  DATA AIPCS(18) /   -.0000000000147475E0 /
  DATA AIPCS(19) /    .0000000000060851E0 /
  DATA AIPCS(20) /   -.0000000000025552E0 /
  DATA AIPCS(21) /    .0000000000010906E0 /
  DATA AIPCS(22) /   -.0000000000004725E0 /
  DATA AIPCS(23) /    .0000000000002076E0 /
  DATA AIPCS(24) /   -.0000000000000924E0 /
  DATA AIPCS(25) /    .0000000000000417E0 /
  DATA AIPCS(26) /   -.0000000000000190E0 /
  DATA AIPCS(27) /    .0000000000000087E0 /
  DATA AIPCS(28) /   -.0000000000000040E0 /
  DATA AIPCS(29) /    .0000000000000019E0 /
  DATA AIPCS(30) /   -.0000000000000009E0 /
  DATA AIPCS(31) /    .0000000000000004E0 /
  DATA AIPCS(32) /   -.0000000000000002E0 /
  DATA AIPCS(33) /    .0000000000000001E0 /
  DATA AIPCS(34) /   -.0000000000000000E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  AIE
  if (FIRST) THEN
     ETA = 0.1*R1MACH(3)
     NAIF  = INITS (AIFCS, 9, ETA)
     NAIG  = INITS (AIGCS, 8, ETA)
     NAIP  = INITS (AIPCS, 34, ETA)
!
     X3SML = ETA**0.3333
     X32SML = 1.3104*X3SML**2
     XBIG = R1MACH(2)**0.6666
  end if
  FIRST = .FALSE.
!
  if (X >= (-1.0)) go to 20
  call R9AIMP (X, XM, THETA)
  AIE = XM * COS(THETA)
  return
!
 20   if (X > 1.0) go to 30
  Z = 0.0
  if (ABS(X) > X3SML) Z = X**3
  AIE = 0.375 + (CSEVL (Z, AIFCS, NAIF) - X*(0.25 + &
    CSEVL (Z, AIGCS, NAIG)) )
  if (X > X32SML) AIE = AIE * EXP(2.0*X*SQRT(X)/3.0)
  return
!
 30   SQRTX = SQRT(X)
  Z = -1.0
  if (X < XBIG) Z = 2.0/(X*SQRTX) - 1.0
  AIE = (.28125 + CSEVL (Z, AIPCS, NAIP))/SQRT(SQRTX)
  return
!
end
