function BESI1E (X)
!
!! BESI1E computes the exponentially scaled modified (hyperbolic) ...
!  Bessel function of the first kind of order one.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      SINGLE PRECISION (BESI1E-S, DBSI1E-D)
!***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
!             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
!             ORDER ONE, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BESI1E(X) calculates the exponentially scaled modified (hyperbolic)
! Bessel function of the first kind of order one for real argument X;
! i.e., EXP(-ABS(X))*I1(X).
!
! Series for BI1        on the interval  0.          to  9.00000D+00
!                                        with weighted error   2.40E-17
!                                         log weighted error  16.62
!                               significant figures required  16.23
!                                    decimal places required  17.14
!
! Series for AI1        on the interval  1.25000D-01 to  3.33333D-01
!                                        with weighted error   6.98E-17
!                                         log weighted error  16.16
!                               significant figures required  14.53
!                                    decimal places required  16.82
!
! Series for AI12       on the interval  0.          to  1.25000D-01
!                                        with weighted error   3.55E-17
!                                         log weighted error  16.45
!                               significant figures required  14.69
!                                    decimal places required  17.12
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890210  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920618  Removed space from variable names.  (RWC, WRB)
!***END PROLOGUE  BESI1E
  DIMENSION BI1CS(11), AI1CS(21), AI12CS(22)
  LOGICAL FIRST
  SAVE BI1CS, AI1CS, AI12CS, NTI1, NTAI1, NTAI12, XMIN, XSML, FIRST
  DATA BI1CS( 1) /   -.001971713261099859E0 /
  DATA BI1CS( 2) /    .40734887667546481E0 /
  DATA BI1CS( 3) /    .034838994299959456E0 /
  DATA BI1CS( 4) /    .001545394556300123E0 /
  DATA BI1CS( 5) /    .000041888521098377E0 /
  DATA BI1CS( 6) /    .000000764902676483E0 /
  DATA BI1CS( 7) /    .000000010042493924E0 /
  DATA BI1CS( 8) /    .000000000099322077E0 /
  DATA BI1CS( 9) /    .000000000000766380E0 /
  DATA BI1CS(10) /    .000000000000004741E0 /
  DATA BI1CS(11) /    .000000000000000024E0 /
  DATA AI1CS( 1) /   -.02846744181881479E0 /
  DATA AI1CS( 2) /   -.01922953231443221E0 /
  DATA AI1CS( 3) /   -.00061151858579437E0 /
  DATA AI1CS( 4) /   -.00002069971253350E0 /
  DATA AI1CS( 5) /    .00000858561914581E0 /
  DATA AI1CS( 6) /    .00000104949824671E0 /
  DATA AI1CS( 7) /   -.00000029183389184E0 /
  DATA AI1CS( 8) /   -.00000001559378146E0 /
  DATA AI1CS( 9) /    .00000001318012367E0 /
  DATA AI1CS(10) /   -.00000000144842341E0 /
  DATA AI1CS(11) /   -.00000000029085122E0 /
  DATA AI1CS(12) /    .00000000012663889E0 /
  DATA AI1CS(13) /   -.00000000001664947E0 /
  DATA AI1CS(14) /   -.00000000000166665E0 /
  DATA AI1CS(15) /    .00000000000124260E0 /
  DATA AI1CS(16) /   -.00000000000027315E0 /
  DATA AI1CS(17) /    .00000000000002023E0 /
  DATA AI1CS(18) /    .00000000000000730E0 /
  DATA AI1CS(19) /   -.00000000000000333E0 /
  DATA AI1CS(20) /    .00000000000000071E0 /
  DATA AI1CS(21) /   -.00000000000000006E0 /
  DATA AI12CS( 1) /    .02857623501828014E0 /
  DATA AI12CS( 2) /   -.00976109749136147E0 /
  DATA AI12CS( 3) /   -.00011058893876263E0 /
  DATA AI12CS( 4) /   -.00000388256480887E0 /
  DATA AI12CS( 5) /   -.00000025122362377E0 /
  DATA AI12CS( 6) /   -.00000002631468847E0 /
  DATA AI12CS( 7) /   -.00000000383538039E0 /
  DATA AI12CS( 8) /   -.00000000055897433E0 /
  DATA AI12CS( 9) /   -.00000000001897495E0 /
  DATA AI12CS(10) /    .00000000003252602E0 /
  DATA AI12CS(11) /    .00000000001412580E0 /
  DATA AI12CS(12) /    .00000000000203564E0 /
  DATA AI12CS(13) /   -.00000000000071985E0 /
  DATA AI12CS(14) /   -.00000000000040836E0 /
  DATA AI12CS(15) /   -.00000000000002101E0 /
  DATA AI12CS(16) /    .00000000000004273E0 /
  DATA AI12CS(17) /    .00000000000001041E0 /
  DATA AI12CS(18) /   -.00000000000000382E0 /
  DATA AI12CS(19) /   -.00000000000000186E0 /
  DATA AI12CS(20) /    .00000000000000033E0 /
  DATA AI12CS(21) /    .00000000000000028E0 /
  DATA AI12CS(22) /   -.00000000000000003E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BESI1E
  if (FIRST) THEN
     NTI1 = INITS (BI1CS, 11, 0.1*R1MACH(3))
     NTAI1 = INITS (AI1CS, 21, 0.1*R1MACH(3))
     NTAI12 = INITS (AI12CS, 22, 0.1*R1MACH(3))
!
     XMIN = 2.0*R1MACH(1)
     XSML = SQRT (4.5*R1MACH(3))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 3.0) go to 20
!
  BESI1E = 0.0
  if (Y == 0.0)  return
!
  if (Y  <=  XMIN) call XERMSG ('SLATEC', 'BESI1E', &
     'ABS(X) SO SMALL I1 UNDERFLOWS', 1, 1)
  if (Y > XMIN) BESI1E = 0.5*X
  if (Y > XSML) BESI1E = X * (.875 + CSEVL(Y*Y/4.5-1., BI1CS,NTI1))
  BESI1E = EXP(-Y) * BESI1E
  return
!
 20   if (Y <= 8.) BESI1E = (.375 + CSEVL ((48./Y-11.)/5., AI1CS, NTAI1) &
    ) / SQRT(Y)
  if (Y > 8.) BESI1E = (.375 + CSEVL (16./Y-1.0, AI12CS, NTAI12)) &
    / SQRT(Y)
  BESI1E = SIGN (BESI1E, X)
!
  return
end
