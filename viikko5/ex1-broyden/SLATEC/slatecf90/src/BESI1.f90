function BESI1 (X)
!
!! BESI1 computes the modified (hyperbolic) Bessel function of the ...
!  first kind of order one.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      SINGLE PRECISION (BESI1-S, DBESI1-D)
!***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BESI1(X) calculates the modified (hyperbolic) Bessel function
! of the first kind of order one for real argument X.
!
! Series for BI1        on the interval  0.          to  9.00000D+00
!                                        with weighted error   2.40E-17
!                                         log weighted error  16.62
!                               significant figures required  16.23
!                                    decimal places required  17.14
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  BESI1E, CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  BESI1
  DIMENSION BI1CS(11)
  LOGICAL FIRST
  SAVE BI1CS, NTI1, XMIN, XSML, XMAX, FIRST
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
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BESI1
  if (FIRST) THEN
     NTI1 = INITS (BI1CS, 11, 0.1*R1MACH(3))
     XMIN = 2.0*R1MACH(1)
     XSML = SQRT (4.5*R1MACH(3))
     XMAX = LOG (R1MACH(2))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 3.0) go to 20
!
  BESI1 = 0.0
  if (Y == 0.0)  return
!
  if (Y  <=  XMIN) call XERMSG ('SLATEC', 'BESI1', &
     'ABS(X) SO SMALL I1 UNDERFLOWS', 1, 1)
  if (Y > XMIN) BESI1 = 0.5*X
  if (Y > XSML) BESI1 = X * (.875 + CSEVL(Y*Y/4.5-1., BI1CS, NTI1))
  return
!
 20   if (Y  >  XMAX) call XERMSG ('SLATEC', 'BESI1', &
     'ABS(X) SO BIG I1 OVERFLOWS', 2, 2)
!
  BESI1 = EXP(Y) * BESI1E(X)
!
  return
end
