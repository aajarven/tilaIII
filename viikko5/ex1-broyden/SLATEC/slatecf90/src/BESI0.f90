function BESI0 (X)
!
!! BESI0 computes the hyperbolic Bessel function of the first kind ...
!  of order zero.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      SINGLE PRECISION (BESI0-S, DBESI0-D)
!***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BESI0(X) computes the modified (hyperbolic) Bessel function
! of the first kind of order zero and real argument X.
!
! Series for BI0        on the interval  0.          to  9.00000D+00
!                                        with weighted error   2.46E-18
!                                         log weighted error  17.61
!                               significant figures required  17.90
!                                    decimal places required  18.15
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  BESI0E, CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  BESI0
  DIMENSION BI0CS(12)
  LOGICAL FIRST
  SAVE BI0CS, NTI0, XSML, XMAX, FIRST
  DATA BI0CS( 1) /   -.07660547252839144951E0 /
  DATA BI0CS( 2) /   1.927337953993808270E0 /
  DATA BI0CS( 3) /    .2282644586920301339E0 /
  DATA BI0CS( 4) /    .01304891466707290428E0 /
  DATA BI0CS( 5) /    .00043442709008164874E0 /
  DATA BI0CS( 6) /    .00000942265768600193E0 /
  DATA BI0CS( 7) /    .00000014340062895106E0 /
  DATA BI0CS( 8) /    .00000000161384906966E0 /
  DATA BI0CS( 9) /    .00000000001396650044E0 /
  DATA BI0CS(10) /    .00000000000009579451E0 /
  DATA BI0CS(11) /    .00000000000000053339E0 /
  DATA BI0CS(12) /    .00000000000000000245E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BESI0
  if (FIRST) THEN
     NTI0 = INITS (BI0CS, 12, 0.1*R1MACH(3))
     XSML = SQRT (4.5*R1MACH(3))
     XMAX = LOG (R1MACH(2))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 3.0) go to 20
!
  BESI0 = 1.0
  if (Y > XSML) BESI0 = 2.75 + CSEVL (Y*Y/4.5-1.0, BI0CS, NTI0)
  return
!
 20   if (Y  >  XMAX) call XERMSG ('SLATEC', 'BESI0', &
     'ABS(X) SO BIG I0 OVERFLOWS', 1, 2)
!
  BESI0 = EXP(Y) * BESI0E(X)
!
  return
end
