function BESI0E (X)
!
!! BESI0E computes the exponentially scaled modified (hyperbolic) ...
!  Bessel function of the first kind of order zero.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      SINGLE PRECISION (BESI0E-S, DBSI0E-D)
!***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
!             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
!             ORDER ZERO, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BESI0E(X) calculates the exponentially scaled modified (hyperbolic)
! Bessel function of the first kind of order zero for real argument X;
! i.e., EXP(-ABS(X))*I0(X).
!
!
! Series for BI0        on the interval  0.          to  9.00000D+00
!                                        with weighted error   2.46E-18
!                                         log weighted error  17.61
!                               significant figures required  17.90
!                                    decimal places required  18.15
!
!
! Series for AI0        on the interval  1.25000D-01 to  3.33333D-01
!                                        with weighted error   7.87E-17
!                                         log weighted error  16.10
!                               significant figures required  14.69
!                                    decimal places required  16.76
!
!
! Series for AI02       on the interval  0.          to  1.25000D-01
!                                        with weighted error   3.79E-17
!                                         log weighted error  16.42
!                               significant figures required  14.86
!                                    decimal places required  17.09
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, INITS, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890313  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  BESI0E
  DIMENSION BI0CS(12), AI0CS(21), AI02CS(22)
  LOGICAL FIRST
  SAVE BI0CS, AI0CS, AI02CS, NTI0, NTAI0, NTAI02, XSML, FIRST
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
  DATA AI0CS( 1) /    .07575994494023796E0 /
  DATA AI0CS( 2) /    .00759138081082334E0 /
  DATA AI0CS( 3) /    .00041531313389237E0 /
  DATA AI0CS( 4) /    .00001070076463439E0 /
  DATA AI0CS( 5) /   -.00000790117997921E0 /
  DATA AI0CS( 6) /   -.00000078261435014E0 /
  DATA AI0CS( 7) /    .00000027838499429E0 /
  DATA AI0CS( 8) /    .00000000825247260E0 /
  DATA AI0CS( 9) /   -.00000001204463945E0 /
  DATA AI0CS(10) /    .00000000155964859E0 /
  DATA AI0CS(11) /    .00000000022925563E0 /
  DATA AI0CS(12) /   -.00000000011916228E0 /
  DATA AI0CS(13) /    .00000000001757854E0 /
  DATA AI0CS(14) /    .00000000000112822E0 /
  DATA AI0CS(15) /   -.00000000000114684E0 /
  DATA AI0CS(16) /    .00000000000027155E0 /
  DATA AI0CS(17) /   -.00000000000002415E0 /
  DATA AI0CS(18) /   -.00000000000000608E0 /
  DATA AI0CS(19) /    .00000000000000314E0 /
  DATA AI0CS(20) /   -.00000000000000071E0 /
  DATA AI0CS(21) /    .00000000000000007E0 /
  DATA AI02CS( 1) /    .05449041101410882E0 /
  DATA AI02CS( 2) /    .00336911647825569E0 /
  DATA AI02CS( 3) /    .00006889758346918E0 /
  DATA AI02CS( 4) /    .00000289137052082E0 /
  DATA AI02CS( 5) /    .00000020489185893E0 /
  DATA AI02CS( 6) /    .00000002266668991E0 /
  DATA AI02CS( 7) /    .00000000339623203E0 /
  DATA AI02CS( 8) /    .00000000049406022E0 /
  DATA AI02CS( 9) /    .00000000001188914E0 /
  DATA AI02CS(10) /   -.00000000003149915E0 /
  DATA AI02CS(11) /   -.00000000001321580E0 /
  DATA AI02CS(12) /   -.00000000000179419E0 /
  DATA AI02CS(13) /    .00000000000071801E0 /
  DATA AI02CS(14) /    .00000000000038529E0 /
  DATA AI02CS(15) /    .00000000000001539E0 /
  DATA AI02CS(16) /   -.00000000000004151E0 /
  DATA AI02CS(17) /   -.00000000000000954E0 /
  DATA AI02CS(18) /    .00000000000000382E0 /
  DATA AI02CS(19) /    .00000000000000176E0 /
  DATA AI02CS(20) /   -.00000000000000034E0 /
  DATA AI02CS(21) /   -.00000000000000027E0 /
  DATA AI02CS(22) /    .00000000000000003E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BESI0E
  if (FIRST) THEN
     NTI0 = INITS (BI0CS, 12, 0.1*R1MACH(3))
     NTAI0 = INITS (AI0CS, 21, 0.1*R1MACH(3))
     NTAI02 = INITS (AI02CS, 22, 0.1*R1MACH(3))
     XSML = SQRT (4.5*R1MACH(3))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 3.0) go to 20
!
  BESI0E = 1.0 - X
  if (Y > XSML) BESI0E = EXP(-Y) * ( 2.75 + &
    CSEVL (Y*Y/4.5-1.0, BI0CS, NTI0) )
  return
!
 20   if (Y <= 8.) BESI0E = (.375 + CSEVL ((48./Y-11.)/5., AI0CS, NTAI0) &
    ) / SQRT(Y)
  if (Y > 8.) BESI0E = (.375 + CSEVL (16./Y-1., AI02CS, NTAI02)) &
    / SQRT(Y)
!
  return
end
