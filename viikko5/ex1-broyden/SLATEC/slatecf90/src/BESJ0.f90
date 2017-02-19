function BESJ0 (X)
!
!! BESJ0 computes the Bessel function of the first kind of order zero.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10A1
!***TYPE      SINGLE PRECISION (BESJ0-S, DBESJ0-D)
!***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ZERO,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BESJ0(X) calculates the Bessel function of the first kind of
! order zero for real argument X.
!
! Series for BJ0        on the interval  0.          to  1.60000D+01
!                                        with weighted error   7.47E-18
!                                         log weighted error  17.13
!                               significant figures required  16.98
!                                    decimal places required  17.68
!
! Series for BM0        on the interval  0.          to  6.25000D-02
!                                        with weighted error   4.98E-17
!                                         log weighted error  16.30
!                               significant figures required  14.97
!                                    decimal places required  16.96
!
! Series for BTH0       on the interval  0.          to  6.25000D-02
!                                        with weighted error   3.67E-17
!                                         log weighted error  16.44
!                               significant figures required  15.53
!                                    decimal places required  17.13
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
!***END PROLOGUE  BESJ0
  DIMENSION BJ0CS(13), BM0CS(21), BTH0CS(24)
  LOGICAL FIRST
  SAVE BJ0CS, BM0CS, BTH0CS, PI4, NTJ0, NTM0, NTTH0, XSML, XMAX, &
     FIRST
  DATA BJ0CS( 1) /    .100254161968939137E0 /
  DATA BJ0CS( 2) /   -.665223007764405132E0 /
  DATA BJ0CS( 3) /    .248983703498281314E0 /
  DATA BJ0CS( 4) /   -.0332527231700357697E0 /
  DATA BJ0CS( 5) /    .0023114179304694015E0 /
  DATA BJ0CS( 6) /   -.0000991127741995080E0 /
  DATA BJ0CS( 7) /    .0000028916708643998E0 /
  DATA BJ0CS( 8) /   -.0000000612108586630E0 /
  DATA BJ0CS( 9) /    .0000000009838650793E0 /
  DATA BJ0CS(10) /   -.0000000000124235515E0 /
  DATA BJ0CS(11) /    .0000000000001265433E0 /
  DATA BJ0CS(12) /   -.0000000000000010619E0 /
  DATA BJ0CS(13) /    .0000000000000000074E0 /
  DATA BM0CS( 1) /    .09284961637381644E0 /
  DATA BM0CS( 2) /   -.00142987707403484E0 /
  DATA BM0CS( 3) /    .00002830579271257E0 /
  DATA BM0CS( 4) /   -.00000143300611424E0 /
  DATA BM0CS( 5) /    .00000012028628046E0 /
  DATA BM0CS( 6) /   -.00000001397113013E0 /
  DATA BM0CS( 7) /    .00000000204076188E0 /
  DATA BM0CS( 8) /   -.00000000035399669E0 /
  DATA BM0CS( 9) /    .00000000007024759E0 /
  DATA BM0CS(10) /   -.00000000001554107E0 /
  DATA BM0CS(11) /    .00000000000376226E0 /
  DATA BM0CS(12) /   -.00000000000098282E0 /
  DATA BM0CS(13) /    .00000000000027408E0 /
  DATA BM0CS(14) /   -.00000000000008091E0 /
  DATA BM0CS(15) /    .00000000000002511E0 /
  DATA BM0CS(16) /   -.00000000000000814E0 /
  DATA BM0CS(17) /    .00000000000000275E0 /
  DATA BM0CS(18) /   -.00000000000000096E0 /
  DATA BM0CS(19) /    .00000000000000034E0 /
  DATA BM0CS(20) /   -.00000000000000012E0 /
  DATA BM0CS(21) /    .00000000000000004E0 /
  DATA BTH0CS( 1) /   -.24639163774300119E0 /
  DATA BTH0CS( 2) /    .001737098307508963E0 /
  DATA BTH0CS( 3) /   -.000062183633402968E0 /
  DATA BTH0CS( 4) /    .000004368050165742E0 /
  DATA BTH0CS( 5) /   -.000000456093019869E0 /
  DATA BTH0CS( 6) /    .000000062197400101E0 /
  DATA BTH0CS( 7) /   -.000000010300442889E0 /
  DATA BTH0CS( 8) /    .000000001979526776E0 /
  DATA BTH0CS( 9) /   -.000000000428198396E0 /
  DATA BTH0CS(10) /    .000000000102035840E0 /
  DATA BTH0CS(11) /   -.000000000026363898E0 /
  DATA BTH0CS(12) /    .000000000007297935E0 /
  DATA BTH0CS(13) /   -.000000000002144188E0 /
  DATA BTH0CS(14) /    .000000000000663693E0 /
  DATA BTH0CS(15) /   -.000000000000215126E0 /
  DATA BTH0CS(16) /    .000000000000072659E0 /
  DATA BTH0CS(17) /   -.000000000000025465E0 /
  DATA BTH0CS(18) /    .000000000000009229E0 /
  DATA BTH0CS(19) /   -.000000000000003448E0 /
  DATA BTH0CS(20) /    .000000000000001325E0 /
  DATA BTH0CS(21) /   -.000000000000000522E0 /
  DATA BTH0CS(22) /    .000000000000000210E0 /
  DATA BTH0CS(23) /   -.000000000000000087E0 /
  DATA BTH0CS(24) /    .000000000000000036E0 /
  DATA PI4 / 0.78539816339744831E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BESJ0
  if (FIRST) THEN
     NTJ0 = INITS (BJ0CS, 13, 0.1*R1MACH(3))
     NTM0 = INITS (BM0CS, 21, 0.1*R1MACH(3))
     NTTH0 = INITS (BTH0CS, 24, 0.1*R1MACH(3))
!
     XSML = SQRT (8.0*R1MACH(3))
     XMAX = 1.0/R1MACH(4)
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 4.0) go to 20
!
  BESJ0 = 1.0
  if (Y > XSML) BESJ0 = CSEVL (.125*Y*Y-1., BJ0CS, NTJ0)
  return
!
 20   if (Y  >  XMAX) call XERMSG ('SLATEC', 'BESJ0', &
     'NO PRECISION BECAUSE ABS(X) IS TOO BIG', 1, 2)
!
  Z = 32.0/Y**2 - 1.0
  AMPL = (0.75 + CSEVL (Z, BM0CS, NTM0)) / SQRT(Y)
  THETA = Y - PI4 + CSEVL (Z, BTH0CS, NTTH0) / Y
  BESJ0 = AMPL * COS (THETA)
!
  return
end
