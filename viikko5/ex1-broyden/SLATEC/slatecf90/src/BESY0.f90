function BESY0 (X)
!
!! BESY0 computes the Bessel function of the second kind of order zero.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10A1
!***TYPE      SINGLE PRECISION (BESY0-S, DBESY0-D)
!***KEYWORDS  BESSEL FUNCTION, FNLIB, ORDER ZERO, SECOND KIND,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BESY0(X) calculates the Bessel function of the second kind
! of order zero for real argument X.
!
! Series for BY0        on the interval  0.          to  1.60000D+01
!                                        with weighted error   1.20E-17
!                                         log weighted error  16.92
!                               significant figures required  16.15
!                                    decimal places required  17.48
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
!***ROUTINES CALLED  BESJ0, CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  BESY0
  DIMENSION BY0CS(13), BM0CS(21), BTH0CS(24)
  LOGICAL FIRST
  SAVE BY0CS, BM0CS, BTH0CS, TWODPI, PI4, &
   NTY0, NTM0, NTTH0, XSML, XMAX, FIRST
  DATA BY0CS( 1) /   -.011277839392865573E0 /
  DATA BY0CS( 2) /   -.12834523756042035E0 /
  DATA BY0CS( 3) /   -.10437884799794249E0 /
  DATA BY0CS( 4) /    .023662749183969695E0 /
  DATA BY0CS( 5) /   -.002090391647700486E0 /
  DATA BY0CS( 6) /    .000103975453939057E0 /
  DATA BY0CS( 7) /   -.000003369747162423E0 /
  DATA BY0CS( 8) /    .000000077293842676E0 /
  DATA BY0CS( 9) /   -.000000001324976772E0 /
  DATA BY0CS(10) /    .000000000017648232E0 /
  DATA BY0CS(11) /   -.000000000000188105E0 /
  DATA BY0CS(12) /    .000000000000001641E0 /
  DATA BY0CS(13) /   -.000000000000000011E0 /
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
  DATA TWODPI / 0.63661977236758134E0 /
  DATA PI4 / 0.78539816339744831E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BESY0
  if (FIRST) THEN
     NTY0 = INITS (BY0CS, 13, 0.1*R1MACH(3))
     NTM0 = INITS (BM0CS, 21, 0.1*R1MACH(3))
     NTTH0 = INITS (BTH0CS, 24, 0.1*R1MACH(3))

     XSML = SQRT (4.0*R1MACH(3))
     XMAX = 1.0/R1MACH(4)
  end if
  FIRST = .FALSE.

  if (X  <=  0.) call XERMSG ('SLATEC', 'BESY0', &
     'X IS ZERO OR NEGATIVE', 1, 2)
  if (X > 4.0) go to 20

  Y = 0.
  if (X > XSML) Y = X*X
  BESY0 = TWODPI*LOG(0.5*X)*BESJ0(X) + .375 + CSEVL (.125*Y-1., &
    BY0CS, NTY0)
  return
!
 20   if (X  >  XMAX) call XERMSG ('SLATEC', 'BESY0', &
     'NO PRECISION BECAUSE X IS BIG', 2, 2)
!
  Z = 32.0/X**2 - 1.0
  AMPL = (0.75 + CSEVL (Z, BM0CS, NTM0)) / SQRT(X)
  THETA = X - PI4 + CSEVL (Z, BTH0CS, NTTH0) / X
  BESY0 = AMPL * SIN (THETA)
!
  return
end
