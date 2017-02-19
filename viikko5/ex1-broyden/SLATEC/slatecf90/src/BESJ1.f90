function BESJ1 (X)
!
!! BESJ1 computes the Bessel function of the first kind of order one.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10A1
!***TYPE      SINGLE PRECISION (BESJ1-S, DBESJ1-D)
!***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ONE,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BESJ1(X) calculates the Bessel function of the first kind of
! order one for real argument X.
!
! Series for BJ1        on the interval  0.          to  1.60000D+01
!                                        with weighted error   4.48E-17
!                                         log weighted error  16.35
!                               significant figures required  15.77
!                                    decimal places required  16.89
!
! Series for BM1        on the interval  0.          to  6.25000D-02
!                                        with weighted error   5.61E-17
!                                         log weighted error  16.25
!                               significant figures required  14.97
!                                    decimal places required  16.91
!
! Series for BTH1       on the interval  0.          to  6.25000D-02
!                                        with weighted error   4.10E-17
!                                         log weighted error  16.39
!                               significant figures required  15.96
!                                    decimal places required  17.08
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   780601  DATE WRITTEN
!   890210  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  BESJ1
  DIMENSION BJ1CS(12), BM1CS(21), BTH1CS(24)
  LOGICAL FIRST
  SAVE BJ1CS, BM1CS, BTH1CS, PI4, NTJ1, NTM1, NTTH1, &
   XSML, XMIN, XMAX, FIRST
  DATA BJ1CS( 1) /   -.11726141513332787E0 /
  DATA BJ1CS( 2) /   -.25361521830790640E0 /
  DATA BJ1CS( 3) /    .050127080984469569E0 /
  DATA BJ1CS( 4) /   -.004631514809625081E0 /
  DATA BJ1CS( 5) /    .000247996229415914E0 /
  DATA BJ1CS( 6) /   -.000008678948686278E0 /
  DATA BJ1CS( 7) /    .000000214293917143E0 /
  DATA BJ1CS( 8) /   -.000000003936093079E0 /
  DATA BJ1CS( 9) /    .000000000055911823E0 /
  DATA BJ1CS(10) /   -.000000000000632761E0 /
  DATA BJ1CS(11) /    .000000000000005840E0 /
  DATA BJ1CS(12) /   -.000000000000000044E0 /
  DATA BM1CS( 1) /    .1047362510931285E0 /
  DATA BM1CS( 2) /    .00442443893702345E0 /
  DATA BM1CS( 3) /   -.00005661639504035E0 /
  DATA BM1CS( 4) /    .00000231349417339E0 /
  DATA BM1CS( 5) /   -.00000017377182007E0 /
  DATA BM1CS( 6) /    .00000001893209930E0 /
  DATA BM1CS( 7) /   -.00000000265416023E0 /
  DATA BM1CS( 8) /    .00000000044740209E0 /
  DATA BM1CS( 9) /   -.00000000008691795E0 /
  DATA BM1CS(10) /    .00000000001891492E0 /
  DATA BM1CS(11) /   -.00000000000451884E0 /
  DATA BM1CS(12) /    .00000000000116765E0 /
  DATA BM1CS(13) /   -.00000000000032265E0 /
  DATA BM1CS(14) /    .00000000000009450E0 /
  DATA BM1CS(15) /   -.00000000000002913E0 /
  DATA BM1CS(16) /    .00000000000000939E0 /
  DATA BM1CS(17) /   -.00000000000000315E0 /
  DATA BM1CS(18) /    .00000000000000109E0 /
  DATA BM1CS(19) /   -.00000000000000039E0 /
  DATA BM1CS(20) /    .00000000000000014E0 /
  DATA BM1CS(21) /   -.00000000000000005E0 /
  DATA BTH1CS( 1) /    .74060141026313850E0 /
  DATA BTH1CS( 2) /   -.004571755659637690E0 /
  DATA BTH1CS( 3) /    .000119818510964326E0 /
  DATA BTH1CS( 4) /   -.000006964561891648E0 /
  DATA BTH1CS( 5) /    .000000655495621447E0 /
  DATA BTH1CS( 6) /   -.000000084066228945E0 /
  DATA BTH1CS( 7) /    .000000013376886564E0 /
  DATA BTH1CS( 8) /   -.000000002499565654E0 /
  DATA BTH1CS( 9) /    .000000000529495100E0 /
  DATA BTH1CS(10) /   -.000000000124135944E0 /
  DATA BTH1CS(11) /    .000000000031656485E0 /
  DATA BTH1CS(12) /   -.000000000008668640E0 /
  DATA BTH1CS(13) /    .000000000002523758E0 /
  DATA BTH1CS(14) /   -.000000000000775085E0 /
  DATA BTH1CS(15) /    .000000000000249527E0 /
  DATA BTH1CS(16) /   -.000000000000083773E0 /
  DATA BTH1CS(17) /    .000000000000029205E0 /
  DATA BTH1CS(18) /   -.000000000000010534E0 /
  DATA BTH1CS(19) /    .000000000000003919E0 /
  DATA BTH1CS(20) /   -.000000000000001500E0 /
  DATA BTH1CS(21) /    .000000000000000589E0 /
  DATA BTH1CS(22) /   -.000000000000000237E0 /
  DATA BTH1CS(23) /    .000000000000000097E0 /
  DATA BTH1CS(24) /   -.000000000000000040E0 /
  DATA PI4 / 0.78539816339744831E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BESJ1
  if (FIRST) THEN
     NTJ1 = INITS (BJ1CS, 12, 0.1*R1MACH(3))
     NTM1 = INITS (BM1CS, 21, 0.1*R1MACH(3))
     NTTH1 = INITS (BTH1CS, 24, 0.1*R1MACH(3))
!
     XSML = SQRT (8.0*R1MACH(3))
     XMIN = 2.0*R1MACH(1)
     XMAX = 1.0/R1MACH(4)
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 4.0) go to 20
!
  BESJ1 = 0.
  if (Y == 0.0) RETURN
  if (Y  <=  XMIN) call XERMSG ('SLATEC', 'BESJ1', &
     'ABS(X) SO SMALL J1 UNDERFLOWS', 1, 1)
  if (Y > XMIN) BESJ1 = 0.5*X
  if (Y > XSML) BESJ1 = X * (.25 + CSEVL(.125*Y*Y-1., BJ1CS, NTJ1))
  return
!
 20   if (Y  >  XMAX) call XERMSG ('SLATEC', 'BESJ1', &
     'NO PRECISION BECAUSE ABS(X) IS TOO BIG', 2, 2)
  Z = 32.0/Y**2 - 1.0
  AMPL = (0.75 + CSEVL (Z, BM1CS, NTM1)) / SQRT(Y)
  THETA = Y - 3.0*PI4 + CSEVL (Z, BTH1CS, NTTH1) / Y
  BESJ1 = SIGN (AMPL, X) * COS (THETA)
!
  return
end
