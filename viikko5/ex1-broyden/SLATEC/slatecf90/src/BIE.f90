function BIE (X)
!
!! BIE calculates the Bairy function for a negative argument and an
!  exponentially scaled Bairy function for a non-negative argument.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10D
!***TYPE      SINGLE PRECISION (BIE-S, DBIE-D)
!***KEYWORDS  BAIRY FUNCTION, EXPONENTIALLY SCALED, FNLIB,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate BI(X) for X  <=  0  and  BI(X)*EXP(ZETA)  where
! ZETA = 2/3 * X**(3/2)  for X  >=  0.0
!
! Series for BIF        on the interval -1.00000D+00 to  1.00000D+00
!                                        with weighted error   1.88E-19
!                                         log weighted error  18.72
!                               significant figures required  17.74
!                                    decimal places required  19.20
!
! Series for BIG        on the interval -1.00000D+00 to  1.00000D+00
!                                        with weighted error   2.61E-17
!                                         log weighted error  16.58
!                               significant figures required  15.17
!                                    decimal places required  17.03
!
! Series for BIF2       on the interval  1.00000D+00 to  8.00000D+00
!                                        with weighted error   1.11E-17
!                                         log weighted error  16.95
!                        approx significant figures required  16.5
!                                    decimal places required  17.45
!
! Series for BIG2       on the interval  1.00000D+00 to  8.00000D+00
!                                        with weighted error   1.19E-18
!                                         log weighted error  17.92
!                        approx significant figures required  17.2
!                                    decimal places required  18.42
!
! Series for BIP        on the interval  1.25000D-01 to  3.53553D-01
!                                        with weighted error   1.91E-17
!                                         log weighted error  16.72
!                               significant figures required  15.35
!                                    decimal places required  17.41
!
! Series for BIP2       on the interval  0.          to  1.25000D-01
!                                        with weighted error   1.05E-18
!                                         log weighted error  17.98
!                               significant figures required  16.74
!                                    decimal places required  18.71
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, INITS, R1MACH, R9AIMP
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890206  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  BIE
  LOGICAL FIRST
  DIMENSION BIFCS(9), BIGCS(8), BIF2CS(10), BIG2CS(10), BIPCS(24), &
    BIP2CS(29)
  SAVE BIFCS, BIGCS, BIF2CS, BIG2CS, BIPCS, BIP2CS, ATR, BTR, &
   NBIF, NBIG, NBIF2, NBIG2, NBIP, NBIP2, X3SML, X32SML, XBIG, FIRST
  DATA BIFCS( 1) /   -.01673021647198664948E0 /
  DATA BIFCS( 2) /    .1025233583424944561E0 /
  DATA BIFCS( 3) /    .00170830925073815165E0 /
  DATA BIFCS( 4) /    .00001186254546774468E0 /
  DATA BIFCS( 5) /    .00000004493290701779E0 /
  DATA BIFCS( 6) /    .00000000010698207143E0 /
  DATA BIFCS( 7) /    .00000000000017480643E0 /
  DATA BIFCS( 8) /    .00000000000000020810E0 /
  DATA BIFCS( 9) /    .00000000000000000018E0 /
  DATA BIGCS( 1) /    .02246622324857452E0 /
  DATA BIGCS( 2) /    .03736477545301955E0 /
  DATA BIGCS( 3) /    .00044476218957212E0 /
  DATA BIGCS( 4) /    .00000247080756363E0 /
  DATA BIGCS( 5) /    .00000000791913533E0 /
  DATA BIGCS( 6) /    .00000000001649807E0 /
  DATA BIGCS( 7) /    .00000000000002411E0 /
  DATA BIGCS( 8) /    .00000000000000002E0 /
  DATA BIF2CS( 1) /   0.09984572693816041E0 /
  DATA BIF2CS( 2) /    .478624977863005538E0 /
  DATA BIF2CS( 3) /    .0251552119604330118E0 /
  DATA BIF2CS( 4) /    .0005820693885232645E0 /
  DATA BIF2CS( 5) /    .0000074997659644377E0 /
  DATA BIF2CS( 6) /    .0000000613460287034E0 /
  DATA BIF2CS( 7) /    .0000000003462753885E0 /
  DATA BIF2CS( 8) /    .0000000000014288910E0 /
  DATA BIF2CS( 9) /    .0000000000000044962E0 /
  DATA BIF2CS(10) /    .0000000000000000111E0 /
  DATA BIG2CS( 1) /    .033305662145514340E0 /
  DATA BIG2CS( 2) /    .161309215123197068E0 /
  DATA BIG2CS( 3) /    .0063190073096134286E0 /
  DATA BIG2CS( 4) /    .0001187904568162517E0 /
  DATA BIG2CS( 5) /    .0000013045345886200E0 /
  DATA BIG2CS( 6) /    .0000000093741259955E0 /
  DATA BIG2CS( 7) /    .0000000000474580188E0 /
  DATA BIG2CS( 8) /    .0000000000001783107E0 /
  DATA BIG2CS( 9) /    .0000000000000005167E0 /
  DATA BIG2CS(10) /    .0000000000000000011E0 /
  DATA BIPCS( 1) /   -.08322047477943447E0 /
  DATA BIPCS( 2) /    .01146118927371174E0 /
  DATA BIPCS( 3) /    .00042896440718911E0 /
  DATA BIPCS( 4) /   -.00014906639379950E0 /
  DATA BIPCS( 5) /   -.00001307659726787E0 /
  DATA BIPCS( 6) /    .00000632759839610E0 /
  DATA BIPCS( 7) /   -.00000042226696982E0 /
  DATA BIPCS( 8) /   -.00000019147186298E0 /
  DATA BIPCS( 9) /    .00000006453106284E0 /
  DATA BIPCS(10) /   -.00000000784485467E0 /
  DATA BIPCS(11) /   -.00000000096077216E0 /
  DATA BIPCS(12) /    .00000000070004713E0 /
  DATA BIPCS(13) /   -.00000000017731789E0 /
  DATA BIPCS(14) /    .00000000002272089E0 /
  DATA BIPCS(15) /    .00000000000165404E0 /
  DATA BIPCS(16) /   -.00000000000185171E0 /
  DATA BIPCS(17) /    .00000000000059576E0 /
  DATA BIPCS(18) /   -.00000000000012194E0 /
  DATA BIPCS(19) /    .00000000000001334E0 /
  DATA BIPCS(20) /    .00000000000000172E0 /
  DATA BIPCS(21) /   -.00000000000000145E0 /
  DATA BIPCS(22) /    .00000000000000049E0 /
  DATA BIPCS(23) /   -.00000000000000011E0 /
  DATA BIPCS(24) /    .00000000000000001E0 /
  DATA BIP2CS( 1) /   -.113596737585988679E0 /
  DATA BIP2CS( 2) /    .0041381473947881595E0 /
  DATA BIP2CS( 3) /    .0001353470622119332E0 /
  DATA BIP2CS( 4) /    .0000104273166530153E0 /
  DATA BIP2CS( 5) /    .0000013474954767849E0 /
  DATA BIP2CS( 6) /    .0000001696537405438E0 /
  DATA BIP2CS( 7) /   -.0000000100965008656E0 /
  DATA BIP2CS( 8) /   -.0000000167291194937E0 /
  DATA BIP2CS( 9) /   -.0000000045815364485E0 /
  DATA BIP2CS(10) /    .0000000003736681366E0 /
  DATA BIP2CS(11) /    .0000000005766930320E0 /
  DATA BIP2CS(12) /    .0000000000621812650E0 /
  DATA BIP2CS(13) /   -.0000000000632941202E0 /
  DATA BIP2CS(14) /   -.0000000000149150479E0 /
  DATA BIP2CS(15) /    .0000000000078896213E0 /
  DATA BIP2CS(16) /    .0000000000024960513E0 /
  DATA BIP2CS(17) /   -.0000000000012130075E0 /
  DATA BIP2CS(18) /   -.0000000000003740493E0 /
  DATA BIP2CS(19) /    .0000000000002237727E0 /
  DATA BIP2CS(20) /    .0000000000000474902E0 /
  DATA BIP2CS(21) /   -.0000000000000452616E0 /
  DATA BIP2CS(22) /   -.0000000000000030172E0 /
  DATA BIP2CS(23) /    .0000000000000091058E0 /
  DATA BIP2CS(24) /   -.0000000000000009814E0 /
  DATA BIP2CS(25) /   -.0000000000000016429E0 /
  DATA BIP2CS(26) /    .0000000000000005533E0 /
  DATA BIP2CS(27) /    .0000000000000002175E0 /
  DATA BIP2CS(28) /   -.0000000000000001737E0 /
  DATA BIP2CS(29) /   -.0000000000000000010E0 /
  DATA ATR / 8.7506905708484345E0 /
  DATA BTR / -2.093836321356054E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BIE
  if (FIRST) THEN
     ETA = 0.1*R1MACH(3)
     NBIF = INITS (BIFCS, 9, ETA)
     NBIG = INITS (BIGCS, 8, ETA)
     NBIF2 = INITS (BIF2CS, 10, ETA)
     NBIG2 = INITS (BIG2CS, 10, ETA)
     NBIP  = INITS (BIPCS , 24, ETA)
     NBIP2 = INITS (BIP2CS, 29, ETA)
!
     X3SML = ETA**0.3333
     X32SML = 1.3104*X3SML**2
     XBIG = R1MACH(2)**0.6666
  end if
  FIRST = .FALSE.
!
  if (X >= (-1.0)) go to 20
  call R9AIMP (X, XM, THETA)
  BIE = XM * SIN(THETA)
  return
!
 20   if (X > 1.0) go to 30
  Z = 0.0
  if (ABS(X) > X3SML) Z = X**3
  BIE = 0.625 + CSEVL (Z, BIFCS, NBIF) + X*(0.4375 + &
    CSEVL (Z, BIGCS, NBIG))
  if (X > X32SML) BIE = BIE * EXP(-2.0*X*SQRT(X)/3.0)
  return
!
 30   if (X > 2.0) go to 40
  Z = (2.0*X**3 - 9.0) / 7.0
  BIE = EXP(-2.0*X*SQRT(X)/3.0) * (1.125 + CSEVL (Z, BIF2CS, NBIF2) &
    + X*(0.625 + CSEVL (Z, BIG2CS, NBIG2)) )
  return
!
 40   if (X > 4.0) go to 50
  SQRTX = SQRT(X)
  Z = ATR/(X*SQRTX) + BTR
  BIE = (0.625 + CSEVL (Z, BIPCS, NBIP)) / SQRT(SQRTX)
  return
!
 50   SQRTX = SQRT(X)
  Z = -1.0
  if (X < XBIG) Z = 16.0/(X*SQRTX) - 1.0
  BIE = (0.625 + CSEVL (Z, BIP2CS, NBIP2))/SQRT(SQRTX)
  return
!
end
