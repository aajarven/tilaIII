function BI (X)
!
!! BI evaluates the Bairy function (the Airy function of the second kind).
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10D
!***TYPE      SINGLE PRECISION (BI-S, DBI-D)
!***KEYWORDS  BAIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BI(X) calculates the Airy function of the second kind for real
! argument X.
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
!***REFERENCES  (NONE)
!***ROUTINES CALLED  BIE, CSEVL, INITS, R1MACH, R9AIMP, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  BI
  DIMENSION BIFCS(9), BIGCS(8), BIF2CS(10), BIG2CS(10)
  LOGICAL FIRST
  SAVE BIFCS, BIGCS, BIF2CS, BIG2CS, NBIF, NBIG, NBIF2, &
   NBIG2, X3SML, XMAX, FIRST
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
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BI
  if (FIRST) THEN
     ETA = 0.1*R1MACH(3)
     NBIF  = INITS (BIFCS , 9, ETA)
     NBIG  = INITS (BIGCS , 8, ETA)
     NBIF2 = INITS (BIF2CS, 10, ETA)
     NBIG2 = INITS (BIG2CS, 10, ETA)
!
     X3SML = ETA**0.3333
     XMAX = (1.5*LOG(R1MACH(2)))**0.6666
  end if
  FIRST = .FALSE.
!
  if (X >= (-1.0)) go to 20
  call R9AIMP (X, XM, THETA)
  BI = XM * SIN(THETA)
  return
!
 20   if (X > 1.0) go to 30
  Z = 0.0
  if (ABS(X) > X3SML) Z = X**3
  BI = 0.625 + CSEVL (Z, BIFCS, NBIF) + X*(0.4375 + &
    CSEVL (Z, BIGCS, NBIG))
  return
!
 30   if (X > 2.0) go to 40
  Z = (2.0*X**3 - 9.0) / 7.0
  BI = 1.125 + CSEVL (Z, BIF2CS, NBIF2) + X*(0.625 + &
    CSEVL (Z, BIG2CS, NBIG2))
  return
!
 40   if (X  >  XMAX) call XERMSG ('SLATEC', 'BI', &
     'X SO BIG THAT BI OVERFLOWS', 1, 2)
!
  BI = BIE(X) * EXP(2.0*X*SQRT(X)/3.0)
  return
!
end
