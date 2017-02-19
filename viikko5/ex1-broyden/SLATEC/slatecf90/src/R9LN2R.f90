function R9LN2R (X)
!
!! R9LN2R evaluates LOG(1+X) from second order relative accuracy so ...
!            that LOG(1+X) = X - X**2/2 + X**3*R9LN2R(X).
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4B
!***TYPE      SINGLE PRECISION (R9LN2R-S, D9LN2R-D, C9LN2R-C)
!***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM, SECOND ORDER
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate  LOG(1+X)  from 2-nd order with relative error accuracy so
! that    LOG(1+X) = X - X**2/2 + X**3*R9LN2R(X)
!
! Series for LN21       on the interval -6.25000D-01 to  0.
!                                        with weighted error   2.49E-17
!                                         log weighted error  16.60
!                               significant figures required  15.87
!                                    decimal places required  17.31
!
! Series for LN22       on the interval  0.          to  8.12500D-01
!                                        with weighted error   1.42E-17
!                                         log weighted error  16.85
!                               significant figures required  15.95
!                                    decimal places required  17.50
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   780401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  R9LN2R
  REAL LN21CS(26), LN22CS(20)
  LOGICAL FIRST
  SAVE LN21CS, LN22CS, NTLN21, NTLN22, XMIN, XBIG, XMAX, FIRST
  DATA LN21CS( 1) /    .18111962513478810E0 /
  DATA LN21CS( 2) /   -.15627123192872463E0 /
  DATA LN21CS( 3) /    .028676305361557275E0 /
  DATA LN21CS( 4) /   -.005558699655948139E0 /
  DATA LN21CS( 5) /    .001117897665229983E0 /
  DATA LN21CS( 6) /   -.000230805089823279E0 /
  DATA LN21CS( 7) /    .000048598853341100E0 /
  DATA LN21CS( 8) /   -.000010390127388903E0 /
  DATA LN21CS( 9) /    .000002248456370739E0 /
  DATA LN21CS(10) /   -.000000491405927392E0 /
  DATA LN21CS(11) /    .000000108282565070E0 /
  DATA LN21CS(12) /   -.000000024025872763E0 /
  DATA LN21CS(13) /    .000000005362460047E0 /
  DATA LN21CS(14) /   -.000000001202995136E0 /
  DATA LN21CS(15) /    .000000000271078892E0 /
  DATA LN21CS(16) /   -.000000000061323562E0 /
  DATA LN21CS(17) /    .000000000013920858E0 /
  DATA LN21CS(18) /   -.000000000003169930E0 /
  DATA LN21CS(19) /    .000000000000723837E0 /
  DATA LN21CS(20) /   -.000000000000165700E0 /
  DATA LN21CS(21) /    .000000000000038018E0 /
  DATA LN21CS(22) /   -.000000000000008741E0 /
  DATA LN21CS(23) /    .000000000000002013E0 /
  DATA LN21CS(24) /   -.000000000000000464E0 /
  DATA LN21CS(25) /    .000000000000000107E0 /
  DATA LN21CS(26) /   -.000000000000000024E0 /
  DATA LN22CS( 1) /   -.22242532535020461E0 /
  DATA LN22CS( 2) /   -.061047100108078624E0 /
  DATA LN22CS( 3) /    .007427235009750394E0 /
  DATA LN22CS( 4) /   -.000933501826163697E0 /
  DATA LN22CS( 5) /    .000120049907687260E0 /
  DATA LN22CS( 6) /   -.000015704722952820E0 /
  DATA LN22CS( 7) /    .000002081874781051E0 /
  DATA LN22CS( 8) /   -.000000278919557764E0 /
  DATA LN22CS( 9) /    .000000037693558237E0 /
  DATA LN22CS(10) /   -.000000005130902896E0 /
  DATA LN22CS(11) /    .000000000702714117E0 /
  DATA LN22CS(12) /   -.000000000096748595E0 /
  DATA LN22CS(13) /    .000000000013381046E0 /
  DATA LN22CS(14) /   -.000000000001858102E0 /
  DATA LN22CS(15) /    .000000000000258929E0 /
  DATA LN22CS(16) /   -.000000000000036195E0 /
  DATA LN22CS(17) /    .000000000000005074E0 /
  DATA LN22CS(18) /   -.000000000000000713E0 /
  DATA LN22CS(19) /    .000000000000000100E0 /
  DATA LN22CS(20) /   -.000000000000000014E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  R9LN2R
  if (FIRST) THEN
     EPS = R1MACH(3)
     NTLN21 = INITS (LN21CS, 26, 0.1*EPS)
     NTLN22 = INITS (LN22CS, 20, 0.1*EPS)
!
     XMIN = -1.0 + SQRT(R1MACH(4))
     SQEPS = SQRT(EPS)
     TXMAX = 6.0/SQEPS
     XMAX = TXMAX - (EPS*TXMAX**2 - 2.0*LOG(TXMAX)) / &
                                                (2.0*EPS*TXMAX)
     TXBIG = 4.0/SQRT(SQEPS)
     XBIG = TXBIG - (SQEPS*TXBIG**2 - 2.0*LOG(TXBIG)) / &
                                                  (2.*SQEPS*TXBIG)
  end if
  FIRST = .FALSE.
!
  if (X < (-0.625) .OR. X > 0.8125) go to 20
!
  if (X < 0.0) R9LN2R = 0.375 + CSEVL (16.*X/5.+1.0, LN21CS, &
    NTLN21)
  if (X >= 0.0) R9LN2R = 0.375 + CSEVL (32.*X/13.-1.0, LN22CS, &
    NTLN22)
  return
!
 20   if (X  <  XMIN) call XERMSG ('SLATEC', 'R9LN2R', &
     'ANSWER LT HALF PRECISION BECAUSE X IS TOO NEAR -1', 1, 1)
  if (X  >  XMAX) call XERMSG ('SLATEC', 'R9LN2R', &
     'NO PRECISION IN ANSWER BECAUSE X IS TOO BIG', 3, 2)
  if (X  >  XBIG) call XERMSG ('SLATEC', 'R9LN2R', &
     'ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG', 2, 1)
!
  R9LN2R = (LOG(1.0+X) - X*(1.0-0.5*X) ) / X**3
  return
!
end
