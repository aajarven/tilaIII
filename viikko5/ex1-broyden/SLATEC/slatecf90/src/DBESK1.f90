  DOUBLE PRECISION FUNCTION DBESK1 (X)
!
!! DBESK1 computes the modified (hyperbolic) Bessel function of the ...
!            third kind of order one.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      DOUBLE PRECISION (BESK1-S, DBESK1-D)
!***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBESK1(X) calculates the double precision modified (hyperbolic)
! Bessel function of the third kind of order one for double precision
! argument X.  The argument must be large enough that the result does
! not overflow and small enough that the result does not underflow.
!
! Series for BK1        on the interval  0.          to  4.00000E+00
!                                        with weighted error   9.16E-32
!                                         log weighted error  31.04
!                               significant figures required  30.61
!                                    decimal places required  31.64
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DBESI1, DBSK1E, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBESK1
  DOUBLE PRECISION X, BK1CS(16), XMAX, XMAXT, XMIN, XSML, Y, &
    D1MACH, DCSEVL, DBESI1, DBSK1E
  LOGICAL FIRST
  SAVE BK1CS, NTK1, XMIN, XSML, XMAX, FIRST
  DATA BK1CS(  1) / +.25300227338947770532531120868533D-1     /
  DATA BK1CS(  2) / -.35315596077654487566723831691801D+0     /
  DATA BK1CS(  3) / -.12261118082265714823479067930042D+0     /
  DATA BK1CS(  4) / -.69757238596398643501812920296083D-2     /
  DATA BK1CS(  5) / -.17302889575130520630176507368979D-3     /
  DATA BK1CS(  6) / -.24334061415659682349600735030164D-5     /
  DATA BK1CS(  7) / -.22133876307347258558315252545126D-7     /
  DATA BK1CS(  8) / -.14114883926335277610958330212608D-9     /
  DATA BK1CS(  9) / -.66669016941993290060853751264373D-12    /
  DATA BK1CS( 10) / -.24274498505193659339263196864853D-14    /
  DATA BK1CS( 11) / -.70238634793862875971783797120000D-17    /
  DATA BK1CS( 12) / -.16543275155100994675491029333333D-19    /
  DATA BK1CS( 13) / -.32338347459944491991893333333333D-22    /
  DATA BK1CS( 14) / -.53312750529265274999466666666666D-25    /
  DATA BK1CS( 15) / -.75130407162157226666666666666666D-28    /
  DATA BK1CS( 16) / -.91550857176541866666666666666666D-31    /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBESK1
  if (FIRST) THEN
     NTK1 = INITDS (BK1CS, 16, 0.1*REAL(D1MACH(3)))
     XMIN = EXP(MAX(LOG(D1MACH(1)), -LOG(D1MACH(2))) + 0.01D0)
     XSML = SQRT(4.0D0*D1MACH(3))
     XMAXT = -LOG(D1MACH(1))
     XMAX = XMAXT - 0.5D0*XMAXT*LOG(XMAXT)/(XMAXT+0.5D0)
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.D0) call XERMSG ('SLATEC', 'DBESK1', &
     'X IS ZERO OR NEGATIVE', 2, 2)
  if (X > 2.0D0) go to 20
!
  if (X  <  XMIN) call XERMSG ('SLATEC', 'DBESK1', &
     'X SO SMALL K1 OVERFLOWS', 3, 2)
  Y = 0.D0
  if (X > XSML) Y = X*X
  DBESK1 = LOG(0.5D0*X)*DBESI1(X) + (0.75D0 + DCSEVL (.5D0*Y-1.D0, &
    BK1CS, NTK1))/X
  return
!
 20   DBESK1 = 0.D0
  if (X  >  XMAX) call XERMSG ('SLATEC', 'DBESK1', &
     'X SO BIG K1 UNDERFLOWS', 1, 1)
  if (X > XMAX) RETURN
!
  DBESK1 = EXP(-X) * DBSK1E(X)
!
  return
end
