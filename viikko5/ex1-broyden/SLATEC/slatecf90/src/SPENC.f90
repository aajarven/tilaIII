function SPENC (X)
!
!! SPENC computes a form of Spence's integral due to K. Mitchell.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C5
!***TYPE      SINGLE PRECISION (SPENC-S, DSPENC-D)
!***KEYWORDS  FNLIB, SPECIAL FUNCTIONS, SPENCE'S INTEGRAL
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate a form of Spence's function defined by
!        integral from 0 to X of  -LOG(1-Y)/Y  DY.
! For ABS(X)  <=  1, the uniformly convergent expansion
!        SPENC = sum K=1,infinity  X**K / K**2     is valid.
!
! Spence's function can be used to evaluate much more general integral
! forms.  For example,
!        integral from 0 to Z of  LOG(A*X+B)/(C*X+D)  DX  =
!             LOG(ABS(B-A*D/C))*LOG(ABS(A*(C*X+D)/(A*D-B*C)))/C
!             - SPENC (A*(C*Z+D)/(A*D-B*C)) / C.
!
! Ref -- K. Mitchell, Philosophical Magazine, 40, p. 351 (1949).
!        Stegun and Abromowitz, AMS 55, p. 1004.
!
!
! Series for SPEN       on the interval  0.          to  5.00000D-01
!                                        with weighted error   6.82E-17
!                                         log weighted error  16.17
!                               significant figures required  15.22
!                                    decimal places required  16.81
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, INITS, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   780201  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  SPENC
  DIMENSION SPENCS(19)
  LOGICAL FIRST
  SAVE SPENCS, PI26, NSPENC, XBIG, FIRST
  DATA SPENCS( 1) /    .1527365598892406E0 /
  DATA SPENCS( 2) /    .08169658058051014E0 /
  DATA SPENCS( 3) /    .00581415714077873E0 /
  DATA SPENCS( 4) /    .00053716198145415E0 /
  DATA SPENCS( 5) /    .00005724704675185E0 /
  DATA SPENCS( 6) /    .00000667454612164E0 /
  DATA SPENCS( 7) /    .00000082764673397E0 /
  DATA SPENCS( 8) /    .00000010733156730E0 /
  DATA SPENCS( 9) /    .00000001440077294E0 /
  DATA SPENCS(10) /    .00000000198444202E0 /
  DATA SPENCS(11) /    .00000000027940058E0 /
  DATA SPENCS(12) /    .00000000004003991E0 /
  DATA SPENCS(13) /    .00000000000582346E0 /
  DATA SPENCS(14) /    .00000000000085767E0 /
  DATA SPENCS(15) /    .00000000000012768E0 /
  DATA SPENCS(16) /    .00000000000001918E0 /
  DATA SPENCS(17) /    .00000000000000290E0 /
  DATA SPENCS(18) /    .00000000000000044E0 /
  DATA SPENCS(19) /    .00000000000000006E0 /
  DATA PI26 / 1.644934066848226E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  SPENC
  if (FIRST) THEN
     NSPENC = INITS (SPENCS, 19, 0.1*R1MACH(3))
     XBIG = 1.0/R1MACH(3)
  end if
  FIRST = .FALSE.
!
  if (X > 2.0) go to 60
  if (X > 1.0) go to 50
  if (X > 0.5) go to 40
  if (X >= 0.0) go to 30
  if (X > (-1.)) go to 20
!
! HERE if X  <=  -1.0
!
  ALN = LOG(1.0-X)
  SPENC = -PI26 - 0.5*ALN*(2.0*LOG(-X)-ALN)
  if (X > (-XBIG)) SPENC = SPENC &
    + (1.0 + CSEVL (4.0/(1.0-X)-1.0, SPENCS, NSPENC)) / (1.0-X)
  return
!
! -1.0  <  X  <  0.0
!
 20   SPENC = -0.5*LOG(1.0-X)**2 &
    - X*(1.0 + CSEVL (4.0*X/(X-1.0)-1.0, SPENCS, NSPENC)) / (X-1.0)
  return
!
! 0.0  <=  X  <=  0.5
!
 30   SPENC = X*(1.0 + CSEVL (4.0*X-1.0, SPENCS, NSPENC))
  return
!
! 0.5  <  X  <=  1.0
!
 40   SPENC = PI26
  if (X /= 1.0) SPENC = PI26 - LOG(X)*LOG(1.0-X) &
    - (1.0-X)*(1.0 + CSEVL (4.0*(1.0-X)-1.0, SPENCS, NSPENC))
  return
!
! 1.0  <  X  <=  2.0
!
 50   SPENC = PI26 - 0.5*LOG(X)*LOG((X-1.0)**2/X) &
    + (X-1.)*(1.0 + CSEVL (4.0*(X-1.)/X-1.0, SPENCS, NSPENC))/X
  return
!
! X  >  2.0
!
 60   SPENC = 2.0*PI26 - 0.5*LOG(X)**2
  if (X < XBIG) SPENC = SPENC &
    - (1.0 + CSEVL (4.0/X-1.0, SPENCS, NSPENC))/X
  return
!
end
