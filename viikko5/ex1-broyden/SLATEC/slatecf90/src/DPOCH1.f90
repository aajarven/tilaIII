  DOUBLE PRECISION FUNCTION DPOCH1 (A, X)
!
!! DPOCH1 calculates a generalization of Pochhammer's symbol starting ...
!  from first order.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C1, C7A
!***TYPE      DOUBLE PRECISION (POCH1-S, DPOCH1-D)
!***KEYWORDS  FIRST ORDER, FNLIB, POCHHAMMER, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate a double precision generalization of Pochhammer's symbol
! for double precision A and X for special situations that require
! especially accurate values when X is small in
!        POCH1(A,X) = (POCH(A,X)-1)/X
!                   = (GAMMA(A+X)/GAMMA(A) - 1.0)/X .
! This specification is particularly suited for stably computing
! expressions such as
!        (GAMMA(A+X)/GAMMA(A) - GAMMA(B+X)/GAMMA(B))/X
!             = POCH1(A,X) - POCH1(B,X)
! Note that POCH1(A,0.0) = PSI(A)
!
! When ABS(X) is so small that substantial cancellation will occur if
! the straightforward formula is used, we use an expansion due
! to Fields and discussed by Y. L. Luke, The Special Functions and Their
! Approximations, Vol. 1, Academic Press, 1969, page 34.
!
! The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as
!        (A+(X-1)/2)**X * polynomial in (A+(X-1)/2)**(-2) .
! In order to maintain significance in POCH1, we write for positive a
!        (A+(X-1)/2)**X = EXP(X*LOG(A+(X-1)/2)) = EXP(Q)
!                       = 1.0 + Q*EXPREL(Q) .
! Likewise the polynomial is written
!        POLY = 1.0 + X*POLY1(A,X) .
! Thus,
!        POCH1(A,X) = (POCH(A,X) - 1) / X
!                   = EXPREL(Q)*(Q/X + Q*POLY1(A,X)) + POLY1(A,X)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCOT, DEXPRL, DPOCH, DPSI, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  DPOCH1
  DOUBLE PRECISION A, X, ABSA, ABSX, ALNEPS, ALNVAR, B, BERN(20), &
    BINV, BP, GBERN(21), GBK, PI, POLY1, Q, RHO, SINPXX, SINPX2, &
    SQTBIG, TERM, TRIG, VAR, VAR2, D1MACH, DPSI, DEXPRL, DCOT, DPOCH
  LOGICAL FIRST
  EXTERNAL DCOT
  SAVE BERN, PI, SQTBIG, ALNEPS, FIRST
  DATA BERN  (  1) / +.833333333333333333333333333333333D-1    /
  DATA BERN  (  2) / -.138888888888888888888888888888888D-2    /
  DATA BERN  (  3) / +.330687830687830687830687830687830D-4    /
  DATA BERN  (  4) / -.826719576719576719576719576719576D-6    /
  DATA BERN  (  5) / +.208767569878680989792100903212014D-7    /
  DATA BERN  (  6) / -.528419013868749318484768220217955D-9    /
  DATA BERN  (  7) / +.133825365306846788328269809751291D-10   /
  DATA BERN  (  8) / -.338968029632258286683019539124944D-12   /
  DATA BERN  (  9) / +.858606205627784456413590545042562D-14   /
  DATA BERN  ( 10) / -.217486869855806187304151642386591D-15   /
  DATA BERN  ( 11) / +.550900282836022951520265260890225D-17   /
  DATA BERN  ( 12) / -.139544646858125233407076862640635D-18   /
  DATA BERN  ( 13) / +.353470703962946747169322997780379D-20   /
  DATA BERN  ( 14) / -.895351742703754685040261131811274D-22   /
  DATA BERN  ( 15) / +.226795245233768306031095073886816D-23   /
  DATA BERN  ( 16) / -.574472439520264523834847971943400D-24   /
  DATA BERN  ( 17) / +.145517247561486490186626486727132D-26   /
  DATA BERN  ( 18) / -.368599494066531017818178247990866D-28   /
  DATA BERN  ( 19) / +.933673425709504467203255515278562D-30   /
  DATA BERN  ( 20) / -.236502241570062993455963519636983D-31   /
  DATA PI / 3.141592653589793238462643383279503D0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DPOCH1
  if (FIRST) THEN
     SQTBIG = 1.0D0/SQRT(24.0D0*D1MACH(1))
     ALNEPS = LOG(D1MACH(3))
  end if
  FIRST = .FALSE.
!
  if (X == 0.0D0) DPOCH1 = DPSI(A)
  if (X == 0.0D0) RETURN
!
  ABSX = ABS(X)
  ABSA = ABS(A)
  if (ABSX > 0.1D0*ABSA) go to 70
  if (ABSX*LOG(MAX(ABSA,2.0D0)) > 0.1D0) go to 70
!
  BP = A
  if (A < (-0.5D0)) BP = 1.0D0 - A - X
  INCR = 0
  if (BP < 10.0D0) INCR = 11.0D0 - BP
  B = BP + INCR
!
  VAR = B + 0.5D0*(X-1.0D0)
  ALNVAR = LOG(VAR)
  Q = X*ALNVAR
!
  POLY1 = 0.0D0
  if (VAR >= SQTBIG) go to 40
  VAR2 = (1.0D0/VAR)**2
!
  RHO = 0.5D0*(X+1.0D0)
  GBERN(1) = 1.0D0
  GBERN(2) = -RHO/12.0D0
  TERM = VAR2
  POLY1 = GBERN(2)*TERM
!
  NTERMS = -0.5D0*ALNEPS/ALNVAR + 1.0D0
  if (NTERMS  >  20) call XERMSG ('SLATEC', 'DPOCH1', &
     'NTERMS IS TOO BIG, MAYBE D1MACH(3) IS BAD', 1, 2)
  if (NTERMS < 2) go to 40
!
  DO 30 K=2,NTERMS
    GBK = 0.0D0
    DO 20 J=1,K
      NDX = K - J + 1
      GBK = GBK + BERN(NDX)*GBERN(J)
 20     CONTINUE
    GBERN(K+1) = -RHO*GBK/K
!
    TERM = TERM * (2*K-2-X)*(2*K-1-X)*VAR2
    POLY1 = POLY1 + GBERN(K+1)*TERM
 30   CONTINUE
!
 40   POLY1 = (X-1.0D0)*POLY1
  DPOCH1 = DEXPRL(Q)*(ALNVAR+Q*POLY1) + POLY1
!
  if (INCR == 0) go to 60
!
! WE HAVE DPOCH1(B,X), BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION
! TO OBTAIN DPOCH1(BP,X).
!
  DO 50 II=1,INCR
    I = INCR - II
    BINV = 1.0D0/(BP+I)
    DPOCH1 = (DPOCH1 - BINV) / (1.0D0 + X*BINV)
 50   CONTINUE
!
 60   if (BP == A) RETURN
!
! WE HAVE DPOCH1(BP,X), BUT A IS LT -0.5.  WE THEREFORE USE A REFLECTION
! FORMULA TO OBTAIN DPOCH1(A,X).
!
  SINPXX = SIN(PI*X)/X
  SINPX2 = SIN(0.5D0*PI*X)
  TRIG = SINPXX*DCOT(PI*B) - 2.0D0*SINPX2*(SINPX2/X)
!
  DPOCH1 = TRIG + (1.0D0 + X*TRIG)*DPOCH1
  return
!
 70   DPOCH1 = (DPOCH(A,X) - 1.0D0) / X
  return
!
end
