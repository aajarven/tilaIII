function POCH1 (A, X)
!
!! POCH1 calculates a generalization of Pochhammer's symbol starting
!  from first order.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C1, C7A
!***TYPE      SINGLE PRECISION (POCH1-S, DPOCH1-D)
!***KEYWORDS  FIRST ORDER, FNLIB, POCHHAMMER, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate a generalization of Pochhammer's symbol for special
! situations that require especially accurate values when X is small in
!        POCH1(A,X) = (POCH(A,X)-1)/X
!                   = (GAMMA(A+X)/GAMMA(A) - 1.0)/X .
! This specification is particularly suited for stably computing
! expressions such as
!        (GAMMA(A+X)/GAMMA(A) - GAMMA(B+X)/GAMMA(B))/X
!             = POCH1(A,X) - POCH1(B,X)
! Note that POCH1(A,0.0) = PSI(A)
!
! When ABS(X) is so small that substantial cancellation will occur if
! the straightforward formula is used, we  use an expansion due
! to Fields and discussed by Y. L. Luke, The Special Functions and Their
! Approximations, Vol. 1, Academic Press, 1969, page 34.
!
! The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as
!        (A+(X-1)/2)**X * polynomial in (A+(X-1)/2)**(-2) .
! In order to maintain significance in POCH1, we write for positive A
!        (A+(X-1)/2)**X = EXP(X*LOG(A+(X-1)/2)) = EXP(Q)
!                       = 1.0 + Q*EXPREL(Q) .
! Likewise the polynomial is written
!        POLY = 1.0 + X*POLY1(A,X) .
! Thus,
!        POCH1(A,X) = (POCH(A,X) - 1) / X
!                   = EXPREL(Q)*(Q/X + Q*POLY1(A,X)) + POLY1(A,X)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  COT, EXPREL, POCH, PSI, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  POCH1
  DIMENSION BERN(9), GBERN(10)
  LOGICAL FIRST
  EXTERNAL COT
  SAVE BERN, PI, SQTBIG, ALNEPS, FIRST
  DATA BERN( 1) /   .83333333333333333E-01 /
  DATA BERN( 2) /  -.13888888888888889E-02 /
  DATA BERN( 3) /   .33068783068783069E-04 /
  DATA BERN( 4) /  -.82671957671957672E-06 /
  DATA BERN( 5) /   .20876756987868099E-07 /
  DATA BERN( 6) /  -.52841901386874932E-09 /
  DATA BERN( 7) /   .13382536530684679E-10 /
  DATA BERN( 8) /  -.33896802963225829E-12 /
  DATA BERN( 9) /   .85860620562778446E-14 /
  DATA PI / 3.14159265358979324E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  POCH1
  if (FIRST) THEN
     SQTBIG = 1.0/SQRT(24.0*R1MACH(1))
     ALNEPS = LOG(R1MACH(3))
  end if
  FIRST = .FALSE.
!
  if (X == 0.0) POCH1 = PSI(A)
  if (X == 0.0) RETURN
!
  ABSX = ABS(X)
  ABSA = ABS(A)
  if (ABSX > 0.1*ABSA) go to 70
  if (ABSX*LOG(MAX(ABSA,2.0)) > 0.1) go to 70
!
  BP = A
  if (A < (-0.5)) BP = 1.0 - A - X
  INCR = 0
  if (BP < 10.0) INCR = 11.0 - BP
  B = BP + INCR
!
  VAR = B + 0.5*(X-1.0)
  ALNVAR = LOG(VAR)
  Q = X*ALNVAR
!
  POLY1 = 0.0
  if (VAR >= SQTBIG) go to 40
  VAR2 = (1.0/VAR)**2
!
  RHO = 0.5*(X+1.0)
  GBERN(1) = 1.0
  GBERN(2) = -RHO/12.0
  TERM = VAR2
  POLY1 = GBERN(2)*TERM
!
  NTERMS = -0.5*ALNEPS/ALNVAR + 1.0
  if (NTERMS  >  9) call XERMSG ('SLATEC', 'POCH1', &
     'NTERMS IS TOO BIG, MAYBE R1MACH(3) IS BAD', 1, 2)
  if (NTERMS < 2) go to 40
!
  DO 30 K=2,NTERMS
    GBK = 0.0
    DO 20 J=1,K
      NDX = K - J + 1
      GBK = GBK + BERN(NDX)*GBERN(J)
 20     CONTINUE
    GBERN(K+1) = -RHO*GBK/K
!
    TERM = TERM * (2*K-2.-X)*(2*K-1.-X)*VAR2
    POLY1 = POLY1 + GBERN(K+1)*TERM
 30   CONTINUE
!
 40   POLY1 = (X-1.0)*POLY1
  POCH1 = EXPREL(Q)*(ALNVAR + Q*POLY1) + POLY1
!
  if (INCR == 0) go to 60
!
! WE HAVE POCH1(B,X).  BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION
! TO OBTAIN POCH1(BP,X).
!
  DO 50 II=1,INCR
    I = INCR - II
    BINV = 1.0/(BP+I)
    POCH1 = (POCH1-BINV)/(1.0+X*BINV)
 50   CONTINUE
!
 60   if (BP == A) RETURN
!
! WE HAVE POCH1(BP,X), BUT A IS LT -0.5.  WE THEREFORE USE A REFLECTION
! FORMULA TO OBTAIN POCH1(A,X).
!
  SINPXX = SIN(PI*X)/X
  SINPX2 = SIN(0.5*PI*X)
  TRIG = SINPXX*COT(PI*B) - 2.0*SINPX2*(SINPX2/X)
!
  POCH1 = TRIG + (1.0 + X*TRIG) * POCH1
  return
!
 70   POCH1 = (POCH(A,X) - 1.0) / X
  return
!
end
