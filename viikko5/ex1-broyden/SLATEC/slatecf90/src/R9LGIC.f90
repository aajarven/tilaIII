function R9LGIC (A, X, ALX)
!
!! R9LGIC computes the log complementary incomplete Gamma function ...
!            for large X and for A  <=  X.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      SINGLE PRECISION (R9LGIC-S, D9LGIC-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, LARGE X,
!             LOGARITHM, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the log complementary incomplete gamma function for large X
! and for A  <=  X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  R9LGIC
  SAVE EPS
  DATA EPS / 0.0 /
!***FIRST EXECUTABLE STATEMENT  R9LGIC
  if (EPS == 0.0) EPS = 0.5*R1MACH(3)
!
  XPA = X + 1.0 - A
  XMA = X - 1.0 - A
!
  R = 0.0
  P = 1.0
  S = P
  DO 10 K=1,200
    FK = K
    T = FK*(A-FK)*(1.0+R)
    R = -T/((XMA+2.0*FK)*(XPA+2.0*FK)+T)
    P = R*P
    S = S + P
    if (ABS(P) < EPS*S) go to 20
 10   CONTINUE
  call XERMSG ('SLATEC', 'R9LGIC', &
     'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 1, 2)
!
 20   R9LGIC = A*ALX - X + LOG(S/XPA)
!
  return
end
