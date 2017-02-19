FUNCTION GAMIT (A, X)
!
!! GAMIT calculates Tricomi's form of the incomplete Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      SINGLE PRECISION (GAMIT-S, DGAMIT-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
!             SPECIAL FUNCTIONS, TRICOMI
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!   Evaluate Tricomi's incomplete gamma function defined by
!
!   GAMIT = X**(-A)/GAMMA(A) * integral from 0 to X of EXP(-T) *
!             T**(A-1.)
!
!   for A  >  0.0 and by analytic continuation for A  <=  0.0.
!   GAMMA(X) is the complete gamma function of X.
!
!   GAMIT is evaluated for arbitrary real values of A and for non-
!   negative values of X (even though GAMIT is defined for X  <
!   0.0), except that for X = 0 and A  <=  0.0, GAMIT is infinite,
!   which is a fatal error.
!
!   The function and both arguments are REAL.
!
!   A slight deterioration of 2 or 3 digits accuracy will occur when
!   GAMIT is very large or very small in absolute value, because log-
!   arithmic variables are used.  Also, if the parameter  A  is very
!   close to a negative integer (but not a negative integer), there is
!   a loss of accuracy, which is reported if the result is less than
!   half machine precision.
!
!***REFERENCES  W. Gautschi, A computational procedure for incomplete
!                 gamma functions, ACM Transactions on Mathematical
!                 Software 5, 4 (December 1979), pp. 466-481.
!               W. Gautschi, Incomplete gamma functions, Algorithm 542,
!                 ACM Transactions on Mathematical Software 5, 4
!                 (December 1979), pp. 482-489.
!***ROUTINES CALLED  ALGAMS, ALNGAM, GAMR, R1MACH, R9GMIT, R9LGIC,
!                    R9LGIT, XERCLR, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
!***END PROLOGUE  GAMIT
  LOGICAL FIRST
  REAL GAMIT
  SAVE ALNEPS, SQEPS, BOT, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  GAMIT
  if (FIRST) THEN
     ALNEPS = -LOG(R1MACH(3))
     SQEPS = SQRT(R1MACH(4))
     BOT = LOG(R1MACH(1))
  end if
  FIRST = .FALSE.
!
  if (X  <  0.0) call XERMSG ('SLATEC', 'GAMIT', 'X IS NEGATIVE', &
     2, 2)
!
  if (X /= 0.0) ALX = LOG(X)
  SGA = 1.0
  if (A /= 0.0) SGA = SIGN (1.0, A)
  AINTA = AINT (A+0.5*SGA)
  AEPS = A - AINTA
!
  if (X > 0.0) go to 20
  GAMIT = 0.0
  if (AINTA > 0.0 .OR. AEPS /= 0.0) GAMIT = GAMR(A+1.0)
  return
!
 20   if (X > 1.0) go to 40
  if (A >= (-0.5) .OR. AEPS /= 0.0) call ALGAMS (A+1.0, ALGAP1, &
    SGNGAM)
  GAMIT = R9GMIT (A, X, ALGAP1, SGNGAM, ALX)
  return
!
 40   if (A < X) go to 50
  T = R9LGIT (A, X, ALNGAM(A+1.0))
  if (T < BOT) call XERCLR
  GAMIT = EXP(T)
  return
!
 50   ALNG = R9LGIC (A, X, ALX)
!
! EVALUATE GAMIT IN TERMS OF LOG(GAMIC(A,X))
!
  H = 1.0
  if (AEPS == 0.0 .AND. AINTA <= 0.0) go to 60
  call ALGAMS (A+1.0, ALGAP1, SGNGAM)
  T = LOG(ABS(A)) + ALNG - ALGAP1
  if (T > ALNEPS) go to 70
  if (T > (-ALNEPS)) H = 1.0 - SGA*SGNGAM*EXP(T)
  if (ABS(H) > SQEPS) go to 60
  call XERCLR
  call XERMSG ('SLATEC', 'GAMIT', 'RESULT LT HALF PRECISION', 1, 1)
!
 60   T = -A*ALX + LOG(ABS(H))
  if (T < BOT) call XERCLR
  GAMIT = SIGN (EXP(T), H)
  return
!
 70   T = T - A*ALX
  if (T < BOT) call XERCLR
  GAMIT = -SGA*SGNGAM*EXP(T)
  return
!
end
