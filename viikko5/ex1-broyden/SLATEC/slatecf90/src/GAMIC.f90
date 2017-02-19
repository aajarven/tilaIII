FUNCTION GAMIC (A, X)
!
!! GAMIC calculates the complementary incomplete Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      SINGLE PRECISION (GAMIC-S, DGAMIC-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!   Evaluate the complementary incomplete gamma function
!
!   GAMIC = integral from X to infinity of EXP(-T) * T**(A-1.)  .
!
!   GAMIC is evaluated for arbitrary real values of A and for non-
!   negative values of X (even though GAMIC is defined for X  <
!   0.0), except that for X = 0 and A  <=  0.0, GAMIC is undefined.
!
!   GAMIC, A, and X are REAL.
!
!   A slight deterioration of 2 or 3 digits accuracy will occur when
!   GAMIC is very large or very small in absolute value, because log-
!   arithmic variables are used.  Also, if the parameter A is very close
!   to a negative integer (but not a negative integer), there is a loss
!   of accuracy, which is reported if the result is less than half
!   machine precision.
!
!***REFERENCES  W. Gautschi, A computational procedure for incomplete
!                 gamma functions, ACM Transactions on Mathematical
!                 Software 5, 4 (December 1979), pp. 466-481.
!               W. Gautschi, Incomplete gamma functions, Algorithm 542,
!                 ACM Transactions on Mathematical Software 5, 4
!                 (December 1979), pp. 482-489.
!***ROUTINES CALLED  ALGAMS, ALNGAM, R1MACH, R9GMIC, R9GMIT, R9LGIC,
!                    R9LGIT, XERCLR, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
!***END PROLOGUE  GAMIC
  LOGICAL FIRST
  REAL GAMIC
  SAVE EPS, SQEPS, ALNEPS, BOT, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  GAMIC
  if (FIRST) THEN
     EPS = 0.5*R1MACH(3)
     SQEPS = SQRT(R1MACH(4))
     ALNEPS = -LOG(R1MACH(3))
     BOT = LOG(R1MACH(1))
  end if
  FIRST = .FALSE.
!
  if (X  <  0.0) call XERMSG ('SLATEC', 'GAMIC', 'X IS NEGATIVE', &
     2, 2)
!
  if (X > 0.0) go to 20
  if (A  <=  0.0) call XERMSG ('SLATEC', 'GAMIC', &
     'X = 0 AND A LE 0 SO GAMIC IS UNDEFINED', 3, 2)
!
  GAMIC = EXP (ALNGAM(A+1.0) - LOG(A))
  return
!
 20   ALX = LOG(X)
  SGA = 1.0
  if (A /= 0.0) SGA = SIGN (1.0, A)
  MA = A + 0.5*SGA
  AEPS = A - MA
!
  IZERO = 0
  if (X >= 1.0) go to 60
!
  if (A > 0.5 .OR. ABS(AEPS) > 0.001) go to 50
  FM = -MA
  E = 2.0
  if (FM > 1.0) E = 2.0*(FM+2.0)/(FM*FM-1.0)
  E = E - ALX*X**(-0.001)
  if (E*ABS(AEPS) > EPS) go to 50
!
  GAMIC = R9GMIC (A, X, ALX)
  return
!
 50   call ALGAMS (A+1.0, ALGAP1, SGNGAM)
  GSTAR = R9GMIT (A, X, ALGAP1, SGNGAM, ALX)
  if (GSTAR == 0.0) IZERO = 1
  if (GSTAR /= 0.0) ALNGS = LOG (ABS(GSTAR))
  if (GSTAR /= 0.0) SGNGS = SIGN (1.0, GSTAR)
  go to 70
!
 60   if (A < X) GAMIC = EXP (R9LGIC(A, X, ALX))
  if (A < X) RETURN
!
  SGNGAM = 1.0
  ALGAP1 = ALNGAM (A+1.0)
  SGNGS = 1.0
  ALNGS = R9LGIT (A, X, ALGAP1)
!
! EVALUATION OF GAMIC(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
!
 70   H = 1.0
  if (IZERO == 1) go to 80
!
  T = A*ALX + ALNGS
  if (T > ALNEPS) go to 90
  if (T > (-ALNEPS)) H = 1.0 - SGNGS*EXP(T)
!
  if (ABS(H) < SQEPS) call XERCLR
  if (ABS(H)  <  SQEPS) call XERMSG ('SLATEC', 'GAMIC', &
     'RESULT LT HALF PRECISION', 1, 1)
!
 80   SGNG = SIGN (1.0, H) * SGA * SGNGAM
  T = LOG(ABS(H)) + ALGAP1 - LOG(ABS(A))
  if (T < BOT) call XERCLR
  GAMIC = SGNG * EXP(T)
  return
!
 90   SGNG = -SGNGS * SGA * SGNGAM
  T = T + ALGAP1 - LOG(ABS(A))
  if (T < BOT) call XERCLR
  GAMIC = SGNG * EXP(T)
  return
!
end
