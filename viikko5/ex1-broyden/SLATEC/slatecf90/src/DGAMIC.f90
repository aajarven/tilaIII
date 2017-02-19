  DOUBLE PRECISION FUNCTION DGAMIC (A, X)
!
!! DGAMIC calculates the complementary incomplete Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (GAMIC-S, DGAMIC-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!   Evaluate the complementary incomplete Gamma function
!
!   DGAMIC = integral from X to infinity of EXP(-T) * T**(A-1.)  .
!
!   DGAMIC is evaluated for arbitrary real values of A and for non-
!   negative values of X (even though DGAMIC is defined for X  <
!   0.0), except that for X = 0 and A  <=  0.0, DGAMIC is undefined.
!
!   DGAMIC, A, and X are DOUBLE PRECISION.
!
!   A slight deterioration of 2 or 3 digits accuracy will occur when
!   DGAMIC is very large or very small in absolute value, because log-
!   arithmic variables are used.  Also, if the parameter A is very close
!   to a negative INTEGER (but not a negative integer), there is a loss
!   of accuracy, which is reported if the result is less than half
!   machine precision.
!
!***REFERENCES  W. Gautschi, A computational procedure for incomplete
!                 gamma functions, ACM Transactions on Mathematical
!                 Software 5, 4 (December 1979), pp. 466-481.
!               W. Gautschi, Incomplete gamma functions, Algorithm 542,
!                 ACM Transactions on Mathematical Software 5, 4
!                 (December 1979), pp. 482-489.
!***ROUTINES CALLED  D1MACH, D9GMIC, D9GMIT, D9LGIC, D9LGIT, DLGAMS,
!                    DLNGAM, XERCLR, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
!***END PROLOGUE  DGAMIC
  DOUBLE PRECISION A, X, AEPS, AINTA, ALGAP1, ALNEPS, ALNGS, ALX, &
    BOT, E, EPS, GSTAR, H, SGA, SGNG, SGNGAM, SGNGS, SQEPS, T, &
    D1MACH, DLNGAM, D9GMIC, D9GMIT, D9LGIC, D9LGIT
  LOGICAL FIRST
  SAVE EPS, SQEPS, ALNEPS, BOT, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DGAMIC
  if (FIRST) THEN
     EPS = 0.5D0*D1MACH(3)
     SQEPS = SQRT(D1MACH(4))
     ALNEPS = -LOG (D1MACH(3))
     BOT = LOG (D1MACH(1))
  end if
  FIRST = .FALSE.
!
  if (X  <  0.D0) call XERMSG ('SLATEC', 'DGAMIC', 'X IS NEGATIVE' &
     , 2, 2)
!
  if (X > 0.D0) go to 20
  if (A  <=  0.D0) call XERMSG ('SLATEC', 'DGAMIC', &
     'X = 0 AND A LE 0 SO DGAMIC IS UNDEFINED', 3, 2)
!
  DGAMIC = EXP (DLNGAM(A+1.D0) - LOG(A))
  return
!
 20   ALX = LOG (X)
  SGA = 1.0D0
  if (A /= 0.D0) SGA = SIGN (1.0D0, A)
  AINTA = AINT (A + 0.5D0*SGA)
  AEPS = A - AINTA
!
  IZERO = 0
  if (X >= 1.0D0) go to 40
!
  if (A > 0.5D0 .OR. ABS(AEPS) > 0.001D0) go to 30
  E = 2.0D0
  if (-AINTA > 1.D0) E = 2.D0*(-AINTA+2.D0)/(AINTA*AINTA-1.0D0)
  E = E - ALX * X**(-0.001D0)
  if (E*ABS(AEPS) > EPS) go to 30
!
  DGAMIC = D9GMIC (A, X, ALX)
  return
!
 30   call DLGAMS (A+1.0D0, ALGAP1, SGNGAM)
  GSTAR = D9GMIT (A, X, ALGAP1, SGNGAM, ALX)
  if (GSTAR == 0.D0) IZERO = 1
  if (GSTAR /= 0.D0) ALNGS = LOG (ABS(GSTAR))
  if (GSTAR /= 0.D0) SGNGS = SIGN (1.0D0, GSTAR)
  go to 50
!
 40   if (A < X) DGAMIC = EXP (D9LGIC(A, X, ALX))
  if (A < X) RETURN
!
  SGNGAM = 1.0D0
  ALGAP1 = DLNGAM (A+1.0D0)
  SGNGS = 1.0D0
  ALNGS = D9LGIT (A, X, ALGAP1)
!
! EVALUATION OF DGAMIC(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
!
 50   H = 1.D0
  if (IZERO == 1) go to 60
!
  T = A*ALX + ALNGS
  if (T > ALNEPS) go to 70
  if (T > (-ALNEPS)) H = 1.0D0 - SGNGS*EXP(T)
!
  if (ABS(H) < SQEPS) call XERCLR
  if (ABS(H)  <  SQEPS) call XERMSG ('SLATEC', 'DGAMIC', &
     'RESULT LT HALF PRECISION', 1, 1)
!
 60   SGNG = SIGN (1.0D0, H) * SGA * SGNGAM
  T = LOG(ABS(H)) + ALGAP1 - LOG(ABS(A))
  if (T < BOT) call XERCLR
  DGAMIC = SGNG * EXP(T)
  return
!
 70   SGNG = -SGNGS * SGA * SGNGAM
  T = T + ALGAP1 - LOG(ABS(A))
  if (T < BOT) call XERCLR
  DGAMIC = SGNG * EXP(T)
  return
!
end
