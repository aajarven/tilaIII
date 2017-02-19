function R9GMIT (A, X, ALGAP1, SGNGAM, ALX)
!
!! R9GMIT computes Tricomi's incomplete Gamma function for small arguments.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      SINGLE PRECISION (R9GMIT-S, D9GMIT-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
!             SPECIAL FUNCTIONS, TRICOMI
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute Tricomi's incomplete gamma function for small X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  ALNGAM, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  R9GMIT
  SAVE EPS, BOT
  DATA EPS, BOT / 2*0.0 /
!***FIRST EXECUTABLE STATEMENT  R9GMIT
  if (EPS == 0.0) EPS = 0.5*R1MACH(3)
  if (BOT == 0.0) BOT = LOG(R1MACH(1))
!
  if (X  <=  0.0) call XERMSG ('SLATEC', 'R9GMIT', &
     'X SHOULD BE GT 0', 1, 2)
!
  MA = A + 0.5
  if (A < 0.0) MA = A - 0.5
  AEPS = A - MA
!
  AE = A
  if (A < (-0.5)) AE = AEPS
!
  T = 1.0
  TE = AE
  S = T
  DO 20 K=1,200
    FK = K
    TE = -X*TE/FK
    T = TE/(AE+FK)
    S = S + T
    if (ABS(T) < EPS*ABS(S)) go to 30
 20   CONTINUE
  call XERMSG ('SLATEC', 'R9GMIT', &
     'NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SERIES', 2, 2)
!
 30   if (A >= (-0.5)) ALGS = -ALGAP1 + LOG(S)
  if (A >= (-0.5)) go to 60
!
  ALGS = -ALNGAM(1.0+AEPS) + LOG(S)
  S = 1.0
  M = -MA - 1
  if (M == 0) go to 50
  T = 1.0
  DO 40 K=1,M
    T = X*T/(AEPS-M-1+K)
    S = S + T
    if (ABS(T) < EPS*ABS(S)) go to 50
 40   CONTINUE
!
 50   R9GMIT = 0.0
  ALGS = -MA*LOG(X) + ALGS
  if (S == 0.0 .OR. AEPS == 0.0) go to 60
!
  SGNG2 = SGNGAM*SIGN(1.0,S)
  ALG2 = -X - ALGAP1 + LOG(ABS(S))
!
  if (ALG2 > BOT) R9GMIT = SGNG2*EXP(ALG2)
  if (ALGS > BOT) R9GMIT = R9GMIT + EXP(ALGS)
  return
!
 60   R9GMIT = EXP(ALGS)
  return
!
end
