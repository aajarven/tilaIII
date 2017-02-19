function R9GMIC (A, X, ALX)
!
!! R9GMIC computes the complementary incomplete Gamma function for A ...
!            near a negative integer and for small X.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      SINGLE PRECISION (R9GMIC-S, D9GMIC-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the complementary incomplete gamma function for A near
! a negative integer and for small X.
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
!***END PROLOGUE  R9GMIC
  SAVE EULER, EPS, BOT
  DATA EULER / .5772156649015329E0 /
  DATA EPS, BOT / 2*0.0 /
!***FIRST EXECUTABLE STATEMENT  R9GMIC
  if (EPS == 0.0) EPS = 0.5*R1MACH(3)
  if (BOT == 0.0) BOT = LOG(R1MACH(1))
!
  if (A  >  0.0) call XERMSG ('SLATEC', 'R9GMIC', &
     'A MUST BE NEAR A NEGATIVE INTEGER', 2, 2)
  if (X  <=  0.0) call XERMSG ('SLATEC', 'R9GMIC', &
     'X MUST BE GT ZERO', 3, 2)
!
  MA = A - 0.5
  FM = -MA
  M = -MA
!
  TE = 1.0
  T = 1.0
  S = T
  DO 20 K=1,200
    FKP1 = K + 1
    TE = -X*TE/(FM+FKP1)
    T = TE/FKP1
    S = S + T
    if (ABS(T) < EPS*S) go to 30
 20   CONTINUE
  call XERMSG ('SLATEC', 'R9GMIC', &
     'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 4, 2)
!
 30   R9GMIC = -ALX - EULER + X*S/(FM+1.0)
  if (M == 0) RETURN
!
  if (M == 1) R9GMIC = -R9GMIC - 1.0 + 1.0/X
  if (M == 1) RETURN
!
  TE = FM
  T = 1.0
  S = T
  MM1 = M - 1
  DO 40 K=1,MM1
    FK = K
    TE = -X*TE/FK
    T = TE/(FM-FK)
    S = S + T
    if (ABS(T) < EPS*ABS(S)) go to 50
 40   CONTINUE
!
 50   DO 60 K=1,M
    R9GMIC = R9GMIC + 1.0/K
 60   CONTINUE
!
  SGNG = 1.0
  if (MOD(M,2) == 1) SGNG = -1.0
  ALNG = LOG(R9GMIC) - ALNGAM(FM+1.0)
!
  R9GMIC = 0.0
  if (ALNG > BOT) R9GMIC = SGNG*EXP(ALNG)
  if (S /= 0.0) R9GMIC = R9GMIC + SIGN (EXP(-FM*ALX+LOG(ABS(S)/FM)) &
    , S)
!
  if (R9GMIC  ==  0.0 .AND. S  ==  0.0) call XERMSG ('SLATEC', &
     'R9GMIC', 'RESULT UNDERFLOWS', 1, 1)
  return
!
end
