DOUBLE PRECISION FUNCTION D9GMIC (A, X, ALX)
!
!! D9GMIC computes the complementary incomplete Gamma function ...
!            for A near a negative integer and X small.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (R9GMIC-S, D9GMIC-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the complementary incomplete gamma function for A near
! a negative integer and for small X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DLNGAM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  D9GMIC
  DOUBLE PRECISION A, X, ALX, ALNG, BOT, EPS, EULER, FK, FKP1, FM, &
    S, SGNG, T, TE, D1MACH, DLNGAM
  LOGICAL FIRST
  SAVE EULER, EPS, BOT, FIRST
  DATA EULER / 0.57721566490153286060651209008240D0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  D9GMIC
  if (FIRST) THEN
     EPS = 0.5D0*D1MACH(3)
     BOT = LOG (D1MACH(1))
  end if
  FIRST = .FALSE.
!
  if (A  >  0.D0) call XERMSG ('SLATEC', 'D9GMIC', &
     'A MUST BE NEAR A NEGATIVE INTEGER', 2, 2)
  if (X  <=  0.D0) call XERMSG ('SLATEC', 'D9GMIC', &
     'X MUST BE GT ZERO', 3, 2)
!
  M = -(A - 0.5D0)
  FM = M
!
  TE = 1.0D0
  T = 1.0D0
  S = T
  DO 20 K=1,200
    FKP1 = K + 1
    TE = -X*TE/(FM+FKP1)
    T = TE/FKP1
    S = S + T
    if (ABS(T) < EPS*S) go to 30
 20   CONTINUE
  call XERMSG ('SLATEC', 'D9GMIC', &
     'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 4, 2)
!
 30   D9GMIC = -ALX - EULER + X*S/(FM+1.0D0)
  if (M == 0) RETURN
!
  if (M == 1) D9GMIC = -D9GMIC - 1.D0 + 1.D0/X
  if (M == 1) RETURN
!
  TE = FM
  T = 1.D0
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
    D9GMIC = D9GMIC + 1.0D0/K
 60   CONTINUE
!
  SGNG = 1.0D0
  if (MOD(M,2) == 1) SGNG = -1.0D0
  ALNG = LOG(D9GMIC) - DLNGAM(FM+1.D0)
!
  D9GMIC = 0.D0
  if (ALNG > BOT) D9GMIC = SGNG * EXP(ALNG)
  if (S /= 0.D0) D9GMIC = D9GMIC + &
    SIGN (EXP(-FM*ALX+LOG(ABS(S)/FM)), S)
!
  if (D9GMIC  ==  0.D0 .AND. S  ==  0.D0) call XERMSG ('SLATEC', &
     'D9GMIC', 'RESULT UNDERFLOWS', 1, 1)
  return
!
end
