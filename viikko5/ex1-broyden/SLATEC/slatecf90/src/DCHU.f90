FUNCTION DCHU (A, B, X)
!
!! DCHU computes the logarithmic confluent hypergeometric function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C11
!***TYPE      DOUBLE PRECISION (CHU-S, DCHU-D)
!***KEYWORDS  FNLIB, LOGARITHMIC CONFLUENT HYPERGEOMETRIC FUNCTION,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DCHU(A,B,X) calculates the double precision logarithmic confluent
! hypergeometric function U(A,B,X) for double precision arguments
! A, B, and X.
!
! This routine is not valid when 1+A-B is close to zero if X is small.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, D9CHU, DEXPRL, DGAMMA, DGAMR, DPOCH,
!                    DPOCH1, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  DCHU
  DOUBLE PRECISION DCHU
  DOUBLE PRECISION A, B, X, AINTB, ALNX, A0, BEPS, B0, C0, EPS, &
    FACTOR, GAMRI1, GAMRNI, PCH1AI, PCH1I, PI, POCHAI, SUM, T, &
    XEPS1, XI, XI1, XN, XTOEPS,  D1MACH, DPOCH, DGAMMA, DGAMR, &
    DPOCH1, DEXPRL, D9CHU
  EXTERNAL DGAMMA
  SAVE PI, EPS
  DATA PI / 3.141592653589793238462643383279503D0 /
  DATA EPS / 0.0D0 /
!***FIRST EXECUTABLE STATEMENT  DCHU
  if (EPS == 0.0D0) EPS = D1MACH(3)
!
  if (X  ==  0.0D0) call XERMSG ('SLATEC', 'DCHU', &
     'X IS ZERO SO DCHU IS INFINITE', 1, 2)
  if (X  <  0.0D0) call XERMSG ('SLATEC', 'DCHU', &
     'X IS NEGATIVE, USE CCHU', 2, 2)
!
  if (MAX(ABS(A),1.0D0)*MAX(ABS(1.0D0+A-B),1.0D0) <  &
    0.99D0*ABS(X)) go to 120
!
! THE ASCENDING SERIES WILL BE USED, BECAUSE THE DESCENDING RATIONAL
! APPROXIMATION (WHICH IS BASED ON THE ASYMPTOTIC SERIES) IS UNSTABLE.
!
  if (ABS(1.0D0+A-B)  <  SQRT(EPS)) call XERMSG ('SLATEC', 'DCHU', &
     'ALGORITHMIS BAD WHEN 1+A-B IS NEAR ZERO FOR SMALL X', 10, 2)
!
  if (B >= 0.0D0) AINTB = AINT(B+0.5D0)
  if (B < 0.0D0) AINTB = AINT(B-0.5D0)
  BEPS = B - AINTB
  N = AINTB
!
  ALNX = LOG(X)
  XTOEPS = EXP (-BEPS*ALNX)
!
! EVALUATE THE FINITE SUM.     -----------------------------------------
!
  if (N >= 1) go to 40
!
! CONSIDER THE CASE B  <  1.0 FIRST.
!
  SUM = 1.0D0
  if (N == 0) go to 30
!
  T = 1.0D0
  M = -N
  DO 20 I=1,M
    XI1 = I - 1
    T = T*(A+XI1)*X/((B+XI1)*(XI1+1.0D0))
    SUM = SUM + T
 20   CONTINUE
!
 30   SUM = DPOCH(1.0D0+A-B, -A)*SUM
  go to 70
!
! NOW CONSIDER THE CASE B  >=  1.0.
!
 40   SUM = 0.0D0
  M = N - 2
  if (M < 0) go to 70
  T = 1.0D0
  SUM = 1.0D0
  if (M == 0) go to 60
!
  DO 50 I=1,M
    XI = I
    T = T * (A-B+XI)*X/((1.0D0-B+XI)*XI)
    SUM = SUM + T
 50   CONTINUE
!
 60   SUM = DGAMMA(B-1.0D0) * DGAMR(A) * X**(1-N) * XTOEPS * SUM
!
! NEXT EVALUATE THE INFINITE SUM.     ----------------------------------
!
 70   ISTRT = 0
  if (N < 1) ISTRT = 1 - N
  XI = ISTRT
!
  FACTOR = (-1.0D0)**N * DGAMR(1.0D0+A-B) * X**ISTRT
  if (BEPS /= 0.0D0) FACTOR = FACTOR * BEPS*PI/SIN(BEPS*PI)
!
  POCHAI = DPOCH (A, XI)
  GAMRI1 = DGAMR (XI+1.0D0)
  GAMRNI = DGAMR (AINTB+XI)
  B0 = FACTOR * DPOCH(A,XI-BEPS) * GAMRNI * DGAMR(XI+1.0D0-BEPS)
!
  if (ABS(XTOEPS-1.0D0) > 0.5D0) go to 90
!
! X**(-BEPS) IS CLOSE TO 1.0D0, SO WE MUST BE CAREFUL IN EVALUATING THE
! DIFFERENCES.
!
  PCH1AI = DPOCH1 (A+XI, -BEPS)
  PCH1I = DPOCH1 (XI+1.0D0-BEPS, BEPS)
  C0 = FACTOR * POCHAI * GAMRNI * GAMRI1 * ( &
    -DPOCH1(B+XI,-BEPS) + PCH1AI - PCH1I + BEPS*PCH1AI*PCH1I)
!
! XEPS1 = (1.0 - X**(-BEPS))/BEPS = (X**(-BEPS) - 1.0)/(-BEPS)
  XEPS1 = ALNX*DEXPRL(-BEPS*ALNX)
!
  DCHU = SUM + C0 + XEPS1*B0
  XN = N
  DO 80 I=1,1000
    XI = ISTRT + I
    XI1 = ISTRT + I - 1
    B0 = (A+XI1-BEPS)*B0*X/((XN+XI1)*(XI-BEPS))
    C0 = (A+XI1)*C0*X/((B+XI1)*XI) &
      - ((A-1.0D0)*(XN+2.D0*XI-1.0D0) + XI*(XI-BEPS)) * B0 &
      / (XI*(B+XI1)*(A+XI1-BEPS))
    T = C0 + XEPS1*B0
    DCHU = DCHU + T
    if (ABS(T) < EPS*ABS(DCHU)) go to 130
 80   CONTINUE
  call XERMSG ('SLATEC', 'DCHU', &
     'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES', 3, 2)
!
! X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE STRAIGHTFORWARD
! FORMULATION IS STABLE.
!
 90   A0 = FACTOR * POCHAI * DGAMR(B+XI) * GAMRI1 / BEPS
  B0 = XTOEPS * B0 / BEPS
!
  DCHU = SUM + A0 - B0
  DO 100 I=1,1000
    XI = ISTRT + I
    XI1 = ISTRT + I - 1
    A0 = (A+XI1)*A0*X/((B+XI1)*XI)
    B0 = (A+XI1-BEPS)*B0*X/((AINTB+XI1)*(XI-BEPS))
    T = A0 - B0
    DCHU = DCHU + T
    if (ABS(T) < EPS*ABS(DCHU)) go to 130
 100  CONTINUE
  call XERMSG ('SLATEC', 'DCHU', &
     'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES', 3, 2)
!
! USE LUKE-S RATIONAL APPROXIMATION IN THE ASYMPTOTIC REGION.
!
 120  DCHU = X**(-A) * D9CHU(A,B,X)
!
 130  return
end
