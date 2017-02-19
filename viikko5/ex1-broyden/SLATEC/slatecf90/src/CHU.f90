function CHU (A, B, X)
!
!! CHU computes the logarithmic confluent hypergeometric function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C11
!***TYPE      SINGLE PRECISION (CHU-S, DCHU-D)
!***KEYWORDS  FNLIB, LOGARITHMIC CONFLUENT HYPERGEOMETRIC FUNCTION,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CHU computes the logarithmic confluent hypergeometric function,
! U(A,B,X).
!
! Input Parameters:
!       A   real
!       B   real
!       X   real and positive
!
! This routine is not valid when 1+A-B is close to zero if X is small.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  EXPREL, GAMMA, GAMR, POCH, POCH1, R1MACH, R9CHU,
!                    XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  CHU
  EXTERNAL GAMMA
  SAVE PI, EPS
  DATA PI / 3.14159265358979324E0 /
  DATA EPS / 0.0 /
!***FIRST EXECUTABLE STATEMENT  CHU
  if (EPS == 0.0) EPS = R1MACH(3)
!
  if (X  ==  0.0) call XERMSG ('SLATEC', 'CHU', &
     'X IS ZERO SO CHU IS INFINITE', 1, 2)
  if (X  <  0.0) call XERMSG ('SLATEC', 'CHU', &
     'X IS NEGATIVE, USE CCHU', 2, 2)
!
  if (MAX(ABS(A),1.0)*MAX(ABS(1.0+A-B),1.0) < 0.99*ABS(X)) &
    go to 120
!
! THE ASCENDING SERIES WILL BE USED, BECAUSE THE DESCENDING RATIONAL
! APPROXIMATION (WHICH IS BASED ON THE ASYMPTOTIC SERIES) IS UNSTABLE.
!
  if (ABS(1.0+A-B)  <  SQRT(EPS)) call XERMSG ('SLATEC', 'CHU', &
     'ALGORITHM IS BAD WHEN 1+A-B IS NEAR ZERO FOR SMALL X', 10, 2)
!
  AINTB = AINT(B+0.5)
  if (B < 0.0) AINTB = AINT(B-0.5)
  BEPS = B - AINTB
  N = AINTB
!
  ALNX = LOG(X)
  XTOEPS = EXP(-BEPS*ALNX)
!
! EVALUATE THE FINITE SUM.     -----------------------------------------
!
  if (N >= 1) go to 40
!
! CONSIDER THE CASE B  <  1.0 FIRST.
!
  SUM = 1.0
  if (N == 0) go to 30
!
  T = 1.0
  M = -N
  DO 20 I=1,M
    XI1 = I - 1
    T = T*(A+XI1)*X/((B+XI1)*(XI1+1.0))
    SUM = SUM + T
 20   CONTINUE
!
 30   SUM = POCH(1.0+A-B, -A) * SUM
  go to 70
!
! NOW CONSIDER THE CASE B  >=  1.0.
!
 40   SUM = 0.0
  M = N - 2
  if (M < 0) go to 70
  T = 1.0
  SUM = 1.0
  if (M == 0) go to 60
!
  DO 50 I=1,M
    XI = I
    T = T * (A-B+XI)*X/((1.0-B+XI)*XI)
    SUM = SUM + T
 50   CONTINUE
!
 60   SUM = GAMMA(B-1.0) * GAMR(A) * X**(1-N) * XTOEPS * SUM
!
! NOW EVALUATE THE INFINITE SUM.     -----------------------------------
!
 70   ISTRT = 0
  if (N < 1) ISTRT = 1 - N
  XI = ISTRT
!
  FACTOR = (-1.0)**N * GAMR(1.0+A-B) * X**ISTRT
  if (BEPS /= 0.0) FACTOR = FACTOR * BEPS*PI/SIN(BEPS*PI)
!
  POCHAI = POCH (A, XI)
  GAMRI1 = GAMR (XI+1.0)
  GAMRNI = GAMR (AINTB+XI)
  B0 = FACTOR * POCH(A,XI-BEPS) * GAMRNI * GAMR(XI+1.0-BEPS)
!
  if (ABS(XTOEPS-1.0) > 0.5) go to 90
!
! X**(-BEPS) IS CLOSE TO 1.0, SO WE MUST BE CAREFUL IN EVALUATING
! THE DIFFERENCES
!
  PCH1AI = POCH1 (A+XI, -BEPS)
  PCH1I = POCH1 (XI+1.0-BEPS, BEPS)
  C0 = FACTOR * POCHAI * GAMRNI * GAMRI1 * ( &
    -POCH1(B+XI, -BEPS) + PCH1AI - PCH1I + BEPS*PCH1AI*PCH1I )
!
! XEPS1 = (1.0 - X**(-BEPS)) / BEPS
  XEPS1 = ALNX * EXPREL(-BEPS*ALNX)
!
  CHU = SUM + C0 + XEPS1*B0
  XN = N
  DO 80 I=1,1000
    XI = ISTRT + I
    XI1 = ISTRT + I - 1
    B0 = (A+XI1-BEPS)*B0*X/((XN+XI1)*(XI-BEPS))
    C0 = (A+XI1)*C0*X/((B+XI1)*XI) - ((A-1.0)*(XN+2.*XI-1.0) &
      + XI*(XI-BEPS)) * B0/(XI*(B+XI1)*(A+XI1-BEPS))
    T = C0 + XEPS1*B0
    CHU = CHU + T
    if (ABS(T) < EPS*ABS(CHU)) go to 130
 80   CONTINUE
  call XERMSG ('SLATEC', 'CHU', &
     'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES', 3, 2)
!
! X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE STRAIGHTFORWARD
! FORMULATION IS STABLE.
!
 90   A0 = FACTOR * POCHAI * GAMR(B+XI) * GAMRI1 / BEPS
  B0 = XTOEPS*B0/BEPS
!
  CHU = SUM + A0 - B0
  DO 100 I=1,1000
    XI = ISTRT + I
    XI1 = ISTRT + I - 1
    A0 = (A+XI1)*A0*X/((B+XI1)*XI)
    B0 = (A+XI1-BEPS)*B0*X/((AINTB+XI1)*(XI-BEPS))
    T = A0 - B0
    CHU = CHU + T
    if (ABS(T) < EPS*ABS(CHU)) go to 130
 100  CONTINUE
  call XERMSG ('SLATEC', 'CHU', &
     'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES', 3, 2)
!
! USE LUKE-S RATIONAL APPROX IN THE ASYMPTOTIC REGION.
!
 120  CHU = X**(-A) * R9CHU(A, B, X)
!
 130  return
end
